"""Shared DGGRID glue — the batch-subprocess engine behind the
DGGRID-backed systems/ modules (isea4t).

DGGRID (Kevin Sahr's reference DGGS implementation, C++) has no supported
in-process API: the canonical interface is a metafile-driven CLI with
batch operations. That shape fits this pipeline well — we always want
cells by the hundred-thousand — and measured throughput is high (~160k
cells/s whole-earth generation; 1M points -> cells in ~4s; boundaries for
1M specific deep cells in ~6s, each in a single invocation). A fresh
process per operation also gives free immunity to the leak/degradation
issues that in-process FFI bindings can develop.

The binary: set DGGS_COMPARE_DGGRID, or drop a build at .tools/dggrid
(no PyPI wheels exist; conda-forge ships `dggrid`, or build from source
with cmake — see the justfile recipe).

Zone handles are (res, seqnum) tuples: DGGRID SEQNUM addresses are only
unique within a resolution, unlike DGGAL's level-encoding zone ints.
Output parsing uses the AIGEN text format (id + centroid line, one
"lon lat" per vertex, closing repeat, END) — the only output type that
needs neither GDAL nor XML.
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import sparea

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _find_bin():
    """DGGS_COMPARE_DGGRID env var, then .tools/dggrid, then PATH."""
    env = os.environ.get('DGGS_COMPARE_DGGRID')
    if env:
        return Path(env)
    tools = _REPO_ROOT / '.tools' / 'dggrid'
    if tools.exists():
        return tools
    on_path = shutil.which('dggrid')
    if on_path:
        return Path(on_path)
    raise FileNotFoundError(
        'no dggrid binary found: set DGGS_COMPARE_DGGRID, or run '
        '`just install-dggrid` (builds into .tools/)')

# AIGEN emits this many decimal degrees digits; 12 keeps float64 fidelity
# through the text round-trip.
PRECISION = 12


def _run(params, workdir):
    """Write a metafile from `params` and run dggrid on it in `workdir`."""
    meta = Path(workdir) / 'job.meta'
    meta.write_text(''.join(f'{k} {v}\n' for k, v in params.items()))
    proc = subprocess.run([str(_find_bin()), str(meta)],
                          capture_output=True, text=True, cwd=workdir)
    if proc.returncode != 0:
        raise RuntimeError(
            f'dggrid failed (exit {proc.returncode}) for {params!r}:\n'
            f'{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}')


def _parse_aigen(path):
    """Yield (seqnum, [(lat, lng), ...]) open CCW rings from an AIGEN file.

    Winding is uniform within one DGGRID output (measured: always CCW), so
    the sparea orientation test runs on the first ring only — a violation
    of the uniformity assumption would still fail loudly downstream in
    stats.cell_stats' area-over-2pi guard.
    """
    flip = None
    with open(path) as f:
        cur, ring = None, []
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == 'END':
                if cur is not None:
                    if len(ring) > 1 and ring[0] == ring[-1]:
                        ring = ring[:-1]
                    if flip is None:
                        flip = sparea.area(ring, signed=True) < 0
                    yield cur, ring[::-1] if flip else ring
                cur, ring = None, []
            elif cur is None:
                cur = int(t[0])              # "seqnum centroid_lon centroid_lat"
            else:
                ring.append((float(t[1]), float(t[0])))   # lon lat -> (lat, lng)


def _parse_aigen_ids(path):
    """Yield only the seqnums from an AIGEN file (~15x cheaper than the
    full parse when the rings would be discarded)."""
    with open(path) as f:
        expect_header = True
        for line in f:
            if expect_header:
                t = line.split(None, 1)
                if t and t[0] != 'END':
                    yield int(t[0])
                    expect_header = False
            elif line.startswith('END'):
                expect_header = True


class Engine:
    """Batch operations for one DGGRID dggs_type (e.g. 'ISEA4T')."""

    def __init__(self, dggs_type):
        self.dggs_type = dggs_type

    def _base(self, res):
        return {'dggs_type': self.dggs_type, 'dggs_res_spec': res,
                'precision': PRECISION}

    def enumerate_ids(self, res):
        """Every seqnum at `res` (whole-earth generation, ids only)."""
        with tempfile.TemporaryDirectory() as d:
            self._generate(res, d, clip=None)
            yield from _parse_aigen_ids(Path(d) / 'cells.gen')

    def cells_at(self, res, latlngs):
        """seqnum (or None) for each (lat, lng) — one TRANSFORM_POINTS run."""
        with tempfile.TemporaryDirectory() as d:
            pts = Path(d) / 'points.txt'
            with open(pts, 'w') as f:
                for lat, lng in latlngs:
                    f.write(f'{lng:.{PRECISION}f} {lat:.{PRECISION}f}\n')
            _run({**self._base(res),
                  'dggrid_operation': 'TRANSFORM_POINTS',
                  'input_file_name': 'points.txt',
                  'input_address_type': 'GEO',
                  'input_delimiter': '" "',
                  'output_file_name': 'ids.txt',
                  'output_address_type': 'SEQNUM',
                  'output_delimiter': '" "'}, d)
            out = []
            with open(Path(d) / 'ids.txt') as f:
                for line in f:
                    line = line.strip()
                    out.append(int(line) if line else None)
            if len(out) != len(latlngs):
                raise RuntimeError(
                    f'{self.dggs_type} r{res}: {len(latlngs)} points in, '
                    f'{len(out)} ids out')
            return out

    def boundaries(self, res, seqnums, refine=0):
        """{seqnum: open CCW ring} for exactly `seqnums` (clipped
        generation); `refine` densification points per edge."""
        with tempfile.TemporaryDirectory() as d:
            ids = Path(d) / 'ids.txt'
            ids.write_text(''.join(f'{s}\n' for s in seqnums))
            self._generate(res, d, clip='ids.txt', densification=refine)
            rings = dict(_parse_aigen(Path(d) / 'cells.gen'))
        missing = set(seqnums) - rings.keys()
        if missing:
            raise RuntimeError(
                f'{self.dggs_type} r{res}: no boundary returned for '
                f'{len(missing)} of {len(seqnums)} cells '
                f'(e.g. {sorted(missing)[:3]})')
        return rings

    def _generate(self, res, workdir, clip, densification=0):
        params = {**self._base(res),
                  'dggrid_operation': 'GENERATE_GRID',
                  'cell_output_type': 'AIGEN',
                  'cell_output_file_name': 'cells',
                  'densification': densification}
        if clip is None:
            params['clip_subset_type'] = 'WHOLE_EARTH'
        else:
            params.update({'clip_subset_type': 'INPUT_ADDRESS_TYPE',
                           'input_address_type': 'SEQNUM',
                           'clip_region_files': clip})
        _run(params, workdir)
