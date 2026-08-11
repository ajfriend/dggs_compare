"""Shared DGGRID glue — the batch-subprocess engine and the GridImpl
adapter behind the DGGRID-backed scripts here. Underscore-prefixed:
shared script code, not a registry entry, not part of the library (see
_dggal_engine's docstring for the pattern).

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

Zone handles are plain SEQNUM ints — only unique within a resolution
(unlike DGGAL's level-encoding zone ints), which is why every engine call
takes `res`. Output parsing uses the AIGEN text format (id + centroid
line, one "lon lat" per vertex, closing repeat, END) — the only output
type that needs neither GDAL nor XML.
"""

import os
import shutil
import subprocess
import tempfile
from functools import cache
from pathlib import Path

from dggs_compare.stats import is_ccw

_REPO_ROOT = Path(__file__).resolve().parents[2]


@cache
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

# dggal's documented ISEA icosahedron placement: vertex 0 at the authalic
# latitude of arctan(golden ratio) N, 11.20 E, azimuth 0 — dggal's README
# states these exact DGGRID settings. Scripts whose grid must coincide
# cell-for-cell with a dggal implementation pass this as `params`;
# DGGRID's own default differs only in longitude (11.25 E). ONE home: a
# drifted copy would silently de-align a pair and read as a spurious
# rotation in `just cross-impl` — the very artifact the alignment kills.
DGGAL_ORIENTATION = {
    'dggs_vert0_lon': '11.20',
    'dggs_vert0_lat': '58.282525588538994675786',
    'dggs_vert0_azimuth': '0.0',
}


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
    the orientation test runs on the first ring only — a violation of the
    uniformity assumption would still fail loudly downstream in
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
                        flip = not is_ccw(ring)
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
    """Batch operations for one DGGRID dggs_type (e.g. 'ISEA4T').

    `params`: extra metafile lines (DGGRID vocabulary, e.g. dggs_vert0_*)
    merged into every invocation — how a script pins a non-default grid
    placement."""

    def __init__(self, dggs_type, params=None):
        self.dggs_type = dggs_type
        self._params = dict(params or {})
        # These three carry the Engine's own contract (identity, per-call
        # res, the float64 text round-trip) — a param override would be a
        # silent lie about all of them.
        assert not self._params.keys() & {'dggs_type', 'dggs_res_spec',
                                          'precision'}, self._params

    def _base(self, res):
        return {'dggs_type': self.dggs_type, 'dggs_res_spec': res,
                'precision': PRECISION, **self._params}

    def enumerate_ids(self, res):
        """Every seqnum at `res` by whole-earth generation (ids only).

        Not the pipeline path — Adapter.enumerate_cells uses the 1..N
        closed form. This is how that assumption gets re-verified when a
        new dggs_type joins: assert this equals range(1, N+1) at a few
        coarse resolutions."""
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

    def boundaries(self, res, seqnums, samples_per_edge=0):
        """Open CCW vertex lists aligned with `seqnums` (clipped
        generation); `samples_per_edge` extra vertices per edge (DGGRID's
        `densification`)."""
        with tempfile.TemporaryDirectory() as d:
            ids = Path(d) / 'ids.txt'
            ids.write_text(''.join(f'{s}\n' for s in seqnums))
            self._generate(res, d, clip='ids.txt', densification=samples_per_edge)
            rings = dict(_parse_aigen(Path(d) / 'cells.gen'))
        missing = set(seqnums) - rings.keys()
        if missing:
            raise RuntimeError(
                f'{self.dggs_type} r{res}: no boundary returned for '
                f'{len(missing)} of {len(seqnums)} cells '
                f'(e.g. {sorted(missing)[:3]})')
        return [rings[s] for s in seqnums]

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


class Adapter:
    """The GridImpl contract over one DGGRID dggs_type — each script
    passes its identity and count formula (the mirror of
    _dggal_engine.Adapter). `grid` derives from the dggs_type; `impl` is
    always 'dggrid'; `packages` is empty because the binary is not a
    python package (see the module docstring)."""

    impl = 'dggrid'
    packages = ()

    def __init__(self, dggs_type, max_res, num_cells, params=None,
                 pentagons=False):
        self.grid = dggs_type.lower()
        self.max_res = max_res
        self._num_cells = num_cells
        self._pentagons = pentagons
        self._pent_ids = {}          # res -> frozenset of pentagon seqnums
        self._engine = Engine(dggs_type, params)
        if params:
            # Non-default metafile settings change which cells exist, so
            # they are artifact identity — record them in the provenance
            # (runner merges the contract's optional metadata attr).
            self.metadata = {'engine_params': '; '.join(
                f'{k} {v}' for k, v in params.items())}
        # Zero-padded to the max seqnum's width: rows are sorted by cid
        # TEXT for spatially coherent Parquet pages, and variable-width
        # decimals would sort '10' < '9'.
        self._cid_width = len(str(num_cells(max_res)))

    def resolutions(self):
        return range(self.max_res + 1)

    def num_cells(self, res):
        return self._num_cells(res)

    def cells_at(self, res, points):
        return self._engine.cells_at(res, points)

    def cid_strs(self, cells):
        return [f'{c:0{self._cid_width}d}' for c in cells]

    def boundaries(self, res, cells, samples_per_edge=0):
        return self._engine.boundaries(res, cells, samples_per_edge)

    def enumerate_cells(self, res):
        # SEQNUMs are 1..N contiguous ascending, which is also exactly
        # the order whole-earth AIGEN output arrives in (verified for
        # ISEA3H and ISEA4T) — the closed form replaces a whole-earth
        # generation subprocess and its full-geometry text parse.
        return range(1, self._num_cells(res) + 1)

    def irregular(self, res, cells):
        """Contract declaration (interface.py): the 12 pentagons for the
        hexagonal grids (`pentagons=True`); otherwise every cell is
        regular. DGGRID's SEQNUM layout puts the pentagons at
        {1, 2, 2+q, ..., 2+10q} with q = (N-2)/10 (each block's pentagon
        comes first; degenerates to 1..12 at r0) — VERIFIED per
        resolution by checking all 12 boundaries have 5 corners, so a
        numbering change in a future DGGRID fails loudly. (A point-lookup
        alternative misfired at deep odd levels: TRANSFORM_POINTS near an
        icosahedron vertex is knife-edge once pentagons are meters wide.)"""
        if not self._pentagons:
            return [False] * len(cells)
        if res not in self._pent_ids:
            q = (self._num_cells(res) - 2) // 10
            ids = [1, 2] + [2 + k * q for k in range(1, 11)]
            nv = [len(r) for r in self._engine.boundaries(res, ids)]
            assert nv == [5] * 12, (res, ids, nv)
            self._pent_ids[res] = frozenset(ids)
        pent = self._pent_ids[res]
        return [c in pent for c in cells]
