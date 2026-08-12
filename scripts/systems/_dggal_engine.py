"""Shared DGGAL binding glue — the common code behind the DGGAL-backed
scripts here (the `*-dggal.py` listing).

System-specific code lives OUT of the library, next to the scripts that
use it: `uv run scripts/systems/<key>.py` puts this directory on
sys.path, so the scripts import this module directly. The underscore
prefix keeps it out of the registry (the {grid}-{impl}.py listing) and
cannot collide with the `dggal` package itself.

DGGAL (Ecere's Discrete Global Grid Abstraction Library, `pip install dggal`,
BSD-3-Clause) exposes many DGGRSs through a single `DGGRS` API. This module
initializes the DGGAL `Application` once at import and wraps a DGGRS
instance in an `Adapter` speaking the GridImpl contract (interface.py).
Import it lazily (the scripts do) so `import dggs_compare` doesn't load
the engine.

Platform note: dggal 0.0.6's macOS arm64 wheel is still half-broken — the
Python extensions are arm64 but the bundled libecrt/libdggal dylibs are
x86_64 — so on Apple Silicon the project env runs x86_64 under Rosetta (see
the justfile `sync` recipe; one env, no re-exec tricks). Linux wheels are
correct. The guarded dlopen below is belt-and-suspenders for wheels with
missing dylib load commands.
"""

import ctypes
import glob
import importlib.util
import os

import csar
import numpy as np

from dggs_compare.stats import authalic_lat, authalic_rings, orient_ccw


def _geodetic_lat(xi):
    """Geodetic latitude whose authalic latitude is `xi` degrees
    (bisection; the map is monotone; ~1e-13 deg after 60 halvings)."""
    lo, hi = -90.0, 90.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if authalic_lat(mid) < xi:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _unit_mean(ring):
    """Normalized mean of an open (lat, lng)-degree ring's unit vectors,
    as (lat, lng) degrees — for a ring symmetric about an axis, the axis."""
    v = csar.to_vec3(ring, geo='latlng_deg').mean(0)
    v /= np.linalg.norm(v)
    return (float(np.degrees(np.arcsin(v[2]))),
            float(np.degrees(np.arctan2(v[1], v[0]))))


def _preload_native():
    """dlopen libecrt/libdggal RTLD_GLOBAL so flat-namespace symbols resolve."""
    for pkg, stem in (('ecrt', 'libecrt'), ('dggal', 'libdggal')):
        spec = importlib.util.find_spec(pkg)
        if spec is None or not spec.origin:
            continue
        libdir = os.path.join(os.path.dirname(spec.origin), 'lib')
        for lib in sorted(glob.glob(os.path.join(libdir, stem + '.*'))):
            try:
                ctypes.CDLL(lib, mode=ctypes.RTLD_GLOBAL)
            except OSError:
                pass


try:
    import dggal as _dggal  # noqa: F401
except ImportError:
    _preload_native()
    import dggal as _dggal  # noqa: F401

from dggal import *  # noqa: E402,F401,F403  upstream-documented setup pattern

_app = Application(appGlobals=globals())
pydggal_setup(_app)

NULL_ZONE = 0xFFFFFFFFFFFFFFFF   # DGGAL's nullZone sentinel (failed lookup)


def _drain(points):
    """Copy a DGGAL vertex Array out as [(lat, lng), ...], then FREE it.

    The dggal 0.0.6 binding never frees the eC Array a vertex call returns
    (~0.5 KB per corners call, ~3.6 KB per refined call). At a million
    boundary calls per table that OOMs the data-release runners — isea3h's
    refined odd levels were the first to hit it — so every Array is deleted
    here the moment its floats are copied out.
    """
    try:
        return [(float(p.lat), float(p.lon)) for p in points]
    finally:
        # Instance.delete explicitly: on Arrays, `.delete` resolves to eC
        # Container's per-ELEMENT delete and would fail wanting an argument.
        Instance.delete(points)


def latlng_ring(points):
    """DGGAL WGS84 vertex points (`.lat`/`.lon` Degrees) -> [(lat, lng), ...].

    Frees the passed Array (see _drain). Strips a closing repeat if present
    (matches the H3/S2/A5 adapters; handles hexagons and the 12 pentagons).
    """
    ring = _drain(points)
    if len(ring) >= 2 and ring[0] == ring[-1]:
        ring = ring[:-1]
    return ring


# DGGAL is inconsistent about corner winding: the 7H grids and rHEALPix
# come back CCW but ISEA3H/IVEA3H come back CW (dggal 0.0.6), except two
# pentagons per level. Every ring an Adapter hands out is normalized
# (stats.orient_ccw).
class Adapter:
    """Wrap one DGGAL DGGRS as the GridImpl the DGGAL scripts hand to
    `runner.generate`.

    `cls` is the DGGRS class name (e.g. 'ISEA7H'); the artifact key's
    grid half is its lowercase (override with `grid` if they ever
    diverge). DGGAL's coordinates are WGS84 geodetic, so the registry
    scripts pass `to_sphere` (a rings->rings function, e.g.
    `stats.authalic_rings`) — its presence in the script call is the
    record of how the system gets to the sphere. The contract methods
    below return that frame; the per-cell singles after them (their
    building blocks) return the binding's raw geodetic frame.
    """

    impl = 'dggal'
    packages = ('dggal',)

    def __init__(self, cls, to_sphere=None, grid=None, pentagons=False):
        self.name = cls
        self.grid = cls.lower() if grid is None else grid
        self._to_sphere = to_sphere
        self._pentagons = pentagons
        self._pent_centers = None    # geodetic (lat, lng) of the 12 vertices
        self._pent_zones = {}        # res -> frozenset of pentagon zones
        self.dggrs = globals()[cls]()
        # One reused query point: a GeoPoint constructed per zone_at call
        # is never freed by the binding (~1M leaked instances per sampled
        # resolution).
        self._pt = GeoPoint(0.0, 0.0)

    # ----- the GridImpl contract ----------------------------------------
    # ("zone" below is DGGAL's own name for a cell; it stays inside this
    # adapter. Zones encode their level, so `res` is implied by the cell
    # arguments where the contract passes both.)
    def resolutions(self):
        return range(self.dggrs.getMaxDGGRSZoneLevel() + 1)

    def num_cells(self, res):
        return int(self.dggrs.countZones(res))

    def cells_at(self, res, points):
        return [self.zone_at(res, lng, lat) for lat, lng in points]

    def cid_strs(self, cells):
        return [self.cid_str(c) for c in cells]

    def boundaries(self, res, cells, samples_per_edge=0):
        if samples_per_edge:
            rings = [self.refined_boundary(c, samples_per_edge)
                     for c in cells]
        else:
            rings = [self.cell_boundary(c) for c in cells]
        return self._to_sphere(rings) if self._to_sphere else rings

    def enumerate_cells(self, res):
        """Every zone at `res`, whole world."""
        arr = self.dggrs.listZones(res, wholeWorld)
        try:
            yield from arr
        finally:
            Instance.delete(arr)   # the binding never frees the zone Array

    def irregular(self, res, cells):
        """Contract declaration (interface.py): the 12 pentagons — the
        cells centered on the icosahedron vertices — for the hexagonal
        grids (`pentagons=True`); otherwise every cell is regular."""
        if not self._pentagons:
            return [False] * len(cells)
        if res not in self._pent_zones:
            if self._pent_centers is None:
                # The r0 cells ARE the 12 pentagons, each 5-fold symmetric
                # about its icosahedron vertex ON THE SPHERE — so the
                # unit-vector mean of the sphere-frame corners IS the
                # vertex (to numerics), mapped back to the geodetic frame
                # zone_at expects. Exactness matters: a deep-level
                # pentagon is centimeters wide.
                self._pent_centers = []
                for z in self.enumerate_cells(0):
                    ring = authalic_rings([self.cell_boundary(z)])[0]
                    lat, lng = _unit_mean(ring)
                    self._pent_centers.append((_geodetic_lat(lat), lng))
                assert len(self._pent_centers) == 12, self._pent_centers
            zones = self.cells_at(res, self._pent_centers)
            # Self-verifying declaration: all 12 must actually be
            # pentagons (guards vertex-lookup drift at depths where a
            # pentagon is centimeters wide, and future dggal versions).
            assert all(len(self.cell_boundary(z)) == 5 for z in zones), zones
            # int(zone) is the canonical id: listZones yields typed
            # (unhashable) zone objects while zone_at returns plain ints.
            pent = frozenset(int(z) for z in zones)
            assert len(pent) == 12, zones     # 12 lookups, 12 DISTINCT hits
            self._pent_zones[res] = pent
        pent = self._pent_zones[res]
        return [int(c) in pent for c in cells]

    # ----- geometry -----------------------------------------------------
    def cell_boundary(self, zone):
        """Corner vertices of `zone` as [(lat, lng), ...] deg, open ring."""
        ring = _drain(self.dggrs.getZoneWGS84Vertices(zone))
        # DGGAL's corner method collapses 2 of 4 corners to (0, 0) for some
        # rHEALPix polar-cap/equatorial-seam cells (a dggal bug). The collapse
        # shows up as duplicate vertices; fall back to the edge-refined
        # boundary, which traces the real cell correctly.
        if len(set(ring)) < len(ring):
            ring = _drain(self.dggrs.getZoneRefinedWGS84Vertices(zone, 0))
        if len(ring) >= 2 and ring[0] == ring[-1]:
            ring = ring[:-1]
        if len(set(ring)) < 3:
            # Never let garbage geometry reach a table silently (a nullZone's
            # "boundary" is all-(0,0); see zone_at).
            raise ValueError(
                f'{self.name}: degenerate boundary for zone '
                f'{self.cid_str(zone)!r} ({len(ring)} verts)')
        return orient_ccw(ring)

    def refined_boundary(self, zone, refine):
        """Edge-refined [(lat, lng), ...] open ring (`refine` points/edge)."""
        return orient_ccw(
            latlng_ring(self.dggrs.getZoneRefinedWGS84Vertices(zone, refine)))

    # ----- ids / point lookup -------------------------------------------
    def cid_str(self, zone):
        return self.dggrs.getZoneTextID(zone)

    def zone_at(self, level, lng, lat):
        """The zone at `level` containing the (lng, lat) point, or None when
        the engine can't resolve it — DGGAL returns its nullZone sentinel for
        rare singular points at deep levels (observed ~1 per 1M uniform
        draws at isea7h r15, near an icosahedron edge). Callers sampling
        points should skip None and draw again."""
        p = self._pt
        p.lat, p.lon = float(lat), float(lng)
        zone = self.dggrs.getZoneFromWGS84Centroid(level, p)
        return None if zone == NULL_ZONE else zone
