"""Shared DGGAL binding glue — the common code behind the DGGAL-backed
systems/ modules (isea7h, ivea7h, isea3h, ivea3h, rhealpix) and the
live-engine analyses.

DGGAL (Ecere's Discrete Global Grid Abstraction Library, `pip install dggal`,
BSD-3-Clause) exposes many DGGRSs through a single `DGGRS` API. This module
initializes the DGGAL `Application` once at import and wraps a DGGRS instance
in an `Adapter` (count / enumerate / point→cell / boundary / cid / edge
refinement / neighbors). Import it lazily (the systems/ modules do) so
`import dggs_compare` doesn't load the engine.

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

import numpy as np
import csar

from .stats import orient_ccw


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

CHUNK = 50_000              # sampling batch size (keeps memory flat)
NULL_ZONE = 0xFFFFFFFFFFFFFFFF   # DGGAL's nullZone sentinel (failed lookup)


def sample_uniform_lnglat(n, rng):
    """Uniform-on-sphere samples as (lng_deg, lat_deg), shape (n, 2)."""
    lng = 360.0 * rng.random(n) - 180.0
    lat = np.degrees(np.arcsin(2.0 * rng.random(n) - 1.0))  # equal-area in lat
    return np.column_stack([lng, lat])


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
    """Wrap one DGGAL DGGRS. `cls` is the DGGRS class name (e.g. 'ISEA7H')."""

    def __init__(self, cls):
        self.name = cls
        self.dggrs = globals()[cls]()
        # One reused query point: a GeoPoint constructed per zone_at call
        # is never freed by the binding (~1M leaked instances per sampled
        # resolution).
        self._pt = GeoPoint(0.0, 0.0)

    # ----- registry-contract batches ------------------------------------
    # The five DGGAL systems/ modules delegate their batch calls here, so
    # the loop bodies — including the (lat, lng) -> (lng, lat) argument
    # swap — live once.
    def cells_at(self, level, points):
        return [self.zone_at(level, lng, lat) for lat, lng in points]

    def cid_strs(self, zones):
        return [self.cid_str(z) for z in zones]

    def boundaries(self, zones):
        return [self.cell_boundary(z) for z in zones]

    def refined_boundaries(self, zones, refine):
        return [self.refined_boundary(z, refine) for z in zones]

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

    def verts(self, zone):
        """Corner vertices as an (M, 3) unit-vec3 array (corners only)."""
        return csar.to_vec3(self.cell_boundary(zone), geo='latlng_deg')

    # ----- ids / counts -------------------------------------------------
    def cid_str(self, zone):
        return self.dggrs.getZoneTextID(zone)

    def count(self, level):
        return int(self.dggrs.countZones(level))

    def max_level(self):
        return self.dggrs.getMaxDGGRSZoneLevel()

    def level(self, zone):
        return self.dggrs.getZoneLevel(zone)

    # ----- cell streams -------------------------------------------------
    def enumerate(self, level):
        """Every zone at `level`, whole world."""
        arr = self.dggrs.listZones(level, wholeWorld)
        try:
            yield from arr
        finally:
            Instance.delete(arr)   # the binding never frees the zone Array

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

    def sample(self, level, n, rng):
        """`n` zones from uniform-on-sphere points (with repeats; failed
        lookups skipped)."""
        done = 0
        while done < n:
            k = min(CHUNK, n - done)
            for lng, lat in sample_uniform_lnglat(k, rng):
                zone = self.zone_at(level, lng, lat)
                if zone is not None:
                    yield zone
            done += k

    def iter_sample(self, level, n, seed):
        """Yield `(cid_str, verts)` for the distinct cells among `n` samples."""
        rng = np.random.default_rng(seed)
        seen = set()
        for lng, lat in sample_uniform_lnglat(n, rng):
            zone = self.zone_at(level, lng, lat)
            if zone is None or zone in seen:
                continue
            seen.add(zone)
            yield self.cid_str(zone), self.verts(zone)
