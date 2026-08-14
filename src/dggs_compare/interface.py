"""The implementation contract: what a scripts/systems/ script provides.

A script defines one implementation of one grid and hands an instance to
`runner.generate`. Implementations are identified by the pair
`(grid, impl)` — a grid may have several implementations (isea3h via
dggal and via DGGRID), and artifacts are keyed `{grid}-{impl}`.

A cell's boundary is always the same kind of object — a list of vertices,
in order, joined by great-circle edges — and the only knob is the SAMPLING
DENSITY: how many extra vertices are sampled along each native edge.
Density 0 is the minimal vertex set. Bringing cells to the sphere is the
implementation's job, done however it thinks best (the closed-form
authalic map is available for implementations that address WGS84 geodetic
coordinates: `stats.authalic_rings` for a `boundaries()` return value,
`stats.authalic_lat` for bare latitudes); the script is the record of how
its system got to the sphere.

Cell handles are plain hashables whose meaning may be scoped to a
resolution (DGGRID seqnums), so every geometry call carries `res`.
"""

from typing import Iterator, Protocol, Sequence


class GridImpl(Protocol):
    grid: str                   # the grid, e.g. 'isea3h'
    impl: str                   # the implementation, e.g. 'dggal'
    packages: Sequence[str]     # distributions whose versions to record
    # OPTIONAL: extra provenance strings merged into the raw files'
    # metadata (and forwarded to the published tables if named in
    # runner.META_KEYS). For identity-affecting engine settings the
    # script chooses — e.g. a pinned DGGRID orientation — which must be
    # readable from the artifact, not only from script source.
    # metadata: dict[str, str]

    def resolutions(self) -> range:
        """Generatable resolutions (0..finest)."""
        ...

    def num_cells(self, res: int) -> int:
        """Total cells at `res`."""
        ...

    def cells_at(self, res: int, points) -> list:
        """[cell or None, ...] for [(lat, lng), ...] degrees. None = the
        engine couldn't resolve the point (rare deep-level singular
        points); samplers skip and draw again."""
        ...

    def cid_strs(self, cells) -> list[str]:
        """Text ids stored in the tables. Must sort (as text) in a
        spatially coherent order — rows are ordered by them for Parquet
        page locality (fixed-width)."""
        ...

    def boundaries(self, res: int, cells, samples_per_edge: int = 0) -> list:
        """[[(lat, lng), ...] open vertex list, ...] aligned with `cells`,
        degrees, with AT LEAST `samples_per_edge` extra vertices sampled
        along each native edge (binding APIs quantize, so implementations
        may return more)."""
        ...

    def enumerate_cells(self, res: int) -> Iterator:
        """Every cell at `res`."""
        ...
