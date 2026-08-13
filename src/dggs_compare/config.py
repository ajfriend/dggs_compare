"""Pipeline configuration — the single source of truth.

These pin the data-artifact keys: every generator and consumer must agree on
them, so they live here and nowhere else. Changing any of them changes what
the tables contain, which is exactly the trigger for cutting a new data
release (see readme).
"""

import math

# ----- sampling ----------------------------------------------------------
SEED = 0xDECAF   # data-v3+ (data-v1/v2 used 0xC0FFEE); changing the seed is
                 # a data-release trigger — fresh sample points everywhere,
                 # which doubles as a hunt for more DGGAL nullZone examples
# Sampling-seed policy. True (default): a distinct stream per resolution
# (default_rng([SEED, res])), so every resolution draws INDEPENDENT points and
# cross-resolution agreement is honest convergence. False: every resolution
# reuses the SAME points, which reads out a scale-invariant field identically
# across resolutions but partly manufactures that agreement. Deterministic
# either way; changes which cells are cached, not the schema.
PER_RES_SEED = True
# Per-resolution cell budget — EXACTLY this many cells at every sampled
# resolution (coarse resolutions with fewer enumerate every cell). Selection
# is three-regime, by total cell count N at the resolution:
#   N <= N_CELLS                      enumerate all
#   N <= SUBSAMPLE_MAX_RATIO*N_CELLS  enumerate ids, subsample to N_CELLS
#                                     (avoids the coupon-collector blowup of
#                                     point-sampling near the cap)
#   larger                            draw uniform points until N_CELLS
#                                     distinct cells (~1.1x draws out here)
# Nuance, negligible for ~equal-area grids: point-sampling includes a cell
# with probability ~proportional to its area; subsampling is uniform per
# cell.
N_CELLS = 1_000_000
SUBSAMPLE_MAX_RATIO = 4
# Resolutions enumerated EXHAUSTIVELY regardless of the budgets — complete
# coverage in the published tables for these entries (r6 is the 1.18M-cell
# torture test); the runner is the consumer.
FULL_RES = {'ivea7h': (1, 2, 3, 5, 6)}

# ----- per-grid registry metadata ------------------------------------------
# Working ("target") resolution per grid: the finest in actual use, count-
# matched to an H3 r9 cell by scripts/calibrate.py (closed-form, no tables;
# it flags any entry here that disagrees with its pick). Grids appear here
# iff they have an implementation script in scripts/systems/ (the listing
# is the registry; this is its per-grid metadata).
TARGET_RES = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10,
              'rhealpix': 9, 'isea3h': 18, 'ivea3h': 18, 'isea4t': 14,
              'hex9': 9}
# Which implementation's tables the site and analyses read, per grid. A
# grid may have several implementations ({grid}-{impl} artifact keys);
# consumers key by grid and resolve through this dict.
PRIMARY_IMPL = {'h3': 'h3', 's2': 's2', 'a5': 'a5', 'isea7h': 'dggal',
                'ivea7h': 'dggal', 'rhealpix': 'dggal', 'isea3h': 'dggal',
                'ivea3h': 'dggal', 'isea4t': 'dggrid', 'hex9': 'hex9'}
# Expected `irregular` count per resolution, per grid — the check against
# each implementation's DECLARED exceptional cells (the tables' irregular
# column; interface.py). 12 = a hexagonal grid's pentagons; 0 = no
# exceptional class (a5's pentagons are its regular cell; hex9's lattice
# defects are tiling vertices, not cells). The dnc-check gate asserts the
# count exactly at full-enumeration resolutions, so a new hex script that
# forgets its declaration fails loudly instead of publishing an all-False
# column indistinguishable from a legitimate declares-none grid.
EXPECTED_IRREGULAR = {'h3': 12, 's2': 0, 'a5': 0, 'isea7h': 12,
                      'ivea7h': 12, 'rhealpix': 0, 'isea3h': 12,
                      'ivea3h': 12, 'isea4t': 0, 'hex9': 0}
# Plot color (matplotlib cycle index) per system.
SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3',
             'ivea7h': 'C4', 'rhealpix': 'C5', 'isea3h': 'C6',
             'ivea3h': 'C7', 'isea4t': 'C8', 'hex9': 'C9'}
# S2 numbers its resolutions "levels"; everyone else says "r".
RES_PREFIX = {s: 'L' if s == 's2' else 'r' for s in TARGET_RES}

# ----- cell counts: the area-matching currency ----------------------------
# Total cell count at each resolution, per system — exact closed forms. Cells
# partition the sphere, so a cell's average area is exactly 4*pi*R^2 / N(res);
# matching cell COUNTS across systems therefore matches average cell SIZE with
# no tables read (no areas, no medians). Both calibrate (TARGET_RES picks) and
# the globe view below match this way.
CELLS_PER_RES = {
    'h3':       lambda r: 2 + 120 * 7 ** r,               # 122, 842, 5882, 41162, ...
    's2':       lambda r: 6 * 4 ** r,                     # 6, 24, 96, ...
    'a5':       lambda r: 12 if r == 0 else 15 * 4 ** r,  # 12, 60, 240, 960, ...
    'isea7h':   lambda r: 10 * 7 ** r + 2,                # 12, 72, 492, 3432, 24012, ...
    'ivea7h':   lambda r: 10 * 7 ** r + 2,
    'isea3h':   lambda r: 10 * 3 ** r + 2,                # 12, 32, 92, 272, 812, ...
    'ivea3h':   lambda r: 10 * 3 ** r + 2,
    'isea4t':   lambda r: 20 * 4 ** r,                    # 20, 80, 320, 1280, ...
    'rhealpix': lambda r: 6 * 9 ** r,                     # 6, 54, 486, 4374, 39366, ...
    'hex9':     lambda r: 12 * 9 ** r,                    # 12, 108, 972, 8748, ...
}
# THE matching criterion (issue #32): the resolution whose cell count is
# closest to `anchor_n` in log-ratio — the symmetric size mismatch (2x too
# big == 2x too small); a plain count difference would bias coarse. Both
# calibrate (TARGET_RES, the dnc-check gate asserts agreement) and the
# globe view use this. The default candidate range is generous: even the
# slowest-growing family (3^r) crosses any sane anchor well inside it.
def count_match_res(sys, anchor_n, resolutions=range(40)):
    n_of = CELLS_PER_RES[sys]
    return min(resolutions, key=lambda r: abs(math.log(n_of(r) / anchor_n)))


def mean_cell_area(sys, res):
    """Exact mean cell area at `res` in steradians: cells partition the
    sphere, so 4*pi / N(res) — closed form, no tables read. THE
    normalizer for relative-area statistics (survey, web viewer)."""
    return 4.0 * math.pi / CELLS_PER_RES[sys](res)


def sampling_regime(sys, res, n_cells=None):
    """'all' | 'subsam' | 'sample' — the selection regime the runner uses
    for (sys, res) at budget `n_cells` (default N_CELLS; pass a table's
    STAMPED n_cells to classify override-budget builds, where FULL_RES
    is also disabled). The closed-form home of the three-regime boundary
    described at N_CELLS above; runner's selection implements it."""
    n = N_CELLS if n_cells is None else n_cells
    if n == N_CELLS and res in FULL_RES.get(sys, ()):
        return 'all'
    total = CELLS_PER_RES[sys](res)
    if total <= n:
        return 'all'
    if total <= SUBSAMPLE_MAX_RATIO * n:
        return 'subsam'
    return 'sample'


# Every per-system dict a new grid must appear in. The dnc-check release
# gate asserts each registry system is present in all of these, so adding
# a dict here extends the gate for free (FULL_RES stays optional;
# RES_PREFIX is derived).
PER_SYSTEM = {'TARGET_RES': TARGET_RES, 'SYS_COLOR': SYS_COLOR,
              'CELLS_PER_RES': CELLS_PER_RES, 'PRIMARY_IMPL': PRIMARY_IMPL,
              'EXPECTED_IRREGULAR': EXPECTED_IRREGULAR}
# The globe view offers a Resolution menu of H3 anchor levels 0..this, every
# system drawing at the resolution whose cell count is closest to the selected
# anchor's. This sets the FINEST anchor, which is also the menu default
# (r3 ~ 41,162 cells ~ 12,600 km^2/cell). Raising it adds finer menu entries
# whose payloads grow geometrically; the coarser levels are nearly free.
GLOBE_H3_RES = 3

# ----- csar settings (recorded in the tables' metadata) -------------------
GAP_TOL = 1e-6            # csar duality-gap certification threshold
CSAR_METHOD = 'auto'      # solver path ('auto' = csar's recommended method)

# ----- convergence check (pins the convergence artifact's content) ---------
# Density-0 vs dense sampling residuals, measured on K cells per gate level.
CONV_SEED = 0xC0FFEE
CONV_LEVELS = (0, 1, 2, 3, 5, 8, 11)   # incl coarsest + the pentagons
CONV_K = 300
CONV_SAMPLES = 40         # dense-reference sampling density (per edge)
CONV_TOL = 1e-3           # max |dAR| beyond this fails the admission gate
# Grids allowed to exceed CONV_TOL, with the documented reason (the gate
# prints it): isea3h's odd-level density-0 approximation, issue #25.
CONV_EXPECTED_RED = {'isea3h': 'density-0 approximation, issue #25'}

# ----- units ---------------------------------------------------------------
R_KM = 6371.0088          # mean Earth radius; steradian -> km^2 is R^2
SR2KM2 = R_KM * R_KM
