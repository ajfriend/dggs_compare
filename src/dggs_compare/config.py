"""Pipeline configuration — the single source of truth.

These pin the data-artifact keys: every generator and consumer must agree on
them, so they live here and nowhere else. Changing any of them changes what
the tables contain, which is exactly the trigger for cutting a new data
release (see readme).
"""

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
# coverage for the full-globe viewer page (web/globe_full.html renders every
# cell; r6 is the 1.18M-cell torture test). webdata emits that page's
# binaries for exactly these entries, straight from the tables.
FULL_RES = {'ivea7h': (1, 2, 3, 5, 6)}

# ----- per-system registry metadata ---------------------------------------
# Working ("target") resolution per system: the finest in actual use, matched
# to an H3 r9 cell by scripts/calibrate.py. Systems appear here iff they have
# a module in systems/ (the folder is the registry; this is its metadata).
TARGET_RES = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10,
              'rhealpix': 9,
              # confirmed by calibrate on the data-v4 tables (1.210x h3 r9)
              'isea3h': 18, 'ivea3h': 18,
              # confirmed by calibrate on the 1M-cell tables (0.872x h3 r9)
              'isea4t': 14,
              # count-matched (12*9^r vs H3 r9, ratio 0.96); confirm with
              # calibrate once hex9 tables exist
              'hex9': 9}
# Plot color (matplotlib cycle index) per system.
SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3',
             'ivea7h': 'C4', 'rhealpix': 'C5', 'isea3h': 'C6',
             'ivea3h': 'C7', 'isea4t': 'C8', 'hex9': 'C9'}
# S2 numbers its resolutions "levels"; everyone else says "r".
RES_PREFIX = {s: 'L' if s == 's2' else 'r' for s in TARGET_RES}

# ----- globe view: area-matched resolutions -------------------------------
# Total cell count at each resolution, per system — exact closed forms. These
# grids are ~equal-area, so a cell's average area is 4*pi*R^2 / N(res); matching
# cell COUNTS across systems therefore matches average cell SIZE. This is the
# cheap, table-free way to area-match (no reading areas, no medians).
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
# Every per-system dict a new grid must appear in. The dnc-check release
# gate asserts each registry system is present in all of these, so adding
# a dict here extends the gate for free (FULL_RES stays optional;
# RES_PREFIX is derived).
PER_SYSTEM = {'TARGET_RES': TARGET_RES, 'SYS_COLOR': SYS_COLOR,
              'CELLS_PER_RES': CELLS_PER_RES}
# The globe view draws one globe per system, all at a common cell size: this H3
# resolution sets it (r3 ~ 41,162 cells ~ 12,600 km^2/cell), and every other
# system uses the resolution whose cell count is closest to H3's here. Raise for
# finer/heavier globes, lower for coarser/lighter ones.
GLOBE_H3_RES = 3

# ----- solver settings (recorded in the tables' metadata) -----------------
GAP_TOL = 1e-6            # csar duality-gap certification threshold
CSAR_METHOD = 'auto'      # solver path ('auto' = csar's recommended method)

# ----- units ---------------------------------------------------------------
R_KM = 6371.0088          # mean Earth radius; steradian -> km^2 is R^2
SR2KM2 = R_KM * R_KM
