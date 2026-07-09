"""Pipeline configuration — the single source of truth.

These pin the data-artifact keys: every generator and consumer must agree on
them, so they live here and nowhere else. Changing any of them changes what
the tables contain, which is exactly the trigger for cutting a new data
release (see readme).
"""

# ----- sampling ----------------------------------------------------------
SEED = 0xC0FFEE
# Sampling-seed policy. True (default): a distinct stream per resolution
# (default_rng([SEED, res])), so every resolution draws INDEPENDENT points and
# cross-resolution agreement is honest convergence. False: every resolution
# reuses the SAME points, which reads out a scale-invariant field identically
# across resolutions but partly manufactures that agreement. Deterministic
# either way; changes which cells are cached, not the schema.
PER_RES_SEED = True
# Per-resolution cell budget: dense N_BIG up to a system's working (target)
# resolution — where the survey and explorations want lots of cells — and a
# thin N_SMALL beyond it, where only the DNC/calibration tails read. Both are
# capped by the resolution's total cell count (coarse resolutions enumerate
# every cell, exactly).
N_BIG, N_SMALL = 100_000, 25_000

# ----- per-system registry metadata ---------------------------------------
# Working ("target") resolution per system: the finest in actual use, matched
# to an H3 r9 cell by scripts/calibrate.py. Systems appear here iff they have
# a module in systems/ (the folder is the registry; this is its metadata).
TARGET_RES = {'h3': 9, 's2': 15, 'a5': 14, 'isea7h': 10, 'ivea7h': 10,
              'rhealpix': 9}
# Plot color (matplotlib cycle index) per system.
SYS_COLOR = {'h3': 'C0', 's2': 'C1', 'a5': 'C2', 'isea7h': 'C3',
             'ivea7h': 'C4', 'rhealpix': 'C5'}
# S2 numbers its resolutions "levels"; everyone else says "r".
RES_PREFIX = {s: 'L' if s == 's2' else 'r' for s in TARGET_RES}

# ----- solver settings (recorded in the tables' metadata) -----------------
GAP_TOL = 1e-6            # skar duality-gap certification threshold
SKAR_METHOD = 'auto'      # solver path ('auto' = skar's recommended method)

# ----- units ---------------------------------------------------------------
R_KM = 6371.0088          # mean Earth radius; steradian -> km^2 is R^2
SR2KM2 = R_KM * R_KM
