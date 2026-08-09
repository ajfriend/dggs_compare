"""Report every implementation's convergence residuals.

The residuals (max |dAR|, density-0 vs dense vertex sampling, per gate
level) are measured by the metrics stage at build time and stamped into
each final table's metadata; the gate itself (config.CONV_TOL, with the
config.CONV_EXPECTED_RED carve-outs) is enforced there. This script just
reads the stamped values back out of the published tables — a pure
metadata read, no bindings, no measurement.

Run with:  just convergence
No CLI args (project convention).
"""

import json

from dggs_compare import cache, config


def main():
    for (grid, impl), res_list in sorted(cache.available_tables().items()):
        raw = cache.table_metadata(grid, impl, res_list[0]).get(
            'convergence_max_dar')
        if raw is None:
            print(f'{cache.key(grid, impl)}: no convergence metadata '
                  f'(pre-#37 table?)')
            continue
        residuals = {int(r): v for r, v in json.loads(raw).items()}
        worst = max(residuals.values(), default=0.0)
        ladder = '  '.join(f'r{r}={v:.1e}' for r, v in sorted(residuals.items()))
        note = ''
        if worst >= config.CONV_TOL:
            note = f'  [over tolerance: {config.CONV_EXPECTED_RED.get(grid, "UNEXPECTED")}]'
        print(f'{cache.key(grid, impl)}: {ladder}{note}')


if __name__ == '__main__':
    main()
