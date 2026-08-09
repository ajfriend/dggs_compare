"""Build the published tables from the raw geometry in data/raw/.

The binding-free half of the pipeline: reads what the implementation
scripts wrote, computes ar/area with the shared solvers, stamps each
implementation's convergence residuals into the metadata.

Run with:  uv run scripts/metrics.py
No CLI args (project convention).
"""

from dggs_compare import metrics

if __name__ == '__main__':
    metrics.build_all()
