"""Build the published tables from the raw geometry in data/raw/ —
see dggs_compare.metrics.

Run with:  uv run scripts/metrics.py
No CLI args (project convention).
"""

from dggs_compare import metrics

if __name__ == '__main__':
    metrics.build_all()
