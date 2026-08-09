"""Build the web viewer's static data (histograms + globe binaries +
manifest) from the Parquet tables — a column reshape, no measurement.

Run with:  just web-data
No CLI args (project convention) — knobs live in dggs_compare.webdata.
"""

from dggs_compare import webdata

if __name__ == '__main__':
    webdata.build_all()
