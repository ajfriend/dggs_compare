"""Consolidated timing report for one data-release Actions run.

Combines GitHub's per-step durations (env setup, dggrid install, artifact
upload) with the in-process phase lines the pipeline prints — the
runner's `select/cids/bounds = N us/cell` generation lines and metrics'
`solve = N us/cell` lines — into one per-system table plus fleet totals.

Run with:  just timing-report <run-id>      (ids: gh run list)
Runs predating the instrumented log lines report step-level columns only.
Tiny-budget runs under-resolve sub-0.1s phases (the log lines print one
decimal), so e.g. the solve rate reads 0 there; full-budget runs are the
real subject.
"""

import json
import os
import re
import subprocess
from datetime import datetime

RUN = os.environ.get('DGGS_COMPARE_RUN', '')

# The in-process lines (runner._write_res / _write_convergence,
# metrics.build). Job logs prefix each line with "job\tstep\ttimestampZ ".
GEN_PAT = re.compile(
    r'\[\S+_r\d+ \]\s+\w+\s+(?P<n>\d+) cells \(\d+ KiB\) '
    r'\[(?P<wall>[\d.]+)s: select (?P<sel>[\d.]+) cids (?P<cid>[\d.]+) '
    r'bounds (?P<bnd>[\d.]+) =')
CONV_PAT = re.compile(r'convergence pairs written \[(?P<s>[\d.]+)s\]')
MET_PAT = re.compile(
    r'\[\S+ r\d+\s*\]\s+(?P<n>\d+) cells \(DNC \d+\) -> \S+ '
    r'\(\d+ KiB\) \[(?P<wall>[\d.]+)s: solve (?P<slv>[\d.]+) =')

# Workflow step name -> report column (anything else lands in "other").
STEP_COL = {'generate raw cell geometry': 'gen',
            'install dggrid if this implementation needs it': 'install'}
STEP_PREFIX_COL = {'metrics': 'metrics',
                   'Run actions/upload-artifact': 'upload',
                   'Run actions/checkout': 'setup',
                   'Run astral-sh/setup-uv': 'setup',
                   'Run extractions/setup-just': 'setup'}


def sh(*args):
    return subprocess.run(args, check=True, capture_output=True,
                          text=True).stdout


def step_secs(step):
    a, b = step.get('startedAt'), step.get('completedAt')
    if not a or not b:
        return 0.0
    iso = lambda s: datetime.fromisoformat(s.replace('Z', '+00:00'))
    return (iso(b) - iso(a)).total_seconds()


def col_of(name):
    if name in STEP_COL:
        return STEP_COL[name]
    for prefix, col in STEP_PREFIX_COL.items():
        if name.startswith(prefix):
            return col
    return 'other'


def parse_job(job):
    """One build job -> its report row."""
    row = {c: 0.0 for c in ('setup', 'install', 'gen', 'metrics',
                            'upload', 'other')}
    for step in job['steps']:
        row[col_of(step['name'])] += step_secs(step)
    log = sh('gh', 'run', 'view', '--job', str(job['databaseId']), '--log')
    cells = engine = solved = t_solve = conv = 0.0
    for line in log.splitlines():
        if (m := GEN_PAT.search(line)):
            cells += int(m['n'])
            engine += (float(m['sel']) + float(m['cid']) + float(m['bnd']))
        elif (m := MET_PAT.search(line)):
            solved += int(m['n'])
            t_solve += float(m['slv'])
        elif (m := CONV_PAT.search(line)):
            conv += float(m['s'])
    # None = the instrumented lines were absent (a pre-instrumentation
    # run), as opposed to a measured 0.0.
    row.update(cells=cells, engine=engine if cells else None, conv=conv,
               solved=solved, t_solve=t_solve if solved else None)
    return row


def fmt(secs):
    """Minutes; '-' for absent or sub-second (a skipped install, an
    unparsed phase)."""
    return '-' if secs is None or secs < 1 else f'{secs / 60:.1f}'


def rate(secs, n):
    return '-' if not n or secs is None else f'{1e6 * secs / n:.0f}'


def main():
    assert RUN, 'set the run id: just timing-report <run-id>'
    jobs = json.loads(sh('gh', 'run', 'view', RUN, '--json', 'jobs'))['jobs']
    builds = {j['name'][len('build ('):-1]: parse_job(j)
              for j in jobs if j['name'].startswith('build (')}

    hdr = (f'{"system":16} {"setup":>6} {"instl":>6} {"gen":>6} {"engin":>6} '
           f'{"parqt":>6} {"conv":>6} {"metrc":>6} {"solve":>6} {"upld":>6} '
           f'| {"gen":>5} {"solve":>5}')
    print('minutes' + ' ' * (len(hdr) - len('minutes') - len('us/cell')) +
          'us/cell')
    print(hdr)
    def parquet_of(r):
        if r['engine'] is None:
            return None
        return max(0.0, r['gen'] - r['engine'] - r['conv'])

    def row(label, r):
        print(f'{label:16} {fmt(r["setup"]):>6} {fmt(r["install"]):>6} '
              f'{fmt(r["gen"]):>6} {fmt(r["engine"]):>6} '
              f'{fmt(parquet_of(r)):>6} {fmt(r["conv"]):>6} '
              f'{fmt(r["metrics"]):>6} {fmt(r["t_solve"]):>6} '
              f'{fmt(r["upload"]):>6} '
              f'| {rate(r["engine"], r["cells"]):>5} '
              f'{rate(r["t_solve"], r["solved"]):>5}')

    tot = {}
    for key_, r in sorted(builds.items()):
        row(key_, r)
        for k, v in r.items():
            tot[k] = tot.get(k, 0.0) + (v or 0.0)
    if not tot['cells']:
        tot['engine'] = None
    if not tot['solved']:
        tot['t_solve'] = None
    row('TOTAL', tot)

    for j in jobs:
        if j['name'] == 'publish':
            steps = [(s['name'], step_secs(s)) for s in j['steps']
                     if step_secs(s) >= 1]
            print('\npublish: ' + '  '.join(f'{n} {fmt(s)}m'
                                            for n, s in steps))


if __name__ == '__main__':
    main()
