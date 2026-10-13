# -*- coding: utf-8 -*-
"""
check_horizon.py -- the two copies of the horizon labels must agree.

00_parameters.R and 00_parameters.py each declare HORIZON_LABEL, HORIZON_SPAN,
HORIZON_SPAN_PLAIN and HORIZON_BLOCK. Two copies is one copy too many, but R
cannot read a Python module and Python will not source R, so the duplication is
paid for here: this script reads both and fails if any of the four differs.

It also fails if a PRODUCER spells a horizon out itself instead of reading the
constant. The computation is exempt: an assignment of the integer 2095 (the
block label the model actually indexes on) is what the constants exist to leave
alone, and a comment explaining the block is a comment.

Usage:  python tools/check_horizon.py
Exit:   0 if the two agree and no producer hard-codes a label, 1 otherwise.
"""

import io
import os
import re
import sys
from importlib import import_module

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
R_FILE = os.path.join(PKG, '01_estimation', '00_parameters.R')

KEYS = ('HORIZON_LABEL', 'HORIZON_SPAN', 'HORIZON_SPAN_PLAIN', 'HORIZON_BLOCK')

# Producers that must read the constants rather than spell a label out.
PRODUCERS = [
    ('02_GE', 'emit_values.py'),
    ('02_GE', 'emit_values_empirics.py'),
    ('02_GE', 'emit_figures.py'),
    ('02_GE', 'decompose_metric_channels.py'),
    ('02_GE', '02e_north_contrast_band.py'),
    ('04_contrast', 'ct3_weighted_curve_gap.py'),
    ('04_contrast', 'ct4_emit_tab_contribA_grid.R'),
    ('01_estimation', '01d_crh_triangulation.R'),
    ('03_inference', 'step6_emit.R'),
]

# A LABEL spelled out: "2015-2095" / "2015--2095" / "at 2095" / "by 2095".
# NOT matched: `FUT_YEAR <- 2095L`, `{'2095': 15}`, `level_mean_2095`, comments.
LABEL_PAT = re.compile(
    r'(?:20\d\d\s*-{1,2}\s*2095)'          # a span ending at 2095
    r'|(?:\b(?:at|by|to|in)\s+2095\b)',    # a level worded at 2095
    re.IGNORECASE)


def r_constants(path):
    src = io.open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for k in KEYS:
        m = re.search(r'^\s*' + k + r'\s*<-\s*"([^"]*)"', src, re.M)
        if m:
            out[k] = m.group(1)
    return out


def py_constants():
    sys.path.insert(0, os.path.join(PKG, '01_estimation'))
    return {k: getattr(import_module('00_parameters'), k) for k in KEYS}


def strip_comments(text, lang):
    """Drop comment and docstring text, so an explanation of the terminal block
    is not reported as a label. Blanked line by line, never deleted, so the
    line numbers a finding reports still point at the file."""
    if lang == 'py':
        def blank(m):
            return re.sub(r'[^\n]', ' ', m.group(0))
        text = re.sub(r'(?s)("""|\'\'\')(?:\\.|(?!\1).)*\1', blank, text)
    out = []
    for line in text.split('\n'):
        cut = line.find('#')
        out.append(line if cut < 0 else line[:cut])
    return '\n'.join(out)


def main():
    rc = 0

    if not os.path.exists(R_FILE):
        print('[FAIL] 00_parameters.R absent: %s' % R_FILE)
        return 1

    rk, pk = r_constants(R_FILE), py_constants()
    missing = [k for k in KEYS if k not in rk]
    if missing:
        print('[FAIL] 00_parameters.R does not declare: %s' % ', '.join(missing))
        rc = 1
    for k in KEYS:
        if k in rk and rk[k] != pk[k]:
            print('[FAIL] %s: R="%s" vs Python="%s"' % (k, rk[k], pk[k]))
            rc = 1
    if rc == 0:
        print('[OK]   the %d constants agree between R and Python' % len(KEYS))
        for k in KEYS:
            print('         %-18s = %s' % (k, pk[k]))

    n_hits = 0
    for d, f in PRODUCERS:
        p = os.path.join(PKG, d, f)
        if not os.path.exists(p):
            print('[WARN] producer missing: %s/%s' % (d, f))
            continue
        lang = 'py' if f.endswith('.py') else 'r'
        src = io.open(p, encoding='utf-8', errors='replace').read()
        for i, line in enumerate(strip_comments(src, lang).split('\n'), 1):
            if LABEL_PAT.search(line):
                print('[FAIL] %s/%s:%d hard-coded label: %s'
                      % (d, f, i, line.strip()[:90]))
                n_hits += 1
                rc = 1
    if n_hits == 0:
        print('[OK]   no producer hard-codes a horizon label')

    return rc


if __name__ == '__main__':
    sys.exit(main())
