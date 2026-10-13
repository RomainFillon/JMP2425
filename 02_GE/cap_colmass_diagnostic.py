#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
cap_colmass_diagnostic.py -- column-mass diagnostic of lambda_0 (trade) and
mu_0 (migration), for the calibration paragraph of app:ge_robustness.

Reads:  intermediate_outputs/{lambda.h5, mig_flows.h5, indices_drop.csv}
Writes: output/cap_colmass_diagnostic.csv
Emits (through 02d_emit_latex.R -> emit_values.py): CALIBlamMaxColmass,
        CALIBlamNcolsExceed, CALIBmuMaxColmass, CALIBmuNcolsExceed.

--------------------------------------------------------------------------
WHY THIS SCRIPT EXISTS
--------------------------------------------------------------------------
The manuscript states that the degeneracy the column-mass cap repairs is
specific to trade: mu_0, calibrated on the same gravity principle and entering
the same counterfactual, does not exhibit it. Four numbers carry that claim.
Typed by hand they would be four literals; here they are re-derived from the
two matrices the solver itself loads.

--------------------------------------------------------------------------
THE TRAP THIS FILE IS ABOUT: mu_0 AND lambda_0 ARE NOT LOADED THE SAME WAY
--------------------------------------------------------------------------
  main.py l.299  migration0 = load_migration(mig_flows.h5, indices_drop)
  main.py l.304  lambda0    = load_h5(lambda.h5, indices_drop, "lambda", cap)
`load_migration` does not purge values below 1e-15 and forces the diagonal to
one on empty rows; `load_h5` purges and forces nothing. The column-mass probe
lives inside `load_h5`, so it is NEVER run on mu_0 by the pipeline. The two
halves of the manuscript's sentence bear on objects built differently, and that
is reproduced here rather than corrected.

--------------------------------------------------------------------------
NO RE-IMPLEMENTATION, AND THAT IS THE GUARDRAIL
--------------------------------------------------------------------------
An earlier version of this probe (04/08/2026, outside the package) re-coded the
two loading paths by hand and got them wrong -- 3862 / 487 against 2378.713 /
423 -- and only a hard-coded control value caught it. That control cannot
travel into the package: the package regenerates lambda.h5 from its own
producer, so a value certified against another artefact would refuse every
legitimate run.

This version removes the need for it. It IMPORTS `load_h5` and
`load_migration` from 02c_model/data_cleaning.py and calls them. The probe is
the shipped code by construction, so there is nothing left to re-verify: any
divergence between this diagnostic and the solver would be a divergence of the
solver with itself.

What remains checkable is that the two agree in the run, and they do: `load_h5`
prints `lambda colmass AVANT cap: max=..., n_cols>1.5=...` on every GE run, and
those two numbers are the two this script writes for lambda_0. Cross-check
against any logs/run_*/ge_*.log.

LEAD-ONLY on the inputs. Writes output/cap_colmass_diagnostic.csv and nothing
else.

Usage:  python cap_colmass_diagnostic.py [--dry-run]
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))            # .../02_GE/
PKG = os.environ.get('JMP_ROOT') or os.path.dirname(HERE)    # .../replication/
INPUTS = os.path.join(PKG, 'intermediate_outputs')
OUT = os.path.join(PKG, 'output')
CSV = os.path.join(OUT, 'cap_colmass_diagnostic.csv')

## The solver's loaders are imported, not copied. See the header.
sys.path.insert(0, os.path.join(HERE, '02c_model'))
from data_cleaning import load_h5, load_migration   # noqa: E402

DRY = '--dry-run' in sys.argv

## The bound the manuscript names. It is the threshold the two count columns
## are counted against; it is NOT applied here -- this script measures the
## uncapped matrices, which is what the sentence describes.
C_BOUND = 1.5


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    need = ['lambda.h5', 'mig_flows.h5', 'indices_drop.csv']
    missing = [n for n in need if not os.path.exists(os.path.join(INPUTS, n))]
    if missing:
        print('[STOP] intrant(s) absent(s) de intermediate_outputs/ : %s'
              % ', '.join(missing))
        print('       Producteurs : 02a_calibration_at_1deg.r (lambda.h5,')
        print('       indices_drop.csv) et 02b_mig_flows.R (mig_flows.h5).')
        print('       Aucun artefact ecrit.')
        return 2

    ## The md5 are RECORDED, not pinned. Pinning would refuse every run whose
    ## upstream was legitimately regenerated; recording lets a reader tie the
    ## four published numbers to the exact bytes they came from.
    stamps = {n: md5(os.path.join(INPUTS, n)) for n in need}
    print('[md5] intrants')
    for n, h in stamps.items():
        print('  %-18s %s' % (n, h))

    ids = pd.read_csv(os.path.join(INPUTS, 'indices_drop.csv')).values.flatten()
    ids = ids.astype(int)
    ## indices_drop.csv is 1-BASED (written by R); main.py l.278 does the same
    ## subtraction. Feeding the loaders 1-based indices would shift every row
    ## and column by one and still produce plausible numbers.
    idx = ids - 1
    print('\n[load] %d cellules actives (ids 1-based [%d,%d] -> 0-based)'
          % (len(ids), ids.min(), ids.max()))

    print('\n--- lambda_0 : chemin load_h5, SANS cap ---')
    lam_m = load_h5(os.path.join(INPUTS, 'lambda.h5'), idx, 'lambda',
                    cap_colmass=None)
    lam = np.asarray(lam_m.sum(axis=0)).ravel()

    print('\n--- mu_0 : chemin load_migration ---')
    mu_m = load_migration(os.path.join(INPUTS, 'mig_flows.h5'), idx)
    mu = np.asarray(mu_m.sum(axis=0)).ravel()

    rows = [
        {'matrix': 'lambda0', 'object': 'trade shares', 'loader': 'load_h5',
         'capped': 'no', 'max_colmass': round(float(lam.max()), 3),
         'n_cols_gt_1.5': int((lam > 1.5).sum()),
         'n_cols_gt_3': int((lam > 3).sum()), 'n_cells': int(lam.size),
         'source_file': 'lambda.h5', 'source_md5': stamps['lambda.h5']},
        {'matrix': 'mu0', 'object': 'bilateral migration',
         'loader': 'load_migration', 'capped': 'n/a',
         'max_colmass': round(float(mu.max()), 3),
         'n_cols_gt_1.5': int((mu > 1.5).sum()),
         'n_cols_gt_3': int((mu > 3).sum()), 'n_cells': int(mu.size),
         'source_file': 'mig_flows.h5', 'source_md5': stamps['mig_flows.h5']},
    ]
    df = pd.DataFrame(rows)
    print('\n' + df.to_string(index=False))

    ## The claim the paragraph makes, stated as the script measures it. It is
    ## PRINTED, never asserted: this script reports the calibration, it does
    ## not adjudicate it.
    print('\n[lecture] la degenerescence est-elle propre au commerce ?')
    print('  lambda_0 : max %.3f | %d colonnes > %g'
          % (lam.max(), (lam > C_BOUND).sum(), C_BOUND))
    print('  mu_0     : max %.3f | %d colonnes > %g'
          % (mu.max(), (mu > C_BOUND).sum(), C_BOUND))
    ## Console output stays ASCII on purpose: the Windows console of this
    ## project is cp1252, and a warning that raises UnicodeEncodeError is a
    ## warning nobody reads.
    if not (mu.max() < lam.max() and (mu > C_BOUND).sum() < (lam > C_BOUND).sum()):
        print("  [!] mu_0 n'est PAS moins concentre que lambda_0 sur les DEUX")
        print('      statistiques. La phrase du manuscrit -- "la degenerescence')
        print('      est propre au commerce" -- ne tient plus telle quelle sur')
        print('      celle qui echoue. A arbitrer avant publication.')

    if DRY:
        print('\n[dry-run] %s NON ecrit.' % CSV)
        return 0
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(CSV, index=False)
    print('\n[write] %s' % CSV)
    return 0


if __name__ == '__main__':
    sys.exit(main())
