# -*- coding: utf-8 -*-
"""
decompose_metric_channels.py -- splits the northern daily-extra into an
EVALUATION channel (Jensen) and a SHAPE channel (the estimated curve differs).

WHAT IT COMPUTES
----------------
Three GE runs, all at q50, all from the same grid and the same settings:

    sc2 = daily curve    evaluated on the DAILY DISTRIBUTION
    sc3 = annual curve   evaluated at the ANNUAL MEAN
    sc4 = daily curve    evaluated at the ANNUAL MEAN

sc4 is the counterfactual that separates the two effects, by exact telescoping:

    (sc2 - sc3)        =      (sc2 - sc4)       +      (sc4 - sc3)
     daily-extra              EVALUATION               SHAPE
     (the published object)   same daily curve,        same evaluation at the
                              distribution vs mean     mean, different curves
                              -> Jensen                -> aggregating before
                                                          estimating

sc4 is produced by main.py with SCENARIO=3 and ANNUAL_DRF_SOURCE=dailycurve,
which sets _sclabel = 4 and writes welfare_sc4_<q>_dailycurve.csv.

Aggregation: mean over cells at latitude >= 45N, weighted by 2015 population, of
the cumulative dlog U at 2095, the last period.

CROSS-CHECK. The re-derived total must reproduce the northern daily-extra that
emit_values.py computes for the same aggregate. The two paths share no code: this
module builds its own latitude mask, its own population weights and its own
cumulation. Agreement is therefore a check on the implementation, not a tautology,
and it is what separates this decomposition from an ad hoc reconstruction.

WHY THE INPUT CHECK IS BLOCKING
-------------------------------
The three welfare matrices must come from ONE run grid. When they did not -- the
earlier arrangement kept them in three separate directories, one of which held a
welfare_sc3_q50_capNone.csv that was NOT sc4's partner -- the natural pairing,
taking the sc3 sitting beside sc4, gave:

    evaluation = +0.147   shape = -0.241   total = -0.093

a Jensen channel SIXTEEN TIMES too large and larger in absolute value than the
total. One would conclude that Jensen carries the result, which is precisely the
mechanism the paper rules out. A third pairing gave -0.661 / +0.504, with both
signs reversed.

Nothing in a directory tree signals which triplet is the right one, and a comment
protects only those who read it. The guard is therefore mechanical: this module
refuses to run unless the three files are listed in output/ge_runs_manifest.csv
against the SAME input fingerprint and the SAME settings, and unless each file
still hashes to what the manifest recorded.

Reads:  output/welfare_sc{2,3}_q50.csv, output/welfare_sc4_q50_dailycurve.csv,
        output/ge_runs_manifest.csv,
        intermediate_outputs/{indices_drop,identifier_xy_final,population0}.csv
Writes: nothing. compute() returns a dict and prints its report.
Needs:  the GE grid must have run, sc4 included.
"""
import csv
import hashlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))   # .../replication/02_GE
PKG = os.path.dirname(HERE)                         # .../replication

OUT = os.path.join(PKG, 'output')
## LOT 85o : sous le bac a sable de tools/known_bad_emission.py, ce module lit la
## COPIE d'output/ (EMIT_SANDBOX_DIR), comme emit_values.py qui l'importe. Sans
## cela son controle de provenance lisait l'arbre vivant intact et ne voyait pas
## la perturbation de la copie. Lecture seule : ce module n'ecrit rien.
if os.environ.get('EMIT_SANDBOX_DIR'):
    OUT = os.path.normpath(os.path.abspath(os.environ['EMIT_SANDBOX_DIR']))
TMP = os.path.join(PKG, 'intermediate_outputs')

##################################################################
# LE TRIPLET DEPEND DU BRAS D'AMENITE.
#
# Le corps du papier tourne a la borne basse f_a = 0, donc la decomposition
# qu'il cite doit en sortir aussi. Les noms etaient figes sur le bras
# interieur, et \GEnorthChanEval / \GEnorthChanShape ont donc porte une
# decomposition qui sommait a -0,153 a cote d'un daily extra de -0,109 :
# deux bras dans une meme phrase.
#
# Le defaut est INCHANGE (tag vide = bras interieur), et SHARED_KEYS contient
# deja `amenity_spec`, donc un triplet melange reste refuse par le manifeste :
# on parametre, on n'affaiblit pas.
##################################################################
def triplet_for(tag=''):
    """Les trois welfare du bras `tag` ('' = interieur, '_NULLamen' = f_a=0).

    ⚠️ L'ORDRE DES SUFFIXES DE sc4 N'EST PAS CELUI QU'ON SUPPOSE.
    main.py compose suffix = _sc<label><q> + <amenite> + <cap> + <drf>, donc
    le bras vient AVANT « _dailycurve » : welfare_sc4_q50_NULLamen_dailycurve.
    A `tag` vide les deux ordres coincident, ce qui masque l'erreur jusqu'au
    premier bras non trivial.
    """
    return (f'welfare_sc2_q50{tag}.csv',
            f'welfare_sc3_q50{tag}.csv',
            f'welfare_sc4_q50{tag}_dailycurve.csv')


SC2, SC3, SC4 = triplet_for('')
TRIPLET = (SC2, SC3, SC4)

MANIFEST = os.path.join(OUT, 'ge_runs_manifest.csv')

AUX_IDS = os.path.join(TMP, 'indices_drop.csv')
AUX_XY = os.path.join(TMP, 'identifier_xy_final.csv')
AUX_POP = os.path.join(TMP, 'population0.csv')

# What the three runs must share. SCENARIO and QUANTILE are expected to differ:
# that is what the triplet is made of.
#
# inputs_md5, NOT run_id. The requirement is that the three runs read the same
# inputs under the same settings, and inputs_md5 says exactly that. run_id would
# be a proxy that fails in both directions: a resumed grid gets a new run_id
# while consuming identical inputs, and one invocation could in principle span an
# input that changed underneath it. run_id stays in the manifest as information.
SHARED_KEYS = ('inputs_md5', 'cap_colmass', 'amenity_spec')

# Numerical non-regression. Left unset on purpose: the values of the 29 July run
# describe inputs that no longer exist -- ERA5 regenerated, cells renumbered,
# tas_ssp2_* realigned. Pin them by hand, from a run that has been examined, or
# leave None and rely on the cross-check below, which needs no frozen reference.
EXPECTED_VALUES = None
TOL = 1e-3

LAT_MIN = 45.0
CROSS_TOL = 5e-4


def _md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _read_manifest():
    if not os.path.exists(MANIFEST):
        print('[REFUS] output/ge_runs_manifest.csv missing.')
        print('  It is written by the GE loop of run_all.sh, one row per run.')
        print('  Without it, nothing proves that the three welfare files come from the same')
        print('  run with the same settings, and a wrong pairing reverses the')
        print('  conclusion.')
        sys.exit(2)
    with open(MANIFEST, newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh))


def check_inputs(triplet=None):
    """Refuse to continue unless the three welfare files come from one grid.

    This replaces a pin on three fixed md5 values. The danger it guards against
    is the same one: a triplet assembled from runs that do not belong together.
    """
    triplet = triplet or TRIPLET
    rows = {r['file']: r for r in _read_manifest()}
    bad = []
    for name in triplet:
        p = os.path.join(OUT, name)
        if not os.path.exists(p):
            bad.append('  MISSING on disk     %s' % name)
            continue
        if name not in rows:
            bad.append('  MISSING in manifest %s' % name)
            continue
        got = _md5(p)
        if got != rows[name]['md5']:
            bad.append('  MD5 %s\n      manifest %s\n      disk     %s'
                       % (name, rows[name]['md5'], got))
    if bad:
        print('[REFUS] the triplet does not match the manifest:')
        print('\n'.join(bad))
        sys.exit(2)

    for k in SHARED_KEYS:
        vals = {rows[n][k] for n in triplet}
        if len(vals) != 1:
            print('[REFUS] the three runs do not share %s: %s' % (k, sorted(vals)))
            print('  Pairing arms with different settings gives a')
            print('  decomposition whose two channels can change sign.')
            sys.exit(2)

    r = rows[triplet[0]]
    print('[ok] triplet checked: inputs_md5=%s cap_colmass=%s amenity_spec=%s'
          % (r['inputs_md5'][:12], r['cap_colmass'], r['amenity_spec']))
    print('     run_id : %s' % ', '.join(sorted({rows[n]['run_id'] for n in triplet})))
    for n in triplet:
        print('     %-34s %s' % (n, rows[n]['md5'][:12]))


def load_welfare(name):
    d = pd.read_csv(os.path.join(OUT, name))
    d = d.loc[:, ~d.columns.str.contains('^Unnamed')]
    return d.values


def compute(cross_check=None, tag=''):
    """Return {'evaluation', 'shape', 'total'} for the NORTH >= 45N aggregate.

    Entry point for consumers, emit_values.py in particular, which turns these
    into macros. check_inputs() runs here, so importing the function does not
    bypass the provenance check. Raises SystemExit(2) when the triplet does not
    match the manifest, when the shapes disagree, or when the cross-check fails.
    Prints its report: the traceability is part of the deliverable.

    cross_check: the northern daily-extra computed by the caller on the same
    grid. When given, the re-derived total must reproduce it.
    """
    _t = triplet_for(tag)
    check_inputs(_t)

    w2, w3, w4 = (load_welfare(_t[0]), load_welfare(_t[1]), load_welfare(_t[2]))
    if not (w2.shape == w3.shape == w4.shape):
        sys.exit('[REFUS] incompatible shapes: %s %s %s'
                 % (w2.shape, w3.shape, w4.shape))
    T, N = w2.shape
    print('[ok] welfare (%d periods, %d cells)' % (T, N))

    ids = pd.read_csv(AUX_IDS).iloc[:, 0].astype(int).values
    xy = pd.read_csv(AUX_XY)
    xy.columns = [c.strip('"') for c in xy.columns]
    id2y = dict(zip(xy['id'].astype(int), xy['y']))
    pop = pd.read_csv(AUX_POP)
    id2p = dict(zip(pop['id'].astype(int), pop['value']))

    lat = np.array([id2y[i] for i in ids], dtype=float)
    wpop = np.array([id2p[i] for i in ids], dtype=float)
    north = lat >= LAT_MIN

    def cum_last(a, b):
        return np.cumsum(np.log(a / b), axis=0)[-1, :]

    def agg(v, mask):
        return float(np.average(v[mask], weights=wpop[mask]))

    _allcells = np.ones(len(ids), bool)
    got = {
        'evaluation': agg(cum_last(w2, w4), north),   # same daily curve
        'shape':      agg(cum_last(w4, w3), north),   # same evaluation at the mean
        'total':      agg(cum_last(w2, w3), north),   # = the published daily-extra
        # Le monde entier, RENDU et plus seulement imprime. Il coute le meme
        # calcul, et un chiffre imprime sans etre rendu finit recopie a la main.
        'world_evaluation': agg(cum_last(w2, w4), _allcells),
        'world_shape':      agg(cum_last(w4, w3), _allcells),
        'world_total':      agg(cum_last(w2, w3), _allcells),
    }

    print('\nNORTH >= %.0fN: %d cells, %.0f M inhabitants (pop 2015)'
          % (LAT_MIN, north.sum(), wpop[north].sum() / 1e6))
    print('Cumulative dlog U over the whole horizon (last period), population-weighted mean:\n')
    print('  EVALUATION channel (Jensen)   sc2 - sc4 : %+.4f' % got['evaluation'])
    print('  SHAPE channel (aggregation)   sc4 - sc3 : %+.4f' % got['shape'])
    print('  ' + '-' * 46)
    print('  total  (daily-extra)          sc2 - sc3 : %+.4f' % got['total'])
    share = 100 * got['evaluation'] / got['total'] if got['total'] else float('nan')
    print('\n  share of the evaluation channel in the total: %+.1f %%' % share)
    allc = np.ones(len(ids), bool)
    print('  (whole-world report, not claimed: eval %+.4f / shape %+.4f)'
          % (agg(cum_last(w2, w4), allc), agg(cum_last(w4, w3), allc)))

    if cross_check is not None:
        d = abs(got['total'] - cross_check)
        verdict = 'OK' if d <= CROSS_TOL else 'RED'
        print("\n[croise] re-derived total %+.4f vs the caller's Northern daily-extra "
              '%+.4f -- gap %.5f -> %s' % (got['total'], cross_check, d, verdict))
        if d > CROSS_TOL:
            print('  The two paths aggregate the same object without sharing code.')
            print('  A gap signals a divergence of mask, weighting or')
            print('  cumulation, not a result.')
            raise SystemExit(2)

    if EXPECTED_VALUES is None:
        print('\n[note] non-regression NOT ARMED: EXPECTED_VALUES is None.')
    else:
        bad = [k for k, v in EXPECTED_VALUES.items() if abs(got[k] - v) > TOL]
        if bad:
            print('\n[ROUGE] drift on %s (tol %g):' % (', '.join(bad), TOL))
            for k in bad:
                print('    %-11s expected %+.4f got %+.4f'
                      % (k, EXPECTED_VALUES[k], got[k]))
            print('  This is NOT a scientific verdict -- diagnose the pipeline.')
            raise SystemExit(2)
        print('\n[ok] non-regression: 3/3 values reproduced (tol %g).' % TOL)
    return got


def main():
    compute()
    return 0


if __name__ == '__main__':
    sys.exit(main())
