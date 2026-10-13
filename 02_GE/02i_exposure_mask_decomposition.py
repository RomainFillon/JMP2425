# -*- coding: utf-8 -*-
r"""02i_exposure_mask_decomposition.py -- source decomposition of the Northern welfare margins.

THE THREE STATES COMPARED, all at f_a = 0, central draw, CAP_COLMASS=none:
    full                   warming everywhere: the headline runs in output/
    north_warming_only     warming at y >= 45 N only (mask 'north')
    nonnorth_warming_only  warming below 45 N only (mask 'south')
The no-warming run is common to the three (welfare_sc1_q50_NULLamen.csv).

WHAT IS SWITCHED OFF, AND WHAT IS NOT. A mask freezes exposure at 2015 outside
its kept cells (02g_exposure_masks.py), so warming is switched off there. General
equilibrium is not: trade, migration, wages and congestion operate in every run.
The columns therefore name the operation, not a channel.

THE QUANTITIES, population-weighted (population0.csv), latitude rounded to 0.1
degree, cumulative log welfare at TERMINAL_ROW = 15 -- the conventions of
emit_values.py:
    daily_total    cumsum log(W_daily / W_nowarming)
    annual_total   cumsum log(W_annual / W_nowarming)
    daily_extra    cumsum log(W_daily / W_annual)
    amenity_part   cumsum log(amenity_hat_daily / amenity_hat_annual), last row
    income_part    daily_extra - amenity_part   (the object of \GEchanIncome*)
over three regions: y >= 45, 45 <= y < 60, y >= 60. Regions below 45 N are not
produced: their cells are on the frozen side of the north mask.

THE MODEL IS NOT LINEAR. full != north_only + south_only in general; the
difference is the INTERACTION and it is written, never dropped.

BLOCKING CHECKS
    - four masked runs in intermediate_outputs/exposure_masks/runs_manifest.csv,
      md5 on disk equal to the manifest, one base fingerprint, equal to the
      inputs_md5 of the headline runs, same settings as the headline runs;
    - headline runs: md5 on disk equal to output/ge_runs_manifest.csv;
    - masked runs: md5 of both hat-paths equal to the manifest, and drf_md5
      recomputed from output/MC_damage_*.csv equal to the manifest (a
      re-estimation after the masked runs makes the decomposition stale);
    - daily_extra == daily_total - annual_total in every column (1e-12);
    The full column is checked against the body DOWNSTREAM, by emit_values.py,
    at 1e-9 on the unrounded values. This script used to compare it with
    output/values.tex; that file is written by 02d, which runs after this
    step, so the check read an artefact of the previous run. Removed on
    25/09/2026: a guard must read upstream of itself.

Writes: output/exposure_mask_decomposition.csv, output/exposure_mask_cells.csv
Usage:  python 02_GE/02i_exposure_mask_decomposition.py
"""
import os, sys, hashlib
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('JMP_ROOT') or os.path.dirname(HERE)
INP = os.path.join(ROOT, 'intermediate_outputs')
OUT = os.path.join(ROOT, 'output')
MASKS = os.path.join(INP, 'exposure_masks')
TERMINAL_ROW = 15
ATAG = '_NULLamen'
STATES = {'full': OUT,
          'north_warming_only': os.path.join(MASKS, 'north'),
          'nonnorth_warming_only': os.path.join(MASKS, 'south')}
REGIONS = [('north_45', 45, 90), ('band_45_60', 45, 60), ('band_60_90', 60, 90)]


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def read(folder, name):
    a = pd.read_csv(os.path.join(folder, name))
    return a.loc[:, ~a.columns.str.contains('^Unnamed')].values      # (T, N)


# ---- provenance ---------------------------------------------------------------
# dtype=str, keep_default_na=False : amenity_spec=null would otherwise be read as NaN,
# and nunique() ignores NaN.
head = pd.read_csv(os.path.join(OUT, 'ge_runs_manifest.csv'), dtype=str,
                   keep_default_na=False).set_index('file')
hfiles = ['welfare_sc%d_q50%s.csv' % (s, ATAG) for s in (1, 2, 3)]
for f in hfiles:
    if f not in head.index:
        sys.exit('[STOP] %s absent de output/ge_runs_manifest.csv' % f)
    if md5(os.path.join(OUT, f)) != head.loc[f, 'md5']:
        sys.exit('[STOP] %s : md5 disque != manifeste' % f)
hsub = head.loc[hfiles]
for col in ('inputs_md5', 'cap_colmass', 'amenity_spec', 'quantile'):
    if hsub[col].astype(str).nunique() != 1:
        sys.exit('[STOP] runs de tete heterogenes sur %s' % col)
H = hsub.iloc[0]

mman = pd.read_csv(os.path.join(MASKS, 'runs_manifest.csv'), dtype=str, keep_default_na=False)
need = {(m, 'welfare_sc%d_q50%s.csv' % (s, ATAG)) for m in ('north', 'south') for s in (2, 3)}
have = set(zip(mman['mask'], mman['file']))
if need - have:
    sys.exit('[STOP] runs masques absents du manifeste : %s' % sorted(need - have))
for _, r in mman.iterrows():
    if md5(os.path.join(MASKS, r['mask'], r['file'])) != r['md5']:
        sys.exit('[STOP] %s/%s : md5 disque != manifeste' % (r['mask'], r['file']))
    if str(r['inputs_md5_base']) != str(H['inputs_md5']):
        sys.exit('[STOP] %s/%s : empreinte de base %s != runs de tete %s'
                 % (r['mask'], r['file'], str(r['inputs_md5_base'])[:8], str(H['inputs_md5'])[:8]))
    for col in ('cap_colmass', 'amenity_spec', 'quantile'):
        if str(r[col]) != str(H[col]):
            sys.exit('[STOP] %s/%s : %s=%s, runs de tete %s' % (r['mask'], r['file'], col, r[col], H[col]))
    sc = r['file'].split('_')[1]                       # 'sc2' / 'sc3'
    for kind in ('amenity', 'productivity'):
        hp = os.path.join(MASKS, r['mask'], '%s_hat_path_%s_q50%s.csv' % (kind, sc, ATAG))
        if md5(hp) != r[kind + '_md5']:
            sys.exit('[STOP] %s : md5 disque != manifeste' % os.path.relpath(hp, ROOT))
_drf = sorted(f for f in os.listdir(OUT) if f.startswith('MC_damage_') and f.endswith('.csv'))
_drf_md5 = hashlib.md5(('\n'.join('%s:%s' % (f, md5(os.path.join(OUT, f))) for f in _drf)
                        + '\n').encode('ascii')).hexdigest()
if set(mman['drf_md5']) != {_drf_md5}:
    sys.exit('[STOP] les dose-reponses de output/ (drf_md5 %s) ne sont plus celles des '
             'runs masques (%s). Relancer 02g puis 02h.'
             % (_drf_md5[:8], ', '.join(sorted(v[:8] for v in set(mman['drf_md5'])))))
print('[provenance] 3 runs de tete + 4 runs masques ; empreinte de base %s ; cap=%s amenity=%s'
      % (str(H['inputs_md5'])[:8], H['cap_colmass'], H['amenity_spec']))

# ---- cells ---------------------------------------------------------------------
ids = pd.read_csv(os.path.join(INP, 'indices_drop.csv')).iloc[:, 0].astype(int).values
xy = pd.read_csv(os.path.join(INP, 'identifier_xy_final.csv'))
xy.columns = [c.strip().strip('"') for c in xy.columns]
y = xy.set_index(xy['id'].astype(int))['y'].reindex(ids).round(1).values
pop = pd.read_csv(os.path.join(INP, 'population0.csv')).set_index('id')['value'].reindex(ids).values
if np.isnan(y).any() or np.isnan(pop).any():
    sys.exit('[STOP] latitude ou population manquante pour une cellule active')

# ---- per-cell margins ------------------------------------------------------------
W1 = read(OUT, 'welfare_sc1_q50%s.csv' % ATAG)
cum = lambda a, b: np.cumsum(np.log(a / b), axis=0)[TERMINAL_ROW, :]
cells = {}
for st, d in STATES.items():
    W2, W3 = read(d, 'welfare_sc2_q50%s.csv' % ATAG), read(d, 'welfare_sc3_q50%s.csv' % ATAG)
    A2 = read(d, 'amenity_hat_path_sc2_q50%s.csv' % ATAG)
    A3 = read(d, 'amenity_hat_path_sc3_q50%s.csv' % ATAG)
    if not (W1.shape == W2.shape == W3.shape == A2.shape == A3.shape) or W1.shape[1] != len(ids):
        sys.exit('[STOP] %s : formes incoherentes %s' % (st, [W1.shape, W2.shape, W3.shape, A2.shape]))
    extra = cum(W2, W3)
    amen = np.cumsum(np.log(A2 / A3), axis=0)[-1, :]
    cells[st] = {'daily_total': cum(W2, W1), 'annual_total': cum(W3, W1),
                 'daily_extra': extra, 'amenity_part': amen, 'income_part': extra - amen}

QUANT = ['daily_extra', 'income_part', 'amenity_part', 'daily_total', 'annual_total']
rows = []
for reg, lo, hi in REGIONS:
    m = (y >= lo) & (y < hi)
    for q in QUANT:
        v = {st: float(np.average(cells[st][q][m], weights=pop[m])) for st in STATES}
        inter = v['full'] - v['north_warming_only'] - v['nonnorth_warming_only']
        rows.append({'region': reg, 'lat_lo': lo, 'lat_hi': hi, 'n_cells': int(m.sum()),
                     'pop_m': pop[m].sum() / 1e6, 'quantity': q, **v, 'interaction': inter,
                     'share_north_warming_only': 100 * v['north_warming_only'] / v['full'],
                     'share_nonnorth_warming_only': 100 * v['nonnorth_warming_only'] / v['full'],
                     'share_interaction': 100 * inter / v['full']})
dec = pd.DataFrame(rows)

# ---- identity: daily_extra = daily_total - annual_total, in every column ----------
for reg, _, _ in REGIONS:
    s = dec[dec['region'] == reg].set_index('quantity')
    for st in list(STATES) + ['interaction']:
        gap = s.loc['daily_extra', st] - (s.loc['daily_total', st] - s.loc['annual_total', st])
        if abs(gap) > 1e-12:
            sys.exit('[STOP] %s/%s : daily_extra != daily_total - annual_total (%.3e)' % (reg, st, gap))

dec.to_csv(os.path.join(OUT, 'exposure_mask_decomposition.csv'), index=False, float_format='%.10g')
cc = pd.read_csv(os.path.join(MASKS, 'masks_cells.csv'))
cc.to_csv(os.path.join(OUT, 'exposure_mask_cells.csv'), index=False)

pd.set_option('display.width', 200)
print(dec[['region', 'quantity', 'full', 'north_warming_only', 'nonnorth_warming_only',
           'interaction', 'share_north_warming_only', 'share_nonnorth_warming_only',
           'share_interaction']]
      .to_string(index=False, float_format=lambda v: '%+.6f' % v))
print('[ecrit] output/exposure_mask_decomposition.csv, output/exposure_mask_cells.csv')
