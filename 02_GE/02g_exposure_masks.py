# -*- coding: utf-8 -*-
r"""02g_exposure_masks.py -- masked climate exposure for the source decomposition.

WHAT IT BUILDS
    Two input sets for the GE solver, identical to intermediate_outputs/ except
    for the two climate files:

      north   cells at y >= 45 keep their projected exposure; every other cell
              is FROZEN at its 2015 exposure.
      south   the complement: the active cells below 45 N keep their projected
              exposure; every cell at y >= 45 is frozen.

    The two kept sets partition the active cells (indices_drop.csv) exactly.

HOW FREEZING SWITCHES WARMING OFF
    A frozen cell has the 2015 value copied into every later period, in the
    daily file (one column per period) and in the annual file (one row per
    cell and period). The solver computes damage from the CHANGE in exposure
    between periods (02c_model/forward.py), so a frozen cell receives exactly
    zero damage in every period, as in the no-warming baseline. No line of the
    solver changes, and no general-equilibrium channel is switched off: goods
    trade, migration, wages and congestion operate everywhere in both runs. What
    is switched off is warming, outside the North in the first set and inside it
    in the second.

WHERE IT WRITES
    intermediate_outputs/exposure_masks/<mask>/
      - the seven non-climate inputs main.py reads, copied and checked by md5;
      - the two masked climate files;
      - the dose-response tables main.py reads from its output folder
        (MC_damage_*.csv from output/), copied and checked by md5, because the
        solver reads its inputs and its dose-responses from JMP_DATA and JMP_OUT
        and the masked runs use this folder as both.

    It reads intermediate_outputs/ and output/ and writes nothing there.

BLOCKING CHECKS, READ BACK FROM DISK
    - every frozen cell has all its periods equal to 2015, daily and annual;
    - every kept cell is unchanged, daily and annual;
    - the two kept sets are disjoint and their union is the active set;
    - the rounded and unrounded latitude give the same partition;
    - every copied file has the md5 of its source.

Usage:  python 02_GE/02g_exposure_masks.py
"""
import os, shutil, hashlib, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('JMP_ROOT') or os.path.dirname(HERE)
INP = os.path.join(ROOT, 'intermediate_outputs')
OUT = os.path.join(ROOT, 'output')
DEST = os.path.join(INP, 'exposure_masks')
LAT_NORTH = 45.0

# What main.py reads from JMP_DATA (02c_model/main.py and data_cleaning.py).
SOLVER_INPUTS = ['lambda.h5', 'mig_flows.h5', 'indices_drop.csv', 'population0.csv',
                 'wage0.csv', 'data_final_population.csv', 'data_final_productivity.csv']
CLIMATE = ['tas_ssp2_daily.csv', 'tas_ssp2_annual.csv']


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def copy_checked(src, dst):
    shutil.copyfile(src, dst)
    if md5(src) != md5(dst):
        sys.exit('[STOP] copie alteree : %s' % dst)


# ---- the partition ----------------------------------------------------------
xy = pd.read_csv(os.path.join(INP, 'identifier_xy_final.csv'))
xy.columns = [c.strip().strip('"') for c in xy.columns]
lat = dict(zip(xy['id'].astype(int), xy['y'].astype(float)))
active = set(pd.read_csv(os.path.join(INP, 'indices_drop.csv')).iloc[:, 0].astype(int))
north = {i for i in active if round(lat[i], 1) >= LAT_NORTH}
if north != {i for i in active if lat[i] >= LAT_NORTH}:
    sys.exit('[STOP] la latitude arrondie et la latitude brute ne donnent pas la meme partition.')
south = active - north
assert north and south and not (north & south) and (north | south) == active
print('[masks] active cells %d : north (y >= %g) %d | south %d'
      % (len(active), LAT_NORTH, len(north), len(south)))

daily = pd.read_csv(os.path.join(INP, 'tas_ssp2_daily.csv'))
annual = pd.read_csv(os.path.join(INP, 'tas_ssp2_annual.csv'))
periods = sorted([c for c in daily.columns if str(c).isdigit()], key=int)
if periods[0] != '2015':
    sys.exit('[STOP] premiere periode du fichier journalier : %s, attendu 2015' % periods[0])
if not {'id', 'year', 'temp_mean'} <= set(annual.columns):
    sys.exit('[STOP] colonnes annuelles inattendues : %s' % list(annual.columns))
a2015 = annual[annual['year'] == 2015].set_index('id')['temp_mean']


def freeze(keep_ids):
    d = daily.copy()
    frozen = ~d['id'].isin(keep_ids)
    for c in periods[1:]:
        d.loc[frozen, c] = d.loc[frozen, '2015']
    a = annual.copy()
    fa = ~a['id'].isin(keep_ids)
    a.loc[fa, 'temp_mean'] = a.loc[fa, 'id'].map(a2015).values
    return d, a


drf = sorted(f for f in os.listdir(OUT) if f.startswith('MC_damage_') and f.endswith('.csv'))
if not drf:
    sys.exit('[STOP] aucune table MC_damage_*.csv dans output/')

for name, keep in (('north', north), ('south', south)):
    dd = os.path.join(DEST, name)
    os.makedirs(dd, exist_ok=True)
    for f in SOLVER_INPUTS:
        copy_checked(os.path.join(INP, f), os.path.join(dd, f))
    for f in drf:
        copy_checked(os.path.join(OUT, f), os.path.join(dd, f))
    d, a = freeze(keep)
    d.to_csv(os.path.join(dd, 'tas_ssp2_daily.csv'), index=False)
    a.to_csv(os.path.join(dd, 'tas_ssp2_annual.csv'), index=False)

    # ---- blocking checks, read back from disk ------------------------------
    d2 = pd.read_csv(os.path.join(dd, 'tas_ssp2_daily.csv'))
    a2 = pd.read_csv(os.path.join(dd, 'tas_ssp2_annual.csv'))
    if d2.shape != daily.shape or a2.shape != annual.shape:
        sys.exit('[STOP] %s : forme des fichiers climat alteree' % name)
    k = d2['id'].isin(keep).values
    fz = d2.loc[~k, periods].fillna(-9e9).values
    if not (fz == fz[:, [0]]).all():
        sys.exit('[STOP] %s : une cellule gelee a une periode differente de 2015' % name)
    if not np.array_equal(d2.loc[k, periods].fillna(-9e9).values,
                          daily.loc[k, periods].fillna(-9e9).values):
        sys.exit('[STOP] %s : une cellule gardee a ete alteree (journalier)' % name)
    ka = a2['id'].isin(keep).values
    piv = a2[~ka].pivot(index='id', columns='year', values='temp_mean')
    if not np.allclose(piv.values, piv[[2015]].values, equal_nan=True):
        sys.exit('[STOP] %s : annuel, une cellule gelee differe de 2015' % name)
    if not np.array_equal(a2.loc[ka, 'temp_mean'].values, annual.loc[ka, 'temp_mean'].values):
        sys.exit('[STOP] %s : une cellule gardee a ete alteree (annuel)' % name)
    n_keep = len(set(d2.loc[k, 'id']) & active)
    if n_keep != len(keep):
        sys.exit('[STOP] %s : %d cellules actives gardees, attendu %d' % (name, n_keep, len(keep)))
    print('[%s] %d active cells keep their trajectory ; %d copied files checked ; '
          'lambda.h5 %s' % (name, n_keep, len(SOLVER_INPUTS) + len(drf),
                            md5(os.path.join(dd, 'lambda.h5'))[:12]))

pd.DataFrame({'mask': ['north', 'south'],
              'n_cells_kept': [len(north), len(south)],
              'n_cells_active': [len(active)] * 2,
              'lat_threshold': [LAT_NORTH] * 2}).to_csv(
    os.path.join(DEST, 'masks_cells.csv'), index=False)
print('[masks] ecrit : %s' % os.path.relpath(DEST, ROOT))
