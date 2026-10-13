#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""ge_calibration_diagnostic.py -- statistiques de l'etat initial calibre.

POURQUOI CE FICHIER EXISTE
--------------------------
lambda.h5 fait 867 Mo. Le lire dans l'emetteur ajouterait une minute a chaque
emission, pour des statistiques qui ne dependent PAS du run : elles decrivent
les matrices d'entree, pas le bien-etre. Ce script les calcule UNE FOIS et
depose un CSV compagnon ; emit_values.py lit le CSV, jamais la matrice.

C'est la forme de cap_colmass_diagnostic.py, et comme lui il IMPORTE les
chargeurs du solveur -- load_h5, load_migration -- au lieu de relire le HDF5 a
la main. Deux chemins de chargement, deux resultats : c'est le defaut que ce
paquet a deja paye.

CE QU'IL PRODUIT
----------------
  lam_diag_median / mean / share_gt_099
      La diagonale de lambda_0 : la part de la demande de r satisfaite par r.
  clear_median / clear_q05 / clear_q95
      ventes / revenu par cellule, ventes_r = sum_s lambda_0(s,r) * w_s * L_s.
      CALCUL EXACT, par somme de colonne sur la matrice creuse -- et NON
      l'approximation lambda[r,s] ~ lambda[s,r] que 02a l.433 emploie sur un
      echantillon en la declarant « coarse, but fast ».
  mig_senders / mig_receivers
      Les plus gros emetteurs et recepteurs NETS de mu_0, en codes pays.
      Aucun volume : le manuscrit n'en publie pas, et la macro non plus.

CE QU'IL NE PRODUIT PAS, ET POURQUOI
------------------------------------
L'inversion Eaton-Kortum qui tient (w, L, tau) fixes et retrouve Z. Le
manuscrit l'annonce comme un controle et en cite les 5e-95e percentiles. Elle
n'existe NULLE PART dans le paquet : 02a ne calcule qu'un ratio approche sur un
echantillon. L'ecrire serait un calcul de modele neuf, pas une emission.
Les deux bornes restent donc non emises, et le CSV le DIT au lieu de laisser un
trou muet.

Reads:  intermediate_outputs/{lambda.h5, mig_flows.h5, indices_drop.csv,
        wage0.csv, population0.csv, identifier_xy_final.csv, poids_par_cellule.csv,
        population_all.csv}, data_inputs/abel_cohen_migration/bilat_mig.csv (domaine, lot 68)
Writes: output/ge_calib_diagnostic.csv
"""
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
INPUTS = os.path.join(PKG, 'intermediate_outputs')
OUT = os.path.join(PKG, 'output')
MODEL = os.path.join(HERE, '02c_model')

N_TOP = 5


def _load_data_cleaning():
    """Importe les chargeurs du solveur. On ne relit pas le HDF5 a la main."""
    spec = importlib.util.spec_from_file_location(
        'data_cleaning', os.path.join(MODEL, 'data_cleaning.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['data_cleaning'] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    need = ['population_all.csv',
            'lambda.h5', 'mig_flows.h5', 'indices_drop.csv', 'wage0.csv',
            'population0.csv', 'identifier_xy_final.csv', 'poids_par_cellule.csv']
    missing = [f for f in need if not os.path.exists(os.path.join(INPUTS, f))]
    if missing:
        print('[STOP] intrants absents : %s' % ', '.join(missing))
        print('       Producteur : 02a_calibration_at_1deg.r, 02b_mig_flows.R.')
        return 2

    dc = _load_data_cleaning()

    # indices_drop.csv est 1-BASED (ecrit par R) ; main.py l.278 fait la meme
    # soustraction. L'oublier decale toute la matrice, et sur une grille qui
    # tolererait l'indice cela passerait SANS erreur.
    ids = pd.read_csv(os.path.join(INPUTS, 'indices_drop.csv')).iloc[:, 0].astype(int).values
    idx = ids - 1
    n = len(idx)
    print('[grille] %d cellules actives' % n)

    t0 = time.time()
    # cap_colmass=None : output/ge_runs_manifest.csv porte cap_colmass=none sur
    # la totalite des runs. Lire la matrice cappee decrirait un bras qui n'a pas
    # tourne.
    lam = dc.load_h5(os.path.join(INPUTS, 'lambda.h5'), idx, 'lambda', None)
    print('[lambda] lue en %.0f s, forme %s' % (time.time() - t0, lam.shape))

    rows = []

    # ---- diagonale -----------------------------------------------------------
    diag = np.asarray(lam.diagonal()).ravel()
    rows += [('lam_diag_median', '%.4f' % np.median(diag), 'lambda.h5 via load_h5'),
             ('lam_diag_mean', '%.4f' % diag.mean(), 'lambda.h5 via load_h5'),
             ('lam_diag_share_gt_099', '%.1f' % (100.0 * (diag > 0.99).mean()),
              'lambda.h5 via load_h5, part des cellules en %'),
             # Mean import share, unweighted, in per cent: computed here from the
             # unrounded diagonal, not as the complement of the rounded mean.
             ('import_share_mean_pct', '%.6f' % (100.0 * (1.0 - diag.mean())),
              'lambda.h5 via load_h5, 100 x (1 - moyenne de lambda_rr), en %')]
    print('[diag] mediane %.4f | moyenne %.4f | >0.99 %.1f %%'
          % (np.median(diag), diag.mean(), 100.0 * (diag > 0.99).mean()))

    # ---- clearing : ventes / revenu -----------------------------------------
    # lambda_0 est ligne-stochastique sur les SOURCES : lambda[r, s] est la part
    # des achats de r venant de s. Les ventes de s valent donc
    #     ventes_s = sum_r lambda[r, s] * revenu_r,
    # une somme de COLONNE ponderee -- pas la somme de ligne, et pas non plus la
    # symetrie lambda[r,s] ~ lambda[s,r] que 02a emploie.
    wage = pd.read_csv(os.path.join(INPUTS, 'wage0.csv')).set_index('id')['value']
    pop = pd.read_csv(os.path.join(INPUTS, 'population0.csv')).set_index('id')['value']
    inc = (wage.reindex(ids).values * pop.reindex(ids).values)
    if not np.all(np.isfinite(inc)) or (inc <= 0).any():
        print('[STOP] revenu non fini ou non positif sur %d cellules.'
              % int((~np.isfinite(inc)).sum() + (inc <= 0).sum()))
        return 2
    sales = np.asarray(lam.multiply(inc[:, None]).sum(axis=0)).ravel()
    ratio = sales / inc
    rows += [('clear_median', '%.3f' % np.median(ratio), 'lambda.h5 x wage0 x population0'),
             ('clear_q05', 'NON CALCULE', 'inversion Eaton-Kortum absente du paquet'),
             ('clear_q95', 'NON CALCULE', 'inversion Eaton-Kortum absente du paquet')]
    print('[clear] ventes/revenu : mediane %.3f | q05 %.3f | q95 %.3f  '
          '(BRUT, sans inversion)'
          % (np.median(ratio), np.quantile(ratio, .05), np.quantile(ratio, .95)))
    print('        les deux bornes publiees supposent une inversion EK qui '
          'n_existe pas dans le paquet -> NON CALCULE')

    # ---- migration nette par pays -------------------------------------------
    mu = dc.load_migration(os.path.join(INPUTS, 'mig_flows.h5'), idx)
    xy = pd.read_csv(os.path.join(INPUTS, 'identifier_xy_final.csv'))
    xy.columns = [c.strip().strip('"') for c in xy.columns]
    ## Lot 51 : ce bloc prenait l'argmax d'une part normalisee PAR PAYS
    ## (lengths_country.csv, 1995 ; RECAP 3.229). La migration nette d'un pays est
    ## desormais la somme des migrations nettes de ses cellules, chacune ponderee
    ## par sa part de population (poids_par_cellule.csv, periode 2015).
    sys.path.insert(0, os.path.join(PKG, 'tools'))
    from temoins_poids import temoin_part_de_la_cellule
    parts = pd.read_csv(os.path.join(INPUTS, 'poids_par_cellule.csv'), keep_default_na=False)
    parts = parts[parts['GID_0'] != ''].copy()
    temoin_part_de_la_cellule(parts)
    parts['x'] = parts['x'].round(1)
    parts['y'] = parts['y'].round(1)
    cell = pd.DataFrame({'id': ids})
    m = xy.set_index('id')
    cell['x'] = m.reindex(ids)['x'].round(1).values
    cell['y'] = m.reindex(ids)['y'].round(1).values

    # mu_0 est ligne-stochastique : mu[o, d] est la part des partants de o vers d.
    # Les flux se pesent donc par la population d'origine.
    p0 = pop.reindex(ids).values
    outflow = np.asarray(mu.multiply(p0[:, None]).sum(axis=1)).ravel() - p0 * np.asarray(mu.diagonal()).ravel()
    inflow = np.asarray(mu.multiply(p0[:, None]).sum(axis=0)).ravel() - p0 * np.asarray(mu.diagonal()).ravel()
    cell['net'] = inflow - outflow
    cp = cell.merge(parts[['x', 'y', 'GID_0', 'part_de_la_cellule']], on=['x', 'y'], how='inner')
    cp['net'] = cp['net'] * cp['part_de_la_cellule']         # repartition PROPORTIONNELLE
    by = cp.groupby('GID_0')['net'].sum().sort_values()
    senders = list(by.head(N_TOP).index)
    receivers = list(by.tail(N_TOP).index[::-1])
    rows += [('mig_senders', ', '.join(senders), 'mig_flows.h5 via load_migration'),
             ('mig_receivers', ', '.join(receivers), 'mig_flows.h5 via load_migration')]
    print('[mig] emetteurs nets : %s' % ', '.join(senders))
    print('[mig] recepteurs nets: %s' % ', '.join(receivers))

    # ---- domaine du modele (lot 68, 30/09/2026) --------------------------------
    # Ce que la grille GE ne couvre pas. Trois ensembles de cellules HABITEES
    # (poids_par_cellule x population_all, la convention de 02b) : actives
    # (indices_drop), dans la grille mais inactives, hors de la grille
    # (identifier_xy_final n'a pas la cellule : Bagdad, Taipei, Doha...). Un pays
    # est HORS DOMAINE s'il n'a aucune population active, PARTIEL s'il en a hors
    # de la grille. Les flux internationaux perdus sont recalcules depuis la
    # source avec l'allocation de 02b STEP 8 (P_out, P_in proportionnels a
    # pop_zc) : ce sont ceux que 02b STEP 13 retire. Rien ici ne descend de mu_0.
    pa = pd.read_csv(os.path.join(INPUTS, 'population_all.csv'))
    pa.columns = [c.strip().strip('"') for c in pa.columns]
    pw = parts[['x', 'y', 'GID_0', 'part_de_la_cellule']].copy()
    pa['x'] = pa['x'].round(1)
    pa['y'] = pa['y'].round(1)
    xyr = xy[['id', 'x', 'y']].copy()
    xyr['x'] = xyr['x'].round(1)
    xyr['y'] = xyr['y'].round(1)
    dc_ = (pw.merge(pa.rename(columns={'value': 'pop_cell'}), on=['x', 'y'], how='left')
             .merge(xyr, on=['x', 'y'], how='left'))
    dc_['pop_zc'] = dc_['part_de_la_cellule'] * dc_['pop_cell'].fillna(0)
    dc_['etat'] = np.where(dc_['id'].isna(), 'hors_grille',
                           np.where(dc_['id'].isin(set(ids)), 'actif', 'grille_inactif'))
    world = float(dc_['pop_zc'].sum())
    hg = dc_[dc_['etat'] == 'hors_grille']
    byc = dc_.pivot_table(index='GID_0', columns='etat', values='pop_zc', aggfunc='sum', fill_value=0)
    for col in ('actif', 'hors_grille'):
        if col not in byc.columns:
            byc[col] = 0.0
    ctot = byc.sum(axis=1)
    c_out = ctot[byc['actif'] <= 0]
    c_part = byc[(byc['actif'] > 0) & (byc['hors_grille'] > 0)]
    bm = pd.read_csv(os.path.join(PKG, 'data_inputs', 'abel_cohen_migration', 'bilat_mig.csv'))
    bm = bm[(bm['year0'] == 2010) & (bm['orig'] != bm['dest'])]
    bm = bm[bm['orig'].isin(ctot.index) & bm['dest'].isin(ctot.index)]
    a_c = (byc['hors_grille'] / ctot.replace(0, np.nan)).fillna(0.0)
    ao = a_c.reindex(bm['orig']).values
    ad = a_c.reindex(bm['dest']).values
    lost = float((bm['da_pb_closed'].values * (1.0 - (1.0 - ao) * (1.0 - ad))).sum())
    rows += [
        ('dom_cells_out', '%d' % hg[['x', 'y']].drop_duplicates().shape[0],
         'poids_par_cellule : cellules habitees absentes de identifier_xy_final'),
        ('dom_pop_out_M', '%.6f' % (hg['pop_zc'].sum() / 1e6),
         'population de ces cellules (part x population_all), millions'),
        ('dom_pop_out_share_pct', '%.6f' % (100.0 * hg['pop_zc'].sum() / world),
         'idem, en %% de la population de poids_par_cellule x population_all (%.1f M)' % (world / 1e6)),
        ('dom_cty_out_n', '%d' % len(c_out), 'pays sans aucune population active'),
        ('dom_cty_out_pop_M', '%.6f' % (c_out.sum() / 1e6), 'population de ces pays, millions'),
        ('dom_cty_out_list', ' '.join(c_out.sort_values(ascending=False).index), 'ISO3, population decroissante'),
        ('dom_cty_partial_n', '%d' % len(c_part), 'pays avec population active ET population hors grille'),
        ('dom_mig_lost_M', '%.6f' % (lost / 1e6),
         'bilat_mig.csv (year0 2010, da_pb_closed) : flux touchant une cellule hors grille, '
         'allocation de 02b STEP 8 ; = flux retires par 02b STEP 13')]
    print('[domaine] hors grille : %d cellules, %.1f M hab (%.2f %%) | %d pays hors domaine '
          '(%.1f M) | %d partiels | flux internationaux perdus %.2f M'
          % (hg[['x', 'y']].drop_duplicates().shape[0], hg['pop_zc'].sum() / 1e6,
             100.0 * hg['pop_zc'].sum() / world, len(c_out), c_out.sum() / 1e6,
             len(c_part), lost / 1e6))

    df = pd.DataFrame(rows, columns=['quantity', 'value', 'source'])
    dest = os.path.join(OUT, 'ge_calib_diagnostic.csv')
    df.to_csv(dest, index=False)
    print('\n[write] output/ge_calib_diagnostic.csv (%d lignes)' % len(df))
    return 0


if __name__ == '__main__':
    sys.exit(main())
