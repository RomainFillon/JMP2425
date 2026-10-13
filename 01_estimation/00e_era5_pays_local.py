"""00e_era5_pays_local.py -- le climat par pays, reconstruit LOCALEMENT depuis la grille (lot 46-51).

POURQUOI LOCAL. 00_era5_distribution.R (serveur) ecrit deux sorties d'un meme passage :
la grille ERA_tas_distribution_grid_annual.csv (jours par cellule, annee et intervalle de
1 degC) et le fichier par pays ERA_tas_annual.csv. Le second n'est qu'une ponderation de
la premiere : pour un pays et une annee,
    count(valeur) = somme sur les cellules de  poids(cellule, pays, periode) x jours,
la periode GHSL etant la derniere <= annee ; share = count / somme par (pays, annee)
(00_era5_distribution.R, recette identique). Changer les poids ne demande donc pas le
serveur : la frontiere serveur s'arrete a la GRILLE, qui reste un artefact livre.

Poids : poids_par_pays.csv, colonne part_du_pays (00d_poids_population.py) -- population
de l'intersection cellule x pays sur la population du pays, sept periodes GHSL.

CERTIFICATION DU CHEMIN LOCAL : tools/test_era5_chemin_local.py applique agrege() aux
ANCIENS poids et exige de reproduire au bit l'ancien fichier serveur (mesure du
27/09/2026 : 225 346 lignes sur 225 346, ecart 5,6e-15). Elle lit ces deux fichiers par
CHEMIN EXPLICITE (archive) ; la chaine, elle, ne les lit plus (lot 51).
Usage : python 01_estimation/00e_era5_pays_local.py <poids_par_pays.csv> <sortie.csv>
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
I = os.path.join(os.path.dirname(HERE), 'intermediate_outputs')
PERIODES = np.array([1985, 1990, 1995, 2000, 2005, 2010, 2015])


def charge_grille(chemin=None):
    g = pd.read_csv(chemin or os.path.join(I, 'ERA_tas_distribution_grid_annual.csv'))
    g['per'] = PERIODES[np.searchsorted(PERIODES, g['year'].values, side='right') - 1]
    g['x'] = g['x'].round(1); g['y'] = g['y'].round(1)
    return g


def agrege(g, w, col):
    """Climat par pays : grille g ponderee par la colonne `col` de w (x, y, GID_0, year)."""
    w = w[['x', 'y', 'GID_0', 'year', col]].copy()
    w['x'] = w['x'].round(1); w['y'] = w['y'].round(1)
    j = g.merge(w, left_on=['x', 'y', 'per'], right_on=['x', 'y', 'year'], suffixes=('', '_w'))
    j['count_p'] = j[col] * j['count']
    c = j.groupby(['GID_0', 'year', 'value'], as_index=False)['count_p'].sum()
    c = c.rename(columns={'count_p': 'count'})
    c['share'] = c['count'] / c.groupby(['GID_0', 'year'])['count'].transform('sum')
    return c[['GID_0', 'year', 'value', 'count', 'share']].sort_values(['GID_0', 'year', 'value'])


if __name__ == '__main__':
    ## Sans argument (run_all.sh), les chemins canoniques de intermediate_outputs/.
    POIDS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(I, 'poids_par_pays.csv')
    SORTIE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(I, 'ERA_tas_annual.csv')
    poids = pd.read_csv(POIDS, keep_default_na=False)
    poids = poids[poids['GID_0'] != '']
    nouveau = agrege(charge_grille(), poids, 'part_du_pays')
    os.makedirs(os.path.dirname(os.path.abspath(SORTIE)), exist_ok=True)
    tmp = SORTIE + '.tmp'
    nouveau.to_csv(tmp, index=False)
    os.replace(tmp, SORTIE)                    # ecriture atomique
    print('[write] %s : %d lignes, %d pays, annees %d-%d' % (SORTIE, len(nouveau),
          nouveau.GID_0.nunique(), nouveau.year.min(), nouveau.year.max()))
