"""00d_poids_population.py -- les deux poids de population, depuis la matrice (lot 46).

Remplace la repartition PROPORTIONNELLE qui n'existait nulle part dans la chaine : les
anciens poids (lengths_country.csv, population_weight_all.csv) comptaient la population
ENTIERE d'une cellule pour chaque pays (ou polygone) qui la touche (RECAP 3.230-3.232).

Entree  : matrice_population_cellule_pays.csv (00c_matrice_population.py) :
          x, y, GID_0, year, pop -- population de l'INTERSECTION cellule x pays, GHSL 30''
          dans gadm_410.gpkg non simplifie, sept periodes GHSL.
Sorties (dans <repertoire>) :
  poids_par_cellule.csv : x, y, GID_0, pop, part_de_la_cellule   (periode 2015 SEULE)
      part_de_la_cellule = pop / population de la cellule rattachee a un pays
      -> somme 1 PAR CELLULE. Repartit un taux national entre les cellules (02a) et
         agrege un resultat de cellule vers les pays en proportion (aval).
  poids_par_pays.csv    : x, y, GID_0, year, year_end, pop, part_du_pays (7 periodes)
      part_du_pays = pop / population du pays -> somme 1 PAR PAYS et periode.
      Pondere le climat d'un pays ; year_end = fenetre de validite, comme
      lengths_country.csv (00_population_weights.R l.195-207).
Periode 2015 : l'annee de calibration du modele (arbitrage du lot 46, 3.3).
Temoins, AVANT toute ecriture : les deux normalisations (tools/temoins_poids.py), puis
la population contre le monde et la Banque mondiale (tools/temoin_population_pays.R).
Usage : python 01_estimation/00d_poids_population.py <matrice.csv> <repertoire>
"""
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'tools'))
from temoins_poids import temoin_part_de_la_cellule, temoin_part_du_pays  # noqa: E402

## Sans argument (run_all.sh), les chemins canoniques de intermediate_outputs/.
_I = os.path.join(ROOT, 'intermediate_outputs')
MAT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_I, 'matrice_population_cellule_pays.csv')
OUT = sys.argv[2] if len(sys.argv) > 2 else _I
PERIODE_CALIBRATION = 2015
YEAR_LAST = 2019          # derniere annee ERA5 couverte, comme YEAR_LAST de 00_population_weights.R

m = pd.read_csv(MAT, keep_default_na=False)
m = m[(m.GID_0 != '') & (m['pop'] > 0)].copy()
m['x'] = m['x'].round(1); m['y'] = m['y'].round(1)

# ---- par cellule, periode 2015 ---------------------------------------------------
c = m[m.year == PERIODE_CALIBRATION][['x', 'y', 'GID_0', 'pop']].copy()
c['part_de_la_cellule'] = c['pop'] / c.groupby(['x', 'y'])['pop'].transform('sum')
n_c = temoin_part_de_la_cellule(c)
print('[temoin PAR CELLULE] %d cellules, chacune somme a 1' % n_c)

# ---- par pays, sept periodes ------------------------------------------------------
p = m[['x', 'y', 'GID_0', 'year', 'pop']].copy()
p['part_du_pays'] = p['pop'] / p.groupby(['GID_0', 'year'])['pop'].transform('sum')
n_p = temoin_part_du_pays(p)
print('[temoin PAR PAYS] %d pays-periodes, chacun somme a 1' % n_p)
ans = sorted(p.year.unique())
fin = {a: (ans[i + 1] - 1 if i + 1 < len(ans) else YEAR_LAST) for i, a in enumerate(ans)}
p['year_end'] = p['year'].map(fin)

# ---- population contre le monde (R, meme fonction que le test connu-faux) ---------
rs = os.environ.get('RSCRIPT', r'C:\Program Files\R\R-4.3.1\bin\Rscript.exe')
r = subprocess.run([rs, '--vanilla', os.path.join(ROOT, 'tools', 'valide_matrice_population.R'), MAT],
                   cwd=ROOT, capture_output=True, text=True)
print(r.stdout.strip()); print(r.stderr.strip()[-1500:])
if r.returncode != 0:
    sys.exit('[STOP] temoin population : la matrice est refusee, rien n est ecrit.')

os.makedirs(OUT, exist_ok=True)
c.sort_values(['x', 'y', 'GID_0']).to_csv(os.path.join(OUT, 'poids_par_cellule.csv'), index=False)
p[['x', 'y', 'GID_0', 'year', 'year_end', 'pop', 'part_du_pays']].sort_values(
    ['x', 'y', 'year', 'GID_0']).to_csv(os.path.join(OUT, 'poids_par_pays.csv'), index=False)
print('[write] poids_par_cellule.csv (%d lignes) et poids_par_pays.csv (%d lignes) dans %s'
      % (len(c), len(p), OUT))
