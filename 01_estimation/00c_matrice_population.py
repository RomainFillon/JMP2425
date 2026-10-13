"""00c_matrice_population.py -- population de chaque intersection (cellule de 1 degre x pays).

Pays    : gadm_410.gpkg NON simplifie (356 508 polygones, 263 GID_0), rasterise sur la
          grille GHSL 30'' par la regle du CENTRE de pixel (~1 km) -> PARTITION.
REGLE DE RATTACHEMENT (lot 47) : un pixel HABITE dont le centre tombe hors de tout
          polygone est rattache au polygone le plus proche si et seulement si la distance
          de son centre a ce polygone est inferieure au PAS de la trame ; au-dela, c'est
          une anomalie et tools/temoin_rattachement.py refuse avant toute ecriture.
Cellule : decoupage de la chaine (00_population_weights.R l.130-131) : blocs de 120 x 120
          pixels depuis l'origine GHSL, rattaches au centre de degre le plus proche :
          x = -179.5 + col // 120, y = 88.5 - row // 120.
Sortie  : (x, y, GID_0, year, pop) ; GID_0 = '' pour les pixels encore hors pays.
Usage   : python 01_estimation/00c_matrice_population.py <sortie.csv> [<liste_rattaches.csv>]
Versions precedentes : .bak_20260927_lot47 (sans rattachement), md5 e25dd0b7...
"""
import os
import sys
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
from rasterio.features import rasterize
from rasterio.windows import Window

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(R, 'tools'))
from temoin_rattachement import temoin_rattachement  # noqa: E402

## Sans argument (run_all.sh), les chemins canoniques de intermediate_outputs/.
_I = os.path.join(R, 'intermediate_outputs')
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_I, 'matrice_population_cellule_pays.csv')
LISTE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_I, 'pixels_rattaches.csv')
TIF = os.path.join(R, 'data_inputs', 'GHSL_POP', 'tif', 'GHS_POP_E%d_GLOBE_R2023A_4326_30ss_V1_0.tif')
GPKG = os.path.join(R, 'data_inputs', 'gadm36_shp', 'gadm_410.gpkg')
ANS = (1985, 1990, 1995, 2000, 2005, 2010, 2015)
t0 = time.time()

pays = gpd.read_file(GPKG, columns=['GID_0'])
pays = pays[pays.geometry.notna() & pays['GID_0'].notna()].reset_index(drop=True)
codes = sorted(pays['GID_0'].unique())
idx = {c: i + 1 for i, c in enumerate(codes)}
print('polygones : %d ; pays : %d' % (len(pays), len(codes)), flush=True)
with rasterio.open(TIF % 2015) as src:
    W, H, T = src.width, src.height, src.transform
PAS = abs(T.a)
cid = rasterize(((g, idx[c]) for g, c in zip(pays.geometry, pays['GID_0'])),
                out_shape=(H, W), transform=T, fill=0, dtype='uint16', all_touched=False)
print('rasterise en %.0f s ; pixels avec pays : %d' % (time.time() - t0, int((cid > 0).sum())),
      flush=True)

# ---- les pixels habites hors polygone, sur l'une quelconque des sept periodes --------
hors = np.zeros((H, W), dtype=bool)
pop15 = None
for an in ANS:
    with rasterio.open(TIF % an) as src:
        for r0 in range(0, H, 1200):
            h = min(1200, H - r0)
            v = src.read(1, window=Window(0, r0, W, h))
            hors[r0:r0 + h] |= (cid[r0:r0 + h] == 0) & np.isfinite(v) & (v > 0)
rr, cc = np.nonzero(hors)
xs = T.c + (cc + 0.5) * T.a
ys = T.f + (rr + 0.5) * T.e
print('pixels habites hors polygone : %d (%.0f s)' % (len(rr), time.time() - t0), flush=True)

# ---- distance au polygone le plus proche, bornee au pas puis a 1 degre ----------------
arbre = shapely.STRtree(pays.geometry.values)
pts = shapely.points(xs, ys)
dist = np.full(len(pts), np.inf)
proche = np.full(len(pts), -1, dtype=np.int64)
for borne in (PAS, 1.0):
    reste = np.nonzero(~np.isfinite(dist))[0]
    if not len(reste):
        break
    (ip, ig), d = arbre.query_nearest(pts[reste], max_distance=borne, return_distance=True,
                                      all_matches=False)
    dist[reste[ip]] = d
    proche[reste[ip]] = ig
    print('  distance <= %.5f deg : %d pixels (%.0f s)' % (borne, len(ip), time.time() - t0), flush=True)

# ---- temoin AVANT toute ecriture ------------------------------------------------------
pp = np.zeros(len(rr))                    # population 2015 des pixels hors polygone
with rasterio.open(TIF % 2015) as src:
    monde = 0.0
    for r0 in range(0, H, 1200):
        h = min(1200, H - r0)
        v = src.read(1, window=Window(0, r0, W, h))
        v = np.where(np.isfinite(v) & (v > 0), v, 0.0)
        monde += float(v.sum())
        dans = (rr >= r0) & (rr < r0 + h)
        pp[dans] = v[rr[dans] - r0, cc[dans]]
# Conservation des quantites intermediaires (lot 48) : TOUS les pixels hors polygone,
# avec leur distance et leur polygone le plus proche, AVANT le temoin -- qu'il accepte ou
# refuse, la question suivante ne doit pas couter un recalcul.
if LISTE:
    pd.DataFrame({'x': xs, 'y': ys, 'pop2015': pp, 'dist_deg': dist,
                  'GID_0_proche': [pays['GID_0'].iat[g] if g >= 0 else '' for g in proche]}
                 ).to_csv(LISTE.replace('.csv', '_tous.csv'), index=False)
## Anomalies connues (author's decision (a), lot 49) : pixels au-dela du pas, arretes
## par leur position dans la trame. Tout NOUVEAU pixel au-dela du pas fait refuser.
_ac = pd.read_csv(os.path.join(R, 'tools', 'anomalies_rattachement_connues.csv'))
CONNUES = set(zip(_ac['ligne'].astype(int), _ac['colonne'].astype(int)))
temoin_rattachement(xs, ys, pp, dist, PAS, monde, sortie=LISTE,
                    cles=list(zip(rr.astype(int), cc.astype(int))), connues=CONNUES)

# ---- rattachement, puis agregation par (cellule, pays) --------------------------------
ok = dist < PAS
cid[rr[ok], cc[ok]] = np.array([idx[pays['GID_0'].iat[g]] for g in proche[ok]], dtype='uint16')
NP = len(codes) + 1
blocs_col = np.arange(W) // 120
res = []
for an in ANS:
    with rasterio.open(TIF % an) as src:
        for r0 in range(0, H, 120):
            h = min(120, H - r0)
            v = src.read(1, window=Window(0, r0, W, h)).astype('float64')
            v[~np.isfinite(v) | (v < 0)] = 0.0
            c = cid[r0:r0 + h, :].astype('int64')
            cle = (np.broadcast_to(blocs_col, (h, W)).astype('int64') * NP + c).ravel()
            s = np.bincount(cle, weights=v.ravel(), minlength=(W // 120 + 1) * NP)
            nz = np.nonzero(s)[0]
            bc, pi = nz // NP, nz % NP
            k = bc < 360
            y = 88.5 - (r0 // 120)
            for b, p, val in zip(bc[k], pi[k], s[nz][k]):
                res.append((round(-179.5 + b, 1), y, codes[p - 1] if p else '', an, val))
    print('%d fait (%.0f s)' % (an, time.time() - t0), flush=True)
df = pd.DataFrame(res, columns=['x', 'y', 'GID_0', 'year', 'pop'])
df.to_csv(OUT, index=False)
print('ecrit %s : %d lignes, %.0f s' % (OUT, len(df), time.time() - t0))
