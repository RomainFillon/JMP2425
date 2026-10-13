#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
02f_fig_mig_latitude.py -- latitudinal profile of initial net migration.

Reads:  intermediate_outputs/{mig_matrix_final_raw.h5, indices_drop.csv,
        identifier_xy_final.csv}
Writes: output/fig_M1_net_by_latitude.png
        output/migration_net_by_latitude.csv

--------------------------------------------------------------------------
WHY THIS SCRIPT EXISTS (09/09/2026)
--------------------------------------------------------------------------
The manuscript loads `fig_M1_net_by_latitude.png` and no step of the package
produced it. The closing inventory of run_all.sh already said so, classifying
it "replication_package/ -- NON MIGRE": the figure is the external-validation
evidence of the migration block, and the sentence that cites it says the
latitudinal profile "directly motivates the welfare narrative". Removing the
figure would have removed the only geographic check on the migration
calibration, so it was migrated instead.

Ported from code_inari/validation/validate_migration.py (M1 path only,
l.145-158 and l.225-231). What is NOT ported, and why:
  - Fig M2, the mobility map: it needs a Natural Earth shapefile that the
    package does not ship, and the manuscript does not load it.
  - the country-level validation table: it needs an external net-migration
    series (JMP_NET_OBS), which is unset and undocumented.
  - the country attribution through lengths_country.csv: the latitudinal
    profile does not use it. Carrying it would import a dependency for a
    column nothing reads.

--------------------------------------------------------------------------
🚨 ORIENTATION -- THE TRAP OF THIS FILE, AND THE ONE THE ORIGINAL FELL INTO
--------------------------------------------------------------------------
`mig_matrix_final_raw.h5` holds RAW flows in persons, cumulated over five
years, on the FULL 17,048-cell grid. The net flow of a cell is

    net(r) = inflow(r) - outflow(r)

and it changes SIGN under a transpose. A wrongly oriented read produces a
perfectly plausible mirror-image profile, with no error raised anywhere.

`validate_migration.py` l.15-17 declares it keeps "l'orientation NATIVE R
(ligne = origine)" and reads `f[:]` WITHOUT `.T`. Measured on 09/09/2026, that
is wrong for this file, and the transpose IS applied here. hdf5r writes
column-major, so h5py reads the matrix transposed -- which is exactly what
`data_cleaning.load_h5` and `load_migration` document and correct for.

⛔ AND IT IS NOT A MATTER OF CONVENTION, IT IS DECIDED BY A PHYSICAL BOUND. A
cell cannot send more people than it holds. Counting, per orientation, the
cells whose flow exceeds their own population:

    h5py rows     624 cells over population   (worst: 447,254 leaving a cell of 7)
    h5py columns    0 cells over population

and the asymmetry survives every population cut: at pop > 100,000 it is still
3 against 0. The bounded quantity is therefore the COLUMN sum, so h5py columns
are origins, h5py rows are destinations, and `.T` restores rows = origin.

CLOSURE IS CHECKED, BUT IT PROVES LESS THAN IT LOOKS. These flows are internal,
so the net sums to zero -- but `sum(colsums) - sum(rowsums)` is identically
`sum(S) - sum(S)` for ANY submatrix, so it cannot detect a mis-indexed subset,
and it is invariant to the transpose that matters here. It is kept as an
ARITHMETIC check only, and it earned its place: it caught a real defect, the
float32 accumulator. numpy sums a float32 array in float32, and pairwise
summation over 17,048 terms drifts to ~1.4e-05 relative -- so the accumulation
is forced to float64 below. The orientation guard is the one that carries the
scientific claim.

⚠️ indices_drop.csv is 1-BASED (written by R); main.py l.278 does the same
subtraction.
"""
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))            # .../02_GE/
PKG = os.environ.get('JMP_ROOT') or os.path.dirname(HERE)    # .../replication/
INPUTS = os.path.join(PKG, 'intermediate_outputs')
OUT = os.path.join(PKG, 'output')

## Width of the latitude bands, in degrees. 5 is the convention of the original
## exhibit and of the caption; the bar height below tracks it.
LAT_BIN = 5.0
## Tolerance of the closure check, as a share of gross flow. The matrix is
## float32 on 17,048^2 entries, so exact cancellation is not available.
CLOSURE_TOL = 1e-6


def main():
    need = ['mig_matrix_final_raw.h5', 'indices_drop.csv',
            'identifier_xy_final.csv']
    missing = [n for n in need if not os.path.exists(os.path.join(INPUTS, n))]
    if missing:
        print('[STOP] intrant(s) absent(s) de intermediate_outputs/ : %s'
              % ', '.join(missing))
        print('       Producteur : 02b_mig_flows.R (matrice brute et ids),')
        print('       02a_calibration_at_1deg.r (identifier_xy_final.csv).')
        print('       fig_M1_net_by_latitude.png NON produite.')
        return 2

    import h5py

    if not os.path.exists(os.path.join(INPUTS, 'population0.csv')):
        print('[STOP] population0.csv absent : la garde d orientation ne peut')
        print('       pas tourner, et sans elle le signe de chaque barre est')
        print('       indecidable. Producteur : 02a_calibration_at_1deg.r.')
        return 2

    ids = pd.read_csv(os.path.join(INPUTS, 'indices_drop.csv')).values.flatten()
    ids = ids.astype(int)
    pos = ids - 1                      # 1-based (R) -> 0-based
    xy = pd.read_csv(os.path.join(INPUTS, 'identifier_xy_final.csv'))
    xy.columns = [c.strip().strip('"') for c in xy.columns]
    pix = pd.DataFrame({'id': ids}).merge(xy[['id', 'x', 'y']], on='id', how='left')
    if pix['y'].isna().any():
        print('[STOP] %d cellules actives sans coordonnee dans '
              'identifier_xy_final.csv.' % int(pix['y'].isna().sum()))
        print('       Voir le piege de renumerotation documente au README.')
        return 2
    print('[load] %d cellules actives | latitudes [%.1f, %.1f]'
          % (len(pix), pix['y'].min(), pix['y'].max()))

    pop = pd.read_csv(os.path.join(INPUTS, 'population0.csv'))
    P = pd.DataFrame({'id': ids}).merge(pop, on='id', how='left')['value'] \
        .values.astype(float)

    with h5py.File(os.path.join(INPUTS, 'mig_matrix_final_raw.h5'), 'r') as f:
        key = 'mig_flows' if 'mig_flows' in f else list(f.keys())[0]
        print('[load] %s %s %s' % (key, f[key].shape, f[key].dtype))
        M = f[key][:]
    M = np.nan_to_num(M, nan=0.0)
    M = M[pos][:, pos]

    ## dtype=float64 EXPLICITE. numpy somme un tableau float32 en float32, et la
    ## sommation par paires sur 17 048 termes derive a ~1,4e-05 en relatif. Sans
    ## ce dtype le controle de fermeture echoue sur de la precision seule.
    s_rows = M.sum(axis=1, dtype=np.float64)
    s_cols = M.sum(axis=0, dtype=np.float64)
    gross = float(M.sum(dtype=np.float64))
    del M

    ## ---- GARDE D'ORIENTATION : la borne physique --------------------------
    ## Une cellule ne peut pas emettre plus d'habitants qu'elle n'en compte. Le
    ## sortant est donc borne par la population, l'entrant ne l'est pas. C'est
    ## ce qui distingue les deux orientations, et rien d'autre ici ne le fait.
    v_rows = int((s_rows > P).sum())
    v_cols = int((s_cols > P).sum())
    print('[orientation] cellules dont le flux depasse leur population :')
    print('    lignes h5py   %5d' % v_rows)
    print('    colonnes h5py %5d' % v_cols)
    if v_cols < v_rows:
        outflow, inflow, sens = s_cols, s_rows, 'colonnes = origines (.T requis)'
    elif v_rows < v_cols:
        outflow, inflow, sens = s_rows, s_cols, 'lignes = origines (lecture native)'
    else:
        print('[STOP] les deux orientations violent la borne autant de fois.')
        print('       Le sens des flux est indecidable sur cet artefact, et le')
        print('       signe de chaque barre en depend. Aucune figure.')
        return 2
    print('    => %s' % sens)
    net = inflow - outflow

    ## ---- CONTROLE ARITHMETIQUE : fermeture ---------------------------------
    ## ⚠️ CE CONTROLE PROUVE MOINS QU'IL N'EN A L'AIR, et le dire ici evite qu'on
    ## s'y fie : sur flux internes le net somme a zero, mais l'egalite est une
    ## TAUTOLOGIE pour tout sous-ensemble et elle est invariante a la
    ## transposition. Elle ne valide pas la selection des cellules, et surtout
    ## pas l'orientation -- c'est la garde ci-dessus qui le fait. Gardee parce
    ## qu'elle attrape ce dont elle est capable : la derive d'accumulation.
    resid = abs(float(net.sum())) / max(gross, 1.0)
    print('[arithmetique] flux brut %.2f M pers | net total / brut = %.2e'
          % (gross / 1e6, resid))
    if resid > CLOSURE_TOL:
        print('[STOP] le net ne somme pas a zero a %.0e pres.' % CLOSURE_TOL)
        print('       Cause attendue : accumulation en float32. Aucune figure.')
        return 2
    print('               == FERME ==')

    ## Bandes de latitude, etiquetees par leur CENTRE.
    pix['net'] = net
    pix['lat_bin'] = np.floor(pix['y'] / LAT_BIN) * LAT_BIN + LAT_BIN / 2.0
    prof = pix.groupby('lat_bin', as_index=False)['net'].sum().sort_values('lat_bin')
    os.makedirs(OUT, exist_ok=True)
    prof.to_csv(os.path.join(OUT, 'migration_net_by_latitude.csv'), index=False)

    ## Ce que la phrase du manuscrit affirme, IMPRIME : gains nets entre +25 et
    ## +45 degres, pertes sous les tropiques. Le lecteur du log voit donc le
    ## chiffre que la figure montre, sans avoir a ouvrir le PDF.
    band = prof[(prof.lat_bin >= 25) & (prof.lat_bin <= 45)]['net'].sum()
    trop = prof[prof.lat_bin.abs() <= 23.5]['net'].sum()
    peak = prof.loc[prof['net'].idxmax(), 'lat_bin']
    print('[lecture] +25 a +45 N : %+.2f M | tropiques (|lat|<=23,5) : %+.2f M'
          % (band / 1e6, trop / 1e6))
    print('          bande de gain maximal : %+.1f degres' % peak)

    fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    ax.barh(prof['lat_bin'], prof['net'] / 1e6,
            height=LAT_BIN * 0.8, color='#b03a2e')
    ax.axvline(0, color='k', lw=0.8)
    ax.axhspan(-23.5, 23.5, color='grey', alpha=0.15)
    ax.set_xlabel('Net migration (millions, 5-yr)')
    ax.set_ylabel('Latitude (°)')
    ax.set_title('Initial net migration by latitude band')
    ax.grid(alpha=0.3, axis='x')
    f_out = os.path.join(OUT, 'fig_M1_net_by_latitude.png')
    fig.savefig(f_out, dpi=130)
    plt.close(fig)
    print('[write] %s' % f_out)
    print('[write] %s' % os.path.join(OUT, 'migration_net_by_latitude.csv'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
