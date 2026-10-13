# -*- coding: utf-8 -*-
"""
02e_north_contrast_band.py -- Northern welfare contrast (sc2 - sc3) at the
quantiles of the damage function, computed from THIS package's own GE grid.

WHAT IT COMPUTES
----------------
One object, evaluated at each quantile that the grid ran:

    contrast = population-weighted mean, over cells with y >= 45N,
               of the 2095 cumulative sum of dlog(U_sc2 / U_sc3)

    sc2 = daily curve evaluated on the DAILY DISTRIBUTION
    sc3 = annual curve evaluated at the ANNUAL MEAN

Same geography, same weights, same cumulation as the "daily-extra" of
tab_ge_north and the "total" of decompose_metric_channels.py.

WHY THIS FILE REPLACES replication_package/ge_analysis/ge_north_contrast_band.py
--------------------------------------------------------------------------------
The former producer lived OUTSIDE the package and read nine `welfare_outer*`
dumps. Those dumps came from an OUTER FIXED-POINT LOOP over utility that the
package model NO LONGER HAS: 02c_model/main.py initialises `iter_outer` and
`i_outer` (l. ~432) and never increments or tests them again. Its only loop is
the local population fixed point. The solver is a single forward sweep over t,
which is why a run takes ~350 s where the old one took hours.

CONSEQUENCE, MEASURED ON 09/09/2026 AND ACTED ON:
  - the point estimate and its bounds ARE derivable here, from welfare_sc2_*
    and welfare_sc3_* -- no dump, no re-run;
  - the three macros that described the OSCILLATION of the outer loop
    (\\GEbandQnfiveNneg, \\GEbandQnfiveLeast, \\GEbandOsc) describe a solver
    this package does not contain. They are NOT re-derived here, and the
    sentence that used them has been removed from the manuscript. Reproducing
    them would have meant re-importing the old solver, which would move every
    GE number in the paper -- a re-specification, not a repair.

DUMP_OUTER IS A DEAD FLAG. main.py reads it and prints it; nothing acts on it.
Do not reach for it expecting welfare_outer* to appear.

THREE DECIMALS ARE REPORTED, TWO ARE PUBLISHED. The old header forbade a third
decimal because the outer loop did not converge. That reason is gone with the
loop, but the replacement reason is not weaker: what this paper establishes about
the contrast is its SIGN, and it reports the magnitude to an order of magnitude.
The CSV carries the rounding the manuscript uses; the log prints four decimals so
that a drift is visible before it reaches a macro.

GUARDRAIL POLARITY. The hard stop belongs to NUMERICAL DRIFT -- a missing input,
a shape mismatch, a grid whose runs do not belong together. It NEVER belongs to a
scientific conclusion: if the sign of the contrast changed, this script must
DISPLAY it, not die on it.

Reads:  output/welfare_sc2_{q05,q50,q95}.csv, output/welfare_sc3_{q05,q50,q95}.csv
        output/ge_runs_manifest.csv
        intermediate_outputs/{indices_drop,identifier_xy_final,population0}.csv
Writes: output/ge_north_contrast_band.csv
"""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
OUT = os.path.join(PKG, 'output')
TMP = os.path.join(PKG, 'intermediate_outputs')

LAT_MIN = 45.0
CUM_LABEL = '2095'          # the last row of the cumulated path
QUANTILES = ('q05', 'q50', 'q95')


def _read_matrix(path):
    d = pd.read_csv(path)
    d = d.loc[:, ~d.columns.str.contains('^Unnamed')]
    return d.values


def geography():
    """North mask and population weights, on the GE cell ordering.

    indices_drop.csv gives the cells the model kept, IN ORDER. The welfare
    matrices are (T-1, N) on that same ordering, so the mask must be built from
    it and never from a filter on identifier_xy_final.csv, whose row order is a
    different thing entirely.
    """
    ids = pd.read_csv(os.path.join(TMP, 'indices_drop.csv')).iloc[:, 0].astype(int).values
    xy = pd.read_csv(os.path.join(TMP, 'identifier_xy_final.csv'))
    xy.columns = [c.strip('"') for c in xy.columns]
    pop = pd.read_csv(os.path.join(TMP, 'population0.csv'))
    id2y = dict(zip(xy['id'].astype(int), xy['y']))
    id2p = dict(zip(pop['id'].astype(int), pop['value']))
    missing = [i for i in ids if i not in id2y or i not in id2p]
    if missing:
        sys.exit(f"[STOP] {len(missing)} identifiants de indices_drop.csv absents de "
                 f"identifier_xy_final.csv ou population0.csv (ex. {missing[:5]}). "
                 "Les trois fichiers ne viennent pas du meme millesime de calibration.")
    lat = np.array([id2y[i] for i in ids], dtype=float)
    wpop = np.array([id2p[i] for i in ids], dtype=float)
    return lat >= LAT_MIN, wpop


def contrast(f2, f3, north, wpop):
    w2, w3 = _read_matrix(f2), _read_matrix(f3)
    if w2.shape != w3.shape:
        sys.exit(f"[STOP] formes incompatibles : {os.path.basename(f2)} {w2.shape} "
                 f"vs {os.path.basename(f3)} {w3.shape}.")
    if w2.shape[1] != north.size:
        sys.exit(f"[STOP] {w2.shape[1]} colonnes contre {north.size} cellules gardees. "
                 "Le welfare et indices_drop.csv ne decrivent pas la meme grille.")
    if not (np.all(w2 > 0) and np.all(w3 > 0)):
        sys.exit("[STOP] welfare non strictement positif : le log n'est pas defini.")
    v = np.cumsum(np.log(w2 / w3), axis=0)[-1, :]
    return float(np.average(v[north], weights=wpop[north]))


def check_one_grid(needed):
    """The triplet must come from ONE grid under ONE set of settings.

    decompose_metric_channels.py refuses to run without this check because a
    triplet assembled from runs that do not belong together reverses the sign of
    both channels. The same exposure applies here.
    """
    man = os.path.join(OUT, 'ge_runs_manifest.csv')
    if not os.path.exists(man):
        sys.exit("[STOP] output/ge_runs_manifest.csv absent. Sans lui rien ne prouve "
                 "que sc2 et sc3 viennent du meme run : c'est exactement la faute "
                 "que l'ancien producteur pinglait au md5, fichier par fichier.")
    m = pd.read_csv(man).set_index('file')
    absent = [f for f in needed if f not in m.index]
    if absent:
        sys.exit("[STOP] absents du manifeste : " + ', '.join(absent))
    sub = m.loc[list(needed)]
    ##########################################################################
    # `run_id` RETIRE DE LA CLE D'UNICITE. Mesure du 21/09/2026.
    #
    # Cette garde exigeait un `run_id` unique sur les six fichiers. Sous le
    # lanceur actuel, `run_one_ge.sh` estampille UN run_id PAR RUN :
    # `one_sc2_q05_20260921_114403`, `one_sc3_q05_20260921_115221`, ... Mesure
    # sur le manifeste du 21/09 : 78 runs, 78 run_id distincts. La condition
    # `nunique() == 1` sur six fichiers differents etait donc STRUCTURELLEMENT
    # INSATISFAISABLE -- pas une garde qui mord a tort, une garde qui ne peut
    # JAMAIS etre verte.
    #
    # Elle etait satisfaisable avant le 10/09, quand un pilote de lot posait un
    # run_id commun (`run_20260910_092031` sur les six). La convention a change
    # sous elle et personne ne l'a vu, parce que le producteur n'a plus tourne
    # depuis.
    #
    # LA CLE JUSTE EST L'EMPREINTE D'INTRANTS. `inputs_md5` hache SCALE_PROD
    # depuis le 20/09 : deux conventions ne peuvent plus partager une empreinte,
    # et c'est exactement ce que `assert_one_arm` (emit_values.py l.128) verifie
    # -- lui imprime `HETEROGENE:` sur les run_id et continue. L'en-tete de
    # `assert_one_arm` affirme d'ailleurs employer « la meme cle que
    # 02e_north_contrast_band.py » : la premisse etait fausse, elle devient
    # vraie ici.
    #
    # Les run_id restent IMPRIMES. Une heterogeneite n'est pas une panne, mais
    # elle doit rester lisible : c'est ce qui permet de dire de quelle passe
    # vient chaque colonne quand on relit un log.
    ##########################################################################
    # ⛔ `fillna` AVANT DE COMPTER, ET CE N'EST PAS UNE PRECAUTION DE STYLE.
    # `amenity_spec` vaut la chaine « null » pour le bras a f_a = 0, et
    # pd.read_csv la convertit en NaN (elle est dans les na_values par defaut).
    # `nunique()` IGNORE les NaN : une bande melangeant le bras interieur
    # (`baseline`) et le bras du corps (`null`) rendait donc nunique() == 1 et
    # passait au VERT. C'est exactement le melange que ce producteur existe
    # pour refuser. Trouve par mutation le 21/09/2026 ; `assert_one_arm`
    # (emit_values.py l.160) portait deja le fillna, pas celui-ci.
    for col in ('inputs_md5', 'cap_colmass', 'amenity_spec'):
        vals = sub[col].fillna('(vide=null)').astype(str)
        if vals.nunique() != 1:
            sys.exit(f"[STOP] les runs de la bande ne partagent pas {col} : "
                     f"{sorted(set(vals))}. Un contraste entre deux millesimes "
                     "n'est pas un contraste.")
    _runs = sorted(set(sub['run_id']))
    print(f"[ok] un seul millesime : inputs_md5={sub['inputs_md5'].iloc[0][:8]}, "
          f"cap={sub['cap_colmass'].iloc[0]}, amenity={sub['amenity_spec'].iloc[0]}")
    print("[ok] run_id : " + (_runs[0] if len(_runs) == 1
                              else f"{len(_runs)} passes -> " + ', '.join(_runs)))


def main():
    pairs = {q: (f'welfare_sc2_{q}.csv', f'welfare_sc3_{q}.csv') for q in QUANTILES}
    have = {q: fs for q, fs in pairs.items()
            if all(os.path.exists(os.path.join(OUT, f)) for f in fs)}
    if 'q50' not in have or 'q95' not in have:
        sys.exit("[STOP] q50 et q95 sont l'un et l'autre necessaires : le manuscrit "
                 "les cite comme l'estimation interieure et sa borne. "
                 "Presents : " + (', '.join(sorted(have)) or 'aucun'))
    # A mechanism that can do nothing must say whether it acted.
    for q in QUANTILES:
        if q not in have:
            print(f"[!] {q} : welfare absent, quantile NON rapporte.")

    needed = [f for fs in have.values() for f in fs]
    check_one_grid(needed)

    north, wpop = geography()
    print(f"[ok] NORD >= {LAT_MIN:.0f}N : {north.sum()} cellules, "
          f"{wpop[north].sum() / 1e6:.0f} M habitants (population de base)")

    rows = []
    for q in sorted(have):
        f2, f3 = (os.path.join(OUT, f) for f in have[q])
        c = contrast(f2, f3, north, wpop)
        print(f"  {q} : contraste = {c:+.4f}  (publie a 2 decimales : {c:+.2f})")
        rows.append({'quantity': f'corrected_{q}', 'value': f'{c:+.2f}'.replace('+', ''),
                     'value_raw': f'{c:.6f}',
                     'source': f'{have[q][0]} vs {have[q][1]}, cumul {CUM_LABEL}, '
                               f'moyenne ponderee population y>={LAT_MIN:.0f}N'})

    # The SIGN is reported as an observation, never asserted. A guardrail that
    # encodes a conclusion turns a finding into a crash.
    signs = {r['quantity']: float(r['value_raw']) for r in rows}
    if all(v < 0 for v in signs.values()):
        print("[obs] signe negatif a tous les quantiles rapportes, "
              "conforme au resultat publie.")
    else:
        pos = [k for k, v in signs.items() if v >= 0]
        print(f"[obs] SIGNE POSITIF OU NUL sur {pos}. Ce n'est pas une panne de "
              "pipeline, c'est un fait a examiner avant publication.")

    df = pd.DataFrame(rows)
    dest = os.path.join(OUT, 'ge_north_contrast_band.csv')
    df.to_csv(dest, index=False)
    print(f"[write] output/ge_north_contrast_band.csv ({len(df)} lignes)")


if __name__ == '__main__':
    main()
