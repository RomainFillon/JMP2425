"""gardes_lot74.py -- les predicats des gardes posees au lot 74 (01/10/2026).

Fonctions PURES : elles disent si une condition tient, et pourquoi pas. Les
emetteurs (02_GE/emit_values.py, 02_GE/emit_values_empirics.py) et l'etape
01c_graine2 de run_all.sh les appellent ; tools/test_gardes_lot74.py les exerce
sur des connus-faux. Une garde dont le predicat vit dans un script non importable
ne se teste qu'en rejouant la chaine ; ici, chacune se teste en une seconde.

Chaque predicat rend (ok, detail) : ok booleen, detail une chaine qui nomme les
valeurs comparees, pour le message de refus ou la ligne du mode rapport.
"""
import os
import re

import pandas as pd

# ---------------------------------------------------------------- (11) zero signe
_RE_ZERO_SIGNE = re.compile('^' + re.escape('$') + r'?\s*[+\-−]0(\.0+)?\s*' + re.escape('$') + '?$')


def zero_signe(val):
    """Vrai si la valeur s'imprime comme un zero signe (« -0.000 », « +0.000 »,
    « $-0.000$ »). Un zero sans signe (« 0.000 ») n'est pas vise."""
    return bool(_RE_ZERO_SIGNE.match(str(val).strip()))


def decimales_non_nulles(v, d0=3, dmax=6):
    """Formate v avec un signe, a d0 decimales, et en ajoute tant que le chiffre
    imprime est nul (dmax au plus). Meme regle que fb() d'emit_values.py, sans
    les dollars. Ne s'applique qu'a une valeur dont le signe est l'information."""
    for d in range(d0, dmax + 1):
        t = f"{v:+.{d}f}"
        if float(t) != 0:
            return t
    return f"{v:+.{dmax}f}"


# ---------------------------------------------------------------- (7) B.2
def garde_b2(identifies, n_id):
    """Supplement B.2 nomme le schema exporter x year comme le seul identifie."""
    ok = ('exporter' in set(identifies)) and int(n_id) == 1
    return ok, 'identifies = {%s}, EMPohSchemesId = %s' % (', '.join(sorted(identifies)), n_id)


# ---------------------------------------------------------------- (8) graines
def garde_graines(q95_base, q95_noy, move_max):
    """B.2 : « the largest such move ..., within which the baseline's percentile
    lies and the other's does not ».
    LOT 85o (07/10/2026) : la clause « marge du baseline de l'ordre de son propre
    ecart entre graines (0,5x a 2x) » est RETIREE. Elle gardait la phrase du §2.3
    « crosses zero by a margin of the order of the difference between two
    independent seeds », retiree du .tex au lot 81 (rapport_lot85n §1a)."""
    a = q95_base <= move_max < q95_noy
    return a, ('q95 baseline %.6f <= max gap %.6f < q95 noy %.6f: %s'
               % (q95_base, move_max, q95_noy, a))


# ---------------------------------------------------------------- (9) amenite
def garde_amenite_ge(annual_q50, extra_q50):
    """§3.4 et supplement : a l'estimation interieure le bras annuel projette
    une perte pour le Nord, et le supplement journalier est negatif."""
    ok = annual_q50 < 0 and extra_q50 < 0
    return ok, 'GEnorthAnnualQfifty %.6f < 0, GEnorthExtraQfifty %.6f < 0' % (annual_q50, extra_q50)


def garde_pnonlin_annuel(pmin):
    """Supplement : « carries no non-linearity the data support »."""
    return pmin > 0.05, 'min pval_nonlin annual %.6f > 0.05' % pmin


# ---------------------------------------------------------------- (10) arctique
def garde_seuil_arctique(v, all_neg):
    """§3.2 : « the daily total is negative at every latitude threshold of the grid,
    and within rounding of zero above 66.5°N ».
      - « every latitude threshold » : all_neg (= \\GEnorthLatAllNeg) est vrai ;
      - « within rounding of zero » : |v| < 0,0005 (a trois decimales), v = total
        journalier central >= 66,5N.
    LOT 85q (07/10/2026) : la clause all_neg est ajoutee avec la phrase re-ecrite (la
    version precedente citait \\GEnorthLatLastNeg, « up to and including »)."""
    a = all_neg is True
    b = abs(v) < 0.0005
    return a and b, ('every threshold negative (GEnorthLatAllNeg): %s; '
                     '|central daily total >= 66.5N| = %.6f < 0.0005: %s'
                     % (bool(all_neg), abs(v), b))


# ---------------------------------------------------------------- (18) graine 2
def comparer_numerique(ref_dir, new_dir, fichiers, tol=1e-8):
    """Compare des CSV cellule par cellule : memes colonnes, memes lignes, memes
    textes, nombres egaux a `tol` pres (absolu). Rend (ok, [lignes de detail])."""
    det, ok = [], True
    for f in fichiers:
        pr, pn = os.path.join(ref_dir, f), os.path.join(new_dir, f)
        if not (os.path.exists(pr) and os.path.exists(pn)):
            det.append('%s: MISSING (%s)' % (f, pr if not os.path.exists(pr) else pn))
            ok = False
            continue
        a, b = pd.read_csv(pr), pd.read_csv(pn)
        if list(a.columns) != list(b.columns) or a.shape != b.shape:
            det.append('%s: different SHAPE %s / %s' % (f, a.shape, b.shape))
            ok = False
            continue
        emax = 0.0
        # Lot 85j : une ligne par colonne fautive. Jusqu'au 85i `bad` etait ecrase
        # a chaque colonne et seule la DERNIERE etait nommee (le 05/10 : slope_q95
        # citee, slope_q05 et slope_q50 passees sous silence alors qu'elles bougeaient aussi).
        bad = []
        for c in a.columns:
            if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]) \
                    and not pd.api.types.is_bool_dtype(a[c]):
                e = float((a[c] - b[c]).abs().max()) if len(a) else 0.0
                if not (e <= tol):
                    bad.append('%s: column %s gap %r > %g' % (f, c, e, tol))
                emax = max(emax, e)
            elif not a[c].astype(str).equals(b[c].astype(str)):
                bad.append('%s: column %s differs (text or boolean)' % (f, c))
        if bad:
            det.extend(bad)
            ok = False
        else:
            det.append('%s: identical to within %g (max gap %.2e)' % (f, tol, emax))
    return ok, det


if __name__ == '__main__':
    # Usage (etape 01c_graine2) : python tools/gardes_lot74.py comparer REF NEW f1 f2 ...
    import sys
    if len(sys.argv) >= 5 and sys.argv[1] == 'comparer':
        ok, det = comparer_numerique(sys.argv[2], sys.argv[3], sys.argv[4:])
        for d in det:
            print('  [graine2] ' + d)
        if not ok:
            print('[STOP] 01c_graine2: the output does not reproduce the reference to within 1e-8.')
        sys.exit(0 if ok else 1)
    print(__doc__)
    sys.exit(2)
