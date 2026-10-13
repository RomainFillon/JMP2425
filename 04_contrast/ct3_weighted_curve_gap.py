# -*- coding: utf-8 -*-
r"""
weighted_curve_gap_north.py — ventilation par tranche de temperature de l'ecart
entre les deux fonctions de dommage, PONDERE PAR LE CHANGEMENT D'EXPOSITION
nordique, et non brut.

POURQUOI CE SCRIPT EXISTE
-------------------------
Le manuscrit (Section 4, DRF journaliere) affirmait que la difference entre les
fonctions que les deux metriques recuperent « is concentrated above 20-25 C ».
Mesure sur les deux courbes publiees, l'ecart BRUT dit le contraire : sur le
support commun aux deux bras, 83 % de l'ecart absolu vit a T <= 10 C et 8 % a
T >= 20 C, avec un maximum de 2,00 log-point au bord froid contre 0,09 a 25 C.

Mais l'ecart brut n'est pas ce qui porte le contraste. Un desaccord de 2 log-
points a -8 C multiplie par un changement de comptage quasi nul ne pese rien ;
0,09 a 25 C multiplie par une forte hausse du nombre de jours peut dominer. Ce
script mesure l'objet pondere, qui est le seul dont la phrase puisse parler.

CE QUE CE SCRIPT CALCULE
------------------------
Geographie et poids : REPRIS A L'IDENTIQUE de ct1_delta_north_grid.R
sect. 3 (cellules 1 deg de latitude y >= 45N, poids population 2005 de
amen_prod.csv, base 2015 et futur 2095 tous deux issus de CMIP6_tas.csv).

Deux mesures signees d'exposition, chacune sur SON axe, toutes deux de masse
totale nulle :

    ds(T)   = sum_c omega_c [ c_fut(T ; c) - c_base(T ; c) ]      (axe : bin journalier)
    dg(Tbar) = mesure a atomes : +omega_c en Tbar_c^fut, -omega_c en Tbar_c^base
                                                                  (axe : moyenne annuelle)

Les increments nordiques de Test 1 s'ecrivent alors exactement

    Delta^daily  = sum_T  ds(T)    f_d(T)
    Delta^annual = sum_Tb dg(Tbar) f_a(Tbar)

et le contraste se telescope EXACTEMENT en deux termes :

    Delta^daily - Delta^annual
        = sum_T ds(T) [ f_d(T) - f_a(T) ]                    <- terme A : ECART DE COURBE
        + [ sum_T ds(T) f_a(T) - sum_Tb dg(Tbar) f_a(Tbar) ] <- terme B : ECART DE MESURE

Le terme A est l'objet demande : l'ecart entre les deux fonctions pondere par le
changement d'exposition journaliere. Il est VENTILABLE par tranche de T, et c'est
sa ventilation qui dit si la phrase du manuscrit est exacte ou non.

Le terme B n'est pas un residu de calcul : c'est le fait que les deux bras
n'integrent meme pas contre la meme mesure d'exposition (une distribution de
jours d'un cote, une moyenne annuelle par cellule de l'autre). Il est rapporte,
jamais tu : le rapporter est ce qui rend la decomposition exacte plutot
qu'approchee.

INVARIANCE DE NORMALISATION (verifiee par assertion)
----------------------------------------------------
ds et dg sont de masse totale nulle, donc A et B sont invariants a l'ajout d'une
constante a f_d ou a f_a. Les courbes publiees sont normalisees a leur propre
maximum ; le resultat ne depend pas de ce choix. L'assertion le prouve au lieu de
l'affirmer.

CE QUE CE SCRIPT N'EST PAS
--------------------------
Ce n'est pas un re-run de ct1. ct1 evalue w'beta specification par
specification dans la base d'ESTIMATION ; ici on prend les deux courbes PUBLIEES
(AIC-ponderees, grille entiere), qui sont exactement les objets dont la phrase du
manuscrit parle et ceux que climate_variability.r consomme. Le contraste de
Test 1 sur la specification de reference (spline4 wl0.5, prod) vaut -0,1099 ; il
est imprime comme repere d'ordre de grandeur, PAS comme assertion : les deux
objets ne vivent pas dans la meme base et n'ont aucune raison de coincider au
chiffre pres.

Ce n'est pas non plus une mesure sur le contraste GE. Le GE realloue, et la
transposition n'est pas monotone (contrat v5, sect. 0 du 31/07). Ce script mesure
l'incrementau niveau de la fonction de dommage, ce qui est le niveau auquel la
phrase auditee est ecrite.

SORTIE
------
output/pe_weighted_gap_north.csv         une ligne par (source de courbe, tranche de T)
output/pe_weighted_gap_north_scalars.csv les dix scalaires des macros \EMPgap*
Nom descriptif, sans conclusion (contrat v5, sect. 0, corollaire de nommage).

Migre au paquet de replication — 2026-09-04. Voir l'en-tete de
ct1_delta_north_grid.R pour le dossier complet du changement de garde-fou.

Romain Fillon — 2026-08-06
"""
import hashlib
import os
import sys
import time

# 🚨 STDOUT EN UTF-8, POSE LE 2026-09-07 — SINON L'ECHEC EST UN FAUX ECHEC.
# Sous Windows, stdout redirige vers un fichier ou un pipe (ce que fait
# run_all.sh) prend la codepage ANSI, cp1252 sur cette machine. Un caractere
# hors de cette page fait lever UnicodeEncodeError a print(), et ce script en
# imprime quatre lignes non-ASCII dont un U+26A0 (mesure du 07/09).
# ⛔ CE QUI REND LA CHOSE VICIEUSE : le print fautif vient APRES les deux
# ecritures. Les CSV etaient corrects sur le disque, les chiffres justes, et le
# script rendait quand meme 1 -- que l'orchestrateur classe en echec d'etape.
# Un diagnostic pose sur le calcul n'aurait rien trouve.
# errors non force : on veut que l'encodage marche, pas qu'il masque.
for _flux in (sys.stdout, sys.stderr):
    try:
        _flux.reconfigure(encoding="utf-8")
    except (AttributeError, OSError):   # flux non reconfigurable : on n'insiste pas
        pass

import numpy as np
import pandas as pd

# Racine deduite de l'emplacement de CE fichier : il vit dans <racine>/04_contrast/.
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("JMP_ROOT") or os.path.dirname(HERE)
REPO = REPO.replace("\\", "/")

# Les constantes de grille viennent de 00_parameters.py, le jumeau Python de
# 00_parameters.R. Ce script en lisait UNE en dur -- la specification de
# reference du pre-engagement (P4) -- et c'etait la seconde copie d'une
# designation qui doit etre unique : ct2 portait la premiere, en R. Depuis le
# 18/09/2026 les deux la lisent.
sys.path.insert(0, os.path.join(REPO, "01_estimation"))
from importlib import import_module as _import_module   # noqa: E402
_PARAMS = _import_module("00_parameters")

# --- constantes reprises de ct1_delta_north_grid.R -------------------------
BASE_YEAR = 2015
FUT_YEAR = 2095
NORTH_LAT = 45.0
POP_YEAR = 2005

CMIP6 = "intermediate_outputs/CMIP6_tas.csv"
# amen_prod.csv est ECRIT PAR LE PAQUET (01d_crh_triangulation.R l.332). Hors du
# paquet ce script lisait l'exemplaire de l'archive CRH ; ici on lit le producteur.
AMEN = "intermediate_outputs/crh/amen_prod.csv"
CURVE_D = "output/damage_prod_bin.csv"
CURVE_A = "output/damage_prod_annual_bin.csv"
MC_D = "output/MC_damage_prod_daily_quantiles.csv"
MC_A = "output/MC_damage_prod_annual_quantiles.csv"
# Les bornes de support et les reperes du Test 1 ne sont JAMAIS retapes ici :
# ils sont DERIVES de la sortie de ct1. Depuis la migration du 2026-09-04, ct1
# est une etape du paquet, donc ce fichier a un producteur dans le paquet.
T1 = "output/contribA/test1_delta_north.csv"

# --- LA PROVENANCE, ET CE QU'ELLE REMPLACE  (migration 2026-09-04) -----------
# Hors du paquet, ces sept entrees etaient EPINGLEES a des md5 en dur, et un
# ecart etait un refus (sys.exit(2)). Le danger vise etait precis, et il faut le
# nommer avant de changer de garde : l'APPARIEMENT CROISE. Le terme A evalue la
# courbe ANNUELLE a des temperatures de BIN JOURNALIER, et la repartition entre A
# et B depend des bornes de support, qui viennent de ct1. Combiner des courbes
# d'un millesime avec des bornes d'un autre donne une decomposition qui
# TELESCOPE EXACTEMENT — B absorbe tout — et qu'aucun controle de somme ne peut
# donc detecter.
#
# ⛔ DES md5 GELES NE CONVIENNENT PLUS ICI. Les six entrees sont desormais
# produites par des etapes numerotees du meme run_all.sh ; tout rejeu legitime
# ferait echouer l'epinglage, qui deviendrait le garde-fou qu'on desarme.
#
# CE QUI REMPLACE : l'EMPREINTE DE PROVENANCE que ct1 inscrit en colonne de
# test1_delta_north.csv. On la recalcule sur les six memes entrees et on refuse
# si elle differe. Le refus porte donc sur le danger vise, et pas sur un proxy.
#
# ⚠ CE QUE L'EMPREINTE NE VOIT PAS : deux artefacts du meme run estimes sous des
# SPECIFICATIONS differentes — par exemple une courbe estimee sans controles
# gravitaires et des bornes issues d'une grille qui les porte. Cela se traite en
# amont, en ecrivant les controles a chaque site d'appel (correctif de 01a du
# 2026-09-04), pas ici.
FP_INPUTS = [
    "intermediate_outputs/BACI_raw_pairs_annual.rds",
    "intermediate_outputs/climate_long_annual.rds",
    "intermediate_outputs/step2_migration.rds",
    "intermediate_outputs/climate_long.rds",
    CMIP6,
    AMEN,
]
# Les noms COURTS servent au jeton : ct1 les forme en R a partir de basename-like
# cles, et les deux implementations doivent tomber sur la meme chaine.
FP_KEYS = {
    "intermediate_outputs/BACI_raw_pairs_annual.rds": "BACI_raw_pairs_annual.rds",
    "intermediate_outputs/climate_long_annual.rds":   "climate_long_annual.rds",
    "intermediate_outputs/step2_migration.rds":       "step2_migration.rds",
    "intermediate_outputs/climate_long.rds":          "climate_long.rds",
    "intermediate_outputs/CMIP6_tas.csv":             "CMIP6_tas.csv",
    "intermediate_outputs/crh/amen_prod.csv":         "amen_prod.csv",
}
# Les quatre artefacts de 01a que ce script consomme. Ils ne sont pas dans
# l'empreinte (ct1 ne les lit pas) mais leur ABSENCE doit etre nommee.
ARTEFACTS_01A = [CURVE_D, CURVE_A, MC_D, MC_A]

T1_REF_SPEC = _PARAMS.GRID_REF_SPEC   # 00_parameters.py, pre-engagement (P4)

# Tranches de temperature de la ventilation (bornes superieures incluses).
BANDS = [(-np.inf, 0), (0, 5), (5, 10), (10, 15), (15, 20),
         (20, 25), (25, 30), (30, np.inf)]

# --- LIMITE CONNUE DU TERME A, ET SON CONTROLE ----------------------------
# Le terme A evalue la courbe ANNUELLE a des temperatures de BIN JOURNALIER.
# Or f_a est une fonction de la moyenne annuelle : l'evaluer a -30 C n'a pas de
# sens empirique, et sur la grille publiee elle y est ecretee a plat. Le
# telescopage reste exact (B absorbe le choix), mais la REPARTITION entre A et B
# depend de cette extension arbitraire. On mesure donc A aussi en restreignant au
# support d'estimation du bras annuel, et on verifie que la conclusion ne bouge
# pas. Bornes DERIVEES de output/contribA/test1_delta_north.csv, colonnes
# support_lo/support_hi du bras annuel (min et max sur les specifications de
# la grille). Ces bornes ne dependent que de la WINSORISATION, pas de la forme :
# le passage de la grille de 24 a 15 specifications le 18/09/2026 retire des
# formes, jamais un seuil, et laisse donc les deux enveloppes inchangees.


def annual_support():
    """(large, etroit) : les deux enveloppes du support d'estimation annuel."""
    d = pd.read_csv(os.path.join(REPO, T1))
    a = d[(d.channel == "prod") & (d.arm == "annual")]
    return {"le plus large":  (float(a.support_lo.min()), float(a.support_hi.max())),
            "le plus etroit": (float(a.support_lo.max()), float(a.support_hi.min()))}


def t1_reference():
    """Increments nordiques Test 1 sur la specification de reference."""
    d = pd.read_csv(os.path.join(REPO, T1))
    r = d[(d.channel == "prod") & (d.spec_label == T1_REF_SPEC)]
    return (float(r[r.arm == "daily"].delta_north.iloc[0]),
            float(r[r.arm == "annual"].delta_north.iloc[0]))


ANNUAL_SUPPORT = None   # renseigne par main(), depuis l'artefact

# --- non-regression numerique : VALEURS, jamais une conclusion -------------
# Renseignees apres le premier run. On asserte les deux termes du telescopage et
# la part de A au-dessus de 20 C ; un ecart signale une derive de pipeline, il ne
# prejuge pas du sens du resultat.


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fingerprint():
    """Le jeton de provenance, forme EXACTEMENT comme ct1 le forme en R :
    lignes "cle=md5" triees par cle, jointes par un saut de ligne, sans saut
    final. Les deux implementations doivent tomber sur la meme valeur, sinon le
    refus se declencherait sur une difference d'encodage et non de provenance."""
    parts = {FP_KEYS[rel]: _md5(os.path.join(REPO, rel)) for rel in FP_INPUTS}
    txt = "\n".join(f"{k}={parts[k]}" for k in sorted(parts))
    return hashlib.md5(txt.encode("utf-8")).hexdigest()


def check_inputs():
    """Presence, puis coherence de provenance avec ct1. Deux refus distincts :
    une entree absente est une etape qui n'a pas tourne ; une empreinte qui
    diverge est un appariement croise."""
    besoin = FP_INPUTS + ARTEFACTS_01A + [T1]
    absents = [r for r in besoin if not os.path.exists(os.path.join(REPO, r))]
    if absents:
        print("[REFUS] entree(s) absente(s) :")
        for r in absents:
            qui = ("ct1_delta_north_grid.R" if r == T1
                   else "01a_productivity_drf.R" if r in ARTEFACTS_01A
                   else "voir run_all.sh")
            print(f"  {r}")
            print(f"      produit par : {qui}")
        sys.exit(2)

    fp = fingerprint()
    t1 = pd.read_csv(os.path.join(REPO, T1), nrows=1)
    if "inputs_fingerprint" not in t1.columns:
        print("[REFUS] test1_delta_north.csv ne porte pas inputs_fingerprint :")
        print("  artefact anterieur a la migration du 2026-09-04. Relancer ct1")
        print("  plutot que de desarmer ce controle.")
        sys.exit(2)
    t1_fp = str(t1["inputs_fingerprint"].iloc[0])
    if t1_fp != fp:
        print("[REFUS] APPARIEMENT CROISE.")
        print(f"  ct1 a tourne sur l'empreinte {t1_fp}")
        print(f"  ct3 lit aujourd'hui l'empreinte {fp}")
        print("  Les bornes de support et les courbes ne viendraient pas du meme")
        print("  etat d'entree. La decomposition telescoperait exactement malgre")
        print("  tout : B absorbe le residu, donc AUCUN controle de somme ne peut")
        print("  voir l'erreur. Relancer ct1 sur l'etat courant.")
        sys.exit(2)
    print(f"[ok] provenance : empreinte {fp}, identique a celle de ct1.")
    return fp


def load_north_exposure():
    """Reproduit ct1_delta_north_grid.R sect. 3. Retourne (ds, muw, ncell).

    ds  : Series indexee par temperature entiere, masse-jour nordique signee.
    muw : DataFrame (x, y, omega, mu_base, mu_fut).
    """
    path = os.path.join(REPO, CMIP6)
    keep = []
    for chunk in pd.read_csv(path, chunksize=4_000_000,
                             dtype={"x": "float32", "y": "float32",
                                    "year": "int32", "value": "int32",
                                    "count": "float64"}):
        m = (chunk["y"] >= NORTH_LAT) & chunk["year"].isin([BASE_YEAR, FUT_YEAR])
        if m.any():
            keep.append(chunk.loc[m])
    cm = pd.concat(keep, ignore_index=True)
    del keep

    # part de jours dans le bin, a l'interieur de la cellule-annee
    tot = cm.groupby(["x", "y", "year"], sort=False)["count"].transform("sum")
    cm["c"] = cm["count"] / tot

    pop = pd.read_csv(os.path.join(REPO, AMEN), usecols=["year", "lat", "lon", "pop_x"])
    pop = pop.loc[pop["year"] == POP_YEAR, ["lat", "lon", "pop_x"]]
    pop = pop.rename(columns={"lon": "x", "lat": "y", "pop_x": "pop"})
    pop["pop"] = pd.to_numeric(pop["pop"], errors="coerce")
    pop["x"] = pop["x"].astype("float32")
    pop["y"] = pop["y"].astype("float32")
    cm = cm.merge(pop, on=["x", "y"], how="inner")

    # une cellule ne compte que si elle a SES DEUX annees
    ny = cm.groupby(["x", "y"], sort=False)["year"].nunique()
    ok2 = ny[ny == 2].index
    cm = cm.set_index(["x", "y"]).loc[ok2].reset_index()
    ncell = len(ok2)

    wcell = cm[["x", "y", "pop"]].drop_duplicates()
    wcell["omega"] = wcell["pop"] / wcell["pop"].sum()
    assert abs(wcell["omega"].sum() - 1) < 1e-12
    cm = cm.merge(wcell[["x", "y", "omega"]], on=["x", "y"])

    # DAILY : masse-jour signee par temperature
    cm["sgn"] = np.where(cm["year"] == FUT_YEAR, 1.0, -1.0)
    ds = cm.assign(m=cm["omega"] * cm["sgn"] * cm["c"]).groupby("value")["m"].sum().sort_index()

    # ANNUAL : moyenne annuelle par cellule
    mu = cm.assign(vt=cm["value"] * cm["c"]).groupby(["x", "y", "year"])["vt"].sum().rename("mu")
    muw = mu.reset_index().pivot(index=["x", "y"], columns="year", values="mu").reset_index()
    muw = muw.rename(columns={BASE_YEAR: "mu_base", FUT_YEAR: "mu_fut"})
    muw = muw.merge(wcell[["x", "y", "omega"]], on=["x", "y"])
    muw = muw[np.isfinite(muw["mu_base"]) & np.isfinite(muw["mu_fut"])]

    print(f"  cellules nordiques retenues (base ET futur) : {ncell}")
    print(f"  somme des masses signees ds (doit etre ~0)  : {ds.sum():.3e}")
    print(f"  Tbar nordique : base [{muw.mu_base.min():.1f}, {muw.mu_base.max():.1f}] "
          f" futur [{muw.mu_fut.min():.1f}, {muw.mu_fut.max():.1f}] degC")
    return ds, muw, ncell


def read_curve(rel, col):
    """Lit (temperature, valeur) dans un artefact de courbe.

    Le controle de colonne est explicite : un KeyError pandas ne nomme ni
    l'artefact ni ce qu'il porte, et c'est exactement l'erreur qui a arrete
    ct3 le 07/09 -- il demandait "damage", nom des artefacts hors paquet,
    quand 01a ecrit "f_z_raw". Une assertion de FORME (nombre de lignes,
    bornes) ne l'aurait pas vue : seule une assertion de COLONNE l'attrape.
    """
    d = pd.read_csv(os.path.join(REPO, rel))
    for need in ("temperature_C", col):
        if need not in d.columns:
            raise SystemExit(
                f"[STOP] {rel} : colonne {need!r} absente.\n"
                f"  colonnes presentes : {list(d.columns)}\n"
                f"  Ce n'est pas un probleme de chemin : l'artefact existe mais\n"
                f"  ne porte pas ce que ce script attend. Verifier le producteur.")
    return d["temperature_C"].to_numpy(float), d[col].to_numpy(float)


def evaluate(ds, muw, xd, fd, xa, fa, label, rows):
    """Telescopage exact + ventilation du terme A. Ajoute les lignes a `rows`."""
    T = ds.index.to_numpy(float)
    w = ds.to_numpy(float)

    fd_T = np.interp(T, xd, fd)
    fa_T = np.interp(T, xa, fa)

    # Delta^daily et Delta^annual
    delta_daily = float(w @ fd_T)
    fa_fut = np.interp(muw["mu_fut"].to_numpy(float), xa, fa)
    fa_base = np.interp(muw["mu_base"].to_numpy(float), xa, fa)
    om = muw["omega"].to_numpy(float)
    delta_annual = float(om @ (fa_fut - fa_base))

    gap = fd_T - fa_T
    A = float(w @ gap)
    B = float(w @ fa_T) - delta_annual
    contrast = delta_daily - delta_annual

    # --- l'exactitude du telescopage se PROUVE, elle ne s'affirme pas -----
    assert abs((A + B) - contrast) < 1e-12, f"telescopage rompu : {A+B} vs {contrast}"

    # --- invariance de normalisation : ajouter des constantes ne change rien -
    A2 = float(w @ ((fd_T + 3.7) - (fa_T - 1.9)))
    B2 = float(w @ (fa_T - 1.9)) - float(om @ ((fa_fut - 1.9) - (fa_base - 1.9)))
    assert abs(A2 - A) < 1e-12 and abs(B2 - B) < 1e-12, "invariance de normalisation violee"

    print(f"\n=== {label} ===")
    print(f"  Delta^daily  = sum ds . f_d      : {delta_daily:+.5f}")
    print(f"  Delta^annual = sum dg . f_a      : {delta_annual:+.5f}")
    print(f"  contraste (daily - annual)       : {contrast:+.5f}")
    print(f"    terme A  ECART DE COURBE pondere : {A:+.5f}  ({100*A/contrast:5.1f} % du contraste)")
    print(f"    terme B  ECART DE MESURE         : {B:+.5f}  ({100*B/contrast:5.1f} % du contraste)")

    tot_abs = float(np.abs(w * gap).sum())
    print(f"\n  Ventilation du terme A par tranche de temperature :")
    print(f"  {'tranche':>12} {'ds (masse)':>12} {'ecart moyen':>12} "
          f"{'contrib A':>11} {'part |A|':>9}")
    above20 = 0.0
    for lo, hi in BANDS:
        m = (T > lo) & (T <= hi)
        if not m.any():
            continue
        contrib = float(w[m] @ gap[m])
        dsm = float(w[m].sum())
        gm = float(np.average(gap[m], weights=np.abs(w[m]))) if np.abs(w[m]).sum() > 0 else np.nan
        share = 100 * float(np.abs(w[m] * gap[m]).sum()) / tot_abs if tot_abs else np.nan
        lbl = (f"<= {hi:g}" if lo == -np.inf else
               f"> {lo:g}" if hi == np.inf else f"{lo:g}-{hi:g}")
        print(f"  {lbl:>12} {dsm:+12.5f} {gm:+12.4f} {contrib:+11.5f} {share:8.1f} %")
        if lo >= 20:
            above20 += abs(float((w[m] * gap[m]).sum()))
        rows.append({"curve_source": label, "band_lo": lo, "band_hi": hi,
                     "exposure_change": dsm, "mean_curve_gap": gm,
                     "contrib_A": contrib, "share_abs_A_pct": share})
    share20 = 100 * above20 / tot_abs if tot_abs else np.nan
    print(f"  -> part de |A| a T > 20 C : {share20:.1f} %")
    rows.append({"curve_source": label, "band_lo": np.nan, "band_hi": np.nan,
                 "exposure_change": np.nan, "mean_curve_gap": np.nan,
                 "contrib_A": A, "share_abs_A_pct": share20})

    # --- controle de support : la conclusion depend-elle de l'extrapolation ? --
    print("  Controle — A restreint au support d'estimation du bras annuel :")
    share20_supp = {}
    for lbl, (lo, hi) in ANNUAL_SUPPORT.items():
        m = (T >= lo) & (T <= hi)
        A_in = float(w[m] @ gap[m])
        tin = float(np.abs(w[m] * gap[m]).sum())
        m20 = m & (T > 20)
        s20 = 100 * float(np.abs(w[m20] * gap[m20]).sum()) / tin if tin else np.nan
        share20_supp[lbl] = s20
        print(f"    support {lbl:<14} [{lo:+.1f} ; {hi:+.1f}] : "
              f"A_dans = {A_in:+.5f} ({100*A_in/A:3.0f} % de A) — "
              f"part de |A| a T > 20 C : {s20:.1f} %")
        rows.append({"curve_source": f"{label} / support annuel {lbl}",
                     "band_lo": lo, "band_hi": hi, "exposure_change": np.nan,
                     "mean_curve_gap": np.nan, "contrib_A": A_in,
                     "share_abs_A_pct": s20})
    below0 = 100 * float(np.abs(w[T <= 0] * gap[T <= 0]).sum()) / tot_abs
    return {"A": A, "B": B, "contrast": contrast, "share20": share20,
            "share_below0": below0, "share20_supp": share20_supp,
            "delta_daily": delta_daily, "delta_annual": delta_annual}


def cold_side_ratio(T, w, FD, FA, xd, fd, xa, fa, rows):
    """Combien de fois plus l'annuel valorise-t-il la MEME perte de jours froids ?

    Deux mesures, parce qu'aucune ne suffit seule :

    (a) PONDEREE — [somme ds.f_a] / [somme ds.f_d] sur la bande froide. C'est
        l'objet du papier (meme changement d'exposition, deux courbes), mais une
        contribution de BANDE n'est pas invariante a la normalisation : somme(ds)
        n'est pas nulle sur une bande. Les deux courbes etant normalisees a leur
        propre optimum, la lecture est « benefice attribue, relatif au mieux
        atteignable » — la meme convention que tab:spec_invariance.

    (b) SANS NORMALISATION — la montee de chaque courbe entre deux temperatures
        situees DANS le support des deux bras. Une difference de niveaux ne
        depend d'aucune normalisation. Sert de controle a (a).

    /!\\ La restriction au support d'estimation du bras annuel AUGMENTE le
    rapport (elle retire les bins ou f_a est ecretee a plat). La borne
    CONTRAIGNANTE de l'argument est donc la plus BASSE, pas la plus haute :
    c'est elle qui doit etre citee comme plancher.
    """
    print("\n=== De combien l'annuel survalorise-t-il la perte de jours froids ? ===")
    print("  (a) pondere par le changement d'exposition, bande T <= 0 C")
    out = {}
    variants = [("aucune restriction", -np.inf, np.inf)]
    variants += [(f"support annuel {k}", lo, hi) for k, (lo, hi) in ANNUAL_SUPPORT.items()]
    for lab, lo, hi in variants:
        m = (T <= 0) & (T >= lo) & (T <= hi)
        a = float(w[m] @ FA[m])
        d = float(w[m] @ FD[m])
        print(f"      {lab:<32} annuel {a:+.5f}  journalier {d:+.5f}  -> x{a/d:.2f}")
        out[lab] = a / d
        rows.append({"curve_source": f"ratio froid pondere / {lab}",
                     "band_lo": lo, "band_hi": 0.0, "exposure_change": float(w[m].sum()),
                     "mean_curve_gap": np.nan, "contrib_A": a - d,
                     "share_abs_A_pct": 100 * a / d})

    print("  (b) sans normalisation : montee de chaque courbe, bornes dans le support")
    for a0, a1 in ((-8, 0), (-8, 2), (-5, 5)):
        ra = float(np.interp(a1, xa, fa) - np.interp(a0, xa, fa))
        rd = float(np.interp(a1, xd, fd) - np.interp(a0, xd, fd))
        print(f"      {a0:+g} -> {a1:+g} C{'':<24} annuel {ra:+.4f}  journalier {rd:+.4f}"
              f"  -> x{ra/rd:.2f}")
        out[f"montee {a0}->{a1}"] = ra / rd
        rows.append({"curve_source": f"ratio froid sans normalisation / {a0} a {a1} C",
                     "band_lo": a0, "band_hi": a1, "exposure_change": np.nan,
                     "mean_curve_gap": np.nan, "contrib_A": ra - rd,
                     "share_abs_A_pct": 100 * ra / rd})
    lo_bound = min(out.values())
    print(f"  -> PLANCHER sur toutes les variantes : x{lo_bound:.2f} "
          f"(extremite contraignante, contrat v5 sect. 0)")
    out["plancher"] = lo_bound
    return out


def main():
    global ANNUAL_SUPPORT, FP
    t0 = time.perf_counter()
    FP = check_inputs()
    ANNUAL_SUPPORT = annual_support()
    print("[ok] bornes de support annuel derivees de " + T1 + " : "
          + " / ".join(f"{k} [{lo:+.2f} ; {hi:+.2f}]"
                       for k, (lo, hi) in ANNUAL_SUPPORT.items()))
    print("\n--- geographie nordique (miroir de ct1_delta_north_grid.R sect. 3) ---")
    ds, muw, _ = load_north_exposure()

    rows = []
    # "f_z_raw" et non "damage" : c'est le nom que 01a donne a la courbe
    # normalisee a max = 0 (01a_productivity_drf.R l.686). "damage" etait le
    # nom des artefacts hors paquet.
    xd, fd = read_curve(CURVE_D, "f_z_raw")
    xa, fa = read_curve(CURVE_A, "f_z_raw")
    main_res = evaluate(ds, muw, xd, fd, xa, fa,
                        "courbes publiees (AIC-ponderees)", rows)

    T = ds.index.to_numpy(float)
    ratios = cold_side_ratio(T, ds.to_numpy(float),
                             np.interp(T, xd, fd), np.interp(T, xa, fa),
                             xd, fd, xa, fa, rows)

    xdm, fdm = read_curve(MC_D, "q50")
    xam, fam = read_curve(MC_A, "q50")
    evaluate(ds, muw, xdm, fdm, xam, fam,
             "MC q50 (prod daily corrige, non consomme par le GE)", rows)

    t1d, t1a = t1_reference()
    print(f"\n[repere, NON assertif] Test 1 {T1_REF_SPEC} (base d'estimation) : "
          f"daily {t1d:+.5f} / annual {t1a:+.5f} / contraste {t1d - t1a:+.5f}")

    # --- scalaires destines aux macros du manuscrit -------------------------
    # Une macro se cable sur une COLONNE, jamais sur un rapport de session.
    T = ds.index.to_numpy(float)
    w = ds.to_numpy(float)
    scal = [
        ("A_share_of_contrast_pct", 100 * main_res["A"] / main_res["contrast"],
         "part du contraste portee par l'ecart de courbe"),
        ("share_absA_below0_pct", main_res["share_below0"],
         "part de |A| sous 0 degC"),
        ("share_absA_above20_pct", main_res["share20"],
         "part de |A| au-dessus de 20 degC"),
        ("share_absA_above20_supp_min_pct",
         min(main_res["share20_supp"].values()),
         "idem, restreint au support annuel — borne basse"),
        ("share_absA_above20_supp_max_pct",
         max(main_res["share20_supp"].values()),
         "idem — borne haute"),
        ("expo_lost_below0_pp", 100 * abs(float(w[T <= 0].sum())),
         "masse-jour nordique perdue sous 0 degC, en points"),
        ("expo_gained_above20_pp", 100 * float(w[T > 20].sum()),
         "masse-jour nordique gagnee au-dessus de 20 degC, en points"),
        ("cold_ratio_floor", ratios["plancher"],
         "sur-valorisation annuelle du froid — PLANCHER, extremite contraignante"),
        ("cold_ratio_supp_min", min(ratios[k] for k in ratios
                                    if k.startswith("support")),
         "idem, pondere et restreint au support — borne basse"),
        ("cold_ratio_supp_max", max(ratios[k] for k in ratios
                                    if k.startswith("support")),
         "idem — borne haute"),
    ]
    out_s = os.path.join(REPO, "output", "pe_weighted_gap_north_scalars.csv")
    df_s = pd.DataFrame(scal, columns=["key", "value", "description"])
    # La provenance voyage AVEC les scalaires : un consommateur peut prouver que
    # ces dix nombres et les bornes de ct1 viennent du meme etat d'entree.
    df_s["inputs_fingerprint"] = FP
    df_s.to_csv(out_s, index=False)
    print(f"[write] output/pe_weighted_gap_north_scalars.csv ({len(scal)} scalaires)")

    out = os.path.join(REPO, "output", "pe_weighted_gap_north.csv")
    df_r = pd.DataFrame(rows)
    df_r["inputs_fingerprint"] = FP
    df_r.to_csv(out, index=False)
    print(f"\n[write] output/pe_weighted_gap_north.csv ({len(rows)} lignes)")

    # ---- TEMOIN DE REPRODUCTION, RAPPORTE ET NON BLOQUANT ------------------
    # Voir reference_repro_contrast/README.md pour pourquoi ce n'est plus un
    # refus. En deux lignes : les courbes ont ete reestimees, donc un refus ici
    # serait le garde-fou qu'on desarme ; et le danger que le refus visait
    # (l'appariement croise) est traite par check_inputs(), qui lui s'arrete.
    wit = os.path.join(HERE, "reference_repro_contrast",
                       "weighted_curve_gap_20260806.csv")
    got = {"A_total": main_res["A"], "B_total": main_res["B"],
           "A_share_above20_pct": main_res["share20"],
           "A_share_above20_supp_large_pct":
               main_res["share20_supp"]["le plus large"],
           "A_share_above20_supp_etroit_pct":
               main_res["share20_supp"]["le plus etroit"],
           "cold_ratio_supp_large": ratios["support annuel le plus large"],
           "cold_ratio_supp_etroit": ratios["support annuel le plus etroit"],
           "cold_ratio_plancher": ratios["plancher"]}
    if not os.path.exists(wit):
        print(f"\n[temoin] ABSENT : {os.path.relpath(wit, REPO)} — ecart non mesure.")
    else:
        ref = pd.read_csv(wit).set_index("key")["value"].to_dict()
        print("\n--- ecart au temoin gele du 06/08/2026 (RAPPORTE, non bloquant) ---")
        ecarts = []
        for k in sorted(ref):
            if k not in got:
                print(f"  {k:34s} absent de cette execution")
                continue
            d = got[k] - ref[k]
            ecarts.append(abs(d))
            print(f"  {k:34s} temoin {ref[k]:+9.5f}   obtenu {got[k]:+9.5f}"
                  f"   ecart {d:+9.5f}")
        if ecarts:
            print(f"\n  ecart absolu maximal : {max(ecarts):.5f}")
            print("  ⚠ Un ecart non nul n'est pas un echec : il mesure de combien la")
            print("    reestimation a deplace la decomposition. Les dix macros \\EMPgap*")
            print("    du manuscrit (l.567-581) en descendent — ce paragraphe est a")
            print("    relire, et le temoin n'est PAS a rafraichir pour le faire taire.")
    print("")
    print(f"Temps total: {(time.perf_counter() - t0) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
