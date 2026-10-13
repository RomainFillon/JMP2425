#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
emit_values_empirics.py -- writes output/values_empirics.tex.

Emits the LaTeX macros carrying the purely EMPIRICAL numbers of the paper:
Stage 1 diagnostics, the Oster bound, eta, Delta-rho, the CRH triangulation,
the Akaike weights and the within R-squared. Every macro is prefixed \EMP.
No \GE macro is emitted here; those belong to emit_values.py.

Reads:  output/*.csv produced by the estimation steps, and a few rendered
        tables when no machine-readable source exists
Writes: output/values_empirics.tex
Exit:   1 when at least one macro mismatches its expected value, 0 otherwise

How a value gets here:
  - it is PARSED from a generated artefact, a .csv or a rendered .tex, whenever
    a machine source exists;
  - a value with no machine source is emitted as a LITERAL, marked [LITERAL]
    in the report, so that it is visible rather than silent.

The script prints a SUBSTITUTION REPORT: macro by macro, the value it emits
against the value expected. A MISMATCH means the macro is NOT emitted. A macro
that quietly replaced a hard-coded number by a different one would move a
published figure without anyone reading a diff.
"""
import hashlib
import math
import os
import re
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))          # .../02_GE/
PKG  = os.path.dirname(HERE)                               # .../replication/
REPO = PKG
INPUTS     = os.path.join(PKG, 'intermediate_outputs')
_RUNS_ENV  = os.environ.get('RUNS_DIR')
RUNS_DIR   = os.path.join(PKG, 'output')
VALUES_TEX = os.path.join(PKG, 'output', 'values.tex')
OUT = os.path.join(PKG, 'output')

# La GRILLE DE CONTRASTE vient de 00_parameters.py, le jumeau Python de
# 00_parameters.R. Elle porte la decision d'auteur du 18/09/2026 -- quinze
# specifications -- et elle ne descend d'aucun run : c'est ce qui en fait un
# tiers utilisable pour confronter contribA_grid.csv, qui, lui, descend de ct1.
sys.path.insert(0, os.path.join(PKG, '01_estimation'))
from importlib import import_module as _import_module   # noqa: E402
_PARAMS = _import_module('00_parameters')


def read(path):
    with open(path, encoding='utf-8') as fh:
        return fh.read()

# --dry-run parses, reports, and does NOT write output/values_empirics.tex.
# Emission is a single global pass run after validation, so the report can be
# read without producing an artefact.
DRY = '--dry-run' in sys.argv

## MODE RAPPORT (lot 53, 28/09/2026). Pendant une REPRISE DELIBEREE declaree par
## `tools/reference_phase.py --declarer-reprise`, un ecart a la valeur attendue est
## RAPPORTE (statut ECART-REPRISE, macro emise, ligne dans ecarts_reprise_empirics.tsv)
## au lieu de bloquer. Aucun argument ni variable d'environnement n'ouvre ce mode.
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'tools'))
from reprise import reprise_declaree, bandeau   # noqa: E402
REPRISE = reprise_declaree()
if REPRISE:
    print(bandeau(REPRISE, 'emit_values_empirics.py'))


def garde_phrase(ok, garde, attendu, obtenu, message):
    """Lot 57 (author's decision, 28/09/2026) : une garde qui verifie une PHRASE du papier (et
    non un calcul ou la structure d'un fichier). Hors reprise : elle ARRETE, comme
    avant. En reprise declaree : son refus est une INFORMATION pour l'auteur -- il
    est rapporte avec ses deux valeurs (ecarts_reprise_gardes.tsv) et l'emission
    continue. Les gardes de structure et de coherence entre deux calculs restent
    bloquantes et n'appellent pas cette fonction."""
    if ok:
        return
    if not REPRISE:
        raise SystemExit(message)
    if DRY:
        print('[REPORT MODE] (dry-run, not written) guard "%s" -- expected %s, obtained %s'
              % (garde, attendu, obtenu))
    else:
        from reprise import rapporte   # noqa: E402
        rapporte('emit_values_empirics.py', garde, attendu, obtenu)

macros   = []   # (name, emitted value, value expected in the tex, source, status)
def add(name, value, expected, source, literal=False, regen=False):
    """expected = valeur en dur historique dans JMP_V9.tex (garde-fou de
    substitution). regen=True : la valeur PEUT legitimement avoir bouge depuis ce
    hard-code (regle (b) du runbook) -> la source fait foi, la macro est emise, et
    l'ecart au hard-code est REPORTE au lieu de bloquer.

    ⛔ UN MISMATCH N'EMET RIEN. La macro sort en `% \\NOM NOT EMITTED` et l'appel
    du .tex reste `Undefined control sequence`. Un MISMATCH ne se ferme donc PAS
    en rejouant la chaine : le rejeu reproduit le meme artefact, donc la meme
    valeur, donc le meme MISMATCH. Il se ferme en ARBITRANT la derive, puis en
    portant la valeur mesuree dans `expected`.

    ------------------------------------------------------------------------
    REBASAGE DU 09/09/2026 -- 25 `expected` deplaces, et pourquoi
    ------------------------------------------------------------------------
    Les `expected` etaient les litteraux de JMP_V9. L'ecart mesure la derive du
    paquet depuis, pas un bug. Ont ete rebasees UNIQUEMENT les macros dont
    l'artefact producteur est FRAIS -- 01a (run 14:21-14:35, posterieur au
    correctif GLORIA), 01b (aucune lecture GLORIA, donc intact), et la part de
    01e qui descend de 01b ou du climat brut.

    ⛔ RIEN N'A ETE REBASE sur un artefact PERIME : 01c, 01d, la chaine CRH,
    ct1-ct4, la part productivite de 01e et les onze runs GE datent d'avant le
    correctif. Les rebaser aurait fige une valeur que le prochain rejeu
    deplacera -- c'est-a-dire le meme travail deux fois.

    ⛔ ET RIEN N'A ETE REBASE quand la derive DEPLACE L'ENONCE. Cinq macros
    sont restees en MISMATCH volontaire, en attente d'arbitrage :
    EMPamenPjointLo, EMPamenAnnualPjointLo, EMPamenPnonlinMin (un p-value
    traverse ou frole le seuil de 5 %), EMPetaDailyAmen et EMPetaAnnualAmen
    (le contraste annuel/journalier du canal amenite s'annule).

    ⚠️ La cause portee en commentaire a chaque site est celle du CANAL et de
    son producteur. Elle n'est PAS une decomposition mesuree : departager les
    six causes connues demanderait des rejeux differentiels qui n'ont pas ete
    faits.
    """
    if regen:
        status = 'REGEN' if value == expected else f'REGEN (stale hard-code: {expected})'
    elif value == expected:
        status = 'MATCH'
    elif REPRISE and not literal:
        status = f'ECART-REPRISE (expected: {expected})'
    else:
        status = 'MISMATCH'
    if literal:
        status += ' [LITERAL]'
    macros.append((name, value, expected, source, status))

# ---------------------------------------------------------------- Stage 1 ----
# The producer CSV is read first, with the regex on the rendered table as a
# fallback, and the script SAYS which one it used. The CSV is keyed by a stable
# identifier; the regex matches DISPLAY LABELS and breaks as soon as the table
# changes shape. A fallback that fires silently is indistinguishable from a
# source that still works.
def _stage1(name):
    _p = os.path.join(OUT, 'stage1_%s_diagnostics.csv' % name)
    if os.path.exists(_p):
        return pd.read_csv(_p).set_index('key')['value'], 'stage1_%s_diagnostics.csv' % name
    print('[stage1] %s: CSV missing, falling back on the REGEX of the rendered table.' % name)
    return None, 'tab_stage1_%s.tex (REGEX, fallback)' % name

def _thousands(x):
    return format(int(round(float(x))), ',').replace(',', '{,}')


# regen=True cible par NOM, pas par canal. La boucle ci-dessous sert le commerce
# ET la migration ; un regen pose sur la boucle entiere aurait relache le
# garde-fou du canal commerce, qui n'a pas ete re-specifie.
#
# MOTIF (09/09/2026) : l'intrant de migration a change de MILLESIME. 01b, 02a et
# 02b lisaient bilat_mig_sex.csv, la serie PAR SEXE d'Abel & Cohen ; ce n'est pas
# une decomposition de la serie totale, `da_pb_closed` etant une comptabilite
# demographique en systeme FERME, et fermer par sexe puis sommer ne donne pas la
# fermeture du total. Les trois scripts lisent desormais bilat_mig.csv, flux
# TOTAUX, Figshare 10.6084/m9.figshare.7731233, VERSION 8 (IMS2024/WPP2024,
# publiee le 26/03/2025).
#
# ET LA SUBSTITUTION N'EST PAS LA CAUSE DE L'ECART. Mesure des trois etats :
#     archive 62.45  |  par sexe 54.88  |  flux totaux v8 53.69
# La bascule par-sexe -> totaux ne deplace que 1.2 point ; le millesime en
# deplace 7.6. `nobs` reste a 146 004 : MEME ECHANTILLON, valeurs revisees.
# La version 8 est ASSUMEE -- c'est ce que le DOI rend aujourd'hui, et epingler
# une version archivee d'un jeu regulierement revise vieillirait mal.
_REGEN_STAGE1 = {'EMPmigWaldF'}

for _nm, _tex, _macros in (
        ('trade', 'tab_stage1_trade.tex',
         # Rebases le 09/09/2026 sur le run 01a POSTERIEUR au correctif GLORIA
         # (14:21-14:35). Cause : le panel commerce -- repli sectoriel retire,
         # `sum()` du legacy corrige, et les deux pertes d'appariement GLORIA.
         # L'ancien couple etait 114.81 / 3{,}806{,}447 ; les 75 489 observations
         # perdues sont celles que la note de bas de page de la Section 4 chiffre
         # desormais (11,87 % des lignes).
         (('EMPtradeWaldF',     'wald_f',    '115.16',        lambda v: '%.2f' % float(v)),
          ('EMPtradeNobs',      'nobs',      '3{,}730{,}958', _thousands),
          ('EMPtradeNsectors',  'n_sectors', '20',            lambda v: str(int(float(v)))))),
        ('migration', 'tab_stage1_migration.tex',
         (('EMPmigWaldF',       'wald_f',    '62.45',         lambda v: '%.2f' % float(v)),
          ('EMPmigNobs',        'nobs',      '146{,}004',     _thousands)))):
    _d, _src = _stage1(_nm)
    _t = read(os.path.join(OUT, _tex))
    for _macro, _key, _exp, _fmt in _macros:
        if _d is not None and _key in _d.index:
            add(_macro, _fmt(_d[_key]), _exp, _src,
                regen=_macro in _REGEN_STAGE1)
        else:
            _rx = {'wald_f':    r'Wald \$F\$ \(gravity controls\) & ([\d.]+)',
                   'nobs':      r'Observations & ([\d,]+)',
                   'n_sectors': r'N\. sectors \(ISIC Rev\.3\) & (\d+)'}[_key]
            _m = re.search(_rx, _t)
            if _m is None:
                print('[stage1] %s: the REGEX of %s no longer matches -- macro NOT emitted.'
                      % (_macro, _tex))
                continue
            _v = _m.group(1)
            add(_macro, _v.replace(',', '{,}') if _key == 'nobs' else _v, _exp, _src)

# ------------------------------------------------- perte sectorielle GLORIA --
# Les trois grandeurs de la note de bas de page de la Section 4. Elles n'ont
# JAMAIS ete des litteraux du manuscrit : la note ne nommait qu'une cause
# d'exclusion sur trois et ne chiffrait pas la perte. Elles entrent donc sans
# valeur attendue historique -- `expected` porte la valeur MESUREE au run du
# 09/09/2026, ce qui fait du garde-fou une non-regression et non une
# confrontation a un hard-code qui n'a jamais existe.
_GL = os.path.join(OUT, 'trade_gloria_loss.csv')
if not os.path.exists(_GL):
    # Un mecanisme qui ne fait rien doit dire qu'il n'a rien fait. Aucune
    # valeur de repli : l'absence doit sortir en "Undefined control sequence"
    # a la compilation, pas en nombre plausible.
    print('[gloria] trade_gloria_loss.csv missing -> 3 EMPtradeGloriaLoss* macros '
          'NOT emitted. Producer: 01a_productivity_drf.R, section 3bis.')
else:
    _gl = pd.read_csv(_GL).iloc[0]
    _SRC_GL = 'trade_gloria_loss.csv (01a_productivity_drf.R, section 3bis)'
    add('EMPtradeGloriaLossObs',   f"{float(_gl['pct_rows']):.2f}",  '11.87', _SRC_GL)
    add('EMPtradeGloriaLossVal',   f"{float(_gl['pct_value']):.2f}", '9.08',  _SRC_GL)
    add('EMPtradeGloriaLossPairs', str(int(_gl['pairs_lost'])),      '191',   _SRC_GL)
    print('[gloria] %d divisions without a denominator (%s): %.2f %% of the rows, '
          '%.2f %% of the value, %d pairs'
          % (int(_gl['n_divisions_dropped']), _gl['divisions'],
             float(_gl['pct_rows']), float(_gl['pct_value']), int(_gl['pairs_lost'])))

# ---------------------------------------------------------------- Oster ------
# Read from the producer CSV rather than by regex on the rendered table. An
# emitter that parses a FORMATTED table depends on its layout: merging two
# columns is enough to break the match, and the failure surfaces as an
# AttributeError on a group that no longer exists. The CSV does not move when
# the table changes columns.
#
# LA CELLULE LUE EST LA REFERENCE DECLAREE, ET ELLE EST DERIVEE (09/09/2026).
# Les deux macros Oster visaient `spline4_wl0.5wu0.5`, tandis que le manuscrit
# les appelle « the reference specification ». La reference du paquet est
# declaree dans 00_parameters.R -- REF_SPEC / REF_DEG / REF_WLOWER / REF_WUPPER,
# soit polynomial degre 3 a 1 % -- et c'est une AUTRE cellule. Nommer
# « reference » une cellule qui ne l'est pas est le mislabel que ce paquet
# retire partout ailleurs.
#
# L'etiquette est RECONSTRUITE depuis 00_parameters.R plutot que retapee : un
# tag en dur ici redeviendrait faux au premier changement de reference, et
# silencieusement. Le motif suit celui de 01c l.331,
#   sprintf("%s%d_wl%swu%s", spec, degree, w$lower, w$upper)
# ⚠️ avec le piege deja documente plus bas : R rend 1.0 comme "1" et 0.5 comme
# "0.5", d'ou le formatage conditionnel.
def _ref_tag():
    _p = os.path.join(PKG, '01_estimation', '00_parameters.R')
    _s = read(_p)
    _g = {}
    for _k in ('REF_SPEC', 'REF_DEG', 'REF_WLOWER', 'REF_WUPPER'):
        _m = re.search(r'^\s*%s\s*<-\s*"?([A-Za-z0-9._]+?)"?L?\s*(?:#.*)?$' % _k,
                       _s, re.M)
        if _m is None:
            raise SystemExit(
                "[STOP] %s not found in 00_parameters.R.\n"
                "  The reference specification cannot be derived, and typing it\n"
                "  here would engrave the mislabel this block removes." % _k)
        _g[_k] = _m.group(1)

    def _w(v):
        f = float(v)
        return str(int(f)) if f == int(f) else ('%g' % f)
    return '%s%d_wl%swu%s' % (_g['REF_SPEC'], int(float(_g['REF_DEG'])),
                              _w(_g['REF_WLOWER']), _w(_g['REF_WUPPER']))


_REF_TAG = _ref_tag()
_od = pd.read_csv(os.path.join(OUT, 'oster_delta_full.csv')).set_index('tag')
if _REF_TAG not in _od.index:
    raise SystemExit(
        "[STOP] oster_delta_full.csv: reference label %r missing.\n"
        "  labels present: %s\n"
        "  Producer: 01c_robustness.R (section R1). The reference comes from\n"
        "  00_parameters.R (REF_SPEC/REF_DEG/REF_WLOWER/REF_WUPPER)."
        % (_REF_TAG, sorted(_od.index)))
print('[oster] reference specification: %s (00_parameters.R)' % _REF_TAG)
##############################################################################
# OSTER -- LA CELLULE LA MOINS FAVORABLE, ET D3 REECRITE (lot 63, 29/09/2026).
#
# Jusqu'au lot 62, \EMPosterDelta, \EMPosterDeltaSpec et \EMPosterMoveRel lisaient
# la REFERENCE de 00_parameters.R, et D3 arretait si la reference cessait de porter
# le plus petit delta* positif. Or le papier ne dit jamais que la reference gagne :
# §2.3 dit « the smallest positive selection ratio », le supplement « the least
# favourable entry on the grid ». D3 gardait une coincidence de grille, qui s'est
# defaite le 29/09 (reference 14,05, rang 3 ; minimum 2,17, polynomial3_wl2.5wu2.5).
#
# AUTHOR'S DECISION : les macros lisent la cellule la moins favorable, et D3
# garde ce qui porte le resultat -- que ce minimum depasse UN -- plus quatre
# coherences macros/table. Detail et tests : 02_GE/garde_oster.py,
# tools/test_garde_oster.py. La reference de 00_parameters.R reste celle des
# diagnostics de l'annexe C ; son existence dans la table est verifiee plus haut.
##############################################################################
sys.path.insert(0, HERE)
import garde_oster as _go   # noqa: E402
_OSTER = _go.emissions(_od)
print('[oster] least favourable cell: %s (delta* = %s; reference %s = %.4f)'
      % (_OSTER['tag'], _OSTER['EMPosterDelta'], _REF_TAG,
         float(_od.loc[_REF_TAG, 'delta_term1'])))
_d3_err, _d3_monde = _go.verifie(_od, _OSTER)
if _d3_err:   # coherences macros/table : BLOQUANTES, meme en reprise
    raise SystemExit("[STOP] D3 (Oster): macros and table disagree:\n  "
                     + "\n  ".join(_d3_err))
garde_phrase(_d3_monde is None,
             'D3 §2.3: the smallest positive delta* exceeds the threshold of one',
             '> 1', _d3_monde or 'ok',
             "[STOP] §2.3 (D3): %s. The sentence \"against a threshold of one\" "
             "no longer holds." % _d3_monde)
# `expected` = valeurs de l'archive (la reference etait alors le minimum) : en
# reprise, l'ecart est rapporte et la macro emise.
add('EMPosterDelta', _OSTER['EMPosterDelta'], '2.2',
    'oster_delta_full.csv (%s = smallest delta_term1 > 0)' % _OSTER['tag'])
add('EMPosterDeltaNegN', _OSTER['EMPosterDeltaNegN'], '4',
    'oster_delta_full.csv ((delta_term1 < 0).sum(), STRICTLY negative)')
add('EMPosterDeltaSpec', _OSTER['EMPosterDeltaSpec'],
    'the cubic polynomial winsorised at $2.5\\%$',
    'oster_delta_full.csv, row %s: spec/degree/wval_* put into words' % _OSTER['tag'])
# \EMPosterDeltaRef (lot 65, author's decision) : le delta* de la specification de REFERENCE
# de 00_parameters.R (REF_SPEC/REF_DEG/REF_W*), a cote de la cellule la moins
# favorable. La phrase du §2.3 dira les deux ; celle-ci lit une cellule NOMMEE,
# pas un argmin. `expected` = archive (8,1) : en reprise, l'ecart est rapporte.
add('EMPosterDeltaRef', '%.1f' % float(_od.loc[_REF_TAG, 'delta_term1']), '14.0',
    'oster_delta_full.csv (%s = reference of 00_parameters.R, delta_term1)' % _REF_TAG)
# Retenu pour la garde de concordance de \EMPgridNspec, plus bas : la phrase
# de §2.3 emploie \EMPgridNspec comme DENOMINATEUR du compte Oster.
_OSTER_NSPEC = len(_od)

# ------------------------------------------------- tests joints / p-values ---
# The former table parser is gone: the macros it fed were retired, and the
# table it read now holds non-linearity p-values for two arms rather than the
# joint p-values it looked for. No regex could satisfy it.
#
# SOURCE = the xtable outputs `*_robustness.tex`, NOT the `tab_*_formatted.tex`
# files. The four formatted tables have no producer in this package: they are
# laid out by hand. Reading them would certify a number that a human typed,
# which is circular, not provenance. The generated ones are prod_robustness.tex,
# amenity_robustness.tex, amenity_annual_robustness.tex and
# prod_annual_robustness.tex.

# Productivity, daily distribution. The two statements in the text are MAXIMA
# over the specification grid, so the maxima are what gets emitted.
# The two macros that parsed the rendered table have been retired. They read the
# daily arm only, from values xtable had already rounded, for sentences covering
# BOTH arms.
# They emit MAXIMA, not thresholds. A maximum stays true whatever threshold the
# reader applies, whereas a claim of unanimity depends on the threshold chosen:
# linearity is rejected in 24 cells out of 24 at 5 per cent, but in 19 at 1 per
# cent and 17 at 0.1 per cent. A unanimity that holds at one threshold and fails
# at another is manufactured by the choice of threshold unless it carries it.
##################################################################################
# LA SOURCE EST LA GRILLE DE QUINZE, PAS CELLE DE DOUZE (18/09/2026, lot 10).
#
# prod_joint_tests_{daily,annual}.csv sont les DOUZE specifications que 01a
# estime -- SPEC_CONFIGS, la grille du corps. Les phrases que ces macros
# alimentent parlent de « the grid » au sens du contraste, qui en compte
# QUINZE. On lit donc prod_joint_tests_grid_*.csv, emis par ct4 sur
# GRID_SPECS, et ct4 y a deja verifie valeur par valeur que les douze
# specifications communes reproduisent celles de 01a.
#
# LE COMPTE EST VERIFIE CONTRE 00_parameters.py, pas contre lui-meme. Compter
# les lignes du CSV recu rendrait la macro vraie de ce que ct4 a ecrit, jamais
# de ce que l'auteur a decide : c'est exactement le defaut ferme dans ct4 le
# 18/09 pour \EMPgridProdNspec.
##################################################################################
_gfd = os.path.join(OUT, 'prod_joint_tests_grid_daily.csv')
_gfa = os.path.join(OUT, 'prod_joint_tests_grid_annual.csv')
for _f in (_gfd, _gfa):
    if not os.path.exists(_f):
        sys.exit("[STOP] %s missing.\n"
                 "  <- 04_contrast/ct4_emit_tab_contribA_grid.R, same pass\n"
                 "     as tab_robustness_prod_joint.tex. Without it the grid\n"
                 "     macros have no source, and the emitter does NOT fall back\n"
                 "     on the twelve of 01a: that would publish a subset under\n"
                 "     the name of the grid." % os.path.relpath(_f, REPO))
_cd = pd.read_csv(_gfd)
_ca = pd.read_csv(_gfa)
for _nm, _df in (('daily', _cd), ('annual', _ca)):
    if set(_df.model_name) != set(_PARAMS.GRID_LABELS):
        sys.exit("[STOP] prod_joint_tests_grid_%s.csv does not carry GRID_SPECS.\n"
                 "  extra   : %s\n  missing : %s"
                 % (_nm,
                    ', '.join(sorted(set(_df.model_name) - set(_PARAMS.GRID_LABELS))) or '-',
                    ', '.join(sorted(set(_PARAMS.GRID_LABELS) - set(_df.model_name))) or '-'))
_alljoint  = list(_cd.p_joint)  + list(_ca.p_joint)
_allnonlin = list(_cd.p_nonlin) + list(_ca.p_nonlin)
_ncells    = len(_alljoint)
if _ncells != 2 * _PARAMS.GRID_NSPEC:
    sys.exit("[STOP] %d cells for %d specifications x 2 metrics."
             % (_ncells, _PARAMS.GRID_NSPEC))
_ceil = lambda x, d: math.ceil(x * 10**d) / 10**d      # an UPPER bound, never a rounding
add('EMPprodPjointMax',  f"{_ceil(max(_alljoint), 4):.4f}",  '0.0001',
    'prod_joint_tests_grid_{daily,annual}.csv, upper bound of the max over the two arms')
add('EMPprodPnonlinMax', f"{_ceil(max(_allnonlin), 3):.3f}", '0.066',
    'same, non-linearity (0.039 on the grid of 12: the three Poly 4 of the daily '
    'arm carry the weakest non-linearity of the grid)')
add('EMPgridNcells', str(_ncells), '30',
    'prod_joint_tests_grid_*.csv: GRID_NSPEC x 2 exposure metrics '
    '(24 on the grid of 12, until 18/09/2026)')
##################################################################################
# LA GRILLE D'AMENITE N'EST PAS CELLE DE PRODUCTIVITE, ET ELLE NE LE DEVIENT PAS.
#
# \EMPgridNcells etait appelee par DEUX passages du supplement : celui de la
# productivite (Table de robustesse, 15 specifications depuis le lot 10) et
# celui de l'amenite, qui lit amenity_joint_tests.csv -- douze specifications,
# vingt-quatre cellules, inchangees. Une seule macro pour deux grilles etait
# juste tant que les deux avaient la meme taille ; elle ne l'est plus.
# La macro d'amenite est donc EMISE SEPAREMENT. ⚠️ Les deux phrases du
# supplement qui l'appellent doivent etre repointees par l'auteur -- elles
# disent aujourd'hui \EMPgridNcells la ou \EMPamenGridNcells est vraie.
##################################################################################
# N du stage 3, LU DANS LE MEME DATA FRAME que la note de
# tab_robustness_prod_joint, qui l'imprime a l'identique. Un N recopie d'une
# note dans une macro est le geste que ce paquet s'interdit.
# Le nom ne peut pas s'ecrire « Stage3 » : LaTeX refuse le chiffre dans un nom
# de commande, et \EMPtradeNobsStage3 etait lu \EMPtradeNobsStage + « 3 ».
_nobs3 = sorted(set(_cd.nobs) | set(_ca.nobs))
## Lot 57 : RESTE BLOQUANTE. Elle verifie une phrase ET une coherence de calcul :
## toutes les cellules de la grille doivent etre estimees sur le meme echantillon.
if len(_nobs3) != 1:
    raise SystemExit(
        "[STOP] nobs not constant over the stage-3 grid: %s\n"
        "  The note of tab_robustness_prod_joint states \"N ... in every cell,\n"
        "  identical across the two metrics\". It would be false." % _nobs3)
add('EMPtradeNobsStageThree', _thousands(_nobs3[0]), '3{,}511{,}478',
    'prod_joint_tests_{daily,annual}.csv, column nobs (same object as the note '
    'of tab_robustness_prod_joint)')

##################################################################################
# L'EXCLUSION FAUTE DE CLIMAT (lot 82, 02/10/2026). Entre l'echantillon de la
# premiere etape (\EMPtradeNobs) et celui de l'estimation (\EMPtradeNobsStageThree),
# les lignes perdues ne viennent PAS des exclusions sectorielles de B.1 (deja
# faites avant la premiere etape) : ce sont les flux dont l'exportateur ou
# l'importateur n'a aucun releve climatique ERA5. C'est le filtre ecrit
# stage3_climate_filter() (lib_drf.R), que 01a applique depuis le lot 41.
# Calcule ICI sur les deux intrants de ce filtre (BACI_raw_pairs_annual.rds et
# climate_long_annual.rds, ecrits par 01a), et confronte a la difference des deux
# effectifs que les macros ci-dessus lisent dans DEUX autres artefacts : l'egalite
# est BLOQUANTE (coherence entre deux calculs).
# ⚠ 01c ecarte les memes lignes par la gestion des NA de feglm, sans le dire
#   (registre, lot 41) : « 01c tombe sur le meme effectif » reste a verifier apres
#   course.
##################################################################################
import pyreadr   # noqa: E402
_bx = pyreadr.read_r(os.path.join(INPUTS, 'BACI_raw_pairs_annual.rds'))[None]
_cx = pyreadr.read_r(os.path.join(INPUTS, 'climate_long_annual.rds'))[None]
_clim_iso = set(_cx['GID_0'].astype(str))
_pays_ech = set(_bx['orig'].astype(str)) | set(_bx['dest'].astype(str))
_sans_clim = sorted(_pays_ech - _clim_iso)
_lignes_sans = int((~(_bx['orig'].astype(str).isin(_clim_iso)
                      & _bx['dest'].astype(str).isin(_clim_iso))).sum())
_n1 = int(float(_stage1('trade')[0]['nobs']))
import gardes_lot82 as _g82   # noqa: E402
_ok_ef, _det_ef = _g82.garde_effectifs(int(len(_bx)), _n1, _lignes_sans, int(_nobs3[0]),
                                       len(_sans_clim))
if not _ok_ef:
    sys.exit("[STOP] exclusion for lack of climate: %s. The three artefacts no longer "
             "describe the same sample." % _det_ef)
print('[effectifs] %d rows, %d countries without ERA5 climate (%s) -> %d'
      % (_lignes_sans, len(_sans_clim), ' '.join(_sans_clim), int(_nobs3[0])))
_SRC_NOCLIM = ('BACI_raw_pairs_annual.rds x climate_long_annual.rds (inputs of '
               'stage3_climate_filter, 01a); blocking equality with EMPtradeNobs - '
               'EMPtradeNobsStageThree')
add('EMPtradeNoClimCountries', str(len(_sans_clim)), '21',
    _SRC_NOCLIM + ': countries of the sample without an ERA5 record')
add('EMPtradeNoClimObs', _thousands(_lignes_sans), '219{,}480',
    _SRC_NOCLIM + ': rows whose exporter or importer has no ERA5 record')
add('EMPprodNnonlinRejOnePct', str(sum(x < 0.01 for x in _allnonlin)), '22',
    'prod_joint_tests_grid_*.csv: cells rejecting linearity at the '
    '1 % level; 19/24 on the grid of 12')

##################################################################################
# « IN ALL CELLS » EST UNE AFFIRMATION D'UNANIMITE, DONC UN COMPTE.
#
# Le manuscrit portait « linearity is rejected in all \EMPgridNcells cells at
# 5 % ». C'etait vrai des 24 cellules de la grille de douze et FAUX des 30 :
# les trois polynomes de degre quatre du bras journalier sortent a p = 0,053,
# 0,066 et 0,082. Le mot « all » etait le seul litteral de la phrase, et c'est
# lui qui a menti -- exactement le mecanisme de \EMPdiffExclNaiveExceptions,
# ou « the three degree-two polynomials » decrivait six exceptions.
#
# ⇒ Le compte de rejets, la LISTE des exceptions et leurs bornes de p sont
# desormais trois macros. La liste est CONSTRUITE par ct4 avec le constructeur
# qui refuse de sur-generaliser : si les trois Poly 4 ne sont plus exceptees
# ensemble, il rend l'enumeration au lieu de la famille (mutation verifiee le
# 18/09/2026 : deux Poly 4 sur trois rendent « X7 and X8 », pas « the three »).
##################################################################################
_pjcp = os.path.join(OUT, 'prod_joint_tests_grid_counts.csv')
if not os.path.exists(_pjcp):
    sys.exit("[STOP] prod_joint_tests_grid_counts.csv missing from output/\n"
             "  <- 04_contrast/ct4_emit_tab_contribA_grid.R, same pass as\n"
             "     tab_robustness_prod_joint.tex.")
_pjc = pd.read_csv(_pjcp).set_index('macro')['value']
SRC_PJC = ('prod_joint_tests_grid_counts.csv <- ct4_emit_tab_contribA_grid.R, '
           'same pass as the note of tab_robustness_prod_joint')


def _pj(name):
    if name not in _pjc.index:
        sys.exit("[STOP] %s missing from prod_joint_tests_grid_counts.csv "
                 "(present: %s)." % (name, ', '.join(map(str, _pjc.index))))
    return str(_pjc.loc[name])


## LE COMPTE EST RECALCULE ICI ET COMPARE A CELUI DE ct4. Les deux lisent la
## meme source, donc ils doivent coincider ; s'ils divergent, c'est que la
## note de la table et la phrase du corps ne decrivent plus la meme grille --
## le defaut meme que ce bloc existe pour fermer.
_n05_local = sum(x < 0.05 for x in _allnonlin)
if _pj('EMPprodNnonlinRejFivePct') != str(_n05_local):
    sys.exit("[STOP] rejections at 5 %%: ct4 reports %s, this pass counts %d.\n"
             "  The note of the table and the sentence of the body would describe\n"
             "  two different grids." % (_pj('EMPprodNnonlinRejFivePct'),
                                         _n05_local))
add('EMPprodNnonlinRejFivePct', _pj('EMPprodNnonlinRejFivePct'), '27', SRC_PJC)
add('EMPprodNnonlinExceptions', _pj('EMPprodNnonlinExceptions'),
    'the three degree-four polynomials under the daily distribution',
    SRC_PJC, regen=True)
add('EMPprodPnonlinExcLo', _pj('EMPprodPnonlinExcLo'), '0.050', SRC_PJC)
add('EMPprodPnonlinExcHi', _pj('EMPprodPnonlinExcHi'), '0.065', SRC_PJC)

# Amenity, ANNUAL arm. The joint p-values are emitted further down, in the same
# block as the F statistics and from the SAME artefact. Taking two columns of one
# test from two different sources is what allows them to drift apart.

# Amenity, DAILY arm: both the F statistics and the p-values come from the
# producer CSV, which is preferred to the rendered table because it carries the
# F, which the table does not.
_jt = pd.read_csv(os.path.join(OUT, 'amenity_joint_tests.csv'))
if len(_jt) != 12:
    sys.exit(f"[amenite] amenity_joint_tests.csv: {len(_jt)} rows, 12 expected.")
# Le compte de cellules du canal AMENITE : douze specifications x deux metriques.
# Il vaut 24, et il valait deja 24 quand \EMPgridNcells valait 24 -- ce qui est
# precisement pourquoi la confusion avait pu tenir sans se voir.
add('EMPamenGridNcells', str(2 * len(_jt)), '24',
    'amenity_joint_tests.csv: 12 specifications x 2 metrics; the amenity '
    'channel stays on SPEC_CONFIGS, not on GRID_SPECS')
add('EMPamenFjointLo', f"{_jt.fval_joint.min():.2f}", '1.72', 'amenity_joint_tests.csv (min joint F, daily)')
add('EMPamenFjointHi', f"{_jt.fval_joint.max():.2f}", '2.60', 'same, max')
# THREE decimals. The table note prints %.3f; at %.2f the same joint p would
# read "0.20" in the body and "0.196" in the note of the same paragraph, and the
# reader would not know which one governs.
add('EMPamenPjointLo', f"{_jt.pval_joint.min():.3f}", '0.074', 'same, min joint p')
add('EMPamenPjointHi', f"{_jt.pval_joint.max():.3f}", '0.161', 'same, max')

# OPEN DEBT, deliberately not turned into a macro. The joint F of the ANNUAL arm
# has no artefact: the persistence written by the estimation step covers the
# daily arm only, the annual arm going through a code path that exposes no object
# at the top level. It therefore remains a BARE LITERAL in the manuscript, and it
# is deliberately absent from LITERAL_KEEP, so that the literal guardrail keeps
# reporting it as a residue. Closing it means writing the annual object at its
# own site and re-running the estimation. Do not invent a source for it.

# ---------------------------------------------------------------- CRH --------
# Source = the inference artefact, PUBLISHED variant, which is the country
# clustering, the more conservative of the two. The LaTeX table is no longer
# parsed: a table parser is fragile, and it broke as soon as the table gained
# interval columns. It failed loudly rather than emitting a wrong number, and it
# was replaced rather than patched.
#
# The path below points inside the package. An earlier version read the artefact
# from a working directory excluded by .gitignore, so a clean clone produced a
# tree where this emitter died on FileNotFoundError: the manuscript compiled and
# the chain that makes its numbers did not.
# It lives in output/ because it is a computed output, not a chosen input.
CRH_CSV = os.path.join(REPO, 'output', 'crh_inference', 'crh_inference_deliverable.csv')
_crh_all = pd.read_csv(CRH_CSV)



# --- choice of standard-error variant ----------------------------------------
# `country` is the PUBLISHED variant, wider on every interval; `conley` is the
# original CRH structure, cited only where the manuscript reports what the
# pre-committed rule costs. The variant is an EXPLICIT argument and never a
# silent default: reading the wrong one breaks nothing, it merely returns
# narrower intervals, which is exactly the kind of error that has to be made
# noisy.
_PUB = {'country': 'PUBLISH', 'conley': 'reference'}

def _crh(channel, unit, obj, variant='country'):
    """Renvoie (point, ci_lo, ci_hi, se) depuis le livrable d'inference."""
    r = _crh_all[(_crh_all.channel == channel) & (_crh_all.unit == unit)
                 & (_crh_all.object == obj) & (_crh_all.vcov_variant == variant)]
    if len(r) != 1:
        sys.exit(f"[CRH] {channel}/{unit}/{obj}/{variant}: {len(r)} row(s), 1 expected "
                 f"— {CRH_CSV} has changed structure.")
    r = r.iloc[0]
    if not str(r.publish).startswith(_PUB[variant]):
        sys.exit(f"[CRH] {channel}/{unit}/{obj}/{variant}: column publish = "
                 f"'{r.publish}', expected a prefix '{_PUB[variant]}' — the "
                 f"mapping variant<->publication status has moved.")
    return (float(r.point), float(r.ci95_lo), float(r.ci95_hi),
            float(r.se) if pd.notna(r.se) else None)

def _sg(x, nd=3):
    """Signe explicite. Les estimations en points de log sont des quantites
    SIGNEES : le signe appartient a la valeur, pas a la phrase. Un '+' code en
    dur dans le .tex devient un faux le jour ou la valeur change de signe —
    c'est precisement ce qui est arrive a la France (contrat §5)."""
    return f"{x:+.{nd}f}"

def _ab(x, nd=3):
    return f"{x:.{nd}f}"

# --------------------------------------------------------------------------
# REBASAGE DU 10/09/2026 -- les 14 `expected` de la chaine CRH
# --------------------------------------------------------------------------
# Les quatorze macros marquees ci-dessous sortaient en MISMATCH au TROISIEME
# chiffre. Cause : la chaine CRH a ete rejouee sur les intrants d'exposition
# re-emis. Les quatre artefacts producteurs datent tous du 10/09
# (crh_inference_deliverable, crh_inference_raw, crh_bootstrap_raw,
# crh_period_sensitivity), et le deplacement d'entree s'est propage a travers
# la chaine sans s'amortir ni s'amplifier.
#
# LE DEPLACEMENT EST ABSOLU, PAS RELATIF. Valeur brute mesuree, et ecart au
# litteral ancien, en points de log :
#   NorthMean        +0.558741 (+0.0027)   NorthMeanHi      +0.929131 (+0.0061)
#   NorthDaily       +0.116158 (+0.0022)   NorthDailyLo     -0.184076 (+0.0029)
#   NorthDailyHi     +0.414218 (+0.0022)   NorthDiff        +0.442583 (+0.0016)
#   NorthDiffLo      +0.234744 (+0.0017)   FraDailyLo       -0.450068 (+0.0009)
#   NorthDiffBoot    +0.442583 (+0.0016)   NorthDiffBootHi  +0.608393 (+0.0014)
#   NorthDiffDecTwo  +0.404997 (+0.0020)   NorthDiffSEDecTwo 0.098887 (+0.0009)
#   AmenNorthDiffHi  +0.059636 (+0.0006)   AmenDiffConleyHi +0.051471 (-0.0005)
# Les quatorze tiennent dans [-0.0006, +0.0062], mediane +0.0017.
#
# ⚠️ NE PAS RESUMER CET ECART PAR UN POURCENTAGE UNIQUE. Les memes quatorze
# s'etalent, en relatif, de -1.6 % a +1.9 % : un decalage absolu quasi constant
# divise par des niveaux allant de 0.05 a 0.93 ne rend pas un pourcentage
# constant. Le « +0.30 % sur mean - daily » ne decrit QUE la difference
# (+0.359 % mesure) ; le porter sur les niveaux serait un faux -- le niveau
# journalier bouge de +1.9 %, six fois plus.
#
# CE QUI NE BOUGE PAS, et c'est ce qui autorise le rebasage : le signe des
# quatorze ; l'ordre mean > daily ; l'exclusion de zero par l'intervalle de la
# difference (+0.235 ; +0.647) ; l'inclusion de zero par celui du niveau
# journalier (-0.184 ; +0.414) ; et la comptabilite, qui se referme au
# millionieme -- 0.558741 - 0.116158 = 0.442583. Aucun enonce du manuscrit ne
# repose sur le troisieme chiffre de ces quantites.
#
# --- productivity, Northern aggregate above 45N, published variant -----------
_p, _lo, _hi, _se = _crh('prod', 'NORTH>=45N', 'level_mean_2095')
add('EMPcrhNorthMean',   _sg(_p),  '+0.559', 'crh_inference_deliverable.csv (prod/NORTH/mean, country)')
## LOT 85o : bornes MC re-posees sur la course Cholesky (step4_mc ; rapport_lot85n §2d).
add('EMPcrhNorthMeanLo', _sg(_lo), '+0.186', 'same, ci95_lo')
add('EMPcrhNorthMeanHi', _sg(_hi), '+0.943', 'same, ci95_hi')

_p, _lo, _hi, _se = _crh('prod', 'NORTH>=45N', 'level_daily_2095')
add('EMPcrhNorthDaily',   _sg(_p),  '+0.116', 'crh_inference_deliverable.csv (prod/NORTH/daily, country)')
add('EMPcrhNorthDailyLo', _sg(_lo), '-0.184', 'same, ci95_lo')
add('EMPcrhNorthDailyHi', _sg(_hi), '+0.414', 'same, ci95_hi')

_p, _lo, _hi, _se = _crh('prod', 'NORTH>=45N', 'difference_mean_minus_daily')
add('EMPcrhNorthDiff',   _sg(_p),  '+0.443', 'crh_inference_deliverable.csv (prod/NORTH/difference, country)')
add('EMPcrhNorthDiffLo', _sg(_lo), '+0.238', 'same, ci95_lo')
add('EMPcrhNorthDiffHi', _sg(_hi), '+0.651', 'same, ci95_hi')

# ---- LES COMPTES « six cases out of six » / « one case out of six » (lot 81) ----
# §2.5 les ecrivait en litteraux. Unites : les six de la table du supplement E
# (cinq pays et l'agregat Nord), canal productivite, intervalles groupes par pays
# (variante publiee). « Significativement positif » = ci95_lo > 0. Emis en toutes
# lettres : le rendu du papier ne change pas. Garde G10 : tools/gardes_lot81.py.
import gardes_lot81 as _g81   # noqa: E402
_cu = _crh_all[(_crh_all.channel == 'prod') & (_crh_all.vcov_variant == 'country')
               & (_crh_all.object.isin(['level_mean_2095', 'level_daily_2095']))]
if (_cu.groupby('object').size() != _cu['unit'].nunique()).any() or _cu['unit'].nunique() == 0 \
        or not all(str(p).startswith(_PUB['country']) for p in _cu['publish']):
    sys.exit('[CRH] counts of §2.5: incomplete deliverable or unexpected publication status.')
_cn = int(_cu['unit'].nunique())
_ca = int((_cu.loc[_cu.object == 'level_mean_2095', 'ci95_lo'] > 0).sum())
_cd = int((_cu.loc[_cu.object == 'level_daily_2095', 'ci95_lo'] > 0).sum())
_ok, _det = _g81.garde_g10(_ca, _cd, _cn)
print('[CRH] G10 %s: %s' % ('ok' if _ok else 'REFUSED', _det))
garde_phrase(_ok, '§2.5 CRH counts six/six, one/six', 'true', _det,
             "[STOP] §2.5 (G10): %s. The sentence \"significantly positive in ... cases out of ...; "
             "... in ... case out of ...\" is false." % _det)
_SRC_CU = ('crh_inference_deliverable.csv (prod, country, %d units): count of ci95_lo > 0'
           % _cn)
add('EMPcrhSigPosAnnual', _g81.MOTS[_ca], 'six', _SRC_CU + ', level_mean_2095, in words')
add('EMPcrhSigPosDaily',  _g81.MOTS[_cd], 'one', _SRC_CU + ', level_daily_2095, in words')
add('EMPcrhNunits',       _g81.MOTS[_cn], 'six', 'crh_inference_deliverable.csv (prod, country): '
    'units of the table of supplement E, in words')

# ---- L'ATTENUATION RELATIVE : la seule quantite comparable entre pipelines --
# Posee le 20/09/2026. Les NIVEAUX des deux pipelines ne se comparent pas --
# horizon, scenario d'emissions, et projection en equilibre partiel sans leur
# innovation endogene. La PART DU GAIN ANNUEL QUE LE JOURNALIER RETIRE, elle,
# est sans dimension et elle se compare.
#
# POURQUOI CE RAPPORT ET PAS SON INVERSE. Le rapport journalier/annuel a pour
# numerateur un objet que le papier declare NON SIGNE : aucune des 15
# specifications n'exclut zero, et le journalier nordique de CRH vaut +0,116
# avec un IC [-0,184 ; +0,414]. Un rapport dont le numerateur enjambe zero n'a
# pas d'intervalle exploitable. Son COMPLEMENT -- difference / annuel -- a pour
# numerateur la DIFFERENCE, le mieux identifie des trois objets (exclut zero en
# 15/15 chez nous, [+0,235 ; +0,647] chez eux). Meme information, meme
# invariance a l'echelle, deux termes signes.
#
# ⚠️ AU-DELA DE 100 % le bras journalier est NEGATIF : la metrique ne retire pas
# seulement le gain, elle le retourne. Le manuscrit le dit en clause, sans quoi
# le mot « share » contredirait le chiffre qu'il introduit.
_p, _lo, _hi, _se = _crh('prod', 'NORTH>=45N', 'attenuation_ratio')
add('EMPcrhAtten',   _ab(100 * _p, 1),  '79.2',
    'crh_inference_deliverable.csv (prod/NORTH/attenuation_ratio, country): PUBLISHED variant, vcov clustered by country; the boot_country row of the same file is marked CONTROL ONLY and is not this one')
add('EMPcrhAttenLo', _ab(100 * _lo, 1), '47.8', 'same, ci95_lo')
add('EMPcrhAttenHi', _ab(100 * _hi, 1), '181.6', 'same, ci95_hi')

# Notre cote, dans NOTRE idiome : l'ETENDUE SUR LA GRILLE, comme le fait la
# Section 2.5 partout. Les deux cotes ne sont donc PAS apparies, et c'est
# voulu : apparier deux constructions d'intervalle differentes serait comparer
# deux CONSTRUCTIONS en croyant comparer deux RESULTATS -- la faute des
# conventions d'echelle, transposee aux methodes d'inference. Chaque cote est
# rendu dans l'idiome de sa propre section, et aucun appariement n'est
# revendique.
_atn = pd.read_csv(os.path.join(OUT, 'contribA', 'tache1_diff_bras.csv'))
_atn = _atn[_atn.channel == 'prod']
if len(_atn) == 0:
    sys.exit('[STOP] tache1_diff_bras.csv: no row of the prod channel.')
_att = 100.0 * _atn['diff'] / _atn['delta_annual']
add('EMPattenNorthMin', _ab(_att.min(), 1), '83.6',
    'tache1_diff_bras.csv (prod): min of diff/delta_annual over the grid')
add('EMPattenNorthMax', _ab(_att.max(), 1), '104.2', 'same, max')

# France is the only negative daily point of the productivity panel. Its sign is
# NOT identified, so the macro carries the sign explicitly.
_p, _lo, _hi, _se = _crh('prod', 'FRA', 'level_daily_2095')
add('EMPcrhFraDaily',   _sg(_p),  '-0.187', 'crh_inference_deliverable.csv (prod/FRA/daily, country)')
add('EMPcrhFraDailyLo', _sg(_lo), '-0.453', 'same, ci95_lo')
add('EMPcrhFraDailyHi', _sg(_hi), '+0.071', 'same, ci95_hi')

# Correlation between the two arms. The argument is that drawing independently
# from two separate covariance matrices would give a wrong variance for the
# difference. Reported across the six units of the published variant.
_ca = pd.read_csv(os.path.join(REPO, 'output', 'crh_inference', 'crh_inference_raw.csv'))
_ca = _ca[(_ca.channel == 'prod') & (_ca.vcov == 'country')].cor_arms.dropna()
add('EMPcrhCorArmsLo', _ab(_ca.min(), 2), '0.67', 'crh_inference_raw.csv (cor_arms, prod/country, min)')
add('EMPcrhCorArmsHi', _ab(_ca.max(), 2), '0.89', 'same, max')

# --- bootstrap control, independent of the analytic route --------------------
_bt = pd.read_csv(os.path.join(REPO, 'output', 'crh_inference', 'crh_bootstrap_raw.csv'))
_bt = _bt[(_bt.channel == 'prod') & (_bt.unit == 'NORTH>=45N')]
if len(_bt) != 1:
    sys.exit(f"[CRH] bootstrap prod/NORTH: {len(_bt)} row(s), 1 expected.")
_bt = _bt.iloc[0]
add('EMPcrhNorthDiffBoot',   _sg(_bt.diff_point), '+0.443', 'crh_bootstrap_raw.csv (prod/NORTH)')
add('EMPcrhNorthDiffBootLo', _sg(_bt.diff_lo),    '+0.118', 'same, diff_lo')
add('EMPcrhNorthDiffBootHi', _sg(_bt.diff_hi),    '+0.608', 'same, diff_hi')

# --- period sensitivity -------------------------------------------------------
# P_all is primary and P_dec2 is the robustness arm. P_no90 must not be cited,
# because its windows overlap, so it deliberately has no macro.
_ps = pd.read_csv(os.path.join(REPO, 'output', 'crh_inference', 'crh_period_sensitivity.csv'))
def _per(channel, unit, spec_prefix, col):
    r = _ps[(_ps.channel == channel) & (_ps.unit == unit)
            & (_ps.spec.str.startswith(spec_prefix))]
    if len(r) != 1:
        sys.exit(f"[CRH] periods {channel}/{unit}/{spec_prefix}: {len(r)} row(s), 1 expected.")
    return float(r[col].iloc[0])

add('EMPcrhNorthDiffDecTwo', _sg(_per('prod', 'NORTH>=45N', 'P_dec2', 'diff')), '+0.405',
    'crh_period_sensitivity.csv (prod/NORTH, P_dec2)')
add('EMPcrhNorthDiffSE',       _ab(_per('prod', 'NORTH>=45N', 'P_all',  'se_diff')), '0.105',
    'crh_period_sensitivity.csv (prod/NORTH, P_all, se_diff)')
add('EMPcrhNorthDiffSEDecTwo', _ab(_per('prod', 'NORTH>=45N', 'P_dec2', 'se_diff')), '0.099',
    'same, P_dec2')

# --- amenity channel: inconclusive, and reported as such ---------------------
_p, _lo, _hi, _se = _crh('amen', 'NORTH>=45N', 'difference_mean_minus_daily')
add('EMPcrhAmenNorthDiff',   _sg(_p),  '+0.028', 'crh_inference_deliverable.csv (amen/NORTH/difference, country)')
add('EMPcrhAmenNorthDiffLo', _sg(_lo), '-0.003', 'same, ci95_lo')
add('EMPcrhAmenNorthDiffHi', _sg(_hi), '+0.060', 'same, ci95_hi')

# The conley variant serves the one sentence that states what the pre-committed
# rule costs. It is the only site in the manuscript allowed to cite it.
_p, _lo, _hi, _se = _crh('amen', 'NORTH>=45N', 'difference_mean_minus_daily', variant='conley')
add('EMPcrhAmenNorthDiffConleyLo', _sg(_lo), '+0.005', 'crh_inference_deliverable.csv (amen/NORTH/difference, conley)')
add('EMPcrhAmenNorthDiffConleyHi', _sg(_hi), '+0.051', 'same, ci95_hi')

add('EMPcrhAmenFraDaily',       _sg(_per('amen', 'FRA', 'P_all',  'daily')), '-0.004',
    'crh_period_sensitivity.csv (amen/FRA, P_all, daily)')
add('EMPcrhAmenFraDailyDecTwo', _sg(_per('amen', 'FRA', 'P_dec2', 'daily')), '+0.002',
    'same, P_dec2')

# --- first-moment coherence, the guardrail of the contrast -------------------
_mc = pd.read_csv(os.path.join(OUT, 'crh_metric_coherence.csv')).iloc[0]
add('EMPcrhMetricCor',   _ab(_mc.cor, 3),           '0.996', 'crh_metric_coherence.csv (cor)')
add('EMPcrhMetricSlope', _ab(_mc.lm_slope, 2),      '0.97',  'crh_metric_coherence.csv (lm_slope)')
add('EMPcrhMetricBias',  _sg(_mc.bias_era_crh, 2),  '-0.03', 'crh_metric_coherence.csv (bias_era_crh)')

# --- bliss points, re-estimated on the CRH design -----------------------------
# The values 1.6 and 8.6 are Cruz and Rossi-Hansberg own published figures, so
# they stay LITERALS with a citation. Turning them into macros would claim that
# this package recomputes them. Only our values are emitted here.
# REBASAGE DU 10/09/2026, meme rejeu CRH que ci-dessus. Bruts mesures :
#   prod 1.90691676 (etait 1.8)   amen 8.55440104 (etait 8.9)
# Les deux points de bonheur bougent en sens OPPOSES, ce qui exclut un decalage
# commun de la courbe et designe une re-estimation, pas un offset.
#
# ⚠️ CONSEQUENCE SUR L'ENONCE, ET ELLE EST FAVORABLE D'UN COTE, DEFAVORABLE DE
# L'AUTRE. La phrase (Section CRH) dit « reproduces their published bliss points
# on both channels », puis « I obtain X pour la productivite contre leur 1.6, et
# Y pour les amenites contre leur 8.6 ».
#   - amenites : 8.9 -> 8.6, soit EXACTEMENT le chiffre publie par Cruz et
#     Rossi-Hansberg. L'ecart de portage disparait.
#   - productivite : 1.8 -> 1.9, l'ecart a leur 1.6 passe de 0.2 a 0.3 C.
# Aucun des deux ne renverse la phrase : elle n'a jamais chiffre l'ecart, elle
# affirme la reproduction. Mais le second sens est le sens qui degrade, et il
# est ecrit ici pour qu'il ne se lise pas seulement dans le premier.
for _n, _f, _exp in [('EMPcrhOptProd', 'crh_design_optimum_prod.csv',    '1.9'),
                     ('EMPcrhOptAmen', 'crh_design_optimum_amen.csv',    '8.6')]:
    add(_n, _ab(pd.read_csv(os.path.join(OUT, _f)).T_opt.iloc[0], 1), _exp, f'{_f} (T_opt)')

# --- peaks of the smoothed delta curve, and the support edge -----------------
# These are the peaks of the SMOOTHED delta(T) in the CRH design. They are
# neither the bliss points above nor the argmax of the fitted response reported
# in the paper. The CSV names its object explicitly so that no one has to guess
# which notion of optimum a number belongs to.
_pk = pd.read_csv(os.path.join(OUT, 'crh_daily_response_peaks.csv')).set_index('channel')
add('EMPcrhPeakProd', _ab(_pk.loc['prod', 'peak_T'], 1), '12.2', 'crh_daily_response_peaks.csv (peak_T, prod)')
add('EMPcrhPeakAmen', _ab(_pk.loc['amen', 'peak_T'], 1), '13.5', 'same, amen')

##################################################################################
# LES TROIS CHIFFRES QUE LA LEGENDE DE fig:crh_response LIT (18/09/2026, lot 10b)
#
# Ils sortent de crh_mean_response_optima.csv, ecrit par 01d_fig_crh_responses.R
# DANS LE MEME PASSAGE que la figure et depuis les memes fichiers : la legende ne
# peut donc pas citer un zero que la figure ne trace pas.
#
# ⚠️ \EMPcrhMeanOpt N'EST PAS \EMPcrhOptProd, et la confusion des deux est le
# defaut que ce lot ferme.
#   \EMPcrhOptProd = 1,9 degC  : racine de la logistique de run_crh_baseline
#     (01d, C2), la specification d'ORIGINE de Cruz et Rossi-Hansberg, dont le
#     regresseur est une moyenne decennale de JANVIER au nord / JUILLET au sud.
#     Leur papier l'appelle « optimal winter temperature » (p. 20). Cite dans la
#     PROSE, pour dire que la replication reproduit leurs points de felicite.
#   \EMPcrhMeanOpt = 20,0 degC : zero de la pente du bras MOYENNE ANNUELLE
#     (01d, C4(b)), celui que le manuscrit compare au journalier et que la
#     projection C6 consomme. Cite dans la LEGENDE, parce que c'est ce que la
#     colonne de droite trace.
# Les deux sont justes, chacun de son objet ; un seul nom pour les deux serait
# faux de l'un des deux.
#
# \EMPcrhDailyRefBin ferme un litteral qui serait apparu dans la legende : le
# bin de reference du bras journalier, dont le coefficient est nul par
# construction et par rapport auquel delta(T) se lit.
##################################################################################
_mop = os.path.join(OUT, 'crh_mean_response_optima.csv')
if not os.path.exists(_mop):
    sys.exit("[STOP] crh_mean_response_optima.csv missing: "
             "01_estimation/01d_fig_crh_responses.R has not run.")
_mo = pd.read_csv(_mop).set_index('channel')
add('EMPcrhMeanOpt', _ab(_mo.loc['prod', 'T_opt'], 1), '20.0',
    'crh_mean_response_optima.csv (01d_fig_crh_responses.R, same pass as the '
    'figure): zero of the slope of the annual-mean arm, productivity channel. '
    'NOT to be confused with EMPcrhOptProd, the bliss point of the original CRH '
    'specification, on a winter-month temperature')
add('EMPcrhMeanOptAmen', _ab(_mo.loc['amen', 'T_opt'], 1), '21.7',
    'same, amenity channel')
add('EMPcrhDailyRefBin', _ab(_mo.loc['prod', 'daily_ref_bin_T'], 1), '25.9',
    'same: centre of the reference bin of the daily arm, whose coefficient '
    'is zero by construction. delta(T) is read relative to it')
# The support edge is an OUTPUT of the rule "hottest bin carrying at least one
# per cent of the day mass", not a chosen threshold. The one per cent is the
# input, and it stays a literal.
_sb = pd.read_csv(os.path.join(PKG, 'intermediate_outputs', 'crh', 'crh_support_bounds.csv'))
add('EMPcrhSupportHi', _ab(_sb.k_max.iloc[0], 1), '40.0',
    'crh_support_bounds.csv (k_max, pre-registered 1 % rule)')

# --- CRH grid counts ----------------------------------------------------------
# LITERALS. The count can be checked against the per-cell projection file, but
# that file is not versioned, so wiring it here would make this emitter die on a
# clean clone. Closing the provenance means writing the count into an artefact
# under output/.
add('EMPcrhGridCells',   '17{,}048', '17{,}048', 'LITERAL — cells of the 1-degree CRH grid; '
    'to be emitted from crh_daily_projection_cells.csv once deposited', literal=True)
add('EMPcrhCellPeriods', '68{,}192', '68{,}192', 'LITERAL — 17,048 x 4 decadal cross-sections = '
    'cell-periods; the manuscript called them "cells", a mislabel', literal=True)

##################################################################################
# L'EFFET DU CLAMP SUR LA QUEUE CHAUDE — LITTERAL FERME LE 18/09/2026
#
# \EMPcrhClampMove etait un LITTERAL, et son commentaire disait pourquoi :
# « la paire de runs 24/06 a ete regeneree, re-derivation par run requise ».
# La re-derivation n'avait pas besoin d'un run. Le clamp n'entre pas dans
# l'estimation, seulement dans l'EVALUATION de delta(T) sur les distributions
# CMIP6 : 01d_crh_clamp_sensitivity.R le refait sur les archives, apres avoir
# prouve qu'il reproduit le bras non clampe de 01d a 1e-16.
#
# 🚨 LA VALEUR MESUREE N'EST PAS LE LITTERAL : 3,8e-8 contre 3,9e-6 ecrits,
# soit un facteur cent. regen=True, donc la SOURCE fait foi et l'ecart est
# REPORTE plutot que bloquant -- c'est la regle (b) du runbook, et c'est le
# cas qu'elle vise : un chiffre dont la provenance etait declaree ouverte.
# L'enonce du supplement ne change pas de sens, il se renforce : le
# deplacement est cent fois plus petit que ce qui etait publie.
#
# MAGNITUDE POUR L'UNE, SIGNE POUR L'AUTRE, et ce n'est pas une inadvertance :
# voir l'en-tete de 01d_crh_clamp_sensitivity.R. La phrase du supplement dit
# « moves the Northern aggregate BY x log points » ; celle du clamp severe
# porte le deplacement de la DIFFERENCE, dont le signe est tout le contenu.
#
# LES DEUX SEUILS VIENNENT DE 00_parameters, ET LA CONCORDANCE R/PYTHON EST
# VERIFIEE ICI. Le CSV est ecrit par le cote R ; on y CHERCHE la valeur que le
# cote Python declare. Si les deux jumeaux divergeaient, la ligne serait
# introuvable et l'emetteur s'arreterait -- ce qui est le seul controle croise
# possible sans faire tourner R depuis Python.
##################################################################################
_CL_SUP = float(_PARAMS.CRH_CLAMP_SUPPORT)
_CL_SEV = float(_PARAMS.CRH_CLAMP_SEVERE)
_cs = pd.read_csv(os.path.join(OUT, 'crh_clamp_sensitivity.csv'))
_dm = pd.read_csv(os.path.join(OUT, 'crh_north_daymass.csv'))


def _sci(x, signed=False):
    """Notation scientifique LaTeX, au format deja publie : 3.9\times10^{-6}."""
    if x == 0:
        return '0'
    e = int(math.floor(math.log10(abs(x))))
    m = round(x / 10.0 ** e, 1)
    if abs(m) >= 10:                      # l arrondi peut franchir la decade
        m, e = m / 10.0, e + 1
    head = ('%+.1f' % m) if signed else ('%.1f' % abs(m))
    return head + (r'\times10^{%d}' % e)


def _clamp_row(ch, hi):
    r = _cs[(_cs.channel == ch) & (_cs.clamp_hi == hi)]
    if len(r) != 1:
        sys.exit("[STOP] crh_clamp_sensitivity.csv: %d row(s) for channel "
                 "%s at threshold %g.\n"
                 "  The threshold is declared by 00_parameters.py; the CSV is "
                 "written under 00_parameters.R.\n"
                 "  A missing row is a divergence between the two twins, "
                 "not a reading problem." % (len(r), ch, hi))
    return r.iloc[0]


_c40 = _clamp_row('prod', _CL_SUP)
_c30 = _clamp_row('prod', _CL_SEV)

add('EMPcrhClampMove', _sci(abs(_c40.move_daily)), r'3.9\times10^{-6}',
    'crh_clamp_sensitivity.csv (01d_crh_clamp_sensitivity.R): MAGNITUDE of the '
    'move of the daily Northern aggregate when delta(T) is held flat '
    'beyond %g degC' % _CL_SUP, regen=True)
add('EMPcrhClampThirtyMove', _sci(_c30.move_diff, signed=True),
    r'+8.2\times10^{-4}',
    'same, clamp at %g degC: SIGNED move of the Northern DIFFERENCE '
    '(decadal mean minus daily). Positive = cutting the hot tail '
    'SEPARATES the two metrics instead of bringing them together' % _CL_SEV)

_d30 = _dm[_dm.threshold_C == _CL_SEV]
if len(_d30) != 1:
    sys.exit("[STOP] crh_north_daymass.csv: %d row(s) at threshold %g."
             % (len(_d30), _CL_SEV))
add('EMPcrhNorthDayMassAboveThirty', '%.2f' % float(_d30.share_pct.iloc[0]),
    '0.28',
    'crh_north_daymass.csv: share of the day mass of cells >=45N '
    'above %g degC at the horizon, population-weighted' % _CL_SEV)
## Lot 82 : le « fewer than 0.1% » du supplement E etait un litteral de resultat.
## Meme fichier, meme objet, au seuil du support (CRH_CLAMP_SUPPORT) : la part de la
## masse de jours projetee des cellules >= 45N au-dela du dernier bin tenu.
_d40 = _dm[_dm.threshold_C == _CL_SUP]
if len(_d40) != 1:
    sys.exit("[STOP] crh_north_daymass.csv: %d row(s) at threshold %g."
             % (len(_d40), _CL_SUP))
add('EMPcrhNorthDayMassAboveSupport', _sci(float(_d40.share_pct.iloc[0])),
    r'1.1\times10^{-5}',
    'crh_north_daymass.csv: share (%%) of the day mass of cells >=45N '
    'beyond %g degC (support) at the horizon, population-weighted' % _CL_SUP)
## Lot 82 : « p<10^{-16} » de C.1 etait un litteral. La source que le texte decrit
## (le test de Wald des controles de gravite, premiere etape migration) ecrit la valeur.
_m1 = _stage1('migration')[0]
if _m1 is None or 'wald_p' not in _m1.index:
    sys.exit("[STOP] stage1_migration_diagnostics.csv: key wald_p missing.")
add('EMPmigWaldP', _sci(float(_m1['wald_p'])), r'6.8\times10^{-56}',
    'stage1_migration_diagnostics.csv (01b), wald_p of the joint test of the controls')

# ---------------------------------------------------------------- eta --------
# ---------------------------------------------------------------- eta --------
# Read from output/, where the producer writes. A duplicate of this file used to
# live in a directory git does not track; an identical name invites the
# assumption that two copies are identical, and nothing guarantees it tomorrow.
d = pd.read_csv(os.path.join(REPO, 'output', 'table_eta.csv'))
amen = d[(d.metric == 'Damage index (daily bins)') & (d.channel == 'Amenity')]
prod = d[(d.metric == 'Damage index (daily bins)') & (d.channel == 'Productivity')]
add('EMPetaDailyAmen', f"{amen.eta_pct.iloc[0]:.1f}", '45.6', 'table_eta.csv')
# REBASAGE DU 10/09/2026 : 18.2 -> 20.0 (brut 0.20005021). Voir le bloc annuel
# ci-dessous : les deux macros PRODUCTIVITE bougent, les deux macros AMENITE ne
# bougent pas, ce qui localise la cause.
add('EMPetaDailyProd', f"{prod.eta_pct.iloc[0]:.1f}", '14.4', 'table_eta.csv')

# The ANNUAL counterparts are macros as well, not just the daily ones. The same
# sentence reading the same CSV under two regimes is what manufactures a drift:
# regenerating the CSV would have moved one pair and left the other in place.
# The annual-versus-daily contrast IS the argument of the subsection, so the four
# numbers must move together or not at all.
amen_a = d[(d.metric == 'Damage index (annual mean)') & (d.channel == 'Amenity')]
prod_a = d[(d.metric == 'Damage index (annual mean)') & (d.channel == 'Productivity')]
add('EMPetaAnnualAmen', f"{amen_a.eta_pct.iloc[0]:.1f}", '39.3', 'table_eta.csv')
# REBASAGE DU 10/09/2026 : 8.1 -> 9.5 (brut 0.09483841).
#
# 🔑 LA CAUSE EST LOCALISEE PAR LES DEUX MACROS QUI N'ONT PAS BOUGE. Au meme
# rejeu, EMPetaDailyAmen (31.0) et EMPetaAnnualAmen (29.3) sont sorties en
# MATCH. Le deplacement est donc confine au canal PRODUCTIVITE, c'est-a-dire au
# canal que 01a a re-emis au correctif GLORIA ; il n'y a pas de derive generale
# de table_eta.csv. Le commentaire ci-dessus exigeait que les quatre bougent
# « ensemble ou pas du tout » : la condition est respectee au niveau du CANAL,
# qui est l'unite que le rejeu a touchee.
#
# ⚠️ PHRASES RELUES, LES DEUX TIENNENT, ET AUCUNE N'A ETE REECRITE :
#   Section 3 : « blind to about a fifth ... the remaining four fifths ».
#     18.2 % arrondissait deja au cinquieme ; 20.0 % l'est exactement. La
#     fraction ecrite en toutes lettres devient juste au lieu d'approchee.
#   Section 3, contraste : « a far smaller common component ... (9.5 % contre
#     20.0 %) ». Le rapport passe de 2.25 a 2.11 -- le facteur deux tient, et
#     AUCUN rapport n'est cite en dur dans la phrase, qui se contente de
#     juxtaposer les deux macros. Rien a corriger.
#   Section 3, amenites : « almost the same figure (29.3 contre 31.0) » et
#     « nearer a third » -- macros inchangees, phrase intacte.
add('EMPetaAnnualProd', f"{prod_a.eta_pct.iloc[0]:.1f}", '9.7',  'table_eta.csv')

# Benchmark of the first row of tab:eta: eta on the RAW annual mean temperature,
# before any damage transformation. It lived in the table and its note without
# ever being cited by the prose, which invoked "a plausible order" against an
# external reference carrying no number, while the only numerical yardstick
# available sat in the same exhibit.
raw_t = d[d.metric == 'Annual mean T']
add('EMPetaAnnualMeanT', f"{raw_t.eta_pct.iloc[0]:.1f}", '43.8', 'table_eta.csv')

# ------------------------------------------------- ecart de courbe pondere ---
# ---------------------------------------------- weighted curve gap ------------
# Where the two functions differ, measured and weighted by the change in
# Northern exposure rather than asserted as a temperature range.
#
# Naming convention: a macro carrying one end of a range SAYS so in its name
# (Floor, Min, Max). Without that mark, a macro carrying a ceiling can be read
# as carrying a floor.
# LITTERAUX MIS A JOUR LE 2026-09-07 (6 des 10). Cause : la re-estimation de
# 01a avec les controles gravitaires, absents de ses huit sites d'appel
# jusqu'au 04/09. Le deplacement est concentre au BORD FROID.
#   ShareOfContrast 76.6->73.1 | ShareAboveTwenty 8->7
#   SuppMin 10->8 | SuppMax 11->9 | ColdRatioSuppMin 3.2->3.0 | Max 3.5->3.3
# ⚠️ Les quatre autres N'ONT PAS bouge, et c'est le controle : ExpoLostBelowZero
# et ExpoGainedAboveTwenty se calculent sur la seule exposition CMIP6, sans
# aucune courbe -- leur immobilite prouve que le deplacement vient des courbes.
# ⛔ Le temoin gele du 06/08 n'est pas rafraichi : l'ecart reste la mesure.
g = pd.read_csv(os.path.join(REPO, 'output', 'pe_weighted_gap_north_scalars.csv'))
g = dict(zip(g.key, g.value))
SRC_G = 'pe_weighted_gap_north_scalars.csv (04_contrast/ct3_weighted_curve_gap.py)'
add('EMPgapShareOfContrast',    f"{g['A_share_of_contrast_pct']:.1f}", '70.1', SRC_G)
add('EMPgapShareBelowZero',     f"{g['share_absA_below0_pct']:.0f}",   '76',   SRC_G)
add('EMPgapShareAboveTwenty',   f"{g['share_absA_above20_pct']:.0f}",  '8',    SRC_G)
add('EMPgapShareAboveTwentySuppMin',
    f"{g['share_absA_above20_supp_min_pct']:.0f}", '9', SRC_G)
add('EMPgapShareAboveTwentySuppMax',
    f"{g['share_absA_above20_supp_max_pct']:.0f}", '10', SRC_G)
add('EMPgapExpoLostBelowZero',  f"{g['expo_lost_below0_pp']:.1f}",     '5.1',  SRC_G)
add('EMPgapExpoGainedAboveTwenty',
    f"{g['expo_gained_above20_pp']:.1f}", '5.8', SRC_G)
add('EMPgapColdRatioFloor',     f"{g['cold_ratio_floor']:.1f}",        '1.8',  SRC_G)
add('EMPgapColdRatioSuppMin',   f"{g['cold_ratio_supp_min']:.1f}",     '3.0',  SRC_G)
add('EMPgapColdRatioSuppMax',   f"{g['cold_ratio_supp_max']:.1f}",     '3.4',  SRC_G)

# ---------------------------------------------------------------- Delta-rho --

# the canonical artefact of the current climate vintage.
d = pd.read_csv(os.path.join(OUT, 'rho_combined_by_country.csv'))
dr = d['delta_rho_a'].dropna()
add('EMPdrhoMedian',    f"{dr.median():.2f}",      '-0.31', 'rho_combined_by_country.csv (delta_rho_a)')
add('EMPdrhoMean',      f"{dr.mean():.2f}",        '-0.33', 'rho_combined_by_country.csv (delta_rho_a)')
add('EMPdrhoMedianAbs', f"{abs(dr.median()):.2f}", '0.31',  'same, absolute value (intro)')

# rho_hist is the historical LEVEL, not to be confused with the Delta-rho above,
# which is the CHANGE. Two objects, two artefacts.
d = pd.read_csv(os.path.join(OUT, 'rho_hist_by_country.csv'))
_rh = d['rho_a_daily'].dropna()
# SIGNED format (:+.2f) on all four. Writing the sign by hand in the .tex makes
# it a silent falsehood the day a maximum turns negative. The signed format keeps
# the macro right in both directions, and costs nothing while the sign holds.
add('EMProhistMedian', f"{_rh.median():+.2f}", '-0.65',
    'rho_hist_by_country.csv (rho_a_daily)')
add('EMProhistMean',   f"{_rh.mean():+.2f}",   '-0.60',
    'rho_hist_by_country.csv (rho_a_daily)')
add('EMProhistMin',    f"{_rh.min():+.2f}",    '-0.93',
    'rho_hist_by_country.csv (rho_a_daily)')
add('EMProhistMax',    f"{_rh.max():+.2f}",    '+0.59',
    'rho_hist_by_country.csv (rho_a_daily)')

# ---------------------------------------------------------------- AIC --------
# Deux regimes (16/07) : prod = selection decisive (delta AIC enormes),
# amenite = poids reellement disperses. Source machine = aic_weights_*.csv
# Emitted by the Monte Carlo of the R scripts, with a literal fallback until the
# CSV files are re-emitted by a run.
def _aic_from_csv(fname, mode):
    p = os.path.join(OUT, fname)
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    if mode == 'delta':
        return str(int(df.loc[df.delta_aic > 0, 'delta_aic'].min()))
    return f"{df.weight.max():.2f}"

# The AIC delta macros of the productivity channel have been retired along with
# the sentences that used them: at this sample size an Akaike weight of one is a
# property of N, not evidence about functional form. Emitting a macro nobody
# calls is the mirror image of the defect being repaired. The amenity weight
# stays: there the weights really are dispersed, which is what makes model
# averaging legitimate in that channel.
for _name, _f, _mode, _exp, _lit in [
    ('EMPamenAicMaxW',        'aic_weights_Amenity_Daily.csv',       'maxw',  '0.18', 'log MC 16/07 : 0.179')]:
    _v = _aic_from_csv(_f, _mode)
    if _v is None:
        add(_name, _exp, _exp, f'LITERAL — {_lit}; will be MATCHed on {_f} at the next R run', literal=True)
    else:
        add(_name, _v, _exp, _f)

# ---- AIC gap BETWEEN THE TWO ARMS -------------------------------------------
# A measurement, not an argument. The two arms share the source file, the sample,
# the dependent variable, the fixed effects, the gravity controls and the
# clustering: only the exposure metric on the right-hand side changes. The AIC is
# therefore comparable across the two tables.
#   - at an IDENTICAL cell, the number of parameters k is the same on both sides,
#     so the gap is exactly -2 times the log-likelihood difference. That is the
#     statistic to quote first.
#
# Sign: positive means the ANNUAL arm has the lower AIC.
#
# What these macros do NOT say is that the AIC selects anything. PPML is a
# quasi-likelihood, so the criterion is invalid here in BOTH directions. They
# state the gap; they never invoke it as proof.
# LA CELLULE EST LA REFERENCE DU PAPIER, ET ELLE EST LUE (18/09/2026).
# Elle etait retapee ici, en clair : une troisieme copie de la designation
# que ct2 et ct3 portaient deja chacun de leur cote. GRID_REF_SPEC est le
# pre-engagement (P4) de 00_parameters.py. La cellule existe dans les deux
# tables d'Akaike -- c'est la grille du CORPS, 12 specifications, qui
# contient la spline df=4 a 0,5 % -- et l'absence est traitee plus bas.
_AIC_CELL = _PARAMS.GRID_REF_SPEC
_ad = os.path.join(OUT, 'aic_weights_Productivity_Daily.csv')
_aa = os.path.join(OUT, 'aic_weights_Productivity_Annual.csv')
if os.path.exists(_ad) and os.path.exists(_aa):
    _dd = pd.read_csv(_ad).set_index('model_name')
    _da = pd.read_csv(_aa).set_index('model_name')
    _gap_cell = _dd.loc[_AIC_CELL, 'aic'] - _da.loc[_AIC_CELL, 'aic']
    _gap_best = _dd.aic.min() - _da.aic.min()
    # Rebases le 09/09/2026 sur le run 01a post-correctif GLORIA. Cause : le
    # panel commerce. Les deux ecarts se resserrent de ~3 % (56{,}508 ->
    # 54{,}811 et 58{,}399 -> 56{,}837), ce qui ne touche ni le SENS de l'ecart
    # -- le bras annuel garde l'AIC la plus basse -- ni l'ordre de grandeur sur
    # lequel la phrase s'appuie.
    add('EMPaicGapArmsCell', f"{_gap_cell:,.0f}".replace(',', '{,}'), '59{,}530',
        f'aic_weights_Productivity_{{Daily,Annual}}.csv, cell {_AIC_CELL} '
        '(same k -> -2 delta logL)')
    add('EMPaicGapArmsBest', f"{_gap_best:,.0f}".replace(',', '{,}'), '61{,}064',
        'aic_weights_Productivity_{Daily,Annual}.csv, min(aic) of each arm '
        '(different k: spline4 daily vs spline3 annual)')
else:
    add('EMPaicGapArmsCell', '56{,}508', '56{,}508',
        'LITERAL — aic_weights_Productivity_*.csv missing', literal=True)
    add('EMPaicGapArmsBest', '58{,}399', '58{,}399',
        'LITERAL — aic_weights_Productivity_*.csv missing', literal=True)

# ---- GE band of contribution B ----------------------------------------------
# The Northern contrast at the quantiles of the damage function, now produced
# INSIDE the package by 02_GE/02e_north_contrast_band.py from this grid's own
# welfare matrices. The producer checks that the runs share a run_id and an
# inputs_md5, which is what the previous, out-of-package producer bought with a
# 46-entry md5 table over four different run vintages.
#
# THREE MACROS WERE RETIRED HERE ON 09/09/2026, AND NOT BECAUSE THEY WERE
# INCONVENIENT. \GEbandQnfiveNneg, \GEbandQnfiveLeast and \GEbandOsc reported the
# number of OUTER iterations, the least negative among them, and their
# oscillation over iterations four to nine. 02c_model/main.py HAS NO OUTER LOOP:
# `iter_outer` and `i_outer` are initialised around l. 432 and never incremented
# or tested; the solver is a single forward sweep over t. The three macros
# described a solver this package does not contain, and the sentence they served
# defended against an artefact of where the solver was stopped -- there is no
# longer a stopping point to defend against. The sentence went with them.
#
# Restoring them would have meant re-importing the old outer loop, which moves
# EVERY GE number in the paper. That is a re-specification, not a repair, and it
# is not this emitter's call to make.
#
# TWO DECIMALS, never three. The old reason -- the outer loop does not converge
# -- died with the loop. The standing reason does not: this paper establishes the
# SIGN of the contrast and reports its magnitude to an order of magnitude.
_bp = os.path.join(OUT, 'ge_north_contrast_band.csv')
if os.path.exists(_bp):
    # dtype=str from the start: pandas would infer float64 and re-render the
    # rounding. The CSV already carries the exact text wanted. One formatting,
    # done by the producer.
    _b = pd.read_csv(_bp, dtype=str).set_index('quantity')['value']
    # regen=True on the three: the hard-codes were established on the
    # out-of-package runs of August, and these come from the package grid. The
    # gap is REPORTED, not blocking -- rule (b) of the runbook.
    for _nm, _k, _exp in (('GEbandQfifty', 'corrected_q50', '-0.13'),
                          ('GEbandQnfive', 'corrected_q95', '-0.06'),
                          ('GEbandQofive', 'corrected_q05', '-0.17')):
        if _k in _b.index:
            add(_nm, _b[_k], _exp, 'ge_north_contrast_band.csv (%s)' % _k, regen=True)
        else:
            # A mechanism that can do nothing must say whether it acted.
            print('[!] %s missing from ge_north_contrast_band.csv: %s NOT emitted.'
                  % (_k, _nm))
else:
    print('[!] ge_north_contrast_band.csv missing (02e_north_contrast_band.py): '
          'GEband* macros NOT emitted.')

# --- appendix OH: the denominator and numerator of the insularity sentence ----
# The count is EMITTED, and from the criterion the code ACTUALLY applies: the
# `straddle` column, that is `hi >= 0`. Do not key on the word "identified" in
# the table: that word names a CAUSE, the collapse of identifying variation,
# where the code tests a SYMPTOM. The two coincide on the current grid, and
# nothing guarantees they will after a re-estimation.
# ---- grille de contraste : taille et gain de variance ------------------------
# Meme fichier et meme colonne que la note de tab_contrast_grid : `se_ratio`,
# soit SE(difference) / SE(independance). Deux arms estimes sur les memes
# cellules ont des coefficients correles, et la variance de la difference porte
# cette covariance ; le rapport dit de combien le calcul d'independance est
# conservateur.
# ---- absorption des effets fixes produit x annee ----------------------------
# max de r2_pyr : la part de variance des termes climatiques que les effets
# produit x annee absorbent. Meme CSV que tab_fe_absorption, donc la macro ne
# peut pas contredire la table qui l'entoure.
_fap = os.path.join(OUT, 'fe_absorption.csv')
if not os.path.exists(_fap):
    raise SystemExit(
        "[STOP] fe_absorption.csv missing from output/\n"
        "  <- 01a_productivity_drf.R l.528")
_fa = pd.read_csv(_fap)
_col = 'r2_pyr_raw' if 'r2_pyr_raw' in _fa.columns else 'r2_pyr'
add('EMPprodYearAbsorpMax', f"{100 * _fa[_col].max():.2f}", '0.41',
    f'fe_absorption.csv, max of {_col} (x100)')

# ---- parametres declares du dispositif --------------------------------------
# Lus dans config_params.csv, depose par 02d depuis 00_parameters.R. Ce ne sont
# pas des resultats : ce sont les reglages, et ils doivent venir du fichier qui
# les declare.
_cfp = os.path.join(OUT, 'config_params.csv')
if not os.path.exists(_cfp):
    raise SystemExit(
        "[STOP] config_params.csv missing from output/\n"
        "  <- 02d_emit_latex.R. It carries MC_N_SIM, which governs the caption of\n"
        "  the fan charts. Hard-coded in the .tex it would become false at the\n"
        "  first change, silently.")
_cf = pd.read_csv(_cfp).set_index('param')['value']
add('EMPnMC', str(int(float(_cf.loc['MC_N_SIM']))), '300',
    'config_params.csv <- 00_parameters.R, MC_N_SIM (draws of the fan charts)')

_cgp = os.path.join(OUT, 'contribA_grid.csv')
if not os.path.exists(_cgp):
    raise SystemExit(
        "[STOP] contribA_grid.csv missing from output/\n"
        "  <- 04_contrast/ct4_emit_tab_contribA_grid.R\n"
        "  It carries the GRID_SPECS grid and the se_ratio column.")
_cg = pd.read_csv(_cgp)
if 'se_ratio' not in _cg.columns:
    raise SystemExit(
        "[STOP] contribA_grid.csv: column 'se_ratio' missing.\n"
        f"  columns: {list(_cg.columns)}\n"
        "  Producer: ct4_emit_tab_contribA_grid.R l.202.")

# LA GRILLE RECUE CONTRE LA GRILLE DECIDEE (18/09/2026).
# \EMPnSpecs et \EMPgridProdNspec comptaient les lignes du CSV recu, donc ce que
# ct1 avait tourne -- jamais ce que l'auteur avait decide. Une grille amputee en
# amont serait passee dans la table, dans la figure et dans les deux macros sans
# qu'aucune des trois ne bronche. 00_parameters.py est le tiers : il ne descend
# d'aucun run. On compare les ETIQUETTES, pas seulement le compte : quinze lignes
# dont une serait NS 5 au lieu de Poly 4 donneraient le bon nombre.
_grid_recue = set(_cg['spec_label'].astype(str))
_grid_voulue = set(_PARAMS.GRID_LABELS)
if _grid_recue != _grid_voulue:
    raise SystemExit(
        "[STOP] contribA_grid.csv does not carry GRID_SPECS.\n"
        "  received %d | GRID_SPECS %d\n"
        "  extra   : %s\n"
        "  missing : %s\n"
        "  Re-run ct1, ct2 then ct4 after any change of grid."
        % (len(_grid_recue), len(_grid_voulue),
           ', '.join(sorted(_grid_recue - _grid_voulue)) or '(none)',
           ', '.join(sorted(_grid_voulue - _grid_recue)) or '(none)'))
print('[grille] contribA_grid.csv = GRID_SPECS of 00_parameters.py (%d specifications)'
      % len(_grid_voulue))

add('EMPnSpecs', str(len(_cg)), '15',
    'contribA_grid.csv (rows = specifications of the contrast grid), '
    'checked against GRID_SPECS of 00_parameters.py')
add('EMPdiffSEshareLo', f"{100 * _cg['se_ratio'].min():.0f}", '53',
    'contribA_grid.csv, min of se_diff/se_diff_naive over the grid')
add('EMPdiffSEshareHi', f"{100 * _cg['se_ratio'].max():.0f}", '74',
    'contribA_grid.csv, max of se_diff/se_diff_naive over the grid')

# ---- les comptes de la grille, PRODUITS PAR L'EXHIBIT LUI-MEME ---------------
# Ces six nombres sont des enonces SUR LA FIGURE et SUR LA TABLE : « tant
# d'intervalles excluent zero ». Ils etaient recomptes ici, sur les CSV AMONT
# (test1_delta_north.csv, tache1_diff_bras.csv). Deux chemins de calcul pour un
# meme nombre, c'est deux nombres qui peuvent diverger sans que rien le dise.
# ct4 les calcule desormais sur `contribA_grid.csv` -- l'objet dont il trace les
# segments et imprime les lignes -- et les ASSERTE contre les segments traces.
# Ici on ne fait que les porter.
_cgcp = os.path.join(OUT, 'contribA_grid_counts.csv')
if not os.path.exists(_cgcp):
    raise SystemExit(
        "[STOP] contribA_grid_counts.csv missing from output/\n"
        "  <- 04_contrast/ct4_emit_tab_contribA_grid.R, same pass as\n"
        "     tab_contrast_grid.tex and fig_contrast_grid.pdf.")
_cgc = pd.read_csv(_cgcp).set_index('macro')['value']
SRC_GRIDC = ('contribA_grid_counts.csv <- ct4_emit_tab_contribA_grid.R, '
             'count asserted against the segments of fig_contrast_grid.pdf')


def _gc(name):
    """Le compte emis par ct4. Une macro absente est une etape qui n'a pas
    tourne, pas une valeur par defaut : on s'arrete."""
    if name not in _cgc.index:
        raise SystemExit(
            "[STOP] %s missing from contribA_grid_counts.csv (present: %s).\n"
            "  Re-run 04_contrast/ct4_emit_tab_contribA_grid.R."
            % (name, ', '.join(map(str, _cgc.index))))
    return str(int(_cgc.loc[name]))


# ---- la PHRASE des exceptions naives ----------------------------------------
# Le manuscrit la portait en dur, et elle etait FAUSSE (« the three degree-two
# polynomials » pour six exceptions). ct4 la derive de la grille ; on la porte.
def _gs(name):
    """Une valeur de chaine emise par ct4."""
    if name not in _cgc.index:
        raise SystemExit(
            "[STOP] %s missing from contribA_grid_counts.csv.\n"
            "  Re-run 04_contrast/ct4_emit_tab_contribA_grid.R." % name)
    return str(_cgc.loc[name])


add('EMPdiffExclNaiveExceptions', _gs('EMPdiffExclNaiveExceptions'),
    'the three degree-two polynomials', SRC_GRIDC, regen=True)

# ---- LA SPECIFICATION TRACEE PAR LES EVENTAILS -------------------------------
# ⚠️ IL Y EN A DEUX, ET CE N'EST PAS LA MEME. mc_run() (lib_drf.R l.303) tire la
# SPECIFICATION selon les poids d'Akaike avant de tirer beta. Ces poids sont
# DEGENERES -- 1,000 sur une seule specification dans chaque bras -- mais pas sur
# la MEME : le bras journalier tire spline4_wl0.5wu0.5, le bras annuel
# spline3_wl0.5wu0.5 (Delta AIC 2026 pour spline4, poids 0).
# Une macro unique \EMPfigDrfSpec ne peut donc pas etre vraie des deux panneaux.
# Les poids sont lus dans les artefacts que 01a archive lui-meme ; rien n'est
# recompte ici, et aucune estimation n'est relancee.
_AICF = {'daily': 'aic_weights_Productivity_Daily.csv',
         'annual': 'aic_weights_Productivity_Annual.csv'}
_drfspec = {}
for _arm, _fn in _AICF.items():
    _fp = os.path.join(OUT, _fn)
    if not os.path.exists(_fp):
        raise SystemExit(
            "[STOP] %s missing from output/\n"
            "  <- 01_estimation/01a_productivity_drf.R" % _fn)
    _aw = pd.read_csv(_fp).sort_values('weight', ascending=False)
    _top, _w = str(_aw.iloc[0]['model_name']), float(_aw.iloc[0]['weight'])
    # Un poids non degenere rendrait la phrase fausse : le Monte Carlo
    # melangerait alors plusieurs specifications et aucune ne serait « la »
    # specification tracee. On s'arrete plutot que de nommer la premiere.
    garde_phrase(_w >= 0.999,
        'arm %s: "the specification traced" (maximal Akaike weight >= 0.999)' % _arm,
        '>= 0.999', '%s weight %.4f' % (_top, _w),
        "[STOP] arm %s: maximal Akaike weight %.4f < 0.999.\n"
        "  The Monte Carlo mixes several specifications; \"the\n"
        "  specification traced\" no longer has a referent. A decision is needed."
        % (_arm, _w))
    _m = re.match(r'^(polynomial|spline)(\d+)_wl([0-9.]+)wu([0-9.]+)$', _top)
    if _m is None:
        raise SystemExit("[STOP] unrecognised specification label: %s" % _top)
    _form, _df, _wl, _wu = _m.group(1), int(_m.group(2)), _m.group(3), _m.group(4)
    if _wl != _wu:
        raise SystemExit("[STOP] asymmetric winsorisation: %s" % _top)
    _drfspec[_arm] = dict(
        form=_form, df=_df, wl=_wl,
        phrase=("%s with %d degrees of freedom, winsorised at %s\\%%"
                % ('a natural spline' if _form == 'spline' else
                   'a polynomial', _df, _wl)))

add('EMPfigDrfSpecDaily', _drfspec['daily']['phrase'],
    'a natural spline with 4 degrees of freedom, winsorised at 0.5\\%',
    'aic_weights_Productivity_Daily.csv <- 01a (weight degenerate at 1.000)',
    regen=True)
add('EMPfigDrfSpecAnnual', _drfspec['annual']['phrase'],
    'a natural spline with 3 degrees of freedom, winsorised at 0.5\\%',
    'aic_weights_Productivity_Annual.csv <- 01a (weight degenerate at 1.000)',
    regen=True)

# La chaine combinee, pour une legende qui couvre les DEUX panneaux.
if (_drfspec['daily']['form'] == _drfspec['annual']['form']
        and _drfspec['daily']['wl'] == _drfspec['annual']['wl']):
    _combo = ("natural splines winsorised at %s\\%%, with %d degrees of freedom "
              "in (a) and %d in (b)"
              % (_drfspec['daily']['wl'], _drfspec['daily']['df'],
                 _drfspec['annual']['df']))
else:
    _combo = ("%s in (a) and %s in (b)"
              % (_drfspec['daily']['phrase'], _drfspec['annual']['phrase']))
add('EMPfigDrfSpec', _combo,
    'natural splines winsorised at 0.5\\%, with 4 degrees of freedom in (a) and 3 in (b)',
    'aic_weights_Productivity_{Daily,Annual}.csv <- 01a', regen=True)

# ---- EMPohSurvivingShare : MESUREE LE 16/09/2026 -----------------------------
# Definition arbitree : moyenne sur les trois termes de (1 - R2) de la
# projection du differentiel climatique bilateral sur les EF exportateur x
# annee -- la construction de tab_fe_absorption, avec un autre jeu d'EF.
# Producteur : 01_estimation/01f_oh_surviving_share.R, qui REFUSE d'ecrire s'il
# ne retrouve pas les R2 paire x produit deja publies (ecart mesure 4,7e-05) :
# le chiffre neuf n'est publiable que si le chiffre connu revient.
#
# ⚠️ MESURE 49,9 %, CONTRE « roughly 69 % » DANS V10. Les deux ne sont pas
# reconcilies, et rien ici ne presume lequel est juste : le 69 % n'a ni
# producteur ni artefact, c'est pourquoi il a fallu le recalculer.
_ohs = os.path.join(OUT, 'oh_surviving_share.csv')
if not os.path.exists(_ohs):
    raise SystemExit(
        "[STOP] oh_surviving_share.csv missing from output/\n"
        "  <- 01_estimation/01f_oh_surviving_share.R")
_oh_df = pd.read_csv(_ohs)
_oh_tot = _oh_df[_oh_df['term'].isna()]
if len(_oh_tot) != 1:
    raise SystemExit(
        "[STOP] oh_surviving_share.csv: %d total row(s), 1 expected."
        % len(_oh_tot))
add('EMPohSurvivingShare', f"{float(_oh_tot.iloc[0]['surviving_share_pct']):.0f}",
    '69', 'oh_surviving_share.csv <- 01f_oh_surviving_share.R '
          '(tab_fe_absorption anchor checked)', regen=True)


_ohp = os.path.join(OUT, 'oh_hotslope_schemes.csv')
if os.path.exists(_ohp):
    _oh = pd.read_csv(_ohp)
    _n_all = len(_oh)
    # `identified` et non `straddle` (07/09) : c'est le nom que 01c_robustness.R
    # l.412 donne a la colonne, et sa polarite est INVERSE de straddle.
    if 'identified' not in _oh.columns:
        raise SystemExit(
            f"[STOP] {_ohp}: column 'identified' missing.\n"
            f"  columns present: {list(_oh.columns)}\n"
            f"  Producer: 01c_robustness.R l.412.")
    if _oh['identified'].dtype != bool:
        raise SystemExit(
            f"[STOP] {_ohp}: column 'identified' of type "
            f"{_oh['identified'].dtype}, boolean expected.\n"
            f"  Read as strings, .astype(bool) would return TRUE everywhere and the count\n"
            f"  would be wrong WITHOUT an error. Refusing rather than counting wrong.")
    # ------------------------------------------------------------------
    # UN SCHEME, PAS UNE LIGNE
    # ------------------------------------------------------------------
    # Defaut ferme le 17/09/2026. \EMPohSchemesN valait 5 et \EMPohSchemesId 3,
    # comptes en LIGNES de oh_hotslope_schemes.csv. Or la tache C reporte son
    # propre bras de controle : la ligne (C, withy) est byte-identique a la
    # ligne (A, baseline) -- memes q05/q50/q95 -- parce que c'est le MEME
    # scheme, sous deux etiquettes de tache. Le CSV decrit donc cinq lignes et
    # QUATRE schemes, et la table que la phrase accompagne en montre quatre :
    # Baseline, Exporter x year, Importer x year, Without income.
    #
    # Le compte identifie tombe de 3 a 2 pour la meme raison (baseline compte
    # deux fois), et c'est ce que la table montre : deux colonnes portent une
    # bande entre crochets, deux portent « --- ».
    #
    # La deduplication porte sur le TRIPLET de bande, pas sur les noms : deux
    # lignes qui rendent les trois memes quantiles decrivent la meme estimation,
    # quel que soit le nom de tache qui les a produites. Une table de
    # correspondance de noms serait une hypothese de plus a maintenir.
    _sch = _oh.drop_duplicates(subset=['slope_q50', 'slope_q05', 'slope_q95'])
    _n_sch = len(_sch)
    _n_id  = int(_sch['identified'].astype(bool).sum())

    # ------------------------------------------------------------------
    # ET ON LE CONFRONTE AUX COLONNES DE LA TABLE
    # ------------------------------------------------------------------
    # Le compte est un enonce SUR LA TABLE (« across the N schemes »). Le
    # deriver du CSV seul, meme deduplique, laisserait revenir l'ecart au
    # premier scheme ajoute d'un cote sans l'autre. On lit donc l'en-tete de
    # tab_oh_bounds.tex et on exige l'egalite -- la meme polarite que les
    # assertions de segments de ct4.
    _ohtp = os.path.join(OUT, 'tab_oh_bounds.tex')
    if not os.path.exists(_ohtp):
        raise SystemExit(
            "[STOP] tab_oh_bounds.tex missing: the count of schemes cannot "
            "be checked against the table it describes.\n"
            "  Producer: 01c_robustness.R (R4ter).")
    with open(_ohtp, encoding='utf-8') as _fh:
        _oht = _fh.read()
    _hdr = [l for l in _oht.split(chr(10))
            if l.strip().startswith('&') and 'Baseline' in l]
    if len(_hdr) != 1:
        raise SystemExit(
            "[STOP] tab_oh_bounds.tex: %d header line(s) recognised, one "
            "expected.\n  The count of columns cannot be established, so "
            "\\EMPohSchemesN cannot be checked." % len(_hdr))
    _ncol = len([c for c in _hdr[0].split('&')[1:]
                 if c.strip().strip(chr(92)).strip()])
    if _ncol != _n_sch:
        raise SystemExit(
            "[STOP] OH schemes: %d in oh_hotslope_schemes.csv after "
            "deduplication, %d columns in tab_oh_bounds.tex.\n"
            "  The sentence of C.2 would count something other than what the reader "
            "sees.\n  columns read: %s" % (_n_sch, _ncol, _hdr[0].strip()))
    print('[oh] %d rows -> %d distinct schemes, %d identified; '
          'table: %d columns (they agree)'
          % (_n_all, _n_sch, _n_id, _ncol))

    add('EMPohSchemesN',  str(_n_sch), '4',
        'oh_hotslope_schemes.csv, DISTINCT schemes (band triplet), '
        'count asserted against the columns of tab_oh_bounds.tex', regen=True)
    add('EMPohSchemesId', str(_n_id),  '2',
        'oh_hotslope_schemes.csv, distinct schemes with identified == TRUE '
        '(90% MC band excluding zero)', regen=True)
    ##########################################################################
    # LE COMPTE QUI PERD SA VARIATION -- CALCULE DEPUIS LE 20/09, EMIS LE 21.
    #
    # `_oh_lost` existe plus bas et ne servait qu'a un print et a une garde.
    # La phrase d'ouverture du Supplement B.2 dit « N schemes, <k> lose their
    # identifying variation ». Ecrire ce <k> en toutes lettres serait un
    # LITTERAL : c'est la regle posee le 20/09 a propos de \EMPosterDeltaNegN
    # -- un compte ecrit en mots est un compte que rien ne verifie.
    #
    # ⚠️ Il n'est PAS pose comme `_n_sch - _n_id`. Les deux comptes portent des
    # unites differentes : `_n_sch` deduplique sur le TRIPLET DE BANDE, tandis
    # que la perte de variation est une propriete du SCHEMA. Aujourd'hui la
    # soustraction donnerait le bon nombre par coincidence -- 4 - 2 = 2 -- et
    # elle cesserait d'etre juste le jour ou un doublon d'etiquetage tombe du
    # cote non identifie. On compte donc l'objet que le mot nomme.
    ##########################################################################
    ## CORRIGE (lot 61). Le compte etait pris sur les LIGNES de _oh pendant que
    ## EMPohSchemesN et EMPohSchemesId le sont sur les SCHEMES dedupliques (_sch) :
    ## le 29/09, withy (re-etiquetage de baseline, 01c l.557) est tombe du cote non
    ## identifie et le CSV a donne N = 4, Id = 1, Lost = 4. Les trois comptes portent
    ## desormais la meme unite, et le temoin l'exige.
    _n_lost = int((~_sch['identified'].astype(bool)).sum())
    if _n_id + _n_lost != _n_sch:
        raise SystemExit(
            "[STOP] OH schemes: identified %d + lost %d != schemes %d. The three "
            "counts do not carry the same unit." % (_n_id, _n_lost, _n_sch))
    add('EMPohSchemesLost', str(_n_lost), '2',
        'oh_hotslope_schemes.csv, DISTINCT schemes (band triplet) with identified == FALSE; guard: '
        "exporter must be the only identified scheme",
        regen=True)

    ##########################################################################
    # LA GARDE DU BLOC OH -- ET POURQUOI LES TROIS MACROS PROPOSEES SONT
    # RETIREES. Pose puis corrige le 20/09/2026.
    #
    # CE QUE J'AVAIS PROPOSE, ET QUI ETAIT FAUX. J'avais lu dans
    # oh_hotslope_schemes.csv « 5 schemas nommes, dont 3 identifies, qui ne
    # rendent que 2 bandes distinctes », et propose de publier les trois
    # comptes. Verifie en source : `withy` N'EST PAS UN SCHEMA DISTINCT.
    #   00_parameters.R l.417 : TRADE_STAGE3_CONTROLS = c(GRAVITY_CONTROLS,
    #                           "GDP_diff") -- le revenu est DEJA dans la
    #                           specification de reference ;
    #   01c l.557           : le bras 1 d'un test de RETRAIT prend
    #                           TRADE_STAGE3_CONTROLS tel quel.
    # Donc `withy` = meme effets fixes, memes controles, meme echantillon que
    # `baseline` : c'est la reference RE-ETIQUETEE en bras « with income » du
    # test de retrait. La bande identique AU BIT est une IDENTITE, pas une
    # coincidence, et publier « trois schemas qui n'en rendent que deux »
    # presenterait un doublon d'etiquetage comme un resultat.
    #
    # ⇒ LES TROIS MACROS SONT RETIREES. Les deux qui existaient valent juste :
    # EMPohSchemesN = 4 et EMPohSchemesId = 2 comptent les schemas REELLEMENT
    # distincts (baseline, exporter, importer, noy), ce que `emit_tab_oh_bounds`
    # sait deja puisqu'il compose exactement ces quatre colonnes et omet withy.
    #
    # CE QUI RESTE FRAGILE, ET CE QUE LA GARDE COUVRE. Le compte passe par une
    # deduplication sur le TRIPLET DE BANDE, pas sur le nom. Aujourd'hui les
    # deux coincident parce que l'unique doublon est un re-etiquetage. Si deux
    # schemas REELLEMENT differents rendaient un jour la meme bande, la
    # deduplication les fusionnerait en silence et le compte deviendrait faux.
    # On ARRETE donc si l'ensemble des coincidences n'est plus exactement
    # {baseline, withy}, ou si les schemas qui perdent leur variation ne sont
    # plus exactement {exporter, importer} -- les deux proprietes que la phrase
    # du Supplement B.2 nomme.
    ##########################################################################
    _oh_lost = sorted(_sch.loc[~_sch['identified'].astype(bool), 'scheme'].astype(str))   # lot 61 : schemes, pas lignes
    ## Lot 74 : l'attente « exporter, importer » encodait l'ancien monde (avant la
    ## reprise, exporter perdait sa bande). B.2 et §2.3 nomment desormais exporter
    ## x year comme le SEUL schema identifie : c'est ce que la garde verifie.
    import gardes_lot74 as _g74   # noqa: E402
    _oh_id = sorted(_sch.loc[_sch['identified'].astype(bool), 'scheme'].astype(str))
    _ok, _det = _g74.garde_b2(_oh_id, _n_id)
    garde_phrase(_ok,
        'Supplement B.2 and §2.3: exporter x year is the only identified scheme',
        'identified = {exporter}, EMPohSchemesId = 1', _det,
        "[STOP] OH schemes: %s. B.2 and §2.3 name exporter x year as the only "
        "scheme whose band excludes zero. Re-read the sentences before re-running."
        % _det)
    _dup = (_oh.groupby(['slope_q50', 'slope_q05', 'slope_q95'])['scheme']
            .apply(lambda s: frozenset(map(str, s))))
    _dup = sorted(set(d for d in _dup if len(d) > 1), key=sorted)
    ## Lot 57 : RESTE BLOQUANTE. Deux schemas reellement differents a bandes
    ## identiques signalent un defaut de CALCUL (la degenerescence du 07/09), pas
    ## seulement une phrase fausse.
    if _dup != [frozenset({'baseline', 'withy'})]:
        sys.exit(
            "[STOP] OH schemes: the only expected band coincidence is "
            "{baseline, withy}, which is a relabelling (01c l.557). "
            "Measured: %s. A coincidence between two REALLY different "
            "schemes would make EMPohSchemesN silently wrong."
            % ([sorted(d) for d in _dup] or 'none'))
    print("  [OH] %d distinct scheme(s) | %d identified | losing: %s | "
          "labelling duplicate: baseline == withy"
          % (_n_sch, _n_id, ', '.join(_oh_lost)))
    # Fourchette des pentes IDENTIFIEES. Meme colonne `identified` que le compte
    # ci-dessus : une macro qui citerait un autre sous-ensemble contredirait le
    # nombre imprime a cote d'elle dans la meme phrase.
    _id = _sch[_sch['identified'].astype(bool)]
    if not len(_id):
        print('[!] no identified scheme -> EMPohSlopeLo/Hi NOT emitted.')
    else:
        ######################################################################
        # LES BORNES VIENNENT DES ESTIMATIONS PONCTUELLES, PAS DU MONTE CARLO.
        # Arbitrage du 22/09/2026, tranche par une mesure.
        #
        # LE FAIT. `slope_q50` est la MEDIANE DE 3000 TIRAGES (01c l.357-365) :
        # elle herite du flux aleatoire, donc de la rotation de `psd_sqrt` sur
        # une `Sigma` mal conditionnee -- meme graine, meme flux, base tournee,
        # donc un AUTRE tirage de la MEME loi. Mesure du 22/09, rapport des
        # valeurs d'avant / apres la division par theta = 6,5 :
        #
        #   schema      identifie   beta ponctuel   q50 (median MC)
        #   baseline    oui         6.500000        6.50000
        #   withy       oui         6.500000        6.50000
        #   noy         oui         6.500000        6.56718   <-- il bouge
        #   exporter    non         6.500000        6.54831
        #   importer    non         6.500000        6.67022
        #
        # Le coefficient divise EXACTEMENT par theta sur les CINQ schemas, ecart
        # maximal 5,1e-11 : l'estimation est saine, et toute l'instabilite est
        # dans les tirages.
        #
        # ⇒ ON RAPPORTE LA PENTE ESTIMEE, PAS LA MEDIANE D'UNE SIMULATION DE LA
        #   PENTE. C'est l'objet que la phrase nomme -- « the slope lies between
        #   X and Y » -- il est exactement covariant par construction, et il ne
        #   bougera pas le jour ou `nearPD` tombera (liste post-soumission 13).
        #
        # ⛔ ET LA BANDE RESTE AU MONTE CARLO. `identified` ci-dessus en vient
        #   toujours : « whose band excludes zero » est une affirmation sur la
        #   LOI, pas sur le point, et c'est la que le Monte Carlo est a sa
        #   place. Les deux phrases du manuscrit font deja ce partage.
        ######################################################################
        _effp = os.path.join(OUT, 'oh_scheme_effects.csv')
        if not os.path.exists(_effp):
            sys.exit("[STOP] %s missing: EMPohSlopeLo/Hi now read the "
                     "POINT slope, which comes from it." % _effp)
        _effd = pd.read_csv(_effp)
        _wide = _effd.pivot_table(index='scheme', columns='temp_C',
                                  values='rel_effect_pct')
        # MEME FORMULE QUE 01c `slope_band()` l.360-362 : (y(30) - y(20)) / 10,
        # a ceci pres qu'elle s'applique a la courbe AJUSTEE et non a un tirage.
        for _t in (20, 30):
            if _t not in _wide.columns:
                sys.exit("[STOP] oh_scheme_effects.csv: temperature %d missing "
                         "from REF_T; the hot slope 20->30 cannot be "
                         "rebuilt. Columns: %s"
                         % (_t, sorted(_wide.columns)))
        _pente = (_wide[30] - _wide[20]) / 10.0
        _abs_id = sorted(set(_id['scheme']) - set(_pente.index))
        if _abs_id:
            sys.exit("[STOP] identified scheme(s) missing from "
                     "oh_scheme_effects.csv: %s. The two CSV files come from the same "
                     "run; a gap signals a mix of vintages."
                     % _abs_id)
        # Controle de coherence : le point et la mediane doivent au moins
        # s'accorder sur le SIGNE. S'ils divergent, ce n'est plus un ecart de
        # tirage, et la phrase « stays negative » ne tient plus.
        for _r in _id.itertuples():
            ## Lot 57 : RESTE BLOQUANTE -- coherence entre deux calculs (point
            ## estime et mediane Monte-Carlo), pas seulement une phrase.
            if _pente[_r.scheme] * float(_r.slope_q50) <= 0:
                sys.exit("[STOP] scheme %r: point slope %.6f and MC median "
                         "%.6f of opposite signs. The disagreement is no longer a "
                         "draw gap." % (_r.scheme, _pente[_r.scheme],
                                               float(_r.slope_q50)))
        # Lo = la plus NEGATIVE, Hi = la moins negative : l'intervalle se lit
        # [Lo, Hi] dans le manuscrit, donc il doit etre ordonne.
        _sel = _pente[list(_id['scheme'])]
        _lo, _hi = float(_sel.min()), float(_sel.max())
        ######################################################################
        # GARDE DE MILLESIME, AUTO-LEVEE. Posee le 21/09/2026.
        #
        # `oh_hotslope_schemes.csv` est le DERNIER artefact d'echelle d'avant la
        # division par theta : son md5 est identique a l'instantane du 19/09, et
        # `slope_q50` est une pente du canal PRODUCTIVITE, en points de log par
        # degre. Elle doit donc diviser par 6,5, et elle ne l'a pas fait --
        # parce que son producteur, `01c_robustness.R`, n'a pas tourne depuis
        # (8 h de calcul, 24 Go libres exiges).
        #
        # Publier -5,75 / -5,49 aujourd'hui serait imprimer une magnitude
        # perimee d'un facteur 6,5 dans une phrase par ailleurs correcte. On
        # emet donc \PROVISOIRE : le manuscrit rend « ?? » avec un warning
        # LaTeX, le controle 3 de check_macros nomme les deux macros en attente,
        # et AUCUNE valeur perimee n'atteint le PDF.
        #
        # ⚠️ ET ELLE SE LEVE TOUTE SEULE. Le declencheur est le md5 du fichier,
        # pas un drapeau a la main : le jour ou `01c` rejoue, l'empreinte change
        # et les deux macros s'emettent sans que personne ait a se souvenir de
        # retirer quoi que ce soit. Un drapeau manuel, lui, se serait oublie --
        # c'est la difference entre une garde et un pense-bete.
        ######################################################################
        _OH_MD5_AVANT_THETA = 'd86fd55169a3830bc936ab470d14c593'
        _oh_md5 = hashlib.md5(
            open(os.path.join(OUT, 'oh_hotslope_schemes.csv'), 'rb').read()).hexdigest()
        if _oh_md5 == _OH_MD5_AVANT_THETA:
            print("  [OH] vintage from BEFORE the division (md5 %s) -> "
                  "EMPohSlopeLo/Hi as PROVISOIRE. Lifted when 01c is re-run."
                  % _oh_md5[:8])
            add('EMPohSlopeLo', '\\PROVISOIRE', '\\PROVISOIRE',
                'oh_hotslope_schemes.csv AWAITING a re-run of 01c (pre-theta)')
            add('EMPohSlopeHi', '\\PROVISOIRE', '\\PROVISOIRE',
                'oh_hotslope_schemes.csv AWAITING a re-run of 01c (pre-theta)')
        else:
            # DEUX ACQUITTEMENTS LE 22/09/2026, ET ILS DISENT LA MEME CHOSE :
            # une epingle n'est une garde que si elle porte une MESURE.
            #
            # ① L'epingle valait `-0.88`. C'etait une PREDICTION -- -5.75/6.5 =
            #    -0.8850, arrondie a la main vers le bas -- pas une mesure. Le
            #    rejeu a rendu -0.885024336, qui s'arrondit a -0.89
            #    (|-0.885024 + 0.89| = 0.004976 < 0.005024) : MISMATCH, macro
            #    non emise, `Undefined control sequence` au manuscrit. La bonne
            #    polarite.
            # ② Puis la SOURCE a change le meme jour -- pente ponctuelle au lieu
            #    de mediane MC, voir le bloc ci-dessus -- et les bornes sont
            #    passees a -0.90 / -0.85. Second MISMATCH, second acquittement.
            #    Les deux epingles portent desormais des valeurs MESUREES :
            #    -0.898879 et -0.851839, arrondies a deux decimales.
            add('EMPohSlopeLo', f"{_lo:.2f}", '-1.05',
                'oh_scheme_effects.csv: POINT slope (y30-y20)/10, min over '
                'the schemes for which oh_hotslope_schemes.csv says identified == TRUE')
            add('EMPohSlopeHi', f"{_hi:.2f}", '-1.05',
                'oh_scheme_effects.csv: POINT slope (y30-y20)/10, max over '
                'the schemes for which oh_hotslope_schemes.csv says identified == TRUE')
            # \EMPohSlopeExporter (E1, lot 74) : la pente ponctuelle du schema que
            # B.2 et §2.3 NOMMENT. \EMPohSlopeLo/Hi sont un min/max sur un ensemble
            # de taille variable ; la phrase nomme un schema, elle lit sa ligne. La
            # garde B.2 ci-dessus assure qu'exporter est identifie.
            if 'exporter' not in _pente.index:
                sys.exit("[STOP] oh_scheme_effects.csv: no exporter row.")
            add('EMPohSlopeExporter', f"{float(_pente['exporter']):.2f}", '-1.05',
                'oh_scheme_effects.csv: POINT slope (y30-y20)/10 of the exporter scheme')
    # Nombre de tirages du bloc OH. Il gouverne la legende de fig_oh_asym_fe et
    # la note du texte, et il est passe de 300 a 3000 le 07/09 : en dur dans le
    # .tex il redeviendrait faux sans bruit. La colonne est CONSTANTE par
    # construction (01c la pose depuis MC_N_SIM_OH) -- on le VERIFIE avant de
    # lire une ligne, sinon lire la premiere serait lire au hasard.
    if 'n_sim' not in _oh.columns:
        print('[!] oh_hotslope_schemes.csv: column n_sim missing '
              '-> EMPohMCdraws NOT emitted (re-run of 01c required).')
    elif _oh['n_sim'].nunique() != 1:
        raise SystemExit(
            f"[STOP] {_ohp}: n_sim takes {_oh['n_sim'].nunique()} values "
            f"({sorted(_oh['n_sim'].unique())}).\n"
            f"  It is a parameter of the block; it cannot vary by scheme.\n"
            f"  Producer: 01c_robustness.R (R4bis), from MC_N_SIM_OH.")
    else:
        # LOT 75 : separateur de milliers « 3,000 », comme \GEncells (« 13,986 »).
        add('EMPohMCdraws', _thousands(_oh['n_sim'].iloc[0]), '3{,}000',
            'oh_hotslope_schemes.csv (n_sim = MC_N_SIM_OH of 00_parameters.R)')
    ##########################################################################
    # LA SECONDE GRAINE (lot 68, 30/09/2026)
    #
    # La phrase de B.2 disait depuis le 07/09 que deux graines independantes
    # s'accordent sur chaque schema. Elle n'avait pas de source dans le paquet :
    # la mesure du 07/09 portait sur un autre millesime. L'etape 01c_graine2 de
    # run_all.sh rejoue 01c avec la seconde graine (JMP_MC_SEED_PROD_DAILY) dans
    # output_01c_graine2/, et ces macros en descendent : le q95 de la pente chaude
    # de chaque schema DISTINCT sous chacune des deux graines, et le nombre de
    # schemas dont le verdict d'identification est le meme sous les deux.
    # Un ecart entre deux graines est UNE realisation, pas une dispersion : aucune
    # macro ici ne s'appelle erreur standard.
    #
    # Deux gardes BLOQUANTES (coherence entre deux calculs) : les memes schemas
    # et le meme n_sim de part et d'autre, et des estimations ponctuelles
    # identiques -- elles ne dependent pas de la graine, donc un ecart dirait que
    # les deux runs n'ont pas vu les memes donnees.
    ##########################################################################
    _oh2d = os.path.join(PKG, 'output_01c_graine2')
    _oh2p = os.path.join(_oh2d, 'oh_hotslope_schemes.csv')
    if not os.path.exists(_oh2p):
        print('[!] %s missing: EMPohSeed* NOT emitted. Step 01c_graine2 of '
              'run_all.sh.' % os.path.relpath(_oh2p, PKG))
    else:
        _oh2 = pd.read_csv(_oh2p)
        _k1 = sorted(zip(_oh['task'].astype(str), _oh['scheme'].astype(str)))
        _k2 = sorted(zip(_oh2['task'].astype(str), _oh2['scheme'].astype(str)))
        if _k1 != _k2 or _oh2['identified'].dtype != bool:
            raise SystemExit("[STOP] second seed: schemes %s against %s, or column "
                             "identified not boolean." % (_k2, _k1))
        if sorted(_oh2['n_sim'].unique()) != sorted(_oh['n_sim'].unique()):
            raise SystemExit("[STOP] second seed: n_sim %s against %s."
                             % (sorted(_oh2['n_sim'].unique()), sorted(_oh['n_sim'].unique())))
        _e1 = pd.read_csv(os.path.join(OUT, 'oh_scheme_effects.csv'))
        _e2 = pd.read_csv(os.path.join(_oh2d, 'oh_scheme_effects.csv'))
        _ee = _e1.merge(_e2, on=['task', 'scheme', 'temp_C'], suffixes=('_1', '_2'), how='outer')
        _dmax = float((_ee['rel_effect_pct_1'] - _ee['rel_effect_pct_2']).abs().max())
        if len(_ee) != len(_e1) or not (_dmax <= 1e-8):
            raise SystemExit("[STOP] second seed: the point estimates differ "
                             "(max gap %r over %d rows). They do not depend on the "
                             "seed: the two runs did not see the same data."
                             % (_dmax, len(_ee)))
        _NOM = {'baseline': 'Baseline', 'exporter': 'Exporter',
                'importer': 'Importer', 'noy': 'Noy'}
        # Valeurs mesurees le 30/09/2026 (graine 44 : output/ du 28/09 ; seconde
        # graine 20260930 : output_01c_graine2_20260930/).
        # LOT 75 : SIGNE EXPLICITE (« +0.016 » a cote de « -0.100 » dans B.2), comme
        # la famille n3 signee d'emit_values.py. Valeurs inchangees, format seul.
        # LOT 85o (07/10/2026) : _EXP1 et _EXP2 re-poses sur la course Cholesky
        # (logs/run_20261005_215131, rapport_lot85n §3a). Anciennes valeurs :
        # _EXP1 +0.016/-0.100/+1.192/+0.080 ; _EXP2 +0.032/-0.056/+1.223/+0.071.
        _EXP1 = {'baseline': '+0.001', 'exporter': '-0.036', 'importer': '+1.212', 'noy': '+0.061'}
        # LOT 85h (04/10/2026) : seconde graine re-posee sur la course sous renv (85f) ;
        # importer +1.244 -> +1.223 (bibliotheque utilisateur : +1.244, rapport_lot85g).
        _EXP2 = {'baseline': '+0.027', 'exporter': '-0.023', 'importer': '+1.253', 'noy': '+0.095'}
        _s2 = _oh2.set_index(['task', 'scheme'])
        _agree = 0
        for _r in _sch.itertuples():
            _sc = str(_r.scheme)
            if _sc not in _NOM:
                raise SystemExit("[STOP] second seed: distinct scheme %r without a macro "
                                 "name (known: %s)." % (_sc, sorted(_NOM)))
            _r2 = _s2.loc[(_r.task, _r.scheme)]
            add('EMPohSeedOneQnfive' + _NOM[_sc], '%+.3f' % float(_r.slope_q95), _EXP1[_sc],
                'oh_hotslope_schemes.csv, slope_q95 of scheme %s (seed MC_SEED prod_daily)' % _sc)
            add('EMPohSeedTwoQnfive' + _NOM[_sc], '%+.3f' % float(_r2['slope_q95']), _EXP2[_sc],
                'output_01c_graine2/oh_hotslope_schemes.csv, slope_q95 of scheme %s '
                '(second seed, step 01c_graine2)' % _sc)
            _agree += int(bool(_r.identified) == bool(_r2['identified']))
        add('EMPohSeedAgreeN', str(_agree), '4',
            'DISTINCT schemes whose identified is the same under the two seeds '
            '(output/ and output_01c_graine2/oh_hotslope_schemes.csv)')
        ######################################################################
        # \EMPohSeedMoveMax (lot 74) ET LA GARDE DES GRAINES.
        # DEFINITION (author's decision, lot 74-C) : le MAXIMUM, sur les QUATRE schemas
        # distincts (baseline, exporter, importer, noy), de |q95 graine 1 - q95
        # graine 2|. Le 01/10 il est atteint par importer (0,052 ; exporter 0,044).
        # LOT 85h : sous renv, importer 0,031 ; le maximum passe a exporter (0,044).
        # LOT 85o : apres le tirage Cholesky, importer 0,041 porte de nouveau le maximum.
        # B.2 le dit : « the largest such move across the four schemes is ... ».
        # Garde de phrase : q95(baseline) <= ecart max < q95(noy) (« within which the
        # baseline's percentile lies and the other's does not »). La clause « marge du
        # baseline de l'ordre de son ecart » est retiree au lot 85o : sa phrase (§2.3)
        # a quitte le .tex au lot 81.
        ######################################################################
        _CITES = ['baseline', 'exporter', 'importer', 'noy']
        _q1 = {str(r.scheme): float(r.slope_q95) for r in _sch.itertuples()}
        _q2 = {str(r.scheme): float(_s2.loc[(r.task, r.scheme)]['slope_q95']) for r in _sch.itertuples()}
        if not set(_CITES) <= set(_q1):
            raise SystemExit("[STOP] seeds: cited schemes %s missing from %s." % (_CITES, sorted(_q1)))
        _mv = {s: abs(_q1[s] - _q2[s]) for s in _CITES}
        _mvmax = max(_mv.values())
        if set(_CITES) != set(_q1):
            raise SystemExit("[STOP] seeds: distinct schemes %s, four expected %s."
                             % (sorted(_q1), _CITES))
        add('EMPohSeedMoveMax', '%.3f' % _mvmax, '0.041',
            'max |q95 seed 1 - q95 seed 2| over the four distinct schemes (%s); reached by %s'
            % (', '.join(_CITES), max(_mv, key=_mv.get)))
        _ok, _det = _g74.garde_graines(_q1['baseline'], _q1['noy'], _mvmax)
        garde_phrase(_ok, 'B.2: the baseline band within the gap between seeds, not the noy band',
                     'q95 base <= max gap < q95 noy', _det,
                     "[STOP] seeds: %s. The sentence of B.2 on the gap between two "
                     "seeds no longer holds." % _det)
        print('  [OH] second seed: %d/%d distinct schemes with the same verdict; '
              'point estimates equal to within %.1e' % (_agree, _n_sch, _dmax))
else:
    print('[!] oh_hotslope_schemes.csv missing: EMPohSchemes* macros NOT emitted. '
          'Re-run robustness_oh.R (task INS) -- it re-reads the three '
          'oh_task*_objects.rds and re-estimates nothing.')

# Size of the specification grid, and therefore the DENOMINATOR of every
# invariance sentence. It is never written "all twelve" by hand: the grid has
# already changed once, and a typed denominator survives such a change without
# saying so.
_g = [len(pd.read_csv(os.path.join(OUT, f))) for f in
      ('aic_weights_Productivity_Daily.csv', 'aic_weights_Productivity_Annual.csv',
       'pe_spec_ranking.csv') if os.path.exists(os.path.join(OUT, f))]
# QUATRIEME SOURCE, ajoutee le 22/09/2026 avec D3. La phrase de §2.3 ecrit
# desormais « \EMPosterDeltaNegN of the \EMPgridNspec specifications return a
# negative delta* » : le denominateur y compte la grille OSTER. S'il etait lu
# sur trois fichiers qui ne la contiennent pas, la phrase attribuerait a une
# grille un compte etabli sur une autre -- le defaut de §3.179, mais sur le
# denominateur. On verse donc oster_delta_full.csv dans la concordance : une
# divergence devient bruyante au lieu d'etre attribuee en silence.
_g.append(_OSTER_NSPEC)
if _g and len(set(_g)) == 1:
    add('EMPgridNspec', str(_g[0]), '12',
        'aic_weights_Productivity_{Daily,Annual}.csv, pe_spec_ranking.csv and '
        'oster_delta_full.csv, four sources in agreement')
else:
    add('EMPgridNspec', '12', '12',
        'LITERAL -- sources disagree or are missing: %s' % _g, literal=True)

##########################################################################
# LA PHRASE DE §2.2 AFFIRME UNE RELATION ENTRE DEUX COMPTES INDEPENDANTS.
# Posee le 22/09/2026, et c'est « la prose affirme, la garde verifie ».
#
# §2.2 ecrit : « the identification diagnostics rest on the forms of the
# published grid; the contrast ADDS a fourth-degree polynomial, so they
# report \EMPgridNspec specifications where it reports \EMPnSpecs ».
#
# ⛔ OR LES DEUX COMPTES NE VIENNENT PAS DU MEME RECENSEMENT, et c'est
#   voulu : chacun compte l'objet que son mot nomme. \EMPnSpecs lit les
#   lignes de contribA_grid.csv ; \EMPgridNspec lit la concordance de
#   quatre artefacts de diagnostic. Leur RELATION, elle, n'etait affirmee
#   que par la phrase. Si la grille bougeait de sorte que la relation casse,
#   chaque macro suivrait sa propre source et la phrase deviendrait fausse
#   EN SILENCE -- exactement le mode de defaillance de ce dossier.
#
# ⇒ ON VERIFIE LES DEUX MOITIES DE L'AFFIRMATION :
#   1. les formes des diagnostics sont un SOUS-ENSEMBLE STRICT de celles du
#      contraste ;
#   2. ce que le contraste AJOUTE est exactement le polynome de degre 4.
##########################################################################
try:
    _f_diag = set(zip(_od['spec'], _od['degree'].astype(int)))
    _cgf = _cg[['form', 'df']].drop_duplicates()
    _f_ctr = set(zip(_cgf['form'], _cgf['df'].astype(int)))
except Exception as _e:
    sys.exit('[STOP] additive relation of the two grids: columns '
             'unreadable (%s). A check that cannot conclude does not '
             'pass.' % _e)
if not _f_diag < _f_ctr:
    sys.exit(
        "[STOP] §2.2 states that the diagnostics rest on the forms of\n"
        "  the published grid and that the contrast ADDS to them. That is no longer\n"
        "  a strict subset.\n"
        "    diagnostics: %s\n    contrast   : %s"
        % (sorted(_f_diag), sorted(_f_ctr)))
_ajout = _f_ctr - _f_diag
# ⛔ ON NE NOMME PAS LA FORME ATTENDUE ICI -- ce serait une entree de grille en
#   dur, et le controle 9 de check_macros la refuse a juste titre. On verifie la
#   PROPRIETE que la phrase affirme -- « le contraste AJOUTE une forme », donc
#   exactement UNE -- et on EMET son nom. §2.2 cite la macro : le jour ou la
#   forme ajoutee change, la phrase suit, et le jour ou il y en a deux, on
#   s'arrete au lieu de publier un singulier faux.
if len(_ajout) != 1:
    sys.exit(
        "[STOP] §2.2 states that the contrast adds ONE form to the diagnostics.\n"
        "  Measured: %d form(s) added, %s.\n"
        "  Correct the sentence, or the grid."
        % (len(_ajout), sorted(_ajout) or '(none)'))
_a_spec, _a_deg = sorted(_ajout)[0]
_ORD = {2: 'second', 3: 'third', 4: 'fourth', 5: 'fifth', 6: 'sixth'}
if _a_deg not in _ORD:
    sys.exit("[STOP] degree %d has no ordinal: add the word rather than "
             "let the macro print a bare digit." % _a_deg)
if _a_spec == 'polynomial':
    _a_mots = 'a %s-degree polynomial' % _ORD[_a_deg]
elif _a_spec == 'spline':
    _a_mots = 'a natural spline with %s degrees of freedom' % _ORD[_a_deg]
else:
    sys.exit("[STOP] unknown form family: %r" % _a_spec)
add('EMPgridAddedForm', _a_mots, 'a fourth-degree polynomial',
    'contribA_grid.csv minus oster_delta_full.csv: the form that the contrast '
    'adds to the diagnostics, put into words')
print('[grilles] diagnostics %d forms, contrast %d forms; '
      'measured addition: %s -- the sentence of §2.2 holds'
      % (len(_f_diag), len(_f_ctr), sorted(_ajout)))
##########################################################################
# GARDE B6 (lot 78) -- §2.2 dit desormais que la forme ajoutee « carries no
# Akaike weight in either arm ». Lue dans la sortie de ct1 (poids d'Akaike sur
# la grille du contraste, par bras). Seuil 1e-6 dans les deux bras. Mesure du
# 02/10/2026 : poids nuls dans les deux bras, ecart d'AIC minimal 3 454 points
# (journalier) et 4 152 (annuel). Connu-faux : tools/test_gardes_lot78.py.
##########################################################################
import gardes_lot78 as _g78   # noqa: E402
_t1aic = pd.read_csv(os.path.join(OUT, 'contribA', 'test1_delta_north.csv'))
_t1aic = _t1aic[(_t1aic.channel == 'prod') & (_t1aic.status == 'OK')]
_w_ajout = {b: list(_t1aic.loc[(_t1aic.arm == b) & (_t1aic.form == _a_spec)
                               & (_t1aic.df.astype(int) == _a_deg), 'aic_weight'].astype(float))
            for b in ('daily', 'annual')}
_ok_b6, _det_b6 = _g78.garde_b6(_w_ajout)
print('[grilles] B6: %s' % _det_b6)
garde_phrase(_ok_b6, 'B6 Akaike weight of the added form', '< 1e-6', _det_b6,
             "[STOP] §2.2 (B6): %s. The sentence \"that form carries no Akaike weight in "
             "either arm\" is false." % _det_b6)

# La taille de grille est RETENUE pour le bloc Finlande, qui cite desormais
# \EMPgridNspec comme denominateur d'une fourchette. Un denominateur affiche par
# une macro et une fourchette calculee sur un autre nombre de specifications
# diraient deux choses differentes dans la meme phrase. None = grille non
# etablie : la garde en aval se desarme au lieu de comparer a un nombre invente.
_GRID_N = _g[0] if (_g and len(set(_g)) == 1) else None

# ----------------------------------------------------------------------
# LA QUEUE CHAUDE A TRAVERS LA GRILLE -- SUR LA GRILLE QUE LE TEXTE NOMME
# ----------------------------------------------------------------------
# Defaut ferme le 17/09/2026. \EMPspecInvarShare valait 91, mesure sur
# pe_spec_ranking.csv, qui porte DOUZE specifications : Poly 2-3 et NS 3-4,
# chacune a trois winsorisations. La phrase qui la citait dit « across the
# grid », et la grille que \EMPgridProdNspec compte n'etait pas celle-la.
# Le correctif n'a pas ete de changer le nombre : il a ete de changer la
# SOURCE. Ces cinq macros descendent de contribA_grid.csv, c'est-a-dire de la
# grille que le texte nomme, quelle qu'elle soit.
#
# 🔑 ET C'EST BIEN LA SOURCE QUI COMPTAIT, LE 18/09 LE MONTRE. La grille passe
# de 24 a 15 specifications (decision d'auteur : les NS 5-6-7 oscillent dans la
# partie dense du support). \EMPspecInvarShare revient a 91 -- la valeur meme
# que le defaut du 17/09 avait ecartee. Ce n'est pas un retour en arriere :
#   - le 17/09, 91 etait mesure sur DOUZE specifications sous une phrase qui en
#     annoncait vingt-quatre. Le nombre etait juste, la phrase le donnait pour
#     autre chose ;
#   - le 18/09, 91 est mesure sur les QUINZE que la phrase annonce.
# Une valeur identique sous une provenance differente reste une autre valeur.
#
# Mesure, colonne fz30 de contribA_grid.csv (courbe normalisee a max = 0,
# en unites log ; x100 = points de log) :
#   15 lignes sur 15 de signe NEGATIF     -> le SIGNE est invariant
#   de -1.233545 (NS 4 [0.5%]) a -1.118996 (NS 3 [2.5%])
#   rapport d'amplitude 1,1 ; la plus plate vaut 91 % de la plus raide
#   (grille large, pour comparaison : de -1.554382 a -0.218316, rapport 7,1,
#    la plus plate a 14 % -- les deux extremes etaient NS 5 et NS 6)
#
# Ce que l'enonce peut dire, et ce qu'il ne pouvait pas dire avant : sur la
# grille a 15, « estimated stably across the grid » est vrai du SIGNE ET de la
# MAGNITUDE. Les deux macros de bornes restent, pour que la phrase dise
# laquelle des deux elle affirme plutot que de le laisser deviner.
#
_cg30 = pd.read_csv(os.path.join(OUT, 'contribA_grid.csv'))
if 'fz30' not in _cg30.columns:
    raise SystemExit(
        "[STOP] contribA_grid.csv: column 'fz30' missing.\n"
        f"  columns: {list(_cg30.columns)}\n"
        "  Producer: 04_contrast/ct4_emit_tab_contribA_grid.R.")
_fz30 = _cg30['fz30'].astype(float)
if _fz30.isna().any():
    raise SystemExit("[STOP] contribA_grid.csv: fz30 has NA values; the "
                     "sign count would bear on a silent subset.")
# Le SIGNE, d'abord : c'est lui que la phrase affirme desormais.
_n_neg30 = int((_fz30 < 0).sum())
if _n_neg30 != len(_fz30):
    print("[!] hot tail: %d specification(s) out of %d with a POSITIVE sign at "
          "30 degC -- the statement of sign invariance no longer holds."
          % (len(_fz30) - _n_neg30, len(_fz30)))
add('EMPspecInvarNneg', str(_n_neg30), '15',
    'contribA_grid.csv, fz30 < 0 over the whole grid (GRID_SPECS)', regen=True)
# Les deux BORNES, signees, en points de log : la phrase les lit « from ... to ».
add('EMPspecInvarMin', f"{100 * _fz30.min():.1f}", '-123.4',
    'contribA_grid.csv, min of fz30 over the whole grid (log points)',
    regen=True)
add('EMPspecInvarMax', f"{100 * _fz30.max():.1f}", '-111.9',
    'contribA_grid.csv, max of fz30 over the whole grid (log points)',
    regen=True)
# Le RAPPORT d'amplitude, pour que le texte n'ecrive pas « threefold » a vue.
add('EMPspecInvarFold', f"{_fz30.abs().max() / _fz30.abs().min():.1f}", '1.1',
    'contribA_grid.csv, |fz30| max / |fz30| min over the whole grid',
    regen=True)
# La part conservee par la plus PLATE : meme definition qu'avant, autre grille.
_share = _fz30.abs().min() / _fz30.abs().max()
# REBASAGE DU 10/09/2026 : 89 -> 91. Brut mesure : plus plate -1.11899602,
# plus raide -1.23354522, rapport 0.90714. L'invariance se RENFORCE.
#
# ⚠️ TROIS SITES CONSOMMATEURS RELUS -- deux dans JMP_V10.tex, un dans le
# corps genere de tab_spec_invariance.tex. AUCUN NE FIGE DE SEUIL : les trois
# ecrivent « the smallest is/retaining X% of the largest », c'est-a-dire la
# valeur elle-meme, sans « au moins », sans « plus de », sans nombre rond en
# regard. Une valeur qui monte ne peut donc rien y invalider. Rien a reecrire.
add('EMPspecInvarShare', f"{100*_share:.0f}", '91',
    'contribA_grid.csv, |fz30| min / |fz30| max over the whole grid',
    regen=True)

# ---- hot-tail drop a 30C (site l.~757 ; string validee 16/07 soir) ----------
# Definition: the value of the daily dose-response at 30 C MINUS the mean of the
# cold plateau over 5 to 10 C. The unit is log points; the level equivalent is
# 1 - exp(drop), RELATIVE TO THE SAME PLATEAU, not to the optimum.
d = pd.read_csv(os.path.join(OUT, 'damage_prod_bin.csv'))
# 'f_z_raw' et non 'damage' (07/09) : c'est le nom que 01a_productivity_drf.R
# l.686 donne a la courbe normalisee a max = 0. Meme objet, autre nom.
if 'f_z_raw' not in d.columns:
    raise SystemExit(
        f"[STOP] damage_prod_bin.csv: column 'f_z_raw' missing.\n"
        f"  columns present: {list(d.columns)}\n"
        f"  Producer: 01a_productivity_drf.R l.686.")
_plateau = d.loc[(d.temperature_C >= 5) & (d.temperature_C <= 10), 'f_z_raw'].mean()
_v30  = d.loc[d.temperature_C == 30, 'f_z_raw'].iloc[0]
_drop = _v30 - _plateau
add('EMPhotDropThirty',    f"{abs(100*_drop):.0f}",            '19',
    'damage_prod_bin.csv (v30 - mean of the 5-10C plateau)')
add('EMPhotDropThirtyPct', f"{100*(1-math.exp(_drop)):.0f}",   '17',
    'same, level effect 1-exp(drop), relative to the plateau')

# ---- Oster : mouvement relatif spec AIC + beta_tilde du poly2[1%] -----------
d = pd.read_csv(os.path.join(OUT, 'oster_delta_full.csv'))
# MEME CELLULE QUE \EMPosterDelta, et pour la meme raison : la phrase du
# manuscrit lie les deux (« the reference specification, delta* = ..., whose
# coefficient moves by ... »). Les faire pointer deux cellules differentes
# ferait dire a la phrase autre chose que ce qu'elle affirme.
_r = d[d.tag == _OSTER['tag_moverel']].iloc[0]   # lot 63 : meme cellule que \EMPosterDelta (garde_oster, point 5)
# `expected` = la nouvelle cellule, meme raison qu'a \EMPosterDelta ci-dessus :
# un MISMATCH n'est pas emis. Ancienne cellule spline4_wl0.5wu0.5 = 5,2.
add('EMPosterMoveRel', f"{100*abs((_r.b_dot - _r.b_tilde)/_r.b_tilde):.1f}", '89.7',
    'oster_delta_full.csv (%s, |b_dot-b_tilde|/|b_tilde|)' % _OSTER['tag_moverel'])
# \EMPosterMoveRelRef (lot 74) : le meme mouvement relatif, sur la specification de
# REFERENCE de 00_parameters.R -- la cellule de \EMPosterDeltaRef. La phrase de §2.3
# oppose les deux mouvements pour dire pourquoi le minimum est petit : delta* porte
# le facteur b_tilde / (b_dot - b_tilde) (01c_robustness.R, run_oster), donc
# l'inverse de ce mouvement.
_rr = d[d.tag == _REF_TAG].iloc[0]
add('EMPosterMoveRelRef', f"{100*abs((_rr.b_dot - _rr.b_tilde)/_rr.b_tilde):.1f}", '12.9',
    'oster_delta_full.csv (%s = reference of 00_parameters.R, |b_dot-b_tilde|/|b_tilde|)' % _REF_TAG)
# ---- \EMPosterPolyTwoBtilde RETIREE le 10/09/2026 --------------------------
# La macro lisait oster_delta_full.csv (polynomial2_wl1wu1, b_tilde) et sortait
# en MISMATCH : +0.005 attendu, -0.126 mesure -- un changement de SIGNE.
#
# Elle n'est pas rebasee, elle part : elle n'a plus de consommateur. Mesure du
# 10/09 : `grep -c EMPosterPolyTwoBtilde JMP_V10.tex` rend 0, contre 1 dans les
# trois sauvegardes anterieures (bak_20260909_avant_lot1, bak_20260909_1613,
# chain_20260907). La phrase de l'entree quasi nulle a ete supprimee au lot
# precedent ; la macro, elle, est restee.
#
# Une valeur calculee puis jamais lue est indiscernable d'une valeur inventee :
# aucune phrase ne la confronte, donc rien n'aurait signale ce changement de
# signe. C'est le motif que ce paquet ferme partout ailleurs. Le bloc part en
# ENTIER -- l'etiquette, son garde-fou et la lecture -- parce que retirer le seul
# `add` laisserait trois lectures orphelines qui ne servent plus personne.
#
# Producteur inchange : 01c_robustness.R (section R1). La cellule redevient
# citable en re-declarant ce bloc, le jour ou une phrase la reclame.

# ---- amenity MDE block -- 6 macros, computed on the unrounded chain ---------
# sigma_eps is on the LOG-RATIO scale, the standard deviation of the residuals
# of log(y/mu) over strictly positive flows. That scale is the one consistent
# with a coefficient on the log scale of the PPML; a deviance or Pearson residual
# would give a minimum detectable slope about a hundred times too small, and a
# flatteringly low ratio.
# N_eff is the number of origin-destination pairs, that is the size of the pair
# fixed effect on the STAGE 3 sample. The power rests on the number of clusters,
# not on the larger flow count of stage 1.
# Largest admissible slope: EMPhotDropThirty / (30 - 7.5), where 7.5 is the
# midpoint of the plateau.
_mde   = pd.read_csv(os.path.join(OUT, 'mde_amenity_inputs.csv'))
_best  = _mde.loc[_mde.weight.idxmax()]
_fs    = pd.read_csv(os.path.join(OUT, 'mde_amenity_fixef_sizes.csv'))
_neff  = int(_fs.loc[_fs.fixef == 'fe_orig_dest', 'size'].iloc[0])
_sig   = float(_best.sd_logratio)
_OMEGA = float(_PARAMS.OMEGA)   # jumeau de Omega, 00_parameters.R (lot 38)
# sigma_{Delta T} : ecart-type INTRA-PAIRE du differentiel bilateral de
# Temperature variability over quinquennial windows. It is a MEASURED quantity,
# not a setting, but no artefact of this package carries it yet; it is therefore
# a literal, marked as such.
# manuscript at four sites. Two independent copies of the same measured number
# are free to diverge; there is now only one.
_SIGMA_DT = 0.15
add('EMPsigmaDeltaT', f"{_SIGMA_DT:.2f}", '0.15',
    'LITERAL — within-pair SD of the bilateral temperature differential; '
    'to be emitted from drf_amenity_final.r', literal=True)
_gmin  = 2.80 * _OMEGA * _sig / (_SIGMA_DT * math.sqrt(_neff)) * 100   # log pts/degC
_slope = abs(100 * _drop) / (30 - 7.5)
_ratio = _gmin / _slope
_ci    = 1.96 * _OMEGA * _sig / (_SIGMA_DT * math.sqrt(_neff)) * 100
add('EMPmdeSigmaEps',  f"{_sig:.2f}",  '1.62', 'mde_amenity_inputs.csv (sd_logratio, max-weight spec)')
add('EMPmdeNeff',      f"{_neff:,}".replace(',', '{,}'), '13{,}929',
    'mde_amenity_fixef_sizes.csv (fe_orig_dest)')
add('EMPmdeGprimeMin', f"{_gmin:.1f}", '12.8', 'eq:mde_amenity : 2.80*Omega*sigma/(0.15*sqrt(Neff))')
add('EMPmdeSlopeMax',  f"{_slope:.1f}", '0.8', 'EMPhotDropThirty / (30-7.5)')
add('EMPmdeRatio',     f"{_ratio:.1f}", '15.6', 'gmin/slope_max, unrounded chain')
add('EMPmdeCI',        f"{_ci:.1f}",    '9.0', '1.96*Omega*sigma/(0.15*sqrt(Neff))')

# ---- de combien f_a chapeau depasse la restriction d'echelle ---------------
# POSE LE 20/09/2026, ET CE N'EST PAS UNE CONSEQUENCE DE LA DIVISION : C'EST CE
# QU'ELLE A REVELE.
# §2.4 pose |f_a| <= |f_z chapeau| -- une restriction en unites de BIEN-ETRE, car
# u = amenite x revenu et l'eq. (2) fait des deux des decalages en log. Elle exige
# donc f_z en unites de productivite FONDAMENTALE. Tant que f_z portait theta en
# trop, la borne etait gonflee de 6,5, et c'est la SEULE raison pour laquelle le
# point estime passait : 0,91 avant, 5,93 apres.
# Le numerateur est la courbe d'amenite QUE LE GE LIT pour amenity=baseline
# (MC_damage_amenity_daily_quantiles.csv, colonne q50), pas la ponderee d'Akaike
# (damage_amenity_bin.csv, qui donne 5,10) : la macro doit decrire le bras qui
# tourne, pas un voisin.
# Meme construction que EMPmdeSlopeMax des deux cotes -- plateau moyen sur
# [5, 10] degC, chute jusqu'a 30 degC -- donc le (30 - 7.5) se simplifie et le
# rapport est celui des chutes.
_amen_curve = pd.read_csv(os.path.join(OUT, 'MC_damage_amenity_daily_quantiles.csv'))
_a_pl = _amen_curve.loc[(_amen_curve.temperature_C >= 5)
                        & (_amen_curve.temperature_C <= 10), 'q50'].mean()
_a_30 = float(np.interp(30, _amen_curve.temperature_C, _amen_curve.q50))
_a_drop = _a_30 - _a_pl
if _drop == 0:
    raise SystemExit('[STOP] zero drop of f_z: the scale ratio has no meaning.')
_excess = abs(_a_drop) / abs(_drop)
##############################################################################
# ⛔ CE QUE CETTE MACRO EST, ET LE MOT QUE LA PROSE NE DOIT PAS EMPLOYER.
# Ajoute le 21/09/2026, decision d'auteur.
#
# C'est un rapport de CHUTES EN NIVEAU : f(30 degC) - moyenne(f sur [5,10] degC),
# pour l'amenite sur la productivite, sur LA MEME PLAGE. Ce n'est PAS un rapport
# de pentes. La mesure du 20/09 (RECAP §3.154) donne les deux :
#     en pente  : 6,02      en niveau : 5,87  ->  la macro vaut 5,9
# Les deux disent la meme chose ; ils ne se valent pas au dixieme.
#
# POURQUOI LE NIVEAU EST LA BONNE CONVENTION ICI. La restriction de §2.4,
# |f_a| <= |f_z|, porte sur une MAGNITUDE de reponse, pas sur une derivee. Le
# rapport de niveaux est donc la mesure la plus fidele de ce que la restriction
# enonce, et il evite d'introduire une seconde convention pour 0,15 d'ecart.
#
# ⇒ LA PROSE ECRIT « exceeds ... by a factor of \EMPamenScaleExcess, measured
#   over the same range ». Elle N'ECRIT PAS « slope », et elle n'ajoute PAS
#   « per degC » : les deux transformeraient un rapport de niveaux en une
#   affirmation sur des pentes, que cette valeur ne porte pas. C'est la famille
#   du controle 10 -- un nom, ou un mot, qui decrit autre chose que ce que la
#   valeur contient.
##############################################################################
add('EMPamenScaleExcess', f'{_excess:.1f}', '10.3',   # lot 85o : 9.3 -> 10.3 (q50 MC Cholesky)
    'ratio of LEVEL DROPS over [5,10]->30 degC: f_a (MC q50) / f_z '
    '(damage_prod_bin). NOT a slope -- the prose must not write "slope".')

# ------------------------------------------------------- within-R2 amenity ---
# The daily and annual arms are reported separately: a single range mixing them
# The former regex spanned the bottom of the daily block and the top of the
# annual one, so it matched no table at all.
#
# HISTORY, KEPT BECAUSE IT NAMES THE FAILURE MODE. These four values were once
# emitted as [LITERAL]: no script of the package wrote the two files they came
# from, and what sat in output/ was a hand transcription from a console log.
# Reading a transcription would have given it the authority of a macro, which is
# circular rather than traceable. A MACRO WITH A FALSE PROVENANCE IS MORE
# DANGEROUS THAN AN ACKNOWLEDGED LITERAL: its form certifies a traceability that
# does not exist.
#
# CLOSED ON 09/09/2026, BOTH ARMS. tab_exogeneity_mrt_migration.tex and
# tab_exogeneity_mrt_annual_migration.tex are now both generated by
# 01b_amenity_drf.R, and their sixth column is the within R-squared per
# specification. Fixed effects are origin x destination and year -- the PAIR
# fixed effects -- on both arms.
def _within_r2(fname):
    """within-R2 par specification depuis une table de placebo pre-tendance.
    Retourne (liste des 12 valeurs, valeur de la spec de REFERENCE ns 4 df 0,5 %).
    Refuse de rendre autre chose que 12 lignes : un format qui change doit crier,
    pas emettre une etendue tronquee."""
    txt = read(os.path.join(OUT, fname))
    # `\$?` ajoute le 07/09 : la colonne winsorisation est desormais en mode
    # mathematique -- "& $[0.5\\%, 0.5\\%]$ &" et non "& [0.5\\%". Optionnel,
    # donc les deux conventions passent. Mesure : 0 ligne avant, 12 apres.
    rows = re.findall(r'^(Polynomial|Natural spline)\s*&\s*(\d+)\s*&\s*\$?\[([\d.]+)\\%.*?&'
                      r'\s*[\d.]+\s*&\s*[\d.]+\s*&\s*([\d.]+)\s*&', txt, re.M)
    if len(rows) != 12:
        raise SystemExit(f"{fname}: {len(rows)} rows parsed, 12 expected — "
                         "the format has changed; not emitting silently.")
    ref = [float(r[3]) for r in rows
           if r[0] == 'Natural spline' and r[1] == '4' and r[2] == '0.5']
    if len(ref) != 1:
        raise SystemExit(f"{fname}: reference spec (ns 4 df, 0.5 %) not found "
                         f"or ambiguous ({len(ref)} rows).")
    return [float(r[3]) for r in rows], ref[0]


# regen=True, meme motif que EMPmigWaldF : millesime v8 d'Abel & Cohen. Ces deux
# bornes ne portent aucune affirmation de seuil -- la phrase qui les cite dit que
# le within-R2 est ELEVE, ce que 0,729--0,733 dit aussi bien que 0,748--0,750.
_wr, _ = _within_r2('tab_exogeneity_mrt_migration.tex')
add('EMPamenWithinRdailyLo', f"{min(_wr):.3f}", '0.748',
    'tab_exogeneity_mrt_migration.tex (within-R2 min over 12 specs)', regen=True)
add('EMPamenWithinRdailyHi', f"{max(_wr):.3f}", '0.750',
    'tab_exogeneity_mrt_migration.tex (within-R2 max over 12 specs)', regen=True)

# EMPidRsqPair is deliberately not emitted. It had a real source, but its only
# consumer was a sentence whose second term is not reproducible. Emitting a macro
# with no consumer is the mirror image of a hard-coded number: the orphan
# detector would count it as unused. The format check on the reference
# specification stays active in _within_r2().

# The annual table exists on disk and gives the interval reported in the
# manuscript, so this value has a machine source.
def _placebo_p(fname):
    """p-values du test de leads, une table de placebo pre-tendance.
    Meme discipline que _within_r2 : refuse de rendre autre chose que 12 lignes."""
    txt = read(os.path.join(OUT, fname))
    # `\$?` ajoute le 07/09, meme motif que _within_r2 : la winsorisation est
    # passee en mode mathematique. Optionnel, donc les deux conventions passent.
    rows = re.findall(r'^(?:Polynomial|Natural spline)\s*&\s*\d+\s*&\s*\$?\[[\d.]+'
                      r'.*?&\s*[\d.]+\s*&\s*([\d.]+)\s*&', txt, re.M)
    if len(rows) != 12:
        raise SystemExit(f"{fname}: {len(rows)} rows parsed, 12 expected — "
                         "the format has changed; not emitting silently.")
    return [float(x) for x in rows]


# Controle prealable, pose le 07/09 : dire CE QUI manque et QUI devrait le
# produire, plutot que de laisser tomber un FileNotFoundError sur un chemin.
_PLAC_TABLES = {
    'tab_exogeneity_mrt_migration.tex':        '01b_amenity_drf.R',
    # 09/09/2026 : le bras ANNUEL a desormais un producteur. emit_placebo_table()
    # est UNE fonction parametree par is_annual, appelee deux fois ; les deux
    # appels ne different que par cet argument, les controles et le signe etant
    # dans le corps. La limite 4 du README est levee.
    'tab_exogeneity_mrt_annual_migration.tex': '01b_amenity_drf.R',
}
_plac_abs = [f for f in _PLAC_TABLES if not os.path.exists(os.path.join(OUT, f))]
if _plac_abs:
    raise SystemExit(
        "[STOP] placebo table(s) missing from output/:\n"
        + "\n".join(f"  {f}  <- {_PLAC_TABLES[f]}" for f in _plac_abs)
        + "\n  EMPplaceboP{min,max} bear on 24 cells (12 daily +\n"
          "  12 annual). Computing them on the 12 available would change a\n"
          "  published macro under cover of a repair: refused.")
# ----------------------------------------------------------------------
# LE MINIMUM PORTE SUR LA TABLE QUE LA PHRASE MONTRE
# ----------------------------------------------------------------------
# Defaut ferme le 17/09/2026. \EMPplaceboPmin valait 0.212, minimum sur
# VINGT-QUATRE cellules : les douze de tab_exogeneity_mrt_migration ET les
# douze de son jumeau annuel. Or le manuscrit n'\input QUE la premiere
# (D.2), et D.2 ecrit, a cote d'elle : « Table 11 substitutes one-period-ahead
# differentials [...] the smallest p-value on the leads is \EMPplaceboPmin ».
# Le lecteur qui verifiait trouvait 0.5553 dans la colonne, et une phrase
# disant 0.212 juste au-dessus. La table annuelle existe, elle est produite,
# mais elle n'est PAS montree : son minimum ne peut pas etre celui qu'une
# phrase attribue a une table visible.
#
# Le bras annuel n'est pas perdu : il part dans une macro qui le nomme. Le
# jour ou le manuscrit \input la table annuelle, la phrase aura son chiffre.
_plac_daily  = _placebo_p('tab_exogeneity_mrt_migration.tex')
_plac_annual = _placebo_p('tab_exogeneity_mrt_annual_migration.tex')
_plac = _plac_daily + _plac_annual
add('EMPplaceboPminAnnual', f"{min(_plac_annual):.4f}", '0.2118',
    'tab_exogeneity_mrt_annual_migration.tex (min p of the leads, annual arm, '
    'table NOT shown by the manuscript)', regen=True)
add('EMPplaceboPminBoth', f"{min(_plac):.3f}", '0.212',
    'tab_exogeneity_mrt{,_annual}_migration.tex (min p over the 24 cells '
    'of BOTH arms)', regen=True)
# regen=True. Le minimum passe de 0,064 a 0,212 : la phrase qui le cite affirme
# que le placebo n'est PAS rejete au seuil de 5 %, et 0,212 le dit plus largement
# que 0,064, qui en etait a un cheveu. Le sens de la phrase est renforce, pas
# retourne -- c'est la condition pour qu'un regen soit une reparation.
add('EMPplaceboPmin', f"{min(_plac_daily):.4f}", '0.5553',
    'tab_exogeneity_mrt_migration.tex (min p of the leads over the 12 cells of '
    'THE table that D.2 shows; same data frame as its note)', regen=True)

# Les jumelles annuelles des deux bornes ci-dessus, produites par la table que
# 01b emet depuis le 09/09. Meme motif, meme absence d'affirmation de seuil.
_wra, _ = _within_r2('tab_exogeneity_mrt_annual_migration.tex')
add('EMPamenWithinRannualLo', f"{min(_wra):.3f}", '0.746',
    'tab_exogeneity_mrt_annual_migration.tex (within-R2 min over 12 specs)', regen=True)
add('EMPamenWithinRannualHi', f"{max(_wra):.3f}", '0.751',
    'tab_exogeneity_mrt_annual_migration.tex (within-R2 max over 12 specs)', regen=True)

# ---- ANNUAL joint F ----------------------------------------------------------
# Read from amenity_joint_tests_annual.csv.
# The recomputed bounds differ from the value the manuscript carried at the upper
# end. `regen=True` means the source governs and the gap to the expected value is
# REPORTED rather than blocking. The sentence claims indistinguishability from
# zero, which the p-values carry, so the upper bound crossing one changes no
# argument.
_jta = pd.read_csv(os.path.join(OUT, 'amenity_joint_tests_annual.csv'))
if len(_jta) != 12:
    sys.exit(f"[amenite] amenity_joint_tests_annual.csv: {len(_jta)} rows, 12 expected.")
add('EMPamenAnnualFjointLo', f"{_jta.fval_joint.min():.2f}", '0.30',
    'amenity_joint_tests_annual.csv (min joint F over 12 specs)', regen=True)
add('EMPamenAnnualFjointHi', f"{_jta.fval_joint.max():.2f}", '0.98',
    'amenity_joint_tests_annual.csv (max joint F over 12 specs)', regen=True)
add('EMPamenAnnualPjointLo', f"{_jta.pval_joint.min():.3f}", '0.138',
    'amenity_joint_tests_annual.csv (min joint p over 12 specs)')
add('EMPamenAnnualPjointHi', f"{_jta.pval_joint.max():.3f}", '0.413',
    'same, max')
# \EMPamenAnnualPnonlin{Lo,Hi} (lot 74) : le test des seuls termes de degre >= 2 du
# bras annuel (colonne pval_nonlin), pas le test joint de tous les termes ci-dessus.
# La phrase du supplement (amenite interieure) dit « no non-linearity the data
# support » : garde de phrase, le minimum doit rester au-dessus de 5 %.
add('EMPamenAnnualPnonlinLo', f"{_jta.pval_nonlin.min():.3f}", '0.077',
    'amenity_joint_tests_annual.csv (min non-lin. p over the 12 annual specs)')
add('EMPamenAnnualPnonlinHi', f"{_jta.pval_nonlin.max():.3f}", '0.648',
    'same, max')
import gardes_lot74 as _g74   # noqa: E402
_ok, _det = _g74.garde_pnonlin_annuel(float(_jta.pval_nonlin.min()))
garde_phrase(_ok, 'Supplement (interior amenity): "no non-linearity the data support"',
             '> 0.05', _det,
             "[STOP] annual amenity: %s. The sentence \"carries no non-linearity the "
             "data support\" no longer holds." % _det)

# Minimum of the NON-LINEARITY p-value over BOTH arms, across all cells. It
# replaces a threshold count, which covered the joint test only while the table
# also displays the non-linearity column, whose weakest cell sits well above the
# conventional level. A minimum names its object and lets the reader apply a
# threshold.
add('EMPamenPnonlinMin', f"{min(_jt.pval_nonlin.min(), _jta.pval_nonlin.min()):.3f}",
    '0.027', 'amenity_joint_tests{,_annual}.csv (min non-lin. p over the 24 cells)')

# ---- sensibilite a la ponderation de tab:dbar_cross : RETIREE (11/09/2026) --
# Les quatre macros EMPdbarMADequal / RMSequal / MADpop / RMSpop chiffraient la
# divulgation de la section F. La section F est retiree du manuscrit et
# tab_dbar_cross avec elle ; les quatre macros n'ont plus aucun appel (mesure :
# 0 dans JMP_V11.tex).
#
# La garde SystemExit qui exigeait dbar_cross_weighting.csv est retiree AVEC
# elles : elle aurait arrete l'emetteur entier pour un artefact que plus
# personne ne lit. Une garde qui protege un consommateur disparu ne protege
# rien, elle bloque.
#
# Le CALCUL survit dans 01e_diagnostics.R section 10, ou il porte le controle
# positif qui compare le bras equipondere a `res`. Si la divulgation revient au
# manuscrit, les macros se re-emettent d'ici, depuis le meme CSV.

# ---------------------------- absorption FE / controles (annexe identification) ----
# The four absorption macros read CALCULATED values: the producer CSV when the
# estimation has been re-run, otherwise the BODY of the tables, never their note.
# A table note is a string typed by hand in the producing script, and it can
# disagree with the column right above it.
#
# The names carry which end of the range they hold (Min, Max). Without that mark,
# a macro holding a MAXIMUM can be used in a sentence that says "more than",
# which is a falsehood that correcting the number alone does not repair.

def _absorption():
    """Renvoie (r2_pair[], r2_pyr[], sd_within[], absorbed[]) et la source utilisee.
    Priorite au CSV producteur ; repli sur le CORPS des tables .tex.
    Rapporte TOUJOURS le chemin emprunte (contrat §0, quatrieme regle : un
    mecanisme qui peut ne rien faire doit dire s'il a agi)."""
    fe_csv, ct_csv = (os.path.join(OUT, f) for f in ('fe_absorption.csv',
                                                     'ctrl_absorption.csv'))
    if os.path.exists(fe_csv) and os.path.exists(ct_csv):
        fe, ct = pd.read_csv(fe_csv), pd.read_csv(ct_csv)
        return (fe.r2_pair_raw.tolist(), fe.r2_pyr_raw.tolist(),
                fe.sd_within_ratio.tolist(), ct.r2_absorbed_raw.tolist(),
                'fe_absorption.csv / ctrl_absorption.csv (producer CSV)')

    row = r'\$\\Delta_\d\$ & ([\d.]+) & ([\d.]+) & ([\d.]+)'
    fe = re.findall(row, read(os.path.join(OUT, 'tab_fe_absorption.tex')))
    ct = re.findall(row, read(os.path.join(OUT, 'tab_ctrl_absorption.tex')))
    assert fe and ct, 'table body not found — has the format of the .tex changed?'
    r2_pair = [float(a) for a, _, _ in fe]
    r2_pyr  = [float(b) for _, b, _ in fe]
    # sd_within re-derived from an r2 rounded to four decimals: VALID FOR THE
    # MAXIMUM ONLY. At that rounding the smallest term is undetermined within a
    # tenth of a point, so its first decimal is not recoverable here. The
    # producer CSV above removes this limit.
    sd_within = [math.sqrt(1 - v) for v in r2_pair]
    absorbed  = [float(c) for _, _, c in ct]
    return (r2_pair, r2_pyr, sd_within, absorbed,
            'BODY of tab_fe_absorption.tex / tab_ctrl_absorption.tex '
            '(fallback: producer CSV missing, drf_prod_final.r not re-run)')

_r2pair, _r2pyr, _sdw, _absorbed, _absorc = _absorption()
print(f'[absorption] source = {_absorc}')

add('EMPfeAbsorptionMin', f'{100 * min(_r2pair):.1f}', '97.7', _absorc)
add('EMPfeAbsorptionMax', f'{100 * max(_r2pair):.1f}', '99.7', _absorc)
# REBASAGE DU 10/09/2026 : 0.44 -> 0.42. Brut mesure max(r2_pyr_raw) =
# 0.00421357 (fe_absorption.csv, 01a_productivity_drf.R l.528, ecrit le 09/09
# 23:27). Cause : le panel commerce de 01a a change d'echantillon au correctif
# GLORIA -- meme rejeu que EMPtradeNobs, 3 730 958 contre 3 806 447.
#
# 🔑 POURQUOI CELLE-CI SEULE A BOUGE dans une famille de six. Les cinq soeurs
# (feAbsorptionMin/Max 97.8/99.7, withinPairSDMax 14.9, ctrlAbsorpLo/Hi 6.4/14.8)
# sont sorties en MATCH au meme rejeu. Elles ne sont pas plus stables : elles
# sont arrondies plus grossierement. EMPfePyrMax est la seule de la famille a
# deux decimales sur une grandeur de l'ordre de 0.4 %, donc la seule dont la
# resolution attrape un deplacement d'echantillon de cette taille.
#
# ⚠️ ZERO APPELANT dans JMP_V10.tex (mesure le 10/09), et zero dans les trois
# sauvegardes : contrairement a EMPosterPolyTwoBtilde, cette macro n'a JAMAIS eu
# de phrase consommatrice. Elle est rebasee et non retiree parce que le retrait
# n'a pas ete arbitre ; le rebasage ne le prejuge pas.
add('EMPfePyrMax',        f'{100 * max(_r2pyr):.2f}',  '0.41', _absorc)
# MAX only: it is the binding end of the argument, that within-pair variation is
# small, and the only one robust to the rounding in fallback mode.
add('EMPwithinPairSDMax', f'{100 * max(_sdw):.1f}',    '15.2', _absorc)
add('EMPctrlAbsorpLo',    f'{100 * min(_absorbed):.1f}', '4.9', _absorc)
add('EMPctrlAbsorpHi',    f'{100 * max(_absorbed):.1f}', '14.6', _absorc)

# --- DEPRECATED ALIASES --------------------------------------------------------
# Two earlier macro names are still called by the manuscript. Removing them would
# add undefined control sequences, so they are emitted, and emitted on the
# BINDING end of their range, which makes the existing sentences true as written.
# Do not reuse them in a new sentence: their names do not say which end they
# carry, and that is the defect they came from.
add('EMPfeAbsorption', f'{100 * min(_r2pair):.1f}', '99.7',
    _absorc + ' — DEPRECATED ALIAS (min, not max: the sentence says "more than")',
    regen=True)
add('EMPwithinPairSD', f'{100 * max(_sdw):.1f}', '5.6',
    _absorc + ' — DEPRECATED ALIAS (max; the former 5.6 was Delta_1 silently generalised '
    'to three terms ranging from 5.5 to 14.9)',
    regen=True)
# Four macros of the cross-sectional identification passage are deliberately NOT
# emitted, and should not be reinstated.
#
# The reason, so that no one restores them:
#   - the specification that produced them, with SEPARATE origin and destination
#     fixed effects, exists in no live script; the live one is a PAIR fixed
#     effect, which is a different object;
#   - that cross-sectional design was removed because it FAILS the pre-trend
#     placebo, which the manuscript did not report;
#   - their values appear in no tracked artefact. Their only home was a prose
#     appendix written by hand, so reading them would have been circular.
#
# This was therefore not an execution debt: no re-run could close it. An
# unverifiable robustness argument is worse than a stale number, because its form
# claims a check that never happened.
#
# Do not reintroduce them unless the cross-sectional specification is recoded AND
# its placebo failure is reported in the text.
add('EMPdispThirty',    '8',    '8',    'LITERAL — cross-country dispersion 2030 (percentage points); '
    'to be emitted from damage_projection.r', literal=True)
add('EMPdispHundred',   '10',   '10',   'LITERAL — cross-country dispersion 2100; '
    'to be emitted from damage_projection.r', literal=True)

# NOTE on the hot slope: NO macro is emitted. The available artefact gives a
# value at 30 C under a max-zero normalisation, which is a DIFFERENT object from
# a drop relative to the cold plateau. No machine source carries the second, so
# it is a structural mismatch and not a rounding one. It must be resolved before
# any macro is created for it.
# NOTE on the external comparator: a LITERATURE value with a citation, so it
# stays hard-coded.

# =============================================================================
# CONTRIBUTION A -- 30 macros of the specification grid. FROZEN, not active.
#
# The block below is written and ready, and it sits here because this is where it
# would run: after every add() and before the report. Uncommenting it is an
# execution, not an edit.
#
# Commented code that is WAITING looks exactly like commented code that was
# ABANDONED. This one is waiting.
#
# Why it is frozen, and the reason is not technical:
#   - its sources live in a working directory that is neither under output/ nor
#     tracked. Pointing the emitter at that directory is the one-line fix, and it
#     is the wrong one: thirty PUBLISHED macros would depend on an unversioned
#     directory. The freeze holds until the producers are versioned.
#   - the text that consumes these macros is not inserted in the manuscript.
#     Emitting macros no sentence calls creates the orphan the rest of this file
#     works to avoid.
#
# On the pairing of `expected` with `regen=True`: these macros have no hard-coded
# history in the manuscript, so passing the measured value as `expected` with
# regen=False would compare a value to itself, which checks nothing. With
# regen=True the same `expected` becomes a dated BASELINE: the source governs,
# and the gap is reported rather than blocking.
#
# NO entry in the read-cache registry: that registry covers the .rds caches
# written by one emitter and never rewritten by a second producer. The problem it
# addresses does not arise here, and adding an entry would be a decorative guard.
#
# What the sentence asserts must be what the macro carries. A macro whose name
# does not name its object lets the two drift apart, which is why the counts
# below are emitted as DISTINCT quantities rather than one count reused: a
# unanimity that broke must be visible in the rendered text.
#
# ---- PARTIAL UNFREEZE: only the macros the text actually cites ----------------
# The other macros of this block stay COMMENTED below, each with its reason.
# Emitting them would create orphans, which is the defect that cost the removal
# of an earlier block of macros, and which this file warns about elsewhere.
#
# Signed zero: the bounds print with an explicit sign. Both sit far enough above
# the rounding threshold that no "-0.000" can appear, and the guardrail below
# checks it on every emission rather than assuming it.
CONTRIB_A = os.path.join(OUT, 'contribA')
SRC_T1 = 'contribA/test1_delta_north.csv (04_contrast/ct1_delta_north_grid.R)'
SRC_DB = 'contribA/tache1_diff_bras.csv (04_contrast/ct2_between_arm_test.R)'

_g = pd.read_csv(os.path.join(CONTRIB_A, 'test1_delta_north.csv'))
_g = _g[_g.channel == 'prod']
_ga = _g[_g.arm == 'annual']
_gd = _g[_g.arm == 'daily']


def _nz(x, dec=3):
    """Un zero signe est un signe sans magnitude : sous le dernier chiffre
    affichable on rend un zero NON signe."""
    return ('0.' + '0' * dec) if abs(x) < 0.5 * 10 ** (-dec) else f"{x:+.{dec}f}"


add('EMPgridProdNspec',   _gc('EMPgridProdNspec'),   '15',     SRC_GRIDC, regen=True)
add('EMPgridAnnualExcl',  _gc('EMPgridAnnualExcl'),  '12',     SRC_GRIDC, regen=True)
add('EMPgridDailyExcl',   _gc('EMPgridDailyExcl'),   '0',      SRC_GRIDC, regen=True)
add('EMPgridAnnualMin',   _nz(_ga.delta_north.min()),            '+0.073', SRC_T1, regen=True)
add('EMPgridAnnualMax',   _nz(_ga.delta_north.max()),            '+0.154', SRC_T1, regen=True)
add('EMPgridDailyMin',    _nz(_gd.delta_north.min()),            '-0.008', SRC_T1, regen=True)
add('EMPgridDailyMax',    _nz(_gd.delta_north.max()),            '+0.015', SRC_T1, regen=True)

# ---- the four macros of section D -------------------------------------------
# They carry what distinguishes "no signal" from "a weak signal not detected":
# the points are CENTRAL (median, count of negatives) and the intervals are WIDE
# (an envelope). The envelope must never be presented as a tight bound: it is of
# the order of a tenth of a log point, and a tight bound would make the paper say
# that the North does not move.
add('EMPgridDailyNneg',    _gc('EMPgridDailyNneg'),  '2',      SRC_GRIDC, regen=True)
add('EMPgridDailyMedian',  _nz(_gd.delta_north.median()),         '+0.006', SRC_T1, regen=True)
add('EMPgridDailyCIloMin', _nz(_gd.ci_lo.min(), 4),               '-0.0700', SRC_T1, regen=True)
add('EMPgridDailyCIhiMax', _nz(_gd.ci_hi.max(), 4),               '+0.0934', SRC_T1, regen=True)

# guardrail: no signed zero among the values just emitted
for _n, _v in (('AnnualMin', _ga.delta_north.min()), ('AnnualMax', _ga.delta_north.max()),
               ('DailyMin', _gd.delta_north.min()),  ('DailyMax', _gd.delta_north.max()),
               ('DailyMedian', _gd.delta_north.median())):
    assert not _nz(_v).startswith('-0.000'), f"[FATAL] signed zero on {_n}"
# add('EMPgridDailyAbsMax', f"{_gd.delta_north.abs().max():.3f}",  '0.036',  SRC_T1, regen=True)
# add('EMPgridSeAnnualLo',  f"{_ga.se.min():.3f}",                 '0.034',  SRC_T1, regen=True)
# add('EMPgridSeAnnualHi',  f"{_ga.se.max():.3f}",                 '0.053',  SRC_T1, regen=True)
# add('EMPgridSeDailyLo',   f"{_gd.se.min():.3f}",                 '0.024',  SRC_T1, regen=True)
# add('EMPgridSeDailyHi',   f"{_gd.se.max():.3f}",                 '0.052',  SRC_T1, regen=True)
#
# # hors-support : masse-jour nordique ecretee, par winsorisation (bras daily)
# _oos = _gd.groupby('winsor_lo').oos_mass_fut.first().sort_index()
# add('EMPoosLo',  f"{100*_oos.loc[0.5]:.2f}", '0.51', SRC_T1, regen=True)
# add('EMPoosMid', f"{100*_oos.loc[1.0]:.2f}", '1.18', SRC_T1, regen=True)
# add('EMPoosHi',  f"{100*_oos.loc[2.5]:.2f}", '3.62', SRC_T1, regen=True)
# # deplacement de delta a travers les 3 winsorisations : MAX par (forme, df),
# # pas la mediane. « au plus 0,012 » est une borne ; la mediane (0,009) n'en
# # est pas une, et l'ecrire sans le dire laisse croire au contraire.
# _sp = _gd.groupby(['form', 'df']).delta_north.agg(lambda s: s.max() - s.min())
# add('EMPoosShift', f"{_sp.max():.3f}", '0.012', SRC_T1, regen=True)
#
_d = pd.read_csv(os.path.join(CONTRIB_A, 'tache1_diff_bras.csv'))
_d = _d[_d.channel == 'prod']
# TWO counts, and they DIFFER: one with the covariance and one without. Emitting
#    them separately is the substance of the argument. A single claim that "the
#    sign holds under both" would be a disguised unanimity.
add('EMPdiffExcl',      _gc('EMPdiffExcl'),      '15',   SRC_GRIDC, regen=True)
add('EMPdiffExclNaive', _gc('EMPdiffExclNaive'), '12',   SRC_GRIDC, regen=True)
add('EMPdiffRhoLo',     f"{_d.rho.min():.2f}",                               '0.49', SRC_DB, regen=True)
add('EMPdiffRhoHi',     f"{_d.rho.max():.2f}",                               '0.74', SRC_DB, regen=True)
garde_phrase(_d.rho.min() > 0, 'note: every rho positive', '> 0',
             'rho min %.4f' % _d.rho.min(),
             "[FATAL] a negative rho: the meaning of the note is reversed")

# ---------------------------------------------------------------------------
# LA GRILLE LARGE : ARCHIVEE, PLUS CITEE, DONC PLUS EMISE (18/09/2026)
# ---------------------------------------------------------------------------
# Trois macros ont vecu ici quelques heures -- \EMPgridWideNspec,
# \EMPdiffExclWide, \EMPdiffExclNaiveWide -- pour une phrase du supplement sur
# la grille a vingt-quatre. Decision d'auteur du meme jour : la phrase est
# retiree, les macros avec elle.
#
# 🔑 LE MOTIF EST CELUI QUI AVAIT ETE SIGNALE EN LES ECRIVANT. Elles etaient les
# SEULES macros du manuscrit dont la source ne soit pas une etape de la chaine
# mais un artefact archive, output/contribA_grid_wide.csv, qu'aucun producteur
# ne sait plus emettre depuis que ct1, ct2, ct4 et ct6 lisent GRID_SPECS. Un
# chiffre publie dont la source ne se rejoue pas est un chiffre qu'on ne peut
# pas defendre deux fois.
#
# ⛔ NE PAS LES REMETTRE SANS REMETTRE UN PRODUCTEUR. L'archive reste sur
# disque (avec contribA_grid_wide.PROVENANCE.txt) comme trace de la grille
# large ; elle n'a plus de consommateur, et c'est l'etat voulu. Les faire
# revivre demande de reposer GRID_FORMS a huit formes et de rejouer ct1 et ct2,
# pas de relire le CSV.
#
# # ANNEXE seulement. ⛔ NE PAS emettre de macro pour l'ECART entre les deux
# #    optima ni pour leur rapport de dispersion : le contrat interdit de
# #    chiffrer la comparaison, et un rapport EST la comparaison sous un autre
# #    nom. Les deux ecarts-types sont emis separement, et c'est tout.
# _a = pd.read_csv(os.path.join(CONTRIB_A, 'tache3_argmax_boot.csv')).set_index('arm')
# _w = pd.read_csv(os.path.join(CONTRIB_A, 'tache3_argmax_boot_draws.csv'))
# _w = _w[_w.ok.astype(str).str.upper().isin(('TRUE', '1'))]
# add('EMPargmaxAnnual',    f"{_a.loc['annual', 'argmax_point']:.1f}",          '11.9', SRC_AX, regen=True)
# add('EMPargmaxDaily',     f"{_a.loc['daily',  'argmax_point']:.1f}",          '7.8',  SRC_AX, regen=True)
# add('EMPargmaxB',         f"{int(_a.loc['daily', 'B'])}",                     '200',  SRC_AX, regen=True)
# add('EMPargmaxBoundFrac', f"{int(_a.frac_draws_at_bound.max() * _a.B.max())}", '0',   SRC_AX, regen=True)
# add('EMPargmaxSdAnnual',  f"{_w[_w.arm == 'annual'].argmax.std(ddof=0):.2f}", '0.76', SRC_AX, regen=True)
# add('EMPargmaxSdDaily',   f"{_w[_w.arm == 'daily' ].argmax.std(ddof=0):.2f}", '1.78', SRC_AX, regen=True)
#
# # CRH : AUCUN fichier a verser — crh_inference_deliverable.csv est deja lu par
# # cet emetteur (CRH_CSV, l. 159) et sa copie output/ est SUIVIE par git.
# # Verifie le 07/08 : md5 a1f3740d identique a _work_crh_inference/.
# _c = _crh_all[(_crh_all.channel == 'prod')
#               & (_crh_all.publish.str.startswith('PUBLISH'))
#               & (_crh_all.vcov_variant == 'country')]
# _cnt = lambda o: int((_c[_c.object == o].ci95_lo > 0).sum())
# add('EMPcrhUnits',      f"{_c[_c.object == 'level_mean_2095'].unit.nunique()}", '6', CRH_CSV, regen=True)
# add('EMPcrhMeanCount',  f"{_cnt('level_mean_2095')}",                          '6', CRH_CSV, regen=True)
# add('EMPcrhDailyCount', f"{_cnt('level_daily_2095')}",                         '1', CRH_CSV, regen=True)
# add('EMPcrhDiffCount',  f"{_cnt('difference_mean_minus_daily')}",              '6', CRH_CSV, regen=True)
# =============================================================================

# ========================= FINLAND -- Section 5 illustration =================
# The nine literals of the Finland paragraph. They were verified but not emitted
# until the columns carrying them existed.
#
# The historical arm and the projected arm give values that agree to the rounding
# on one figure and differ on another. Reading the wrong arm would therefore be
# invisible on the first and visible on the second. Hence the check on the base
# column itself, never on the file name.
#
# No regen=True here. These nine values are EXISTING literals of the manuscript,
# and a macro must emit exactly what was hard-coded. A gap must BLOCK the
# emission and be seen, not be reported as a new baseline.
_FIN_DP = os.path.join(OUT, 'damage_projection_country.csv')
_FIN_PE = os.path.join(OUT, 'pe_signtest_by_country.csv')
_FIN_BASE, _FIN_PROJ = 2015, 2095
_fin_skip = None

if not (os.path.exists(_FIN_DP) and os.path.exists(_FIN_PE)):
    _fin_skip = 'artefact(s) missing'
else:
    _dp = pd.read_csv(_FIN_DP)
    _pe = pd.read_csv(_FIN_PE)
    _miss = ([c for c in ('days_ge20', 'days_ge25', 'days_ltm10') if c not in _dp.columns]
             + [c for c in ('base_source_cmip', 'daily_change_cmip', 'annual_change_cmip')
                if c not in _pe.columns])
    if _miss:
        _fin_skip = 'missing columns %s -- re-run 1 and/or 2 not done' % _miss
    else:
        _dp = _dp[_dp.GID_0 == 'FIN'].set_index('year')
        _pe = _pe[(_pe.GID_0 == 'FIN') & (_pe.horizon == _FIN_PROJ)]
        if _FIN_BASE not in _dp.index or _FIN_PROJ not in _dp.index:
            _fin_skip = 'year %d or %d missing for FIN' % (_FIN_BASE, _FIN_PROJ)
        elif len(_pe) == 0:
            _fin_skip = 'no FIN row at horizon %d' % _FIN_PROJ
        elif set(_pe.base_source_cmip.unique()) != {_FIN_BASE}:
            # The base is READ; it is not inferred from the producing script.
            _fin_skip = 'CMIP6 base = %s, %d expected' % (
                sorted(_pe.base_source_cmip.unique()), _FIN_BASE)
        elif _pe.annual_change_cmip.round(9).nunique() != 1:
            # The ANNUAL arm does not depend on the daily specification. If it
            # varies, it is not the quantity the text describes.
            _fin_skip = 'annual_change_cmip not invariant across specs (%d values)' % (
                _pe.annual_change_cmip.round(9).nunique())
        elif abs(float(_pe.aic_weight.sum()) - 1.0) > 1e-9:
            _fin_skip = 'AIC weights sum to %.6f, 1 expected' % _pe.aic_weight.sum()
        elif _GRID_N is not None and len(_pe) != _GRID_N:
            # La phrase cite \EMPgridNspec comme denominateur de la fourchette
            # calculee ici. Les deux doivent porter sur la MEME grille, sans quoi
            # la phrase annonce douze specifications et en resume un autre nombre.
            _fin_skip = 'FIN grid has %d specs, %d announced by EMPgridNspec' % (
                len(_pe), _GRID_N)

if _fin_skip:
    print('[FIN] macros NOT emitted: %s' % _fin_skip)
else:
    _SRC_DP = ('damage_projection_country.csv (01e_diagnostics.R section 8, '
               'CMIP6 base %d)' % _FIN_BASE)
    _SRC_PE = ('pe_signtest_by_country.csv (01e_diagnostics.R section 8bis, '
               'CMIP6 arm %d)' % _FIN_BASE)
    _ann = float(_pe.annual_change_cmip.iloc[0])
    _day = float((_pe.aic_weight * _pe.daily_change_cmip).sum())
    print('[FIN] AIC weights: max=%.3f over %d specs (weighted daily arm)'
          % (_pe.aic_weight.max(), len(_pe)))
    add('EMPfinTmeanBase',   f"{_dp.loc[_FIN_BASE, 'T_mean']:.1f}",     '5.3',    _SRC_DP)
    add('EMPfinTmeanProj',   f"{_dp.loc[_FIN_PROJ, 'T_mean']:.1f}",     '7.2',    _SRC_DP)
    # REBASAGE DU 10/09/2026 : +0.215 -> +0.216. Brut mesure 0.21614795, un seul
    # chiffre de dernier rang. Cause : rejeu de 01e section 8bis sur le climat
    # re-emis. Le bras ANNUEL est invariant a travers les 12 specifications (le
    # garde-fou `annual_change_cmip non invariant` ci-dessus le verifie), donc ce
    # deplacement vient de l'exposition, pas d'un choix de specification.
    add('EMPfinGainAnnual',  f"{_ann:+.3f}",                            '+0.036', _SRC_PE)
    # ------------------------------------------------------------------
    # REBASAGE ARBITRE DU 10/09/2026 : +0.027 -> +0.049
    # ------------------------------------------------------------------
    # CE QUE LA MACRO MESURE : la somme ponderee par AIC de daily_change_cmip
    # sur la grille, pour FIN a l'horizon 2095, base CMIP6 2015.
    #
    # ⚠️ LE VECTEUR DE POIDS EST DEGENERE, et il faut le savoir pour lire
    # l'ecart. Mesure du 10/09 : 1.00 sur spline4_wl0.5wu0.5, 0.00 sur les onze
    # autres. La « moyenne ponderee » est donc la specification argmax-AIC
    # SEULE -- ce que le manuscrit dit deja de ses poids d'Akaike en deux
    # endroits (note de la Figure DRF, legende de fig:oh_asym_fe).
    #
    # L'ECART N'EST PAS UNE DERIVE DE DONNEES, C'EST UN BASCULEMENT D'ARGMAX.
    # Les douze specifications s'etalent de 0.025619 a 0.048657 pour FIN :
    # l'ancien +0.027 etait au BAS de cette grille, le nouveau +0.048657 en est
    # le HAUT exact. Le bras ANNUEL, lui, ne bouge que d'un dernier chiffre
    # (+0.215 -> +0.216) et reste invariant a travers les douze. Ce n'est donc
    # pas l'exposition qui a double le bras journalier : c'est la forme
    # fonctionnelle retenue qui a change.
    #
    # CONSEQUENCE SUR L'ENONCE, ET C'EST LA RAISON DES TROIS MACROS CI-DESSOUS.
    # La phrase codait « about an eighth » en dur. 0.027/0.215 = 0.126 valait
    # bien un huitieme ; 0.049/0.216 = 0.225 n'en vaut plus un. Et la fraction
    # publiee n'avait JAMAIS ete invariante a la grille : le rapport y va de
    # 11.9 % a 22.5 %. Une fraction en dur posee sur une valeur selectionnee par
    # AIC redevient fausse au prochain basculement, en silence. Elle est donc
    # remplacee par des macros, et la fourchette est divulguee.
    _ratio = _day / _ann
    _r_lo  = float(_pe.daily_change_cmip.min()) / _ann
    _r_hi  = float(_pe.daily_change_cmip.max()) / _ann

    # ---- gardes Finlande (a) et (b) : RETIREES (lot 74, 01/10/2026) ----
    # (a) gardait « the Akaike-selected form is the one that makes the gap
    # smallest », (b) « the daily arm stays a fraction of the annual-mean one ».
    # Aucune des deux phrases n'existe plus dans JMP_V11.tex, ni dans le texte
    # vivant ni en commentaire (releve du lot 73). Une garde de phrase sans phrase
    # ne protege rien et rapporte un ecart qui ne decrit aucun texte. Si l'une des
    # deux phrases revient, sa garde revient avec elle (.bak_20261001_lot74).

    add('EMPfinGainDaily',   f"{_day:+.3f}",                            '+0.009', _SRC_PE)
    ##########################################################################
    # \EMPfinGainFactor -- LE FACTEUR ENTRE LES DEUX GAINS FINLANDAIS.
    # Posee le 22/09/2026, et elle nait d'un piege evite de justesse.
    #
    # LA PHRASE de §3.3 dit desormais que les deux fonctions laissent la
    # Finlande gagnante et « disagree by a factor of ... on how much ». La
    # macro qu'on allait y mettre etait \EMPfinGainRatio -- et c'est FAUX :
    # elle est declaree `range` dans values_pairs.tsv, c'est le MIN/MAX du
    # rapport SUR LA GRILLE (22.5, egal a la borne haute par construction), et
    # son propre commentaire la dit « NON CITEE ». Le lecteur qui divise les
    # deux nombres imprimes a cote trouve ~4.7, pas 22.5.
    #
    # ⇒ ON EMET LE QUOTIENT DES DEUX GRANDEURS QUE LA PHRASE IMPRIME, calcule
    #   AVANT arrondi -- la regle de la note de convention de §2.3. C'est le
    #   seul objet dont la phrase parle.
    ##########################################################################
    if _day == 0:
        sys.exit("[STOP] Finland: zero daily gain, the factor diverges.")
    # ACQUITTE au premier passage : l'epingle 4.7 etait une PREDICTION, faite
    # en divisant les deux nombres ARRONDIS (+0.033 / +0.007). La mesure sur
    # les valeurs non arrondies rend 4.4 -- et c'est precisement le phenomene
    # que la note de convention de §2.3 annonce. Cinquieme rapport de cette
    # famille, et le seul cree aujourd'hui.
    add('EMPfinGainFactor', f"{_ann / _day:.1f}", '4.2',
        _SRC_PE + ': EMPfinGainAnnual / EMPfinGainDaily, before rounding')
    # `expected` = la valeur MESUREE au run du 10/09 : ces trois macros sont
    # NEUVES, aucune n'a jamais eu de litteral dans le manuscrit. Le garde-fou
    # est donc une non-regression, pas une confrontation a un hard-code.
    ##########################################################################
    # `Min`/`Max`, PAS `Lo`/`Hi`, ET C'EST LE CORRECTIF DU 20/09/2026.
    # _r_lo et _r_hi sont le MIN et le MAX sur les douze specifications, pas
    # des bornes de confiance. Et le point estime EGALE LA BORNE HAUTE PAR
    # CONSTRUCTION : la garde ci-dessus exige que le bras journalier retenu
    # par AIC soit le maximum de la grille. `22.5 [11.9 ; 22.5]` se lisait
    # donc comme un intervalle degenere qui n'existe pas.
    # ⛔ ET ELLE RESTE INEMPLOYEE : c'est une attenuation au niveau PAYS, et le
    # niveau pays n'est identifie d'aucun cote -- CRH marque la France
    # FRAGILE a 162,8 %, et le registre notait le 28/07 « fourchette pays
    # 44-79 % NON identifiee ». L'agregat, jamais les pays.
    ##########################################################################
    add('EMPfinGainRatio',   f"{100*_ratio:.1f}",                       '23.7', _SRC_PE)
    add('EMPfinGainRatioMin', f"{100*_r_lo:.1f}",                        '13.2', _SRC_PE)
    add('EMPfinGainRatioMax', f"{100*_r_hi:.1f}",                        '24.1', _SRC_PE)
    add('EMPfinDaysHotBase', f"{_dp.loc[_FIN_BASE, 'days_ge20']:.1f}",  '14.3',   _SRC_DP)
    add('EMPfinDaysHotProj', f"{_dp.loc[_FIN_PROJ, 'days_ge20']:.1f}",  '30.6',   _SRC_DP)
    add('EMPfinDaysVHotBase', f"{_dp.loc[_FIN_BASE, 'days_ge25']:.1f}", '0.5',    _SRC_DP)
    add('EMPfinDaysVHotProj', f"{_dp.loc[_FIN_PROJ, 'days_ge25']:.1f}", '2.0',    _SRC_DP)
    add('EMPfinDaysColdProj', f"{_dp.loc[_FIN_PROJ, 'days_ltm10']:.0f}", '10',    _SRC_DP)

##################################################################################
# LE BIAIS DE BINNING (lot 78, 02/10/2026) -- sorties de ct8 et ct9.
#
# Jones, Moscona, Olken, von Dessauer (2026) : en panel a effets fixes de lieu,
# les comptes de jours extremes portent des tendances de rechauffement correlees
# au climat de base, et la regression peut fabriquer un bras froid. ct8
# re-estime les specifications de reference sous quatre structures d'effets fixes
# (base = la PPML publiee, between, pp, pp5) ; ct9 mesure les ingredients du biais
# et la part intra-cellule du design CRH (supplement E).
#
# Toutes les macros sont NEUVES : `expected` = la mesure du bac a sable du lot 77
# (out_v2/*.csv, out/*.csv), donc une non-regression, pas un hard-code.
# Gardes : G1 (coherence entre deux calculs : base == ct1/ct2) BLOQUANTE ; G2-G7,
# gardes de PHRASE, par garde_phrase() et tools/gardes_lot78.py (connu-faux :
# tools/test_gardes_lot78.py). Une garde de phrase qui refuse dit que la phrase
# est fausse : on la reecrit, on ne relache pas la garde.
##################################################################################
import gardes_lot78 as _g78   # noqa: E402
_DCB = os.path.join(OUT, 'contribA')
_BIN = {}
for _k, _f in (('V', 'binning_variants.csv'), ('C', 'binning_contrasts.csv'),
               ('K8', 'ct8_controls.csv'), ('M', 'binning_mechanism.csv'),
               ('W', 'crh_within_cell.csv'), ('K9', 'ct9_controls.csv')):
    _p = os.path.join(_DCB, _f)
    if not os.path.exists(_p):
        sys.exit("[STOP] binning: %s missing from output/contribA/\n"
                 "  <- 04_contrast/%s (step %s of run_all.sh)"
                 % (_f, 'ct8_binning_variants.R' if _k in ('V', 'C', 'K8') else
                    'ct9_binning_mechanism.R', 'ct8' if _k in ('V', 'C', 'K8') else 'ct9'))
    _BIN[_k] = pd.read_csv(_p)
_bV, _bC, _bM, _bW = _BIN['V'], _BIN['C'], _BIN['M'], _BIN['W']
_BVARS = ['base', 'between', 'pp', 'pp5']
_BFITS = {'daily_ref', 'annual_ref', 'annual_ge'}
## ---- structure (bloquant) ----
for _k in ('K8', 'K9'):
    if (_BIN[_k]['verdict'] != 'OK').any() or len(_BIN[_k]) == 0:
        sys.exit("[STOP] binning: controls of %s not all OK:\n%s"
                 % (_k, _BIN[_k].to_string()))
if sorted(_bV['variant'].unique()) != sorted(_BVARS) or \
        any(set(_bV.loc[_bV.variant == v, 'fit']) != _BFITS for v in _BVARS) or \
        len(_bV) != len(_BVARS) * len(_BFITS) or (_bV['status'] != 'OK').any():
    sys.exit("[STOP] binning_variants.csv: %d rows, variants %s; %d x %d expected, all OK."
             % (len(_bV), sorted(_bV['variant'].unique()), len(_BVARS), len(_BFITS)))
if len(_bC) != 2 * len(_BVARS) or set(_bC['contrast']) != {'ref_ct2', 'ref_ge'}:
    sys.exit("[STOP] binning_contrasts.csv: %d rows, %d expected." % (len(_bC), 2 * len(_BVARS)))
if _bV['inputs_empreinte'].nunique() != 1 or _bC['inputs_empreinte'].nunique() != 1 or \
        _bV['inputs_empreinte'].iloc[0] != _bC['inputs_empreinte'].iloc[0]:
    sys.exit("[STOP] binning: fits and contrasts do not carry the same input fingerprint.")
_bv = _bV.set_index(['variant', 'fit'])
_bc = _bC[_bC.contrast == 'ref_ct2'].set_index('variant')

## ---- G1 : base == ct1/ct2, a l'emission (BLOQUANT, coherence entre deux calculs) ----
_t1b = pd.read_csv(os.path.join(_DCB, 'test1_delta_north.csv'))
_t1b = _t1b[(_t1b.channel == 'prod') & (_t1b.status == 'OK')]
_t2b = pd.read_csv(os.path.join(_DCB, 'tache1_diff_bras.csv'))
_t2b = _t2b[(_t2b.channel == 'prod') & (_t2b.status == 'OK')]
_g1 = []
for _fit in sorted(_BFITS):
    _r = _bv.loc[('base', _fit)]
    _p = _t1b[(_t1b.arm == _r['arm']) & (_t1b.spec_label == _r['spec_label'])]
    if len(_p) != 1:
        sys.exit("[STOP] G1: %s / %s missing or multiple in test1_delta_north.csv"
                 % (_r['arm'], _r['spec_label']))
    _g1 += [abs(_r['delta_north'] - _p['delta_north'].iloc[0]), abs(_r['se'] - _p['se'].iloc[0])]
_rb = _bc.loc['base']
_pb = _t2b[_t2b.spec_label == _PARAMS.GRID_REF_SPEC]
if len(_pb) != 1:
    sys.exit("[STOP] G1: reference %s missing from tache1_diff_bras.csv" % _PARAMS.GRID_REF_SPEC)
_g1 += [abs(_rb['diff'] - _pb['diff'].iloc[0]), abs(_rb['se_diff'] - _pb['se_diff'].iloc[0])]
if not all(math.isfinite(x) for x in _g1) or max(_g1) > 1e-12:
    sys.exit("[STOP] G1: the base variant does not fall back on the published values (max gap %.3e > 1e-12).\n"
             "  binning_variants.csv / binning_contrasts.csv against test1_delta_north.csv and "
             "tache1_diff_bras.csv." % max(_g1))
print('[binning] G1: base == ct1/ct2, max gap %.1e' % max(_g1))

## ---- les grandeurs ----
_bD = {v: _bc.loc[v] for v in _BVARS}
_bCold = {v: _bv.loc[(v, 'daily_ref')] for v in _BVARS}
_shift_d = (_bD['between']['diff'] - _bD['base']['diff']) / _bD['base']['se_diff']
_shift_c = (_bCold['between']['f_m10'] - _bCold['base']['f_m10']) / _bCold['base']['f_m10_se']
_bm = _bM.set_index('quantity')['value']
for _q in ('part_intra_ref_min', 'part_intra_ref_max', 'ytrend_cor_dTbar',
           'warming_trend_cor_tbase_popw'):
    if _q not in _bm.index or not math.isfinite(_bm[_q]):
        sys.exit("[STOP] binning_mechanism.csv: %s missing or not finite." % _q)
_bins = [(float(re.search(r'lower edge (-?(?:[Ii]nf|[0-9.]+))\)', d).group(1)), float(v))
         for q, v, d in zip(_bM['quantity'], _bM['value'], _bM['definition'])
         if q.startswith('bin_trend_cor_tbase|')]
_bw = _bW.set_index(['channel', 'arm'])
try:
    _crh_p = 100 * float(_bw.loc[('prod', 'daily (s*)'), 'intra_resid_med'])
    _crh_a = 100 * float(_bw.loc[('amen', 'daily (s*)'), 'intra_resid_med'])
    _crh_tp = 100 * float(_bw.loc[('prod', 'mean_daily'), 'intra_resid_med'])
    _crh_ta = 100 * float(_bw.loc[('amen', 'mean_daily'), 'intra_resid_med'])
except KeyError as _e:
    sys.exit("[STOP] crh_within_cell.csv: row missing (%s)." % _e)

## ---- gardes de phrase G2-G7 ----
_G = [('G2', _g78.garde_g2(_shift_d, _shift_c),
       "B1/B2/B3: \"the design does not use that variation\" -- the between variant moves "
       "D or the cold arm by half a standard error or more."),
      ('G3', _g78.garde_g3(_bCold['pp']['f_m10_lo'], _bD['pp']['ci_lo'], _bD['pp']['ci_hi'],
                           _bD['base']['diff'], _bv.loc[('pp', 'daily_ref'), 'delta_north'],
                           _bv.loc[('pp', 'annual_ref'), 'delta_north'],
                           _bv.loc[('base', 'daily_ref'), 'delta_north'],
                           _bv.loc[('base', 'annual_ref'), 'delta_north']),
       "B2/B3: \"positive and significant\", \"reverses the contrast\", \"turns negative "
       "under both arms\"."),
      ('G4', _g78.garde_g4(_bD['pp5']['ci_lo'], _bD['pp5']['ci_hi'],
                           _bCold['pp5']['f_m10'], _bCold['pp']['f_m10']),
       "B2/B3: \"which includes zero\", \"the cold arm falls to\"."),
      ('G5', _g78.garde_g5(100 * _bm['part_intra_ref_max'], _bm['ytrend_cor_dTbar']),
       "B2/B3: the ingredients of the bias are no longer small."),
      ('G6', _g78.garde_g6(_crh_p, _crh_a, _crh_tp, _crh_ta),
       "B4: the description of supplement E (within-cell share) no longer holds."),
      ('G7', _g78.garde_g7(_bm['warming_trend_cor_tbase_popw'], _bins),
       "B3: \"faster in cold countries ... positively ... in both tails and not in the centre\".")]
## G8 (lot 79) : supplement E, « The design has no cell fixed effect ». Le code des
## estimateurs CRH (formules d'effets fixes) ET les donnees (part intra-cellule < 1).
import gardes_lot79 as _g79   # noqa: E402
import glob as _glob   # noqa: E402
_crh_src = [os.path.join(PKG, '01_estimation', '01d_crh_triangulation.R')] + \
    sorted(_glob.glob(os.path.join(PKG, '03_inference', 'step*.R')))
_crh_txt = {os.path.relpath(p, PKG): read(p) for p in _crh_src}
_G.append(('G8', _g79.garde_g8(_crh_txt, [float(x) for x in _bW['intra_resid_max']]),
           "S4: \"The design has no cell fixed effect: identification is across cells within "
           "a sub-region and period\"."))
for _nom, (_ok, _det), _phr in _G:
    print('[binning] %s %s: %s' % (_nom, 'ok' if _ok else 'REFUSED', _det))
    garde_phrase(_ok, 'binning %s' % _nom, 'true', _det,
                 "[STOP] %s refuses: %s\n  %s\n  The sentence is false: rewrite it, "
                 "do not relax the guard." % (_nom, _det, _phr))

## ---- les macros ----
_SRCV = 'contribA/binning_variants.csv <- 04_contrast/ct8_binning_variants.R'
_SRCC = 'contribA/binning_contrasts.csv <- 04_contrast/ct8_binning_variants.R'
_SRCM = 'contribA/binning_mechanism.csv <- 04_contrast/ct9_binning_mechanism.R'
_SRCW = 'contribA/crh_within_cell.csv <- 04_contrast/ct9_binning_mechanism.R'
for _v, _nm, _e, _elo, _ehi in (('base', 'EMPbinDbase', '+0.0169', None, None),
                                ('between', 'EMPbinDbetween', '+0.0172', None, None),
                                ('pp', 'EMPbinDpp', '-0.0108', '-0.0167', '-0.0049'),
                                ('pp5', 'EMPbinDppFive', '-0.0023', '-0.0079', '+0.0032')):
    _s = '%s, ref_ct2 (%s), diff' % (_SRCC, _v)
    add(_nm, _sg(_bD[_v]['diff'], 4), _e, _s)
    if _elo is not None:
        add(_nm + 'Lo', _sg(_bD[_v]['ci_lo'], 4), _elo, _s.replace('diff', 'ci_lo'))
        add(_nm + 'Hi', _sg(_bD[_v]['ci_hi'], 4), _ehi, _s.replace('diff', 'ci_hi'))
add('EMPbinShiftDbetween', _sg(_shift_d, 2), '+0.05',
    _SRCC + ': (diff between - diff base) / se_diff base, ref_ct2')
add('EMPbinShiftColdBetween', _sg(_shift_c, 2), '-0.24',
    _SRCV + ': (f_m10 between - f_m10 base) / f_m10_se base, daily_ref')
for _v, _nm, _e, _elo, _ehi in (('base', 'EMPbinColdBase', '-0.178', '-0.396', '+0.040'),
                                ('pp', 'EMPbinColdPp', '+0.310', '+0.143', '+0.478'),
                                ('pp5', 'EMPbinColdPpFive', '+0.190', '+0.087', '+0.292')):
    _s = '%s, %s / daily_ref, f(-10) - f(15)' % (_SRCV, _v)
    add(_nm, _sg(_bCold[_v]['f_m10'], 3), _e, _s + ' (f_m10)')
    add(_nm + 'Lo', _sg(_bCold[_v]['f_m10_lo'], 3), _elo, _s + ' (f_m10_lo)')
    add(_nm + 'Hi', _sg(_bCold[_v]['f_m10_hi'], 3), _ehi, _s + ' (f_m10_hi)')
add('EMPbinNorthPpDaily', _sg(_bv.loc[('pp', 'daily_ref'), 'delta_north'], 4), '-0.0157',
    _SRCV + ', pp / daily_ref, delta_north')
add('EMPbinNorthPpAnnual', _sg(_bv.loc[('pp', 'annual_ref'), 'delta_north'], 4), '-0.0265',
    _SRCV + ', pp / annual_ref, delta_north')
add('EMPbinWithinPairMin', '%.1f' % (100 * _bm['part_intra_ref_min']), '0.3',
    _SRCM + ', part_intra_ref_min x 100')
add('EMPbinWithinPairMax', '%.1f' % (100 * _bm['part_intra_ref_max']), '1.9',
    _SRCM + ', part_intra_ref_max x 100')
add('EMPbinTrendCorY', '%.2f' % _bm['ytrend_cor_dTbar'], '0.12', _SRCM + ', ytrend_cor_dTbar')
add('EMPcrhWithinCellProd', '%.0f' % _crh_p, '6', _SRCW + ', prod / daily (s*), intra_resid_med x 100')
add('EMPcrhWithinCellAmen', '%.0f' % _crh_a, '25', _SRCW + ', amen / daily (s*), intra_resid_med x 100')
add('EMPcrhWithinCellTmeanProd', '%.1f' % _crh_tp, '0.6', _SRCW + ', prod / mean_daily, intra_resid_med x 100')
add('EMPcrhWithinCellTmeanAmen', '%.1f' % _crh_ta, '2.4', _SRCW + ', amen / mean_daily, intra_resid_med x 100')

## ---- la table : output/tab_binning.tex ----
## Les etiquettes de ligne viennent de la colonne fe_structure du CSV (ecrite par
## ct8) : aucune variante n'est nommee ni designee ici.
def _bfmt(x, lo, hi, nd):
    return '$%s$ & ' % _sg(x, nd), '$[%s, %s]$ & ' % (_sg(lo, nd), _sg(hi, nd))

_trows = []
for _v in _BVARS:
    _lab = str(_bv.loc[(_v, 'daily_ref'), 'fe_structure']).replace(' x ', r'$\times$')
    _lab = _lab[0].upper() + _lab[1:]
    _cells = [_bfmt(_bv.loc[(_v, 'daily_ref'), 'delta_north'], _bv.loc[(_v, 'daily_ref'), 'ci_lo'],
                    _bv.loc[(_v, 'daily_ref'), 'ci_hi'], 4),
              _bfmt(_bv.loc[(_v, 'annual_ref'), 'delta_north'], _bv.loc[(_v, 'annual_ref'), 'ci_lo'],
                    _bv.loc[(_v, 'annual_ref'), 'ci_hi'], 4),
              _bfmt(_bD[_v]['diff'], _bD[_v]['ci_lo'], _bD[_v]['ci_hi'], 4),
              _bfmt(_bCold[_v]['f_m10'], _bCold[_v]['f_m10_lo'], _bCold[_v]['f_m10_hi'], 3)]
    _trows.append(_lab + ' & ' + ''.join(c[0] for c in _cells)[:-3] + r' \\')
    _trows.append(' & ' + ''.join(c[1] for c in _cells)[:-3] + r' \\[2pt]')
_tref = int(_bv.loc[('base', 'daily_ref'), 'T_ref'])
_tab = [
    '%% GENERATED by 02_GE/emit_values_empirics.py -- DO NOT EDIT BY HAND.',
    '%% Source: contribA/binning_variants.csv, binning_contrasts.csv <- 04_contrast/ct8_binning_variants.R',
    r'\begin{table}[H]\centering\footnotesize\setlength{\tabcolsep}{3pt}',
    r'\caption{Binning bias: the productivity dose-response under four fixed-effect structures.}',
    r'\label{tab:binning}',
    r'\begin{tabular}{>{\raggedright\arraybackslash}p{3.3cm}cccc}',
    r'\toprule',
    r'Fixed effects & $\Delta_{\mathrm{North}}$, daily & $\Delta_{\mathrm{North}}$, annual '
    r'& Contrast $D$ & $f(-10)-f(%d)$, daily \\' % _tref,
    r'\midrule',
] + _trows + [
    r'\bottomrule',
    r'\multicolumn{5}{p{\dimexpr\linewidth-2\tabcolsep\relax}}{\footnotesize',
    r'\textit{Note.} Each row re-estimates the two reference specifications of',
    r'Section~\ref{sec:contrast} under the structure named in the first column; with',
    r"pair-mean differentials, each pair's bilateral temperature differentials are",
    r'replaced by their means over the estimation window.',
    r'$\Delta_{\mathrm{North}}$ is the projected change in log productivity north of',
    r'$45^\circ$N; $D$ is the annual minus the daily change, with the covariance',
    r'between arms; the last column is the daily response at $-10^\circ$C relative to',
    r'$%d^\circ$C. 95\%% intervals in brackets, standard errors clustered by exporter' % _tref,
    r'and importer.',
    r'}',
    r'\end{tabular}',
    r'\end{table}',
]

# ---------------------------------------------------------------- rapport ----
print(f"{'macro':<20} {'emitted':>14} {'expected':>14}  status")
print('-' * 78)
bad = 0
for name, val, exp, src, status in macros:
    print(f"\\{name:<19} {val:>14} {exp:>14}  {status}")
    if status.startswith('MISMATCH'):
        bad += 1
        print(f"{'':20} !! source: {src}")
print('-' * 78)

## Mode rapport : les ecarts forment le releve (valeur attendue -> valeur emise).
_ecarts = [(n, e, v, s) for n, v, e, s, st in macros if st.startswith('ECART-REPRISE')]
if REPRISE:
    print(f"[REPORT MODE] {len(_ecarts)} difference(s) reported out of "
          f"{sum(1 for m in macros if not m[4].startswith('REGEN'))} macros compared")
    if not DRY:
        _pe = os.path.join(OUT, 'ecarts_reprise_empirics.tsv')
        with open(_pe + '.tmp', 'w', encoding='utf-8', newline='\n') as _fh:
            _fh.write('# %s\nmacro\tattendu\temis\tsource\n' % bandeau(REPRISE, 'emit_values_empirics.py'))
            for n, e, v, s in _ecarts:
                _fh.write(f'{n}\t{e}\t{v}\t{s}\n')
        os.replace(_pe + '.tmp', _pe)
        print(f"[write] output/ecarts_reprise_empirics.tsv ({len(_ecarts)} rows)")

# ---------------------------------------------------------------- emission ---
lines = [
    '% ====================================================================',
    '% GENERATED by 02_GE/emit_values_empirics.py',
    '% DO NOT EDIT BY HAND. Purely EMPIRICAL figures (no \\GE*).',
    '% Macros marked LITERAL carry a provenance still to be closed (see the script).',
    '% ====================================================================',
    '',
]
import gardes_lot74 as _g74   # noqa: E402
for name, val, exp, src, status in macros:
    # LOT 74 : garde d'EMISSION, bloquante meme en reprise -- aucun zero signe.
    if not status.startswith('MISMATCH') and _g74.zero_signe(val):
        raise SystemExit('[STOP] %s would print %r: a signed zero.' % (name, val))
    if status.startswith('MISMATCH'):
        lines.append(f"% \\{name} NOT EMITTED -- MISMATCH ({val} vs {exp}) -- {src}")
    else:
        lines.append(f"\\newcommand{{\\{name}}}{{{val}}}   % {src}")

##################################################################################
# LA NATURE DE CHAQUE PAIRE : INTERVALLE OU ETENDUE. Posee le 20/09/2026.
#
# LE DEFAUT QU'ELLE FERME. Le suffixe `Lo`/`Hi` porte DEUX conventions dans ce
# paquet : une borne de confiance (`ci95_lo/hi`) et un min/max sur un ensemble.
# Releve du 20/09 : 39 paires en tout, dont 14 citees dans le flux vivant.
# Le `.tex` les rend identiquement, et ce qui les distingue aujourd'hui n'est
# pas le NOM mais la TOURNURE -- `[a,b]` pour un intervalle, « A to B » pour une
# etendue. Une discipline qui vit dans la prose est une discipline qu'un futur
# auteur ne peut pas voir.
#
# ⛔ ET ELLE A DEJA MORDU. `\EMPohSlopeLo/Hi` sont le min et le max de la pente
# PONCTUELLE sur les schemas identifies -- une ETENDUE -- et le manuscrit les rend
# (le 20/09, quand la garde a mordu, la source etait `slope_q50` ; elle est
# passee a la pente ponctuelle le 22/09, la NATURE de la paire est inchangee)
# `$[\EMPohSlopeLo,\EMPohSlopeHi]$ per $^\circ$C`, c'est-a-dire comme un
# intervalle de confiance sur la pente chaude. Ce n'en est pas un.
#
# CE QUI EST DECLARE ICI, ET CE QUI NE L'EST PAS. Les paires CITEES dans le flux
# vivant, parce que c'est la que le lecteur peut etre trompe. Une paire emise et
# jamais citee ne peut induire personne en erreur ; le controle 10 ne la reclame
# donc pas, et le renommage des huit autres attend le lot post-soumission.
#
# `interval` : les deux bornes encadrent UNE quantite, avec un niveau de
#              confiance. Rendu attendu : $[\Lo,\Hi]$.
# `range`    : les deux bornes sont le min et le max d'un ENSEMBLE de valeurs
#              (grille de specifications, schemas, tirages). Rendu attendu :
#              « \Lo to \Hi », « \Lo--\Hi », « between \Lo and \Hi ».
##################################################################################
PAIR_NATURE = {
    # --- intervalles de confiance -------------------------------------------
    'EMPcrhNorthDiff':    ('interval', 'ci95 clustered by country (crh_inference_deliverable)'),
    'EMPcrhAtten':        ('interval', 'ci95 clustered by country, attenuation_ratio'),
    # --- etendues sur un ensemble -------------------------------------------
    'EMPamenPjoint':      ('range', 'min/max of the joint p, daily arm'),
    'EMPamenAnnualPjoint':('range', 'min/max of the joint p, annual arm'),
    'EMPamenAnnualPnonlin':('range', 'min/max of the non-linear p (degree >= 2), annual arm'),
    'EMPcrhCorArms':      ('range', 'min/max of cor_arms'),
    'EMPctrlAbsorp':      ('range', 'min/max of the absorbed share'),
    'EMPdiffRho':         ('range', 'min/max of rho over the grid'),
    'EMPdiffSEshare':     ('range', 'min/max of se_ratio'),
    'EMPfeAbsorption':    ('range', 'min/max of the pair x product R2'),
    'EMPgridAnnual':      ('range', 'min/max of delta_north, annual arm'),
    'EMPgridDaily':       ('range', 'min/max of delta_north, daily arm'),
    'EMPohSlope':         ('range', 'min/max of slope_q50 over the identified schemes'),
    'EMPprodPnonlinExc':  ('range', 'min/max of the non-linearity p, exceptions'),
    'EMPspecInvar':       ('range', 'min/max of fz30 over the grid'),
    'EMPattenNorth':      ('range', 'min/max of diff/annual over the grid'),
    'EMPfinGainRatio':    ('range', 'min/max of the ratio over the grid -- NOT CITED'),
    'GEnorthDailyRq':     ('range', 'quantiles 5 and 95 of the arm B outputs'),
    # --- lot 78, biais de binning (ct8, ct9) ----------------------------------
    'EMPbinDpp':          ('interval', 'ci95 of the contrast D, variant pp (binning_contrasts.csv)'),
    'EMPbinDppFive':      ('interval', 'ci95 of the contrast D, variant pp5 (binning_contrasts.csv)'),
    'EMPbinColdBase':     ('interval', 'ci95 of f(-10)-f(15), daily, base (binning_variants.csv)'),
    'EMPbinColdPp':       ('interval', 'ci95 of f(-10)-f(15), daily, pp (binning_variants.csv)'),
    'EMPbinColdPpFive':   ('interval', 'ci95 of f(-10)-f(15), daily, pp5 (binning_variants.csv)'),
    'EMPbinWithinPair':   ('range', 'min/max of the within-pair share over the reference regressors'),
}

# --- le sidecar que check_macros controle 10 lit -------------------------
# Meme dispositif que values_provenance.tsv : l'emetteur DECLARE, le controle
# VERIFIE. Un fichier, pas un commentaire.
_pairs_txt = 'base\tnature\tsource\n' + ''.join(
    '%s\t%s\t%s\n' % (b, n, src)
    for b, (n, src) in sorted(PAIR_NATURE.items()))
if not DRY:
    ## Lot 86d : toutes les ecritures de cet emetteur sont en LF explicite.
    with open(os.path.join(OUT, 'values_pairs.tsv'), 'w', encoding='utf-8', newline='\n') as _fh:
        _fh.write(_pairs_txt)
    print('[write] output/values_pairs.tsv: %d pairs declared' % len(PAIR_NATURE))

dest = os.path.join(OUT, 'values_empirics.tex')
if DRY:
    print(f"[dry-run] NOTHING written. {len(macros)} macros parsed, {bad} mismatch.")
    print(f"[dry-run] target that would have been written: {dest}")
else:
    with open(dest, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"[write] output/values_empirics.tex ({len(macros)} macros, {bad} mismatch)")
    ## lot 78 : la table du biais de binning, ecriture atomique
    _tb = os.path.join(OUT, 'tab_binning.tex')
    with open(_tb + '.tmp', 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(_tab) + '\n')
    os.replace(_tb + '.tmp', _tb)
    print(f"[write] output/tab_binning.tex ({len(_BVARS)} variant rows)")
# CODE DE SORTIE DISTINCT POUR LE VERDICT (2026-09-07).
# Ce script sortait en 1 sur un mismatch de macro -- le MEME code que ses
# neuf gardes dures (sys.exit(str) l.163/211/215/266/280/710 et
# raise SystemExit(str) l.654/659/685) et que toute exception non
# rattrapee, Python rendant 1 dans les trois cas. Un appelant ne pouvait
# donc PAS distinguer le verdict d'un plantage : un FileNotFoundError se
# lisait << MISMATCH >> et le run etait declare reussi.
#   0 = toutes les macros concordent
#   1 = garde dure ou exception -- le script n'a pas rendu son verdict
#   3 = verdict rendu, au moins une macro en non-concordance
sys.exit(3 if bad else 0)
