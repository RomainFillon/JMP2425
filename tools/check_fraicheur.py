#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""check_fraicheur.py -- un artefact publie est-il plus vieux que ses intrants ?

LA REGLE, ET ELLE EST MECANIQUE :

    ⇒ UN ARTEFACT PLUS VIEUX QUE LA DONNEE QU'IL LIT EST PERIME,
      QUOI QUE DISE SON CONTENU.

Elle ne demande ni de lire une couche texte, ni d'inferer une magnitude
attendue. Elle a tranche les Figures 3 et 4 le 22/09/2026 la ou le raisonnement
sur les echelles ne donnait qu'une presomption : les deux figures dataient du
18/09 17:58 et les trois grilles qu'elles lisent des 20 et 21/09.

⛔ POURQUOI LA TABLE EST TENUE A LA MAIN. `audit_paquet.py` construit le graphe
   automatiquement, mais il porte encore 83 indetermines et cinq aveuglements
   connus (variables de chemin, f-strings, etapes declarees sautees,
   dependances intra-etape, troncature). Onze figures et une poignee de tables
   se recensent en une lecture : ici la main est plus sure que le graphe, et
   c'est le seul endroit de ce depot ou je l'ecris.

⛔ ET UN ARTEFACT NON RECENSE NE PASSE PAS AU VERT. Il est rapporte comme
   « non recense », rc=2 -- un controle qui ne peut pas conclure ne conclut pas.

🚨 LIMITE MESUREE AU PREMIER PASSAGE, ET ELLE EST DANS UN SEUL SENS. Le mtime
est un PROXY du contenu. Ce controle a signale `fig_contrast_grid.pdf` et
`tab_contrast_grid.tex` comme perimes parce que `contribA_grid.csv` portait une
date plus recente -- or le CSV avait ete RESTAURE apres un test de mutation, et
son md5 est identique au bit (`7bad475a` des deux cotes). `cp` deplace
l'horodatage sans toucher au contenu.

⇒ IL SUR-RAPPORTE DANS UN SENS : une restauration, un `touch`, un `cp`
  produisent une alerte sans defaut. **Verifier le md5 de l'intrant signale
  avant de regenerer.** C'est la distinction « horodatage seul / contenu
  modifie » que `cmp_rerun.py` fait depuis le 19/09, et qu'il reste a porter
  ici -- liste post-soumission.

🚨 ET IL SOUS-RAPPORTE DANS L'AUTRE, SOUS UNE CONDITION QU'IL FAUT ECRIRE.
« Un intrant modifie a forcement un mtime plus recent » suppose que TOUTE
ecriture met l'horodatage a jour. C'est vrai d'une ecriture ordinaire et FAUX
de toute operation qui preserve les horodatages :

    cp -p / cp --preserve=timestamps · rsync -t / rsync -a
    tar -x et unzip (qui preservent PAR DEFAUT) · touch -t

Sous l'une d'elles, un contenu MODIFIE arrive avec un mtime ANCIEN, et un
perime reel echappe a la garde. C'est la direction qui compte, et ce depot est
exactement la ou elle mord : la procedure courante est `.bak_<date>` plus un
instantane hors arbre, et une restauration depuis un instantane du 20 avec une
copie qui preserve rendrait fraiche une figure du 22 qui ne l'est pas.

⇒ LA PREMISSE N'EST PAS UNE PROPRIETE DU SYSTEME DE FICHIERS, C'EST UNE
  DISCIPLINE DU DEPOT : **une restauration ne preserve jamais les
  horodatages** -- `cp` nu, jamais `cp -p`. Ecrite dans tools/README.md a cote
  de la procedure de sauvegarde. Sous cette regle, et seulement sous elle,
  l'asymetrie ci-dessus est exacte.

Usage:  python tools/check_fraicheur.py
Sortie: rc=0 tout frais, rc=1 au moins un perime, rc=2 non concluant.
"""
import io
import os
import sys

ROOT = os.environ.get('JMP_FRAICHEUR_ROOT') or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'output')

# ----------------------------------------------------------------------------
# artefact publie -> (producteur, [intrants])
# Les chemins sont relatifs a la racine du paquet.
# ----------------------------------------------------------------------------
# Les runs du CORPS (f_a = 0), ceux que lisent emit_figures.py et emit_values.py.
_BODY = ['output/welfare_sc%d_q50_NULLamen.csv' % s for s in (1, 2, 3)]
# Les soixante tirages de bande, que emit_values.py lit pour fig_ge_tiers.csv.
_BAND = ['output/welfare_sc2_rb%02d_NULLamen.csv' % k for k in range(1, 61)]
# La geographie que emit_figures.py lit pour agreger par pays.
_GEO = ['intermediate_outputs/%s.csv' % f for f in
        ('indices_drop', 'identifier_xy_final', 'population0', 'poids_par_cellule')]

TABLE = {
    # --- les trois figures d'equilibre general -----------------------------
    # 🚨 CORRIGE LE 25/09/2026. Cette table listait welfare_sc*_q50.csv, les
    # runs du bras INTERIEUR, alors que emit_figures.py (l.100-125) lit le bras
    # du corps, welfare_sc*_q50_NULLamen.csv, et que la figure des quartiles ne
    # lit aucun run : elle lit fig_ge_tiers.csv, ecrit par emit_values.py. La
    # garde n'a donc jamais rien certifie sur les Figures 3 et 4, et elle se
    # taisait comme une garde qui fonctionne. Les intrants ci-dessous sont ceux
    # que les producteurs lisent, relus dans leur source.
    'output/fig_ge_tiers.pdf': (
        '02_GE/emit_figures.py',
        ['output/fig_ge_tiers.csv']),
    'output/fig_ge_tiers.csv': (
        '02_GE/emit_values.py',
        _BODY + _BAND),
    'output/map_welfare_2100_country.pdf': (
        '02_GE/emit_figures.py',
        _BODY + _GEO),
    'output/ge_welfare_trajectories.pdf': (
        '02_GE/emit_figures.py',
        _BODY + _GEO),
    # --- la grille de contraste --------------------------------------------
    'output/fig_contrast_grid.pdf': (
        '04_contrast/ct4_emit_tab_contribA_grid.R',
        ['output/contribA_grid.csv']),
    'output/tab_contrast_grid.tex': (
        '04_contrast/ct4_emit_tab_contribA_grid.R',
        ['output/contribA_grid.csv']),
    'output/tab_robustness_prod_joint.tex': (
        '04_contrast/ct4_emit_tab_contribA_grid.R',
        ['output/prod_joint_tests_grid_daily.csv',
         'output/prod_joint_tests_grid_annual.csv']),
    # --- le bloc OH ---------------------------------------------------------
    'output/tab_oh_bounds.tex': (
        '01_estimation/01c_robustness.R',
        ['output/oh_scheme_effects.csv', 'output/oh_hotslope_schemes.csv']),
    # --- l'absorption (le cas du matin) ------------------------------------
    'output/tab_ctrl_absorption.tex': (
        '01_estimation/01a_productivity_drf.R',
        ['output/ctrl_absorption.csv']),
    'output/tab_fe_absorption.tex': (
        '01_estimation/01a_productivity_drf.R',
        ['output/fe_absorption.csv']),
    # --- les eventails Monte Carlo -----------------------------------------
    'output/MC_prod_daily_fan_paper.pdf': (
        '01_estimation/01a_productivity_drf.R',
        ['output/MC_damage_prod_daily_quantiles.csv']),
    'output/MC_prod_annual_fan_paper.pdf': (
        '01_estimation/01a_productivity_drf.R',
        ['output/MC_damage_prod_annual_quantiles.csv']),
    'output/MC_amenity_daily_fan_paper.pdf': (
        '01_estimation/01b_amenity_drf.R',
        ['output/MC_damage_amenity_daily_quantiles.csv']),
    'output/MC_amenity_annual_fan_paper.pdf': (
        '01_estimation/01b_amenity_drf.R',
        ['output/MC_damage_amenity_annual_quantiles.csv']),
    # --- les quatre qui manquaient au premier passage ----------------------
    # Aucune n'etait « sans producteur » : c'etait une omission de PERIMETRE de
    # ma part, pas une trouvaille. Les quatre producteurs se trouvent par un
    # grep, et les recenser vaut mieux que de les declarer hors perimetre.
    'output/fig_crh_responses.pdf': (
        '01_estimation/01d_fig_crh_responses.R',
        ['intermediate_outputs/crh/crh_support_bounds.csv',
         'intermediate_outputs/crh/amen_prod.csv']),
    'output/fig_M1_net_by_latitude.png': (
        '02_GE/02f_fig_mig_latitude.py',
        ['intermediate_outputs/mig_matrix_final_raw.h5']),
    'output/myplot_migration.pdf': (
        '02_GE/02a_calibration_at_1deg.r',
        ['intermediate_outputs/identifier_xy_final.csv',
         'intermediate_outputs/poids_par_cellule.csv']),
    'output/no_activity_map_smooth_light.pdf': (
        '02_GE/02a_calibration_at_1deg.r',
        ['intermediate_outputs/identifier_xy_final.csv',
         'intermediate_outputs/poids_par_cellule.csv']),
}

# Artefacts publies dont les intrants ne sont PAS recenses. Ils ne passent pas
# au vert : ils sont NON CONCLUANTS, et la liste est la pour retrecir.
##############################################################################
# ⛔ UNE GARDE QUI RESTE A rc=2 EST UNE GARDE QU'ON ARRETE DE LIRE.
#
# Au premier passage ce controle rendait rc=2 PAR CONSTRUCTION, sur quatre
# figures non recensees. Un feu qui n'est jamais vert n'informe plus : au bout
# de trois jours, rc=2 se lit « normal ». C'est la meme logique qui a fait jeter
# le balayage d'octets a 53 faux positifs et le controle 4 a 80 entrees.
#
# Les quatre sont desormais RECENSEES -- aucune n'etait sans producteur, c'etait
# une omission de perimetre. Ce dictionnaire reste, vide, pour que la sortie
# dise EXPLICITEMENT combien d'artefacts sont ignores deliberement, sur la meme
# ligne que le verdict : la regle de troncature du 22/09.
##############################################################################
NON_RECENSES = {}


def mt(p):
    q = os.path.join(ROOT, p)
    return os.path.getmtime(q) if os.path.exists(q) else None


def main():
    perimes, absents, ok = [], [], []
    for art, (prod, intrants) in sorted(TABLE.items()):
        ta = mt(art)
        if ta is None:
            absents.append('%s: ARTEFACT MISSING' % art)
            continue
        plus_recents = []
        for src in intrants:
            ts = mt(src)
            if ts is None:
                absents.append('%s: input missing -> %s' % (art, src))
                continue
            if ts > ta:
                plus_recents.append((src, ts))
        if plus_recents:
            perimes.append((art, prod, ta, plus_recents))
        elif not any(a.startswith(art) for a in absents):
            ok.append(art)

    import time

    def d(x):
        return time.strftime('%d/%m %H:%M', time.localtime(x))

    print('[fraicheur] %d artefact(s) listed, %d not listed'
          % (len(TABLE), len(NON_RECENSES)))
    print('\n1. ARTEFACTS OLDER THAN THEIR INPUTS: %d' % len(perimes))
    for art, prod, ta, srcs in perimes:
        print('   STALE     %s (%s)' % (art, d(ta)))
        print('             producer: %s' % prod)
        for s, ts in srcs:
            print('             NEWER input: %s (%s)' % (s, d(ts)))
    if not perimes:
        print('   none')

    print('\n2. INCONCLUSIVE: %d' % (len(absents) + len(NON_RECENSES)))
    for a in absents:
        print('   MISSING   ' + a)
    for a, pourquoi in sorted(NON_RECENSES.items()):
        print('   NOT LISTED   %-44s %s' % (a, pourquoi))

    print('\n3. FRESH: %d' % len(ok))
    print('\nVERDICT: %s (%d artefact(s) listed, %d deliberately '
          'ignored)'
          % ('%d stale artefact(s)' % len(perimes) if perimes
             else ('INCONCLUSIVE on %d artefact(s)' % len(absents)
                   if absents
                   else 'every listed artefact is newer than its inputs'),
             len(TABLE), len(NON_RECENSES)))
    if perimes:
        return 1
    return 2 if (absents or NON_RECENSES) else 0


if __name__ == '__main__':
    sys.exit(main())
