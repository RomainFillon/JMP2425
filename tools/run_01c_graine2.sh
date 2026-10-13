#!/usr/bin/env bash
# run_01c_graine2.sh -- l'etape 01c_graine2 de run_all.sh (lot 74, 01/10/2026).
#
# 01c_robustness.R a la seconde graine (JMP_OUTPUT_DIR et JMP_MC_SEED_PROD_DAILY
# poses par run_all.sh), PUIS la garde de reproductibilite : la sortie doit
# reproduire la reference output_01c_graine2_20260930/ A 1e-8 PRES, cellule par
# cellule (tools/gardes_lot74.py comparer).
# Lot 85f (04/10/2026) : la garde lit la copie a l'octet de ces trois CSV dans
# 01_estimation/reference_repro_01c_graine2/ (dans le perimetre ; output_*/ ne
# l'est pas : un clone rendait ABSENT puis STOP).
#
# POURQUOI DANS L'ETAPE ET PAS APRES ELLE. Le 30/09, la verification etait une
# boucle `cmp` tapee a la main apres `--until 01c_graine2` : octet a octet, elle a
# rendu « differ » pour deux fichiers dont l'ecart max etait 1,1e-12 (ordre de
# sommation flottant) -- bruit publie inchange. Un critere d'octets refuse ce que
# le papier ne voit pas ; le critere est donc NUMERIQUE, au seuil de la garde
# d'emit_values_empirics.py sur les estimations ponctuelles (1e-8). Et il vit
# dans l'etape, parce que arret_si_until sort juste apres elle : une garde posee
# a la ligne suivante de run_all.sh ne tournerait jamais sous --until.
#
# Les trois CSV comparees portent tout ce que la table tab_oh_bounds et les
# macros EMPohSeed* lisent ; tab_oh_bounds.tex en est une fonction (01c_reemit).
set -u
: "${RSCRIPT:?RSCRIPT non pose}" "${PYTHON:?PYTHON non pose}" "${JMP_ROOT:?JMP_ROOT non pose}"
"$RSCRIPT" "$JMP_ROOT/01_estimation/01c_robustness.R" || exit $?
echo "--- garde de reproductibilite (1e-8) contre 01_estimation/reference_repro_01c_graine2/ ---"
"$PYTHON" "$JMP_ROOT/tools/gardes_lot74.py" comparer \
  "$JMP_ROOT/01_estimation/reference_repro_01c_graine2" "$JMP_OUTPUT_DIR" \
  oh_hotslope_schemes.csv oh_scheme_effects.csv oster_delta_full.csv
