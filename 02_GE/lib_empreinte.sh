# lib_empreinte.sh -- l'empreinte des intrants du GE, cote run_all.sh (lot 53, 28/09/2026).
#
# MEME RECETTE que 02_GE/run_one_ge.sh (l.99-114) et que fingerprint() de
# 02_GE/02h_exposure_mask_runs.py : le md5 des md5 des neuf fichiers, un par
# ligne, SUIVI de la ligne SCALE_CONVENTION lue dans param_theta.
#
# Pourquoi ce fichier : run_all.sh hachait les neuf fichiers SANS la ligne
# SCALE_CONVENTION. 02h compare l'empreinte vivante (neuf + convention) a
# celle que run_all.sh ecrit au manifeste : l'arret en `maskruns` etait certain,
# APRES les treize heures de la grille GE. Teste par tools/test_empreinte_maskruns.py.
#
# Usage : source "$JMP_ROOT/02_GE/lib_empreinte.sh" ; INPUTS_MD5="$(empreinte_intrants)"

convention_echelle() {
  # theta lu a sa source unique, comme run_one_ge.sh l.87-98 : une occurrence, aucun repli.
  local f="$JMP_ROOT/02_GE/02c_model/specification_study.csv" hit theta
  [ -f "$f" ] || { echo "[FATAL] specification_study.csv introuvable : $f" >&2; return 3; }
  hit="$(grep -E '^[[:space:]]*param_theta[[:space:]]*;' "$f" || true)"
  if [ "$(printf '%s' "$hit" | grep -c .)" != "1" ]; then
    echo "[FATAL] cle 'param_theta' absente ou multiple dans $f. Aucun repli." >&2; return 3
  fi
  theta="$(printf '%s' "$hit" | cut -d';' -f2 | tr -d '[:space:]')"
  case "$theta" in
    ''|*[!0-9.]*) echo "[FATAL] param_theta illisible : '$hit'" >&2; return 3 ;;
  esac
  printf 'SCALE_PROD=1/THETA_GE;THETA_GE=%s' "$theta"
}

empreinte_intrants() {
  # Les lignes sont assemblees AVANT le hachage : dans un `{ ... } | md5sum`, un
  # refus a l'interieur du bloc laisserait md5sum rendre une empreinte quand meme.
  local sc f lignes=""
  sc="$(convention_echelle)" || return 3
  for f in lambda.h5 mig_flows.h5 indices_drop.csv identifier_xy_final.csv \
           population0.csv tas_ssp2_daily.csv tas_ssp2_annual.csv \
           data_final_population.csv data_final_productivity.csv ; do
    [ -f "$JMP_ROOT/intermediate_outputs/$f" ] || { echo "[FATAL] intrant GE absent : $f" >&2; return 3; }
    lignes="$lignes$(md5sum "$JMP_ROOT/intermediate_outputs/$f" | cut -d' ' -f1)
"
  done
  printf '%s%s\n' "$lignes" "$sc" | md5sum | cut -d' ' -f1
}
