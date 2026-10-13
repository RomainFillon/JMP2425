#!/usr/bin/env bash
# =============================================================================
# run_one_ge.sh -- runs ONE GE run and records it in output/ge_runs_manifest.csv.
#
# The grid in run_all.sh runs all eight and rewrites the manifest from scratch.
# This script exists for the other case: adding one run to a grid whose other
# runs are already on disk and already current, without re-running work that
# would reproduce itself byte for byte.
#
# A run launched by hand and absent from the manifest is exactly what the
# manifest exists to prevent: decompose_metric_channels.py would refuse the
# triplet, and the refusal would be right. So the manifest line is written here,
# by the same recipe the grid uses, and only after a run that returned 0.
#
# The manifest is APPENDED to, and any existing line for the same output file is
# replaced, so that re-running one scenario does not leave two rows for it.
#
# USAGE
#   SCENARIO=3 QUANTILE=q50 ANNUAL_DRF_SOURCE=dailycurve \
#   CAP_COLMASS=none AMENITY_SPEC=baseline ./02_GE/run_one_ge.sh
#
# Reads:  02_GE/02c_model/main.py and everything it reads.
# Writes: output/welfare_<tag>.csv and the two hat-paths, one log under
#         logs/one_<tag>_<timestamp>/, one line in output/ge_runs_manifest.csv.
# =============================================================================
set -euo pipefail

JMP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
command -v cygpath >/dev/null 2>&1 && JMP_ROOT="$(cygpath -m "$JMP_ROOT")"
export JMP_ROOT

PYTHON="${PYTHON:-python}"
SCENARIO="${SCENARIO:?SCENARIO must be set}"
QUANTILE="${QUANTILE:?QUANTILE must be set}"
ANNUAL_DRF_SOURCE="${ANNUAL_DRF_SOURCE:-annual}"
CAP_COLMASS="${CAP_COLMASS:-none}"
AMENITY_SPEC="${AMENITY_SPEC:-baseline}"
export CAP_COLMASS AMENITY_SPEC

# The output name carries the amenity tag too, and main.py appends it after the
# scenario and quantile. A run whose file lands under a name the emitter does not
# look for is a run for nothing, so the two rules must stay in step.
case "$AMENITY_SPEC" in
  null)     ATAG="_NULLamen" ;;
  bell)     ATAG="_BELLamen" ;;
  bell_cap) ATAG="_BELLCAPamen" ;;
  *)        ATAG="" ;;
esac
if [ "$ANNUAL_DRF_SOURCE" = "dailycurve" ]; then
  # 🚨 L'ORDRE DES SUFFIXES SUIT main.py, IL NE LE DEVINE PAS.
  # main.py l.429 : suffix = _sc<label><q> + <amenite> + <cap> + <drf>, donc
  # l'etiquette d'amenite vient AVANT « _dailycurve ». Ce script construisait
  # « _dailycurve${ATAG} », l'ordre inverse. Avec AMENITY_SPEC=baseline ATAG est
  # vide et les deux coincident -- c'est pourquoi le defaut n'a jamais eu
  # d'effet visible. Au premier sc4 a amenity=null (16/09/2026), le run a
  # REUSSI, ecrit welfare_sc4_q50_NULLamen_dailycurve.csv, et ce script a
  # cherche welfare_sc4_q50_dailycurve_NULLamen.csv : introuvable, donc AUCUNE
  # ligne de manifeste pour un run qui avait pourtant tourne onze minutes.
  TAG="sc4_${QUANTILE}${ATAG}_dailycurve"
else
  TAG="sc${SCENARIO}_${QUANTILE}${ATAG}"
fi

GE_MANIFEST="$JMP_ROOT/output/ge_runs_manifest.csv"
HDR="file,md5,run_id,scenario,quantile,cap_colmass,amenity_spec,annual_drf_source,inputs_md5,seconds,rc,timestamp"

# The same fingerprint recipe as the grid. Keep the two in step: a different
# input list would make the two writers disagree about what "same inputs" means,
# and the guard would pass on a triplet it should refuse.
#
# -----------------------------------------------------------------------------
# THE CONVERSION CONSTANT IS PART OF THE FINGERPRINT. Added 20/09/2026.
# -----------------------------------------------------------------------------
# Until this date the recipe hashed nine CALIBRATION files and nothing else. It
# did not hash the dose-responses, and it did not hash the scale the estimation
# was run under. So the sixty runs of 16/09 -- made when the productivity channel
# still carried theta twice -- and the runs of 20/09 -- made after the division --
# came out with the SAME signature, 930ff8e0. Two runs under two different
# conventions were indistinguishable to every consumer of the manifest, and only
# the timestamp separated them. `run_nuit_reprise.sh` had to work around that
# with a date filter; a date filter is a workaround, the fingerprint is the fix.
#
# theta is READ AT ITS SINGLE SOURCE, the solver's own configuration file, by the
# same recipe as 00_parameters.R l.136-147: the key `param_theta`, one occurrence,
# no fallback. A conversion constant that degrades quietly is the failure mode
# being closed here, so an unreadable key STOPS the run rather than hashing a
# placeholder.
THETA_FILE="$JMP_ROOT/02_GE/02c_model/specification_study.csv"
[ -f "$THETA_FILE" ] || { echo "[FATAL] specification_study.csv introuvable : $THETA_FILE"; exit 3; }
THETA_HIT="$(grep -E '^[[:space:]]*param_theta[[:space:]]*;' "$THETA_FILE" || true)"
if [ "$(printf '%s' "$THETA_HIT" | grep -c .)" != "1" ]; then
  echo "[FATAL] cle 'param_theta' absente ou multiple dans $THETA_FILE. Aucun repli."
  exit 3
fi
THETA_GE="$(printf '%s' "$THETA_HIT" | cut -d';' -f2 | tr -d '[:space:]')"
case "$THETA_GE" in
  ''|*[!0-9.]*) echo "[FATAL] param_theta illisible : '$THETA_HIT'"; exit 3 ;;
esac
# La CHAINE de convention, pas la valeur seule : deux conventions differentes qui
# partageraient un theta (par exemple scale=1 et scale=1/theta) doivent sortir
# deux empreintes. Ce qui est hache est ce que le canal productivite applique.
SCALE_CONVENTION="SCALE_PROD=1/THETA_GE;THETA_GE=$THETA_GE"
INPUTS_MD5="$( { for f in \
      "$JMP_ROOT/intermediate_outputs/lambda.h5" \
      "$JMP_ROOT/intermediate_outputs/mig_flows.h5" \
      "$JMP_ROOT/intermediate_outputs/indices_drop.csv" \
      "$JMP_ROOT/intermediate_outputs/identifier_xy_final.csv" \
      "$JMP_ROOT/intermediate_outputs/population0.csv" \
      "$JMP_ROOT/intermediate_outputs/tas_ssp2_daily.csv" \
      "$JMP_ROOT/intermediate_outputs/tas_ssp2_annual.csv" \
      "$JMP_ROOT/intermediate_outputs/data_final_population.csv" \
      "$JMP_ROOT/intermediate_outputs/data_final_productivity.csv" ; do
    md5sum "$f" | cut -d' ' -f1
  done ; printf '%s\n' "$SCALE_CONVENTION" ; } | md5sum | cut -d' ' -f1)"

RUN_ID="one_${TAG}_$(date '+%Y%m%d_%H%M%S')"
LOGDIR="$JMP_ROOT/logs/$RUN_ID"
mkdir -p "$LOGDIR"

echo "[one] $TAG"
echo "     SCENARIO=$SCENARIO QUANTILE=$QUANTILE ANNUAL_DRF_SOURCE=$ANNUAL_DRF_SOURCE"
echo "     CAP_COLMASS=$CAP_COLMASS AMENITY_SPEC=$AMENITY_SPEC"
echo "     inputs_md5=$INPUTS_MD5"
echo "     debut : $(date '+%H:%M:%S')"

t0=$(date +%s)
set +e
SCENARIO="$SCENARIO" QUANTILE="$QUANTILE" ANNUAL_DRF_SOURCE="$ANNUAL_DRF_SOURCE" \
  "$PYTHON" -u "$JMP_ROOT/02_GE/02c_model/main.py" 2>&1 | tee "$LOGDIR/ge_${TAG}.log"
rc=${PIPESTATUS[0]}
set -e
t1=$(date +%s); d=$((t1 - t0))

if [ "$rc" -ne 0 ]; then
  echo "⛔ ECHEC $TAG (code $rc) apres $((d/60)) min $((d%60)) s. Log : $LOGDIR/ge_${TAG}.log"
  echo "   Aucune ligne ecrite au manifeste : un run echoue ne laisse pas de provenance."
  exit "$rc"
fi

WF="welfare_${TAG}.csv"
if [ ! -f "$JMP_ROOT/output/$WF" ]; then
  echo "⛔ $WF introuvable apres un run a 0. Rien n'est ecrit au manifeste." >&2
  exit 1
fi

[ -f "$GE_MANIFEST" ] || echo "$HDR" > "$GE_MANIFEST"
# 🚨 COURSE LECTURE-MODIFICATION-ECRITURE, MESUREE LE 16/09/2026.
# Ce bloc recopiait le manifeste sans sa propre ligne, remettait la copie en
# place, puis ajoutait la ligne. Deux runs qui finissent ensemble lisent donc
# le MEME etat, et le second `mv` ecrase la ligne que le premier venait
# d'ajouter. Sur soixante runs a deux processus, une ligne a ete perdue --
# rb41 -- alors que son log disait « inscrit au manifeste ». Le script
# affirmait une chose que le fichier ne portait pas.
# Le verrou serialise la section critique. `flock` est dans Git-Bash ; a
# defaut, on degrade vers un verrou par repertoire (mkdir est atomique).
_MANIFEST_LOCK="${GE_MANIFEST}.lock"
_maj_manifeste () {
  TMPM="$(mktemp)"
  grep -v "^${WF}," "$GE_MANIFEST" > "$TMPM" || true
  echo "$WF,$(md5sum "$JMP_ROOT/output/$WF" | cut -d' ' -f1),$RUN_ID,$SCENARIO,$QUANTILE,$CAP_COLMASS,$AMENITY_SPEC,$ANNUAL_DRF_SOURCE,$INPUTS_MD5,$d,$rc,$(date '+%Y-%m-%dT%H:%M:%S')" \
    >> "$TMPM"
  mv "$TMPM" "$GE_MANIFEST"        # un seul remplacement, atomique
}
if command -v flock >/dev/null 2>&1; then
  ( flock 9; _maj_manifeste ) 9>"$_MANIFEST_LOCK"
else
  for _i in $(seq 1 600); do
    mkdir "${_MANIFEST_LOCK}.d" 2>/dev/null && break
    sleep 1
  done
  _maj_manifeste
  rmdir "${_MANIFEST_LOCK}.d" 2>/dev/null || true
fi
if ! grep -q "^${WF}," "$GE_MANIFEST"; then
  echo "⛔ $WF ABSENT du manifeste apres ecriture : la ligne a ete perdue." >&2
  exit 1
fi

echo "     fin   : $(date '+%H:%M:%S')  duree $((d/60)) min $((d%60)) s"
echo "[ok] $WF inscrit au manifeste."
