#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- runs the whole replication package, end to end.
#
# Executes the estimation steps, then the general-equilibrium solver, then the
# emitters that produce the numbers and tables the manuscript reads. Each step
# writes its own log under logs/run_<timestamp>/.
#
# Two steps are NOT called because their inputs cannot be shipped; their outputs
# are supplied instead and checked in section 5. They are documented in place.
#
# The run ends with an inventory of what the chain did not produce, measured
# against the targets JMP_V11.tex actually loads (JMP_V10.tex is frozen).
#
# USAGE
#   ./run_all.sh                 # everything, in order
#   ./run_all.sh --from 02a      # resume at the named step
#   ./run_all.sh --list          # print the step names and exit
#   RSCRIPT=... PYTHON=... ./run_all.sh
#
# COST -- MEASURED WALL-CLOCK TIME, declared before you start (lot 34, 25/09/2026)
#   ge          78 GE runs, one solver at a time, ~10 min each (max 13 min):
#               787 min of solving, 14 h 07 of wall-clock including waits for
#               14 GB of free memory (logs/lot25_20260923_2355.etat). The band
#               alone (60 runs) is about 10 h of it.
#   02a         about 1 h 22                   (README, Timings)
#   01a         about 28 min                   (README, Timings)
#   ct1 + ct2   under 38 min (measured on the former 24-spec grid, a bound)
#   01f 2 min 01 s | ct7 11 s | ct5 6 s | 01a_reemit, 01c_reemit about 1 s each
#               (measured 25/09/2026)
#   ct8         5 min 28 s (lot 77 sandbox, 01/10/2026) | ct9 not yet measured
#   Every other step is listed with its measured time in README, Timings.
#   Plan for a full day on one machine. ./run_all.sh --from <step> resumes.
# =============================================================================

# No `|| true` and no `set +e` outside an explicit capture of the return code.
# A step that breaks stops the chain: continuing would build outputs on a missing
# input, and the closing inventory would report them as produced.
set -euo pipefail

FROM=""; LIST=0; UNTIL=""
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM="${2:-}"; shift 2 ;;
    ## --until <etape> (lot 51) : arret PROPRE apres l'etape nommee, pour les points
    ## d'arret de la reprise ; reprendre ensuite par --from <etape suivante>.
    --until) UNTIL="${2:-}"; shift 2 ;;
    --list) LIST=1; shift ;;
    *) echo "⛔ unknown option: $1" >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# 0. ROOT -- derived from the location of THIS file, never from the cwd
# ---------------------------------------------------------------------------
JMP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# Git Bash returns /c/Users/..., which R and Python cannot open. cygpath -m
# returns C:/Users/... with forward slashes, which all three accept.
command -v cygpath >/dev/null 2>&1 && JMP_ROOT="$(cygpath -m "$JMP_ROOT")"
export JMP_ROOT

if [ ! -d "$JMP_ROOT/01_estimation" ] || [ ! -d "$JMP_ROOT/data_inputs" ]; then
  echo "⛔ $JMP_ROOT does not look like the package root." >&2
  echo "   01_estimation/ and data_inputs/ are expected there." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 1. INTERPRETERS
# ---------------------------------------------------------------------------
RSCRIPT="${RSCRIPT:-}"
if [ -z "$RSCRIPT" ]; then
  if command -v Rscript >/dev/null 2>&1; then RSCRIPT="$(command -v Rscript)"
  else for c in "/c/Program Files/R"/R-*/bin/x64/Rscript.exe; do [ -x "$c" ] && RSCRIPT="$c"; done
  fi
fi
[ -n "$RSCRIPT" ] || { echo "⛔ Rscript not found. Set RSCRIPT=..." >&2; exit 1; }

PYTHON="${PYTHON:-}"
if [ -z "$PYTHON" ]; then
  for c in python python3 py; do command -v "$c" >/dev/null 2>&1 && { PYTHON="$c"; break; }; done
fi
[ -n "$PYTHON" ] || { echo "⛔ Python not found. Set PYTHON=..." >&2; exit 1; }

# ---------------------------------------------------------------------------
# 2. ENVIRONMENT -- SET EXPLICITLY, EVEN WHEN THE VALUE IS THE DEFAULT
#
# What is not set runs on a default, and a default appears in no log. Setting
# every variable here makes the log say under which configuration an output was
# produced.
# ---------------------------------------------------------------------------
export ERA5_DIR="${ERA5_DIR:-$JMP_ROOT/data_inputs/ERA}"
export ERA5_OUT="${ERA5_OUT:-$JMP_ROOT/intermediate_outputs}"
export OGR2OGR="${OGR2OGR:-C:/OSGeo4W/bin/ogr2ogr.exe}"   # binaire GDAL, ni R ni Python

# -- GE arm -----------------------------------------------------------------
# CAP_COLMASS is a SPECIFICATION CHOICE, not a technical setting: it caps the
# column mass of lambda0, which changes the counterfactual. Leaving it implicit
# would publish a cap that no one asked for.
export CAP_COLMASS="${CAP_COLMASS:-none}"
export AMENITY_SPEC="${AMENITY_SPEC:-baseline}"
export ANNUAL_DRF_SOURCE="${ANNUAL_DRF_SOURCE:-annual}"
export DUMP_OUTER="${DUMP_OUTER:-0}"

# RUNS_DIR tells emit_values.py to read the run grid from output/, where the
# solver writes. Its PRESENCE is what matters, not its value.
export RUNS_DIR="${RUNS_DIR:-output}"
export CAP_TAG="${CAP_TAG:-}"

# ---------------------------------------------------------------------------
# 3. THE STEPS -- named, so that --from can resume at one
# ---------------------------------------------------------------------------
ETAPES=(matrice poids climat_local 01a 01a_reemit 01f 01b 01c 01c_graine2 01c_reemit 01d crh1 crh2 crh3 crh4 crh5 crh6 crh7 01d_bis crhclamp crhfig 01e ct1 ct2 ct8 ct9 ct3 ct4 ct6 ct5 ct7 02a 02b mig_cross colmass calibdiag miglat ge masks maskruns maskdec band 02d releve_pa1 targets)

if [ "$LIST" = 1 ]; then
  printf "Steps (for --from and --until):\n"; for e in "${ETAPES[@]}"; do echo "  $e"; done; exit 0
fi

if [ -n "$FROM" ]; then
  ok=0; for e in "${ETAPES[@]}"; do [ "$e" = "$FROM" ] && ok=1; done
  [ "$ok" = 1 ] || { echo "⛔ unknown step: $FROM  (see --list)" >&2; exit 2; }
fi

## --until VALIDE (lot 54, 28/09/2026). Sans cette garde, un nom mal tape, `targets`
## (qui n'est pas une etape) ou une etape anterieure a --from n'etait JAMAIS atteint :
## la course allait au bout en silence, point d'arret compris.
valide_until () {
  [ -z "$UNTIL" ] && return 0
  local i=0 iu=-1 if_=0 e
  for e in "${ETAPES[@]}"; do
    [ "$e" = "$UNTIL" ] && iu=$i
    [ -n "$FROM" ] && [ "$e" = "$FROM" ] && if_=$i
    i=$((i + 1))
  done
  if [ "$iu" -lt 0 ] || [ "$UNTIL" = "targets" ]; then
    echo "⛔ --until: unknown step, or a step the run cannot stop after: $UNTIL  (see --list)" >&2; exit 2
  fi
  if [ "$iu" -lt "$if_" ]; then
    echo "⛔ --until $UNTIL comes before --from $FROM: the stop point would never be reached." >&2
    exit 2
  fi
}
valide_until

## L'ARRET, UNE SEULE FOIS (lot 54). Appele par etape() ET par le bloc GE, qui ne
## passe pas par etape() : `--until ge` n'arretait pas la course.
arret_si_until () {   # arret_si_until <nom>
  if [ -n "$UNTIL" ] && [ "$1" = "$UNTIL" ]; then
    echo "🛑 STOP POINT: --until $UNTIL reached; steps done:"
    for l in "${FAITES[@]}"; do echo "  ✔ $l"; done
    exit 0
  fi
}

ACTIF=1; [ -n "$FROM" ] && ACTIF=0
actif () {   # actif <nom> -> 0 si l'etape doit tourner
  [ "$1" = "$FROM" ] && ACTIF=1
  [ "$ACTIF" = 1 ]
}

## Lot 85 (03/10/2026) : le repertoire de la course est cree ICI, apres --list et la
## validation de --from/--until. Cree avant (l.110 jusqu'au lot 84), il laissait un
## logs/run_* vide a chaque --list et a chaque option refusee.
LOGDIR="$JMP_ROOT/logs/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"

# ---------------------------------------------------------------------------
# 4. RUN IDENTITY -- R, Python and package versions, written to the log
#
# The package pins nothing: there is no renv.lock and no requirements.txt. This
# log is therefore the only record of the environment that produced the outputs,
# and without it two runs cannot be compared.
# ---------------------------------------------------------------------------
{
  echo "=== RUN IDENTITY — $(date '+%Y-%m-%d %H:%M:%S') ==="
  echo "root    : $JMP_ROOT"
  echo "Rscript : $RSCRIPT"
  echo "Python  : $PYTHON"
  [ -n "$FROM" ] && echo "resume  : --from $FROM"
  # Lot 53 : le journal inscrit la REPRISE DELIBEREE (mode rapport des gardes gelees).
  if [ -f "$JMP_ROOT/tools/REPRISE_DELIBEREE.tsv" ]; then
    echo "DECLARED RESUMPTION (tools/REPRISE_DELIBEREE.tsv) -- non-regression guards in REPORT MODE:"
    grep -v '^#' "$JMP_ROOT/tools/REPRISE_DELIBEREE.tsv" | sed 's/^/   /'
  else
    echo "declared resumption: none -- non-regression guards are BLOCKING"
  fi
  echo
  "$RSCRIPT" --vanilla -e 'cat(R.version.string, "|", R.version$platform, "\n\n"); print(sessionInfo())'
  echo
  echo "--- R packages ---"
  "$RSCRIPT" --vanilla -e '
    p <- c("data.table","dplyr","tidyr","fixest","MASS","lfe","readxl","lubridate",
           "ggplot2","sf","terra","exactextractr","hdf5r","R.matlab","rnaturalearth",
           "countrycode","scales","Matrix","doParallel","foreach","splines")
    for (x in sort(p)) cat(sprintf("  %-16s %s\n", x,
      tryCatch(as.character(packageVersion(x)), error=function(e) "MISSING")))'
  echo
  echo "--- Python ---"
  "$PYTHON" -c "import sys; print(' ', sys.version.replace(chr(10),' '))"
  "$PYTHON" - <<'PY'
import importlib
for m in ('numpy','pandas','h5py','scipy','rasterio','pyreadr'):
    try:
        mod = importlib.import_module(m)
        print(f"  {m:12s} {getattr(mod,'__version__','?')}")
    except Exception:
        print(f"  {m:12s} MISSING")
PY
  echo
  echo "--- ogr2ogr (GDAL, system dependency) ---"
  if [ -x "$OGR2OGR" ] || command -v "$OGR2OGR" >/dev/null 2>&1; then
    "$OGR2OGR" --version 2>&1 | head -1
  else
    echo "  MISSING: $OGR2OGR  (no longer required by the chain)"
  fi
  echo
  echo "--- environment variables set ---"
  for v in JMP_ROOT ERA5_DIR ERA5_OUT OGR2OGR CAP_COLMASS AMENITY_SPEC \
           ANNUAL_DRF_SOURCE DUMP_OUTER RUNS_DIR CAP_TAG; do
    eval "printf '  %-20s %s\n' $v \"\$$v\""
  done
} 2>&1 | tee "$LOGDIR/00_session.log"

echo
echo "============================================================"
echo " FULL CHAIN — logs: $LOGDIR"
echo "============================================================"
echo

# ---------------------------------------------------------------------------
# 5. INPUT CHECK -- up front, before the first step
#
# Without it, a missing input surfaces only at the fread that needs it, under a
# message that names neither the file, nor the step that should have produced
# it, nor why it cannot be regenerated locally.
# ---------------------------------------------------------------------------
## Lot 84 (03/10/2026) : d'abord tools/check_data.py -- presence ET md5 de chaque intrant
## de DATA.md (environ 30 Go lus). CHECK_DATA=presence le reduit a la presence (rapide).
echo "--- inputs of DATA.md: tools/check_data.py ---"
set +e
if [ "${CHECK_DATA:-md5}" = "presence" ]; then
  "$PYTHON" "$JMP_ROOT/tools/check_data.py" --presence
else
  "$PYTHON" "$JMP_ROOT/tools/check_data.py"
fi
rc_data=$?
set -e
if [ "$rc_data" -ne 0 ]; then
  echo "⛔ tools/check_data.py: code $rc_data. STOPPED BEFORE THE FIRST STEP (see DATA.md)."
  exit 1
fi
echo
MANQUE=0
## Lot 85 (03/10/2026) : la taille d'un FICHIER present est lue sur le disque, jamais
## ecrite a la main -- les tailles de tas_ssp2_{daily,annual}.csv affichees ici etaient
## restees celles d'avant la regeneration du 02/09. <volume> ne decrit plus que le reste.
verifier () {   # verifier <chemin> <volume> <producteur> <lecteurs>
  if [ -e "$JMP_ROOT/$1" ]; then
    local taille=""
    [ -f "$JMP_ROOT/$1" ] && taille="$(stat -c %s "$JMP_ROOT/$1") o"
    printf "  ok       %-50s %s %s\n" "$1" "$taille" "$2"
  else
    printf "  MISSING  %-50s %s\n" "$1" "$2"
    printf "           produced by : %s\n" "$3"
    printf "           read by     : %s\n" "$4"
    MANQUE=$((MANQUE + 1))
  fi
}

echo "--- outputs of the SERVER steps (supplied; cannot be regenerated locally) ---"
## Lot 51 : ERA_tas_annual.csv n'est plus FOURNI -- l'etape climat_local le produit
## depuis la grille ci-dessous et poids_par_pays.csv. La frontiere serveur s'arrete a
## la GRILLE (passage des 97 Go de .nc aux cellules).
verifier intermediate_outputs/ERA_tas_distribution_grid_annual.csv "" \
         "00_era5_distribution.R  [SERVER]"            "climat_local, 01d"
verifier intermediate_outputs/CMIP6_tas.csv                        "" \
         "00_temperature_projections_JMP.R  [SERVER]"  "01d, 01e"
verifier intermediate_outputs/tas_ssp2_daily.csv                   "" \
         "00_temperature_projections_JMP.R  [SERVER]"  "02c_model/main.py"
verifier intermediate_outputs/tas_ssp2_annual.csv                  "" \
         "00_temperature_projections_JMP.R  [SERVER]"  "02c_model/main.py"

echo
echo "--- public inputs ---"
verifier data_inputs/BACI_CEPII_trade/data              "12 GB, 28 files"     "CEPII"   "01a"
verifier data_inputs/BACI_CEPII_trade/WITS_Concordance_H0_to_I3.CSV "549,733 B — UPPER-CASE .CSV" "WITS" "01a"
verifier data_inputs/GLORIA/data_gloria/97              "3.0 GB, 25 files"    "GLORIA"  "01a"
verifier data_inputs/gravity_CEPII/Gravity_V202211.csv  "1.2 GB"              "CEPII"   "01a, 01b"
verifier data_inputs/GHSL_POP/tif                       "2.5 GB, 7 tif"       "JRC GHSL" "matrice (00c), 02a"
verifier data_inputs/gadm36_shp/gadm_simplified.shp     "244 MB"              "GADM 3.6 simplified" "00_era5, 02a"
verifier data_inputs/gadm36_shp/gadm_410.gpkg           "2.6 GB"              "GADM 4.1 (unzipped from gadm_410-gpkg.zip)" "matrice (00c)"
verifier tools/anomalies_rattachement_connues.csv        "31,253 pixels"       "construit_anomalies_connues.py (author's tool)" "matrice (00c)"
verifier data_inputs/EGGW_Replication_Package.zip       "2.2 GB — NESTED ARCHIVE" "Zenodo 10.5281/zenodo.7681424" "01d"
verifier data_inputs/abel_cohen_migration/bilat_mig.csv "9.3 MB" "Abel & Cohen, Figshare 10.6084/m9.figshare.7731233" "01b, 02a, 02b"
verifier data_inputs/Kummu/rast_gdpTot_1990_2022_30arcmin.tif "20 MB"         "Kummu et al." "02a"
verifier data_inputs/Kummu/raster_netMgr_2010_2014_5yrSum.tif "20 MB"         "Kummu et al." "02b"

if [ "$MANQUE" -gt 0 ]; then
  echo
  echo "⛔ $MANQUE input(s) missing. STOPPED BEFORE THE FIRST STEP."
  echo "   Nothing has run: a chain launched on a missing input produces"
  echo "   outputs that look valid."
  exit 1
fi
echo
echo "  -> required inputs present."
echo

# ---------------------------------------------------------------------------
# 6. THE STEP RUNNER
# ---------------------------------------------------------------------------
declare -a FAITES=() SAUTEES=()

etape () {   # etape <nom> <interpreteur> <script> <produit>
  local nom="$1" interp="$2" script="$3" quoi="$4"
  if ! actif "$nom"; then
    SAUTEES+=("$nom (before --from $FROM)")
    printf "  ⏭  %-14s skipped (--from %s)\n" "$nom" "$FROM"
    return 0
  fi
  local t0 t1 d rc
  echo "------------------------------------------------------------"
  echo "[$nom] $(basename "$script")"
  echo "     output  : $quoi"
  echo "     start   : $(date '+%H:%M:%S')"
  t0=$(date +%s)
  set +e
  "$interp" "$JMP_ROOT/$script" 2>&1 | tee "$LOGDIR/${nom}.log"
  rc=${PIPESTATUS[0]}
  set -e
  t1=$(date +%s); d=$((t1 - t0))
  if [ "$rc" -ne 0 ]; then
    echo "     FAILED  : $(date '+%H:%M:%S')  after $((d/60)) min $((d%60)) s  (code $rc)"
    echo
    echo "⛔ STOPPED at step '$nom'. Full log: $LOGDIR/${nom}.log"
    echo "   After the fix, resume with:  ./run_all.sh --from $nom"
    exit "$rc"
  fi
  echo "     end     : $(date '+%H:%M:%S')  duration $((d/60)) min $((d%60)) s"
  echo
  FAITES+=("$(printf '%-14s %3d min %02d s' "$nom" $((d/60)) $((d%60)))")
  arret_si_until "$nom"
}

# 00_parameters.R is not a step: it is the parameter library, sourced by every
# script, and it derives the root from its own location.

## Lot 51 (28/09/2026) : `pop_weights` (00_population_weights.R, ancien
## lengths_country.csv) est RETIRE. Il comptait la population ENTIERE de chaque
## cellule pour chaque pays qui la touche (+23 % au monde, RECAP 3.231). Trois
## etapes le remplacent, chacune gardee par ses temoins AVANT toute ecriture :
etape matrice "$PYTHON" 01_estimation/00c_matrice_population.py \
      "matrice_population_cellule_pays.csv — GHSL 30'' x unsimplified gadm_410, one-step attachment rule (~50 min)"
etape poids "$PYTHON" 01_estimation/00d_poids_population.py \
      "poids_par_cellule.csv (part_de_la_cellule) + poids_par_pays.csv (part_du_pays)"
etape climat_local "$PYTHON" 01_estimation/00e_era5_pays_local.py \
      "ERA_tas_annual.csv — shipped ERA5 grid, weighted by part_du_pays"

# ---------------------------------------------------------------------------
#  [SERVER STEP -- 00_era5_distribution.R]   NOT CALLED
#
#  WHY   It reads about 15.6 GB of daily ERA5 netCDF at 0.25 degrees, which are
#        not shipped: data_inputs/ERA/ is empty and the README points to the
#        Copernicus CDS. The run takes roughly an hour on eight cores.
#
#  ITS OUTPUTS ARE SUPPLIED, AND CHECKED IN SECTION 5:
#        intermediate_outputs/ERA_tas_annual.csv
#        intermediate_outputs/ERA_tas_distribution_grid_annual.csv
#        That check is what replaces the step. Without it the chain would fail
#        at the fread in 01a without naming the cause.
#
#  To re-run it on the server, with the netCDF files in place:
#    ERA5_DIR=<netcdf dir> ERA5_OUT=$JMP_ROOT/intermediate_outputs #    Rscript 01_estimation/00_era5_distribution.R
# ---------------------------------------------------------------------------

etape 01a "$RSCRIPT" 01_estimation/01a_productivity_drf.R \
      "productivity damage functions — BACI_raw_pairs_annual.rds, damage_prod*_bin.csv, 5 .tex, 3 .pdf"
# 01a_reemit et 01f : entres le 25/09/2026 (lot 34). L'audit de graphe
# (tools/audit_graphe.py) les trouvait appeles par personne ; 01f produit
# oh_surviving_share.csv, sans lequel emit_values_empirics.py s'arrete en 02d.
# Tous deux lisent des sorties de 01a et rien d'autre de la chaine.
etape 01a_reemit "$RSCRIPT" 01_estimation/01a_reemit_tab_absorption.R \
      "tab_fe_absorption.tex, tab_ctrl_absorption.tex re-emitted from their CSV files"
etape 01f "$RSCRIPT" 01_estimation/01f_oh_surviving_share.R \
      "oh_surviving_share.csv — share of the climate variation that survives (69%)"
etape 01b "$RSCRIPT" 01_estimation/01b_amenity_drf.R \
      "amenity damage functions — step2_migration.rds, damage_amenity*_bin.csv, 4 .tex, 2 .pdf"
etape 01c "$RSCRIPT" 01_estimation/01c_robustness.R \
      "robustness — tab_oster_delta.tex, tab_oh_*.tex"
# ---------------------------------------------------------------------------
# 01c_graine2 (lot 68, 30/09/2026) -- 01c REJOUE A UNE SECONDE GRAINE, a cote.
# La phrase de B.2 sur l'accord de deux graines independantes n'avait pas de
# source dans le paquet. Cette etape la fournit : meme script, meme donnees, la
# graine du bras productivite journalier surchargee (JMP_MC_SEED_PROD_DAILY),
# TOUTES les sorties dans output_01c_graine2/ (JMP_OUTPUT_DIR ; 00_parameters.R
# refuse la graine sans la redirection). output/ et intermediate_outputs/ ne
# sont pas ecrits. emit_values_empirics.py en tire EMPohSeed* en 02d.
# Reference : output_01c_graine2_20260930/ (rejeu manuel du lot 67), copiee a
# l'octet dans 01_estimation/reference_repro_01c_graine2/ (lot 85f) : la garde lit la copie.
# ---------------------------------------------------------------------------
SEED_OH_2=20260930
_old_outdir="${JMP_OUTPUT_DIR-}"; _old_seed="${JMP_MC_SEED_PROD_DAILY-}"
if actif 01c_graine2; then mkdir -p "$JMP_ROOT/output_01c_graine2"; fi
export JMP_OUTPUT_DIR="$JMP_ROOT/output_01c_graine2" JMP_MC_SEED_PROD_DAILY="$SEED_OH_2"
# Lot 74 : l'etape passe par tools/run_01c_graine2.sh, qui enchaine 01c et la garde
# de reproductibilite NUMERIQUE (1e-8) contre reference_repro_01c_graine2/ -- elle
# remplace la boucle `cmp` octet a octet tapee a la main le 30/09.
export RSCRIPT PYTHON JMP_ROOT
etape 01c_graine2 bash tools/run_01c_graine2.sh \
      "01c at seed $SEED_OH_2 — output_01c_graine2/oh_hotslope_schemes.csv (+ tables, bases), check at 1e-8"
if [ -n "$_old_outdir" ]; then export JMP_OUTPUT_DIR="$_old_outdir"; else unset JMP_OUTPUT_DIR; fi
if [ -n "$_old_seed" ]; then export JMP_MC_SEED_PROD_DAILY="$_old_seed"; else unset JMP_MC_SEED_PROD_DAILY; fi
etape 01c_reemit "$RSCRIPT" 01_estimation/01c_reemit_tab_oh_bounds.R \
      "tab_oh_bounds.tex re-emitted from oh_scheme_effects.csv + oh_hotslope_schemes.csv"

# 01d loads trmult_reduced.mat (2.1 GB) from the CRH archive, so this step is
# markedly slower than the others.
etape 01d "$RSCRIPT" 01_estimation/01d_crh_triangulation.R \
      "CRH triangulation — intermediate_outputs/crh/*, fig_crh_*.pdf, tab_crh_northern.tex"

# ---------------------------------------------------------------------------
# 6bis. THE CRH INFERENCE CHAIN -- 03_inference/step1..step7
#
# It takes the outputs 01d has just written into intermediate_outputs/crh/ and
# returns the 95 per cent intervals the manuscript quotes. The analytical
# interval is the published one: step3 builds the covariance, step4 draws from
# it. step5 is a bootstrap WITNESS, not the number.
#
# THE AMENITY BOOTSTRAP IS NOT RUN, and this is a decision, not an omission. It
# was measured at 4.84 s per iteration, so 80.7 minutes for B = 1000, against a
# budget of 75 minutes fixed BEFORE the measurement. Raising a budget once the
# overrun is known is choosing the threshold from the result. It also feeds no
# macro: the emitter reads the productivity channel at NORTH>=45N only.
#
# 01d RUNS TWICE, and the loop is real rather than an ordering quirk:
#   pass 1  01d writes intermediate_outputs/crh/*, the inputs of step1;
#   chain   step1..step7 write output/crh_inference/*;
#   pass 2  01d reads crh_inference_deliverable.csv and only then can write
#           tab_crh_northern.tex, which carries the intervals.
# On pass 1 the table is skipped with a message; the script does not fail.
#
# COST, measured end to end on 3 September: step1 16 s, step2 7 s, step3 4 min 10,
# step4 3 s, step5 14 min 09, step6 under a second, step7 2 min 38 -- 21 min 29 for
# the seven. step5, the bootstrap, dominates. An earlier estimate of two hours came
# from consecutive file timestamps in an interactive session, which measure the gap
# between two gestures of the operator rather than any runtime.
#
# THE THREE REPRODUCTION CHECKS WILL REPORT A NON-ZERO GAP, and that is what
# they are for. They compare the chain against frozen witnesses in
# 01_estimation/reference_repro_crh/, which no step here writes. Measured on the
# current tree: about 8.8e-02 on the productivity bins, about 2.2e-04 on the
# amenity ones. Pointing them at intermediate_outputs/crh/ would compare 01d
# with itself and turn all three green having verified nothing.
# ---------------------------------------------------------------------------

etape crh1 "$RSCRIPT" 03_inference/step1_gate.R \
      "step1_frame.rds — three reproduction checks per channel"
etape crh2 "$RSCRIPT" 03_inference/step2_maps.R \
      "step2_maps.rds — cell-by-bin maps"
etape crh3 "$RSCRIPT" 03_inference/step3_vcov.R \
      "step3_vcov.rds — Conley and country covariances"
etape crh4 "$RSCRIPT" 03_inference/step4_mc.R \
      "crh_inference_raw.csv — analytical intervals, the published numbers"
etape crh5 "$RSCRIPT" 03_inference/step5_boot.R \
      "crh_bootstrap_raw.csv — bootstrap witness, productivity channel"
etape crh6 "$RSCRIPT" 03_inference/step6_emit.R \
      "crh_inference_deliverable.csv — read by 01d (second pass) and emit_values_empirics.py"
etape crh7 "$RSCRIPT" 03_inference/step7_period_sensitivity.R \
      "crh_period_sensitivity.csv — sensitivity to the periods"

# Second pass. Everything before tab_crh_northern.tex is recomputed identically;
# the table is the only thing this pass adds, and it could not exist before the
# chain had run.
etape 01d_bis "$RSCRIPT" 01_estimation/01d_crh_triangulation.R \
      "tab_crh_northern.tex — second pass, with the intervals of the chain"

# LES DEUX EVALUATIONS SUR ARCHIVES CRH (lot 10, 18/09/2026).
# Aucune des deux n'estime quoi que ce soit : elles LISENT les CSV que 01d vient
# d'ecrire dans intermediate_outputs/crh/ et les evaluent. Elles passent donc
# apres 01d et peuvent passer avant tout le reste ; elles doivent passer avant
# 02d, qui lance emit_values_empirics.py, seul consommateur de leurs sorties.
# Chacune porte son propre controle de non-regression contre un artefact de
# 01d -- l'agregat nordique non clampe pour l'une, crh_daily_response_peaks.csv
# pour l'autre -- de sorte qu'une evaluation qui aurait derive de son pipeline
# d'origine ne produit rien plutot que de produire autre chose.
etape crhclamp "$RSCRIPT" 01_estimation/01d_crh_clamp_sensitivity.R \
      "crh_clamp_sensitivity.csv, crh_north_daymass.csv — 3 macros EMPcrhClamp*"
etape crhfig "$RSCRIPT" 01_estimation/01d_fig_crh_responses.R \
      "fig_crh_responses.pdf (4 panels) + crh_daily_response_peaks.csv"

etape 01e "$RSCRIPT" 01_estimation/01e_diagnostics.R \
      "diagnostics — tab_dbar_cross.tex, map_*.pdf, table_eta.csv"

# ---------------------------------------------------------------------------
# 6ter. THE ANNUAL/DAILY CONTRAST CHAIN -- 04_contrast/ct1..ct4
#
# LE GRAPHE REEL, mesure le 07/09 -- il n'est PAS une chaine lineaire :
#     ct1 --> ct3            (ct3 ne lit AUCUNE sortie de ct2)
#     ct1 --> ct2 --> ct4    (ct4 refuse si l'un des deux manque)
# ct1 ecrit output/contribA/test1_delta_north.csv ; ct2 y lit les scores, ct3
# les bornes de support du bras annuel et l'increment de reference, ct4 la
# grille entiere PLUS tache1_diff_bras.csv de ct2.
# L'ordre numerique ci-dessous respecte les deux branches, donc il convient a
# l'orchestrateur. Mais pour un rejeu a la main dont on attend le resultat de
# ct3 -- les dix macros \EMPgap* -- ct3 se lance DES QUE ct1 est fini, sans
# attendre ct2, qui coute autant que ct1 (meme grille, deux bras par cellule).
#
# ct3 recalcule l empreinte de provenance des six intrants et REFUSE si elle
# differe de celle que ct1 a inscrite en colonne. Le danger vise est
# l appariement croise : combiner des courbes d un millesime avec des bornes
# d un autre donne une decomposition qui telescope exactement, qu aucun
# controle de somme ne detecterait.
#
# COUT MESURE le 2026-09-07, machine locale, SUR LA GRILLE LARGE (24 specs) :
#   ct1  16 min 46 s   96 modeles (8 formes x 3 winsor. x 2 canaux x 2 bras)
#   ct2  21 min 20 s   meme grille, mais les DEUX bras par cellule
#   ct3      20 s      ct4       1 s
# ⚠ LA GRILLE A CHANGE LE 2026-09-18 : GRID_SPECS (01_estimation/00_parameters.R)
# porte QUINZE specifications, soit 60 modeles pour ct1 et 30 paires pour ct2.
# Les durees ci-dessus sont donc des MAJORANTS ; le compte, lui, n'est ecrit
# nulle part ici -- ct1 et ct2 lisent la constante et l'annoncent au demarrage.
# ⚠ Ce n'est PAS ct1 qui domine, contrairement a ce qu'on attendrait du nombre
# de modeles : ct2 estime autant et calcule en plus la covariance entre bras.
# ---------------------------------------------------------------------------
etape ct1 "$RSCRIPT" 04_contrast/ct1_delta_north_grid.R \
      "GRID_SPECS grid x 4 (channel x arm) — output/contribA/test1_delta_north.csv"
etape ct2 "$RSCRIPT" 04_contrast/ct2_between_arm_test.R \
      "difference between arms — output/contribA/tache1_diff_bras.csv"
# ct8 et ct9 (lot 78, 02/10/2026) : le biais de binning de Jones, Moscona, Olken,
# von Dessauer (2026), transpose du bac a sable du lot 77 (_sandbox_binning_20261001).
# ct8 re-estime les specifications de reference sous quatre structures d'effets
# fixes (base, between, pp, pp5) ; son temoin BLOQUANT exige que base retombe sur
# test1_delta_north.csv (ct1) et tache1_diff_bras.csv (ct2) a 1e-12, d'ou sa place
# APRES ct2. ct9 mesure les ingredients du biais et la part intra-cellule du design
# CRH ; il lit l'effectif de base dans la sortie de ct8, d'ou sa place APRES ct8.
# Duree : ct8 5 min 28 s au bac a sable ; ct9 non mesure (quelques minutes attendues).
etape ct8 "$RSCRIPT" 04_contrast/ct8_binning_variants.R \
      "contribA/binning_variants.csv, binning_contrasts.csv, ct8_controls.csv — four fixed-effect structures"
etape ct9 "$RSCRIPT" 04_contrast/ct9_binning_mechanism.R \
      "contribA/binning_mechanism.csv, crh_within_cell.csv, ct9_controls.csv — ingredients of the bias, CRH within-cell share"
etape ct3 "$PYTHON"  04_contrast/ct3_weighted_curve_gap.py \
      "A/B decomposition — pe_weighted_gap_north{,_scalars}.csv (10 EMPgap macros)"
etape ct4 "$RSCRIPT" 04_contrast/ct4_emit_tab_contribA_grid.R \
      "tab_contrast_grid.tex, tab_robustness_prod_joint.tex, fig_contrast_grid.pdf"
# ct6 : la Figure 1 du manuscrit. Elle vient APRES ct1 (grille en fond) et apres
# ct4 (le compte de specifications qu'elle asserte). L'appel correspondant a ete
# RETIRE de 01a le 17/09/2026 : un fichier, un producteur.
etape ct6 "$RSCRIPT" 04_contrast/ct6_fig_drf_paper.R \
      "Figure 1 — MC_prod_{daily,annual}_fan_paper.pdf, grid in the background"
# ct5 et ct7 : entres le 25/09/2026 (lot 34).
# ct5 CLASSE les tirages Monte-Carlo par le changement nord et ecrit les colonnes
# rq05/rq50/rq95 et rb01..rbN dans MC_damage_prod_daily_quantiles.csv et
# MC_damage_amenity_daily_quantiles_NULL.csv, plus ge_runs_band_armB.csv. Il doit
# donc passer AVANT 'ge' : les runs de bande lisent ces colonnes, et emit_values.py
# s'arrete sans ge_runs_band_armB.csv. ARMB_N vient de 00_parameters.R (ci-dessous).
# ⚠ ct5 ne reproduit ses sorties qu'AU DERNIER BIT (ecart <= 1e-16 mesure le
# 25/09) : sur un arbre deja peuple, le relancer seul change le md5 des
# dose-reponses et 02i refuse (drf_md5). Le relancer, c'est relancer tout l'aval.
ARMB_N="$(grep -E '^ARMB_N[[:space:]]*<-' "$JMP_ROOT/01_estimation/00_parameters.R" \
          | head -1 | sed -E 's/.*<-[[:space:]]*([0-9]+)L?.*/\1/')"
case "$ARMB_N" in
  ''|*[!0-9]*) echo "⛔ ARMB_N unreadable in 01_estimation/00_parameters.R" >&2; exit 1 ;;
esac
export ARMB_N_PILOTE="$ARMB_N"
etape ct5 "$RSCRIPT" 04_contrast/ct5_rank_mc_draws.R \
      "rq*/rb* columns of the MC_damage_* tables, ge_runs_band_armB.csv (ARMB_N=$ARMB_N)"
etape ct7 "$RSCRIPT" 04_contrast/ct7_cold_tail_variants.R \
      "contribA/ct7_cold_tail_variants.csv, ct7_controls.csv — cold-tail variants"
etape 02a "$RSCRIPT" 02_GE/02a_calibration_at_1deg.r \
      "one-degree calibration — lambda.h5, population0/wage0/productivity0, 2 .pdf"
etape 02b "$RSCRIPT" 02_GE/02b_mig_flows.R \
      "mig_flows.h5 — migration matrix 17,048^2, written in blocks of 200 rows"

# Lot 85x (11/10/2026) : T1, la masse transfrontaliere du bloc interne de mu_0. Lit
# les RDS et le temoin que 02b vient d'ecrire, et la masse STEP 11 de son journal ;
# ecrit intermediate_outputs/mig_internal_cross.csv, que 02d (emit_values.py) lit
# pour \GEmigInternalCrossM et \GEmigInternalCrossShare. Hors chaine jusqu'au 85w :
# le clone de test les a emises a \PROVISOIRE.
etape mig_cross bash tools/run_mig_cross.sh \
      "mig_internal_cross.csv — GEmigInternalCrossM, GEmigInternalCrossShare"

# Le diagnostic de colmass lit les DEUX matrices ci-dessus et rien d'autre : il
# passe donc apres 02b et peut passer avant 'ge'. Il ne depend d'aucun welfare,
# et ses quatre macros decrivent les INTRANTS, pas un run.
etape colmass "$PYTHON" 02_GE/cap_colmass_diagnostic.py \
      "cap_colmass_diagnostic.csv — 4 CALIB*Colmass macros (app:ge_robustness)"

# Statistics of the calibrated initial state: diagonal of lambda0 (and the mean
# import share), sales/income, largest net migration senders and receivers. It
# reads the input matrices only, so it passes here, before 'ge'. It was run by
# hand and called by no step until 25/09/2026, so a clean replication never
# regenerated ge_calib_diagnostic.csv, which feeds GElamDiag*, GEclear*,
# GEimportShareMean and GEmig*.
etape calibdiag "$PYTHON" 02_GE/ge_calibration_diagnostic.py \
      "ge_calib_diagnostic.csv — GElamDiag*, GEimportShareMean, GEclear*, GEmig*"

# Profil latitudinal de la migration nette. Lit la matrice BRUTE ecrite par 02b
# et ne depend pas du solveur : passe donc avant 'ge'.
etape miglat "$PYTHON" 02_GE/02f_fig_mig_latitude.py \
      "fig_M1_net_by_latitude.png + migration_net_by_latitude.csv"

# ---------------------------------------------------------------------------
#  [SERVER STEP -- 00_temperature_projections_JMP.R]   NOT CALLED
#
#  WHY   It is the single producer of every future temperature used downstream.
#        Its input is ISIMIP3b bias-adjusted W5E5 netCDF, SSP2-4.5, daily tas at
#        0.5 degrees, four models over 2015-2100. Those files are not shipped:
#        data_inputs/CMIP6/ is empty.
#
#  ITS OUTPUTS ARE SUPPLIED, AND CHECKED IN SECTION 5:
#        intermediate_outputs/CMIP6_tas.csv
#        intermediate_outputs/tas_ssp2_daily.csv
#        intermediate_outputs/tas_ssp2_annual.csv
#
#  The outputs are keyed by the cell id of identifier_xy_final.csv. Re-running
#  02a changes that numbering, so these three files must be regenerated with it;
#  otherwise every cell receives another cell temperature and nothing reports it.
#
#  To re-run:  Rscript 02_GE/00_temperature_projections_JMP.R #                <input_dir> <identifier_xy_final.csv> <out_dir>
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 7. THE GE GRID -- 78 runs since 25/09/2026 (eleven before; see GE_GRID below)
#
# Scenario 1 is the no-climate-change baseline: main.py forces both damage
# vectors to one, and the quantile never enters the computation. sc1_q05,
# sc1_q50 and sc1_q95 would therefore be byte-identical. Running all three would
# report nine runs where only seven carry information.
#
# THE EIGHTH RUN IS sc4, AND IT EXISTS FOR THE DECOMPOSITION. It is scenario 3
# read off the DAILY curve (ANNUAL_DRF_SOURCE=dailycurve), which main.py labels
# sc4: the daily dose-response evaluated at the annual mean. Without it the
# evaluation/shape decomposition has no counterfactual, and the two macros
# GEnorthChanEval and GEnorthChanShape cannot be emitted. Only q50 is run: the
# decomposition bears on the headline arm alone.
#
# RUNS NINE TO ELEVEN ARE THE g = 0 AMENITY ARM, the Manski lower bound. They
# re-solve sc1, sc2 and sc3 with AMENITY_SPEC=null, that is with the amenity
# dose-response set to zero everywhere, and their outputs carry a _NULLamen
# suffix. Without them emit_values.py leaves HAVE_NULL false: seven macros stay
# \PROVISOIRE and tab_ge_sensitivity.tex and tab_ge_north_sens.tex are not
# written at all. The question they answer is whether the Northern reversal
# survives switching the amenity channel off.
#
# THREE RUNS, NOT ONE. The emitter reads the suffix on all three scenarios before
# it forms a margin, so a single sc2 under the tag would leave HAVE_NULL false and
# change nothing. sc1 carries no climate damage and its result does not depend on
# the amenity table, but the file must exist under the tagged name.
#
# EACH RUN APPENDS A LINE TO output/ge_runs_manifest.csv. That manifest is what
# lets a consumer prove three welfare matrices come from ONE grid under ONE set
# of settings. decompose_metric_channels.py refuses to run without it, because a
# triplet assembled from runs that do not belong together reverses the sign of
# both channels.
# ---------------------------------------------------------------------------
if actif ge; then
  GE_MAIN="$JMP_ROOT/02_GE/02c_model/main.py"
  GE_MANIFEST="$JMP_ROOT/output/ge_runs_manifest.csv"
  RUN_ID="$(basename "$LOGDIR")"
  # INPUTS_MD5 fingerprints the input state the grid is about to consume. It is
  # what a consumer must compare, not the run identifier: a resumed grid gets a
  # new RUN_ID while reading exactly the same files, and two runs that share
  # these inputs and these settings are comparable whether or not one invocation
  # produced both.
  # Lot 53 : MEME RECETTE que run_one_ge.sh et 02h (neuf fichiers + la ligne
  # SCALE_CONVENTION). Sans cette ligne, 02h refusait en `maskruns` -- apres les
  # treize heures de la grille. Voir 02_GE/lib_empreinte.sh.
  source "$JMP_ROOT/02_GE/lib_empreinte.sh"
  INPUTS_MD5="$(empreinte_intrants)" || { echo "⛔ cannot fingerprint the GE inputs"; exit 3; }
  echo "file,md5,run_id,scenario,quantile,cap_colmass,amenity_spec,annual_drf_source,inputs_md5,seconds,rc,timestamp" \
    > "$GE_MANIFEST"
  # scenario | quantile | ANNUAL_DRF_SOURCE | AMENITY_SPEC for this run
  #
  # LA GRILLE COMPLETE depuis le 25/09/2026 (lot 34) : les onze runs d'origine,
  # plus les 67 que seuls des pilotes de lot produisaient (run_lot25.sh) et
  # qu'aucun reproducteur n'aurait eus :
  #   sc4 q50 dailycurve a f_a = 0   \GEnorthChanEval, \GEnorthChanShape (§3.3)
  #   sc2/sc3 x rq05/rq50/rq95       les tirages classes (bande classee)
  #   sc2 x rb01..rbN, f_a = 0       la bande : Tableau 15, Figure 4, les parts
  #                                  de tirages du Tableau 1
  # N = ARMB_N, lu dans 00_parameters.R plus haut. 78 runs a N = 60.
  GE_GRID=("1 q50 annual baseline" "2 q05 annual baseline"
           "2 q50 annual baseline" "2 q95 annual baseline"
           "3 q05 annual baseline" "3 q50 annual baseline"
           "3 q95 annual baseline" "3 q50 dailycurve baseline"
           "1 q50 annual null" "2 q50 annual null" "3 q50 annual null"
           "3 q50 dailycurve null")
  for Q in rq05 rq50 rq95; do for SC in 2 3; do GE_GRID+=("$SC $Q annual null"); done; done
  for N in $(seq -w 1 "$ARMB_N"); do GE_GRID+=("2 rb$N annual null"); done
  echo "------------------------------------------------------------"
  echo "[ge] 02c_model/main.py — ${#GE_GRID[@]} runs (≈ 10 min each, measured)"
  echo "     CAP_COLMASS=$CAP_COLMASS  AMENITY_SPEC=$AMENITY_SPEC  ANNUAL_DRF_SOURCE=$ANNUAL_DRF_SOURCE"
  echo "     manifest : output/ge_runs_manifest.csv  (run_id=$RUN_ID)"
  t0g=$(date +%s); nge=0
  for combo in "${GE_GRID[@]}"; do
    set -- $combo; SC="$1"; Q="$2"; ADRF="$3"; AMEN="$4"
    case "$AMEN" in null) ATAG="_NULLamen" ;; bell) ATAG="_BELLamen" ;;
                    bell_cap) ATAG="_BELLCAPamen" ;; *) ATAG="" ;; esac
    # ⚠ L'ORDRE DES SUFFIXES EST CELUI DE main.py : _sc<label>_<q> + amenite +
    # cap + _dailycurve. L'ancien TAG mettait _dailycurve AVANT l'amenite : juste
    # a f_a interieur (amenite vide), faux a f_a = 0 -- le manifeste aurait
    # cherche un fichier qui n'existe pas et n'aurait pas ecrit la ligne.
    if [ "$ADRF" = "dailycurve" ]; then TAG="sc4_${Q}${ATAG}_dailycurve"
    else TAG="sc${SC}_${Q}${ATAG}"; fi
    echo "     ▶ $TAG  $(date '+%H:%M:%S')"
    t0=$(date +%s)
    set +e
    SCENARIO="$SC" QUANTILE="$Q" ANNUAL_DRF_SOURCE="$ADRF" AMENITY_SPEC="$AMEN" \
      "$PYTHON" -u "$GE_MAIN" \
      2>&1 | tee "$LOGDIR/ge_${TAG}.log"     # ONE log PER RUN, never appended
    rc=${PIPESTATUS[0]}
    set -e
    t1=$(date +%s); d=$((t1 - t0))
    if [ "$rc" -ne 0 ]; then
      echo "     FAILED $TAG (code $rc) after $((d/60)) min $((d%60)) s"
      echo "⛔ STOPPED. Log: $LOGDIR/ge_${TAG}.log"
      echo "   Resume with:  ./run_all.sh --from ge"
      exit "$rc"
    fi
    # The manifest line is written only after a run that returned 0, so a failed
    # run leaves no provenance behind for a consumer to trust.
    WF="welfare_${TAG}.csv"
    if [ -f "$JMP_ROOT/output/$WF" ]; then
      echo "$WF,$(md5sum "$JMP_ROOT/output/$WF" | cut -d' ' -f1),$RUN_ID,$SC,$Q,$CAP_COLMASS,$AMEN,$ADRF,$INPUTS_MD5,$d,$rc,$(date '+%Y-%m-%dT%H:%M:%S')" \
        >> "$GE_MANIFEST"
    else
      echo "     ⚠ $WF not found: no line in the manifest."
    fi
    echo "       ok  $((d/60)) min $((d%60)) s"
    nge=$((nge+1))
  done
  t1g=$(date +%s); dg=$((t1g - t0g))
  FAITES+=("$(printf '%-14s %3d min %02d s  (%d runs)' "ge" $((dg/60)) $((dg%60)) "$nge")")
  echo
  arret_si_until ge
else
  SAUTEES+=("ge (before --from $FROM)"); echo "  ⏭  ge             skipped (--from $FROM)"
fi

# ---------------------------------------------------------------------------
# 7bis. THE SOURCE DECOMPOSITION -- three steps, four GE runs
#
# Where does the Northern daily extra come from: from warming in the North, or
# from warming elsewhere reaching the North through goods and migration? Two
# masked climate files answer it without touching the solver. Outside its kept
# cells, a mask copies the 2015 exposure into every later period; damage is
# computed on the change between periods, so warming there is exactly zero.
# General equilibrium is not switched off: trade, migration, wages and
# congestion operate in every run.
#
#   masks     02g  builds intermediate_outputs/exposure_masks/{north,south}/
#   maskruns  02h  daily and annual arm per mask, f_a = 0, q50, CAP none: FOUR
#                  runs, one draw, no band. ~10 min each, 14 GB free required.
#                  Own manifest (exposure_masks/runs_manifest.csv): the headline
#                  manifest is indexed by file name and these runs share the
#                  headline names. 02h refuses to run unless the inputs, outside
#                  the two climate files, fingerprint exactly as the headline runs.
#   maskdec   02i  output/exposure_mask_decomposition.csv, interaction included.
#
# They come AFTER 'ge' (maskruns reads the dose-response tables, maskdec the
# headline no-warming, daily and annual runs) and BEFORE '02d' (the emitter
# reads the decomposition).
# ---------------------------------------------------------------------------
etape masks "$PYTHON" 02_GE/02g_exposure_masks.py \
      "intermediate_outputs/exposure_masks/{north,south}/ — masked exposures"
etape maskruns "$PYTHON" 02_GE/02h_exposure_mask_runs.py \
      "4 masked GE runs + exposure_masks/runs_manifest.csv"
etape maskdec "$PYTHON" 02_GE/02i_exposure_mask_decomposition.py \
      "exposure_mask_decomposition.csv, exposure_mask_cells.csv"

# La bande nordique consomme les welfare de la grille ci-dessus et doit donc
# passer APRES 'ge' et AVANT '02d' : emit_values_empirics.py lit son CSV.
etape band "$PYTHON"  02_GE/02e_north_contrast_band.py \
      "ge_north_contrast_band.csv — Northern contrast at the quantiles of the damage functions"

etape 02d "$RSCRIPT" 02_GE/02d_emit_latex.R \
      "tab_eta.tex; values.tex, values_empirics.tex, tab_ge*.tex, tab_binning.tex, Figures 3-4 (the three emitters)"

# ---------------------------------------------------------------------------
# releve_pa1 (lot 57, 28/09/2026) -- le point d'arret 1 devenu RELEVE, NON BLOQUANT.
# Author's decision : la course va au bout ; l'etape ecrit
# logs/releve_pa1_<date>.txt (EMP avant/apres, paires Lo/Hi, tests joints,
# effectifs, neuf controles) et rend la main, QUEL QUE SOIT son resultat.
# Lot 85x (11/10/2026) : deplace d'apres ct7 a apres 02d. Il relance
# emit_values_empirics.py, qui lit output/config_params.csv, et les neuf controles,
# qui lisent output/values.tex : deux sorties que SEUL 02d ecrit. Sur un clone neuf,
# place apres ct7, il rendait « [STOP] config_params.csv absent de output/ ».
# ---------------------------------------------------------------------------
if actif releve_pa1; then
  echo "------------------------------------------------------------"
  echo "[releve_pa1] tools/etape_releve_pa1.py  (NON-BLOCKING)"
  echo "     start   : $(date '+%H:%M:%S')"
  t0r=$(date +%s)
  set +e
  "$PYTHON" "$JMP_ROOT/tools/etape_releve_pa1.py" 2>&1 | tee "$LOGDIR/releve_pa1.log"
  rcr=${PIPESTATUS[0]}
  set -e
  dr=$(( $(date +%s) - t0r ))
  echo "     end     : $(date '+%H:%M:%S')  (code $rcr, non-blocking)"
  FAITES+=("$(printf '%-14s %3d min %02d s  (code %d, non-blocking)' "releve_pa1" $((dr/60)) $((dr%60)) "$rcr")")
  arret_si_until releve_pa1
else
  SAUTEES+=("releve_pa1 (before --from $FROM)"); echo "  ⏭  releve_pa1     skipped (--from $FROM)"
fi

# ---------------------------------------------------------------------------
# LES EMETTEURS NE SONT PAS DES ETAPES : 02d LES LANCE.
#
# 02d_emit_latex.R lance les TROIS emetteurs (emit_values_empirics.py,
# emit_values.py, emit_figures.py) et s'arrete dur si l'un echoue. Deux
# etapes 'emit_values' et 'emit_empirics' etaient declarees ICI, apres 02d :
# ces deux-la tournaient donc DEUX FOIS a chaque run complet (emit_figures.py,
# lui, ne tournait qu'une fois -- l'asymetrie disait deja que le doublon
# n'etait pas voulu).
#
# Et la seconde execution portait une branche FAUSSE : elle lisait le code 1
# comme << au moins une macro en MISMATCH >>, alors que les neuf gardes dures
# de emit_values_empirics.py et toute exception non rattrapee rendent 1 elles
# aussi. Sans 02d en amont, un FileNotFoundError aurait ete classe MISMATCH
# et le run declare reussi. La distinction se fait desormais par un code de
# sortie distinct (3 = verdict), lu par 02d.
#
# Reprise : ces emetteurs se rejouent par  ./run_all.sh --from 02d
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 8. CLOSING INVENTORY
# ---------------------------------------------------------------------------
echo "============================================================"
echo " SUMMARY — $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo
echo "RAN (${#FAITES[@]})"
for l in "${FAITES[@]:-}"; do [ -n "$l" ] && echo "  ✅ $l"; done
echo
if [ "${#SAUTEES[@]}" -gt 0 ]; then
  echo "SKIPPED BY --from (${#SAUTEES[@]})"
  for l in "${SAUTEES[@]}"; do echo "  ⏭  $l"; done
  echo
fi
cat <<'EOF'
SKIPPED BY CONSTRUCTION
  ⏭  00_era5_distribution.R          SERVER — the ERA5 GRID is supplied and checked at
                                     the start; the climate BY COUNTRY is rebuilt locally
                                     (step climat_local)
  ⏭  00_temperature_projections_JMP  SERVER — outputs supplied, checked at the start
EOF
echo

# ---------------------------------------------------------------------------
# 9. WHAT THE CHAIN DID NOT PRODUCE -- checked against JMP_V11.tex, the live
#    manuscript. It pointed at the frozen JMP_V10.tex until 25/09/2026 and
#    reported eight "missing" targets that were V10's, not V11's.
#
# The list of missing targets is measured on the run that just happened, not
# inferred from reading the manuscript.
# ---------------------------------------------------------------------------
TEX="$JMP_ROOT/JMP_V11.tex"
echo "------------------------------------------------------------"
echo "TARGETS OF THE MANUSCRIPT NOT PRODUCED"
echo "------------------------------------------------------------"
if [ ! -f "$TEX" ]; then
  echo "  ⚠ $TEX missing: cannot compare."
else
  # grep returns 1 when it finds nothing and 2 on a real error. Under `set -e`
  # a `|| true` would swallow both, and an unreadable .tex would pass for a .tex
  # with no targets. Codes 0 and 1 are accepted here; 2 stops the run.
  extraire () {   # extraire <motif> <sed> <destination>
    local rc
    set +e
    grep -n "$1" "$TEX" | grep -v '^[0-9]*:[[:space:]]*%' | sed "$2" > "$3"
    rc=${PIPESTATUS[0]}
    set -e
    if [ "$rc" -gt 1 ]; then
      echo "  ⛔ grep failed (code $rc) on $TEX — comparison abandoned." >&2
      return 1
    fi
    return 0
  }
  # Live \input{...} targets; commented lines are excluded
  extraire '\\input{' 's/^\([0-9]*\):.*\\input{\([^}]*\)}.*/\1 \2/' "$LOGDIR/_inputs.txt"
  # Live \includegraphics{...} targets
  extraire '\\includegraphics' 's/^\([0-9]*\):.*\\includegraphics\(\[[^]]*\]\)\?{\([^}]*\)}.*/\1 \3/' "$LOGDIR/_figs.txt"

  NI=$(wc -l < "$LOGDIR/_inputs.txt"); NF=$(wc -l < "$LOGDIR/_figs.txt")
  echo "  live targets: $NI \\input + $NF \\includegraphics = $((NI+NF))"
  echo

  # producteur <target> -> the step and the script that WRITE that target, printed
  # next to each missing target so that a replicator knows where to resume.
  # Lot 86b (12/10/2026) : refaite sur les ecritures reelles (writeLines, fwrite,
  # pdf(), savefig, open(..., 'w') de chaque script), et non plus sur les noms. Les
  # motifs specifiques passent avant les motifs larges. Elle portait trois erreurs :
  # tab_robustness_prod_joint attribuee a 01a (ct4 l'ecrit), fig_ge_tiers et
  # fig_contrast_grid / tab_contrast_grid absentes, re-emissions de 01a/01c ignorees.
  producteur () {
    case "$1" in
      *values_empirics.tex)          echo "step 02d (emit_values_empirics.py)" ;;
      *values.tex)                   echo "step 02d (emit_values.py)" ;;
      *tab_ge.tex|*tab_ge_*)         echo "step 02d (emit_values.py)" ;;
      *tab_binning*)                 echo "step 02d (emit_values_empirics.py, from the outputs of ct8)" ;;
      *tab_eta*)                     echo "step 02d (02d_emit_latex.R)" ;;
      *fig_ge_tiers*|*ge_welfare_trajectories*|*map_welfare_*) echo "step 02d (emit_figures.py)" ;;
      *tab_crh_northern*)            echo "step 01d_bis (01d_crh_triangulation.R, second pass; needs crh_inference_deliverable.csv)" ;;
      *fig_crh_responses*)           echo "step crhfig (01d_fig_crh_responses.R)" ;;
      *fig_drf_crh_*)                echo "step 01d (01d_crh_triangulation.R)" ;;
      *MC_prod_daily_fan_paper*|*MC_prod_annual_fan_paper*) echo "step ct6 (ct6_fig_drf_paper.R)" ;;
      *MC_amenity_daily_fan_paper*|*MC_amenity_annual_fan_paper*) echo "step 01b (01b_amenity_drf.R)" ;;
      *tab_robustness_prod_joint*|*tab_contrast_grid*|*fig_contrast_grid*) echo "step ct4 (ct4_emit_tab_contribA_grid.R)" ;;
      *tab_fe_absorption*|*tab_ctrl_absorption*) echo "step 01a_reemit (01a_reemit_tab_absorption.R)" ;;
      *tab_oh_bounds*)               echo "step 01c_reemit (01c_reemit_tab_oh_bounds.R)" ;;
      *tab_oster*|*tab_oh_*|*fig_oh_asym_fe*) echo "step 01c (01c_robustness.R)" ;;
      *fig_oh_*)                     echo "none: no script of the package writes this fig_oh_*" ;;
      *tab_stage1_migration*|*exogeneity_mrt*|*robustness_amenity*|*MC_amenity_*) echo "step 01b (01b_amenity_drf.R)" ;;
      *tab_stage1_trade*|*stability_controls*|*appendix_identification*|*MC_prod_*) echo "step 01a (01a_productivity_drf.R)" ;;
      *tab_dbar_cross*|*dbar_cross_weighting*|*map_nonstationarity*|*map_rho_hist*) echo "step 01e (01e_diagnostics.R)" ;;
      *myplot_migration*|*no_activity_map*) echo "step 02a (02a_calibration_at_1deg.r)" ;;
      *fig_M1_net_by_latitude*|*migration_net_by_latitude*) echo "step miglat (02f_fig_mig_latitude.py)" ;;
      *) echo "unknown producer" ;;
    esac
  }

  NMISS=0
  while read -r ln p; do
    [ -z "${p:-}" ] && continue
    if [ ! -f "$JMP_ROOT/$p" ]; then
      printf "  ⛔ %-46s (l.%-5s) <- %s\n" "$p" "$ln" "$(producteur "$p")"
      NMISS=$((NMISS+1))
    fi
  done < "$LOGDIR/_inputs.txt"

  while read -r ln p; do
    [ -z "${p:-}" ] && continue
    found=""
    for c in "$p" "output/$p" "$p.pdf" "output/$p.pdf" "$p.png" "output/$p.png"; do
      [ -f "$JMP_ROOT/$c" ] && { found=1; break; }
    done
    if [ -z "$found" ]; then
      printf "  ⛔ %-46s (l.%-5s) <- %s\n" "$p" "$ln" "$(producteur "$p")"
      NMISS=$((NMISS+1))
    fi
  done < "$LOGDIR/_figs.txt"

  echo
  if [ "$NMISS" -eq 0 ]; then
    echo "  ✅ every live target is present."
  else
    echo "  $NMISS target(s) missing out of $((NI+NF))."
    echo "  ⇒ The manuscript will not compile while they are missing."
  fi
fi

echo
echo "Logs of this run: $LOGDIR"
echo "============================================================"
