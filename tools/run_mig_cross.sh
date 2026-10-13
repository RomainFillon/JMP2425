#!/usr/bin/env bash
# run_mig_cross.sh -- l'etape mig_cross de run_all.sh (lot 85x, 11/10/2026).
#
# Appelle tools/mig_internal_cross.R (T1, lots 70-72), qui ecrit
# intermediate_outputs/mig_internal_cross.csv, la source de \GEmigInternalCrossM et
# \GEmigInternalCrossShare (emit_values.py, bloc E4).
#
# POURQUOI UNE ETAPE. Jusqu'au lot 85w, T1 se lancait a la main et le CSV de l'atelier
# datait du 01/10. Le clone de test du 10/10 ne l'avait pas : emit_values.py a emis
# les deux macros a \PROVISOIRE avec un simple [warn], et le PDF imprimait
# « ?? million ... ??% ».
#
# L'ARGUMENT --masse-step11. T1 recalcule STEP 10-11 de 02b et s'assure contre la
# masse que 02b a JOURNALISEE (assertion A1), puis contre mig_witness_by_country.csv
# (assertion A2, connu-faux Syrie). La masse est lue ici, dans le 02b.log le plus
# recent sous logs/run_*/ : celui de la course en cours quand la chaine passe par
# 02b, celui de la derniere course sinon (--from mig_cross). Un journal ANTERIEUR a
# migration_matrix.rds ne decrit pas la matrice du disque : refus. Un journal
# perime mais posterieur ne passerait pas A1 (masse recalculee != journal) ni A2.
set -u
: "${RSCRIPT:?RSCRIPT not set}" "${JMP_ROOT:?JMP_ROOT not set}"
MM="$JMP_ROOT/intermediate_outputs/migration_matrix.rds"
[ -f "$MM" ] || { echo "[STOP] $MM missing: run 02b first."; exit 1; }
LOG02B=$(ls -t "$JMP_ROOT"/logs/run_*/02b.log 2>/dev/null | head -1)
[ -n "$LOG02B" ] || { echo "[STOP] no logs/run_*/02b.log: run 02b first."; exit 1; }
if [ "$LOG02B" -ot "$MM" ]; then
  echo "[STOP] $LOG02B is older than migration_matrix.rds: it does not describe the matrix on disk."
  exit 1
fi
M=$(grep -a -A3 'STEP 11:' "$LOG02B" | grep -a 'total internal flows' | head -1 \
    | tr -d '\r' | sed -E 's/.*total internal flows:[[:space:]]*([0-9.]+)[[:space:]]*M.*/\1/')
case "$M" in
  ''|*[!0-9.]*) echo "[STOP] STEP 11 mass unreadable in $LOG02B: '$M'"; exit 1 ;;
esac
echo "[mig_cross] STEP 11 mass read in $LOG02B: $M M"
"$RSCRIPT" "$JMP_ROOT/tools/mig_internal_cross.R" --masse-step11 "$M"
