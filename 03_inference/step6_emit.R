# step6_emit.R -- reshape the inference results into the deliverable table.
#
# Emits one row per (channel, vcov variant, unit, object), carrying the point
# estimate, the 95% interval, and a flag saying whether the row is publishable.
#
# Reads:  crh_inference_raw.csv, crh_bootstrap_raw.csv (optional)
# Writes: crh_inference_deliverable.csv
# Needs:  step4_mc.R, and step5_boot.R if the bootstrap rows are wanted.
#
# The publication rule is pre-committed: report the more conservative of the two
# variance estimators. The script therefore checks which variant is uniformly
# wider rather than assuming it, and marks every row accordingly. It emits no
# LaTeX macro; the macro layer reads this table.

suppressMessages({ library(data.table) })
## Paths come from 00_parameters.R, which derives the root from its own location.
## THREE SOURCES, NOT ONE. The variable this script used to call `base` stood for
## three different things, and they must not be repointed together:
##   DIR_CRH    outputs of 01d_crh_triangulation.R -- the inputs of this chain;
##   eggw_csv() public inputs read inside the CRH archive, never extracted;
##   DIR_REF    frozen reproduction witnesses, see below.
.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("[STOP] set JMP_ROOT, or launch with Rscript from the repository root.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_eggw.R"))

## DIR_REF holds the FROZEN reproduction witnesses. They are NOT produced by any
## step of run_all.sh and must never be: the day 01d overwrites them, the checks
## below compare the chain against itself, report a difference of zero, and three
## guards go green without having verified anything.
DIR_REF  <- file.path(.root, "01_estimation", "reference_repro_crh")
DIR_WORK <- file.path(DIR_TMP, "crh_inference")
DIR_INF  <- file.path(DIR_OUT, "crh_inference")
dir.create(DIR_WORK, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_INF,  showWarnings = FALSE, recursive = TRUE)
r <- fread(file.path(DIR_INF, "crh_inference_raw.csv"))

mk <- function(obj, pt, lo, hi, se, flag) data.table(
  channel=r$channel, vcov_variant=r$vcov, unit=r$unit, n_cells=r$n_cells, object=obj,
  point=r[[pt]], ci95_lo=r[[lo]], ci95_hi=r[[hi]],
  se=if (is.null(se)) NA_real_ else r[[se]], flag=flag)

out <- rbindlist(list(
  mk("level_mean_2095",   "mean_point","mean_lo","mean_hi","mean_se", "-"),
  mk("level_daily_2095",  "daily_point","daily_lo","daily_hi","daily_se", "-"),
  mk("difference_mean_minus_daily", "diff_point","diff_lo","diff_hi","diff_se", "-"),
  mk("attenuation_ratio", "ratio_point","ratio_lo","ratio_hi", NULL, r$ratio_flag)))
# Share-type outputs: a point only, no interval. A share of draws is not an
# estimate with a sampling distribution of its own.
sh <- rbindlist(list(
  data.table(channel=r$channel, vcov_variant=r$vcov, unit=r$unit, n_cells=r$n_cells,
             object="share_draws_daily_sign_flips", point=r$p_daily_signflip,
             ci95_lo=NA_real_, ci95_hi=NA_real_, se=NA_real_, flag="-"),
  data.table(channel=r$channel, vcov_variant=r$vcov, unit=r$unit, n_cells=r$n_cells,
             object="share_draws_difference_sign_flips", point=r$p_diff_signflip,
             ci95_lo=NA_real_, ci95_hi=NA_real_, se=NA_real_, flag="-"),
  data.table(channel=r$channel, vcov_variant=r$vcov, unit=r$unit, n_cells=r$n_cells,
             object="prob_denominator_le_zero", point=r$p_denom_le0,
             ci95_lo=NA_real_, ci95_hi=NA_real_, se=NA_real_, flag="-")))
out <- rbindlist(list(out, sh))
setorder(out, channel, unit, object, vcov_variant)

# ---- which variant is uniformly wider? ---------------------------------------
# The pre-committed rule picks the more conservative estimator, so any interval
# where Conley is wider than country clustering is an exception that must be
# named rather than averaged away.
w <- dcast(out[!is.na(ci95_lo)], channel+unit+object ~ vcov_variant,
           value.var = c("ci95_lo","ci95_hi"))
w[, width_conley := ci95_hi_conley - ci95_lo_conley]
w[, width_country := ci95_hi_country - ci95_lo_country]
cat("Intervals where CONLEY is wider than COUNTRY (exceptions to 'country is conservative'):\n")
ex <- w[width_conley > width_country]
if (nrow(ex)) print(ex[, .(channel, unit, object, width_conley, width_country)]) else
  cat("  NONE -- country clustering is uniformly the more conservative variant.\n")
out[, publish := fifelse(vcov_variant == "country", "PUBLISH (more conservative)", "reference (original CRH structure)")]

# ---- bootstrap control route, flagged as NOT for publication -----------------
# The bootstrap is a robustness control on the analytic intervals, not a second
# set of published numbers. Its rows carry that flag so that a later reader
# cannot mistake them for the reported interval.
bf <- file.path(DIR_INF, "crh_bootstrap_raw.csv")
if (file.exists(bf)) {
  b <- fread(bf)
  mkb <- function(obj, pt, lo, hi) data.table(
    channel=b$channel, vcov_variant=b$vcov, unit=b$unit, n_cells=NA_integer_, object=obj,
    point=b[[pt]], ci95_lo=b[[lo]], ci95_hi=b[[hi]], se=NA_real_,
    flag=paste0("control s4.3, B=", b$n_boot),
    publish="CONTROL ONLY (never published, per s4.4 pre-commitment)")
  ob <- rbindlist(list(
    mkb("level_mean_2095","mean_point","mean_lo","mean_hi"),
    mkb("level_daily_2095","daily_point","daily_lo","daily_hi"),
    mkb("difference_mean_minus_daily","diff_point","diff_lo","diff_hi"),
    mkb("attenuation_ratio","ratio_point","ratio_lo","ratio_hi")))
  out <- rbindlist(list(out, ob), use.names=TRUE)
  cat("\nAppended", nrow(ob), "control rows for channels:", paste(unique(b$channel), collapse=","), "\n")
}
setorder(out, channel, unit, object, vcov_variant)

fwrite(out, file.path(DIR_INF, "crh_inference_deliverable.csv"))
cat("\nWritten crh_inference_deliverable.csv :", nrow(out), "rows\n")
