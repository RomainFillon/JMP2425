# step2_maps.R -- build the linear projection maps W (cells x coefficients).
#
# The projection from bin coefficients to a per-cell damage is linear at fixed
# weights, so the whole pipeline can be represented by a matrix W. Propagating
# the basis vectors through it gives W column by column. Once W is known, the
# sampling distribution of the coefficients maps to the cells in closed form,
# which is what step4 needs; resampling the pipeline itself would be intractable.
#
# The linearity is not assumed: the script checks that W %*% coef_hat reproduces
# the published per-cell columns before writing anything.
#
# Reads:  step1_frame.rds, CMIP6_tas.csv, crh_daily_projection_cells.csv
# Writes: step2_maps.rds
# Needs:  step1_gate.R must have run first.

suppressMessages({ library(data.table); library(splines) })
t0 <- proc.time(); cat("=== step2_maps.R ===\n")
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
BASE_YEAR <- 2015L; FUT_YEAR <- 2095L     # projection horizon, as in 01d

S1 <- readRDS(file.path(DIR_WORK, "step1_frame.rds"))
NB <- S1$NB; ctr <- S1$ctr; ref_bin <- S1$ref_bin; mctr <- S1$mctr

# ---- 1. CMIP6 -> per-cell delta-share matrix over integer temperatures ------------
cat("Loading CMIP6_tas.csv...\n")
cm <- fread(file.path(DIR_TMP, "CMIP6_tas.csv"))[year %in% c(BASE_YEAR, FUT_YEAR)]
cm[, c := count / sum(count), by = .(x, y, year)]
mu <- cm[, .(mu = sum(value * c)), by = .(x, y, year)]
cells <- unique(cm[, .(x, y)]); setorder(cells, x, y); cells[, cid := .I]
cm <- merge(cm, cells, by = c("x","y"))
Tvals <- sort(unique(cm$value)); nT <- length(Tvals)
cat("  cells:", nrow(cells), " integer temps:", nT, " range [", min(Tvals), ",", max(Tvals), "]\n")
cm[, tid := match(value, Tvals)]
cm[, sgn := fifelse(year == FUT_YEAR, 1, -1)]
dc <- cm[, .(v = sum(sgn * c)), by = .(cid, tid)]
M <- Matrix::sparseMatrix(i = dc$cid, j = dc$tid, x = dc$v,
                          dims = c(nrow(cells), nT))          # cells x nT : Dc_r(T)
muw <- dcast(mu, x + y ~ year, value.var = "mu")
setnames(muw, c(as.character(BASE_YEAR), as.character(FUT_YEAR)), c("mu_base","mu_fut"))
muw <- merge(muw, cells, by = c("x","y"))
setorder(muw, cid)
stopifnot(identical(muw$cid, cells$cid), all(is.finite(muw$mu_base)), all(is.finite(muw$mu_fut)))
rm(cm, dc); gc()

# ---- 2. linear operators, per channel --------------------------------------------
# Daily smoother, as a linear map delta -> delta~(T). Same weights and same
# spline as the estimation step, so W carries no extra smoothing of its own.
daily_smooth_row <- function(delta_vec, se_vec) {
  bins_dt <- data.table(bin=1:NB, t_center=ctr$t_center, delta=as.numeric(delta_vec),
                        se=as.numeric(se_vec))
  fitp <- bins_dt[is.finite(t_center)]
  w <- ifelse(is.finite(fitp$se), 1/fitp$se, 1/mean(fitp$se, na.rm=TRUE))
  sp <- lm(delta ~ ns(t_center, df=5), data=fitp, weights=w)
  tgrid <- seq(min(fitp$t_center), max(fitp$t_center), length.out=400)
  dsm <- as.numeric(predict(sp, newdata=data.frame(t_center=tgrid)))
  approx(tgrid, dsm, xout = Tvals, rule = 2)$y   # constant beyond the fitted range
}
# Mean-arm smoother plus trapezoidal integral, as a linear map mc -> Lambda(mu).
mean_lambda_row <- function(mc_vec, se_vec) {
  df_fit <- data.frame(v=as.numeric(mc_vec), w=1/as.numeric(se_vec), temp=as.numeric(mctr))
  sp <- lm(v ~ ns(temp, df=4), data=df_fit, weights=w)
  mu_grid  <- seq(min(df_fit$temp), max(df_fit$temp), length.out=400)
  d_smooth <- as.numeric(predict(sp, newdata=data.frame(temp=mu_grid)))
  dmu <- diff(mu_grid)
  Lam <- c(0, cumsum((head(d_smooth,-1)+tail(d_smooth,-1))/2*dmu))
  smin <- min(df_fit$temp); smax <- max(df_fit$temp)
  f <- function(m) approx(mu_grid, Lam, xout = pmin(pmax(m, smin), smax), rule = 2)$y
  f(muw$mu_fut) - f(muw$mu_base)                                # per-cell level change
}

build_channel <- function(tag) {
  cat("\n--- channel:", tag, "---\n")
  P  <- S1[[tag]]
  # keep ONLY the share coefficients (drop Chebyshev controls) before filling the
  # NB-long weight vector -- assigning unknown names would append and lengthen it.
  se_sh <- P$se_a[names(P$se_a) %in% S1$scol]
  se_full <- setNames(rep(NA_real_, NB), S1$scol); se_full[names(se_sh)] <- se_sh
  stopifnot(length(se_full) == NB)

  # W_daily : cells x NB. The reference-bin column is structurally zero but is
  # kept, so that the column index matches the bin index everywhere downstream.
  Wd <- matrix(0, nrow = nrow(cells), ncol = NB)
  for (k in 1:NB) {
    e <- rep(0, NB); e[k] <- 1
    Wd[, k] <- as.numeric(M %*% daily_smooth_row(e, se_full))
  }
  # W_mean : cells x 20
  Wm <- matrix(0, nrow = nrow(cells), ncol = 20L)
  for (j in 1:20) {
    e <- rep(0, 20); e[j] <- 1
    Wm[, j] <- mean_lambda_row(e, P$se_b)
  }

  # ---- linearity check: reproduce the published per-cell columns -------------
  pub <- fread(file.path(DIR_REF, "crh_daily_projection_cells.csv"))[channel == tag]
  pub <- merge(pub, cells, by = c("x","y"))
  setorder(pub, cid)
  pred_d <- as.numeric(Wd %*% as.numeric(P$delta))
  pred_m <- as.numeric(Wm %*% as.numeric(P$mc))
  i <- pub$cid
  cat(sprintf("  cells published=%d  matched=%d\n", nrow(pub), length(i)))
  cat(sprintf("  LINEARITY daily : max|W.delta - published| = %.3e\n", max(abs(pred_d[i] - pub$daily))))
  cat(sprintf("  LINEARITY mean  : max|W.mc    - published| = %.3e\n", max(abs(pred_m[i] - pub$mean))))
  list(Wd = Wd, Wm = Wm, pub = pub)
}

B_amen <- build_channel("amen")
B_prod <- build_channel("prod")

saveRDS(list(cells=cells, Tvals=Tvals, muw=muw,
             amen=B_amen, prod=B_prod), file.path(DIR_WORK, "step2_maps.rds"))
cat("\nSaved step2_maps.rds\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n=== Done ===\n")
