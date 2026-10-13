# step4_mc.R -- inference on the aggregated quantities.
#
# Draws from the joint sampling distribution of (delta_daily, delta_mean) and
# pushes each draw through the projection map W and the population-weighted
# aggregation. Because both steps are linear, the draw can be applied to the
# coefficients rather than re-running the pipeline.
#
# Produces, per channel and per unit: the two levels at 2095, their difference,
# the attenuation ratio, and the share of draws in which a sign flips.
#
# Reads:  step1_frame.rds, step2_maps.rds, step3_vcov.rds,
#         tab_crh_northern_reversal_daily.csv (reproduction reference)
# Writes: crh_inference_raw.csv
# Needs:  steps 1 to 3 must have run first.
#
# The attenuation ratio is a ratio of two estimated quantities, so its
# distribution is not symmetric and its denominator can approach zero. The share
# of draws with a non-positive denominator is reported alongside it, because a
# percentile interval on such a ratio is only readable once that share is known.

suppressMessages({ library(data.table) })
t0 <- proc.time(); cat("=== step4_mc.R ===\n")
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
NDRAW <- 20000L
set.seed(20260728L)          # fixed in advance, so the draw cannot be reselected
PROBS <- c(0.025, 0.5, 0.975)

S1 <- readRDS(file.path(DIR_WORK, "step1_frame.rds"))
S2 <- readRDS(file.path(DIR_WORK, "step2_maps.rds"))
S3 <- readRDS(file.path(DIR_WORK, "step3_vcov.rds"))
UNITS <- c("RUS","CAN","FIN","DEU","FRA","NORTH>=45N")

## The root of V. Lot 85l (06/10/2026), the author's decision: negative
## eigenvalues are truncated to zero, V+ = U diag(max(lambda, 0)) U' -- a
## reconstruction that does not depend on the sign of the eigenvectors -- and the
## root is the pivoted Cholesky factor of V+ (R'R = V+). No jitter. It replaces
## the spectral root Z diag(sqrt(lambda+)) U', which jumps when an eigenvector
## flips sign under a 1e-13 perturbation (lot 85k). Same normals, same order.
## One V is not PSD (conley, two eigenvalues near -2.5e-05): V+ then has rank
## k - 2, and chol(pivot = TRUE) warns that it is rank-deficient -- expected, and
## the reconstruction check below is what decides. tools/test_step4_cholesky.R.
psd_chol_root <- function(V) {
  e  <- eigen((V + t(V))/2, symmetric = TRUE)
  Vp <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  Vp <- (Vp + t(Vp))/2
  R  <- withCallingHandlers(chol(Vp, pivot = TRUE), warning = function(w)
          if (grepl("rank-deficient", conditionMessage(w))) invokeRestart("muffleWarning"))
  k  <- ncol(Vp); r <- attr(R, "rank")
  if (r < k) R[(r + 1):k, (r + 1):k] <- 0
  R  <- R[, order(attr(R, "pivot")), drop = FALSE]
  err <- max(abs(crossprod(R) - Vp)) / max(abs(Vp))
  if (!is.finite(err) || err > 1e-10)
    stop(sprintf("[STOP] step4 : R'R does not reconstruct V+ (relative error %.3e).", err))
  list(R = R, rank = r, min_eig = min(e$values), n_neg = sum(e$values < 0),
       tot_neg = sum(pmin(e$values, 0)),
       min_eig_after = min(eigen(Vp, symmetric = TRUE, only.values = TRUE)$values))
}

draw_mvn <- function(V, n) {
  rt <- psd_chol_root(V)
  list(Z = matrix(rnorm(n * ncol(V)), n, ncol(V)) %*% rt$R,
       n_neg = rt$n_neg, min_eig = rt$min_eig, tot_neg = rt$tot_neg,
       min_eig_after = rt$min_eig_after, rank = rt$rank)
}

out <- list()
for (tag in c("prod","amen")) {
  cat("\n=========== channel:", tag, "===========\n")
  P  <- S1[[tag]]; B <- S2[[tag]]; R <- S3[[tag]]
  Wd <- B$Wd; Wm <- B$Wm; pub <- B$pub          # pub is ordered by cid, has pop + iso3
  stopifnot(identical(pub$cid, S2$cells$cid))

  # Point coefficient vectors, as refitted and checked in step1.
  d_hat_full <- as.numeric(P$delta)                    # length NB, ref bin = 0
  m_hat      <- as.numeric(P$mc)                       # length 20
  # Index of the estimated bins inside the NB-long vector. The reference bin is
  # excluded: it is normalised to zero and has no sampling variation.
  ia <- match(R$names_a, S1$scol)
  stopifnot(!any(is.na(ia)), length(ia) == length(R$names_a))
  cat(sprintf("  coef check daily: max|step1 - stacked| = %.3e\n",
              max(abs(d_hat_full[ia] - as.numeric(R$coef_a)))))
  cat(sprintf("  coef check mean : max|step1 - stacked| = %.3e\n",
              max(abs(m_hat        - as.numeric(R$coef_b)))))

  # ---- aggregation weights per unit --------------------------------------------
  # Population weights, fixed at the base year. Letting them move with the
  # scenario would mix a composition change into the damage estimate.
  A <- sapply(UNITS, function(u) {
    sel <- if (u == "NORTH>=45N") pub$y >= 45 else pub$iso3 == u
    sel[is.na(sel)] <- FALSE
    w <- ifelse(sel, pub$pop, 0); w[is.na(w)] <- 0
    w / sum(w)
  })                                                    # cells x units
  ncell <- sapply(UNITS, function(u) {
    sel <- if (u == "NORTH>=45N") pub$y >= 45 else pub$iso3 == u
    sum(sel, na.rm = TRUE) })

  # Each aggregate is a linear functional of the coefficients, so its variance
  # follows in closed form from V and no resampling of the pipeline is needed.
  Vd <- t(Wd) %*% A        # NB x units   (daily)
  Vm <- t(Wm) %*% A        # 20 x units   (mean)
  lvl_d <- as.numeric(t(Vd) %*% d_hat_full)
  lvl_m <- as.numeric(t(Vm) %*% m_hat)

  # ---- reproduction check against the frozen reference aggregates ------------
  # A FROZEN WITNESS, not a published table. It comes from a reference run of our
  # own reproduction of the CRH design, and it must never be repointed at DIR_CRH:
  # the chain would then be compared with itself, the difference would be zero by
  # construction, and this check would go green having verified nothing.
  #
  # A NON-ZERO DIFFERENCE IS EXPECTED HERE, and it is the object of the check. The
  # inputs of this chain are regenerated by 01d, whose productivity arm departs
  # from the reference run; the amenity arm reproduces it to about 1e-05. The size
  # of the gap is what is reported, and it does not stop the chain.
  tabpub <- fread(file.path(DIR_REF, "tab_crh_northern_reversal_daily.csv"))
  tp <- tabpub[channel == tag][match(UNITS, iso3)]
  cat(sprintf("  REPRO vs reference, daily: max|diff| = %.3e\n", max(abs(lvl_d - tp$daily))))
  cat(sprintf("  REPRO vs reference, mean : max|diff| = %.3e\n", max(abs(lvl_m - tp$mean))))
  cat(sprintf("  REPRO n_cells: %s\n", paste(ncell == tp$n_cells, collapse=" ")))

  for (vtag in c("conley","country")) {
    V <- if (vtag == "conley") R$V_conley else R$V_country
    dr <- draw_mvn(V, NDRAW)
    cat(sprintf("\n  [%s] vcov dim=%d  min eigenvalue=%.3e  #negative=%d (clipped)\n",
                vtag, ncol(V), dr$min_eig, dr$n_neg))
    cat(sprintf("  [%s] after truncation: min eigenvalue=%.3e  rank=%d  (pivoted Cholesky)\n",
                vtag, dr$min_eig_after, dr$rank))
    pa <- length(R$names_a)
    Da <- dr$Z[, 1:pa, drop=FALSE]                       # draws of (delta - delta_hat)
    Db <- dr$Z[, (pa+1):ncol(V), drop=FALSE]             # draws of (mc - mc_hat)

    # Levels per unit, one row per draw
    Ld <- Da %*% Vd[ia, , drop=FALSE]; Ld <- sweep(Ld, 2, lvl_d, "+")
    Lm <- Db %*% Vm;                   Lm <- sweep(Lm, 2, lvl_m, "+")
    Diff <- Lm - Ld
    Ratio <- 1 - Ld / Lm

    # Analytic standard error of the levels. The map is linear, so this is exact
    # given V and serves as a check on the draw-based interval.
    se_d <- sqrt(pmax(0, diag(t(Vd[ia,,drop=FALSE]) %*% V[1:pa,1:pa] %*% Vd[ia,,drop=FALSE])))
    se_m <- sqrt(pmax(0, diag(t(Vm) %*% V[(pa+1):ncol(V),(pa+1):ncol(V)] %*% Vm)))
    Gd <- rbind(Vd[ia,,drop=FALSE]*-1, Vm)               # functional of (mean - daily)
    se_diff <- sqrt(pmax(0, diag(t(Gd) %*% V %*% Gd)))
    # What the interval would be if the cross-arm covariance were ignored. Kept
    # to show how much the joint treatment of the two arms actually matters.
    se_diff_indep <- sqrt(se_m^2 + se_d^2)
    cov_arms <- diag(t(Vd[ia,,drop=FALSE]) %*% V[1:pa,(pa+1):ncol(V)] %*% Vm)
    cor_arms <- cov_arms / (se_d * se_m)

    q <- function(X) t(apply(X, 2, quantile, probs = PROBS, na.rm = TRUE))
    qd <- q(Ld); qm <- q(Lm); qf <- q(Diff); qr <- q(Ratio)

    for (k in seq_along(UNITS)) {
      denom_ci_zero <- (qm[k,1] <= 0 & qm[k,3] >= 0)
      p_den_le0 <- mean(Lm[,k] <= 0)
      flag <- if (denom_ci_zero) "NON-REPORTABLE (denominator CI covers 0)" else
              if (p_den_le0 > 0.01) "FRAGILE (P(mean<=0)>1%)" else "ok"
      out[[length(out)+1]] <- data.table(
        channel = tag, vcov = vtag, unit = UNITS[k], n_cells = ncell[k],
        mean_point = lvl_m[k], mean_lo = qm[k,1], mean_hi = qm[k,3], mean_se = se_m[k],
        daily_point = lvl_d[k], daily_lo = qd[k,1], daily_hi = qd[k,3], daily_se = se_d[k],
        diff_point = lvl_m[k]-lvl_d[k], diff_lo = qf[k,1], diff_hi = qf[k,3], diff_se = se_diff[k],
        ratio_point = 1 - lvl_d[k]/lvl_m[k], ratio_lo = qr[k,1], ratio_hi = qr[k,3],
        diff_se_indep = se_diff_indep[k], se_inflation_if_indep = se_diff_indep[k]/se_diff[k],
        cor_arms = cor_arms[k],
        ratio_flag = flag, p_denom_le0 = p_den_le0,
        p_daily_signflip = mean(sign(Ld[,k]) != sign(lvl_d[k])),
        p_diff_signflip  = mean(sign(Diff[,k]) != sign(lvl_m[k]-lvl_d[k])),
        n_neg_eig = dr$n_neg, min_eig = dr$min_eig)
    }
  }
}
res <- rbindlist(out)
fwrite(res, file.path(DIR_INF, "crh_inference_raw.csv"))

cat("\n\n================ RESULTS (95% intervals, ", NDRAW, " draws) ================\n", sep="")
for (tg in c("prod","amen")) for (vt in c("conley","country")) {
  cat(sprintf("\n--- channel=%s  vcov=%s ---\n", tg, vt))
  r <- res[channel==tg & vcov==vt]
  print(r[, .(unit,
              mean = sprintf("%+.3f [%+.3f,%+.3f]", mean_point, mean_lo, mean_hi),
              daily = sprintf("%+.3f [%+.3f,%+.3f]", daily_point, daily_lo, daily_hi),
              diff = sprintf("%+.3f [%+.3f,%+.3f]", diff_point, diff_lo, diff_hi),
              atten = sprintf("%.1f%% [%.1f,%.1f]", 100*ratio_point, 100*ratio_lo, 100*ratio_hi),
              flag = ratio_flag, flipD = round(p_daily_signflip,3), flipDiff = round(p_diff_signflip,3),
              rho = round(cor_arms,3), infl = round(se_inflation_if_indep,2))],
        nrow = 20)
}
cat("\nSaved crh_inference_raw.csv\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n=== Done ===\n")
