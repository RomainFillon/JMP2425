# step5_boot.R -- control route: pairs cluster bootstrap.
#
# Resamples countries with replacement and RE-ESTIMATES both arms at each draw,
# then reports percentile intervals. Unlike step4, which propagates an analytic
# covariance through a linear map, this route makes no distributional assumption
# about the coefficients.
#
# Reads:  step1_frame.rds, step2_maps.rds, step3_vcov.rds
# Writes: crh_bootstrap_raw.csv
# Needs:  steps 1 to 3 must have run first.
#
# This is a control, not a second set of published numbers: the reported figure
# comes from step4 whatever this shows. The rule was fixed in advance so that the
# narrower of the two intervals cannot be chosen after seeing both.
#
# It targets the country-clustered variant only. A Conley analogue would need a
# resampling scheme for overlapping spatial neighbourhoods, which has no standard
# definition.
#
# The cost of five iterations is measured first, and the run proceeds only if the
# extrapolated total is within budget. A bootstrap that re-estimates two fixed
# effect models per draw is expensive enough that this has to be checked rather
# than assumed.

suppressMessages({ library(data.table); library(fixest) })
t0 <- proc.time(); cat("=== step5_boot.R ===\n")
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
BUDGET_MIN <- 75      # fixed ex ante: run only if the extrapolated cost fits this
B_TARGET   <- 1000L
set.seed(20260728L)

S1 <- readRDS(file.path(DIR_WORK, "step1_frame.rds")); S2 <- readRDS(file.path(DIR_WORK, "step2_maps.rds"))
S3 <- readRDS(file.path(DIR_WORK, "step3_vcov.rds"))
dat <- S1$dat; use_s <- S1$use_s; che_vars <- S1$che_vars; mb <- paste0("mb_",1:20)
UNITS <- c("RUS","CAN","FIN","DEU","FRA","NORTH>=45N")

suppressMessages({ library(sf); library(rnaturalearth) }); sf_use_s2(FALSE)
cl <- unique(dat[, .(lon,lat)]); pts <- st_as_sf(cl, coords=c("lon","lat"), crs=4326)
world <- ne_countries(scale="medium", returnclass="sf")
ic <- c("iso_a3_eh","adm0_a3","iso_a3","ISO_A3"); ic <- ic[ic %in% names(world)][1]
world <- world[, ic]; names(world)[1] <- "iso3"
cl[, iso3 := suppressMessages(st_join(pts, world, join=st_within))$iso3]
if (any(is.na(cl$iso3)))
  cl[is.na(iso3), iso3 := world$iso3[suppressMessages(st_nearest_feature(pts[is.na(cl$iso3),], world))]]
dat <- merge(dat, cl, by=c("lon","lat"), all.x=TRUE)

run_channel <- function(lhs, extra_rhs, fe_block, subnal_col, tag) {
  cat("\n--- channel:", tag, "---\n")
  d <- copy(dat)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year-1985)/5 * max(subnal) + subnal]
  keep <- c(lhs, use_s, mb, extra_rhs, fe_block, "fe", "iso3")
  d <- d[, ..keep][complete.cases(d[, ..keep])]
  fml_a <- as.formula(paste(lhs,"~",paste(c(use_s, extra_rhs),collapse=" + "),"|",fe_block,"+ fe"))
  fml_b <- as.formula(paste(lhs,"~",paste(c(mb,    extra_rhs),collapse=" + "),"|",fe_block,"+ fe"))
  ctry <- unique(d$iso3); G <- length(ctry)
  setkey(d, iso3)

  P <- S1[[tag]]; B <- S2[[tag]]; pub <- B$pub
  ia <- match(use_s, S1$scol)
  A <- sapply(UNITS, function(u) {
    sel <- if (u=="NORTH>=45N") pub$y>=45 else pub$iso3==u; sel[is.na(sel)] <- FALSE
    w <- ifelse(sel, pub$pop, 0); w[is.na(w)] <- 0; w/sum(w) })
  Vd <- t(B$Wd) %*% A; Vm <- t(B$Wm) %*% A

  one <- function() {
    g  <- sample(ctry, G, replace = TRUE)
    db <- d[.(g), allow.cartesian = TRUE]
    fa <- tryCatch(feols(fml_a, data=db, vcov="iid", notes=FALSE), error=function(e) NULL)
    fb <- tryCatch(feols(fml_b, data=db, vcov="iid", notes=FALSE), error=function(e) NULL)
    if (is.null(fa) || is.null(fb)) return(NULL)
    ca <- coef(fa)[use_s]; cb <- coef(fb)[mb]
    if (anyNA(ca) || anyNA(cb)) return(NULL)
    dfull <- rep(0, S1$NB); dfull[ia] <- ca
    c(as.numeric(t(Vd) %*% dfull), as.numeric(t(Vm) %*% cb))
  }

  tt <- system.time(for (i in 1:5) one())["elapsed"]
  per <- tt/5; est <- per * B_TARGET / 60
  cat(sprintf("  cost: %.2f s/iteration -> B=%d would take %.1f min (budget %d min)\n",
              per, B_TARGET, est, BUDGET_MIN))
  if (est > BUDGET_MIN) { cat("  -> PROHIBITIVE, not executed (control is not blocking).\n"); return(NULL) }

  R <- matrix(NA_real_, B_TARGET, 2*length(UNITS)); nfail <- 0L
  for (b in 1:B_TARGET) {
    v <- one(); if (is.null(v)) { nfail <- nfail + 1L; next }; R[b,] <- v
    if (b %% 100 == 0) cat("   .", b, "\n")
  }
  R <- R[complete.cases(R), , drop=FALSE]
  cat(sprintf("  usable draws: %d/%d (failed %d)\n", nrow(R), B_TARGET, nfail))
  nu <- length(UNITS)
  Ld <- R[, 1:nu, drop=FALSE]; Lm <- R[, (nu+1):(2*nu), drop=FALSE]
  lvl_d <- as.numeric(t(Vd) %*% as.numeric(P$delta)); lvl_m <- as.numeric(t(Vm) %*% as.numeric(P$mc))
  q <- function(X) t(apply(X, 2, quantile, probs=c(0.025,0.975), na.rm=TRUE))
  qd <- q(Ld); qm <- q(Lm); qf <- q(Lm-Ld); qr <- q(1-Ld/Lm)
  data.table(channel=tag, vcov="boot_country", unit=UNITS, n_boot=nrow(R),
             mean_point=lvl_m, mean_lo=qm[,1], mean_hi=qm[,2],
             daily_point=lvl_d, daily_lo=qd[,1], daily_hi=qd[,2],
             diff_point=lvl_m-lvl_d, diff_lo=qf[,1], diff_hi=qf[,2],
             ratio_point=1-lvl_d/lvl_m, ratio_lo=qr[,1], ratio_hi=qr[,2],
             p_daily_signflip=colMeans(sign(Ld)!=rep(sign(lvl_d),each=nrow(Ld))),
             p_diff_signflip =colMeans(sign(Lm-Ld)!=rep(sign(lvl_m-lvl_d),each=nrow(Ld))))
}

r1 <- run_channel("log_prod_phi", che_vars,     "land_type", "subnal_reg",  "prod")
r2 <- run_channel("log_amen_HDI", character(0), "grid_2c",   "subnal_ctry", "amen")
res <- rbindlist(list(r1, r2))
if (nrow(res)) {
  fwrite(res, file.path(DIR_INF, "crh_bootstrap_raw.csv"))
  print(res[, .(channel, unit,
                mean=sprintf("%+.3f [%+.3f,%+.3f]", mean_point, mean_lo, mean_hi),
                daily=sprintf("%+.3f [%+.3f,%+.3f]", daily_point, daily_lo, daily_hi),
                diff=sprintf("%+.3f [%+.3f,%+.3f]", diff_point, diff_lo, diff_hi),
                atten=sprintf("%.1f%% [%.1f,%.1f]", 100*ratio_point, 100*ratio_lo, 100*ratio_hi),
                flipD=round(p_daily_signflip,3), flipDiff=round(p_diff_signflip,3))], nrow=30)
}
cat("\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n=== Done ===\n")
