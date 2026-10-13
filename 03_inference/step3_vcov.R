# step3_vcov.R -- joint clustered covariance of (delta_daily, delta_mean).
#
# The published quantity is a difference between the two arms, so its variance
# needs their cross covariance. The two arms are separate regressions on the SAME
# observations, and that cross term appears in neither marginal vcov. Stacking
# them into one block-diagonal system, with arm-specific regressors and
# arm-specific fixed effects, makes OLS reproduce the two separate fits exactly
# while returning the full joint sandwich.
#
# Two variance estimators are computed, and both are kept:
#   (i)  Conley(550 km), the spatial HAC structure of the original design;
#   (ii) clustering by country, which is coarser and matches the object of
#        interest, a population-weighted aggregate across cells within a country.
# The publication rule picks the more conservative of the two; step6 checks which
# one that is rather than assuming it.
#
# Reads:  step1_frame.rds, lengths_country.csv
# Writes: step3_vcov.rds
# Needs:  step1_gate.R must have run first.

suppressMessages({ library(data.table); library(fixest) })
t0 <- proc.time(); cat("=== step3_vcov.R ===\n")
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
S1 <- readRDS(file.path(DIR_WORK, "step1_frame.rds"))
dat <- S1$dat; scol <- S1$scol; use_s <- S1$use_s; che_vars <- S1$che_vars
mb <- paste0("mb_", 1:20)

# ---- 0. country assignment for the estimation cells --------------------------
# Needed only by the country-clustered variant. A cell is attached to the
# country holding the largest share of it, so each cell enters exactly one
# cluster.
cat("Assigning countries to estimation cells...\n")
suppressMessages({ library(sf); library(rnaturalearth) }); sf_use_s2(FALSE)
cl <- unique(dat[, .(lon, lat)])
pts <- st_as_sf(cl, coords = c("lon","lat"), crs = 4326)
world <- ne_countries(scale = "medium", returnclass = "sf")
iso_col <- c("iso_a3_eh","adm0_a3","iso_a3","ISO_A3"); iso_col <- iso_col[iso_col %in% names(world)][1]
world <- world[, iso_col]; names(world)[1] <- "iso3"
j <- suppressMessages(st_join(pts, world, join = st_within))
cl[, iso3 := j$iso3]
nmiss <- sum(is.na(cl$iso3))
if (nmiss) {  # coastal / small-island cells falling outside any polygon -> nearest country
  idx <- suppressMessages(st_nearest_feature(pts[is.na(cl$iso3), ], world))
  cl[is.na(iso3), iso3 := world$iso3[idx]]
}
cat(sprintf("  cells=%d  countries=%d  assigned by nearest-feature=%d\n",
            nrow(cl), uniqueN(cl$iso3), nmiss))
dat <- merge(dat, cl, by = c("lon","lat"), all.x = TRUE)
stopifnot(!any(is.na(dat$iso3)))

# ---- 1. stacked system, per channel ------------------------------------------------
build_stack <- function(lhs, extra_rhs, fe_block, subnal_col, tag) {
  cat("\n--- channel:", tag, "---\n")
  d <- copy(dat)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year-1985)/5 * max(subnal) + subnal]
  keep <- c(lhs, use_s, mb, extra_rhs, fe_block, "fe", "lat", "lon", "iso3")
  d <- d[, ..keep]
  d <- d[complete.cases(d)]                      # identical sample for both arms
  cat("  stacked-base N:", nrow(d), "\n")

  dA <- copy(d)[, arm := "A"]; dB <- copy(d)[, arm := "B"]
  ds <- rbind(dA, dB); rm(dA, dB)
  for (v in use_s)     set(ds, j = paste0("A_", v), value = ds[[v]] * (ds$arm == "A"))
  for (v in mb)        set(ds, j = paste0("B_", v), value = ds[[v]] * (ds$arm == "B"))
  for (v in extra_rhs) {
    set(ds, j = paste0("A_", v), value = ds[[v]] * (ds$arm == "A"))
    set(ds, j = paste0("B_", v), value = ds[[v]] * (ds$arm == "B"))
  }
  rhs <- c(paste0("A_", use_s), paste0("B_", mb))
  if (length(extra_rhs)) rhs <- c(rhs, paste0("A_", extra_rhs), paste0("B_", extra_rhs))
  fml <- as.formula(paste(lhs, "~", paste(rhs, collapse=" + "),
                          "|", paste0(fe_block, "^arm"), "+ fe^arm"))

  fit <- feols(fml, data = ds, vcov = "iid", fixef.tol = 1e-11)     # coefficients only; vcov reset below
  # Separate per-arm fits, kept to check that the stacked system reproduces them
  fa <- feols(as.formula(paste(lhs,"~",paste(c(use_s, extra_rhs),collapse=" + "),"|",fe_block,"+ fe")),
              data = d, vcov = conley(550), fixef.tol = 1e-11)
  fb <- feols(as.formula(paste(lhs,"~",paste(c(mb,    extra_rhs),collapse=" + "),"|",fe_block,"+ fe")),
              data = d, vcov = conley(550), fixef.tol = 1e-11)
  cf <- coef(fit)
  cat(sprintf("  coef match daily arm: max|diff|=%.3e\n",
              max(abs(cf[paste0("A_",use_s)] - coef(fa)[use_s]))))
  cat(sprintf("  coef match mean  arm: max|diff|=%.3e\n",
              max(abs(cf[paste0("B_",mb)]    - coef(fb)[mb]))))

  cat("  computing joint conley(550) ...\n")
  Vc <- vcov(fit, vcov = conley(550))
  cat("  computing joint cluster(country) ...\n")
  Vk <- vcov(fit, cluster = ~iso3)

  # The diagonal blocks of the joint sandwich must match the marginal variances
  # of the separate fits. If they do not, the stacking changed the estimate
  # rather than only adding the cross term.
  ra <- sqrt(diag(Vc)[paste0("A_",use_s)]) / se(fa)[use_s]
  rb <- sqrt(diag(Vc)[paste0("B_",mb)])    / se(fb)[mb]
  cat(sprintf("  own-block SE ratio (stacked/marginal) daily: min=%.6f max=%.6f\n", min(ra), max(ra)))
  cat(sprintf("  own-block SE ratio (stacked/marginal) mean : min=%.6f max=%.6f\n", min(rb), max(rb)))

  ia <- paste0("A_", use_s); ib <- paste0("B_", mb)
  list(names_a = use_s, names_b = mb,
       coef_a = cf[ia], coef_b = cf[ib],
       V_conley  = Vc[c(ia, ib), c(ia, ib)],
       V_country = Vk[c(ia, ib), c(ia, ib)],
       se_ratio_conley = c(range(ra), range(rb)),
       n_obs = nrow(d), n_country = uniqueN(d$iso3))
}

R_amen <- build_stack("log_amen_HDI", character(0), "grid_2c",   "subnal_ctry", "amen")
R_prod <- build_stack("log_prod_phi", che_vars,     "land_type", "subnal_reg",  "prod")

saveRDS(list(amen = R_amen, prod = R_prod), file.path(DIR_WORK, "step3_vcov.rds"))
cat("\nSaved step3_vcov.rds\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n=== Done ===\n")
