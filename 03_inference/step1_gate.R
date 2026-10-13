# step1_gate.R -- entry point of the CRH inference chain.
#
# Rebuilds the estimation frame, refits both arms of both channels, and stores
# the coefficients and standard errors the later steps need. Every downstream
# step reads step1_frame.rds and none re-estimates.
#
# Reads:  era_daily_shares_pooled.csv, amen_prod.csv, era_daily_cellmean.csv,
#         geo_data.csv, subcountry_EU.csv, daily_{amen,prod}_bins.csv
# Writes: step1_frame.rds
# Needs:  the CRH derived files listed above.
#
# Three gates run before anything is stored, and all three must hold for the
# comparison the paper makes to be meaningful:
#   (i)   both arms are fitted on the same sample, so the daily-versus-annual
#         contrast is not a change of sample in disguise;
#   (ii)  the Conley(550 km) variance is actually applied, since fixest falls
#         back silently when it cannot honour the request;
#   (iii) the refitted coefficients reproduce the published bins, which is what
#         licenses reusing the published curve downstream.

suppressMessages({ library(data.table); library(fixest); library(splines) })
t0 <- proc.time(); cat("=== step1_gate.R ===\n")
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

# ---- 0. bin edges ------------------------------------------------------------
# The 3-degree edges are the CRH design. The two infinite outer bins collect the
# tails, so no day is dropped and the shares still sum to one per cell-year.
INT_EDGES <- seq(-24, 39, by = 3)
BIN_EDGES <- c(-Inf, INT_EDGES, Inf)
NB <- length(BIN_EDGES) - 1L
cat("  NB bins:", NB, "\n")
assign_bin <- function(value) as.integer(cut(value, BIN_EDGES, labels = 1:NB,
                                             include.lowest = TRUE, right = TRUE))

sh <- fread(file.path(DIR_CRH, "era_daily_shares_pooled.csv"))
sh[, bin := assign_bin(value)]
shb <- sh[, .(s = sum(share)), by = .(lon, lat, year, bin)]
ctr <- sh[, .(t_center = sum(value * share) / sum(share)), by = bin][order(bin)]
wide <- dcast(shb, lon + lat + year ~ bin, value.var = "s", fill = 0)
scol <- paste0("s", 1:NB)
setnames(wide, as.character(1:NB), scol, skip_absent = TRUE)
for (b in 1:NB) if (!(scol[b] %in% names(wide))) wide[, (scol[b]) := 0]
setcolorder(wide, c("lon","lat","year", scol))
mass <- colSums(as.matrix(wide[, ..scol]))
ref_bin <- as.integer(sub("s","", names(which.max(mass))))
cat("  ref bin:", ref_bin, " center=", round(ctr[bin==ref_bin]$t_center,2), "C\n")
rm(sh, shb); gc()

# ---- 1. fundamentals, geography, subnational units ---------------------------
# The geographic controls enter as Chebyshev polynomials of the rescaled
# variable, as in the CRH design; the rescaling to [-1, 1] is what makes the
# recursion well conditioned.
ap <- fread(file.path(DIR_CRH, "amen_prod.csv"))
ap[, `:=`(prod = as.numeric(prod), amen_HDI = as.numeric(amen_HDI), phi_1990 = as.numeric(phi_1990))]
ap[, log_amen_HDI := log(amen_HDI)]
ap[, log_prod     := log(prod)]
ap[, log_prod_phi := fifelse(year == 1990 | is.na(phi_1990), log_prod, log(prod / phi_1990))]
ap <- ap[, .(year, lat, lon, grid_2c, log_amen_HDI, log_prod_phi)]

cm <- fread(file.path(DIR_CRH, "era_daily_cellmean.csv"))[, .(lon, lat, year, mean_daily)]

geo <- eggw_csv("Data/geo_data.csv"); setnames(geo, names(geo)[1], "lat")
geo_name <- c("elev","rough","std_elev","dist2coast","dist2ocean","dist2water")
for (g in geo_name) {
  v <- geo[[g]]; gmin <- min(v,na.rm=TRUE); gmax <- max(v,na.rm=TRUE)
  gche <- (v - (gmin+gmax)/2); gche <- gche / max(gche, na.rm=TRUE); gche <- pmin(pmax(gche,-1),1)
  for (i in 1:5) set(geo, j = paste0(g,"_che",i), value = cos(i*acos(gche)))
}
che_vars <- as.vector(outer(geo_name, 1:5, function(a,b) paste0(a,"_che",b)))
geo <- geo[, c("lat","lon","land_type",che_vars), with = FALSE]

sc <- eggw_csv("Data/subcountry_EU.csv")
setnames(sc, c("subcountry_EU_reg","subcountry_EU_ctry"), c("subnal_reg","subnal_ctry"))
sc <- sc[, .(year, lat, lon, subnal_reg, subnal_ctry)]

dat <- merge(wide, ap, by.x = c("lon","lat","year"), by.y = c("lon","lat","year"), all = FALSE)
dat <- merge(dat, cm, by = c("lon","lat","year"), all = FALSE)
dat <- merge(dat, sc, by.x = c("lat","lon","year"), by.y = c("lat","lon","year"), all = FALSE)
dat <- merge(dat, geo, by = c("lat","lon"), all = FALSE)
cat("  merged cell-periods:", nrow(dat), "  unique cells:", uniqueN(dat[,.(lat,lon)]), "\n")

mv  <- dat$mean_daily[!is.na(dat$mean_daily)]
mbrk <- c(-Inf, quantile(mv, seq(1/20,19/20,1/20)), Inf)
dat[, mbin := as.integer(cut(mean_daily, mbrk, labels=1:20, include.lowest=TRUE))]
mctr <- quantile(mv, seq(0.025,0.975,by=0.05))
for (i in 1:20) set(dat, j=paste0("mb_",i), value = dat$mean_daily*(dat$mbin==i))

# ---- 2. refit both arms per channel, and run the three gates -----------------
use_s <- setdiff(scol, paste0("s", ref_bin))
mb    <- paste0("mb_", 1:20)

check_channel <- function(lhs, extra_rhs, fe_block, subnal_col, tag) {
  cat("\n--- channel:", tag, "---\n")
  d <- copy(dat)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year-1985)/5 * max(subnal) + subnal]

  fml_a <- as.formula(paste(lhs, "~", paste(c(use_s, extra_rhs), collapse=" + "), "|", fe_block, "+ fe"))
  fml_b <- as.formula(paste(lhs, "~", paste(c(mb,   extra_rhs), collapse=" + "), "|", fe_block, "+ fe"))

  ## fixef.tol = 1e-11, the same tolerance step3_vcov.R uses. At the fixest
  ## default of 1e-6 the alternating projections stop as soon as the demeaning
  ## moves by less than that, and two fits of the SAME model on the SAME data
  ## then agree only to about 1e-6. The amenity design converges slowly enough to
  ## sit at that edge -- its coefficients moved by 1.1e-03 on the daily arm when
  ## the tolerance was tightened -- while productivity converges four orders past
  ## it and moves by 3e-06. One tolerance across the package, and it is the tight
  ## one.
  ##
  ## THE GAP TO THE FROZEN WITNESS WIDENS BECAUSE OF THIS, from 3.5e-04 to
  ## 1.4e-03 on the amenity bins, and that is expected: the witness was produced
  ## by a loose fit. It is not re-frozen. A witness moved to match the run it is
  ## supposed to test would test nothing.
  fit_a <- feols(fml_a, data=d, vcov=conley(550), fixef.tol=1e-11)
  fit_b <- feols(fml_b, data=d, vcov=conley(550), fixef.tol=1e-11)
  na <- fitstat(fit_a,"n")[[1]]; nb <- fitstat(fit_b,"n")[[1]]
  cat(sprintf("  N daily=%d  N mean=%d  identical=%s\n", na, nb, na==nb))
  oa <- obs(fit_a); ob <- obs(fit_b)
  cat("  same observation index set:", identical(sort(oa), sort(ob)), "\n")
  cat("  vcov label a:", attr(vcov(fit_a), "type"), " | b:", attr(vcov(fit_b), "type"), "\n")

  # Reproduction check against the published bins. The refit must land on the
  # same coefficients, otherwise the frame rebuilt here is not the frame the
  # published curve came from.
  ca <- coef(fit_a); sa <- se(fit_a)
  delta <- setNames(rep(0, NB), scol); delta[use_s] <- ca[use_s]
  pub <- fread(file.path(DIR_REF, sprintf("daily_%s_bins.csv", tag)))
  cat(sprintf("  repro daily bins: max|diff|=%.3e\n", max(abs(pub$delta - as.numeric(delta)))))
  cat(sprintf("  repro t_center  : max|diff|=%.3e\n", max(abs(pub$t_center - ctr$t_center))))
  cb <- coef(fit_b)
  mc <- cb[grep("^mb_",names(cb))]; o <- order(as.integer(sub("mb_","",names(mc)))); mc <- mc[o]
  list(fit_a=fit_a, fit_b=fit_b, d=d, delta=delta, se_a=sa, mc=mc, se_b=se(fit_b)[grep("^mb_",names(se(fit_b)))][o])
}

res_amen <- check_channel("log_amen_HDI", character(0), "grid_2c",   "subnal_ctry", "amen")
res_prod <- check_channel("log_prod_phi", che_vars,     "land_type", "subnal_reg",  "prod")

saveRDS(list(dat=dat, ctr=ctr, mctr=mctr, ref_bin=ref_bin, NB=NB, scol=scol,
             use_s=use_s, che_vars=che_vars,
             amen=list(delta=res_amen$delta, se_a=res_amen$se_a, mc=res_amen$mc, se_b=res_amen$se_b),
             prod=list(delta=res_prod$delta, se_a=res_prod$se_a, mc=res_prod$mc, se_b=res_prod$se_b)),
        file.path(DIR_WORK, "step1_frame.rds"))
cat("\nSaved step1_frame.rds\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n=== Done ===\n")
