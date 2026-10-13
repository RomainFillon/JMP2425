# step7_period_sensitivity.R -- sensitivity of the estimate to the period set.
#
# The 1990 cross-section covers a single year, whereas the other periods cover
# five. This script asks whether that asymmetry drives the estimate, by refitting
# both arms on three period sets:
#   P_all  = {1990, 1995, 2000, 2005}  the primary specification
#   P_no90 = {1995, 2000, 2005}        drops the one-year block
#   P_dec2 = {2000, 2005}              the two periods with full decadal coverage
#
# Every sample-dependent object is rebuilt inside each period set: the reference
# bin, the mean-bin quantile breaks and centres, and the smoothing weights.
# Reusing the full-sample versions would hold fixed quantities that are part of
# what changes, and the sensitivity would be understated.
#
# Reads:  step1_frame.rds, step2_maps.rds, CMIP6_tas.csv
# Writes: crh_period_sensitivity.csv
# Needs:  steps 1 and 2 must have run first.

suppressMessages({ library(data.table); library(fixest); library(splines) })
t0 <- proc.time(); cat("=== step7_period_sensitivity.R ===\n")
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
BASE_YEAR <- 2015L; FUT_YEAR <- 2095L
UNITS <- c("RUS","CAN","FIN","DEU","FRA","NORTH>=45N")

S1 <- readRDS(file.path(DIR_WORK, "step1_frame.rds")); S2 <- readRDS(file.path(DIR_WORK, "step2_maps.rds"))
dat0 <- S1$dat; scol <- S1$scol; NB <- S1$NB; ctr <- S1$ctr; che_vars <- S1$che_vars
cells <- S2$cells; muw <- S2$muw

# Country assignment for the estimation cells, by the same rule as step3.
suppressMessages({ library(sf); library(rnaturalearth) }); sf_use_s2(FALSE)
cl <- unique(dat0[, .(lon,lat)]); pts <- st_as_sf(cl, coords=c("lon","lat"), crs=4326)
world <- ne_countries(scale="medium", returnclass="sf")
icol <- c("iso_a3_eh","adm0_a3","iso_a3","ISO_A3"); icol <- icol[icol %in% names(world)][1]
world <- world[, icol]; names(world)[1] <- "iso3"
cl[, iso3 := suppressMessages(st_join(pts, world, join=st_within))$iso3]
if (any(is.na(cl$iso3)))
  cl[is.na(iso3), iso3 := world$iso3[suppressMessages(st_nearest_feature(pts[is.na(cl$iso3),], world))]]
dat0 <- merge(dat0, cl, by=c("lon","lat"), all.x=TRUE)

# ---- rebuild the cell x integer-temperature delta-share matrix ---------------
cat("Loading CMIP6_tas.csv...\n")
cmip <- fread(file.path(DIR_TMP, "CMIP6_tas.csv"))[year %in% c(BASE_YEAR, FUT_YEAR)]
cmip[, c := count/sum(count), by=.(x,y,year)]
cmip <- merge(cmip, cells, by=c("x","y"))
Tvals <- S2$Tvals; cmip[, tid := match(value, Tvals)]
cmip[, sgn := fifelse(year==FUT_YEAR, 1, -1)]
dc <- cmip[, .(v=sum(sgn*c)), by=.(cid,tid)]
M <- Matrix::sparseMatrix(i=dc$cid, j=dc$tid, x=dc$v, dims=c(nrow(cells), length(Tvals)))
rm(cmip, dc); gc()

pub_cells <- S2$prod$pub                     # cells x (pop, iso3, y) -- identical across channels
A <- sapply(UNITS, function(u) {
  sel <- if (u=="NORTH>=45N") pub_cells$y>=45 else pub_cells$iso3==u; sel[is.na(sel)] <- FALSE
  w <- ifelse(sel, pub_cells$pop, 0); w[is.na(w)] <- 0; w/sum(w) })

daily_map <- function(se_vec, ref_bin) {          # delta -> per-cell daily change, as a matrix
  W <- matrix(0, nrow(cells), NB)
  for (k in 1:NB) {
    e <- rep(0, NB)
    b <- data.table(t_center=ctr$t_center, delta={e[k] <- 1; e}, se=as.numeric(se_vec))
    f <- b[is.finite(t_center)]
    w <- ifelse(is.finite(f$se), 1/f$se, 1/mean(f$se, na.rm=TRUE))
    sp <- lm(delta ~ ns(t_center, df=5), data=f, weights=w)
    tg <- seq(min(f$t_center), max(f$t_center), length.out=400)
    W[,k] <- as.numeric(M %*% approx(tg, as.numeric(predict(sp, newdata=data.frame(t_center=tg))),
                                     xout=Tvals, rule=2)$y)
  }
  W
}
mean_map <- function(se_vec, mctr) {              # mc -> per-cell mean-arm level change
  W <- matrix(0, nrow(cells), 20L)
  for (j in 1:20) {
    e <- rep(0,20); e[j] <- 1
    df <- data.frame(v=e, w=1/as.numeric(se_vec), temp=as.numeric(mctr))
    sp <- lm(v ~ ns(temp, df=4), data=df, weights=w)
    mg <- seq(min(df$temp), max(df$temp), length.out=400)
    ds <- as.numeric(predict(sp, newdata=data.frame(temp=mg)))
    Lam <- c(0, cumsum((head(ds,-1)+tail(ds,-1))/2*diff(mg)))
    f <- function(m) approx(mg, Lam, xout=pmin(pmax(m, min(df$temp)), max(df$temp)), rule=2)$y
    W[,j] <- f(muw$mu_fut) - f(muw$mu_base)
  }
  W
}

run <- function(tag, lhs, extra_rhs, fe_block, subnal_col, KEEP, label, want_vcov) {
  d <- copy(dat0)[year %in% KEEP]
  # Sample-dependent objects, rebuilt inside this period set rather than
  # inherited from the full sample.
  ref_bin <- as.integer(sub("s","", names(which.max(colSums(as.matrix(d[, ..scol]))))))
  mv <- d$mean_daily[!is.na(d$mean_daily)]
  mbrk <- c(-Inf, quantile(mv, seq(1/20,19/20,1/20)), Inf)
  mctr <- quantile(mv, seq(0.025,0.975,by=0.05))
  d[, mbin := as.integer(cut(mean_daily, mbrk, labels=1:20, include.lowest=TRUE))]
  for (i in 1:20) set(d, j=paste0("mb_",i), value = d$mean_daily*(d$mbin==i))
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year-1985)/5 * max(subnal) + subnal]
  us <- setdiff(scol, paste0("s", ref_bin)); mb <- paste0("mb_",1:20)

  fa <- feols(as.formula(paste(lhs,"~",paste(c(us,extra_rhs),collapse=" + "),"|",fe_block,"+ fe")),
              data=d, vcov=conley(550), notes=FALSE)
  fb <- feols(as.formula(paste(lhs,"~",paste(c(mb,extra_rhs),collapse=" + "),"|",fe_block,"+ fe")),
              data=d, vcov=conley(550), notes=FALSE)
  dfull <- setNames(rep(0,NB), scol); dfull[us] <- coef(fa)[us]
  sea <- setNames(rep(NA_real_,NB), scol); sa <- se(fa)[us]; sea[names(sa)] <- sa
  mc <- coef(fb)[mb]; seb <- se(fb)[mb]

  Wd <- daily_map(sea, ref_bin); Wm <- mean_map(seb, mctr)
  lvl_d <- as.numeric(t(t(Wd) %*% A) %*% as.numeric(dfull))
  lvl_m <- as.numeric(t(t(Wm) %*% A) %*% as.numeric(mc))

  ci <- NULL
  if (want_vcov) {   # joint country-clustered vcov on this sample, for the difference
    keep <- c(lhs, us, mb, extra_rhs, fe_block, "fe", "iso3")
    dd <- d[, ..keep][complete.cases(d[, ..keep])]
    dA <- copy(dd)[, arm:="A"]; dB <- copy(dd)[, arm:="B"]; ds <- rbind(dA,dB)
    for (v in us) set(ds, j=paste0("A_",v), value=ds[[v]]*(ds$arm=="A"))
    for (v in mb) set(ds, j=paste0("B_",v), value=ds[[v]]*(ds$arm=="B"))
    for (v in extra_rhs) { set(ds, j=paste0("A_",v), value=ds[[v]]*(ds$arm=="A"))
                           set(ds, j=paste0("B_",v), value=ds[[v]]*(ds$arm=="B")) }
    rhs <- c(paste0("A_",us), paste0("B_",mb))
    if (length(extra_rhs)) rhs <- c(rhs, paste0("A_",extra_rhs), paste0("B_",extra_rhs))
    ft <- feols(as.formula(paste(lhs,"~",paste(rhs,collapse=" + "),"|",paste0(fe_block,"^arm"),"+ fe^arm")),
                data=ds, vcov="iid", fixef.tol=1e-11, notes=FALSE)
    V <- vcov(ft, cluster=~iso3)[c(paste0("A_",us), paste0("B_",mb)), c(paste0("A_",us), paste0("B_",mb))]
    ia <- match(us, scol); pa <- length(us)
    G <- rbind(-(t(Wd) %*% A)[ia,,drop=FALSE], t(Wm) %*% A)
    ci <- sqrt(pmax(0, diag(t(G) %*% V %*% G)))
  }
  data.table(channel=tag, spec=label, n_obs=fitstat(fa,"n")[[1]], ref_bin=ref_bin, unit=UNITS,
             mean=lvl_m, daily=lvl_d, diff=lvl_m-lvl_d,
             atten=1-lvl_d/lvl_m, se_diff=if (is.null(ci)) NA_real_ else ci)
}

SPECS <- list(list(c(1990,1995,2000,2005), "P_all (published primary)"),
              list(c(1995,2000,2005),      "P_no90 (drop the 1-year block)"),
              list(c(2000,2005),           "P_dec2 (two full clean decades)"))
out <- list()
for (s in SPECS) {
  cat("\n--- prod |", s[[2]], "---\n")
  out[[length(out)+1]] <- run("prod","log_prod_phi",che_vars,"land_type","subnal_reg", s[[1]], s[[2]], TRUE)
  cat("--- amen |", s[[2]], "---\n")
  out[[length(out)+1]] <- run("amen","log_amen_HDI",character(0),"grid_2c","subnal_ctry", s[[1]], s[[2]], FALSE)
}
res <- rbindlist(out)
fwrite(res, file.path(DIR_INF, "crh_period_sensitivity.csv"))

for (tg in c("prod","amen")) {
  cat(sprintf("\n================ channel = %s ================\n", tg))
  w <- dcast(res[channel==tg], unit ~ spec, value.var=c("mean","daily","diff"))
  print(w, nrow=20)
}
cat("\n--- prod: shift vs the published SE of the difference ---\n")
ref <- res[channel=="prod" & spec=="P_all (published primary)"]
for (sp in c("P_no90 (drop the 1-year block)","P_dec2 (two full clean decades)")) {
  a <- res[channel=="prod" & spec==sp]
  cat(sprintf("  %s\n", sp))
  for (k in seq_along(UNITS))
    cat(sprintf("    %-11s diff %+.3f -> %+.3f  (shift %+.3f = %.2f x published SE %.3f); se_diff now %.3f\n",
                UNITS[k], ref$diff[k], a$diff[k], a$diff[k]-ref$diff[k],
                abs(a$diff[k]-ref$diff[k])/ref$se_diff[k], ref$se_diff[k], a$se_diff[k]))
}
cat("\nn_obs by spec:\n"); print(unique(res[, .(channel, spec, n_obs, ref_bin)]))
cat("\nSaved crh_period_sensitivity.csv\nTotal time:", round((proc.time()-t0)["elapsed"]/60,1), "min\n")
