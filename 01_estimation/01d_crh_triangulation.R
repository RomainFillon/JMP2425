##################################################################################
# 01d_crh_triangulation.R -- replication and triangulation of Cruz and
#                           Rossi-Hansberg.
#
# WHAT THIS APPENDIX DOES, IN TWO STAGES
#   A. IT REPRODUCES. The Cruz and Rossi-Hansberg replication package is in
#      MATLAB; sections C1 and C2 are its R port. The spatial model is inverted
#      to recover amenities and productivities per cell, then THEIR damage
#      function is re-estimated on THEIR metric. The control is their own
#      published result: 1.6 degrees C for productivity, 8.6 for amenities.
#   B. IT CHANGES ONE THING. Sections C3 to C7 redo the estimation with the
#      distribution of days in place of the annual mean. Same cells, same
#      fundamentals, same fixed effects, same Conley standard errors. The result
#      is then projected onto CMIP6 to see whether the northern gain survives.
#
# THE ORDER CARRIES THE ARGUMENT. Reproduce first, move one thing second.
# Without the first stage the second says nothing: a gap could come from the
# port rather than from the metric.
#
# Reads:
#   EXTERNAL, from the published Cruz and Rossi-Hansberg package, via the zip:
#     H0_areal.mat, pop_Gecon.mat, wage_Gecon.mat, trmult_reduced.mat (2.1 GB),
#     CO2_EDGAR.mat, C.csv, CO2_country.csv, clean_country.csv, HDI_GDPpc.csv,
#     temp_data.csv, geo_data.csv, subcountry_EU.csv
#   NOT REGENERABLE LOCALLY, produced by the server pipeline:
#     ERA_tas_distribution_grid_annual.csv, CMIP6_tas.csv
#   PRODUCED BY A SEPARATE CHAIN:
#     crh_inference_deliverable.csv, which carries the bootstrap intervals
#   A referee can download the EXTERNAL inputs; nothing can be done about the
#   server ones. Presenting the two alike would mislead by uniformity.
#
# Writes:
#   .tex   tab_crh_northern.tex
#   .pdf   fig_drf_crh_annual_vs_daily_{prod,amenity}.pdf
#          (fig_crh_daily_response.pdf est remplacee par fig_crh_responses.pdf,
#           produite par 01d_fig_crh_responses.R)
#   .csv   crh_design_optimum_{prod,amen}.csv, amen_prod.csv,
#          {amen,prod}_coeff_10y_1h_20b_550d.csv,
#          tab_crh_northern_reversal_daily{,_dec2}.csv,
#          crh_daily_projection_cells.csv,
#          (crh_daily_response_peaks.csv est passe a 01d_fig_crh_responses.R)
#          daily_{amen,prod}_{bins,curve}.csv, mean_{amen,prod}_{curve,coeff}.csv,
#          era_daily_{shares_pooled,cellmean}.csv, crh_metric_guardrail.csv,
#          crh_support_bounds.csv, crh_outofsupport_mass.csv
#
# MATLAB is NOT required: C1 and C2 port the code to R and read the data files of
# the published package directly.
##################################################################################

suppressMessages({
  library(data.table); library(fixest); library(splines)
  library(hdf5r);      library(R.matlab)
})

## The root is derived from the location of this script, not from the working
## directory and not from a required environment variable: neither of those is
## reliable across the ways R can be launched. JMP_ROOT wins when it is set.
.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .self <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                    error = function(e) NA_character_)
  if (is.na(.self)) {
    .a <- commandArgs(trailingOnly = FALSE)
    .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
    if (length(.f)) .self <- normalizePath(.f[1], winslash = "/")
  }
  if (is.na(.self))
    stop("Impossible de deduire la racine. Poser JMP_ROOT, ou lancer par\n",
         "  Rscript depuis la racine du depot.")
  .root <- dirname(dirname(.self))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))

## ---- Lot 57 (28/09/2026) : OBSERVABILITY, not a correction -----------------
## 01d_bis died twice in C1, period 2, with code 1 and no message. Two tools so
## that the next death says more:
##   1. every console write is flushed at once. `cat` is masked by a version
##      that flushes after each call, and a task callback flushes after each
##      top-level expression (the printed gc() tables). Nothing else changes.
##   2. mem_releve() prints, with the clock, the R heap (sum of gc() Mb used),
##      the process RSS, and the free memory of the machine (ps package). It is
##      wrapped in tryCatch: an instrument must never be what kills the run.
cat <- function(...) {
  base::cat(...)
  flush(stdout()); flush(stderr())
  invisible(NULL)
}
invisible(addTaskCallback(function(...) { flush(stdout()); TRUE }))
mem_releve <- function(etiquette) {
  ligne <- tryCatch({
    tas <- sum(gc(verbose = FALSE)[, 2])
    sm  <- ps::ps_system_memory()
    rss <- ps::ps_memory_info()[["rss"]]
    sprintf("  [mem %s] %s | tas R %.0f Mo | RSS %.0f Mo | machine libre %.1f / %.1f Go",
            format(Sys.time(), "%H:%M:%S"), etiquette, tas, rss / 2^20,
            sm$avail / 2^30, sm$total / 2^30)
  }, error = function(e) sprintf("  [mem %s] %s | releve impossible : %s",
                                 format(Sys.time(), "%H:%M:%S"), etiquette,
                                 conditionMessage(e)))
  writeLines(ligne)
  flush(stdout())
}
mem_releve("début du script")


if (!dir.exists(DIR_CRH)) dir.create(DIR_CRH, recursive = TRUE)
t0 <- proc.time()

##################################################################################
# 0. PARAMETERS OF THE CRUZ AND ROSSI-HANSBERG MODEL
#
# Taken from their calibration, not estimated here.
# theta and Omega are the SAME SYMBOLS as in the trade and migration channels of
# this paper, with values of their own. Everything is therefore prefixed, so
# that they cannot be confused with the parameters of 00_parameters.R.
##################################################################################

CRH <- list(
  lambda  = 0.32,   gamma1 = 0.319,  gamma2 = 0.99246,
  eta     = 1,      mu     = 0.8,    nu     = 0.15,
  ksi     = 125,    alpha  = 0.06,   theta  = 6.5,
  epsilon = 1.6,    Omega  = 0.5,    psi_util = 0.045)

CRH_YEARS       <- c(1990L, 1995L, 2000L, 2005L)
CRH_LOGIS_STEEP <- 0.1   # logistic slope, fixed as in damage_function.do

## The years 1985-1989 are EXCLUDED, and that is not a sample choice.
## The server binned those years HOURLY rather than daily, so a cell-year holds
## about 24 times as many values. The MEANS are preserved, since the mean of
## hourly values is the mean of daily ones, so a guardrail on the mean would
## PASS. But the SHAPE of the distribution is inflated by the diurnal cycle,
## which is fatal for a specification whose object IS the shape.
CRH_PERIODS <- list(`1990` = 1990L, `1995` = 1990:1995,
                    `2000` = 1991:2000, `2005` = 1996:2005)
CRH_PERIODS_DEC2 <- c(2000L, 2005L)   # the two full and clean decades

## Interior bins of 3 degrees C from -24 to +39: the weighted day mass lives in
## [-26, +30], and the hot tail is resolved up to +39 for the projection.
CRH_INT_EDGES <- seq(-24, 39, by = 3)
CRH_BIN_EDGES <- c(-Inf, CRH_INT_EDGES, Inf)
CRH_NB        <- length(CRH_BIN_EDGES) - 1L

## CRH_BASE_YEAR, CRH_FUT_YEAR, CRH_COLD et CRH_NORTH_LAT sont declares dans
## 00_parameters.R depuis le 18/09/2026 : ils ont trois lecteurs et non plus un.
## On VERIFIE qu'ils sont poses plutot que de les redeclarer -- une redeclaration
## ici redeviendrait la copie qu'on vient de retirer.
stopifnot(exists("CRH_BASE_YEAR"), exists("CRH_FUT_YEAR"),
          exists("CRH_COLD"), exists("CRH_NORTH_LAT"))
CRH_CONLEY_KM <- 550

## PRE-REGISTERED rule for the well-populated support: the edge bin is the last
## one whose count of cell-periods with positive day mass reaches
## max(1 per cent of the cell-periods, 30). It is fixed before being applied,
## which is what makes it a criterion rather than an adjustment.
## ONE DEMEANING TOLERANCE ACROSS THE PACKAGE: fixef.tol = 1e-11 on every feols
## here, as in 03_inference/. At the fixest default of 1e-6 the alternating
## projections stop as soon as the demeaning moves by less than that, so two
## fits of the SAME model on the SAME data agree only to about 1e-6. The
## amenity design converges slowly enough to sit at that edge; productivity
## converges four orders past it. The provenance guard of section C7 compares
## exactly two such fits, so a loose tolerance made it unsatisfiable.

CRH_SUPPORT_MIN_N   <- 30L
CRH_SUPPORT_MIN_PCT <- 0.01

##################################################################################
# 0bis. THE INPUTS, AND THEIR TWO ORIGINS
#
# The Cruz and Rossi-Hansberg package is PUBLIC AND CITABLE (Zenodo, DOI
# 10.5281/zenodo.7681424). It is not copied into this repository: the CSV files
# are read from INSIDE the archive, and only the .mat files are extracted,
# because hdf5r and R.matlab open a path and cannot read from a connection.
#
# The two climate files, by contrast, come from the SERVER and can be regenerated
# by no one else. Treating them like the public package would mislead by
# uniformity: a referee can download one and not the others.
##################################################################################

source(file.path(.jmp_root, "01_estimation/lib_eggw.R"))
require_eggw()

SERVER_IN <- c(
  era_grid = file.path(DIR_TMP, "ERA_tas_distribution_grid_annual.csv"),
  cmip6    = file.path(DIR_TMP, "CMIP6_tas.csv"))
miss <- SERVER_IN[!file.exists(SERVER_IN)]
if (length(miss))
  stop("⛔ Intrant(s) serveur absent(s) — non régénérables localement :\n  ",
       paste(sprintf("%-9s %s", names(miss), miss), collapse = "\n  "),
       "\n  Produits sur Inari par 01_era5_distribution.R (mode grid) et la ",
       "chaîne CMIP6.")

##################################################################################
# ================= A. ON REPRODUIT CRUZ & ROSSI-HANSBERG ======================
##################################################################################

##################################################################################
# C1. INVERSION OF THE SPATIAL MODEL -- R port of fundamentals.m
#
# Recovers, per cell and per period, the amenity and the productivity that
# rationalise the observed population and wages. It is a nested fixed point:
# amenities, then productivity, then population, until the implied population
# reproduces the observed one.
#
# trmult_reduced.mat is 2.1 GB and takes about a minute to load.
##################################################################################

cat("\n=== C1 : inversion du modèle spatial (portage de fundamentals.m) ===\n")

H0   <- read_mat_h5(eggw_mat("H0_areal.mat"),   "H0",   c(180, 360))
popf <- read_mat_h5(eggw_mat("pop_Gecon.mat"),  "pop",  c(180, 360, 4))
wagf <- read_mat_h5(eggw_mat("wage_Gecon.mat"), "wage", c(180, 360, 4))

mask <- H0 > 0
H <- as.vector(H0[mask]); n <- length(H)
cat(sprintf("  cellules terrestres : n = %d\n", n))

Cc  <- as.vector(as.matrix(eggw_csv("Data/C.csv", header = FALSE)))[mask]
n_C <- length(unique(Cc))
pop_cell  <- sapply(1:4, function(t) as.vector(popf[, , t])[mask])
wage_cell <- sapply(1:4, function(t) as.vector(wagf[, , t])[mask])
wbar <- mean(wage_cell[, 3])

## ---- energy: fossil and clean shares, and the price per cell ---------------
tCO2_toe <- 2.8466
fossil0  <- 25.284967 / tCO2_toe
clean0_g <- 1.3179
GDP0 <- sum(wage_cell[, 3] * pop_cell[, 3])
pf0  <- 72.99e9 * tCO2_toe
pc0  <- 1.15 * 76.34e9
eps  <- CRH$epsilon
fossil_share <- (1 + (pc0 / pf0) * (clean0_g / fossil0)^(1 / eps))^(-1)
pe0 <- (fossil_share^eps * pf0^(1 - eps) +
        (1 - fossil_share)^eps * pc0^(1 - eps))^(1 / (1 - eps))
energy0 <- (fossil_share * fossil0^((eps - 1) / eps) +
            (1 - fossil_share) * clean0_g^((eps - 1) / eps))^(eps / (eps - 1))
chi <- 1 - (pe0 * energy0 * (CRH$mu + CRH$gamma1 / CRH$ksi) / GDP0) / CRH$mu

## CO2_EDGAR is a v7 .mat, not v7.3: readMat, never hdf5r.
co2f  <- readMat(eggw_mat("CO2_EDGAR.mat"))[["CO2.EDGAR"]]
fos_c <- as.matrix(eggw_csv("Data/CO2_country.csv",   header = FALSE))
cln_c <- as.matrix(eggw_csv("Data/clean_country.csv", header = FALSE))
fossil_cell <- clean_cell <- matrix(0, n, 4)

## TWO SOURCES OF NA, and each breaks the `if` differently.
##   co2f : the EDGAR raster carries NA on cells with no data. Without na.rm,
##          sum() propagates them and `if (NA > 0)` stops.
##   Cc   : if a country code is NA for a cell, `Cc == i` returns NA and the
##          logical mask becomes invalid. which() drops them, `==` does not.
## Treating the EDGAR NA as ZERO emissions is the right move: a cell with no data
## emits nothing that is known, and the proportional split that follows would
## give it a zero share anyway.
cat(sprintf("  co2f : %d valeurs non finies | Cc : %d NA sur %d cellules\n",
            sum(!is.finite(co2f)), sum(is.na(Cc)), length(Cc)))

for (t in 1:4) {
  ct <- as.vector(co2f[, , t])[mask]
  ct[!is.finite(ct)] <- 0
  for (i in 1:n_C) {
    idx <- which(Cc == i)
    if (!length(idx)) next
    s <- sum(ct[idx])
    if (s > 0) {
      fossil_cell[idx, t] <- ct[idx] * fos_c[i, t] / s
      clean_cell[idx, t]  <- ct[idx] * cln_c[i, t] / s
    }
  }
}

energy_cell <- (fossil_share * fossil_cell^((eps - 1) / eps) +
                (1 - fossil_share) * clean_cell^((eps - 1) / eps))^(eps / (eps - 1))
safe <- function(x) ifelse(is.nan(x) | is.infinite(x) | x <= 0, 1e-20, x)
kk <- CRH$mu * (1 - chi) / (CRH$mu + CRH$gamma1 / CRH$ksi)
pf_cell <- kk * fossil_share * (pop_cell / H) / safe(energy_cell) *
           safe(energy_cell / fossil_cell)^(1 / eps)
pc_cell <- kk * (1 - fossil_share) * (pop_cell / H) / safe(energy_cell) *
           safe(energy_cell / clean_cell)^(1 / eps)
pe_cell <- (fossil_share^eps * safe(pf_cell)^(1 - eps) +
            (1 - fossil_share)^eps * safe(pc_cell)^(1 - eps))^(1 / (1 - eps))

## ---- Constantes du point fixe ---------------------------------------------------
th <- CRH$theta; lam <- CRH$lambda
denom <- 1 + 2 * th
squ   <- CRH$alpha - 1 + (lam + CRH$gamma1 / CRH$ksi - (1 - CRH$mu)) * th
flatL <- lam - squ / denom
flatR <- 1 - lam * th + (1 + th) * squ / denom
const_N <- wbar^((1 + 2 * th) / th)

cat("  chargement de trmult_reduced.mat (2,1 Go, ~1 min)...\n")
mem_releve("avant chargement de trmult")
ftr <- H5File$new(eggw_mat("trmult_reduced.mat"), mode = "r")
trmult <- ftr[["trmult_reduced"]][1:n, 1:n]; ftr$close_all()

## ---- nested fixed point, four periods --------------------------------------
amen_util_H0 <- prod_H0 <- realgdp_H0 <- matrix(0, n, 4)
amen_util_i <- rep(1, n)
mem_releve("trmult chargée, avant la boucle C1 : 4 périodes (per in 1:4)")
for (per in 1:4) {
  mem_releve(sprintf("période %d/4 (%d) début", per, CRH_YEARS[per]))
  pop0 <- pop_cell[, per]; w0 <- wage_cell[, per]; pe0c <- pe_cell[, per]
  lbar <- sum(pop0); pd <- pop0 / H; l_i <- pd
  NL <- const_N * w0^(-1) * pd^lam
  NR <- w0^(1 + th) * pd^(1 - lam * th) * H
  LL_aux <- H^(-1 / (flatL * denom)) *
            pe0c^(-th * (1 - chi) * CRH$mu / (flatL * denom))
  LR_aux <- H^(th / denom) * pe0c^(-th * (1 + th) * (1 - chi) * CRH$mu / denom)
  prod_aux <- const_N^(-th) * H * w0^(1 + 2 * th) * pd^(-squ) *
              pe0c^(th * (1 - chi) * CRH$mu)
  step <- 0.4; err_o <- Inf
  while (err_o >= 1e-12) {
    step <- step / 2
    err1 <- Inf
    while (err1 >= 1e-15) {
      af <- NL * as.vector(trmult %*% (NR * amen_util_i^(th - step)))^(-1 / th)
      af[!is.finite(af)] <- amen_util_i[!is.finite(af)]
      err1 <- sum((af - amen_util_i)^2)
      amen_util_i <- 0.5 * af + 0.5 * amen_util_i
    }
    prod_v <- prod_aux * amen_util_i^th
    LL <- amen_util_i^((1 + th) / (flatL * denom)) *
          prod_v^(1 / (flatL * denom)) * LL_aux
    LR <- amen_util_i^(th^2 / denom) * prod_v^((1 + th) / denom) * LR_aux
    err2 <- Inf
    while (err2 >= 1e-15) {
      lf <- LL * as.vector(trmult %*% (LR * l_i^flatR))^(1 / (flatL * th))
      lf[!is.finite(lf)] <- l_i[!is.finite(lf)]
      err2 <- sum(((H * l_i - H * lf) / sum(H * lf) * lbar)^2)
      l_i <- 0.5 * lf + 0.5 * l_i
    }
    l_i <- l_i * lbar / sum(H * l_i)
    err_o <- sum((H * l_i - H * pd)^2)
  }
  amen_util_H0[, per] <- amen_util_i
  prod_H0[, per]      <- prod_v
  realgdp_H0[, per]   <- amen_util_i^(-1) * pd^lam
  cat(sprintf("  période %d (%d) : convergée\n", per, CRH_YEARS[per]))
  mem_releve(sprintf("période %d/4 (%d) fin", per, CRH_YEARS[per]))
}

## ---- forward simulation 1990-2005: the knowledge stock phi -----------------
hdi_data <- as.data.frame(eggw_csv("Data/HDI_GDPpc.csv"))
HDI_all <- hdi_data[, 4]; year_all <- hdi_data[, 1]
HDI_mat <- sapply(CRH_YEARS, function(y) HDI_all[year_all == y][1:n])
u_HDI <- exp(HDI_mat^3 / CRH$psi_util)

phi_sim <- matrix(1, n, 16); prod_sim <- matrix(0, n, 16)
prod_sim[, 1] <- prod_H0[, 1]
for (t in 1:15) {
  phi_sim[, t + 1] <- phi_sim[, t]
  if (t < 15) {
    avg <- mean(prod_sim[, t])
    prod_sim[, t + 1] <- CRH$eta * prod_sim[, t]^CRH$gamma2 *
                         avg^(1 - CRH$gamma2) * phi_sim[, t]^(CRH$gamma1 * th)
  }
}
phi_summ <- phi_sim[, c(5, 10, 15)]

## ---- Assemblage d'amen_prod.csv --------------------------------------------------
lat_v <- matrix(rep(seq(89.5, -89.5, by = -1), 360), 180, 360)[mask]
lon_v <- matrix(rep(seq(-179.5, 179.5, by = 1), each = 180), 180, 360)[mask]
g2 <- kronecker(matrix(1:(90 * 180), 90, 180), matrix(1, 2, 2)); g2[H0 == 0] <- NA
HDI0 <- HDI_all[year_all == 2000][1:n]

amen_prod <- data.table(
  year = rep(CRH_YEARS, each = n), lat = rep(lat_v, 4), lon = rep(lon_v, 4),
  amen_HDI = as.vector(u_HDI * amen_util_H0), prod = as.vector(prod_H0),
  phi_1990 = c(rep(NA_real_, n), as.vector(phi_summ)),
  grid_2c = rep(g2[mask], 4),
  subcountry = rep(match(HDI0, sort(unique(HDI0))), 4),
  pop_x = as.vector(pop_cell), wage_x = as.vector(wage_cell),
  realgdp_x = as.vector(realgdp_H0), size_cell_x = rep(H, 4))
fwrite(amen_prod, file.path(DIR_CRH, "amen_prod.csv"))
cat(sprintf("  amen_prod.csv : %d lignes\n", nrow(amen_prod)))
rm(trmult); gc(verbose = FALSE)

##################################################################################
# C2. THE CRH DAMAGE FUNCTION, RE-ESTIMATED
#
# Twenty quantile bins of the decadal hemispheric mean, smoothed by a logistic
# with a fixed slope. The optimum is the root of that logistic.
#
# THIS IS THE CONTROL ON THE PORT. Cruz and Rossi-Hansberg publish 1.6 degrees C
# for productivity and 8.6 for amenities. If the optima here depart far from
# those, the port is wrong and the triangulation says nothing. Both values are
# EMITTED: they previously existed only in a console log, and a comparison with
# only one side verifiable is not publishable.
##################################################################################

cat("\n=== C2 : la DRF de CRH, réestimée ===\n")

## Hemispheric temperature: January in the north, July in the south, on a ten
## year rolling mean. The annual mean of the same object serves as the reference
## for the metric guardrail of C3.
td <- eggw_csv("Data/temp_data.csv")
td[, cell := .GRP, by = .(lat, lon)]; setorder(td, cell, year)
mcols <- grep("^temp_1y_", names(td), value = TRUE)
td[, temp_1y_1a := rowMeans(.SD, na.rm = TRUE), .SDcols = mcols]
td[, temp_10y_1m := frollmean(temp_1y_1m, 10, align = "right"), by = cell]
td[, temp_10y_7m := frollmean(temp_1y_7m, 10, align = "right"), by = cell]
td[, temp_10y_1a := frollmean(temp_1y_1a, 10, align = "right"), by = cell]
td[, temp_10y_1h := fifelse(lat >= 0, temp_10y_1m, temp_10y_7m)]
temp_crh <- td[year %in% CRH_YEARS, .(year, lat, lon, temp_10y_1h, temp_10y_1a)]
rm(td); gc(verbose = FALSE)

geo <- eggw_csv("Data/geo_data.csv"); setnames(geo, names(geo)[1], "lat")
subc <- eggw_csv("Data/subcountry_EU.csv")
setnames(subc, c("subcountry_EU_reg", "subcountry_EU_ctry"),
         c("subnal_reg", "subnal_ctry"), skip_absent = TRUE)

## Chebyshev polynomials of order 5 on six geographic variables: the control of
## the CRH productivity specification, taken as is.
GEO_VARS <- c("elev", "rough", "std_elev", "dist2coast", "dist2ocean", "dist2water")
for (g in GEO_VARS) {
  v <- geo[[g]]; r <- range(v, na.rm = TRUE)
  z <- v - mean(r); z <- pmin(pmax(z / max(z, na.rm = TRUE), -1), 1)
  for (i in 1:5) set(geo, j = paste0(g, "_che", i), value = cos(i * acos(z)))
}
CHE <- as.vector(outer(GEO_VARS, 1:5, function(a, b) paste0(a, "_che", b)))
geo <- geo[, c("lat", "lon", "land_type", CHE), with = FALSE]

crh_base <- merge(copy(amen_prod), temp_crh, by = c("year", "lat", "lon"))
crh_base <- merge(crh_base, subc[, .(year, lat, lon, subnal_reg, subnal_ctry)],
                  by = c("year", "lat", "lon"))
crh_base <- merge(crh_base, geo, by = c("lat", "lon"))
crh_base[, `:=`(log_amen_HDI = log(as.numeric(amen_HDI)),
                log_prod_phi = fifelse(year == 1990 | is.na(phi_1990),
                                       log(as.numeric(prod)),
                                       log(as.numeric(prod) / as.numeric(phi_1990))))]
cat(sprintf("  panel CRH : %d cellules-périodes\n", nrow(crh_base)))
if (nrow(crh_base) == 0L)
  stop("⛔ panel CRH vide après jointures : vérifier les clés lat/lon/year.")

fit_logistic <- function(v, se, temp) {
  d <- data.frame(v = v, se = se, temp = temp)
  fn <- function(t, lo, hi, ctr) lo + (hi - lo) /
        (1 + exp(CRH_LOGIS_STEEP * (t - ctr)))
  st <- list(lo = min(v, na.rm = TRUE), hi = max(v, na.rm = TRUE), ctr = 3.75)
  fit <- tryCatch(
    nls(v ~ fn(temp, lo, hi, ctr), data = d, start = st, weights = 1 / se,
        control = nls.control(tol = 1e-5, maxiter = 200)),
    error = function(e) {
      cat("    NLS tol=1e-5 a échoué, reprise à tol=1e-2 :", conditionMessage(e), "\n")
      nls(v ~ fn(temp, lo, hi, ctr), data = d, start = st, weights = 1 / se,
          control = nls.control(tol = 1e-2, maxiter = 500)) })
  p <- coef(fit)
  list(p = p, T_opt = uniroot(function(t) fn(t, p[["lo"]], p[["hi"]], p[["ctr"]]),
                              interval = c(-50, 50), extendInt = "yes")$root)
}

run_crh_baseline <- function(lhs, extra, fe_block, subnal_col, tag, published) {
  cat(sprintf("  -- canal %s --\n", tag))
  d <- copy(crh_base)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year - 1985) / 5 * max(subnal) + subnal]
  tv <- d$temp_10y_1h[!is.na(d$temp_10y_1h)]
  brk <- c(-Inf, quantile(tv, seq(1/20, 19/20, 1/20)), Inf)
  d[, bin20 := as.integer(cut(temp_10y_1h, brk, labels = 1:20, include.lowest = TRUE))]
  ctr <- quantile(tv, seq(0.025, 0.975, by = 0.05))
  for (i in 1:20) set(d, j = paste0("tb_", i), value = d$temp_10y_1h * (d$bin20 == i))

  f <- as.formula(paste(lhs, "~",
        paste(c(paste0("tb_", 1:20), extra), collapse = " + "), "|", fe_block, "+ fe"))
  fit <- tryCatch(feols(f, data = d, vcov = conley(CRH_CONLEY_KM), fixef.tol = 1e-11),
                  error = function(e) {
                    cat("    Conley a échoué, repli sur cluster :",
                        conditionMessage(e), "\n")
                    feols(f, data = d, cluster = ~ get(fe_block), fixef.tol = 1e-11) })
  cb <- coef(fit)[grep("^tb_", names(coef(fit)))]
  sb <- se(fit)[grep("^tb_", names(se(fit)))]
  o <- order(as.integer(sub("tb_", "", names(cb))))
  L <- fit_logistic(as.numeric(cb[o]), as.numeric(sb[o]), as.numeric(ctr))

  cat(sprintf("    optimum obtenu %.2f °C | CRH publient %.1f °C | écart %+.2f\n",
              L$T_opt, published, L$T_opt - published))
  fwrite(data.table(channel = tag, T_opt = L$T_opt, crh_published = published,
                    source = "05_crh_triangulation.R :: run_crh_baseline"),
         file.path(DIR_OUT, sprintf("crh_design_optimum_%s.csv", tag)))
  short <- ifelse(tag == "prod", "prod", "amen")
  fwrite(data.table(coeff = c("logis_min", "logis_max", "logis_center", "logis_steep"),
                    value = c(L$p[["lo"]], L$p[["hi"]], L$p[["ctr"]], CRH_LOGIS_STEEP)),
         file.path(DIR_CRH, sprintf("%s_coeff_10y_1h_20b_550d.csv", short)))
  fwrite(data.table(bin = 1:20, temp = as.numeric(ctr),
                    v = as.numeric(cb[o]), se = as.numeric(sb[o])),
         file.path(DIR_CRH, sprintf("%s_bins_10y_1h_20b.csv", short)))
  invisible(L$T_opt)
}

invisible(run_crh_baseline("log_amen_HDI", character(0), "grid_2c",
                           "subnal_ctry", "amen", 8.6))
invisible(run_crh_baseline("log_prod_phi", CHE, "land_type",
                           "subnal_reg", "prod", 1.6))

##################################################################################
# ============= B. ONE THING CHANGES: THE METRIC ===============================
##################################################################################

##################################################################################
# C3. POOLING THE DAYS, AND THE METRIC COHERENCE GUARDRAIL
#
# THE GUARDRAIL IS THE CONDITION OF THE EXERCISE. Per period, the mean
# reconstructed from the day shares must track the decadal mean of CRH. A nearly
# constant offset is acceptable, being the same systematic bias as the annual
# non-regression; a wide dispersion would signal an aggregation problem upstream,
# and the whole bins-against-bins exercise would fall with it.
##################################################################################

cat("\n=== C3 : mise en commun des jours + garde-fou ===\n")

era <- fread(SERVER_IN[["era_grid"]])          # x, y, year, value, count
ap  <- copy(amen_prod)
setnames(era, c("x", "y"), c("lon", "lat"), skip_absent = TRUE)

## Restricted to the support of the fundamentals: the exercise only has meaning
## on the cells where CRH inverted their model.
cells_ok <- unique(ap[, .(lon, lat)])
era <- merge(era, cells_ok, by = c("lon", "lat"))

pool <- rbindlist(lapply(names(CRH_PERIODS), function(p) {
  e <- era[year %in% CRH_PERIODS[[p]]]
  e[, .(count = sum(count)), by = .(lon, lat, value)][, period := as.integer(p)][]
}))
pool[, share := count / sum(count), by = .(lon, lat, period)]
setnames(pool, "period", "year")
shares <- pool[, .(lon, lat, year, value, share)]
cellmean <- pool[, .(mean_daily = sum(value * share), ndays = sum(count)),
                 by = .(lon, lat, year)]
fwrite(shares,   file.path(DIR_CRH, "era_daily_shares_pooled.csv"))
fwrite(cellmean, file.path(DIR_CRH, "era_daily_cellmean.csv"))

{
  crh_t <- temp_crh[, .(lon, lat, year, temp_10y_1a)]
  g <- merge(cellmean, crh_t[, .(lon, lat, year, temp_10y_1a)],
             by = c("lon", "lat", "year"))
  if (nrow(g) == 0L) stop("garde-fou : aucune cellule-période appariée à temp_data.csv.")
  bias <- mean(g$mean_daily - g$temp_10y_1a, na.rm = TRUE)
  sdev <- sd(g$mean_daily - g$temp_10y_1a, na.rm = TRUE)
  rho  <- cor(g$mean_daily, g$temp_10y_1a, use = "complete.obs")
  slp  <- coef(lm(mean_daily ~ temp_10y_1a, data = g))[2]
  cat(sprintf("  cohérence métrique : r = %.4f | pente = %.4f | biais = %+.3f °C | sd = %.3f\n",
              rho, slp, bias, sdev))
  ## It DIAGNOSES, it does not stop. Putting a numeric threshold here would
  ## substitute an automatic criterion for a judgement the original left to the
  ## reader.
  cat("  (un biais quasi constant est attendu ; une dispersion large invaliderait\n")
  cat("   l'exercice « on ne change QUE la métrique » — à lire, pas à ignorer.)\n")
  fwrite(data.table(cor = rho, slope = slp, bias = bias, sd = sdev),
         file.path(DIR_CRH, "crh_metric_guardrail.csv"))

  ## The same four numbers, under the names and at the address the macro emitter
  ## reads. They used to be produced outside the package, by a script that also
  ## rewrote era_daily_{shares_pooled,cellmean}.csv -- files this script already
  ## writes twenty lines above. Migrating that script whole would have given two
  ## producers for one artefact, so only the coherence figures come across, and
  ## they are not recomputed: they are the ones printed just above.
  fwrite(data.table(n = nrow(g), cor = rho, bias_era_crh = bias,
                    sd_diff = sdev, lm_slope = unname(slp)),
         file.path(DIR_OUT, "crh_metric_coherence.csv"))
  cat(sprintf("  [write] crh_metric_coherence.csv (n = %d cellules-periodes)\n", nrow(g)))
}

##################################################################################
# C4. RE-ESTIMATION -- TWO METRICS, ONE SAMPLE
#
#   (a) DAILY: log(x) ~ sum_k delta_k * s_k(r,t), plus Chebyshev for productivity
#       where s_k is the share of the cell-period days falling in bin k.
#       The REFERENCE bin, the one with maximum day mass, is dropped from the
#       right-hand side: the shares sum to one and would be collinear with the
#       constant of the fixed effects. The projection is invariant to the choice
#       of reference, since it enters only through differences.
#   (b) MEAN: log(x) ~ sum_i [mean * 1{bin_i}], the exact form of CRH, except
#       that the scalar is the annual mean OF THE DAILY DISTRIBUTION, hence in
#       THE SAME UNIT as (a). Smoothed into a level curve by a natural spline and
#       a trapezoidal integral.
#
# Fixed effects and controls, identical across channels:
#   amenity       log_amen_HDI | grid_2c   + year x subcountry, no geography
#   productivity  log_prod_phi | land_type + year x subregion, Chebyshev order 5
##################################################################################

cat("\n=== C4 : ré-estimation bins-contre-bins ===\n")

assign_bin <- function(v) as.integer(cut(v, CRH_BIN_EDGES, labels = 1:CRH_NB,
                                         include.lowest = TRUE, right = TRUE))

build_panel <- function(keep_periods = NULL) {
  sh <- copy(shares)
  sh[, bin := assign_bin(value)]
  ctr <- sh[, .(t_center = sum(value * share) / sum(share)), by = bin][order(bin)]
  wide <- dcast(sh[, .(s = sum(share)), by = .(lon, lat, year, bin)],
                lon + lat + year ~ bin, value.var = "s", fill = 0)
  scol <- paste0("s", 1:CRH_NB)
  setnames(wide, as.character(1:CRH_NB), scol, skip_absent = TRUE)
  for (b in 1:CRH_NB) if (!(scol[b] %in% names(wide))) wide[, (scol[b]) := 0]
  if (!is.null(keep_periods)) wide <- wide[year %in% keep_periods]

  ## The panel of section C2 is REUSED: same fundamentals, same geography, same
  ## fixed effects. Rebuilding it here would risk letting the sample diverge
  ## between the two stages, which is exactly what the exercise forbids: only the
  ## metric may change.
  d <- merge(wide, crh_base[, c("lon", "lat", "year", "grid_2c", "land_type",
                                "subnal_reg", "subnal_ctry",
                                "log_amen_HDI", "log_prod_phi", CHE),
                            with = FALSE],
             by = c("lon", "lat", "year"))
  d <- merge(d, cellmean[, .(lon, lat, year, mean_daily)],
             by = c("lon", "lat", "year"))
  if (nrow(d) == 0L)
    stop("⛔ panel bins-contre-bins vide : les parts de jours ne s'apparient pas ",
         "au panel de la section C2.")

  ## Reference bin = maximum day mass, DERIVED from the data and not chosen.
  mass <- colSums(as.matrix(d[, ..scol]))
  list(d = d, ctr = ctr, scol = scol,
       ref_bin = as.integer(sub("s", "", names(which.max(mass)))))
}

## Smoothing of the mean response: weighted natural spline, then a trapezoidal
## integral to go from the derivative to the level Lambda(mu).
smooth_mean_drf <- function(v, se, temp) {
  sp <- lm(v ~ ns(temp, df = 4), weights = 1 / se)
  g <- seq(min(temp), max(temp), length.out = 400)
  ds <- as.numeric(predict(sp, newdata = data.frame(temp = g)))
  Lam <- c(0, cumsum((head(ds, -1) + tail(ds, -1)) / 2 * diff(g)))
  cross <- which(diff(sign(ds)) != 0)
  list(mu = g, delta = ds, Lambda = Lam,
       Topt = if (length(cross)) g[cross[1]] else NA_real_, support = range(temp))
}

estimate_channel <- function(PN, lhs, extra_rhs, fe_block, subnal_col, tag, suffix = "") {
  cat(sprintf("  -- canal %s%s --\n", tag, suffix))
  d <- copy(PN$d)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year - 1985) / 5 * max(subnal) + subnal]

  ## (a) DAILY -- the reference bin is dropped.
  use_s <- setdiff(PN$scol, paste0("s", PN$ref_bin))
  f_a <- as.formula(paste(lhs, "~", paste(c(use_s, extra_rhs), collapse = " + "),
                          "|", fe_block, "+ fe"))
  fit_a <- tryCatch(feols(f_a, data = d, vcov = conley(CRH_CONLEY_KM), fixef.tol = 1e-11),
                    error = function(e) {
                      cat("    Conley a échoué, repli sur cluster :",
                          conditionMessage(e), "\n")
                      feols(f_a, data = d, cluster = ~ get(fe_block), fixef.tol = 1e-11) })
  delta <- setNames(rep(0, CRH_NB), PN$scol); delta[use_s] <- coef(fit_a)[use_s]
  sed   <- setNames(rep(NA_real_, CRH_NB), PN$scol); sed[use_s] <- se(fit_a)[use_s]
  bins <- data.table(bin = 1:CRH_NB, t_center = PN$ctr$t_center,
                     delta = as.numeric(delta), se = as.numeric(sed),
                     is_ref = (1:CRH_NB == PN$ref_bin))
  fp <- bins[is.finite(t_center)]
  sp <- lm(delta ~ ns(t_center, df = 5), data = fp,
           weights = ifelse(is.finite(se), 1 / se, 1 / mean(se, na.rm = TRUE)))
  tg <- seq(min(fp$t_center), max(fp$t_center), length.out = 400)
  curve_a <- data.table(t = tg, delta = as.numeric(predict(sp,
                          newdata = data.frame(t_center = tg))))

  ## (b) MEAN -- twenty quantile bins on the daily mean.
  mv <- d$mean_daily[!is.na(d$mean_daily)]
  brk <- c(-Inf, quantile(mv, seq(1/20, 19/20, 1/20)), Inf)
  d[, mbin := as.integer(cut(mean_daily, brk, labels = 1:20, include.lowest = TRUE))]
  mctr <- quantile(mv, seq(0.025, 0.975, by = 0.05))
  for (i in 1:20) set(d, j = paste0("mb_", i), value = d$mean_daily * (d$mbin == i))
  f_b <- as.formula(paste(lhs, "~",
           paste(c(paste0("mb_", 1:20, collapse = " + "), extra_rhs), collapse = " + "),
           "|", fe_block, "+ fe"))
  fit_b <- tryCatch(feols(f_b, data = d, vcov = conley(CRH_CONLEY_KM), fixef.tol = 1e-11),
                    error = function(e) feols(f_b, data = d, cluster = ~ get(fe_block),
                                              fixef.tol = 1e-11))
  mc <- coef(fit_b)[grep("^mb_", names(coef(fit_b)))]
  ms <- se(fit_b)[grep("^mb_", names(se(fit_b)))]
  o  <- order(as.integer(sub("mb_", "", names(mc))))
  sm <- smooth_mean_drf(as.numeric(mc[o]), as.numeric(ms[o]), as.numeric(mctr))

  cat(sprintf("    (a) DAILY N = %d, %d parts utilisées (bin réf. %d = 0)\n",
              fitstat(fit_a, "n")[[1]], length(use_s), PN$ref_bin))
  cat(sprintf("    (b) MEAN  N = %d, optimum = %s, support [%.1f, %.1f]\n",
              fitstat(fit_b, "n")[[1]],
              ifelse(is.na(sm$Topt), "aucun", sprintf("%.2f", sm$Topt)),
              sm$support[1], sm$support[2]))

  ot <- paste0(tag, suffix)
  fwrite(bins,    file.path(DIR_CRH, sprintf("daily_%s_bins.csv", ot)))
  fwrite(curve_a, file.path(DIR_CRH, sprintf("daily_%s_curve.csv", ot)))
  fwrite(data.table(mu = sm$mu, delta = sm$delta, Lambda = sm$Lambda),
         file.path(DIR_CRH, sprintf("mean_%s_curve.csv", ot)))
  fwrite(data.table(param = c("Topt_mean", "mu_support_min", "mu_support_max", "ref_bin"),
                    value = c(sm$Topt, sm$support[1], sm$support[2], PN$ref_bin)),
         file.path(DIR_CRH, sprintf("mean_%s_coeff.csv", ot)))
  invisible(sm$Topt)
}

run_estimation <- function(keep = NULL, suffix = "") {
  PN <- build_panel(keep)
  cat(sprintf("  panel : %d cellules-périodes | bin réf. %d (%.1f °C)\n",
              nrow(PN$d), PN$ref_bin, PN$ctr[bin == PN$ref_bin]$t_center))
  estimate_channel(PN, "log_amen_HDI", character(0), "grid_2c",
                               "subnal_ctry", "amen", suffix)
  estimate_channel(PN, "log_prod_phi", CHE, "land_type",
                               "subnal_reg", "prod", suffix)
}

run_estimation()
run_estimation(CRH_PERIODS_DEC2, "_dec2")

##################################################################################
# C5. WELL-POPULATED SUPPORT AND OUT-OF-SUPPORT MASS
#
# PRE-REGISTERED rule, fixed before being applied, which is what makes it a
# criterion rather than an adjustment: the hot edge is the last bin whose count
# of cell-periods with positive day mass reaches
# max(1 per cent of the cell-periods, 30).
##################################################################################

cat("\n=== C5 : support et masse hors support ===\n")

support_bounds <- function(sh) {
  ncp <- uniqueN(sh[, .(lon, lat, year)])
  thr <- max(CRH_SUPPORT_MIN_N, ceiling(CRH_SUPPORT_MIN_PCT * ncp))
  cnt <- sh[share > 0, .(n = uniqueN(paste(lon, lat, year))), by = value][order(value)]
  ok <- cnt[n >= thr]
  list(k_min = min(ok$value), k_max = max(ok$value), threshold = thr,
       n_cellperiods = ncp, counts = cnt)
}
sb <- support_bounds(shares)
cat(sprintf("  support bien peuplé = [%.1f, %.1f] °C (seuil %d sur %d cellules-périodes)\n",
            sb$k_min, sb$k_max, sb$threshold, sb$n_cellperiods))
fwrite(data.table(k_min = sb$k_min, k_max = sb$k_max, threshold = sb$threshold,
                  n_cellperiods = sb$n_cellperiods),
       file.path(DIR_CRH, "crh_support_bounds.csv"))

##################################################################################
# C6. PROJECTION ONTO CMIP6, AND THE GATE ON THE NORTHERN REVERSAL
#
#   DAILY : dlog = sum_T delta(T) * [share_2095(T) - share_2015(T)]
#   MEAN  : dlog = Lambda(mu_2095) - Lambda(mu_2015),  mu = sum_T T * share(T)
# Both functions are estimated on THE SAME cells and projected onto THE SAME
# CMIP6 distributions. Only the metric differs.
#
# 2095 and not 2100: the last FULL quinquennial block is 2095-2099, and 2100
# would be a degenerate one-year block. This is the convention of the GE, and the
# text labels that horizon 2100 for readability.
##################################################################################

cat("\n=== C6 : projection CMIP6 ===\n")

cmip <- fread(SERVER_IN[["cmip6"]])[year %in% c(CRH_BASE_YEAR, CRH_FUT_YEAR)]
cmip[, c := count / sum(count), by = .(x, y, year)]
mu_cell <- cmip[, .(mu = sum(value * c)), by = .(x, y, year)]

project_channel <- function(tag, suffix = "", reliable_hi = Inf) {
  ot <- paste0(tag, suffix)
  cv <- fread(file.path(DIR_CRH, sprintf("daily_%s_curve.csv", ot)))
  mc <- fread(file.path(DIR_CRH, sprintf("mean_%s_curve.csv", ot)))
  co <- fread(file.path(DIR_CRH, sprintf("mean_%s_coeff.csv", ot)))
  smin <- co[param == "mu_support_min", value]; smax <- co[param == "mu_support_max", value]

  ## delta(T) truncated and then held flat beyond the reliable support. The
  ## default is Inf, that is, the function as estimated. The clamp is a
  ## robustness option: the projected northern day mass beyond the clamp point is
  ## at most a fraction of a per cent, so clamping is immaterial, and it biases
  ## AGAINST the reversal rather than towards it.
  cvt <- cv[t <= reliable_hi]
  d <- copy(cmip)[, dlt := approx(cvt$t, cvt$delta, xout = value, rule = 2)$y]
  gw <- dcast(d[, .(g = sum(dlt * c)), by = .(x, y, year)], x + y ~ year, value.var = "g")
  setnames(gw, as.character(c(CRH_BASE_YEAR, CRH_FUT_YEAR)), c("g_base", "g_fut"))

  Lam <- function(m) approx(mc$mu, mc$Lambda, xout = pmin(pmax(m, smin), smax), rule = 2)$y
  mw <- dcast(mu_cell, x + y ~ year, value.var = "mu")
  setnames(mw, as.character(c(CRH_BASE_YEAR, CRH_FUT_YEAR)), c("mu_base", "mu_fut"))
  mw <- mw[is.finite(mu_base) & is.finite(mu_fut)]
  mw[, `:=`(mean_base = Lam(mu_base), mean_fut = Lam(mu_fut), dmu = mu_fut - mu_base)]

  w <- merge(gw, mw, by = c("x", "y"))
  w[, `:=`(daily = g_fut - g_base, mean = mean_fut - mean_base, channel = tag)]
  w[, .(x, y, channel, dmu, daily, mean)]
}

aggregate_north <- function(res, out_tag) {
  a <- amen_prod[year == 2005, .(lat, lon, pop = as.numeric(pop_x))]
  res <- merge(res, a, by.x = c("x", "y"), by.y = c("lon", "lat"))
  suppressMessages({ library(sf); library(rnaturalearth) }); sf_use_s2(FALSE)
  pts <- st_as_sf(unique(res[, .(x, y)]), coords = c("x", "y"), crs = 4326)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  ## rnaturalearth codes iso_a3 as "-99" for France and Norway; iso_a3_eh
  ## resolves them. France was assigned geometrically but lost its code, which is
  ## a coding problem and not a resolution one.
  ic <- c("iso_a3_eh", "adm0_a3", "iso_a3", "ISO_A3")
  ic <- ic[ic %in% names(world)][1]
  world <- world[, ic]; names(world)[1] <- "iso3"
  j <- st_join(pts, world, join = st_within); cj <- st_coordinates(j)
  res <- merge(res, data.table(x = cj[, 1], y = cj[, 2], iso3 = j$iso3),
               by = c("x", "y"), all.x = TRUE)

  wm <- function(v, p) sum(v * p, na.rm = TRUE) / sum(p, na.rm = TRUE)
  agg <- res[iso3 %in% CRH_COLD,
             .(dmu = wm(dmu, pop), daily = wm(daily, pop), mean = wm(mean, pop),
               n_cells = .N), by = .(iso3, channel)]
  nor <- res[y >= CRH_NORTH_LAT,
             .(iso3 = "NORTH>=45N", dmu = wm(dmu, pop), daily = wm(daily, pop),
               mean = wm(mean, pop), n_cells = .N), by = channel]
  agg <- rbind(agg, nor[, names(agg), with = FALSE])
  ## The gate is a CLASSIFICATION, not a conclusion: the paper reports an erosion
  ## and an asymmetry of significance, not a reversal.
  agg[, verdict := fifelse(mean > 0 & daily < 0, "REVERSAL",
               fifelse(mean > 0 & daily > 0 & daily < mean, "ATTENUATION",
               fifelse(mean > 0 & daily >= mean, "amplified", "other")))]
  agg[, atten_pct := fifelse(mean > 0, round(100 * (1 - daily / mean), 1), NA_real_)]
  setorder(agg, channel, iso3)
  ## The _reversal name of the CSV files is kept DELIBERATELY: they serve as diff
  ## keys against the archived baselines. The .tex, by contrast, carries a name
  ## that describes a content and does not assert a result.
  fwrite(agg, file.path(DIR_CRH, sprintf("tab_crh_northern_reversal_%s.csv", out_tag)))
  fwrite(res, file.path(DIR_CRH, sprintf("crh_daily_projection_cells%s.csv",
                                         sub("^daily", "", out_tag))))
  print(agg[, .(channel, iso3, dmu = round(dmu, 2), mean = round(mean, 4),
                daily = round(daily, 4), atten_pct, verdict)])
  agg
}

agg_all <- aggregate_north(
  rbindlist(lapply(c("amen", "prod"), project_channel)), "daily")
agg_dec2 <- aggregate_north(
  rbindlist(lapply(c("amen", "prod"), project_channel, suffix = "_dec2")), "daily_dec2")

## Mass projected outside the estimation support, per northern country.
{
  a <- amen_prod[year == 2005, .(lat, lon, pop = as.numeric(pop_x))]
  cf <- merge(cmip[year == CRH_FUT_YEAR], a, by.x = c("x", "y"), by.y = c("lon", "lat"))
  oo <- cf[, .(above = sum(c[value > sb$k_max]), below = sum(c[value < sb$k_min])),
           by = .(x, y, pop)]
  cat(sprintf("  masse projetée au-delà de %.1f °C : %.3f %% (pondérée population)\n",
              sb$k_max, 100 * sum(oo$above * oo$pop) / sum(oo$pop)))
  fwrite(oo, file.path(DIR_CRH, "crh_outofsupport_mass.csv"))
}

##################################################################################
# C7. FIGURES ET TABLE
#
# The table requires crh_inference_deliverable.csv, which carries the 95 per cent
# bootstrap intervals under the published variant. Its producer is SEPARATE from
# this script. When the file is missing the figures are produced and the table is
# not: a mechanism that can do nothing must say whether it acted.
##################################################################################

cat("\n=== C7 : figures et table ===\n")

## ---- figure des reponses : PRODUITE AILLEURS (18/09/2026, lot 10) ----------
## fig_crh_daily_response.pdf (deux panneaux, metrique journaliere seule) est
## remplacee par fig_crh_responses.pdf (quatre panneaux, les deux canaux x les
## deux metriques). Producteur : 01_estimation/01d_fig_crh_responses.R, qui
## ecrit aussi crh_daily_response_peaks.csv.
##
## POURQUOI ELLE SORT D'ICI. Les quatre panneaux ne sont que des LECTURES de CSV
## que ce script vient d'ecrire -- daily_*_bins/curve, *_bins_10y_1h_20b,
## *_coeff_10y_1h_20b_550d, crh_design_optimum_*. Rien n'y est re-estime. Les y
## garder obligeait a rejouer une heure d'inversion pour deplacer un trait, et
## ⚠️ ce rejeu n'est PAS neutre : DIR_CRH ne suit pas JMP_OUTPUT_DIR, donc un
## 01d relance sous redirection reecrit quand meme intermediate_outputs/crh/.
## Le nouveau producteur reproduit crh_daily_response_peaks.csv AU BIT
## (verifie le 18/09/2026), ce qui est la preuve que rien n'a bouge au passage.
## L'ancienne figure est conservee dans le snapshot du lot.

## ---- figures: the two metrics, cell by cell --------------------------------
cells <- fread(file.path(DIR_CRH, "crh_daily_projection_cells.csv"))
base_mu <- mu_cell[year == CRH_BASE_YEAR, .(x, y, mu_base = mu)]
for (tag in c("prod", "amen")) {
  d <- merge(cells[channel == tag], base_mu, by = c("x", "y"))
  f <- file.path(DIR_OUT, sprintf("fig_drf_crh_annual_vs_daily_%s.pdf",
                                  ifelse(tag == "amen", "amenity", tag)))
  pdf(f, width = 8, height = 5); par(mar = c(4.2, 4.2, 1.5, 1))
  plot(d$mu_base, d$mean, pch = 16, cex = 0.3, col = "#00000055",
       xlab = "Baseline mean temperature (°C)",
       ## HORIZON_SPAN_PLAIN, not a year spelled out here: the axis shows the
       ## END of the terminal block while CRH_FUT_YEAR keeps indexing 2095L.
       ylab = sprintf("Projected change in log %s, %s",
                      ifelse(tag == "prod", "productivity", "amenity"),
                      HORIZON_SPAN_PLAIN))
  points(d$mu_base, d$daily, pch = 16, cex = 0.3, col = "#d7302755")
  abline(h = 0, lty = 2, col = "grey40")
  legend("topright", bty = "n", pch = 16, col = c("black", "#d73027"),
         legend = c("Mean metric (CRH)", "Daily distribution"))
  dev.off()
  cat(sprintf("  %s\n", basename(f)))
}

## ---- Table ---------------------------------------------------------------------
## ONE DESTINATION. The deliverable is written by 03_inference/step6_emit.R into
## output/crh_inference/, which is also where emit_values_empirics.py reads it.
## This script used to look in intermediate_outputs/, so the same file had two
## homes and the two consumers could see different vintages, or none.
INF <- file.path(DIR_OUT, "crh_inference", "crh_inference_deliverable.csv")
if (!file.exists(INF)) {
  cat("  ⛔ tab_crh_northern.tex NON PRODUITE : crh_inference_deliverable.csv absent.\n")
  cat("     Les intervalles a 95 % viennent de la chaine 03_inference/, qui tourne\n")
  cat("     APRES ce script : cette table demande un second passage de 01d.\n")
} else {
  inf <- fread(INF)[vcov_variant == "country"]
  ord <- c("RUS", "CAN", "FIN", "DEU", "FRA", "NORTH>=45N")
  nm <- c(RUS = "Russia", CAN = "Canada", FIN = "Finland", DEU = "Germany",
          FRA = "France", `NORTH>=45N` = "North $\\geq$45N")
  ci <- function(ch, iso, obj) {
    r <- inf[channel == ch & unit == iso & object == obj]
    if (nrow(r) != 1L) stop("Lecture d'intervalle impossible : ", ch, "/", iso, "/", obj)
    r
  }
  L <- c("%% GENERATED by 01d_crh_triangulation.R -- DO NOT EDIT BY HAND.",
         "\\begin{tabular}{lccc}", "\\toprule",
         ## « Annual mean », pas « Mean » (18/09/2026). La colonne porte le
         ## bras C4(b), estime sur la moyenne annuelle de la distribution
         ## journaliere. « Mean » tout court laissait le lecteur l'attribuer
         ## a la moyenne decennale du mois d'hiver de C2, qui est un AUTRE
         ## objet -- mediane 0,5 degC contre 10,4, et un optimum a 1,9 degC
         ## contre 20,0. Les deux vivent dans ce script ; un en-tete d'une
         ## seule syllabe ne peut pas les distinguer.
         "Region & Annual mean & Daily & Difference \\\\", "\\midrule")
  for (ch in c("prod", "amen")) {
    L <- c(L, sprintf("\\multicolumn{4}{l}{\\emph{%s}} \\\\",
                      ifelse(ch == "prod", "Productivity", "Amenity")))
    for (iso in ord) {
      r1 <- agg_all[channel == ch & iso3 == iso]; if (!nrow(r1)) next
      cm <- ci(ch, iso, "level_mean_2095"); cd <- ci(ch, iso, "level_daily_2095")
      cf <- ci(ch, iso, "difference_mean_minus_daily")
      ## Provenance guardrail: the inference file must carry THE SAME points as
      ## the projection artefact. Two sources that diverge would put an interval
      ## from one run around a point estimate from another.
      for (p in list(c(cm$point, r1$mean), c(cd$point, r1$daily),
                     c(cf$point, r1$mean - r1$daily)))
        if (abs(p[1] - p[2]) > 1e-6)
          stop(sprintf("Point divergent inférence/projection pour %s/%s : %.6f vs %.6f",
                       ch, iso, p[1], p[2]))
      L <- c(L,
        sprintf("%s & $%+.3f$ & $%+.3f$ & $%+.3f$ \\\\", nm[iso],
                cm$point, cd$point, cf$point),
        sprintf(" & {\\footnotesize $[%+.3f,\\,%+.3f]$} & {\\footnotesize $[%+.3f,\\,%+.3f]$} & {\\footnotesize $[%+.3f,\\,%+.3f]$} \\\\",
                cm$ci95_lo, cm$ci95_hi, cd$ci95_lo, cd$ci95_hi,
                cf$ci95_lo, cf$ci95_hi),
        "\\addlinespace[2pt]")
    }
    L <- c(L, "\\addlinespace")
  }
  writeLines(c(L, "\\bottomrule", "\\end{tabular}"),
             file.path(DIR_OUT, "tab_crh_northern.tex"))
  cat("  tab_crh_northern.tex\n")
}

cat("\n=== 05_crh_triangulation.R terminé ===\n")