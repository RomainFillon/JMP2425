##################################################################################
# ct7_cold_tail_variants.R — LA QUEUE FROIDE : le contraste sous trois variantes
#                            de projection, SANS AUCUNE RE-ESTIMATION.
#
# Reads:  intermediate_outputs/contribA/tache1_scores.rds   <- ct2
#         output/contribA/spec_curves_{daily,annual}.csv    <- ct1
#         output/contribA/test1_delta_north.csv             <- ct1  (controle)
#         output/contribA/tache1_diff_bras.csv              <- ct2  (controle)
#         intermediate_outputs/CMIP6_tas.csv, crh/amen_prod.csv
#         intermediate_outputs/climate_long_annual.rds
# Writes: output/contribA/ct7_cold_tail_variants.csv
#         output/tab_cold_tail_variants.tex
# Sourced by: nothing.
#
# CE QUE CE SCRIPT FAIT — ET SURTOUT CE QU'IL NE FAIT PAS
# --------------------------------------------------------
# §2.3 declare la queue froide de la reponse journaliere NON IDENTIFIEE, et §3.3
# attribue a la plage sous 0 degC l'essentiel du contraste nordique. L'objection
# se resume en une phrase : le resultat repose sur la partie de la fonction que le
# papier declare lui-meme fragile. On y repond en RE-PROJETANT, pas en re-estimant.
#
# Delta_north est LINEAIRE en beta :
#     daily   Delta = (sum_b Delta m_b . b(T_b))' beta  =  w' beta
#     annual  Delta = (sum_c omega_c [b(mu_fut) - b(mu_base)])' beta
# Une variante de projection ne touche QUE w. beta, la vcov jointe, les clusters
# et l'echantillon sont ceux du run publie, lus a l'archive. Aucun PPML n'est
# relance : c'est exactement « evaluation des fonctions deja ajustees contre les
# changements d'exposition archives ».
#
# LES TROIS VARIANTES
# -------------------
#   v1  cold_flat_0      reponse journaliere tenue plate SOUS 0 degC a f(0)
#   v2  cold_flat_5      idem SOUS 5 degC a f(5)  -- la borne que le texte nomme
#   v3  common_support   projection restreinte au support COMMUN aux deux bras
# v1 et v2 ne touchent que le bras journalier : c'est sa queue froide qui est en
# cause. v3 contraint les DEUX bras, puisque c'est un support partage.
# Aplatir f sous s revient a evaluer la base en pmax(T, s) : f(max(T,s)) est
# encore lineaire en beta, donc la variance analytique s'ecrit a l'identique.
#
# LE CONTROLE QUI AUTORISE TOUT LE RESTE  (G1, G2, G3 — BLOQUANTS)
# -----------------------------------------------------------------
# Un script qui reconstruit w et V a partir d'une archive doit d'abord REPRODUIRE
# l'etat publie, sinon ses variantes mesurent sa propre derive.
#   G1  le w que je reconstruis == le w archive par ct2            (1e-10)
#   G2  le Delta_north que j'en tire == celui de test1_delta_north.csv
#   G3  la difference et son ecart-type == ceux de tache1_diff_bras.csv
# Les trois portent sur la variante BASELINE, celle qui doit retomber sur le
# publie. Aucune variante n'est ecrite si l'une des trois echoue.
#
# ⚠ beta n'est PAS archive. Il est RECUPERE de la courbe publiee par moindres
# carres sur la base d'estimation (la courbe est exactement base %*% beta - max,
# donc le systeme est consistant et la recuperation exacte). G2 est precisement
# ce qui verifie que cette recuperation est bonne : si beta etait mal recupere,
# le Delta baseline ne retomberait pas sur le publie.
##################################################################################

suppressMessages({
  library(data.table); library(splines); library(tools)
})
options(width = 200)
t0 <- proc.time(); cat("=== ct7_cold_tail_variants.R ===\n")

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "04_contrast", "lib_contrast.R"))

DIR_CONTRIB <- file.path(DIR_OUT, "contribA")
BASE_YEAR <- 2015L; FUT_YEAR <- 2095L; NORTH_LAT <- 45; POP_YEAR <- 2005L
CHAN <- "prod"
TOL_W <- 1e-10; TOL_D <- 1e-8

P <- c(scores  = file.path(DIR_TMP, "contribA", "tache1_scores.rds"),
       cur_d   = file.path(DIR_CONTRIB, "spec_curves_daily.csv"),
       cur_a   = file.path(DIR_CONTRIB, "spec_curves_annual.csv"),
       t1      = file.path(DIR_CONTRIB, "test1_delta_north.csv"),
       diff    = file.path(DIR_CONTRIB, "tache1_diff_bras.csv"),
       cmip    = file.path(DIR_TMP, "CMIP6_tas.csv"),
       pop     = file.path(DIR_CRH, "amen_prod.csv"),
       clim    = file.path(DIR_TMP, "climate_long_annual.rds"))
cat("\n--- entrees ---\n")
miss <- names(P)[!file.exists(P)]
if (length(miss)) stop("[STOP] entree(s) absente(s) : ", paste(miss, collapse = ", "))
for (n in names(P)) cat(sprintf("  ok  %-8s %s\n", n, as.character(tools::md5sum(P[n]))))

STORE <- readRDS(P["scores"])
CUR_D <- fread(P["cur_d"]); CUR_A <- fread(P["cur_a"])
T1    <- fread(P["t1"])[channel == CHAN & status == "OK"]
DIFF  <- fread(P["diff"])[channel == CHAN & status == "OK"]
clim  <- as.data.table(readRDS(P["clim"]))

##################################################################################
# 1. GEOGRAPHIE NORDIQUE — recopie fidele de ct1 §3 (memes fichiers, memes poids)
##################################################################################
cat("\n--- geographie nordique ---\n")
cm <- fread(P["cmip"])[year %in% c(BASE_YEAR, FUT_YEAR) & y >= NORTH_LAT]
cm[, c := count / sum(count), by = .(x, y, year)]
pop <- fread(P["pop"])[year == POP_YEAR, .(lat, lon, pop = as.numeric(pop_x))]
cm <- merge(cm, pop, by.x = c("x", "y"), by.y = c("lon", "lat"), all = FALSE)
ok2 <- cm[, .(ny = uniqueN(year)), by = .(x, y)][ny == 2L, .(x, y)]
cm <- merge(cm, ok2, by = c("x", "y"))
ncell <- nrow(ok2)
wcell <- unique(cm[, .(x, y, pop)]); wcell[, omega := pop / sum(pop)]
stopifnot(abs(sum(wcell$omega) - 1) < 1e-12)
cm <- merge(cm, wcell[, .(x, y, omega)], by = c("x", "y"))
cm[, sgn := fifelse(year == FUT_YEAR, 1, -1)]
dmass <- cm[, .(m = sum(omega * sgn * c)), by = value][order(value)]
share_by_T <- cm[, .(s = sum(omega * c)), by = .(year, value)]
mu  <- cm[, .(mu = sum(value * c)), by = .(x, y, year)]
muw <- dcast(mu, x + y ~ year, value.var = "mu")
setnames(muw, c(as.character(BASE_YEAR), as.character(FUT_YEAR)), c("mu_base", "mu_fut"))
muw <- merge(muw, wcell[, .(x, y, omega)], by = c("x", "y"))
muw <- muw[is.finite(mu_base) & is.finite(mu_fut)]
cat(sprintf("  cellules >= %gN : %d | somme des masses signees : %.2e\n",
            NORTH_LAT, ncell, sum(dmass$m)))
cat(sprintf("  masse-jour nordique sous 0 degC (futur) : %.3f %%\n",
            100 * sum(share_by_T[year == FUT_YEAR & value < 0, s])))
cat(sprintf("  masse-jour nordique sous 5 degC (futur) : %.3f %%\n",
            100 * sum(share_by_T[year == FUT_YEAR & value < 5, s])))
rm(cm); gc(verbose = FALSE)

##################################################################################
# 2. OUTILS
##################################################################################
## beta recupere de la courbe publiee : la courbe est base %*% beta - max, donc
## une regression sur [base, 1] rend beta exactement.
##
## 🚨 LA COLONNE temperature_C EST ARRONDIE, ET S'EN SERVIR COUTE 5,7e-4.
## ct1 ecrit `round(gec$x, 1)` alors que gec$x vaut seq(lo, hi, by = 0.1). Au
## bras JOURNALIER lo est ENTIER, la grille tombe juste et l'arrondi est exact
## (G2 mesure 2,6e-17). Au bras ANNUEL lo est le min des moyennes annuelles,
## donc reel : chaque abscisse est deplacee de jusqu'a 0,05 degC, et beta
## recupere sur ces abscisses-la manque le Delta publie de 5,7e-4 -- assez pour
## faire tomber G2, et c'est bien ce qu'on lui demande.
## On reconstruit donc la grille EXACTE a partir des bornes de la base, et on
## s'en sert pour predire ; la colonne arrondie ne sert plus qu'a verifier
## l'appariement ligne a ligne.
recover_beta <- function(curve, B) {
  x <- seq(B$lo, B$hi, by = 0.1)
  if (length(x) != nrow(curve))
    stop(sprintf("appariement courbe/grille impossible : %d abscisses reconstruites, %d lignes",
                 length(x), nrow(curve)))
  if (max(abs(round(x, 1) - curve$temperature_C)) > 1e-9)
    stop("appariement courbe/grille : les abscisses arrondies ne coincident pas")
  X <- predict(B$basis, x)
  ok <- is.finite(curve$y) & apply(is.finite(X), 1, all)
  fit <- lm.fit(cbind(X[ok, , drop = FALSE], 1), curve$y[ok])
  b <- unname(fit$coefficients[seq_len(ncol(X))])
  list(beta = b, resid = max(abs(fit$residuals)))
}

## V propre et bloc croise, EXACTEMENT comme ct2 (l.641, l.647-650).
Vown <- function(bread, meat, ssc) ssc * (bread %*% meat %*% bread)

## Delta, se, et l'intervalle, pour un couple (w_daily, w_annual) donne.
contrast_for <- function(st, wd, wa, bd, ba) {
  kd <- length(wd); ka <- length(wa)
  Vd <- Vown(st$bread_daily,  st$meat_daily,  st$ssc[["daily"]])
  Va <- Vown(st$bread_annual, st$meat_annual, st$ssc[["annual"]])
  Vx <- st$ssc[["cross"]] * (st$bread_daily %*% st$meat_cross %*% st$bread_annual)
  Vj <- matrix(0, kd + ka, kd + ka)
  Vj[1:kd, 1:kd] <- Vd; Vj[(kd+1):(kd+ka), (kd+1):(kd+ka)] <- Va
  Vj[1:kd, (kd+1):(kd+ka)] <- Vx; Vj[(kd+1):(kd+ka), 1:kd] <- t(Vx)

  u  <- c(-wd, wa)                      # D = annual - daily
  ud <- c(wd, numeric(ka)); ua <- c(numeric(kd), wa)
  var_raw <- as.numeric(t(u) %*% Vj %*% u)
  if (is.finite(var_raw) && var_raw >= 0) {
    Vpub <- Vj; var_pub <- var_raw; psd <- "none"; ntr <- 0L
  } else {
    ei <- eigen((Vj + t(Vj))/2, symmetric = TRUE)
    ntr <- sum(ei$values < 0)
    Vpub <- ei$vectors %*% diag(pmax(ei$values, 0), nrow(Vj)) %*% t(ei$vectors)
    var_pub <- as.numeric(t(u) %*% Vpub %*% u); psd <- "eigen_trunc"
  }
  v_d <- as.numeric(t(ud) %*% Vpub %*% ud)
  v_a <- as.numeric(t(ua) %*% Vpub %*% ua)
  cov_da <- as.numeric(t(ud) %*% Vpub %*% ua)

  dd <- sum(wd * bd); da <- sum(wa * ba); D <- da - dd
  se <- sqrt(var_pub); se_n <- sqrt(v_d + v_a)
  list(delta_daily = dd, delta_annual = da, diff = D,
       se_daily = sqrt(v_d), se_annual = sqrt(v_a),
       se_diff = se, se_diff_naive = se_n,
       rho = cov_da / sqrt(v_d * v_a),
       d_lo = dd - 1.96*sqrt(v_d), d_hi = dd + 1.96*sqrt(v_d),
       a_lo = da - 1.96*sqrt(v_a), a_hi = da + 1.96*sqrt(v_a),
       diff_lo = D - 1.96*se, diff_hi = D + 1.96*se,
       diff_lo_naive = D - 1.96*se_n, diff_hi_naive = D + 1.96*se_n,
       daily_excl = (dd - 1.96*sqrt(v_d)) * (dd + 1.96*sqrt(v_d)) > 0,
       annual_excl = (da - 1.96*sqrt(v_a)) * (da + 1.96*sqrt(v_a)) > 0,
       diff_excl = (D - 1.96*se) * (D + 1.96*se) > 0,
       diff_excl_naive = (D - 1.96*se_n) * (D + 1.96*se_n) > 0,
       psd_fix = psd, n_eigen_trunc = ntr)
}

##################################################################################
# 3. BOUCLE : 15 specifications x 4 variantes, sans une seule re-estimation
##################################################################################
VARIANTS <- c("baseline", "cold_flat_0", "cold_flat_5", "common_support")
OUT <- list(); G <- list()

for (cfg in GRID_FORMS) {
  spec <- cfg$spec; deg <- cfg$degree
  for (wc in GRID_WINSOR) {
    wl <- wc$lower; wu <- wc$upper
    key <- sprintf("%s%d_wl%gwu%g", spec, deg, wl, wu)
    st  <- STORE[[paste(CHAN, key, sep = "|")]]
    if (is.null(st)) stop("[STOP] archive sans entree pour ", key)
    cat(sprintf("\n[%s] ", key)); flush.console()

    Bd <- daily_basis(clim, spec, deg, wl, wu)
    Ba <- annual_basis(clim, spec, deg, wl, wu, mean_climate)

    rb <- recover_beta(CUR_D[model_name == key], Bd)
    ra <- recover_beta(CUR_A[model_name == key], Ba)
    bd <- rb$beta; ba <- ra$beta
    nd <- length(bd); na <- length(ba)
    kd <- length(st$w_daily); ka <- length(st$w_annual)

    ## w pleins : les coefficients climatiques d'abord, les controles a zero.
    ## L'ordre est celui de cc_daily / cc_annual, verifie par G1.
    pad <- function(wc_, k) { z <- numeric(k); z[seq_along(wc_)] <- wc_; z }

    ## supports
    lo_d <- Bd$lo; hi_d <- Bd$hi; lo_a <- Ba$lo; hi_a <- Ba$hi
    cs_lo <- max(lo_d, lo_a); cs_hi <- min(hi_d, hi_a)
    if (cs_lo >= cs_hi) stop("[STOP] support commun vide pour ", key)

    for (v in VARIANTS) {
      ## --- bras journalier : floor sur la temperature avant la base ---------
      fl <- switch(v, baseline = lo_d, cold_flat_0 = max(lo_d, 0),
                   cold_flat_5 = max(lo_d, 5), common_support = cs_lo)
      ch <- if (v == "common_support") cs_hi else hi_d
      Tcl <- pmin(pmax(dmass$value, fl), ch)
      wd_clim <- as.numeric(t(predict(Bd$basis, Tcl)) %*% dmass$m)

      ## --- bras annuel : inchange sauf sous support commun ------------------
      if (v == "common_support") {
        Xf <- predict(Ba$basis, pmin(pmax(muw$mu_fut,  cs_lo), cs_hi))
        Xb <- predict(Ba$basis, pmin(pmax(muw$mu_base, cs_lo), cs_hi))
      } else {
        Xf <- predict(Ba$basis, pmin(pmax(muw$mu_fut,  lo_a), hi_a))
        Xb <- predict(Ba$basis, pmin(pmax(muw$mu_base, lo_a), hi_a))
      }
      wa_clim <- as.numeric(t(Xf - Xb) %*% muw$omega)

      wd <- pad(wd_clim, kd); wa <- pad(wa_clim, ka)

      ## ---- G1 : sur la variante baseline, w doit retomber sur l'archive ----
      if (v == "baseline") {
        g1d <- max(abs(wd - st$w_daily)); g1a <- max(abs(wa - st$w_annual))
        if (!is.finite(g1d) || g1d > TOL_W || !is.finite(g1a) || g1a > TOL_W)
          stop(sprintf(paste0("G1 VIOLEE (%s) : le w reconstruit ne reproduit pas ",
                              "celui de l'archive (daily %.3e, annual %.3e). Les ",
                              "variantes mesureraient ma derive, pas la queue froide."),
                       key, g1d, g1a))
      }

      r <- contrast_for(st, wd, wa, bd, ba)

      ## ---- G2 / G3 : la baseline doit retomber sur le publie ---------------
      if (v == "baseline") {
        pd <- T1[arm == "daily"  & spec_label == key]
        pa <- T1[arm == "annual" & spec_label == key]
        pD <- DIFF[spec_label == key]
        g2d <- abs(r$delta_daily  - pd$delta_north)
        g2a <- abs(r$delta_annual - pa$delta_north)
        g3D <- abs(r$diff - pD$diff); g3s <- abs(r$se_diff - pD$se_diff)
        if (max(g2d, g2a) > TOL_D)
          stop(sprintf("G2 VIOLEE (%s) : Delta baseline != publie (d %.3e, a %.3e)",
                       key, g2d, g2a))
        if (max(g3D, g3s) > TOL_D)
          stop(sprintf("G3 VIOLEE (%s) : diff/se baseline != publie (%.3e, %.3e)",
                       key, g3D, g3s))
        G[[key]] <- data.table(spec_label = key, g1_w_daily = g1d, g1_w_annual = g1a,
                               g2_delta_daily = g2d, g2_delta_annual = g2a,
                               g3_diff = g3D, g3_se_diff = g3s,
                               beta_resid_daily = rb$resid, beta_resid_annual = ra$resid)
        cat("G1/G2/G3 ok ")
      }

      ## masse-jour nordique effectivement aplatie par cette variante
      flat <- sum(share_by_T[year == FUT_YEAR & value < fl, s])

      OUT[[length(OUT)+1L]] <- data.table(
        channel = CHAN, spec_label = key, form = spec, df = deg,
        winsor_lo = wl, winsor_hi = wu, variant = v,
        floor_daily = fl, ceil_daily = ch,
        support_lo_daily = lo_d, support_hi_daily = hi_d,
        support_lo_annual = lo_a, support_hi_annual = hi_a,
        common_lo = cs_lo, common_hi = cs_hi,
        daymass_flattened_fut = flat,
        delta_daily = r$delta_daily, se_daily = r$se_daily,
        daily_lo = r$d_lo, daily_hi = r$d_hi, daily_excl = r$daily_excl,
        delta_annual = r$delta_annual, se_annual = r$se_annual,
        annual_lo = r$a_lo, annual_hi = r$a_hi, annual_excl = r$annual_excl,
        diff = r$diff, se_diff = r$se_diff,
        diff_lo = r$diff_lo, diff_hi = r$diff_hi, diff_excl = r$diff_excl,
        se_diff_naive = r$se_diff_naive,
        diff_lo_naive = r$diff_lo_naive, diff_hi_naive = r$diff_hi_naive,
        diff_excl_naive = r$diff_excl_naive,
        rho = r$rho, psd_fix = r$psd_fix, n_eigen_trunc = r$n_eigen_trunc,
        n_cells_north = ncell)
    }
    cat(sprintf("| base D=%+.4f  v1 D=%+.4f  v2 D=%+.4f  v3 D=%+.4f",
                OUT[[length(OUT)-3L]]$diff, OUT[[length(OUT)-2L]]$diff,
                OUT[[length(OUT)-1L]]$diff, OUT[[length(OUT)]]$diff))
  }
}

R <- rbindlist(OUT)
GG <- rbindlist(G)
stopifnot(nrow(R) == GRID_NSPEC * length(VARIANTS))

##################################################################################
# 4. SORTIE
##################################################################################
fwrite(R,  file.path(DIR_CONTRIB, "ct7_cold_tail_variants.csv"))
fwrite(GG, file.path(DIR_CONTRIB, "ct7_controls.csv"))
cat(sprintf("\n\nEcrit: %s  (%d lignes)\n",
            file.path(DIR_CONTRIB, "ct7_cold_tail_variants.csv"), nrow(R)))

cat("\n================ SYNTHESE PAR VARIANTE ================\n")
for (v in VARIANTS) {
  s <- R[variant == v]
  cat(sprintf("\n%-15s  masse-jour aplatie (futur) : %.3f %%\n",
              v, 100 * mean(s$daymass_flattened_fut)))
  cat(sprintf("   nord daily   : [%+.4f, %+.4f]  mediane %+.4f  | IC excluant 0 : %d/%d\n",
              min(s$delta_daily), max(s$delta_daily), median(s$delta_daily),
              sum(s$daily_excl), nrow(s)))
  cat(sprintf("   nord annual  : [%+.4f, %+.4f]  mediane %+.4f  | IC excluant 0 : %d/%d\n",
              min(s$delta_annual), max(s$delta_annual), median(s$delta_annual),
              sum(s$annual_excl), nrow(s)))
  cat(sprintf("   difference   : [%+.4f, %+.4f]  mediane %+.4f\n",
              min(s$diff), max(s$diff), median(s$diff)))
  cat(sprintf("   diff excluant 0 : %d/%d avec covariance | %d/%d sans (naif)\n",
              sum(s$diff_excl), nrow(s), sum(s$diff_excl_naive), nrow(s)))
  cat(sprintf("   signe de la difference constant : %s\n",
              if (all(s$diff > 0) || all(s$diff < 0)) "OUI" else "NON"))
}
cat(sprintf("\n=== ct7 termine en %.0f s ===\n", (proc.time()-t0)["elapsed"]))
