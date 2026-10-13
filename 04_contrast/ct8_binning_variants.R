##################################################################################
# ct8_binning_variants.R — LE BIAIS DE BINNING : la PPML du canal productivite
#                          sous quatre structures d'effets fixes.
#
# Reads:  intermediate_outputs/BACI_raw_pairs_annual.rds     <- 01a
#         intermediate_outputs/climate_long_annual.rds       <- 01a
#         intermediate_outputs/CMIP6_tas.csv                 [SERVEUR]
#         intermediate_outputs/crh/amen_prod.csv             <- 01d
#         output/contribA/test1_delta_north.csv              <- ct1  (temoin)
#         output/contribA/tache1_diff_bras.csv               <- ct2  (temoin)
# Writes: output/contribA/binning_variants.csv
#         output/contribA/binning_contrasts.csv
#         output/contribA/ct8_controls.csv
#         intermediate_outputs/contribA/ct8_bases_estimation/  (bases, isolees)
# Sourced by: nothing.
#
# ORIGINE. Transpose du bac a sable du lot 77 (JMP/_sandbox_binning_20261001/
# scripts/v2_short_panel.R, rapport_lot77_v2.txt). Les seuls changements sont la
# racine du paquet, les chemins de sortie, l'ordre (base d'abord, temoin
# AVANT les autres variantes) et le temoin rendu bloquant a 1e-12 avec son test
# connu-faux. L'arithmetique des ajustements et des contrastes est inchangee.
#
# L'OBJECTION VISEE. Jones, Moscona, Olken, von Dessauer (2026, NBER w34671) :
# dans un panel a effets fixes de lieu, les comptes de jours extremes portent
# des tendances induites par le rechauffement, correlees au climat de base ; si
# la variable expliquee a aussi une tendance correlee au climat, la regression
# fabrique un U. Remede : effets fixes lieu x periode courte.
#
# LES QUATRE VARIANTES (bras journalier spline4, annuel spline4 et spline3,
# winsorisation 0,5/0,5 — les specifications de reference de ct2 et du GE) :
#   base     produit x annee                            (la PPML publiee de 01a)
#   between  produit x annee, Delta_k remplace par sa moyenne 1995-2019 par
#            paire : toute variation temporelle intra-paire du regresseur, donc
#            le canal du biais, est retiree sans toucher aux effets fixes
#   pp       produit x annee + paire x produit           (le panel de Jones et al.)
#   pp5      produit x annee + paire x produit x 5 ans   (leur remede « short panel »)
# La PPML publiee N'A PAS d'effet fixe paire x produit (01a, fe_str =
# "fe_prod_year") : pp et pp5 sont des contre-factuels, pas des robustesses.
#
# LE TEMOIN QUI AUTORISE LE RESTE (BLOQUANT)
# ------------------------------------------
# base doit retrouver test1_delta_north.csv (Delta_north et se, 3 ajustements)
# et tache1_diff_bras.csv (D et se_diff, spline4) a 1e-12. Il est evalue JUSTE
# APRES base, avant les trois autres variantes : un ecart arrete l'etape sans
# rien ecrire. Mesure au bac a sable : 3,1e-17 et 4,2e-17.
# Connu-faux (execute a chaque course, avant le temoin reel) : la fonction de
# temoin doit REFUSER une base perturbee de 1e-11 et de 1e-3 sur chacune des
# quatre quantites, une base amputee d'une ligne et une base portant un NA.
# Si elle en laisse passer un seul, l'etape s'arrete : un temoin qui ne refuse
# pas n'en est pas un.
#
# DUREE. Bac a sable, 01/10/2026 : 5 min 28 s pour les quatre variantes
# (base 49 s, between 46 s, pp 77 s, pp5 98 s d'ajustements, plus chargement et
# contrastes). Memoire : pic 6-8 Go.
##################################################################################

suppressMessages({ library(data.table); library(splines); library(fixest); library(tools) })
options(width = 200)
t0 <- proc.time(); cat("=== ct8_binning_variants.R ===\n")

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
## Les bases d'estimation que make_climate_terms_wide() enregistre vont dans un
## dossier PROPRE a ct8 : ecrire dans celui de ct1 re-ecrirait ses fichiers sous
## les memes cles (md5 neuf, REGISTRE allonge) pour un usage qui n'est pas le sien.
options(jmp.dir_bases = file.path(DIR_TMP_OUT, "contribA", "ct8_bases_estimation"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "04_contrast", "lib_contrast.R"))

DIR_CONTRIB <- file.path(DIR_OUT, "contribA")
dir.create(DIR_CONTRIB, showWarnings = FALSE, recursive = TRUE)
TOL_TEMOIN <- 1e-12

VARIANTES <- c("base", "between", "pp", "pp5")          # base EN PREMIER (temoin)
FE_LABEL <- c(base = "product x year", between = "product x year, pair-mean differentials",
              pp = "product x year + pair x product",
              pp5 = "product x year + pair x product x 5 years")
## LES TROIS AJUSTEMENTS, SANS UNE ENTREE DE GRILLE EN DUR (controle 9 de
## check_macros). daily_ref et annual_ref : GRID_REF_SPEC dans les deux bras, la
## reference P4 de ct2 et du §2.5. annual_ge : la specification annuelle que le GE
## retient, c'est-a-dire l'argmax du poids d'Akaike du bras annuel dans la sortie
## de ct1 (poids degenere a 1 : on s'arrete s'il ne l'est plus, le choix serait
## alors ambigu). Au 01/10/2026 : spline4/spline4/spline3, winsorisation 0,5.
T1_ALL <- fread(file.path(DIR_OUT, "contribA", "test1_delta_north.csv"))[channel == "prod" & status == "OK"]
grid_entry <- function(lab) {
  z <- Filter(function(s) identical(s$spec_label, lab), GRID_SPECS)
  if (length(z) != 1L) stop("[STOP] specification absente de GRID_SPECS : ", lab)
  z[[1]]
}
.ann <- T1_ALL[arm == "annual"][order(-aic_weight)]
if (!nrow(.ann) || .ann$aic_weight[1] < 0.999)
  stop("[STOP] bras annuel : poids d'Akaike maximal ", if (nrow(.ann)) .ann$aic_weight[1] else NA,
       " < 0,999 -- la specification annuelle du GE n'est plus determinee sans ambiguite.")
fit_of <- function(arm, lab) { z <- grid_entry(lab)
  list(arm = arm, spec = z$spec, deg = z$degree, wl = z$winsor_lo, wu = z$winsor_hi) }
FITS <- list(daily_ref  = fit_of("daily",  GRID_REF_SPEC),
             annual_ref = fit_of("annual", GRID_REF_SPEC),
             annual_ge  = fit_of("annual", .ann$spec_label[1]))
CONTRASTES <- list(              # D = annuel - journalier
  ref_ct2 = c(annual = "annual_ref", daily = "daily_ref"),   # §2.5, reference P4 de ct2
  ref_ge  = c(annual = "annual_ge",  daily = "daily_ref"))   # specs de reference des bras GE
cat(sprintf("[ajustements] daily_ref = annual_ref = %s ; annual_ge = %s\n",
            GRID_REF_SPEC, .ann$spec_label[1]))
T_PTS <- c(-10, -5, 0, 5, 25, 30); T_REF <- 15
t_tag <- function(T) if (T < 0) paste0("m", abs(T)) else as.character(T)

##################################################################################
# 0. ENTREES
##################################################################################
P <- c(baci = file.path(DIR_TMP, "BACI_raw_pairs_annual.rds"),
       clim = file.path(DIR_TMP, "climate_long_annual.rds"),
       cmip = file.path(DIR_TMP, "CMIP6_tas.csv"),
       pop  = file.path(DIR_CRH, "amen_prod.csv"),
       t1   = file.path(DIR_CONTRIB, "test1_delta_north.csv"),
       diff = file.path(DIR_CONTRIB, "tache1_diff_bras.csv"))
cat("\n--- entrees ---\n")
miss <- names(P)[!file.exists(P)]
if (length(miss)) stop("[STOP] entree(s) absente(s) : ", paste(miss, collapse = ", "))
MD5 <- vapply(P, function(p) unname(md5sum(p)), character(1))
for (n in names(P)) cat(sprintf("  ok  %-5s %s  %s\n", n, MD5[[n]], P[[n]]))
## empreinte des QUATRE intrants d'estimation (pas des fichiers temoins)
EMPR <- empreinte(setNames(MD5[c("baci", "clim", "cmip", "pop")],
                           basename(P[c("baci", "clim", "cmip", "pop")])))
cat(sprintf("  empreinte des intrants d'estimation : %s\n", EMPR))

dat  <- readRDS(P[["baci"]]); setDT(dat)
clim <- readRDS(P[["clim"]]); setDT(clim)
cat(sprintf("  BACI : %d lignes | climat : %d lignes\n", nrow(dat), nrow(clim)))
dat[, fe_pp  := paste(orig, dest, product_code_ISICRev3, sep = "_")]
dat[, fe_pp5 := paste(fe_pp, (year - TRADE_YEAR_MIN) %/% 5L, sep = "_")]

##################################################################################
# 1. GEOGRAPHIE NORDIQUE — recopie de ct1 (memes fichiers, memes poids)
##################################################################################
BASE_YEAR <- 2015L; FUT_YEAR <- 2095L; NORTH_LAT <- 45; POP_YEAR <- 2005L
cm <- fread(P[["cmip"]])[year %in% c(BASE_YEAR, FUT_YEAR) & y >= NORTH_LAT]
cm[, c := count / sum(count), by = .(x, y, year)]
pop <- fread(P[["pop"]])[year == POP_YEAR, .(lat, lon, pop = as.numeric(pop_x))]
cm <- merge(cm, pop, by.x = c("x", "y"), by.y = c("lon", "lat"), all = FALSE)
ok2 <- cm[, .(ny = uniqueN(year)), by = .(x, y)][ny == 2L, .(x, y)]
cm <- merge(cm, ok2, by = c("x", "y")); ncell <- nrow(ok2)
wcell <- unique(cm[, .(x, y, pop)]); wcell[, omega := pop / sum(pop)]
cm <- merge(cm, wcell[, .(x, y, omega)], by = c("x", "y"))
cm[, sgn := fifelse(year == FUT_YEAR, 1, -1)]
dmass <- cm[, .(m = sum(omega * sgn * c)), by = value][order(value)]
mu  <- cm[, .(mu = sum(value * c)), by = .(x, y, year)]
muw <- dcast(mu, x + y ~ year, value.var = "mu")
setnames(muw, c(as.character(BASE_YEAR), as.character(FUT_YEAR)), c("mu_base", "mu_fut"))
muw <- merge(muw, wcell[, .(x, y, omega)], by = c("x", "y"))
muw <- muw[is.finite(mu_base) & is.finite(mu_fut)]
cat(sprintf("  Nord : %d cellules\n", ncell)); rm(cm); gc(verbose = FALSE)

##################################################################################
# 2. MACHINERIE CGM — recopie de ct2 (bloc croise entre bras, B1 bloquante)
##################################################################################
rowsel <- function(mod, n_full) {
  r <- seq_len(n_full)
  if (!is.null(mod$obs_selection) && length(mod$obs_selection))
    for (s in mod$obs_selection) r <- r[s]
  r
}
meat_cross <- function(S1, S2, g_or, g_de, g_pair) {
  ag <- function(S, g) rowsum(S, group = g, reorder = TRUE)
  crossprod(ag(S1, g_or), ag(S2, g_or)) + crossprod(ag(S1, g_de), ag(S2, g_de)) -
    crossprod(ag(S1, g_pair), ag(S2, g_pair))
}
ssc_recover <- function(V_raw, V_fixest) {
  c_hat <- sum(V_fixest * V_raw) / sum(V_raw * V_raw)
  list(c = c_hat, rel_norm = max(abs(V_fixest - c_hat * V_raw)) / max(abs(V_fixest)))
}

##################################################################################
# 3. UN AJUSTEMENT
##################################################################################
fe_of <- function(v) switch(v, base = "fe_prod_year", between = "fe_prod_year",
                            pp = c("fe_prod_year", "fe_pp"), pp5 = c("fe_prod_year", "fe_pp5"))

ajuste <- function(v, f) {
  is_ann <- f$arm == "annual"
  tk <- proc.time()
  dp <- prepare_diffs(dat, clim, f$spec, f$deg, f$wl, f$wu, is_ann, scale = SCALE_PROD)
  dp <- stage3_climate_filter(dp, f$spec, f$deg, is_ann)
  dcols <- grep(if (is_ann) "^diff_annual_spline_[1-9]" else "^diff_spline_[1-9]",
                names(dp), value = TRUE)
  if (v == "between") {             # Delta_k -> moyenne 1995-2019 de la paire
    u <- unique(dp[, c("orig", "dest", "year", dcols), with = FALSE])
    u <- u[, lapply(.SD, mean), by = .(orig, dest), .SDcols = dcols]
    dp[, (dcols) := NULL]; dp <- merge(dp, u, by = c("orig", "dest"), sort = FALSE)
  }
  invisible(capture.output(
    r <- run_stage3(dp, f$spec, f$deg, is_ann, controls = TRADE_STAGE3_CONTROLS,
                    fe_str = fe_of(v))))
  mod <- r$model
  rws <- rowsel(mod, nrow(dp))
  keys <- list(or = as.character(dp$orig[rws]), de = as.character(dp$dest[rws]))
  n_in <- nrow(dp); rm(dp); gc(verbose = FALSE)
  cc <- intersect(dcols, names(coef(mod)))
  ## base de Delta_north et des courbes : celle de ct1
  if (!is_ann) {
    B <- daily_basis(clim, f$spec, f$deg, f$wl, f$wu)
    w <- as.numeric(t(predict(B$basis, pmin(pmax(dmass$value, B$lo), B$hi))) %*% dmass$m)
  } else {
    B <- annual_basis(clim, f$spec, f$deg, f$wl, f$wu, mean_climate)
    Xf <- predict(B$basis, pmin(pmax(muw$mu_fut,  B$lo), B$hi))
    Xb <- predict(B$basis, pmin(pmax(muw$mu_base, B$lo), B$hi))
    w  <- as.numeric(t(Xf - Xb) %*% muw$omega)
  }
  ## un coefficient absorbe par les effets fixes ne laisse pas w' beta defini
  if (length(cc) != length(w))
    return(list(ok = FALSE, msg = sprintf("%d/%d coefficients climatiques estimes",
                                          length(cc), length(w)), secs = (proc.time() - tk)[3]))
  wfull <- setNames(numeric(length(coef(mod))), names(coef(mod))); wfull[cc] <- w
  bc <- function(T) predict(B$basis, pmin(pmax(T, B$lo), B$hi))
  list(ok = TRUE, mod = mod, keys = keys, rows = rws, cc = cc, w = w, wfull = wfull,
       B = B, bc = bc, n_in = n_in, wald_p = r$wald_joint$p,
       secs = as.numeric((proc.time() - tk)[3]))
}

##################################################################################
# 4. UNE VARIANTE : trois ajustements, courbes, deux contrastes
##################################################################################
une_variante <- function(v) {
  cat(sprintf("\n######## variante %s (effets fixes : %s) ########\n", v,
              paste(fe_of(v), collapse = " + ")))
  OF <- list(); FT <- list(); OC <- list()
  for (nm in names(FITS)) {
    f <- FITS[[nm]]
    lab <- sprintf("%s%d_wl%gwu%g", f$spec, f$deg, f$wl, f$wu)
    cat(sprintf("  [%s] %s %s ... ", v, f$arm, lab)); flush.console()
    g <- tryCatch(ajuste(v, f), error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    ## un ajustement manquant laisserait une ligne de table vide sous un nom
    ## plein : l'etape s'arrete.
    if (!isTRUE(g$ok)) stop(sprintf("[STOP] ajustement %s / %s / %s : %s", v, f$arm, lab, g$msg))
    mod <- g$mod; V <- vcov(mod)[g$cc, g$cc, drop = FALSE]; beta <- coef(mod)[g$cc]
    dn <- as.numeric(crossprod(g$w, beta)); vv <- delta_var(g$w, V); se <- sqrt(vv$var_psd)
    cat(sprintf("N=%d  Dn=%+.5f se=%.5f  (%.0f s)\n", nobs(mod), dn, se, g$secs))
    row <- data.table(variant = v, fe_structure = FE_LABEL[[v]], fit = nm, arm = f$arm,
      spec_label = lab, status = "OK",
      delta_north = dn, se = se, ci_lo = dn - 1.96 * se, ci_hi = dn + 1.96 * se,
      excludes_zero = (dn - 1.96 * se) * (dn + 1.96 * se) > 0,
      psd_fix = vv$fix, n_obs = nobs(mod), n_in = g$n_in, wald_p = g$wald_p,
      n_fe = paste(mod$fixef_sizes, collapse = "/"), secs = g$secs,
      support_lo = g$B$lo, support_hi = g$B$hi, T_ref = T_REF)
    ## f(T) - f(T_REF), intervalle sur la vcov propre
    for (T in T_PTS) {
      wt <- as.numeric(g$bc(T) - g$bc(T_REF)); vt <- delta_var(wt, V)
      est <- as.numeric(crossprod(wt, beta)); s <- sqrt(vt$var_psd); tg <- t_tag(T)
      row[, (paste0("f_", tg)) := est]
      row[, (paste0("f_", tg, "_se")) := s]
      row[, (paste0("f_", tg, "_lo")) := est - 1.96 * s]
      row[, (paste0("f_", tg, "_hi")) := est + 1.96 * s]
      row[, (paste0("f_", tg, "_out_of_support")) := (T < g$B$lo | T > g$B$hi)]
    }
    row[, inputs_empreinte := EMPR]
    OF[[nm]] <- row; FT[[nm]] <- g; rm(mod); gc(verbose = FALSE)
  }
  ## ---- contrastes : bloc croise CGM (ct2) ----
  for (cn in names(CONTRASTES)) {
    pa <- FT[[CONTRASTES[[cn]][["annual"]]]]; pd <- FT[[CONTRASTES[[cn]][["daily"]]]]
    OWN <- list(); b1 <- 0
    for (a in c("d", "a")) {
      g <- if (a == "d") pd else pa
      S <- as.matrix(g$mod$scores); gp <- paste(g$keys$or, g$keys$de, sep = "\r")
      Br <- solve(g$mod$hessian); Vraw <- Br %*% meat_cross(S, S, g$keys$or, g$keys$de, gp) %*% Br
      rc <- ssc_recover(Vraw, vcov(g$mod)); b1 <- max(b1, rc$rel_norm)
      if (!is.finite(rc$rel_norm) || rc$rel_norm >= 1e-8)
        stop(sprintf("[STOP] B1 violee (%s/%s/%s) : %.3e", v, cn, a, rc$rel_norm))
      OWN[[a]] <- list(V = rc$c * Vraw, bread = Br, ssc = rc$c)
    }
    common <- intersect(pd$rows, pa$rows)
    pos_d <- match(common, pd$rows); pos_a <- match(common, pa$rows)
    Sd <- as.matrix(pd$mod$scores)[pos_d, , drop = FALSE]
    Sa <- as.matrix(pa$mod$scores)[pos_a, , drop = FALSE]
    g_or <- pd$keys$or[pos_d]; g_de <- pd$keys$de[pos_d]
    stopifnot(identical(g_or, pa$keys$or[pos_a]), identical(g_de, pa$keys$de[pos_a]))
    ssc_x <- sqrt(OWN$d$ssc * OWN$a$ssc)
    Vx <- ssc_x * (OWN$d$bread %*% meat_cross(Sd, Sa, g_or, g_de, paste(g_or, g_de, sep = "\r")) %*%
                     OWN$a$bread)
    rm(Sd, Sa)
    kd <- length(coef(pd$mod)); ka <- length(coef(pa$mod))
    Vj <- matrix(0, kd + ka, kd + ka)
    Vj[1:kd, 1:kd] <- OWN$d$V; Vj[(kd + 1):(kd + ka), (kd + 1):(kd + ka)] <- OWN$a$V
    Vj[1:kd, (kd + 1):(kd + ka)] <- Vx; Vj[(kd + 1):(kd + ka), 1:kd] <- t(Vx)
    u <- c(-pd$wfull, pa$wfull); var_raw <- as.numeric(t(u) %*% Vj %*% u)
    if (is.finite(var_raw) && var_raw >= 0) { Vp <- Vj; fix <- "none" } else {
      ei <- eigen((Vj + t(Vj)) / 2, symmetric = TRUE)
      Vp <- ei$vectors %*% diag(pmax(ei$values, 0), nrow(Vj)) %*% t(ei$vectors); fix <- "eigen_trunc" }
    ud <- c(pd$wfull, numeric(ka)); ua <- c(numeric(kd), pa$wfull)
    vd <- as.numeric(t(ud) %*% Vp %*% ud); va <- as.numeric(t(ua) %*% Vp %*% ua)
    cda <- as.numeric(t(ud) %*% Vp %*% ua)
    dd <- as.numeric(crossprod(pd$w, coef(pd$mod)[pd$cc]))
    da <- as.numeric(crossprod(pa$w, coef(pa$mod)[pa$cc]))
    D <- da - dd; seD <- sqrt(as.numeric(t(u) %*% Vp %*% u))
    OC[[cn]] <- data.table(variant = v, fe_structure = FE_LABEL[[v]], contrast = cn,
      annual = CONTRASTES[[cn]][["annual"]], daily = CONTRASTES[[cn]][["daily"]],
      delta_daily = dd, delta_annual = da, diff = D, se_diff = seD,
      ci_lo = D - 1.96 * seD, ci_hi = D + 1.96 * seD,
      excludes_zero = (D - 1.96 * seD) * (D + 1.96 * seD) > 0,
      rho = cda / sqrt(vd * va), se_diff_naive = sqrt(vd + va), psd_fix = fix,
      n_common = length(common), b1_max_relnorm = b1, inputs_empreinte = EMPR)
    cat(sprintf("  contraste %s : D=%+.5f se=%.5f rho=%.3f\n", cn, D, seD, cda / sqrt(vd * va)))
  }
  rm(FT); gc(verbose = FALSE)
  list(fits = rbindlist(OF), contrasts = rbindlist(OC))
}

##################################################################################
# 5. LE TEMOIN, ET SON TEST CONNU-FAUX
##################################################################################
T1   <- fread(P[["t1"]])[channel == "prod"]
DIFF <- fread(P[["diff"]])[channel == "prod"]

## Refuse (stop) si base ne retombe pas sur le publie. Rend la table des ecarts.
temoin <- function(fits, cons, t1, t2, tol = TOL_TEMOIN) {
  fb <- fits[variant == "base"]
  if (nrow(fb) != length(FITS)) stop(sprintf("[TEMOIN] %d ajustements base, %d attendus",
                                             nrow(fb), length(FITS)))
  m1 <- merge(fb[, .(arm, spec_label, delta_north, se)],
              t1[, .(arm, spec_label, dn_pub = delta_north, se_pub = se)],
              by = c("arm", "spec_label"))
  if (nrow(m1) != length(FITS)) stop(sprintf("[TEMOIN] %d ajustements apparies a test1, %d attendus",
                                             nrow(m1), length(FITS)))
  e1 <- max(abs(c(m1$delta_north - m1$dn_pub, m1$se - m1$se_pub)))
  r2 <- cons[variant == "base" & contrast == "ref_ct2"]
  p2 <- t2[spec_label == GRID_REF_SPEC]
  if (nrow(r2) != 1L || nrow(p2) != 1L)
    stop(sprintf("[TEMOIN] contraste ref_ct2 : %d ligne(s) ici, %d dans tache1", nrow(r2), nrow(p2)))
  e2 <- max(abs(c(r2$diff - p2$diff, r2$se_diff - p2$se_diff)))
  if (!is.finite(e1) || e1 > tol)
    stop(sprintf("[TEMOIN] Delta_north/se base vs test1_delta_north.csv : ecart %.3e > %.0e", e1, tol))
  if (!is.finite(e2) || e2 > tol)
    stop(sprintf("[TEMOIN] D/se_diff base vs tache1_diff_bras.csv : ecart %.3e > %.0e", e2, tol))
  data.table(control = c("base Delta_north and se vs test1_delta_north.csv (3 fits)",
                         "base D and se_diff (spline4) vs tache1_diff_bras.csv"),
             n = c(nrow(m1), 1L), max_abs_gap = c(e1, e2), tolerance = tol, verdict = "OK")
}

## Connu-faux : construit une base qui RETOMBE EXACTEMENT sur le publie (lue dans
## les fichiers temoins), verifie qu'elle passe, puis la casse de chaque facon.
connu_faux <- function(t1, t2) {
  fb <- rbindlist(lapply(names(FITS), function(nm) {
    f <- FITS[[nm]]; lab <- sprintf("%s%d_wl%gwu%g", f$spec, f$deg, f$wl, f$wu)
    p <- t1[arm == f$arm & spec_label == lab]
    if (nrow(p) != 1L) stop("[CONNU-FAUX] test1 sans ligne unique pour ", f$arm, " ", lab)
    data.table(variant = "base", arm = f$arm, spec_label = lab, delta_north = p$delta_north, se = p$se)
  }))
  p2 <- t2[spec_label == GRID_REF_SPEC]
  cb <- data.table(variant = "base", contrast = "ref_ct2", diff = p2$diff, se_diff = p2$se_diff)
  invisible(temoin(fb, cb, t1, t2))                       # le vrai doit passer
  cas <- list()
  for (eps in c(1e-11, 1e-3)) {
    for (q in c("delta_north", "se")) {
      x <- copy(fb); x[2, (q) := get(q) + eps]; cas[[sprintf("%s+%g", q, eps)]] <- list(x, cb) }
    for (q in c("diff", "se_diff")) {
      x <- copy(cb); x[1, (q) := get(q) + eps]; cas[[sprintf("%s+%g", q, eps)]] <- list(fb, x) }
  }
  cas[["ligne base retiree"]] <- list(fb[-1], cb)
  cas[["contraste retire"]]   <- list(fb, cb[0])
  x <- copy(fb); x[1, delta_north := NA_real_]; cas[["NA"]] <- list(x, cb)
  res <- vapply(names(cas), function(k) {
    r <- tryCatch({ temoin(cas[[k]][[1]], cas[[k]][[2]], t1, t2); "PASSE" },
                  error = function(e) "REFUSE")
    cat(sprintf("  connu-faux %-22s %s\n", k, r)); r }, character(1))
  if (any(res != "REFUSE"))
    stop("[STOP] le temoin laisse passer : ", paste(names(res)[res != "REFUSE"], collapse = ", "))
  length(res)
}
cat("\n--- connu-faux du temoin ---\n")
N_KB <- connu_faux(T1, DIFF)
cat(sprintf("  %d cas casses, %d refuses\n", N_KB, N_KB))

##################################################################################
# 6. BOUCLE : base, temoin, puis les trois autres
##################################################################################
RES <- list()
RES[["base"]] <- une_variante("base")
cat("\n--- temoin (bloquant) ---\n")
CTRL <- temoin(RES$base$fits, RES$base$contrasts, T1, DIFF)
print(CTRL)
for (v in setdiff(VARIANTES, "base")) RES[[v]] <- une_variante(v)

FITS_OUT <- rbindlist(lapply(RES, `[[`, "fits"), fill = TRUE)
CONS_OUT <- rbindlist(lapply(RES, `[[`, "contrasts"), fill = TRUE)
stopifnot(nrow(FITS_OUT) == length(VARIANTES) * length(FITS),
          nrow(CONS_OUT) == length(VARIANTES) * length(CONTRASTES),
          all(FITS_OUT$status == "OK"))
CTRL <- rbind(CTRL, data.table(control = "known-bad cases refused by the witness",
                               n = N_KB, max_abs_gap = NA_real_, tolerance = TOL_TEMOIN,
                               verdict = "OK"))

##################################################################################
# 7. SORTIE — ecriture atomique (tmp puis rename)
##################################################################################
ecrit <- function(d, nom) {
  p <- file.path(DIR_CONTRIB, nom); tmp <- paste0(p, ".tmp")
  fwrite(d, tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp, p)) stop("[STOP] ecriture impossible de ", p)
  cat(sprintf("Ecrit : %s (%d lignes)\n", p, nrow(d)))
}
ecrit(FITS_OUT, "binning_variants.csv")
ecrit(CONS_OUT, "binning_contrasts.csv")
ecrit(CTRL,     "ct8_controls.csv")

cat("\n=============== Delta_north ===============\n")
print(FITS_OUT[, .(variant, fit, dn = round(delta_north, 5), se = round(se, 5),
                   cold = round(f_m10, 4), cold_lo = round(f_m10_lo, 4), cold_hi = round(f_m10_hi, 4),
                   n_obs, n_fe, secs = round(secs))])
cat("\n=============== contrastes ===============\n")
print(CONS_OUT[, .(variant, contrast, D = round(diff, 5), se = round(se_diff, 5),
                   ci_lo = round(ci_lo, 5), ci_hi = round(ci_hi, 5), excl0 = excludes_zero,
                   rho = round(rho, 3))])
cat(sprintf("\n=== ct8 termine en %.1f min ===\n", (proc.time() - t0)[3] / 60))
