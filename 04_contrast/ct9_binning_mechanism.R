##################################################################################
# ct9_binning_mechanism.R — LES INGREDIENTS DU BIAIS DE BINNING, MESURES SANS
#                           ESTIMATION PPML ; et la part intra-cellule du design CRH.
#
# Reads:  intermediate_outputs/climate_long_annual.rds          <- 01a
#         intermediate_outputs/BACI_raw_pairs_annual.rds        <- 01a
#         intermediate_outputs/poids_par_pays.csv               <- poids (00d)
#         intermediate_outputs/crh_inference/step1_frame.rds    <- crh1 (step1_gate.R)
#         output/fe_absorption.csv                              <- 01a  (temoin)
#         output/contribA/binning_variants.csv                  <- ct8  (effectif)
# Writes: output/contribA/binning_mechanism.csv
#         output/contribA/crh_within_cell.csv
#         output/contribA/ct9_controls.csv
#         intermediate_outputs/contribA/ct9_bases_estimation/  (bases, isolees)
# Sourced by: nothing.
#
# ORIGINE. Bac a sable du lot 77 (JMP/_sandbox_binning_20261001/scripts/),
# restreint a ce que le texte cite :
#   v1_mecanisme.R  §1.1 tendance de rechauffement par pays, §1.2 correlation
#                   des tendances de bins avec le climat de base, §1.3 part
#                   intra-paire x produit des regresseurs (r2_pp seul) et son
#                   temoin contre fe_absorption.csv
#   v1bis_cote_y.R  tendance intra-paire x produit de la variable expliquee et
#                   sa correlation avec dTbar
#   v4_crh_intra_cellule.R  part intra-cellule entre periodes, apres les effets
#                   fixes de CRH, par canal et par bras
# Les copies de lib_drf.R du bac a sable sont remplacees par lib_drf.R lui-meme
# (une source de verite) ; les bases qu'il enregistre vont dans un dossier propre.
#
# DEFINITIONS (celles des phrases du papier)
#   part_intra(x)   = 1 - R2 de x sur les effets paire x produit, au niveau ligne
#                     (paire x produit x annee) ; regresseurs de reference : spline4
#                     journalier (4 termes) et spline3 annuel (3 termes) ; min et max.
#   corr tendance y = correlation, entre paires x produits, de la pente OLS de
#                     log(explained) (retrait produit x annee, >= 15 annees > 0 sur
#                     25) avec dTbar = Tbar_orig - Tbar_dest (moyenne 1995-2019).
#   CRH intra       = var(e intra-cellule) / var(e), e = residu du regresseur apres
#                     les effets fixes et controles de step1_gate.R ; mediane par
#                     canal et bras ; mean_daily = la temperature moyenne.
##################################################################################

suppressMessages({ library(data.table); library(splines); library(fixest); library(tools) })
options(width = 200)
t0 <- proc.time(); cat("=== ct9_binning_mechanism.R ===\n")

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
options(jmp.dir_bases = file.path(DIR_TMP_OUT, "contribA", "ct9_bases_estimation"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
setFixest_nthreads(4)

DIR_CONTRIB <- file.path(DIR_OUT, "contribA")
YMIN <- TRADE_YEAR_MIN; YMAX <- TRADE_YEAR_MAX
TOL_R2 <- 1e-8

P <- c(clim  = file.path(DIR_TMP, "climate_long_annual.rds"),
       baci  = file.path(DIR_TMP, "BACI_raw_pairs_annual.rds"),
       pop   = file.path(DIR_TMP, "poids_par_pays.csv"),
       crh   = file.path(DIR_TMP, "crh_inference", "step1_frame.rds"),
       fea   = file.path(DIR_OUT, "fe_absorption.csv"),
       ct8   = file.path(DIR_CONTRIB, "binning_variants.csv"))
cat("\n--- entrees ---\n")
miss <- names(P)[!file.exists(P)]
if (length(miss)) stop("[STOP] entree(s) absente(s) : ", paste(miss, collapse = ", "))
for (n in names(P)) cat(sprintf("  ok  %-5s %s  %s\n", n, md5sum(P[[n]]), P[[n]]))

clim <- as.data.table(readRDS(P[["clim"]]))
clim <- clim[, .(count = sum(count)), by = .(GID_0, year, value)]
if (min(clim$year) != YMIN || max(clim$year) != YMAX)
  stop(sprintf("[STOP] climat %d-%d, fenetre attendue %d-%d", min(clim$year), max(clim$year), YMIN, YMAX))
B0 <- as.data.table(readRDS(P[["baci"]]))
N_PPML <- unique(fread(P[["ct8"]])[variant == "base", n_obs])
if (length(N_PPML) != 1L) stop("[STOP] binning_variants.csv : effectif base non unique")

slope <- function(y, x) { x <- x - mean(x); sum(x * (y - mean(y))) / sum(x^2) }
wcor  <- function(x, y, w) { w <- w / sum(w); mx <- sum(w * x); my <- sum(w * y)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2)) }
OUT <- list()
add <- function(q, v, n, def) OUT[[length(OUT) + 1L]] <<- data.table(quantity = q, value = v, n = n,
                                                                     definition = def)

##################################################################################
# 1. TENDANCE DE RECHAUFFEMENT 1995-2019 PAR PAYS (v1 §1.1)
##################################################################################
tmean <- clim[, .(tm = sum(value * count) / sum(count)), by = .(GID_0, year)]
pp <- fread(P[["pop"]], na.strings = "", select = c("GID_0", "year", "pop"))
pop <- pp[year == 2005 & !is.na(GID_0) & GID_0 != "", .(pop = sum(pop)), by = GID_0]
pays_ech <- sort(unique(c(B0$orig, B0$dest)))
tr <- tmean[, .(pente_dec = 10 * slope(tm, year), tbase = mean(tm)), by = GID_0]
tr <- merge(tr, pop, by = "GID_0", all.x = TRUE)
te <- tr[GID_0 %in% pays_ech & !is.na(pop)]
add("warming_trend_mean_popw", weighted.mean(te$pente_dec, te$pop), nrow(te),
    "population-weighted mean warming trend 1995-2019, degC per decade, sample countries")
add("warming_trend_cor_tbase_popw", wcor(te$pente_dec, te$tbase, te$pop), nrow(te),
    "population-weighted correlation of the 1995-2019 warming trend with baseline temperature")
add("warming_trend_cor_tbase", cor(te$pente_dec, te$tbase), nrow(te),
    "unweighted correlation of the 1995-2019 warming trend with baseline temperature")

##################################################################################
# 2. TENDANCES DES BINS DE 5 degC ET CLIMAT DE BASE (v1 §1.2)
##################################################################################
BRK <- c(-Inf, seq(-10, 35, 5), Inf)
LAB <- c("<-10", sprintf("[%d,%d)", head(seq(-10, 35, 5), -1), seq(-5, 35, 5)), ">=35")
bl <- clim[, .(GID_0, year, value, count)]
bl[, count := 365 * count / sum(count), by = .(GID_0, year)]
bl[, bin := cut(value, BRK, labels = LAB, right = FALSE)]
bins <- dcast(bl, GID_0 + year ~ bin, value.var = "count", fun.aggregate = sum, fill = 0, drop = FALSE)
bins <- bins[!is.na(GID_0)]
bt <- melt(bins, id.vars = c("GID_0", "year"), variable.name = "bin", value.name = "jours")
bt <- bt[, .(pente = 10 * slope(jours, year)), by = .(GID_0, bin)]
bt <- merge(bt, te[, .(GID_0, tbase, pop)], by = "GID_0")
bc <- bt[, .(n = .N, cor_tbase = if (sd(pente) > 0) cor(pente, tbase) else NA_real_), by = bin]
bc[, lo := c(-Inf, seq(-10, 35, 5))[match(bin, LAB)]]
for (i in seq_len(nrow(bc)))
  add(sprintf("bin_trend_cor_tbase|%s", bc$bin[i]), bc$cor_tbase[i], bc$n[i],
      sprintf("correlation across sample countries of the 1995-2019 trend in days in bin %s with baseline temperature (lower edge %s)",
              bc$bin[i], format(bc$lo[i])))

##################################################################################
# 3. PART INTRA-PAIRE x PRODUIT DES REGRESSEURS (v1 §1.3, r2_pp seul) + TEMOIN
##################################################################################
B <- B0[, .(orig, dest, year, product_code_ISICRev3)]
B[, pp := paste(orig, dest, product_code_ISICRev3, sep = "_")]
ajoute_diff <- function(B, tw, tag) {
  cols <- setdiff(names(tw), c("GID_0", "year"))
  o <- copy(tw); setnames(o, c("GID_0", cols), c("orig", paste0(cols, "_o")))
  d <- copy(tw); setnames(d, c("GID_0", cols), c("dest", paste0(cols, "_d")))
  X <- merge(B[, .(orig, dest, year)], o, by = c("orig", "year"), all.x = TRUE, sort = FALSE)
  X <- merge(X, d, by = c("dest", "year"), all.x = TRUE, sort = FALSE)
  X <- unique(X, by = c("orig", "dest", "year"))
  for (k in cols) X[[paste0(tag, "_", k)]] <- X[[paste0(k, "_o")]] - X[[paste0(k, "_d")]]
  X[, c("orig", "dest", "year", paste0(tag, "_", cols)), with = FALSE]
}
termes <- function(tw, pre, tag) {   # garde GID_0, year et les termes ; renomme t1..tk
  cols <- grep(paste0("^", pre, "[1-9]"), names(tw), value = TRUE)
  x <- tw[, c("GID_0", "year", cols), with = FALSE]; setnames(x, cols, paste0("t", seq_along(cols))); x
}
TW <- list(poly3 = termes(make_climate_terms_wide(clim, "polynomial", 3L, 1, 1), "polynomial_"),
           spl4  = termes(make_climate_terms_wide(clim, "spline", 4L, 0.5, 0.5), "spline_"),
           ann3  = termes(make_annual_terms(clim, "spline", 3L, 0.5, 0.5), "annual_spline_"))
for (nm in names(TW))
  B <- merge(B, ajoute_diff(B, TW[[nm]], nm), by = c("orig", "dest", "year"), sort = FALSE)
B <- B[complete.cases(B[, .(poly3_t1, poly3_t2, poly3_t3)])]
cat(sprintf("[3] %d lignes apres filtre climat\n", nrow(B)))
ycols <- grep("^(poly3|spl4|ann3)_t", names(B), value = TRUE)
m <- feols(as.formula(paste0("c(", paste(ycols, collapse = ","), ") ~ 1 | pp")), data = B,
           notes = FALSE, warn = FALSE)
r2pp <- setNames(vapply(seq_along(ycols), function(i) {
  r <- residuals(m[[i]]); 1 - var(r) / var(B[[ycols[i]]]) }, numeric(1)), ycols)
rm(m); gc(verbose = FALSE)

## temoin : r2_pp poly3 == r2_pair_raw de fe_absorption.csv (01a)
fea <- fread(P[["fea"]])
temoin_r2 <- function(r2, ref, tol = TOL_R2) {
  if (length(r2) != 3L || nrow(ref) != 3L) stop("[TEMOIN] poly3 : 3 termes attendus des deux cotes")
  e <- max(abs(r2 - ref$r2_pair_raw))
  if (!is.finite(e) || e > tol) stop(sprintf("[TEMOIN] r2_pp poly3 vs fe_absorption.csv : %.3e > %.0e", e, tol))
  e
}
r2p3 <- unname(r2pp[c("poly3_t1", "poly3_t2", "poly3_t3")])
## connu-faux : perturbe, ampute, NA -- chacun doit etre refuse
kb <- list(plus1e7 = r2p3 + c(0, 1e-7, 0), plus1e2 = r2p3 + c(1e-2, 0, 0),
           ampute = r2p3[-3], na = replace(r2p3, 1, NA_real_))
kbr <- vapply(kb, function(x) tryCatch({ temoin_r2(x, fea); "PASSE" }, error = function(e) "REFUSE"),
              character(1))
print(kbr)
if (any(kbr != "REFUSE")) stop("[STOP] temoin r2 laisse passer : ", paste(names(kbr)[kbr != "REFUSE"], collapse = ", "))
e_r2 <- temoin_r2(r2p3, fea)
cat(sprintf("[temoin] r2_pp poly3 vs fe_absorption.csv : %.3e\n", e_r2))

for (y in ycols) add(sprintf("part_intra|%s", y), 1 - r2pp[[y]], nrow(B),
                     sprintf("1 - R2 of the bilateral differential %s on pair x product effects", y))
ref <- grep("^(spl4|ann3)_t", ycols, value = TRUE)
add("part_intra_ref_min", min(1 - r2pp[ref]), length(ref),
    "min over the reference regressors (daily spline4, annual spline3) of 1 - R2 on pair x product effects")
add("part_intra_ref_max", max(1 - r2pp[ref]), length(ref),
    "max over the reference regressors (daily spline4, annual spline3) of 1 - R2 on pair x product effects")
rm(B); gc(verbose = FALSE)

##################################################################################
# 4. TENDANCE DE LA VARIABLE EXPLIQUEE ET dTbar (v1bis)
##################################################################################
Y <- B0[orig %in% unique(clim$GID_0) & dest %in% unique(clim$GID_0),
        .(orig, dest, year, k = product_code_ISICRev3, explained)]
if (nrow(Y) != N_PPML)
  stop(sprintf("[STOP] echantillon de l'etape 3 : %d lignes, la PPML de ct8 en a %d", nrow(Y), N_PPML))
Y[, pp := paste(orig, dest, k, sep = "_")]
Y[, py := paste(k, year, sep = "_")]
Y[explained > 0, ly := log(explained)]
Y[explained > 0, ly_t := ly - mean(ly), by = py]
SL <- Y[explained > 0, .(n_pos = .N, pente_log = if (.N >= 15) slope(ly_t, year) else NA_real_),
        by = .(pp, orig, dest)]
SW <- Y[, .(poids = mean(explained)), by = pp]
S <- merge(SL, SW, by = "pp")
Pt <- unique(Y[, .(orig, dest, year)])
o <- copy(tmean); setnames(o, c("GID_0", "tm"), c("orig", "tm_o"))
d <- copy(tmean); setnames(d, c("GID_0", "tm"), c("dest", "tm_d"))
Pt <- merge(merge(Pt, o, by = c("orig", "year")), d, by = c("dest", "year"))
PS <- Pt[, .(dTbar = mean(tm_o - tm_d)), by = .(orig, dest)]
S <- merge(S, PS, by = c("orig", "dest"))[!is.na(pente_log)]
add("ytrend_cor_dTbar", cor(S$pente_log, S$dTbar), nrow(S),
    "correlation across pair x products of the within-pair trend of log(explained) (product x year removed) with the baseline temperature differential")
add("ytrend_cor_dTbar_weighted", wcor(S$pente_log, S$dTbar, S$poids), nrow(S),
    "same, weighted by the pair x product mean of explained")
rm(Y, S, Pt); gc(verbose = FALSE)

##################################################################################
# 5. CRH : PART INTRA-CELLULE ENTRE PERIODES APRES LEURS EFFETS FIXES (v4)
##################################################################################
F <- readRDS(P[["crh"]])
dat <- as.data.table(F$dat)
dat[, cell := .GRP, by = .(lat, lon)]
mesure <- function(lhs_ctrl, fe_block, subnal_col, tag) {
  d <- copy(dat)
  d[, subnal := as.integer(factor(get(subnal_col)))]
  d[, fe := (year - 1985) / 5 * max(subnal) + subnal]
  lhs <- if (tag == "amen") "log_amen_HDI" else "log_prod_phi"
  d <- d[is.finite(get(lhs)) & !is.na(get(fe_block))]
  xs <- c(F$use_s, paste0("mb_", 1:20), "mean_daily")
  xs <- xs[vapply(xs, function(v) var(d[[v]]) > 0, logical(1))]
  rhs <- if (length(lhs_ctrl)) paste(lhs_ctrl, collapse = " + ") else "1"
  m <- feols(as.formula(paste0("c(", paste(xs, collapse = ","), ") ~ ", rhs, " | ",
                               fe_block, " + fe")), data = d, notes = FALSE, warn = FALSE,
             fixef.tol = 1e-11)
  rbindlist(lapply(seq_along(xs), function(i) {
    x <- d[[xs[i]]]; e <- residuals(m[[i]])
    xw <- x - ave(x, d$cell); ew <- e - ave(e, d$cell)
    data.table(channel = tag, regressor = xs[i], n = length(x),
               intra_brut = var(xw) / var(x), intra_resid = var(ew) / var(e),
               intra_total = var(ew) / var(x))
  }))
}
R4 <- rbind(mesure(character(0), "grid_2c", "subnal_ctry", "amen"),
            mesure(F$che_vars, "land_type", "subnal_reg", "prod"))
R4[, arm := fifelse(grepl("^s[0-9]+$", regressor), "daily (s*)",
                    fifelse(grepl("^mb_", regressor), "annual (mb_*)", "mean_daily"))]
CW <- R4[, .(n_regressors = .N, n_obs = n[1], intra_brut_med = median(intra_brut),
             intra_resid_min = min(intra_resid), intra_resid_med = median(intra_resid),
             intra_resid_max = max(intra_resid), intra_total_med = median(intra_total)),
         by = .(channel, arm)]
CW[, fe_structure := fifelse(channel == "prod", "land type + sub-region x period, 30 Chebyshev controls",
                             "2x2 grid blocks + sub-country x period")]

##################################################################################
# 6. SORTIE — ecriture atomique
##################################################################################
ecrit <- function(d, nom) {
  p <- file.path(DIR_CONTRIB, nom); tmp <- paste0(p, ".tmp")
  fwrite(d, tmp)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp, p)) stop("[STOP] ecriture impossible de ", p)
  cat(sprintf("Ecrit : %s (%d lignes)\n", p, nrow(d)))
}
MECH <- rbindlist(OUT)
ecrit(MECH, "binning_mechanism.csv")
ecrit(CW,   "crh_within_cell.csv")
ecrit(data.table(control = c("r2_pp poly3 vs fe_absorption.csv r2_pair_raw",
                             "known-bad r2 cases refused",
                             "stage-3 sample rows == ct8 base n_obs"),
                 value = c(e_r2, length(kbr), N_PPML), verdict = "OK"), "ct9_controls.csv")
print(MECH[!grepl("^part_intra\\|", quantity), .(quantity, value = signif(value, 4), n)])
print(CW)
cat(sprintf("\n=== ct9 termine en %.1f min ===\n", (proc.time() - t0)[3] / 60))
