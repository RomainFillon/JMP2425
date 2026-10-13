##################################################################################
# 01c_robustness.R -- robustness of the damage functions.
#
# FIVE STEPS, in the order of the argument in the paper:
#   R1  selection on unobservables        -> tab_oster_delta.tex
#   R2  instrumenting with latitude       -> tab_iv_latitude.tex
#   R3  asymmetric fixed effects          -> tab_oh_asym_fe.tex,
#                                            tab_oh_hotslope_range.tex
#   R4  controls and income               -> tab_oh_fundamental_controls.tex,
#                                            tab_oh_income_badcontrol.tex
#   R5  exporter cost correction          -> tab_robustness_astier_{prod,amenity}
#
# R3 to R5 all answer the SAME criticism, on the exporter cost term, by three
# different routes: closing the channel with fixed effects (R3), controlling for
# the fundamentals (R4), and purging income of its climate component (R5). Three
# routes agreeing is the argument; any one of them alone would not be.
#
# Reads:  the panels and caches produced by 01a and 01b,
#         output/oh_hdi_country_year.csv when present,
#         intermediate_outputs/{mrt_structural,gdp_purged}.rds for R5
# Writes: tab_oster_delta.tex, oster_delta_full.csv, tab_iv_latitude.tex,
#         tab_oh_*.tex, oh_hotslope_schemes.csv, and the R5 tables
# Needs:  01a and 01b must have run.
#
# R5 additionally needs mrt_structural.rds and gdp_purged.rds, which come from a
# script outside the runbook. When they are absent R5 is SKIPPED and says so.
##################################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(fixest); library(splines)
  library(ggplot2)
})

## Single-threaded estimation (05/10/2026, lot 85j; see rapport_lot85i). With the
## default 4 threads the two-FE schemes (exporter, importer) carried a 1e-13
## run-to-run noise in beta_hat and vcov, which flipped a sign inside
## MASS::mvrnorm and gave the second-seed MC band two exact states. Global
## options: they also govern the feglm() calls made from lib_drf.R.
fixest::setFixest_nthreads(1)
data.table::setDTthreads(1)

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
source(file.path(.root, "01_estimation", "lib_drf.R"))


BACI    <- readRDS(file.path(DIR_TMP, "BACI_raw_pairs_annual.rds"))
climate <- readRDS(file.path(DIR_TMP, "climate_long_annual.rds"))
setDT(BACI); setDT(climate)

GRID <- expand.grid(i = seq_along(SPEC_CONFIGS), j = seq_along(WINSOR_CONFIGS))
grid_iter <- function(f) {
  res <- list()
  for (n in seq_len(nrow(GRID))) {
    cfg <- SPEC_CONFIGS[[GRID$i[n]]]; w <- WINSOR_CONFIGS[[GRID$j[n]]]
    key <- sprintf("%s%d_wl%swu%s", cfg$spec, cfg$degree, w$lower, w$upper)
    res[[key]] <- f(cfg$spec, cfg$degree, w$lower, w$upper, key)
  }
  res
}

##################################################################################
# R1. SELECTION ON UNOBSERVABLES -- Oster (2019) delta*
#
#   delta* = [b_tilde (R_tilde - R_dot)] / [(b_dot - b_tilde)(R_max - R_tilde)]
# with b_dot without the gravity controls, b_tilde with them, and
#   R = 1 - dev(model) / dev(FE only)   McFadden pseudo R-squared, for PPML
#   R_max = min(1, 1.3 * R_tilde)       as recommended by Oster
#
# delta* varies as 1/(b_dot - b_tilde), so a coefficient that BARELY MOVES when
# the controls are added produces an arbitrarily large |delta*| whose sign is not
# informative. Extreme values therefore indicate MAXIMUM STABILITY, not
# fragility. The table note must say so: without it a reader takes a large
# negative number for a disaster.
##################################################################################

cat("\n=== R1 : Oster delta* ===\n")

run_oster <- function(dt, spec, deg) {
  cols <- paste0("diff_", spec, "_", seq_len(deg))
  cols <- cols[cols %in% names(dt)]
  if (!length(cols)) return(NULL)
  fit <- function(rhs) tryCatch(feglm(
    as.formula(paste0("explained ~ ", if (length(rhs)) paste(rhs, collapse = " + ") else "1")),
    data = as.data.frame(dt), family = poisson(),
    fixef = "fe_prod_year", cluster = c("orig", "dest")), error = function(e) NULL)

  m_null <- fit(character(0)); if (is.null(m_null)) return(NULL)
  m_dot  <- fit(cols);         if (is.null(m_dot))  return(NULL)
  m_til  <- fit(c(cols, TRADE_STAGE3_CONTROLS)); if (is.null(m_til)) return(NULL)

  R_dot <- 1 - m_dot$deviance / m_null$deviance
  R_til <- 1 - m_til$deviance / m_null$deviance
  R_max <- min(1, 1.3 * R_til)
  bd <- coef(m_dot)[cols[1]]; bt <- coef(m_til)[cols[1]]
  delta <- if (is.na(bd) || is.na(bt)) NA
           else if (abs(bd - bt) < 1e-10) Inf
           else (bt * (R_til - R_dot)) / ((bd - bt) * (R_max - R_til))
  list(delta = delta, b_dot = bd, b_tilde = bt,
       R_dot = R_dot, R_tilde = R_til, R_max = R_max,
       p_dot = wald(m_dot, cols)$p, p_tilde = wald(m_til, cols)$p,
       spec = spec, degree = deg)
}

res_oster <- grid_iter(function(spec, deg, wl, wu, key) {
  cat(sprintf("  [%s]\n", key))
  dt <- prepare_diffs(BACI, climate, spec, deg, wl, wu, scale = SCALE_PROD)
  ## The differentials themselves are winsorised, at the bounds of the basis.
  for (dc in paste0("diff_", spec, "_", seq_len(deg))) {
    q <- quantile(dt[[dc]], c(wl / 100, 1 - wu / 100), na.rm = TRUE)
    dt <- dt[get(dc) >= q[1] & get(dc) <= q[2]]
  }
  r <- run_oster(dt[is.finite(GDP_diff)], spec, deg)
  if (!is.null(r)) { r$wval_lower <- wl; r$wval_upper <- wu }
  r
})

oster_csv <- rbindlist(lapply(names(res_oster), function(k) {
  r <- res_oster[[k]]; if (is.null(r)) return(NULL)
  data.table(tag = k, spec = r$spec, degree = r$degree,
             wval_lower = r$wval_lower, wval_upper = r$wval_upper,
             delta_term1 = r$delta, b_dot = r$b_dot, b_tilde = r$b_tilde,
             R_dot = r$R_dot, R_tilde = r$R_tilde, R_max = r$R_max,
             p_dot = r$p_dot, p_tilde = r$p_tilde)
}))
fwrite(oster_csv, file.path(DIR_OUT, "oster_delta_full.csv"))

{
  lab <- c(polynomial = "Polynomial", spline = "Natural spline")
  body <- sprintf("%s %d & $[%.1f\\%%,%.1f\\%%]$ & %s & %.4f & %.4f & %.4f & %.4f & %.3f & %.3f \\\\",
    lab[oster_csv$spec], oster_csv$degree, oster_csv$wval_lower, oster_csv$wval_upper,
    ifelse(is.infinite(oster_csv$delta_term1), "$\\infty$",
      ifelse(is.na(oster_csv$delta_term1), "---",
             sprintf("%.2f", oster_csv$delta_term1))),
    oster_csv$b_dot, oster_csv$b_tilde, oster_csv$R_dot, oster_csv$R_tilde,
    oster_csv$p_dot, oster_csv$p_tilde)
  writeLines(c(
    tex_header("01c_robustness.R (R1)", "oster_delta_full.csv"),
    "\\begin{table}[H]", "\\centering\\small",
    "\\caption{Oster (2019) selection-on-unobservables robustness for the",
    "  productivity dose-response function. $\\delta^*$: degree of selection on",
    "  unobservables relative to gravity controls needed to explain away the",
    "  estimated climate effect ($\\beta^*=0$), computed with",
    "  $R_{\\text{max}}=\\min(1,\\,1.3\\times\\tilde{R})$.}",
    "\\label{tab:oster_delta}",
    "\\begin{tabular}{lcccccccc}", "\\toprule",
    paste0("Specification & Winsorization & $\\delta^*$ & $\\dot\\beta_1$ & ",
           "$\\tilde\\beta_1$ & $\\dot{R}$ & $\\tilde{R}$ & $p$ (no ctrl.) & ",
           "$p$ (ctrl.) \\\\"),
    "\\midrule", body, "\\bottomrule",
    tex_note(9, c(
      "\\textit{Note.} $\\delta^*$ is computed for the first basis coefficient",
      ## LOT 74 (01/10/2026, P3) : « $\delta^*>2$ is the conventional threshold »
      ## contredisait le corps (§2.3, « against a threshold of one ») ; un seul seuil.
      "(dominant linear term). $\\delta^*>1$ means unobservables would have to be",
      "more important than the observable gravity controls to explain away the",
      "estimated effect; one is the threshold, as in the text.",
      "Because $\\delta^*$ scales with $1/(\\dot\\beta_1-\\tilde\\beta_1)$, a",
      "coefficient that barely moves when controls are added produces an",
      "arbitrarily large $|\\delta^*|$ of uninformative sign: extreme values in",
      "this column indicate maximal coefficient stability, not fragility.",
      "Within-pseudo-$R^2$ is measured relative to a PPML model with",
      "product$\\times$year fixed effects only.")),
    "\\end{tabular}", "\\end{table}"),
    file.path(DIR_OUT, "tab_oster_delta.tex"))
  cat("[write] tab_oster_delta.tex + oster_delta_full.csv\n")
}

##################################################################################
# R2. INSTRUMENTING WITH LATITUDE
#
# Instrument: Z_rs = |lat_orig| - |lat_dest|, on national centroids. Absolute
# latitude predicts long-run mean temperature strongly and is excluded from the
# trade equation conditional on the controls.
#
# EXPECTED SIGN: Z < 0 means the origin is closer to the equator, hence warmer,
# hence a positive climate differential. The first stage must therefore return a
# NEGATIVE coefficient on Z. A positive sign would invalidate the reading of the
# instrument, not merely its strength.
#
# Estimator: control function (Rivers and Vuong 1988) for IV-PPML.
#   stage 1: linear regression of the climate differential on Z and the
#            controls, with product x year fixed effects;
#   stage 2: PPML with the first-stage residual as a regressor.
#   The endogeneity test is the nullity of the coefficient on that residual.
#
# RESERVATION TO KEEP IN THE TEXT: the exclusion restriction is debatable.
# Latitude correlates with institutions and development, and the gravity
# controls do not absorb all of it. This column bounds a sensitivity; it does not
# provide an alternative identification.
##################################################################################

cat("\n=== R2 : IV latitude ===\n")

## Latitude is NOT in the CEPII gravity database: none of its columns carries it.
## It is computed here from the centroids of the national polygons.
suppressMessages({ library(sf); library(rnaturalearth) }); sf_use_s2(FALSE)
world <- ne_countries(scale = "medium", returnclass = "sf")
## The iso_a3 column codes France and Norway as "-99". Using that column and
## then filtering out "-99" silently DROPS France, one of the five countries the
## paper discusses. iso_a3_eh resolves them.
ic <- c("iso_a3_eh", "adm0_a3", "iso_a3")
ic <- ic[ic %in% names(world)][1]
cent <- suppressWarnings(st_coordinates(st_centroid(st_geometry(world))))
cent_dt <- data.table(iso3 = world[[ic]],
                      lat = cent[, 2])[iso3 != "-99" & !is.na(iso3)]
cat(sprintf("  centroïdes : %d pays (colonne %s) | FRA présent : %s\n",
            nrow(cent_dt), ic, "FRA" %in% cent_dt$iso3))
cent_dt <- data.table(iso3 = world[[ic]],
                      lat = cent[, 2])[iso3 != "-99" & !is.na(iso3)]
## rnaturalearth carries several polygons for some countries, so st_centroid
## returns one per polygon and the iso3 codes come out DUPLICATED. The join would
## then multiply rows without anything saying the sample had changed. One
## latitude per country is kept.
n_dup <- nrow(cent_dt) - uniqueN(cent_dt$iso3)
if (n_dup > 0)
  cat(sprintf("  [!] %d centroide(s) en double, moyenne par pays\n", n_dup))
cent_dt <- cent_dt[, .(lat = mean(lat)), by = iso3]
cat(sprintf("  centroides : %d pays (colonne %s) | FRA present : %s\n",
            nrow(cent_dt), ic, "FRA" %in% cent_dt$iso3))
            
run_iv <- function(dt, spec, deg) {
  d1 <- paste0("diff_", spec, "_1")
  dt <- merge(dt, cent_dt[, .(orig = iso3, la_o = lat)], by = "orig", all.x = TRUE)
  dt <- merge(dt, cent_dt[, .(dest = iso3, la_d = lat)], by = "dest", all.x = TRUE)
  dt[, Z := abs(la_o) - abs(la_d)]
  dt <- dt[is.finite(Z) & is.finite(get(d1))]
  if (nrow(dt) < 1000) return(NULL)

  fs <- feols(as.formula(paste0(d1, " ~ Z + ",
                paste(TRADE_STAGE3_CONTROLS, collapse = " + "), " | fe_prod_year")),
              data = dt, cluster = ~ orig + dest)
  dt[, cf_resid := as.numeric(resid(fs))]
  ## The first-stage F is on the INSTRUMENT ALONE: what is measured is its
  ## strength, not that of the control block.
  F_first <- tryCatch(wald(fs, "Z")$stat, error = function(e) NA_real_)
  beta_Z  <- tryCatch(coef(fs)[["Z"]],    error = function(e) NA_real_)

  cols <- paste0("diff_", spec, "_", seq_len(deg))
  ss <- feglm(as.formula(paste0("explained ~ ",
                paste(c(cols, TRADE_STAGE3_CONTROLS, "cf_resid"), collapse = " + "))),
              data = as.data.frame(dt), family = poisson(),
              fixef = "fe_prod_year", cluster = ~ orig + dest)
  ct <- summary(ss)$coeftable
  in_mod <- cols[cols %in% names(coef(ss))]
  list(F_first = F_first, beta_Z = beta_Z,
       beta_iv = coef(ss)[cols[1]],
       p_joint = if (length(in_mod)) wald(ss, in_mod)$p else NA_real_,
       p_endog = if ("cf_resid" %in% rownames(ct)) ct["cf_resid", 4] else NA_real_,
       nobs = nobs(ss), spec = spec, degree = deg)
}
res_iv <- grid_iter(function(spec, deg, wl, wu, key) {
  cat(sprintf("  [%s]\n", key))
  r <- run_iv(prepare_diffs(BACI, climate, spec, deg, wl, wu,
                            scale = SCALE_PROD), spec, deg)
  if (!is.null(r)) { r$wval_lower <- wl; r$wval_upper <- wu }
  r
})

{
  iv <- rbindlist(lapply(res_iv, function(r) if (is.null(r)) NULL else as.data.table(r)))
  lab <- c(polynomial = "Polynomial", spline = "Natural spline")
  body <- sprintf("%s %d & $[%s\\%%, %s\\%%]$ & %.1f & %.4f & %s & %s \\\\",
    lab[iv$spec], iv$degree, iv$wval_lower, iv$wval_upper, iv$F_first, iv$beta_iv,
    ifelse(iv$p_joint < 0.001, "$<0.001$", sprintf("%.3f", iv$p_joint)),
    ifelse(is.na(iv$p_endog), "---", sprintf("%.3f", iv$p_endog)))
  writeLines(c(
    tex_header("01c_robustness.R (R2)"),
    "\\begin{table}[H]", "\\centering\\small",
    "\\caption{Instrumenting the bilateral temperature differential with the",
    "  difference in absolute latitude.}",
    "\\label{tab:iv_latitude}",
    "\\begin{tabular}{llccc c}", "\\toprule",
    paste0("Functional form & Winsorization & First-stage $F$ & $\\hat\\beta_1$ ",
           "(IV) & $p$ (joint) & $p$ (endogeneity) \\\\"),
    "\\midrule", body, "\\bottomrule",
    tex_note(6, c(
      "\\textit{Note.} Control-function IV-PPML \\citep{rivers1988two}. First stage:",
      "linear regression of the bilateral climate differential on",
      "$Z_{rs}=|\\text{lat}_r|-|\\text{lat}_s|$ and the gravity controls, with",
      "product$\\times$year fixed effects. Second stage: PPML including the",
      "first-stage residual; the endogeneity test is the nullity of its",
      "coefficient. Standard errors clustered by origin and destination.",
      sprintf("$N$ ranges from %s to %s across cells.",
              formatC(min(iv$nobs), format = "d", big.mark = "{,}"),
              formatC(max(iv$nobs), format = "d", big.mark = "{,}")),
      "The exclusion restriction is not innocuous: absolute latitude correlates",
      "with institutions and development, which the gravity controls absorb only",
      "in part. This column bounds a sensitivity; it does not provide an",
      "alternative identification.")),
    "\\end{tabular}", "\\end{table}"),
    file.path(DIR_OUT, "tab_iv_latitude.tex"))
  cat("[write] tab_iv_latitude.tex\n")
}

##################################################################################
# R3. ASYMMETRIC FIXED EFFECTS -- closing the exporter cost channel
#
# Three schemes:
#   baseline  product x year                       bilateral identification
#   exporter  product x year + exporter x year     absorbs the exporter cost
#             term ENTIRELY, so g is identified by importer climate only
#   importer  product x year + importer x year     the mirror image
#
# This is not an ordinary robustness check: each scheme DESTROYS half of the
# identifying variation. The stability of the hot slope across the three brackets
# the bias of the exporter cost term.
#
# The slope is rounded to the unit, and no band is reported when the band covers
# zero: printing a precision on a quantity the table declares unmeasurable would
# contradict the table itself.
##################################################################################

cat("\n=== R3 : effets fixes asymétriques ===\n")

SCHEMES <- list(baseline = "fe_prod_year",
                exporter = c("fe_prod_year", "fe_exp_year"),
                importer = c("fe_prod_year", "fe_imp_year"))
SCHEME_LAB <- c(baseline = "Baseline", exporter = "Exporter$\\times$year",
                importer = "Importer$\\times$year")
REF_T <- c(5, 10, 15, 20, 25, 30)

res_oh_A <- setNames(lapply(names(SCHEMES), function(x) list()), names(SCHEMES))
for (n in seq_len(nrow(GRID))) {
  cfg <- SPEC_CONFIGS[[GRID$i[n]]]; w <- WINSOR_CONFIGS[[GRID$j[n]]]
  key <- sprintf("%s%d_wl%swu%s", cfg$spec, cfg$degree, w$lower, w$upper)
  cat(sprintf("  [%s]\n", key))
  dp <- prepare_diffs(BACI, climate, cfg$spec, cfg$degree, w$lower, w$upper,
                      scale = SCALE_PROD)
  dp[, fe_exp_year := paste(orig, year, sep = "_")]
  dp[, fe_imp_year := paste(dest, year, sep = "_")]
  for (sc in names(SCHEMES)) {
    r <- run_stage3(dp, cfg$spec, cfg$degree, controls = TRADE_STAGE3_CONTROLS,
                    fe_str = SCHEMES[[sc]])
    r$wval_lower <- w$lower; r$wval_upper <- w$upper
    res_oh_A[[sc]][[key]] <- r
  }
  rm(dp); gc(verbose = FALSE)
}

curves_A <- lapply(names(SCHEMES), function(sc) aic_curve(res_oh_A[[sc]], climate))
names(curves_A) <- names(SCHEMES)
mc_A <- lapply(names(SCHEMES), function(sc)
  ## MC_N_SIM_OH et non MC_N_SIM : a 300 tirages le verdict d'identification
  ## est decide par la graine (mesure du 07/09, voir 00_parameters.R).
  mc_run(res_oh_A[[sc]], climate, MC_N_SIM_OH, MC_X_GRID,
         MC_SEED[["prod_daily"]], label = paste0("A:", sc)))
names(mc_A) <- names(SCHEMES)

## Hot slope from 20 to 30 degrees C, with its Monte Carlo band.
slope_band <- function(mc) {
  s <- mc[, .(sl = (approx(x, y, 30, rule = 2)$y -
                    approx(x, y, 20, rule = 2)$y) / 10), by = sim_id]
  quantile(s$sl, c(0.05, 0.5, 0.95), na.rm = TRUE)
}
eff <- sapply(names(SCHEMES), function(sc)
  approx(curves_A[[sc]]$x, curves_A[[sc]]$y_avg, xout = REF_T, rule = 2)$y)
rownames(eff) <- paste0(REF_T, "C")
slp <- sapply(names(SCHEMES), function(sc) slope_band(mc_A[[sc]]))

## No band when it covers zero: the identifying variation has collapsed, and a
## band would put a precision on something unmeasurable.
fmt_slope <- function(v) {
  if (v[1] * v[3] <= 0) sprintf("%.0f \\; [---]", 100 * v[2])
  else sprintf("%.0f \\; [%.0f, %.0f]", 100 * v[2], 100 * v[1], 100 * v[3])
}

{
  L <- c(tex_header("01c_robustness.R (R3)"),
    "\\begin{table}[H]\\centering\\small",
    paste0("\\caption{Stage-1 productivity dose-response under asymmetric fixed ",
           "effects (Astier et al.\\ robustness). Relative effect (\\%, normalised ",
           "to its own maximum) at reference temperatures.}"),
    "\\label{tab:oh_asym_fe}",
    "\\begin{tabular}{lccc}", "\\toprule",
    paste0(" & ", paste(SCHEME_LAB[names(SCHEMES)], collapse = " & "), " \\\\"),
    "\\midrule",
    "\\multicolumn{4}{l}{\\textit{Relative effect on productivity (\\%, max $=0$)}}\\\\")
  for (i in seq_along(REF_T))
    L <- c(L, sprintf("\\quad %d$^\\circ$C & %.2f & %.2f & %.2f \\\\", REF_T[i],
                      100 * eff[i, "baseline"], 100 * eff[i, "exporter"],
                      100 * eff[i, "importer"]))
  L <- c(L, "\\midrule",
    sprintf("Hot-side slope 20$\\to$30$^\\circ$C (\\%%/$^\\circ$C) & %s & %s & %s \\\\",
            fmt_slope(slp[, "baseline"]), fmt_slope(slp[, "exporter"]),
            fmt_slope(slp[, "importer"])),
    "\\bottomrule",
    tex_note(4, c(
      "\\textit{Note.} Exporter$\\times$year fixed effects absorb the exporter",
      "capability term $Z(s)x(s)^{-\\theta}$ in full, so $g$ is identified off",
      "importer climate only; the importer scheme is the mirror. Stability of the",
      "hot-side slope brackets the exporter-side cost bias.",
      sprintf("Each column is a weighted average over the %d specifications with",
              length(res_oh_A[[1]])),
      "Akaike weights. Slopes are rounded to the unit, with the $[Q_{05},Q_{95}]$",
      "Monte-Carlo band in brackets; no band is reported where the band contains",
      "zero (---), since that scheme's identifying climate variation has",
      "collapsed and a band would put a precision on a quantity the table",
      "declares unmeasurable.")),
    "\\end{tabular}", "\\end{table}")
  writeLines(L, file.path(DIR_OUT, "tab_oh_asym_fe.tex"))
  cat("[write] tab_oh_asym_fe.tex\n")
}

##################################################################################
# R3bis. LA FIGURE DU MEME OBJET
#
# fig_oh_asym_fe.pdf. Le manuscrit la charge (l.~1961) et aucune etape ne la
# produisait : elle manquait au flux vivant.
#
# AUCUNE ESTIMATION NEUVE. `curves_A` et `mc_A` sont ceux que la table
# ci-dessus vient de resumer. La figure et la table sortent donc du MEME objet,
# et sont coherentes par construction -- pas par relecture.
#
# LES RUBANS SONT LES QUANTILES DES TIRAGES, pas une bande analytique : c'est ce
# que la legende du manuscrit annonce (Q05-Q95). Ils sont calcules ici sur
# MC_N_SIM_OH tirages ; la legende du .tex doit porter ce nombre et non 300.
##################################################################################

{
  ## ⚠️ SCHEME_LAB porte du LaTeX (« Exporter$\\times$year ») : il est fait pour
  ## la table, et une bande de facette le rendrait litteralement. Les etiquettes
  ## de la figure sont donc les memes en TEXTE BRUT, posees ici et nulle part
  ## ailleurs, dans l'ordre de SCHEMES.
  SCHEME_LAB_PLAIN <- c(baseline = "Baseline",
                        exporter = "Exporter × year",
                        importer = "Importer × year")
  stopifnot(identical(names(SCHEME_LAB_PLAIN), names(SCHEMES)))

  ## Enveloppe par schema : la mediane et les deux quantiles, a x commun.
  band <- rbindlist(lapply(names(SCHEMES), function(sc)
    mc_A[[sc]][, .(y_q05 = quantile(y, 0.05, na.rm = TRUE),
                   y_q50 = quantile(y, 0.50, na.rm = TRUE),
                   y_q95 = quantile(y, 0.95, na.rm = TRUE)),
               by = x][, scheme := SCHEME_LAB_PLAIN[[sc]]][]))
  ## La courbe tracee est la moyenne ponderee AIC, celle de la table -- pas la
  ## mediane des tirages. Les deux different, et confondre les deux ferait lire
  ## la figure comme un objet que la table ne rapporte pas.
  ctr <- rbindlist(lapply(names(SCHEMES), function(sc)
    data.table(x = curves_A[[sc]]$x, y = curves_A[[sc]]$y_avg,
               scheme = SCHEME_LAB_PLAIN[[sc]])))
  ## Ordre des facettes : celui de SCHEMES, donc celui des colonnes de la table.
  lev <- unname(SCHEME_LAB_PLAIN[names(SCHEMES)])
  band[, scheme := factor(scheme, levels = lev)]
  ctr[,  scheme := factor(scheme, levels = lev)]

  p_oh <- ggplot(band, aes(x = x)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(aes(ymin = 100 * y_q05, ymax = 100 * y_q95),
                fill = FAN_COLOR[["prod"]], alpha = 0.18) +
    geom_line(data = ctr, aes(y = 100 * y),
              colour = FAN_COLOR[["prod"]], linewidth = 0.7) +
    facet_wrap(~ scheme, nrow = 1) +
    ## PAS DE SOUS-TITRE. Il disait « AIC-weighted over 12 specifications »,
    ## ce que la legende du manuscrit CONTREDIT : les poids d Akaike sont
    ## degeneres, une seule specification porte le poids un, et la courbe
    ## tracee est celle-la. Un enonce faux grave dans une image ne se corrige
    ## pas a la relecture du .tex -- il survit a toutes les reecritures de
    ## legende. Le nombre de tirages, lui, est deja dans la legende via
    ## \EMPohMCdraws, donc le repeter ici serait un second chemin vers la
    ## meme valeur.
    labs(x = "Temperature (°C)", y = FAN_YLAB[["prod"]]) +
    ## L AXE EST COUPE A -200 %. Les rubans descendent a -310 % dans la seule
    ## facette importer, ou la variation identifiante s est effondree et ou la
    ## table ne rapporte AUCUNE bande. Laisser l axe suivre cette queue ecrase
    ## les trois courbes dans le tiers superieur du cadre et donne du relief a
    ## la facette qui en porte le moins. coord_cartesian ROGNE l affichage sans
    ## retirer d observation : les quantiles restent ceux des tirages.
    coord_cartesian(ylim = c(-200, NA)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          plot.subtitle = element_text(size = 8, colour = "grey30"))

  ggsave(file.path(DIR_OUT, "fig_oh_asym_fe.pdf"), p_oh,
         width = 10, height = 3.6)
  cat(sprintf("[write] fig_oh_asym_fe.pdf (%d schemas, B = %d)\n",
              length(SCHEMES), MC_N_SIM_OH))
  rm(band, ctr, p_oh, lev)
}

## ---- spread of the hot slope across the identified schemes ------------------
## 🚨 EMISSION DEPLACEE APRES R4 (2026-09-07). Le CSV des schemas s'ecrivait ICI,
## donc il ne pouvait porter que la tache A : les bras de R4 sont estimes plus
## bas, et leurs bandes MC -- que run_oh_controls calcule deja pour sa propre
## table -- etaient jetees. Le compte publie (\EMPohSchemesN = 8) vient du
## design complet A+B+C ; un CSV ecrit avant R4 ne pouvait pas le reproduire.
## On accumule ici, on emet apres R4.
sch_rows <- function(m, tache) {
  rbindlist(lapply(colnames(m), function(a) {
    v <- m[, a]
    data.table(task = tache, scheme = a,
               slope_q50 = 100 * v[2], slope_q05 = 100 * v[1],
               slope_q95 = 100 * v[3], identified = v[1] * v[3] > 0)
  }))
}
SCH <- sch_rows(slp, "A")

## L'EFFET AUX TEMPERATURES DE REFERENCE EST ARCHIVE, PLUS JETE (11/09/2026).
## `eff` porte deja les six lignes que la table imprime. Sans CSV, la seule
## trace en etait le .tex rendu -- et reconstruire tab_oh_bounds obligeait soit
## a recopier des valeurs d'un .tex (interdit), soit a relancer une estimation
## de cinquante minutes. On accumule comme SCH, on emet apres R4.
eff_rows <- function(m, tache) {
  rbindlist(lapply(colnames(m), function(a)
    data.table(task = tache, scheme = a, temp_C = REF_T,
               rel_effect_pct = 100 * as.numeric(m[, a]))))
}
EFF <- eff_rows(eff, "A")
cat(sprintf("  [R3] %d schemas accumules (tache A), %d points d'effet\n",
            nrow(SCH), nrow(EFF)))

##################################################################################
# R4. CONTROLS FOR FUNDAMENTALS AND INCOME
#
# Two exercises with the same logic: add the controls the criticism suggests, and
# show that the hot slope does not move.
#   B  fundamentals: human development index, origin and destination;
#   C  income: bilateral GDP enters DIRECTLY, although it is a bad control, being
#      itself a consequence of climate. The column exists to show that even that
#      control does not reverse the slope.
#
# This is a RESTRICTED-SAMPLE comparison: the controls are not available
# everywhere, so the baseline row of R4 is not the baseline row of R3.
##################################################################################

cat("\n=== R4 : contrôles de fondamentaux et revenu ===\n")

hdi_path <- file.path(DIR_OUT, "oh_hdi_country_year.csv")
run_oh_controls <- function(varie, label, output_path, note_extra,
                            mode = c("ajout", "retrait"),
                            noms_bras, entetes, tache) {
  mode <- match.arg(mode)
  ## 🚨 GARDE POSEE LE 2026-09-07 -- ELLE AURAIT ATTRAPE LA DEGENERESCENCE.
  ## Jusqu'ici cette fonction AJOUTAIT toujours le terme. Appelee avec
  ## "GDP_diff", qui figure DEJA dans TRADE_STAGE3_CONTROLS (00_parameters.R
  ## l.170), elle donnait deux bras aux memes controles : la table imprimait six
  ## paires de nombres egaux, sans lever la moindre erreur, et le manuscrit
  ## (l.1953) annoncait pendant ce temps qu'elle RETIRAIT le revenu.
  deja <- varie %in% TRADE_STAGE3_CONTROLS
  if (mode == "ajout" && any(deja))
    stop("run_oh_controls(ajout) : ", paste(varie[deja], collapse = ", "),
         " figure DEJA dans TRADE_STAGE3_CONTROLS. Les deux bras seraient ",
         "identiques et la table ne testerait rien. Utiliser mode = 'retrait'.")
  if (mode == "retrait" && !all(deja))
    stop("run_oh_controls(retrait) : ", paste(varie[!deja], collapse = ", "),
         " ne figure PAS dans TRADE_STAGE3_CONTROLS : il n'y a rien a retirer.")

  ## Bras 1 = la specification de reference ; bras 2 = celle qui varie.
  ctrl_de <- function(a) {
    if (a == noms_bras[1]) TRADE_STAGE3_CONTROLS
    else if (mode == "ajout") c(TRADE_STAGE3_CONTROLS, varie)
    else setdiff(TRADE_STAGE3_CONTROLS, varie)
  }

  res <- setNames(lapply(noms_bras, function(x) list()), noms_bras)
  for (n in seq_len(nrow(GRID))) {
    cfg <- SPEC_CONFIGS[[GRID$i[n]]]; w <- WINSOR_CONFIGS[[GRID$j[n]]]
    key <- sprintf("%s%d_wl%swu%s", cfg$spec, cfg$degree, w$lower, w$upper)
    dp <- prepare_diffs(BACI, climate, cfg$spec, cfg$degree, w$lower, w$upper,
                        scale = SCALE_PROD)
    if (!all(varie %in% names(dp))) return(invisible(NULL))
    ## Les deux bras sur le MEME echantillon : sinon l'ecart melangerait un
    ## effet de controle et un effet de selection.
    dp <- dp[complete.cases(dp[, ..varie])]
    for (a in noms_bras) {
      r <- run_stage3(dp, cfg$spec, cfg$degree, controls = ctrl_de(a))
      r$wval_lower <- w$lower; r$wval_upper <- w$upper
      res[[a]][[key]] <- r
    }
    rm(dp); gc(verbose = FALSE)
  }
  cu <- lapply(res, function(x) aic_curve(x, climate))
  mc <- lapply(names(res), function(a)
    ## MC_N_SIM_OH : meme raison qu'en R3.
    mc_run(res[[a]], climate, MC_N_SIM_OH, MC_X_GRID, MC_SEED[["prod_daily"]],
           label = paste0(label, ":", a)))
  names(mc) <- names(res)
  e <- sapply(names(res), function(a)
    approx(cu[[a]]$x, cu[[a]]$y_avg, xout = REF_T, rule = 2)$y)
  s <- sapply(names(res), function(a) slope_band(mc[[a]]))
  ## Meme geste qu'en R3 : les six points par bras sont accumules au lieu
  ## d'etre jetes apres impression. `<<-` parce que l'accumulateur vit au
  ## niveau du script, comme SCH, et que le contrat de retour de cette
  ## fonction -- les lignes de schemas -- ne doit pas bouger.
  EFF <<- rbind(EFF, eff_rows(e, tache))

  ## CONTROLE DE NON-DEGENERESCENCE. Deux bras qui rendent la meme courbe font
  ## une table qui ne teste rien : on le DIT au lieu de l'imprimer deux fois.
  if (isTRUE(all.equal(as.numeric(e[, 1]), as.numeric(e[, 2]),
                       tolerance = 1e-12)))
    warning("run_oh_controls(", label, ") : les deux bras rendent des effets ",
            "IDENTIQUES -- la table compare une specification a elle-meme.")

  L <- c(tex_header("01c_robustness.R (R4)"),
    "\\begin{table}[H]\\centering\\small",
    sprintf("\\caption{%s}", label),
    ## 🚨 CORRIGE LE 2026-09-07 : le prefixe tab_ est RETIRE avant de former
    ## l'etiquette. Sans cela on gravait \\label{tab:tab_oh_income_badcontrol}
    ## quand le manuscrit demande \\ref{tab:oh_income_badcontrol} -- renvoi non
    ## resolu, et invisible tant que le \\input restait commente.
    sprintf("\\label{tab:%s}",
            sub("^tab_", "", basename(tools::file_path_sans_ext(output_path)))),
    "\\begin{tabular}{lcc}", "\\toprule",
    sprintf(" & %s & %s \\\\", entetes[1], entetes[2]), "\\midrule",
    "\\multicolumn{3}{l}{\\textit{Relative effect on productivity (\\%, max $=0$)}}\\\\")
  for (i in seq_along(REF_T))
    L <- c(L, sprintf("\\quad %d$^\\circ$C & %.2f & %.2f \\\\", REF_T[i],
                      100 * e[i, noms_bras[1]], 100 * e[i, noms_bras[2]]))
  L <- c(L, "\\midrule",
    sprintf("Hot-side slope 20$\\to$30$^\\circ$C (\\%%/$^\\circ$C) & %s & %s \\\\",
            fmt_slope(s[, noms_bras[1]]), fmt_slope(s[, noms_bras[2]])),
    "\\bottomrule", tex_note(3, note_extra), "\\end{tabular}", "\\end{table}")
  writeLines(L, output_path)
  cat(sprintf("[write] %s\n", basename(output_path)))
  invisible(sch_rows(s, tache))
}

if (!isTRUE(OH_TACHE_B_HDI)) {
  ## Lot 56 : sautee PAR DECLARATION (00_parameters.R, OH_TACHE_B_HDI), plus par
  ## absence de fichier. Le lecteur ci-dessous n'a jamais tourne (colonnes attendues
  ## GID_0/hdi, fichier iso3/HDI).
  cat("  ⏭ R4/fundamentals (task B, HDI) SKIPPED BY DECLARATION: OH_TACHE_B_HDI = FALSE",
      "(00_parameters.R).\n")
} else if (file.exists(hdi_path)) {
  hdi <- fread(hdi_path)
  BACI <- merge(BACI, hdi[, .(orig = GID_0, year, hdi_orig = hdi)],
                by = c("orig", "year"), all.x = TRUE)
  BACI <- merge(BACI, hdi[, .(dest = GID_0, year, hdi_dest = hdi)],
                by = c("dest", "year"), all.x = TRUE)
  .r4b <- run_oh_controls(c("hdi_orig", "hdi_dest"),
    "Fundamental controls: human development index, origin and destination",
    file.path(DIR_OUT, "tab_oh_fundamental_controls.tex"),
    c("\\textit{Note.} The comparison is on the RESTRICTED sample for which the",
      "index is available in both countries, so the baseline column here is not",
      "the baseline of Table~\\ref{tab:oh_asym_fe}. Stability of the hot-side",
      "slope is the object of interest, not the level."),
    mode = "ajout", noms_bras = c("baseline", "raw"),
    entetes = c("Baseline", "With HDI"), tache = "B")
  if (!is.null(.r4b)) SCH <- rbind(SCH, .r4b)
} else {
  cat("  ⚠ oh_hdi_country_year.csv absent — R4/fondamentaux SAUTÉE.\n")
}

## 🚨 SENS INVERSE LE 2026-09-07. Cet appel AJOUTAIT GDP_diff a une liste qui le
## contient deja : les deux colonnes etaient identiques ligne pour ligne, et la
## table comparait une specification a elle-meme. Le test que le manuscrit
## decrit (l.1953, "drops the income term to bound the collider path") est le
## RETRAIT du revenu. Bras : withy = controles complets, noy = sans revenu.
.r4c <- run_oh_controls("GDP_diff",
  "Income as a control: a deliberately bad control",
  file.path(DIR_OUT, "tab_oh_income_badcontrol.tex"),
  c("\\textit{Note.} Bilateral income is a bad control: it is itself a",
    "consequence of climate, so conditioning on it removes part of the effect",
    "being estimated. The second column DROPS it, which bounds the collider",
    "path; the comparison is on the sample where income is observed, so both",
    "columns rest on the same observations."),
  mode = "retrait", noms_bras = c("withy", "noy"),
  entetes = c("With income", "Without income"), tache = "C")
if (!is.null(.r4c)) SCH <- rbind(SCH, .r4c)

##################################################################################
# R4bis. LE CSV DES SCHEMAS ET SON ETENDUE -- EMIS ICI, PAS DANS R3
#
# Ecrit apres R4 pour porter TOUTES les taches qui ont tourne. Les bandes des
# bras de R4 sont celles que run_oh_controls calcule pour sa propre table : rien
# n'est re-estime, on cesse simplement de les jeter.
#
# ⚠️ LE COMPTE EST UN RESULTAT, PAS UNE CIBLE. Il vaudra 5 tant que la tache B
# est sautee (A: 3 + C: 2), contre 8 au design complet. On ne comble pas l'ecart
# en comptant autrement.
##################################################################################

##################################################################################
# R4ter. LA TABLE FUSIONNEE DES BORNES OH
#
# Remplace tab_oh_asym_fe, tab_oh_income_badcontrol et tab_oh_hotslope_range :
# trois tables qui montraient le meme objet -- l'effet relatif aux temperatures
# de reference et la pente cote chaud -- sur des sous-ensembles de schemas qui
# se recouvraient. Une seule table, quatre colonnes.
#
# FONCTION PURE de ses deux data.tables, pour etre exercable sans rejouer une
# estimation de cinquante minutes : tools/reemit_from_archives.R l'extrait et
# lui passe oh_scheme_effects.csv et oh_hotslope_schemes.csv.
##################################################################################
emit_tab_oh_bounds <- function(eff_dt, sch_dt, output_path) {
  COLS <- list(list(task = "A", scheme = "baseline", head = "Baseline"),
               list(task = "A", scheme = "exporter", head = "Exporter$\\times$year"),
               list(task = "A", scheme = "importer", head = "Importer$\\times$year"),
               list(task = "C", scheme = "noy",      head = "Without income"))
  ## Une colonne absente est DITE, jamais remplacee par un blanc silencieux.
  keep <- vapply(COLS, function(c0)
    nrow(eff_dt[task == c0$task & scheme == c0$scheme]) > 0L, logical(1))
  if (!any(keep)) {
    cat("  [OMIS] aucune colonne disponible -> tab_oh_bounds.tex NON produite.\n")
    return(invisible(NULL))
  }
  if (!all(keep))
    cat(sprintf("  [!] colonnes absentes de oh_scheme_effects.csv : %s\n",
                paste(vapply(COLS[!keep], function(c0)
                  paste0(c0$task, "/", c0$scheme), character(1)), collapse = ", ")))
  COLS <- COLS[keep]
  temps <- sort(unique(eff_dt$temp_C))
  nc <- length(COLS) + 1L

  cell_eff <- function(c0, t) {
    v <- eff_dt[task == c0$task & scheme == c0$scheme & temp_C == t, rel_effect_pct]
    if (!length(v)) "---" else sprintf("%.2f", v[1])
  }
  ## MEME REGLE QUE LES TABLES REMPLACEES : pas de bande quand elle couvre zero.
  ## Le verdict vient de la colonne `identified`, calculee par R3/R4 ; il n'est
  ## pas refait ici, sinon deux endroits pourraient en decider differemment.
  slope_point <- function(c0) {
    y <- function(t) eff_dt[task == c0$task & scheme == c0$scheme & temp_C == t, rel_effect_pct]
    a <- y(30); b <- y(20)
    if (length(a) != 1L || length(b) != 1L)
      stop("[STOP] pente ponctuelle : effets a 20 et 30 C introuvables pour ",
           c0$task, "/", c0$scheme)
    (a - b) / 10
  }
  cell_slope <- function(c0) {
    r <- sch_dt[task == c0$task & scheme == c0$scheme]
    if (!nrow(r)) return("---")
    ## `return` EXPLICITE sur la branche non identifiee. Sans lui, R evalue le
    ## premier sprintf, le jette, et rend TOUJOURS la bande -- y compris pour un
    ## schema que la table declare non mesurable. Defaut attrape par le test de
    ## mutation ci-dessous, pas par la lecture.
    ## THE POINT IS THE POINT ESTIMATE (lot 34, 25/09/2026): (y30 - y20)/10 from
    ## the rows of this same table, the object the text cites as \EMPohSlopeLo/Hi.
    ## It was slope_q50, the Monte-Carlo median (-0.885 against -0.899): two
    ## decimals on the median would have printed -0.89 beside a -0.90 in the text.
    ## Bounds at two decimals since 25/09 (lot 32); the point follows.
    pt <- slope_point(c0)
    if (!isTRUE(as.logical(r$identified[1])))
      return(sprintf("%.2f \\; [---]", pt))
    sprintf("%.2f \\; [%.2f, %.2f]", pt, r$slope_q05[1], r$slope_q95[1])
  }

  L <- c(tex_header("01c_robustness.R (R4ter)",
                    "oh_scheme_effects.csv + oh_hotslope_schemes.csv"),
    "\\begin{table}[H]\\centering\\small",
    paste0("\\caption{Bounds on the productivity dose-response under asymmetric ",
           "fixed effects and without the income control. Relative effect (\\%, ",
           "normalised to each scheme's own maximum) at reference temperatures, ",
           "and the hot-side slope with its Monte-Carlo band.}"),
    "\\label{tab:oh_bounds}",
    sprintf("\\begin{tabular}{l%s}", paste(rep("c", length(COLS)), collapse = "")),
    "\\toprule",
    paste0(" & ", paste(vapply(COLS, function(c0) c0$head, character(1)),
                        collapse = " & "), " \\\\"),
    "\\midrule",
    sprintf("\\multicolumn{%d}{l}{\\textit{Relative effect on productivity (\\%%, max $=0$)}}\\\\", nc))
  for (t in temps)
    L <- c(L, sprintf("\\quad %g$^\\circ$C & %s \\\\", t,
                      paste(vapply(COLS, cell_eff, character(1), t = t),
                            collapse = " & ")))
  L <- c(L, "\\midrule",
    sprintf("Hot-side slope 20$\\to$30$^\\circ$C (\\%%/$^\\circ$C) & %s \\\\",
            paste(vapply(COLS, cell_slope, character(1)), collapse = " & ")),
    "\\bottomrule",
    tex_note(nc, c(
      "\\textit{Note.} Exporter$\\times$year fixed effects absorb the exporter",
      "capability term in full, so the response is identified off importer",
      "climate only; the importer scheme is the mirror. The last column drops",
      "bilateral income, which is a bad control: it is itself a consequence of",
      "climate. The slope is the point estimate, $(y_{30}-y_{20})/10$ from the",
      ## 22/09/2026 : la note donnait a ses DEUX tirets la raison qui ne vaut
      ## que pour importer x year (F joint < 1).
      ## 25-26/09/2026 (lots 34-35) : "\times" a UN antislash = TABULATION en R ;
      ## toujours deux antislashs ici.
      ## LOT 74 (01/10/2026, P2) : la note disait « under exporter x year effects
      ## the band includes zero » ; c'est desormais la SEULE bande qui exclut zero.
      ## La note dit la regle sans nommer le schema ; la bande imprimee le marque.
      ## Meme texte que 01c_reemit_tab_oh_bounds.R (preuve de non-divergence la-bas).
      "rows above, with the $[Q_{05},Q_{95}]$ Monte-Carlo band in brackets where",
      "the band excludes zero; where it contains zero (---) no band is reported,",
      "since it would not bound the slope. Under importer$\\times$year effects the",
      "identifying climate variation has, in addition, collapsed.")),
    "\\end{tabular}", "\\end{table}")
  writeLines(L, output_path)
  cat(sprintf("  [write] %s (%d colonnes, %d temperatures)\n",
              basename(output_path), length(COLS), length(temps)))
  invisible(L)
}

{
  sch <- SCH
  ## LE NOMBRE DE TIRAGES VOYAGE AVEC LES VERDICTS. La legende de la figure et
  ## la note du texte le citent, et il vient de passer de 300 a 3000 : reste en
  ## dur dans le .tex, il redeviendrait faux au prochain changement, et
  ## silencieusement. Colonne CONSTANTE -- c'est un parametre du bloc, pas une
  ## propriete du schema -- et c'est justement pourquoi l'emetteur peut la lire
  ## sur n'importe quelle ligne apres avoir controle qu'elle ne varie pas.
  sch[, n_sim := MC_N_SIM_OH]
  fwrite(sch, file.path(DIR_OUT, "oh_hotslope_schemes.csv"))
  ## L'effet aux six temperatures, toutes taches confondues. C'est la moitie
  ## haute de tab_oh_bounds ; sans lui elle n'existait que dans le .tex rendu.
  fwrite(EFF, file.path(DIR_OUT, "oh_scheme_effects.csv"))
  emit_tab_oh_bounds(EFF, sch, file.path(DIR_OUT, "tab_oh_bounds.tex"))
  cat(sprintf("  [write] oh_scheme_effects.csv (%d lignes, %d schemas x %d temperatures)\n",
              nrow(EFF), uniqueN(EFF[, .(task, scheme)]), length(REF_T)))
  id <- sch[identified == TRUE]
  writeLines(c(
    tex_header("01c_robustness.R (R4bis)", "oh_hotslope_schemes.csv"),
    "\\begin{table}[H]\\centering\\small",
    "\\caption{Range of the hot-side productivity slope across identified fixed-effect schemes.}",
    "\\label{tab:oh_hotslope_range}",
    "\\begin{tabular}{lc}", "\\toprule",
    "Quantity & Value \\\\", "\\midrule",
    sprintf("Schemes estimated & %d \\\\", nrow(sch)),
    sprintf("Schemes with an identified slope & %d \\\\", nrow(id)),
    sprintf("Hot-side slope range (\\%%/$^\\circ$C) & %s \\\\",
            if (nrow(id) == 0L) "\\textit{no identified scheme}"
            else sprintf("$[%.2f,\\; %.2f]$",
                         min(id$slope_q50), max(id$slope_q50))),
    "\\bottomrule",
    tex_note(2, c(
      "\\textit{Note.} A scheme counts as identified when its Monte-Carlo band on",
      "the $20\\to30^\\circ$C slope excludes zero. Schemes whose band straddles",
      "zero are estimated but not counted: their identifying climate variation",
      "has collapsed under the fixed effects.")),
    "\\end{tabular}", "\\end{table}"),
    file.path(DIR_OUT, "tab_oh_hotslope_range.tex"))
  cat(sprintf("[write] oh_hotslope_schemes.csv + tab_oh_hotslope_range.tex (%d schemas, %d identifies)\n",
              nrow(sch), nrow(id)))
  print(sch[, .(task, scheme, slope_q50, slope_q05, slope_q95, identified)])
}

##################################################################################
# R5. CORRECTION ASTIER ET AL. (2026)
#
# Deux modifications structurelles :
#   (i)  GDP_diff purged of its climate component;
#   (ii) the structural multilateral resistance differential.
#
# On the AMENITY side only (i) applies: the pair fixed effect already absorbs the
# permanent component of the bilateral differential.
#
# These two artefacts come from a script outside the runbook. When they are
# absent the step is skipped and says so: a mechanism that can do nothing must
# say whether it acted.
##################################################################################

cat("\n=== R5 : correction Astier et al. ===\n")

p_mrt <- file.path(DIR_TMP, "mrt_structural.rds")
p_gdp <- file.path(DIR_TMP, "gdp_purged.rds")

if (!file.exists(p_mrt) || !file.exists(p_gdp)) {
  cat(sprintf("  ⛔ R5 SAUTÉE : %s absent.\n",
              paste(basename(c(p_mrt, p_gdp))[!file.exists(c(p_mrt, p_gdp))],
                    collapse = " et ")))
  cat("     Ces artefacts viennent de mrt_compute.r, qui n'est dans aucune étape\n")
  cat("     du runbook. Sans eux, tab_robustness_astier_*.tex ne sont PAS produites.\n")
} else {
  mrt <- as.data.table(readRDS(p_mrt)); gp <- as.data.table(readRDS(p_gdp))

  astier_arms <- function(dt, spec, deg) {
    dt <- merge(dt, gp[, .(orig = GID_0, dest, year, GDP_diff_purged)],
                by = c("orig", "dest", "year"), all.x = TRUE)
    dt <- merge(dt, mrt[, .(orig = GID_0, year, phi_orig = phi)],
                by = c("orig", "year"), all.x = TRUE)
    dt <- merge(dt, mrt[, .(dest = GID_0, year, phi_dest = phi)],
                by = c("dest", "year"), all.x = TRUE)
    dt[, mrt_diff := phi_dest - phi_orig]
    dt <- dt[is.finite(GDP_diff_purged) & is.finite(mrt_diff)]
    base <- setdiff(TRADE_STAGE3_CONTROLS, "GDP_diff")
    list(
      `(1) Baseline`         = run_stage3(dt, spec, deg, controls = c(base, "GDP_diff")),
      `(2) + MRT`            = run_stage3(dt, spec, deg, controls = c(base, "GDP_diff", "mrt_diff")),
      `(3) + purged income`  = run_stage3(dt, spec, deg, controls = c(base, "GDP_diff_purged")),
      `(4) Astier-compatible`= run_stage3(dt, spec, deg, controls = c(base, "GDP_diff_purged", "mrt_diff")))
  }

  emit_astier <- function(arms, output_path, caption, tex_label, note) {
    cn <- names(arms[[1]])
    body <- sapply(cn, function(a) {
      ok <- sapply(arms, function(g) !is.null(g[[a]]))
      p  <- sapply(arms[ok], function(g) g[[a]]$wald_joint$p)
      b1 <- sapply(arms[ok], function(g) {
        cc <- grep("^diff_", names(coef(g[[a]]$model)), value = TRUE)[1]
        coef(g[[a]]$model)[cc] })
      c(sprintf("%.4f", median(b1, na.rm = TRUE)),
        ifelse(min(p, na.rm = TRUE) < 0.001, "$<0.001$",
               sprintf("%.3f", min(p, na.rm = TRUE))))
    })
    writeLines(c(
      tex_header("01c_robustness.R (R5)"),
      "\\begin{table}[H]\\centering\\small",
      sprintf("\\caption{%s}", caption), sprintf("\\label{tab:%s}", tex_label),
      sprintf("\\begin{tabular}{l%s}", paste(rep("c", length(cn)), collapse = "")),
      "\\toprule", paste0(" & ", paste(cn, collapse = " & "), " \\\\"), "\\midrule",
      paste0("Median $\\hat\\beta_1$ across specifications & ",
             paste(body[1, ], collapse = " & "), " \\\\"),
      paste0("Smallest joint $p$-value & ", paste(body[2, ], collapse = " & "), " \\\\"),
      "\\bottomrule", tex_note(length(cn) + 1L, note),
      "\\end{tabular}", "\\end{table}"), output_path)
    cat(sprintf("[write] %s\n", basename(output_path)))
  }

  arms_prod <- grid_iter(function(spec, deg, wl, wu, key) {
    cat(sprintf("  [prod %s]\n", key))
    astier_arms(prepare_diffs(BACI, climate, spec, deg, wl, wu,
                              scale = SCALE_PROD), spec, deg)
  })
  emit_astier(arms_prod, file.path(DIR_OUT, "tab_robustness_astier_prod.tex"),
    "Productivity dose-response under the Astier et al.\\ (2026) correction",
    "robustness_astier_prod",
    c("\\textit{Note.} Column (1) is the published baseline. Column (2) adds the",
      "structural multilateral-resistance differential $\\Phi_d-\\Phi_o$; column (3)",
      "replaces bilateral income by its climate-orthogonal residual; column (4)",
      "applies both. Entries summarise the twelve specifications: the median",
      "first-coefficient estimate and the smallest joint $p$-value.",
      "The exporter-side cost term is the object of the critique; it is bounded",
      "here rather than eliminated."))

  step2 <- readRDS(file.path(DIR_TMP, "step2_migration.rds"))
  clim_a <- readRDS(file.path(DIR_TMP, "climate_long.rds"))
  arms_amen <- grid_iter(function(spec, deg, wl, wu, key) {
    cat(sprintf("  [amenity %s]\n", key))
    dt <- prepare_diffs(step2, clim_a, spec, deg, wl, wu,
                        fe_cols = c("fe_orig_dest", "fe_year"),
                        sign_dest_minus_orig = TRUE, scale = Omega)
    dt <- merge(dt, gp[, .(orig = GID_0, dest, year, GDP_diff_purged)],
                by = c("orig", "dest", "year"), all.x = TRUE)
    dt <- dt[is.finite(GDP_diff_purged)]
    list(`(1) Baseline` = run_stage3(dt, spec, deg,
           controls = AMENITY_STAGE3_CONTROLS, fe_str = c("fe_orig_dest", "fe_year")),
         `(2) Purged income` = run_stage3(dt, spec, deg,
           controls = c("GDP_diff_purged", "populationdiff"),
           fe_str = c("fe_orig_dest", "fe_year")))
  })
  emit_astier(arms_amen, file.path(DIR_OUT, "tab_robustness_amenity_astier.tex"),
    "Amenity dose-response under the Astier et al.\\ (2026) correction",
    "robustness_amenity_astier",
    c("\\textit{Note.} For the amenity channel the origin$\\times$destination fixed",
      "effect already absorbs the permanent component of the bilateral climate",
      "differential, so only the income correction applies: bilateral income is",
      "replaced by its climate-orthogonal residual. The null result of",
      "Section~\\ref{sec:amenity_estimation} is unchanged."))
}



##################################################################################
# R6. SPECIFICATION INVARIANCE -- the quantity the counterfactual integrates
#
# WHAT THIS SECTION REPLACES. The specification used to be justified by AIC
# selection. At this sample size an Akaike weight of one is a property of SAMPLE
# SIZE, not evidence about functional form: the log-likelihood grows fast enough
# with N that the best specification crushes the others whatever the real gap
# between them.
#
# The argument becomes instead: the quantity the counterfactual integrates -- the
# productivity drop at 30 degrees C, relative to the maximum of each curve --
# does not depend on the specification chosen. That is an invariance, not a
# selection.
#
# AND THAT IS WHY NO ROW IS MARKED. An invariance exhibit that names a winner
# puts back, through the layout, the argument the text removes. No AIC marker, no
# bold, no highlighted row.
#
# PRODUIT
#   pe_spec_ranking.csv   -> read by 02d_emit_latex.R for tab_spec_invariance.tex
##################################################################################

cat("\n=== R6 : invariance de spécification ===\n")

## The two temperatures at which the quantity is evaluated. 30 degrees C is the
## one the manuscript cites; 25 is there to show that the ranking does not rest
## on an isolated point of the hot tail.
PE_EVAL_T <- c(30, 25)

res_pe <- grid_iter(function(spec, deg, wl, wu, key) {
  cat(sprintf("  [%s]\n", key))
  dp <- prepare_diffs(BACI, climate, spec, deg, wl, wu, scale = SCALE_PROD)
  r  <- run_stage3(dp, spec, deg, controls = TRADE_STAGE3_CONTROLS)
  ## get_effect_curve reads wval_lower/wval_upper off the object; without them
  ## winsor_bounds receives NULL and seq() fails.
  r$wval_lower <- wl; r$wval_upper <- wu
  rm(dp); gc(verbose = FALSE)
  r
})

## ---- the quantity, specification by specification ---------------------------
## EACH CURVE IS RELATIVE TO ITS OWN MAXIMUM. That is the convention of
## get_effect_curve, and it is what makes the twelve comparable: a common
## normalisation at a fixed temperature would make the ranking depend on the
## anchor point, which is exactly what this exercise sets out to rule out.
pe_rows <- rbindlist(lapply(names(res_pe), function(nm) {
  cc <- get_effect_curve(res_pe[[nm]], climate)
  if (is.null(cc)) {
    ## A mechanism that does nothing must say whether it acted.
    cat(sprintf("  [!] %s : aucune courbe (coefficients absents) — ÉCARTÉE\n", nm))
    return(NULL)
  }
  ## rule = 1, never 2: beyond the estimated support the value does not exist and
  ## NA is returned. Extrapolating by the edge value would put a number where
  ## there is no observation.
  v <- approx(cc$x, cc$y, xout = PE_EVAL_T, rule = 1)$y
  data.table(model_name = nm,
             warm_drop_30 = v[1], warm_drop_25 = v[2],
             support_lo = min(cc$x), support_hi = max(cc$x))
}))

if (!nrow(pe_rows)) {
  cat("  [OMIS] aucune spécification exploitable -> pe_spec_ranking.csv NON produit.\n")
} else {
  n_hors <- sum(is.na(pe_rows$warm_drop_30))
  if (n_hors) {
    ## A support that does not reach 30 degrees C is not an anomaly: it is what
    ## the tightest winsorisation produces. But the specification cannot enter a
    ## table whose column is "at 30 degrees C".
    cat(sprintf("  [!] %d spécification(s) dont le support n'atteint pas %g °C :\n",
                n_hors, PE_EVAL_T[1]))
    print(pe_rows[is.na(warm_drop_30), .(model_name, support_lo, support_hi)])
    pe_rows <- pe_rows[!is.na(warm_drop_30)]
  }

  ## Ordered by the quantity itself, from smallest to largest in absolute value.
  ## NOT by merit: the table does not rank, it shows that the spread is narrow.
  setorder(pe_rows, -warm_drop_30)
  fwrite(pe_rows, file.path(DIR_OUT, "pe_spec_ranking.csv"))

  n_spec <- nrow(pe_rows)
  n_same <- sum(pe_rows$warm_drop_30 < 0)
  share  <- 100 * max(pe_rows$warm_drop_30) / min(pe_rows$warm_drop_30)
  cat(sprintf("  [write] pe_spec_ranking.csv (%d specs)\n", n_spec))
  cat(sprintf("  à %g °C : de %.1f à %.1f points de log | %d/%d de même signe\n",
              PE_EVAL_T[1], 100 * min(pe_rows$warm_drop_30),
              100 * max(pe_rows$warm_drop_30), n_same, n_spec))
  cat(sprintf("  la plus petite en valeur absolue retient %.0f %% de la plus grande\n",
              share))
  ## This percentage feeds the invariance macro. It is EMITTED, not typed into
  ## the .tex: if the invariance degraded, the rendered text would move with it.
  fwrite(data.table(macro = "EMPspecInvarShare", value = round(share, 0)),
         file.path(DIR_OUT, "pe_spec_invariance_share.csv"))
}


cat("\n=== 04_robustness.R terminé ===\n")
