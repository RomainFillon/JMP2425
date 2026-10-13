##################################################################################
# 01b_amenity_drf.R -- damage function for AMENITIES, through migration.
#
# ESTIMATED EQUATION
#   log mu_t(r,s) - log mu_t(r,r)
#     = -(lambda/Omega)[log L_t(s) - log L_t(r)]
#       + (1/Omega)[log y_t(s) - log y_t(r)]
#       + (1/Omega)[f(T_t(s), zeta_a) - f(T_t(r), zeta_a)]
#       + psi_rs + phi_t + eps
#   PPML, with origin x destination and year fixed effects, standard errors
#   clustered on origin and destination.
#   The origin is r and the destination is s, so the difference is dest - orig,
#   the opposite orientation from the trade channel.
#
# Reads:  data_inputs/abel_cohen_migration, WB_pop_country,
#         WB_gdppercapita_country, EMDAT, gravity_CEPII,
#         intermediate_outputs/ERA_tas_distribution_grid_annual.csv
# Writes:
#   .tex   tab_stage1_migration.tex, tab_exogeneity_mrt_migration.tex,
#          tab_exogeneity_mrt_annual_migration.tex,
#          tab_robustness_amenity_formatted.tex,
#          tab_robustness_amenity_annual_formatted.tex
#   .pdf   MC_amenity_{daily,annual}_fan_paper.pdf
#   .csv   MC_damage_amenity_{daily,annual}_quantiles.csv  -> GE inputs
#          stage1_migration_diagnostics.csv                -> stage 1 macros
#          amenity_joint_tests{,_annual}.csv               -> joint test macros
#          aic_weights_Amenity_{Daily,Annual}.csv
#          mde_amenity_{inputs,fixef_sizes}.csv            -> power macros
#   .rds   climate_long.rds                                (read back by the fans)
# Needs:  00_parameters.R and lib_drf.R, sourced below.
##################################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(fixest);     library(splines); library(MASS)
  library(readxl);     library(lubridate); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

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


##################################################################################
# 1. DATA
##################################################################################

## ---- 1.1 bilateral migration flows (Abel and Cohen 2019) --------------------
bilat_mig <- fread(file.path(DIR_IN, "abel_cohen_migration/bilat_mig.csv"))
bilat_mig <- bilat_mig[, .(da_pb_closed = sum(da_pb_closed, na.rm = TRUE)),
                       by = .(orig, dest, year0)]
bilat_mig[, year := recode_year5(year0)]
bilat_mig[, nb_years := uniqueN(year), by = .(orig, dest)]
bilat_mig <- bilat_mig[nb_years == max(nb_years)][, nb_years := NULL]

## ---- 1.2 GDP per capita and population (World Bank) -------------------------
gdp <- load_wb(file.path(DIR_IN,
                 "WB_gdppercapita_country/total/total_GDP_country_constant2015USD.csv"),
               "gdppercapita", period = "quinquennial",
               year_min = 1985L, year_max = 2019L,
               periods_expected = PERIODS_EXPECTED)
gdp[, log_GDP := log(gdppercapita)]

pop <- load_wb(file.path(DIR_IN, "WB_pop_country/population_data.csv"),
               "population", period = "quinquennial",
               year_min = 1985L, year_max = 2019L,
               periods_expected = PERIODS_EXPECTED)
pop[, log_pop := log(population)]

## ---- 1.3 gravity (CEPII) and disasters (EM-DAT) -----------------------------
gravity <- fread(file.path(DIR_IN, "gravity_CEPII/Gravity_V202211.csv"),
                 select = c("year", "iso3_o", "iso3_d", "distw_harmonic", "distcap",
                            "contig", "comlang_off", "sibling_ever",
                            "col_dep_ever", "fta_wto"))
gravity <- gravity[year %in% 1990:2019 & !is.na(comlang_off)]
gravity[, year := recode_year5(year)]
gravity <- gravity[, lapply(.SD, mean, na.rm = TRUE), by = .(iso3_o, iso3_d, year),
                   .SDcols = setdiff(names(gravity), c("iso3_o", "iso3_d", "year"))]
setnames(gravity, c("iso3_o", "iso3_d"), c("orig", "dest"))

emdat_file <- list.files(file.path(DIR_IN, "EMDAT"), pattern = "\\.xlsx$",
                         full.names = TRUE)[1]
EMDAT <- as.data.table(read_excel(emdat_file, sheet = 1)) %>%
  dplyr::select(`Start Year`, ISO, `Start Month`, `Start Day`,
                `End Month`, `End Day`, `Total Affected`) %>%
  setNames(c("year", "orig", "start_month", "start_day",
             "end_month", "end_day", "total")) %>%
  mutate(duration_days = pmax(0, as.numeric(
    make_date(year, end_month, end_day) - make_date(year, start_month, start_day) + 1),
    na.rm = TRUE)) %>%
  group_by(orig, year) %>%
  summarise(nb_events = n(),
            total_duration_days = sum(duration_days, na.rm = TRUE),
            total_affected      = sum(total, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()
EMDAT[, year := recode_year5(year)]
EMDAT <- EMDAT[, lapply(.SD, mean, na.rm = TRUE), by = .(orig, year),
               .SDcols = DISASTER_CONTROLS]

## ---- 1.4 assembling the bilateral panel -------------------------------------
step2 <- copy(bilat_mig)
step2 <- merge(step2, gdp[, .(GID_0, year, log_GDP_orig = log_GDP)],
               by.x = c("orig", "year"), by.y = c("GID_0", "year"), all.x = TRUE)
step2 <- merge(step2, gdp[, .(GID_0, year, log_GDP_dest = log_GDP)],
               by.x = c("dest", "year"), by.y = c("GID_0", "year"), all.x = TRUE)
step2[, GDP_diff := (log_GDP_dest - log_GDP_orig) / Omega]
step2 <- step2[!is.na(GDP_diff)]

step2 <- merge(step2, pop[, .(GID_0, year, log_pop_orig = log_pop, pop_orig = population)],
               by.x = c("orig", "year"), by.y = c("GID_0", "year"), all.x = TRUE)
step2 <- merge(step2, pop[, .(GID_0, year, log_pop_dest = log_pop)],
               by.x = c("dest", "year"), by.y = c("GID_0", "year"), all.x = TRUE)
step2[, populationdiff := -Lambda * (log_pop_dest - log_pop_orig) / Omega]

step2 <- merge(step2, gravity, by = c("orig", "dest", "year"), all.x = TRUE)
step2 <- merge(step2, EMDAT,   by = c("orig", "year"),         all.x = TRUE)
step2[is.na(nb_events), (DISASTER_CONTROLS) := 0]

## Dependent variable: mu(r,s)/mu(r,r), the flow relative to the non-migrants.
step2[, total_out := sum(da_pb_closed, na.rm = TRUE), by = .(orig, year)]
step2[, stayers   := pop_orig - total_out]
step2[, explained := da_pb_closed / stayers]
step2 <- step2[!is.na(explained) & !is.na(GDP_diff) & !is.na(populationdiff)]
saveRDS(step2, file.path(DIR_TMP_OUT, "step2_migration.rds"))

## ---- 1.5 climate: daily distribution aggregated to the country --------------
## The input carries 24 temperature bins per country and year. It is produced
## outside this package, by the server pipeline.
climate <- fread(file.path(DIR_TMP, "ERA_tas_annual.csv"))
setnames(climate, c("GID_0", "year", "value", "count", "share"))
climate <- climate[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
saveRDS(climate, file.path(DIR_TMP_OUT, "climate_long.rds"))

## Guardrail: a quinquennial object has already been written under this name.
{
  yrs <- sort(unique(climate$year))
  if (length(yrs) < 20 || any(diff(yrs) > 1L))
    stop(sprintf("climate_long.rds n'est PAS au pas annuel : %d années, pas max %d.",
                 length(yrs), max(diff(yrs))))
  cat(sprintf("[ok] climat : %d années au pas annuel, %d lignes, %d pays\n",
              length(yrs), nrow(climate), uniqueN(climate$GID_0)))
}

##################################################################################
# 3. STAGE 1 -- GRAVITY WITHOUT TEMPERATURE
#
# Bilateral migration costs m(r,s) are unobserved, so this stage checks that the
# gravity controls proxy them in a jointly informative way.
# Fixed effects: origin x year and destination x year.
##################################################################################

step2_data <- readRDS(file.path(DIR_TMP, "step2_migration.rds"))

run_stage1 <- function(data) {
  dt <- as.data.frame(data)
  dt$fe_orig_year <- paste(dt$orig, dt$year, sep = "_")
  dt$fe_dest_year <- paste(dt$dest, dt$year, sep = "_")
  controls <- c(GRAVITY_CONTROLS, DISASTER_CONTROLS, AMENITY_STAGE3_CONTROLS)
  feglm(as.formula(paste0("explained ~ ", paste(controls, collapse = " + "))),
        data = dt, family = poisson(), cluster = c("orig", "dest"),
        fixef = c("fe_orig_year", "fe_dest_year"))
}

cat("=== Stage 1 ===\n")
s1 <- run_stage1(step2_data)

emit_stage1_table <- function(mod, data, output_path) {
  sq_cor <- round(fitstat(mod, "sq.cor")[[1]], 3)
  ll     <- round(as.numeric(logLik(mod)), 0)
  n_obs  <- nobs(mod)
  gc_in  <- GRAVITY_CONTROLS[GRAVITY_CONTROLS %in% names(coef(mod))]
  w      <- wald(mod, gc_in)
  disp   <- round(sum((mod$y - fitted(mod))^2 / pmax(fitted(mod), 1e-10)) /
                    (length(mod$y) - length(coef(mod))), 3)
  zeros  <- round(mean(data$explained == 0, na.rm = TRUE), 3)

  ## CSV keyed by a STABLE IDENTIFIER, never by a display label. Parsing the
  ## rendered table by regex breaks as soon as the table changes shape.
  fwrite(data.frame(
    key   = c("sq_cor", "loglik", "nobs", "wald_f", "wald_p", "dispersion", "share_zeros"),
    value = c(sq_cor, ll, n_obs, w$stat, w$p, disp, zeros)),
    file.path(DIR_OUT, "stage1_migration_diagnostics.csv"))

  panel_a <- data.frame(
    s = c("Pseudo-$R^2$ (sq. corr.)", "Log-likelihood", "Observations",
          "Wald $F$ (gravity controls)", "$p$-value (Wald)",
          "PPML dispersion", "Share zero flows"),
    v = c(sprintf("%.3f", sq_cor), formatC(ll, format = "d", big.mark = ","),
          formatC(n_obs, format = "d", big.mark = ","), sprintf("%.2f", w$stat),
          sprintf("%.4f", w$p), sprintf("%.3f", disp), sprintf("%.3f", zeros)))

  nm_map <- c(distw_harmonic = "Log distance (harmonic)",
              distcap = "Log distance (capitals)", comlang_off = "Common language",
              col_dep_ever = "Colonial link", fta_wto = "FTA (WTO)",
              nb_events = "N. disaster events",
              total_duration_days = "Total disaster duration (days)",
              total_affected = "Total affected (disaster)",
              GDP_diff = "$(1/\\Omega)\\cdot\\log(y_d/y_o)$",
              populationdiff = "$-(\\lambda/\\Omega)\\cdot\\log(L_d/L_o)$")
  cf <- coef(mod); se <- sqrt(diag(vcov(mod)))
  keep <- names(cf)[names(cf) %in% names(nm_map)]
  panel_b <- data.frame(v = nm_map[keep], c = sprintf("%.4f", cf[keep]),
                        s = sprintf("(%.4f)", se[keep]))

  writeLines(c(
    "%% GENERATED by 01b_amenity_drf.R -- DO NOT EDIT BY HAND.",
    "\\begin{table}[H]", "\\centering\\small",
    "\\caption{Stage 1 gravity regression diagnostics --- Migration}",
    "\\label{tab:stage1_migration}", "\\bigskip",
    "\\textbf{Panel A: Diagnostics}\\\\[4pt]",
    "\\begin{tabular}{lc}", "\\hline\\hline", "Statistic & Value \\\\", "\\hline",
    apply(panel_a, 1, function(r) paste0(r[1], " & ", r[2], " \\\\")),
    "\\hline\\hline",
    "\\multicolumn{2}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    "PPML with origin$\\times$year and destination$\\times$year fixed effects.",
    "Gravity controls: log distance (harmonic and capital-to-capital), common",
    "official language, colonial link, WTO free trade agreement.",
    "Dispersion $=\\sum(y_i-\\hat y_i)^2/\\hat y_i/(N-K)$.}",
    "\\end{tabular}", "\\par\\bigskip",
    "\\textbf{Panel B: Gravity control coefficients}\\\\[4pt]",
    "\\begin{tabular}{lcc}", "\\hline\\hline",
    "Variable & Coefficient & Std. Error \\\\", "\\hline",
    apply(panel_b, 1, function(r) paste0(r[1], " & ", r[2], " & ", r[3], " \\\\")),
    "\\hline\\hline",
    "\\multicolumn{3}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    "Standard errors clustered by origin and destination.}",
    "\\end{tabular}", "\\end{table}"), output_path)
  cat(sprintf("[write] %s\n", basename(output_path)))
}

emit_stage1_table(s1, step2_data, file.path(DIR_OUT, "tab_stage1_migration.tex"))

##################################################################################
# 4. STAGE 3 -- THE DAILY ARM AND THE ANNUAL ARM
#
# WITHIN-PAIR identification: pair and year fixed effects. That is what licenses
# the pre-trend test below.
##################################################################################

step2_mrt <- step2_data[year >= YEAR_MIN_STAGE3]

estimate_grid <- function(is_annual) {
  res <- list()
  for (cfg in SPEC_CONFIGS) for (w in WINSOR_CONFIGS) {
    key <- sprintf("%s%d_wl%swu%s", cfg$spec, cfg$degree, w$lower, w$upper)
    cat(sprintf("  Stage 3 [%s] %s\n", if (is_annual) "annual" else "daily", key))
    dp <- prepare_diffs(step2_mrt, climate, cfg$spec, cfg$degree,
                        w$lower, w$upper, is_annual,
                        fe_cols = c("fe_orig_dest", "fe_year"),
                        sign_dest_minus_orig = TRUE, scale = Omega)
    r <- run_stage3(dp, cfg$spec, cfg$degree, is_annual,
                    controls = AMENITY_STAGE3_CONTROLS,
                    fe_str = c("fe_orig_dest", "fe_year"))
    r$wval_lower <- w$lower; r$wval_upper <- w$upper
    res[[key]] <- r
  }
  res
}

cat("=== Stage 3 : grille de 12 spécifications, deux bras ===\n")
results_daily  <- estimate_grid(is_annual = FALSE)
results_annual <- estimate_grid(is_annual = TRUE)

## The joint tests are extracted DIRECTLY from the models: a value computed and
## then discarded before writing is indistinguishable, once published, from a
## value that was invented.
joint_tests <- function(res) {
  rbindlist(lapply(names(res), function(nm) {
    r <- res[[nm]]
    data.table(model_name = nm, spec = r$spec, degree = r$degree,
               fval_joint  = r$wald_joint$stat,  pval_joint  = r$wald_joint$p,
               fval_nonlin = r$wald_nonlin$stat, pval_nonlin = r$wald_nonlin$p)
  }))
}
jt_daily  <- joint_tests(results_daily)
jt_annual <- joint_tests(results_annual)
fwrite(jt_daily,  file.path(DIR_OUT, "amenity_joint_tests.csv"))
fwrite(jt_annual, file.path(DIR_OUT, "amenity_joint_tests_annual.csv"))
cat(sprintf("[write] amenity_joint_tests{,_annual}.csv — F joint dans [%.2f;%.2f] et [%.2f;%.2f]\n",
            min(jt_daily$fval_joint), max(jt_daily$fval_joint),
            min(jt_annual$fval_joint), max(jt_annual$fval_joint)))

##################################################################################
# 5. PLACEBO -- PRE-TREND ON THE LEADS
#
# The contemporaneous bilateral differential is replaced by the differential ONE
# PERIOD AHEAD. Under the exogeneity condition, future climate does not explain
# present migration.
##################################################################################

## ONE function, TWO arms, and `is_annual` is the ONLY thing that differs between
## the two calls at the bottom. The controls, the direction of the difference
## (dest - orig, the opposite of the trade channel) and the scale (Omega) are
## written in the BODY, not passed by the caller: two callers restating them is
## exactly how one arm ends up estimated without the gravity controls while
## nothing reports it.
emit_placebo_table <- function(output_path, is_annual = FALSE) {
  rows <- rbindlist(lapply(SPEC_CONFIGS, function(cfg) rbindlist(lapply(WINSOR_CONFIGS, function(w) {
    ## Same branch as prepare_diffs(): the annual basis is NOT the daily basis
    ## evaluated at the mean, and its columns are named annual_poly_ /
    ## annual_spline_ -- never annual_polynomial_.
    if (is_annual) {
      cw     <- make_annual_terms(climate, cfg$spec, cfg$degree, w$lower, w$upper)
      prefix <- if (cfg$spec == "polynomial") "annual_poly_" else "annual_spline_"
    } else {
      cw     <- make_climate_terms_wide(climate, cfg$spec, cfg$degree, w$lower, w$upper)
      prefix <- paste0(cfg$spec, "_")
    }
    tc <- paste0(prefix, seq_len(cfg$degree))
    lead <- copy(cw)[order(GID_0, year)]
    lead[, (paste0("lead_", tc)) := lapply(.SD, shift, n = 1, type = "lead"),
         by = GID_0, .SDcols = tc]
    lead <- lead[, c("GID_0", "year", paste0("lead_", tc)), with = FALSE]

    o <- copy(lead); setnames(o, c("GID_0", paste0("lead_", tc)),
                              c("orig", paste0("lead_", tc, "_orig")))
    d <- copy(lead); setnames(d, c("GID_0", paste0("lead_", tc)),
                              c("dest", paste0("lead_", tc, "_dest")))
    dt <- merge(copy(step2_data), o, by = c("orig", "year"), all.x = TRUE)
    dt <- merge(dt, d, by = c("dest", "year"), all.x = TRUE)

    ## The diff prefix carries the arm, as it does in prepare_diffs(). Two arms
    ## writing under one name is the failure the column guardrail of section 6
    ## exists to catch.
    lead_pre  <- if (is_annual) paste0("lead_diff_annual_", cfg$spec, "_")
                 else           paste0("lead_diff_", cfg$spec, "_")
    lead_cols <- paste0(lead_pre, seq_len(cfg$degree))
    for (k in seq_len(cfg$degree))
      dt[[lead_cols[k]]] <- (dt[[paste0("lead_", tc[k], "_dest")]] -
                             dt[[paste0("lead_", tc[k], "_orig")]]) / Omega
    dt[, fe_orig_dest := paste(orig, dest, sep = "_")]
    dt[, fe_year      := as.character(year)]

    lead_cols <- lead_cols[lead_cols %in% names(dt)]
    mod <- tryCatch(feglm(
      as.formula(paste0("explained ~ ",
                        paste(c(lead_cols, AMENITY_STAGE3_CONTROLS), collapse = " + "))),
      data = as.data.frame(dt), family = poisson(),
      cluster = c("orig", "dest"), fixef = c("fe_orig_dest", "fe_year")),
      error = function(e) { message(e$message); NULL })
    if (is.null(mod)) return(NULL)

    lin <- lead_cols[lead_cols %in% names(coef(mod))]
    wl  <- if (length(lin)) wald(mod, lin) else list(stat = NA, p = NA)
    data.table(spec = cfg$spec, deg = cfg$degree,
               w_lower = w$lower, w_upper = w$upper,
               F_lead = wl$stat, p_lead = wl$p,
               within_r2 = 1 - sum(resid(mod)^2) / (var(mod$y) * (length(mod$y) - 1)))
  }))))

  p_min <- min(rows$p_lead, na.rm = TRUE)
  arm_caption <- if (is_annual) "annual mean temperature"
                 else          "daily temperature distribution"
  ## Le bras journalier garde sa note MOT POUR MOT : la table journaliere doit
  ## se reemettre au bit, sinon on ne peut plus distinguer ce que la migration
  ## de l'API a change de ce qu'elle n'a pas touche.
  arm_metric  <- if (is_annual) " annual mean" else ""
  fam <- ifelse(rows$spec == "polynomial", "Polynomial", "Natural spline")
  flag <- ifelse(rows$p_lead > 0.10, "\\checkmark",
          ifelse(rows$p_lead > 0.05, "$\\sim$", "\\texttimes"))
  body <- sprintf("%s & %d & $[%s\\%%, %s\\%%]$ & %.3f & %.4f & %.4f & %s \\\\",
                  fam, rows$deg, rows$w_lower, rows$w_upper,
                  rows$F_lead, rows$p_lead, rows$within_r2, flag)

  writeLines(c(
    "%% GENERATED by 01b_amenity_drf.R -- DO NOT EDIT BY HAND.",
    "\\begin{table}[H]", "\\centering\\small",
    sprintf("\\caption{Placebo pre-trend test --- Migration (%s)}", arm_caption),
    sprintf("\\label{tab:exogeneity_mrt_migration%s}",
            if (is_annual) "_annual" else ""),
    "\\begin{tabular}{llccccc}", "\\hline\\hline",
    paste("Functional form & Degree & Winsorization & $F$-stat (leads) &",
          "$p$-value & Within-$R^2$ & No pre-trend \\\\"),
    "\\hline", body, "\\hline\\hline",
    "\\multicolumn{7}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    sprintf(paste0("\\textit{Note.} One-period-ahead bilateral%s temperature ",
                   "differentials substituted"), arm_metric),
    "for contemporaneous differentials in the Stage~3 PPML regression. Fixed effects:",
    "origin$\\times$destination and year; standard errors clustered by origin and",
    "destination. H$_0$: lead coefficients jointly insignificant.",
    ## LA NOTE CITE LA MACRO, ELLE NE RE-IMPRIME PAS LE NOMBRE.
    ## 17/09/2026 : le nombre imprime ici et la macro \EMPplaceboPmin
    ## etaient deux chemins vers la meme quantite, et ils avaient diverge --
    ## la note disait 0.5553 (min de CETTE table), la macro 0.212 (min sur
    ## les 24 cellules des DEUX bras), a deux centimetres l'une de l'autre
    ## dans D.2. La macro est desormais le minimum de cette table ; la note
    ## l'appelle, donc il n'y a plus deux chemins.
    sprintf("The smallest $p$-value on the leads is $%s$.",
            if (is_annual) "\\EMPplaceboPminAnnual" else "\\EMPplaceboPmin"),
    "$\\checkmark$ = $p>0.10$;\\ $\\sim$ = $0.05<p\\leq0.10$;\\ \\texttimes\\ = $p\\leq0.05$.}",
    "\\end{tabular}", "\\end{table}"), output_path)
  cat(sprintf("[write] %s — min p sur les avances = %.4f\n", basename(output_path), p_min))
  invisible(rows)
}

cat("=== Placebo (pré-tendance) ===\n")
placebo <- emit_placebo_table(
  file.path(DIR_OUT, "tab_exogeneity_mrt_migration.tex"),
  is_annual = FALSE)
placebo_annual <- emit_placebo_table(
  file.path(DIR_OUT, "tab_exogeneity_mrt_annual_migration.tex"),
  is_annual = TRUE)

## COLUMN guardrail, the same discipline as section 6. A SHAPE assertion sees
## NOTHING here: both arms return 12 rows by construction. Only a comparison of
## CONTENT, column against column, catches an arm that was in fact estimated
## twice on the same basis -- the failure that hid for weeks behind a file whose
## name said "withMRT" and whose bytes said "annual".
for (cc in c("F_lead", "p_lead", "within_r2")) {
  if (nrow(placebo_annual) == nrow(placebo) &&
      isTRUE(all.equal(placebo[[cc]], placebo_annual[[cc]])))
    stop(sprintf("emit_placebo_table : colonne %s IDENTIQUE entre les deux bras. ", cc),
         "Le bras annuel a ete estime sur la base journaliere.")
}
if (nrow(placebo_annual) != 12L ||
    !all(is.finite(placebo_annual$p_lead)) ||
    !all(is.finite(placebo_annual$within_r2)))
  stop(sprintf("emit_placebo_table : bras annuel non exploitable (%d lignes, ",
               nrow(placebo_annual)),
       sprintf("%d p non finies, %d within-R2 non finis). ",
               sum(!is.finite(placebo_annual$p_lead)),
               sum(!is.finite(placebo_annual$within_r2))),
       "EMPplaceboPmin porte sur les 24 cellules : ne pas emettre sur 12.")
cat(sprintf(paste("[ok] placebo : deux bras distincts, %d cellules chacun ;",
                  "min p journalier %.4f, annuel %.4f\n"),
            nrow(placebo), min(placebo$p_lead), min(placebo_annual$p_lead)))

##################################################################################
# 6. TABLES DE ROBUSTESSE
#
# COLUMN guardrail. The two tables were once copied from one another, and a
# SHAPE assertion would have seen nothing: both had 12 rows. Only a comparison of
# CONTENT, column against column, detects the duplication.
##################################################################################

emit_robustness_amenity <- function(jt, output_path, arm_label, twin = NULL) {
  if (!is.null(twin)) for (cc in c("fval_joint", "pval_joint", "fval_nonlin", "pval_nonlin"))
    if (nrow(twin) == nrow(jt) && isTRUE(all.equal(jt[[cc]], twin[[cc]])))
      stop(sprintf("emit_robustness_amenity : colonne %s IDENTIQUE entre les deux bras.", cc))

  m <- regmatches(jt$model_name,
                  regexec("^(polynomial|spline)([0-9])_wl([0-9.]+)wu([0-9.]+)$", jt$model_name))
  if (!all(vapply(m, length, 1L) == 5L))
    stop("emit_robustness_amenity : model_name non analysable.")
  fam <- vapply(m, function(x) if (x[2] == "polynomial") "Polynomial" else "Natural spline", "")
  deg <- vapply(m, `[`, "", 3); wlo <- vapply(m, `[`, "", 4); whi <- vapply(m, `[`, "", 5)
  ## Sort on spec, NOT on the display label: "Natural spline" sorts before
  ## "Polynomial" alphabetically, which would reverse the order used in every
  ## other table.
  ord <- order(jt$spec, deg, as.numeric(wlo))

  fmtp <- function(x) ifelse(is.na(x), "---",
                      ifelse(x < 0.001, "$<0.001$", sprintf("%.3f", x)))
  body <- sprintf("%s & %s & $[%s\\%%, %s\\%%]$ & %.3f & %s & %s \\\\",
                  fam[ord], deg[ord], wlo[ord], whi[ord], jt$fval_joint[ord],
                  fmtp(jt$pval_joint[ord]), fmtp(jt$pval_nonlin[ord]))

  ## Two MINIMA, not a threshold count: a minimum names its object and lets the
  ## lecteur appliquer son propre seuil.
  writeLines(c(
    "%% GENERATED by 01b_amenity_drf.R :: emit_robustness_amenity()",
    "%% Source: amenity_joint_tests*.csv. DO NOT EDIT BY HAND.",
    "\\begin{table}[H]", "\\centering\\small",
    sprintf("\\caption{Robustness: Amenity dose-response function (%s)}", arm_label),
    sprintf("\\label{tab:robustness_amenity%s}",
            if (grepl("annual", arm_label, ignore.case = TRUE)) "_annual" else ""),
    "\\begin{tabular}{llcccc}", "\\toprule",
    "Functional form & Deg. & Winsorization & $F$ (joint) & $p$ (joint) & $p$ (non-lin.) \\\\",
    "\\midrule", body, "\\bottomrule",
    "\\multicolumn{6}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    "\\textit{Note.} Stage~3 PPML estimates on bilateral migration flows, with",
    "origin$\\times$destination and year fixed effects. Standard errors clustered by",
    "origin and destination. $F$~(joint) and $p$~(joint): joint significance of all",
    "climate terms; $p$~(non-lin.): joint significance of the degree-$\\geq 2$ terms.",
    sprintf(paste0("Across the %d specifications the smallest joint $p$-value is %.3f; ",
                   "the smallest $p$-value for the degree-$\\geq 2$ terms alone is %.3f."),
            nrow(jt), min(jt$pval_joint, na.rm = TRUE), min(jt$pval_nonlin, na.rm = TRUE)),
    "No within-pair $R^2$ is reported: it is not computed on this specification grid,",
    "and the within-pair $R^2$ quoted elsewhere in this appendix belongs to the",
    "placebo specification. No specification is singled out.}",
    "\\end{tabular}", "\\end{table}"), output_path, useBytes = TRUE)
  cat(sprintf("[write] %s\n", basename(output_path)))
}

emit_robustness_amenity(jt_daily,
  file.path(DIR_OUT, "tab_robustness_amenity_formatted.tex"),
  "daily temperature distribution")
emit_robustness_amenity(jt_annual,
  file.path(DIR_OUT, "tab_robustness_amenity_annual_formatted.tex"),
  "annual mean temperature", twin = jt_daily)

##################################################################################
# 7. MONTE CARLO -- TWO LAYERS OF UNCERTAINTY
#   layer 1: the specification is drawn according to the Akaike weights;
#   couche 2 : beta* ~ N(beta_hat, V_hat).
# Output: the quantiles the GE consumes.
##################################################################################

mc_daily  <- mc_run(results_daily,  climate, MC_N_SIM, MC_X_GRID,
                    MC_SEED[["amenity_daily"]],  FALSE, "Amenity-Daily")
mc_annual <- mc_run(results_annual, climate, MC_N_SIM, MC_X_GRID,
                    MC_SEED[["amenity_annual"]], TRUE,  "Amenity-Annual")
s_daily   <- mc_summarise(mc_daily)
s_annual  <- mc_summarise(mc_annual)

mc_export(s_daily,  file.path(DIR_OUT, "MC_damage_amenity_daily_quantiles.csv"))           
mc_export(s_annual, file.path(DIR_OUT, "MC_damage_amenity_annual_quantiles.csv")) 

## THE g = 0 TABLES, the Manski lower bound on the amenity channel.
##
## They carry no estimate: the amenity response is set to zero everywhere, on the
## same grid and with the same columns as the tables above, so that the GE can be
## re-solved with the channel switched off. The comparison answers one question --
## does the Northern reversal survive g = 0? -- and the two sensitivity tables of
## the manuscript are built from it.
##
## Written here rather than shipped as a fixed file: a constant function needs no
## provenance beyond the grid it is defined on, and a GE input with no producer in
## the package cannot be regenerated by a reader.
null_tab <- data.table(temperature_C = MC_EXPORT_GRID)
for (q in c("q05", "q25", "q50", "q75", "q95", "y_mean")) null_tab[[q]] <- 0
for (fr in c("daily", "annual")) {
  p_null <- file.path(DIR_OUT, sprintf("MC_damage_amenity_%s_quantiles_NULL.csv", fr))
  fwrite(null_tab, p_null)
  cat(sprintf("[write] %s | g = 0 sur %d points\n", basename(p_null), nrow(null_tab)))
}


##################################################################################
# 10. CENTRAL DAMAGE FUNCTIONS -- inputs to 01e_diagnostics.R
#
# Curves averaged with the Akaike weights, on the grid common to all
# specifications. 01e consumes them to build the damage indices.
##################################################################################

for (arm in c("daily", "annual")) {
  res <- if (arm == "daily") results_daily else results_annual
  avg <- aic_curve(res, climate, is_annual = (arm == "annual"))
  setnames(avg, c("x", "y_avg"), c("temperature_C", "zeta_a"))
  f <- file.path(DIR_OUT, sprintf("damage_amenity%s_bin.csv",
                                  if (arm == "daily") "" else "_annual"))
  fwrite(avg, f)
  cat(sprintf("[write] %s (%d points, [%.1f, %.1f] °C)\n",
              basename(f), nrow(avg), min(avg$temperature_C),
              max(avg$temperature_C)))
}

##################################################################################
# 8. PAPER FIGURES -- the two fan charts
#
# Redrawn from the quantile CSV files, with no re-estimation. Three properties:
#   - the support is DERIVED from the artefact by plateau detection, never coded;
#   - the two panels share a window and a scale, so they can be compared by
#     moving the eye;
#   - the shift is SCALAR: the median is zero at the reference temperature, but
#     the band is not pinched there, because each draw was normalised to its own
#     maximum. The axis therefore says "MC median = 0", never "normalised".
##################################################################################

plot_fans("amenity", climate)  

##################################################################################
# 9. INPUTS OF THE POWER CALCULATION (identification appendix)
#
# Emits sigma_eps under three candidate scales, and N_eff, the size of the pair
# fixed effect, for the minimum-detectable-effect macros. The choice of scale
# remains an editorial one and is not made here.
##################################################################################

aic_mde <- aic_weights(results_daily)
mde <- rbindlist(lapply(aic_mde$model_name, function(tag) {
  m  <- results_daily[[tag]]$model
  rr <- tryCatch(residuals(m, type = "response"), error = function(e) NULL)
  if (is.null(rr)) return(NULL)
  mu <- fitted(m); y <- mu + rr; pos <- y > 0 & mu > 0
  data.table(tag = tag, weight = aic_mde$weight[aic_mde$model_name == tag],
             nobs = m$nobs,
             sd_deviance = sd(residuals(m, type = "deviance"), na.rm = TRUE),
             sd_pearson  = sd(residuals(m, type = "pearson"),  na.rm = TRUE),
             sd_logratio = sd(log(y[pos] / mu[pos]), na.rm = TRUE))
}))
fwrite(mde, file.path(DIR_OUT, "mde_amenity_inputs.csv"))

best <- aic_mde$model_name[which.max(aic_mde$weight)]
fs <- results_daily[[best]]$model$fixef_sizes
fwrite(data.table(fixef = names(fs), size = as.integer(fs), best_tag = best),
       file.path(DIR_OUT, "mde_amenity_fixef_sizes.csv"))
cat(sprintf("[write] mde_amenity_{inputs,fixef_sizes}.csv | spec poids-max : %s\n", best))

cat("\n=== 03b_amenity_drf.R terminé ===\n")
