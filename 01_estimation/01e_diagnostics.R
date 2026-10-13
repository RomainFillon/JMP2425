##################################################################################
# 01e_diagnostics.R -- non-stationarity and the common share of damages.
#
# WHAT THIS STEP MEASURES
#   eta^j  the share of the temporal variance of damages that a year fixed
#          effect absorbs, that is, what is COMMON to every country. The lower
#          eta, the less a model with a single world temperature suffices.
#   rho_r  the correlation between a country damage index and world temperature,
#          historical and then projected. Their GAP is the non-stationarity.
#   Dbar   the dispersion across countries around the world mean, at the three
#          horizons of the paper.
#
# THE ORDER OF THE SECTIONS IS NOT THE ORDER OF THE ARGUMENT, DELIBERATELY.
# The common sample -- the countries that have BOTH an ERA baseline and a CMIP6
# projection -- must be known BEFORE the weighting is fixed. Otherwise eta covers
# one set of countries and Dbar another, and the equality that the section 10
# control ASSERTS between its equal-weighted arm and `res` becomes false.
# Hence: load -> indices -> future -> RESTRICT -> weight -> eta.
#
# The three merged scripts shared a convention without sharing code: equal
# weights across countries. One of them declared itself consistent with another
# through a COMMENT. Here the weight vector is defined once, on a single sample,
# and all three sections read it.
#
# Reads:  output/damage_{amenity,prod}_{,annual_}bin.csv, from 01a and 01b
#         output/drf_prod_specs/drf_prod_spec_*.csv, from 01a (section 8bis)
#         intermediate_outputs/ERA_tas_annual.csv, CMIP6_tas.csv,
#         poids_par_pays.csv, poids_par_cellule.csv (lot 51)
# Writes:
#   (.tex  AUCUN depuis le 11/09/2026 : tab_dbar_cross.tex retiree avec la
#          section F du manuscrit. Le calcul `res` demeure : il sert au
#          controle positif de la section 10.)
#   .pdf   map_rho_hist_{amenity,prod}.pdf,
#          map_nonstationarity_{amenity,prod}.pdf
#   .csv   table_eta.csv, table_eta_quinquennial.csv, rho_hist_by_country.csv,
#          rho_combined_by_country.csv, damage_projection_country.csv,
#          pe_signtest_by_country.csv, jensen_gap_country.csv
#          (dbar_cross_weighting.csv n'est plus produit : section 10 retiree
#           le 28/09/2026, lot 54 -- sans consommateur depuis le 11/09)
# Needs:  01a and 01b must have run.

##################################################################################

suppressMessages({
  library(data.table); library(ggplot2); library(scales)
  library(rnaturalearth); library(sf); library(countrycode)
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
    stop("Racine indeterminable. Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(.self))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))



HORIZONS   <- c(2030L, 2050L, 2095L)
BASE_YEARS <- 2011:2015
## Year of the population weights used to aggregate CMIP6 to the country.
## 2015 since lot 51 (28/09/2026) : the calibration year of the model, used
## everywhere a country is attributed or a quantity distributed (author's decision, lot 46).
## It was 1995, "the median of the seven historical periods".
LC_YEAR <- 2015L

t0 <- proc.time()

##################################################################################
# 1. FONCTIONS DE DOMMAGE ET CLIMAT HISTORIQUE
##################################################################################

cat("=== 1. Chargement ===\n")

read_bin <- function(nm, newname) {
  f <- file.path(DIR_OUT, nm)
  if (!file.exists(f)) stop("[STOP] absent : ", f, "\n  Produit par 01a/01b.")
  d <- fread(f); setnames(d, c("temp", newname)); d
}
dam_amenity <- merge(read_bin("damage_amenity_bin.csv",        "zeta_a_daily"),
                     read_bin("damage_amenity_annual_bin.csv", "zeta_a_annual"),
                     by = "temp", all = TRUE)
dam_prod    <- merge(read_bin("damage_prod_bin.csv",           "zeta_z_daily"),
                     read_bin("damage_prod_annual_bin.csv",    "zeta_z_annual"),
                     by = "temp", all = TRUE)

## Extrapolation by the edge value (rule = 2), never linear beyond the support:
## a straight line extended outside the estimated domain invents a slope.
interp_damage <- function(x, tbl, col)
  approx(tbl$temp, tbl[[col]], xout = x, method = "linear", rule = 2)$y

hist_clim <- fread(file.path(DIR_TMP, "ERA_tas_annual.csv"))
if ("value_bin" %in% names(hist_clim)) setnames(hist_clim, "value_bin", "value")
hist_clim[, share := count / sum(count, na.rm = TRUE), by = .(GID_0, year)]
cat(sprintf("  climat : %d lignes | %d pays | %d annees\n",
            nrow(hist_clim), uniqueN(hist_clim$GID_0), uniqueN(hist_clim$year)))

##################################################################################
# 2. INDICES DE DOMMAGE
#
# NO WEIGHTING HERE. It waits for section 4bis, which fixes the sample.
##################################################################################

cat("=== 2. Indices de dommage ===\n")

for (v in c("zeta_a_daily", "zeta_a_annual"))
  hist_clim[, (v) := interp_damage(value, dam_amenity, v)]
for (v in c("zeta_z_daily", "zeta_z_annual"))
  hist_clim[, (v) := interp_damage(value, dam_prod, v)]

damage_idx <- hist_clim[, .(
  D_a_daily  = sum(zeta_a_daily  * share, na.rm = TRUE),
  D_a_annual = sum(zeta_a_annual * share, na.rm = TRUE),
  D_z_daily  = sum(zeta_z_daily  * share, na.rm = TRUE),
  D_z_annual = sum(zeta_z_annual * share, na.rm = TRUE),
  T_mean     = sum(value * share, na.rm = TRUE)), by = .(GID_0, year)]
cat(sprintf("  %d pays avec base ERA (echantillon PROVISOIRE)\n",
            uniqueN(damage_idx$GID_0)))

##################################################################################
# 3. CMIP6 PROJECTION -- aggregation from pixel to country
##################################################################################

cat("=== 3. Projection CMIP6 ===\n")

fut_file <- file.path(DIR_TMP, "CMIP6_tas.csv")
## Lot 51 : poids_par_pays.csv (00d_poids_population.py) remplace lengths_country.csv,
## qui comptait la population ENTIERE de chaque cellule touchee (RECAP 3.231).
lc_file  <- file.path(DIR_TMP, "poids_par_pays.csv")
HAS_FUTURE <- file.exists(fut_file) && file.exists(lc_file)

if (!HAS_FUTURE) {
  ## A mechanism that can do nothing must say whether it acted.
  cat("  [OMIS] partie projetee NON CALCULEE :\n")
  if (!file.exists(fut_file)) cat("    CMIP6_tas.csv absent\n")
  if (!file.exists(lc_file))  cat("    poids_par_pays.csv absent (00d)\n")
  cat("    -> map_nonstationarity_*.pdf NON produites.\n")
} else {
  fut <- fread(fut_file)
  fut[, share := count / sum(count, na.rm = TRUE), by = .(x, y, year)]
  lc <- fread(lc_file, na.strings = NULL)[year == LC_YEAR & GID_0 != "",
                                          .(x, y, GID_0, part_du_pays)]
  setkey(lc, x, y)
  ## Realigned on the one-degree cell centres: CMIP6 and lengths_country do not
  ## use the same coordinate convention.
  fut[, `:=`(x_snap = floor(x) + 0.5, y_snap = floor(y) + 0.5)]
  fj <- merge(fut, lc, by.x = c("x_snap", "y_snap"), by.y = c("x", "y"),
              all.x = FALSE, allow.cartesian = TRUE)
  if (!nrow(fj)) stop("[STOP] aucun pixel CMIP6 apparie a poids_par_pays.")
  ## Renormalised WITHIN the country over the matched pixels: part_du_pays is the
  ## share of the country's population living in the cell (sums to 1 per country).
  fj[, tot := sum(part_du_pays, na.rm = TRUE), by = .(GID_0, year)]
  fj <- fj[tot > 0][, w := part_du_pays / tot]
  fut_clim <- fj[, .(count = sum(share * w, na.rm = TRUE)), by = .(GID_0, year, value)]
  fut_clim[, share := count / sum(count, na.rm = TRUE), by = .(GID_0, year)]
  rm(fut, fj); gc(verbose = FALSE)

  for (v in c("zeta_a_daily", "zeta_a_annual"))
    fut_clim[, (v) := interp_damage(value, dam_amenity, v)]
  for (v in c("zeta_z_daily", "zeta_z_annual"))
    fut_clim[, (v) := interp_damage(value, dam_prod, v)]
  fut_damage <- fut_clim[, .(
    D_a_daily = sum(zeta_a_daily * share, na.rm = TRUE),
    D_z_daily = sum(zeta_z_daily * share, na.rm = TRUE),
    T_mean    = sum(value * share, na.rm = TRUE)), by = .(GID_0, year)]
  cat(sprintf("  %d pays projetes | %d-%d\n", uniqueN(fut_damage$GID_0),
              min(fut_damage$year), max(fut_damage$year)))
}

##################################################################################
# 4. ECHANTILLON COMMUN ET PONDERATION
#
# THESE TWO BLOCKS ARE THE COHERENCE CONDITION OF THIS FILE. Some countries have
# an ERA baseline without a CMIP6 projection, because the climate grid does not
# resolve them, and some the reverse. Computing eta on one set and Dbar on the
# other would break the equality that the section 10 control asserts between
# its equal-weighted arm and `res`, which is now what carries this condition.
#
# EQUAL WEIGHTS ACROSS COUNTRIES, omega_r = 1/n_c. This is NOT a population
# weighting, and the section title used to say "equal population weights"; that
# mislabel propagated as far as the manuscript.
#
# THE AGGREGATION HAS TWO STAGES, and the phrase "pop-weighted" is correct only
# at the first:
#   cells   -> country : POPULATION-WEIGHTED (ERA_tas_annual.csv already carries
#                        population shares; CMIP6 is aggregated in section 3);
#   country -> world   : EQUALLY WEIGHTED (the omega line below).
# Moving to genuine population weights would require doing so in EVERY section of
# this file, and would change the published figures.
##################################################################################

cat("=== 4. Echantillon commun et ponderation ===\n")

if (HAS_FUTURE) {
  common    <- intersect(unique(damage_idx$GID_0), unique(fut_damage$GID_0))
  drop_hist <- setdiff(unique(damage_idx$GID_0), common)
  drop_fut  <- setdiff(unique(fut_damage$GID_0), common)
  if (length(drop_hist))
    cat(sprintf("  [!] %d pays avec base ERA sans projection CMIP6, ecartes :\n      %s\n",
                length(drop_hist), paste(drop_hist, collapse = ", ")))
  if (length(drop_fut))
    cat(sprintf("  [!] %d pays projetes sans base ERA, ecartes\n", length(drop_fut)))
  damage_idx <- damage_idx[GID_0 %in% common]
  fut_damage <- fut_damage[GID_0 %in% common]
  fut_clim   <- fut_clim[GID_0 %in% common]
} else {
  cat("  [!] pas de projection : eta et rho portent sur l'echantillon ERA seul.\n")
}

N_COUNTRIES <- uniqueN(damage_idx$GID_0)
OMEGA <- 1 / N_COUNTRIES            # the convention, set ONCE
damage_idx[, omega := OMEGA]
cat(sprintf("  %d pays | omega = 1/%d | eta, rho et Dbar sur CET ensemble\n",
            N_COUNTRIES, N_COUNTRIES))

damage_idx <- merge(damage_idx,
                    damage_idx[, .(T_global = sum(omega * T_mean) / sum(omega)),
                               by = year], by = "year")

##################################################################################
# 5. ETA -- the common share of the time variance
#
#   eta = Var_t(D_world_t) / Var_{r,t}(D_rt), on within-country demeaned series.
# without them one cannot say whether a lower eta comes from a smaller NUMERATOR
# (less common movement) or a larger DENOMINATOR (more idiosyncratic
# dispersion). The two readings are different, and the ratio alone does not
# separate them.
##################################################################################

cat("=== 5. Eta ===\n")

compute_eta <- function(dt, col) {
  dt <- copy(dt)
  dt[, dm := get(col) - mean(get(col), na.rm = TRUE), by = GID_0]
  w  <- dt[, .(D_world_dm = sum(omega * dm) / sum(omega)), by = year]
  vw <- var(w$D_world_dm, na.rm = TRUE); vt <- var(dt$dm, na.rm = TRUE)
  list(eta = vw / vt, var_world = vw, var_total = vt)
}

## THE ROW ORDER IS STRUCTURAL: the emitter indexes by POSITION, with the
## benchmark on row one and the spacing after rows one and three. Reordering
## here would move the note onto the wrong row there.
ETA_SPECS <- list(
  list(metric = "Annual mean T",              channel = "---",          col = "T_mean"),
  list(metric = "Damage index (daily bins)",  channel = "Amenity",      col = "D_a_daily"),
  list(metric = "Damage index (annual mean)", channel = "Amenity",      col = "D_a_annual"),
  list(metric = "Damage index (daily bins)",  channel = "Productivity", col = "D_z_daily"),
  list(metric = "Damage index (annual mean)", channel = "Productivity", col = "D_z_annual"))

eta_results <- rbindlist(lapply(ETA_SPECS, function(s) {
  r <- compute_eta(damage_idx, s$col)
  data.table(metric = s$metric, channel = s$channel, eta = r$eta,
             var_world = r$var_world, var_total = r$var_total)
}))
eta_results[, eta_pct := round(100 * eta, 1)]
eta_results[, n_countries := N_COUNTRIES]
fwrite(eta_results, file.path(DIR_OUT, "table_eta.csv"))
print(eta_results[, .(metric, channel, eta = round(eta, 3), eta_pct)])

## Eta on quinquennial aggregates, the clause of the table note. This value was
## PRINTED and never written, so it was a published number nobody could check.
## It is now persisted.
{
  q <- copy(damage_idx)
  q[, per := (year %/% 5L) * 5L]
  q5 <- q[, .(T_mean = mean(T_mean), omega = OMEGA), by = .(GID_0, per)]
  setnames(q5, "per", "year")
  r <- compute_eta(q5, "T_mean")
  fwrite(data.table(metric = "Annual mean T (quinquennial)", eta = r$eta,
                    var_world = r$var_world, var_total = r$var_total,
                    n_countries = N_COUNTRIES),
         file.path(DIR_OUT, "table_eta_quinquennial.csv"))
  cat(sprintf("  eta quinquennal (T moyen) : %.3f\n", r$eta))
}

##################################################################################
# 6. RHO -- correlation historique, projetee, et non-stationnarite
##################################################################################

cat("=== 6. Rho ===\n")

rho_hist <- damage_idx[, .(
  rho_a_daily  = cor(D_a_daily,  T_global, use = "complete.obs"),
  rho_a_annual = cor(D_a_annual, T_global, use = "complete.obs"),
  rho_z_daily  = cor(D_z_daily,  T_global, use = "complete.obs"),
  rho_z_annual = cor(D_z_annual, T_global, use = "complete.obs"),
  mean_T       = mean(T_mean, na.rm = TRUE)), by = GID_0]
fwrite(rho_hist, file.path(DIR_OUT, "rho_hist_by_country.csv"))
cat(sprintf("  historique : %d pays | rho_a_daily median %.3f\n",
            nrow(rho_hist), median(rho_hist$rho_a_daily, na.rm = TRUE)))

if (HAS_FUTURE) {
  fut_damage <- merge(fut_damage,
    fut_damage[, .(T_global = mean(T_mean, na.rm = TRUE)), by = year], by = "year")
  rho_proj <- fut_damage[, .(
    rho_a_proj = cor(D_a_daily, T_global, use = "complete.obs"),
    rho_z_proj = cor(D_z_daily, T_global, use = "complete.obs")), by = GID_0]
  rho_combined <- merge(rho_hist, rho_proj, by = "GID_0", all.x = TRUE)
  rho_combined[, `:=`(delta_rho_a = rho_a_proj - rho_a_daily,
                      delta_rho_z = rho_z_proj - rho_z_daily)]
  fwrite(rho_combined, file.path(DIR_OUT, "rho_combined_by_country.csv"))
  cat(sprintf("  projete : %d pays | delta_rho_a median %.3f | mean %.3f\n",
              nrow(rho_proj), median(rho_combined$delta_rho_a, na.rm = TRUE),
              mean(rho_combined$delta_rho_a, na.rm = TRUE)))
} else rho_combined <- copy(rho_hist)

##################################################################################
# 7. CARTES
##################################################################################

cat("=== 7. Cartes ===\n")

world_sf <- ne_countries(scale = "medium", returnclass = "sf")
## iso_a3 codes France and Norway as "-99"; iso_a3_eh resolves them.
ic <- c("iso_a3_eh", "adm0_a3", "iso_a3")
ic <- ic[ic %in% names(world_sf)][1]
world_sf$iso3c_join <- world_sf[[ic]]
world_sf <- world_sf[world_sf$iso3c_join != "-99" & !is.na(world_sf$iso3c_join), ]

make_rho_map <- function(dt, col, title, subtitle, limits = c(-1, 1)) {
  m <- merge(world_sf, dt[, .SD, .SDcols = c("GID_0", col)],
             by.x = "iso3c_join", by.y = "GID_0", all.x = TRUE)
  setnames(m, col, "rho_val")
  ggplot(m) +
    geom_sf(aes(fill = rho_val), colour = "white", linewidth = 0.07) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, limits = limits, na.value = "grey80",
      labels = scales::number_format(accuracy = 0.1),
      guide = guide_colourbar(title = "Correlation", title.position = "top",
                              barwidth = 12, barheight = 0.7)) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
          plot.subtitle = element_text(size = 8.5, colour = "grey40", hjust = 0,
                                       margin = margin(t = 3, b = 6)),
          legend.position = "bottom", legend.title = element_text(size = 8.5),
          legend.text = element_text(size = 8), plot.margin = margin(8, 8, 8, 8))
}

save_map <- function(p, nm) {
  ggsave(file.path(DIR_OUT, nm), plot = p, width = 11, height = 6)
  cat(sprintf("  [write] %s\n", nm))
}

save_map(make_rho_map(rho_hist, "rho_a_daily",
  "Historical correlation: amenity damage index vs global mean temperature",
  paste("rho^hist_r = corr(D^a_t(r), T^global_t), ERA5. Red: damage tracks the",
        "global trend. Blue: idiosyncratic local dynamics.")),
  "map_rho_hist_amenity.pdf")
save_map(make_rho_map(rho_hist, "rho_z_daily",
  "Historical correlation: productivity damage index vs global mean temperature",
  "rho^hist_r = corr(D^z_t(r), T^global_t), ERA5."),
  "map_rho_hist_prod.pdf")

if (HAS_FUTURE) {
  save_map(make_rho_map(rho_combined, "delta_rho_a",
    "Non-stationarity: change in correlation between historical and projected periods",
    paste("rho^proj_r - rho^hist_r, amenity. Red: damages increasingly track the",
          "global trend. Blue: they diverge from it."), limits = c(-0.5, 0.5)),
    "map_nonstationarity_amenity.pdf")
  save_map(make_rho_map(rho_combined, "delta_rho_z",
    "Non-stationarity: change in correlation (productivity damage index)",
    "rho^proj_r - rho^hist_r, productivity.", limits = c(-0.5, 0.5)),
    "map_nonstationarity_prod.pdf")
}

##################################################################################
# 8. PROJECTED DAMAGES BY COUNTRY
#
# Carried over from the damage projection script, reduced to what the
# manuscript consumes.
# ITS CRH AMENITY BRANCH IS NOT REUSED: it applies the logistic of the CRH
# design, which is a different object. Mixing the two would change what the
# damage index denotes, and a published macro descends from it.
# Its three figures are not reused either: the manuscript loads none of them.
##################################################################################

cat("=== 8. Projection par pays ===\n")

if (!HAS_FUTURE) {
  cat("  [OMIS] pas de projection -> damage_projection_country.csv NON produit.\n")
} else {
  ## Baseline: the ERA mean of the last five historical years. The deltas are
  ## measured against it, not against an isolated year.
  base <- damage_idx[year %in% BASE_YEARS,
                     .(D_a_base = mean(D_a_daily, na.rm = TRUE),
                       D_z_base = mean(D_z_daily, na.rm = TRUE),
                       T_base   = mean(T_mean,    na.rm = TRUE)), by = GID_0]
  cat(sprintf("  base ERA %d-%d : %d pays\n",
              min(BASE_YEARS), max(BASE_YEARS), nrow(base)))

  ## The sample was restricted in section 4, so this merge must lose nothing.
  proj <- merge(fut_damage, base, by = "GID_0")
  if (uniqueN(proj$GID_0) != N_COUNTRIES)
    stop(sprintf("[STOP] %d pays apres jointure contre %d attendus. ",
                 uniqueN(proj$GID_0), N_COUNTRIES),
         "La restriction de la section 4 n'a pas tenu.")

  proj[, `:=`(delta_a = D_a_daily - D_a_base,
              delta_z = D_z_daily - D_z_base,
              delta_T = T_mean    - T_base)]
  proj[, delta_total := delta_a + delta_z]
  setnames(proj, c("D_a_daily", "D_z_daily"), c("D_a", "D_z"))

  ## ---- COMPTES DE JOURS PAR SEUIL -------------------------------------------
  ## Ces trois colonnes portent six des neuf macros \EMPfin* du paragraphe
  ## Finlande de la Section 5. Elles etaient produites par le damage_projection.r
  ## HORS PAQUET et par lui seul : la section 8 les avait perdues au passage, si
  ## bien que l'emetteur s'abstenait sur « colonnes absentes » et que les NEUF
  ## macros restaient non definies -- y compris les deux qui ne dependent pas
  ## d'elles, la porte etant commune.
  ##
  ## CALCULEES SUR fut_clim, c'est-a-dire sur la MEME distribution CMIP6
  ## pop-ponderee que T_mean : meme intrant, meme ponderation, meme annee. La
  ## base n'est PAS la baseline ERA 2011-2015 des deltas de dommage, parce que
  ## le texte decrit la moyenne CMIP6 de l'annee (« rises from ... in 2015 »).
  ## `share` est normalisee a 1 par (GID_0, year), donc x 365 = jours par an.
  days_thr <- fut_clim[, .(
    days_ge20  = 365 * sum(share[value >=  20], na.rm = TRUE),
    days_ge25  = 365 * sum(share[value >=  25], na.rm = TRUE),
    days_ltm10 = 365 * sum(share[value <  -10], na.rm = TRUE)), by = .(GID_0, year)]
  proj <- merge(proj, days_thr, by = c("GID_0", "year"), all.x = TRUE)
  ## La jointure ne doit rien perdre : days_thr et proj sortent tous deux de
  ## fut_clim, sur l'echantillon commun deja restreint en section 4.
  if (anyNA(proj$days_ge20))
    stop(sprintf("[STOP] %d lignes sans compte de jours apres jointure. ",
                 sum(is.na(proj$days_ge20))),
         "days_thr et proj descendent pourtant du meme fut_clim.")

  fwrite(proj[, .(GID_0, year, T_mean, delta_T, delta_a, delta_z,
                  delta_total, D_a, D_z, days_ge20, days_ge25, days_ltm10)],
         file.path(DIR_OUT, "damage_projection_country.csv"))
  cat(sprintf("  [write] damage_projection_country.csv (%d lignes, %d pays)\n",
              nrow(proj), uniqueN(proj$GID_0)))

  ################################################################################
  # 8bis. LE PROXY D'EQUILIBRE PARTIEL DU TEST DE SIGNE
  #
  # PRODUIT output/pe_signtest_by_country.csv, qui porte les neuf macros
  # \EMPfin* du paragraphe Finlande de la Section 5.
  #
  # POURQUOI ICI, ET PAS DANS UN SCRIPT A PART (09/09/2026).
  # Le producteur hors paquet, proxy_pe_signtest.r, refaisait POUR SON COMPTE
  # l'assemblage climatique : relecture de CMIP6_tas.csv (462 Mo), re-jointure a
  # lengths_country.csv, re-ponderation par la population, re-calage des
  # coordonnees sur floor(x)+0.5. C'est exactement l'assemblage des sections 1
  # et 3 ci-dessus. Deux implantations du meme calcul derivent, et derivent en
  # silence : on ne le verrait qu'a l'ecart entre deux chiffres publies. La
  # section reutilise donc `hist_clim` et `fut_clim` -- une seule source
  # climatique pour tout le fichier.
  #
  # CE QUE LE PROXY CALCULE. Pour chaque pays froid, chaque specification
  # journaliere et chaque horizon :
  #
  #   Daily_extra_PE = [D^daily(proj) - D^daily(base)]  -  [D^annual(proj) - D^annual(base)]
  #
  # La composante ANNUELLE passe par la DRF annuelle ponderee FIXE, evaluee a la
  # moyenne du pays ; la composante JOURNALIERE integre la courbe de la
  # specification contre la distribution entiere. L'ecart entre les deux est
  # l'objet que le papier rapporte.
  #
  # DEUX BRAS DE BASE, DECLARES, COTE A COTE. Le bras ERA (2011-2015) est celui
  # des deltas de dommage de la section 8 ; le bras CMIP6 (2015) est celui que
  # le manuscrit decrit a l'exemple finlandais (« rises from ... in 2015 »).
  # L'emetteur LIT la colonne base_source_cmip et refuse d'emettre si elle ne
  # vaut pas l'annee attendue : le bras n'est jamais deduit du nom du fichier.
  ################################################################################

  cat("=== 8bis. Proxy PE du test de signe ===\n")

  ## Les cinq economies froides de l'illustration. Le test de signe porte sur
  ## elles parce que c'est la que le renversement se joue.
  PE_COLD  <- c("RUS", "CAN", "FIN", "DEU", "FRA")
  ## Annee de base du bras CMIP6. C'est celle que l'emetteur exige.
  PE_BASE_CMIP <- 2015L

  spec_dir <- file.path(DIR_OUT, "drf_prod_specs")
  spec_f   <- if (dir.exists(spec_dir))
    sort(list.files(spec_dir, pattern = "^drf_prod_spec_\\d+\\.csv$",
                    full.names = TRUE)) else character(0)

  if (!length(spec_f)) {
    ## Un mecanisme qui ne fait rien doit dire qu'il n'a rien fait.
    cat("  [OMIS] output/drf_prod_specs/ vide ou absent",
        " -> pe_signtest_by_country.csv NON produit,\n",
        "         et les 9 macros \\EMPfin* resteront non definies.\n", sep = "")
    cat("         Producteur : 01a_productivity_drf.R, section 9bis.\n")
  } else {
    curves <- rbindlist(lapply(spec_f, fread))
    ## Les poids voyagent AVEC les courbes. L'emetteur controle qu'ils somment
    ## a un ; le controler ici aussi ferait porter le refus par le producteur.
    wsum <- curves[, .(w = aic_weight[1]), by = model_name][, sum(w)]
    if (abs(wsum - 1) > 1e-9)
      stop(sprintf("[STOP] poids AIC des dumps de somme %.9f, 1 attendue. ", wsum),
           "Les douze courbes ne viennent pas du meme appel a aic_weights().")
    cat(sprintf("  %d specifications | poids AIC de somme %.6f\n",
                uniqueN(curves$model_name), wsum))

    manquants <- setdiff(PE_COLD, unique(fut_clim$GID_0))
    if (length(manquants))
      cat(sprintf("  [!] %d pays froids sans projection CMIP6, ecartes : %s\n",
                  length(manquants), paste(manquants, collapse = ", ")))

    ## Interpolation par la valeur de bord hors support (rule = 2), jamais
    ## lineaire : prolonger une droite au-dela du domaine estime inventerait
    ## une pente. C'est la convention de interp_damage ci-dessus.
    integ_daily <- function(dist, cv)
      sum(approx(cv$temperature_C, cv$f_z_raw, xout = dist$value,
                 method = "linear", rule = 2)$y * dist$share, na.rm = TRUE)
    annual_at <- function(Tm) interp_damage(Tm, dam_prod, "zeta_z_annual")

    ## Base ERA : moyenne des parts sur BASE_YEARS, renormalisee. Meme
    ## convention que la base des deltas de la section 8.
    era_base <- hist_clim[GID_0 %in% PE_COLD & year %in% BASE_YEARS,
                          .(share = mean(share, na.rm = TRUE)), by = .(GID_0, value)]
    era_base[, share := share / sum(share, na.rm = TRUE), by = GID_0]

    pe_rows <- list()
    for (cty in intersect(PE_COLD, unique(fut_clim$GID_0))) {
      base_e <- era_base[GID_0 == cty, .(value, share)]
      if (!nrow(base_e)) {
        cat(sprintf("  [!] %s : aucune base ERA, ecarte\n", cty)); next
      }
      Tm_base_e   <- sum(base_e$value * base_e$share, na.rm = TRUE)
      Dz_ann_e    <- annual_at(Tm_base_e)

      base_c <- fut_clim[GID_0 == cty & year == PE_BASE_CMIP, .(value, share)]
      has_c  <- nrow(base_c) > 0L
      Tm_base_c <- if (has_c) sum(base_c$value * base_c$share, na.rm = TRUE) else NA_real_
      Dz_ann_c  <- if (has_c) annual_at(Tm_base_c) else NA_real_
      if (!has_c)
        cat(sprintf("  [!] %s : pas d'annee %d dans CMIP6, bras CMIP6 a NA\n",
                    cty, PE_BASE_CMIP))

      for (h in HORIZONS) {
        pj <- fut_clim[GID_0 == cty & year == h, .(value, share)]
        if (!nrow(pj)) next
        Tm_proj     <- sum(pj$value * pj$share, na.rm = TRUE)
        Dz_ann_proj <- annual_at(Tm_proj)
        for (nm in unique(curves$model_name)) {
          cv <- curves[model_name == nm]
          d_base_e <- integ_daily(base_e, cv)
          d_proj   <- integ_daily(pj, cv)
          d_base_c <- if (has_c) integ_daily(base_c, cv) else NA_real_
          pe_rows[[length(pe_rows) + 1L]] <- data.table(
            GID_0 = cty, horizon = h, spec = nm, model_name = nm,
            aic_weight = cv$aic_weight[1],
            Tmean_base = Tm_base_e, Tmean_proj = Tm_proj,
            daily_change  = d_proj - d_base_e,
            annual_change = Dz_ann_proj - Dz_ann_e,
            daily_extra_PE = (d_proj - d_base_e) - (Dz_ann_proj - Dz_ann_e),
            base_source_cmip    = PE_BASE_CMIP,
            Tmean_base_cmip     = Tm_base_c,
            daily_change_cmip   = d_proj - d_base_c,
            annual_change_cmip  = Dz_ann_proj - Dz_ann_c,
            daily_extra_PE_cmip = (d_proj - d_base_c) - (Dz_ann_proj - Dz_ann_c))
        }
      }
    }

    if (!length(pe_rows)) {
      cat("  [OMIS] aucune ligne -> pe_signtest_by_country.csv NON produit.\n")
    } else {
      pe <- rbindlist(pe_rows)
      fwrite(pe, file.path(DIR_OUT, "pe_signtest_by_country.csv"))
      cat(sprintf("  [write] pe_signtest_by_country.csv (%d lignes | %d pays | %d specs)\n",
                  nrow(pe), uniqueN(pe$GID_0), uniqueN(pe$spec)))
      ## Le bras ANNUEL ne depend pas de la specification journaliere. S'il
      ## variait, ce ne serait pas la quantite que le texte decrit -- et
      ## l'emetteur s'abstiendrait sur ce motif. Le dire ici nomme le
      ## producteur plutot que le consommateur.
      nv <- pe[horizon == max(HORIZONS), uniqueN(round(annual_change_cmip, 9)),
               by = GID_0][, max(V1)]
      cat(sprintf("  controle : bras annuel invariant across specs -> %s\n",
                  if (nv == 1L) "OUI" else sprintf("NON (%d valeurs)", nv)))
      s <- pe[horizon == max(HORIZONS),
              .(extra_min = min(daily_extra_PE), extra_max = max(daily_extra_PE)),
              by = GID_0]
      cat(sprintf("  daily extra PE a %d : %s\n", max(HORIZONS),
                  paste(sprintf("%s [%.3f, %.3f]", s$GID_0, s$extra_min, s$extra_max),
                        collapse = " | ")))
      rm(pe, pe_rows, curves, era_base)
    }
  }

  ## ---- Correction de Jensen, diagnostic -----------------------------------------
  ## A function evaluated at a point is not its mean evaluated over a
  ## distribution. The gap is measured here rather than asserted.
  jensen <- fut_clim[, .(D_dist = sum(zeta_z_daily * share, na.rm = TRUE),
                         T_bar  = sum(value * share, na.rm = TRUE)),
                     by = .(GID_0, year)]
  jensen[, D_point := interp_damage(T_bar, dam_prod, "zeta_z_daily")]
  jensen[, gap := D_dist - D_point]
  fwrite(jensen, file.path(DIR_OUT, "jensen_gap_country.csv"))
  cat(sprintf("  [write] jensen_gap_country.csv | ecart median %.4f | max %.4f\n",
              median(jensen$gap, na.rm = TRUE), max(abs(jensen$gap), na.rm = TRUE)))
}

##################################################################################
# 9. DBAR -- dispersion across countries at the three horizons
#
# THE SAME WEIGHTING AS SECTION 5: omega, defined in section 4 on the common
# sample. One of the merged scripts declared itself consistent with eta through a
# COMMENT and redefined omega on its own side. A declared consistency is not one.
# coherence tenue.
##################################################################################

cat("=== 9. Dbar ===\n")

dp <- file.path(DIR_OUT, "damage_projection_country.csv")
if (!file.exists(dp)) {
  cat("  [OMIS] damage_projection_country.csv absent -> section 9 NON calculee ;\n         le controle positif de la section 10 sera donc INERTE.\n")
} else {
  dam <- fread(dp)
  n_c <- uniqueN(dam$GID_0)
  ## The restriction of section 4 makes this equality true by construction.
  ## Le controle reste : si elle cessait de tenir, la note du tableau mentirait.
  if (n_c != N_COUNTRIES)
    stop(sprintf("[STOP] %d pays dans la projection contre %d pour eta. ",
                 n_c, N_COUNTRIES),
         "eta et Dbar ne porteraient plus sur le meme echantillon, et le ",
         "controle positif de la section 10 comparerait deux objets distincts.")
  dam[, omega := OMEGA]

  res <- dam[year %in% HORIZONS & !is.na(D_z), {
    Dbar <- sum(omega * D_z)
    dc   <- D_z - Dbar
    list(Dbar_z = Dbar,
         sanity_mean    = sum(omega * dc),        # doit valoir ~0
         Dbar_cross_mad = sum(omega * abs(dc)),
         Dbar_cross_sd  = sqrt(sum(omega * dc^2)),
         n_countries    = .N)
  }, by = year]
  setorder(res, year)
  print(res)
  if (max(abs(res$sanity_mean)) > 1e-10)
    cat("  [!] la moyenne des ecarts n'est pas nulle : ponderation incoherente.\n")

  ## ---------------------------------------------------------------------
  ## tab_dbar_cross.tex N'EST PLUS EMISE  (11/09/2026)
  ## ---------------------------------------------------------------------
  ## La section F qui la portait a ete retiree du manuscrit : la table etait
  ## orpheline -- zero \\input, zero \\ref dans JMP_V11.tex, verifie avant
  ## retrait. Un exhibit genere que rien ne consomme finit par etre cite de
  ## memoire, ou re-decouvert et re-integre sans son contexte.
  ##
  ## LE CALCUL `res` RESTE, ET C'EST DELIBERE : le controle positif de la
  ## section 10 compare son bras equipondere a `res` et s'ARRETE sur
  ## divergence. Retirer le calcul avec la table supprimerait ce controle
  ## sans que rien ne le signale -- exactement le geste qu'on se refuse.
  ##
  ## La note retiree portait aussi le seul enonce correct de la convention de
  ## ponderation a deux etages (cellules -> pays : population ; pays -> monde :
  ## EQUIPONDERE, omega_r = 1/n). Si un exhibit reprend Dbar, il doit reprendre
  ## cet enonce : « pop-weighted » est faux au second etage.
  cat(sprintf(paste("  [non emis] tab_dbar_cross.tex : section F retiree du",
                    "manuscrit (%d pays, omega_r = 1/%d). Section 10 retiree le",
                    "28/09/2026.\n"), n_c, n_c))
}

##################################################################################
# 10. DBAR -- SENSIBILITE A LA PONDERATION : RETIREE LE 28/09/2026 (lot 54).
#
# Elle ne produisait que output/dbar_cross_weighting.csv, qu aucun consommateur
# ne lit depuis le 11/09/2026 : ce jour-la, la section F du manuscrit et
# tab_dbar_cross ont ete retirees, et les quatre macros EMPdbarMADequal /
# RMSequal / MADpop / RMSpop avec elles (emit_values_empirics.py, bloc « RETIREE
# (11/09/2026) »). Un calcul dont personne ne lit la sortie ne s execute pas.
# Son retrait supprime aussi un defaut d ORDRE : elle lisait population0.csv et
# identifier_xy_final.csv, que 02a n ecrit que PLUS TARD dans la chaine
# (tools/AUDIT_GRAPHE.md §5c, RECAP 3.240) -- en place, elle lisait le fichier du
# passage precedent. Code d origine : 01e_diagnostics.R.bak_20260928_lot54,
# l.673-809. S il revient, il lira une source produite AVANT 01e.
##################################################################################

cat(sprintf("\n=== 01e_diagnostics.R termine -- %.1f min ===\n",
            (proc.time() - t0)[["elapsed"]] / 60))