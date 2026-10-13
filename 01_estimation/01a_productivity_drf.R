##################################################################################
# 01a_productivity_drf.R -- damage function for PRODUCTIVITY, through trade.
#
# ESTIMATED EQUATION (Caliendo-Parro framework)
#   log(X^k_rs / X^k_rr) = [g(T_s) - g(T_r)] + gamma * theta_k * log(y_s/y_r)
#                          + gravity + product x year FE + eps
#   PPML, standard errors clustered on origin and destination.
#   The exporter is orig (s) and the importer is dest (r), so the difference is
#   orig - dest, the opposite orientation from the migration channel.
#
#   🚨 CORRIGE LE 2026-09-19 -- CE COMMENTAIRE DISAIT L'INVERSE DU CODE.
#   Il ecrivait [g(T_r) - g(T_s)] et log(y_r/y_s), c'est-a-dire importateur
#   moins exportateur, alors que prepare_diffs() (lib_drf.R l.133) calcule
#   `a - b` = orig - dest = exportateur - importateur, et que GDP_diff (l.391
#   ci-dessous) vaut log(gdp_orig/gdp_dest). La phrase qui suivait disait deja
#   « orig - dest » : le commentaire se contredisait sur trois lignes.
#   Aucun effet numerique -- c'est un commentaire -- mais il a servi de source
#   a la lecture de l'orientation, et il l'envoyait du mauvais cote.
#
#   ⚠ `gamma` est ECRIT ICI PARCE QU'IL EST ESTIME. theta_k est cuit dans le
#   regresseur GDP_diff (l.391), mais celui-ci entre par TRADE_STAGE3_CONTROLS
#   avec un coefficient LIBRE : ce n'est ni un offset ni une restriction, c'est
#   un controle. L'ecrire sans gamma laisserait croire que le coefficient est
#   impose a 1.
#
#   ⚠ LE TERME CLIMATIQUE PORTE LE THETA COMMUN, JAMAIS theta_k : les deux appels
#   a prepare_diffs() (reference : `data_ref <-` ; grille : `dp <-`) passent
#   scale = SCALE_PROD = 1/THETA_GE (00_parameters.R:249), et prepare_diffs()
#   divise le differentiel par scale (lib_drf.R:157). Le regresseur vaut donc
#   theta * Delta b, et le coefficient estime est f_z (theta = \GEtheta, lu dans
#   02c_model/specification_study.csv). Corrige le 03/10/2026 (lot 85) : ce
#   commentaire disait « scale = 1 », perime depuis le 19/09/2026.
#   Voir la convention d'unite en tete de lib_drf.R (mc_export) : ce que la
#   courbe exportee represente, et ce que le GE en fait, sont en arbitrage.
#
# IDENTIFICATION
#   Pair x product fixed effects absorb almost all the bilateral climate
#   variation, which destroys the identifying variation rather than cleaning it.
#   They are replaced by product x year fixed effects plus gravity controls, and
#   the estimate exploits variation BETWEEN pairs. The appendix measures how much
#   each set of fixed effects absorbs rather than asserting it.
#
# Reads:  data_inputs/BACI_CEPII_trade, GLORIA, EMDAT, gravity_CEPII,
#         WB_gdppercapita_country, caliendo_parro_2015_thetak,
#         intermediate_outputs/ERA_tas_distribution_grid_annual.csv
# Writes:
#   .tex   tab_stage1_trade.tex, appendix_identification_trade.tex,
#          tab_fe_absorption.tex, tab_ctrl_absorption.tex
#          (tab_robustness_prod_joint.tex est passe a ct4_emit_tab_contribA_grid.R,
#           qui l'emet sur les 15 specifications de GRID_SPECS)
#   .pdf   stability_controls_prod.pdf
#          (MC_prod_{daily,annual}_fan_paper.pdf est passe a ct6_fig_drf_paper.R)
#   .csv   MC_damage_prod_{daily,annual}_quantiles.csv  -> GE inputs
#          damage_prod_bin.csv                          -> GE input
#          fe_absorption.csv, ctrl_absorption.csv       -> absorption macros
#          prod_joint_tests_{daily,annual}.csv          -> ANCRE de non-regression
#                                                          des 12 specifications communes,
#                                                          relue par ct4
#          stage1_trade_diagnostics.csv                 -> stage 1 macros
#          aic_weights_Productivity_{Daily,Annual}.csv
#   .rds   BACI_raw_pairs_annual.rds, climate_long_annual.rds
#          (both read back by several later scripts)
# Needs:  00_parameters.R and lib_drf.R, sourced below.

##################################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(fixest);     library(splines); library(MASS); library(lfe)
  library(readxl);     library(lubridate); library(ggplot2)
})
## The root is derived from the location of this script, not from the working
## directory and not from a required environment variable: neither of those is
## reliable across the ways R can be launched. JMP_ROOT wins when it is set.
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
# 1. BILATERAL TRADE -- BACI, HS92 to ISIC Rev.3 concordance
##################################################################################

BACI <- rbindlist(lapply(
  list.files(file.path(DIR_IN, "BACI_CEPII_trade/data"), full.names = TRUE),
  function(f) { dt <- fread(f, select = 1:5)
                setnames(dt, c("year", "exporter", "importer", "HS92", "value")); dt }))
BACI <- BACI[year >= TRADE_YEAR_MIN & year <= TRADE_YEAR_MAX]
BACI <- BACI[, .(value = sum(value, na.rm = TRUE)),
             by = .(year, exporter, importer, HS92)]
BACI[, HS92 := fifelse(nchar(as.character(HS92)) == 5,
                       paste0("0", as.character(HS92)), as.character(HS92))]

## The .CSV extension is UPPER CASE because that is the name the provider
## delivers, and the package convention is that every input keeps its supplier
## name. Written in lower case this read works on Windows and fails on any Linux
## clone, where case matters.
conc <- fread(file.path(DIR_IN, "BACI_CEPII_trade/WITS_Concordance_H0_to_I3.CSV"),
              colClasses = "character")
setnames(conc, c("HS92", "description_hs92",
                 "product_code_ISICRev3_4digits", "description_rev3"))
conc[, HS92 := fifelse(nchar(HS92) == 5, paste0("0", HS92), HS92)]
BACI <- merge(BACI, conc, by = "HS92", all.x = TRUE)

## Two-digit ISIC division. Integer division handles "0111" as "111" correctly:
## substr() would return "11" instead of "01" for divisions below 10.
BACI[, product_code_ISICRev3 := sprintf("%02d",
  suppressWarnings(as.integer(product_code_ISICRev3_4digits)) %/% 100L)]
BACI[!grepl("^[0-9]{2}$", product_code_ISICRev3), product_code_ISICRev3 := NA_character_]
BACI <- BACI[!is.na(product_code_ISICRev3) & !is.na(value)]
BACI <- BACI[, .(value = sum(value, na.rm = TRUE)),
             by = .(year, exporter, importer, product_code_ISICRev3)]

## ---- sectoral trade elasticities (Caliendo and Parro 2015) -------------------
theta_k <- fread(file.path(DIR_IN, "BACI_CEPII_trade/caliendo_parro_2015_thetak/theta_k.csv"),
                 sep = ";", colClasses = "character", select = 1:4, fill = TRUE)
BACI <- merge(BACI, theta_k, by = "product_code_ISICRev3", all.x = TRUE)
cat(sprintf("[theta_k] %d secteurs appariés | manufacture 15-37 présente : %s\n",
            uniqueN(BACI[!is.na(theta_k)]$product_code_ISICRev3),
            any(BACI$product_code_ISICRev3 %in% sprintf("%02d", 15:37) & !is.na(BACI$theta_k))))
BACI <- BACI[!is.na(theta_k)][, theta_k := as.numeric(theta_k)]

country <- fread(file.path(DIR_IN, "BACI_CEPII_trade/country_codes_V202401.csv"),
                 select = c("country_code", "country_iso3"))
setnames(country, c("code", "iso3"))
BACI <- merge(BACI, country[, .(exporter = code, iso3_o = iso3)], by = "exporter")
BACI <- merge(BACI, country[, .(importer = code, iso3_d = iso3)], by = "importer")

## Entities that are not countries in the sense of the spatial model.
not_countries <- c("ABW","CUW","HKG","FSM","GRL","XKX","XXA","NCL","SMR",
                   "TWN","MAC","PSE","SXM","PYF","PLW","TCA","CYM","SSD")
BACI <- BACI[!(iso3_o %in% not_countries) & !(iso3_d %in% not_countries)]
BACI <- BACI[!is.na(sector)]
BACI <- BACI[, .(year, iso3_o, iso3_d, product_code_ISICRev3, value, theta_k, sector)]

##################################################################################
# 2. GDP PER CAPITA -- annual interpolation
##################################################################################

gdp <- load_wb(file.path(DIR_IN, "WB_gdppercapita_country/total/total_GDP_country_constant2015USD.csv"),
               "gdppercapita", period = "annual",
               year_min = TRADE_YEAR_MIN, year_max = TRADE_YEAR_MAX)

##################################################################################
# 3. GLORIA -- denominator X^k_rr, domestic absorption
#
# Three matching traps, all silent: they produce NA rather than an error.
#   (i)   the labels come as "(COUNTRY) Product name product", so BOTH the
#         prefix and the suffix must be stripped; the concordance table holds
#         the bare name;
#   (ii)  its codes are integers where the trade data carries two digits;
#   (iii) some labels contain double spaces, so both sides are normalised.
##################################################################################

## Label normalisation, applied on BOTH sides of the join.
norm_lab <- function(x) trimws(gsub("\\s+", " ", x))

gloria_map <- fread(file.path(DIR_IN, "BACI_CEPII_trade/gloria_to_isic3.csv"),
                    colClasses = "character")
setnames(gloria_map, c("product", "Code"),
         c("sector_gloria", "product_code_ISICRev3"), skip_absent = TRUE)
gloria_map[, sector_gloria := norm_lab(sector_gloria)]
gloria_map[, product_code_ISICRev3 := sprintf("%02d",
             as.integer(product_code_ISICRev3))]
gloria_map <- unique(gloria_map[, .(sector_gloria, product_code_ISICRev3)])

readme <- read_excel(file.path(DIR_IN, "GLORIA/GLORIA_ReadMe_057.xlsx"),
                     sheet = "Sequential region-sector labels")
## L'annee se lit dans le NOM, et la garde exige exactement les millesimes
## TRADE_YEAR_MIN..TRADE_YEAR_MAX et aucun autre fichier (lib_drf.R,
## gloria_year_files ; lot 34). Avant : annee = position dans le listage.
.gloria <- gloria_year_files(file.path(DIR_IN, "GLORIA/data_gloria/97"),
                             years = TRADE_YEAR_MIN:TRADE_YEAR_MAX)
gloria_files <- .gloria$files
gloria_years <- .gloria$years

rs        <- readme$Sequential_regionSector_labels
fd_lab    <- readme$Sequential_finalDemand_labels[1:984]
ctry_col  <- gsub(".*\\(([A-Z]{3})\\).*", "\\1", fd_lab)
ctry_row  <- gsub(".*\\(([A-Z]{3})\\).*", "\\1", rs)
prod_rows <- which(grepl("product$", rs, ignore.case = TRUE))
## Country prefix stripped, then the trailing "product".
##
## THE STRIP IS ANCHORED ON THE COUNTRY CODE, AND THAT IS NOT COSMETIC. The
## prefix is a country name followed by its three-letter code: "Rest of Americas
## (XAM) ". Requiring \([A-Z]{3}\) anchors the strip on that code.
##
## A greedy ".*\\) " -- what stood here until 09/09/2026 -- runs to the LAST ") "
## in the string, so it swallows the product name whenever the PRODUCT itself
## carries parentheses. Measured: "Growing beverage crops (coffee, tea etc)
## product", the ONLY label of the map in that case, was reduced to the literal
## "product" and matched nothing. Its domestic absorption never entered the
## denominator of ISIC division 01, which came out 0.8 % short.
##
## THE GUARDRAIL HAD BEEN FIRING ALL ALONG and nobody read it: every run printed
## "119/120 libellés appariés" and "non appariés : product". Hence the hard stop
## below -- a label that stops matching must halt the step, not print a line.
sector_lab <- norm_lab(sub("\\s+product$", "",
                           sub(".*\\([A-Z]{3}\\)\\s*", "", rs[prod_rows])))

## Matching check BEFORE reading 25 files: if the labels do not line up, it is
## better to know immediately than after the reads.
{
  ok <- unique(sector_lab) %in% gloria_map$sector_gloria
  cat(sprintf("[GLORIA] %d/%d libellés appariés à la table ISIC\n",
              sum(ok), length(ok)))
  if (!any(ok))
    stop("⛔ AUCUN libellé GLORIA apparié : les étiquettes du ReadMe ne ",
         "correspondent pas à gloria_to_isic3.csv.")
  ## ARRET DUR SUR TOUT LIBELLE NON APPARIE, pose le 09/09/2026.
  ## Les deux pieges d'etiquette de ce jeu sont identifies et traites, et le
  ## couple est COMPLEMENTAIRE -- aucun des deux producteurs n'avait les deux :
  ##   - un DOUBLE ESPACE interne dans "Alcoholic and other  beverages" (div. 15),
  ##     seul libellé du ReadMe dans ce cas : norm_lab() le collapse des DEUX cotes
  ##     de la jointure. Le producteur hors paquet ne le collapsait pas et perdait
  ##     22 % du denominateur de la division 15 ;
  ##   - une PARENTHESE dans le nom du produit (div. 01) : le regex ancre sur
  ##     \([A-Z]{3}\) ci-dessus. Le paquet, avec son ".*\\) " glouton, perdait
  ##     0,8 % du denominateur de la division 01.
  ## 120/120 est donc atteignable et ATTEINT. Un libellé qui cesse de s'apparier
  ## est une DERIVE MECANIQUE d'etiquette, pas un fait scientifique : elle doit
  ## arreter l'etape. Un denominateur de production silencieusement incomplet est
  ## exactement ce qui a coute 0,8 % a la division 01 pendant toute la serie de
  ## runs qui precede.
  if (any(!ok))
    stop(sprintf("⛔ %d/%d libellés GLORIA NON appariés : ", sum(!ok), length(ok)),
         paste(unique(sector_lab)[!ok], collapse = " | "),
         "\n  Un libellé non apparié retire sa production du denominateur X^k(rr) ",
         "de sa division, SANS RIEN DIRE. Verifier le strip du prefixe (ancre sur ",
         "le code pays) et la normalisation des espaces, des DEUX cotes de la ",
         "jointure, avant de relacher cette garde.")
}

final_gloria <- rbindlist(lapply(seq_along(gloria_files), function(i) {
  fd  <- fread(gloria_files[i])
  mat <- as.matrix(fd[prod_rows, lapply(.SD, as.numeric), .SDcols = 1:984])
  cp  <- ctry_row[prod_rows]
  dt <- rbindlist(lapply(unique(cp), function(ct) {
    ri <- which(cp == ct); ci <- which(ctry_col == ct)
    if (!length(ci)) return(NULL)
    data.table(GID_0 = ct, sector_gloria = sector_lab[ri],
               value_gloria = rowSums(mat[ri, ci, drop = FALSE], na.rm = TRUE))
  }))
  dt[, year := gloria_years[i]][]
}))

final_gloria <- merge(final_gloria, gloria_map, by = "sector_gloria", all.x = TRUE)
n_na <- sum(is.na(final_gloria$product_code_ISICRev3))
if (n_na > 0)
  cat(sprintf("⚠ %d/%d lignes GLORIA sans code ISIC — écartées.\n",
              n_na, nrow(final_gloria)))
final_gloria <- final_gloria[!is.na(product_code_ISICRev3),
  .(value_gloria = sum(value_gloria, na.rm = TRUE)),
  by = .(GID_0, year, product_code_ISICRev3)]
cat(sprintf("[GLORIA] %d lignes | %d pays | %d secteurs | %d années\n",
            nrow(final_gloria), uniqueN(final_gloria$GID_0),
            uniqueN(final_gloria$product_code_ISICRev3),
            uniqueN(final_gloria$year)))

##################################################################################
# 4. GRAVITY, DISASTERS, AND THE DEPENDENT VARIABLE
##################################################################################

gravity <- fread(file.path(DIR_IN, "gravity_CEPII/Gravity_V202211.csv"),
                 select = c("year", "iso3_o", "iso3_d", "distw_harmonic", "distcap",
                            "contig", "comlang_off", "sibling_ever",
                            "col_dep_ever", "fta_wto"))
gravity <- gravity[year >= TRADE_YEAR_MIN & year <= TRADE_YEAR_MAX & !is.na(comlang_off)]
setnames(gravity, c("iso3_o", "iso3_d"), c("orig", "dest"))

emdat_file <- list.files(file.path(DIR_IN, "EMDAT"), pattern = "\\.xlsx$",
                         full.names = TRUE)[1]
EMDAT <- as.data.table(read_excel(emdat_file, sheet = 1)) %>%
  dplyr::select(`Start Year`, ISO, `Start Month`, `Start Day`,
                `End Month`, `End Day`, `Total Affected`) %>%
  setNames(c("year", "orig", "sm", "sd", "em", "ed", "total")) %>%
  filter(year >= TRADE_YEAR_MIN, year <= TRADE_YEAR_MAX) %>%
  mutate(duration_days = pmax(0, as.numeric(
    make_date(year, em, ed) - make_date(year, sm, sd) + 1), na.rm = TRUE)) %>%
  group_by(orig, year) %>%
  summarise(nb_events = n(),
            total_duration_days = sum(duration_days, na.rm = TRUE),
            total_affected      = sum(total, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()

## ---- Assemblage ---------------------------------------------------------------
BACI <- merge(BACI, gdp[, .(iso3_o = GID_0, year, gdp_orig = gdppercapita)],
              by = c("iso3_o", "year"), all.x = TRUE)
BACI <- merge(BACI, gdp[, .(iso3_d = GID_0, year, gdp_dest = gdppercapita)],
              by = c("iso3_d", "year"), all.x = TRUE)
BACI <- merge(BACI, final_gloria[, .(iso3_d = GID_0, year, product_code_ISICRev3,
                                     value_gloria)],
              by = c("iso3_d", "year", "product_code_ISICRev3"), all.x = TRUE)

##################################################################################
# 3bis. LA PERTE SECTORIELLE, COMPTEE ET DITE
#
# POURQUOI CE BLOC EXISTE (09/09/2026).
# La nomenclature GLORIA ne porte PAS toutes les divisions ISIC Rev.3 que le
# commerce BACI porte. Les divisions absentes n'obtiennent aucun denominateur de
# production, donc `value_gloria` est NA, et le filtre `!is.na(total_gloria)`
# plus bas les ecarte. C'EST LE COMPORTEMENT VOULU -- mais il etait SILENCIEUX :
# des divisions entieres disparaissaient sans qu'aucune ligne ne le dise.
#
# LE MOTIF DE L'EXCLUSION. La variable expliquee est un RATIO, X^k_rs / X^k_ss :
# le commerce bilateral sur l'absorption domestique de LA MEME division. Replier
# une division absente sur une division voisine -- ce que faisait le producteur
# hors paquet, 18 -> 19 et 30, 33 -> 32 -- met au denominateur la production
# d'une AUTRE division. Le ratio cesse alors de mesurer ce qu'il nomme. On
# prefere perdre la division que publier un ratio dont le denominateur n'est pas
# le sien.
#
# CE BLOC NE FILTRE RIEN. Il mesure et il imprime ; le filtre reste a sa place.
##################################################################################
{
  codes_baci   <- sort(unique(BACI$product_code_ISICRev3))
  codes_gloria <- sort(unique(final_gloria$product_code_ISICRev3))
  sans_denom   <- setdiff(codes_baci, codes_gloria)

  perdu <- BACI[product_code_ISICRev3 %in% sans_denom]
  garde <- BACI[!product_code_ISICRev3 %in% sans_denom]

  ## Paires qui disparaissent A CAUSE de cette perte : celles qui tiennent les
  ## 25 annees completes AVANT le retrait et ne les tiennent plus apres. C'est la
  ## seule definition qui impute la perte a l'exclusion sectorielle plutot qu'aux
  ## filtres voisins -- un comptage brut de paires melangerait les deux.
  pa <- BACI[,  .(n = uniqueN(year)), by = .(iso3_o, iso3_d)][n == TRADE_N_YEARS]
  pb <- garde[, .(n = uniqueN(year)), by = .(iso3_o, iso3_d)][n == TRADE_N_YEARS]
  paires_perdues <- nrow(pa) - nrow(pb)

  cat(sprintf("[gloria] %d divisions ISIC du panel BACI sur %d, soit %s, n'ont AUCUN\n",
              length(sans_denom), length(codes_baci),
              paste(sans_denom, collapse = ", ")))
  cat("         denominateur GLORIA. Elles sont ECARTEES, et non repliees sur une\n")
  cat("         division voisine : le ratio X_rs/X_ss exige le MEME secteur aux\n")
  cat("         deux etages.\n")
  if (length(sans_denom)) {
    cat(sprintf("         perte : %s lignes sur %s (%.2f %%) | %.2f %% du commerce en valeur\n",
                format(nrow(perdu), big.mark = ","), format(nrow(BACI), big.mark = ","),
                100 * nrow(perdu) / nrow(BACI),
                100 * sum(perdu$value, na.rm = TRUE) / sum(BACI$value, na.rm = TRUE)))
    cat(sprintf("         paires perdant de ce fait leurs %d annees completes : %d\n",
                TRADE_N_YEARS, paires_perdues))
    print(perdu[, .(lignes = .N,
                    pct_valeur = round(100 * sum(value, na.rm = TRUE) /
                                       sum(BACI$value, na.rm = TRUE), 3)),
                by = product_code_ISICRev3][order(-pct_valeur)])
  }

  ## ---- CSV COMPAGNON : la perte devient trois macros -------------------------
  ## La note de bas de page de la Section 4 nomme cette exclusion. Un ordre de
  ## grandeur -- 9 % du commerce en valeur -- n'est pas une note anodine, et le
  ## taper en dur serait le defaut que le paquet retire partout. Le CSV est
  ## ecrit MEME quand rien n'est perdu : une note qui ne peut plus se remplir
  ## doit se voir a la compilation, pas s'effacer.
  fwrite(data.table(
    n_divisions_dropped = length(sans_denom),
    divisions           = paste(sans_denom, collapse = ", "),
    n_divisions_baci    = length(codes_baci),
    n_rows_lost         = nrow(perdu),
    n_rows_total        = nrow(BACI),
    pct_rows            = 100 * nrow(perdu) / nrow(BACI),
    pct_value           = 100 * sum(perdu$value, na.rm = TRUE) /
                                sum(BACI$value, na.rm = TRUE),
    pairs_lost          = paires_perdues,
    n_years_required    = TRADE_N_YEARS),
    file.path(DIR_OUT, "trade_gloria_loss.csv"))
  cat("         [write] trade_gloria_loss.csv (3 macros \\EMPtradeGloriaLoss*)\n")

  rm(perdu, garde, pa, pb); invisible(gc(verbose = FALSE))
}

## DISTRIBUTION DE k -- combien de lignes chaque groupe collapse.
## `value_gloria` est CONSTANT dans le groupe : la jointure ci-dessus se fait sur
## une cle unique de final_gloria. D'ou `first()` et non `sum()` a la ligne
## suivante. Le producteur hors paquet ecrivait `sum(value_gloria)`, ce qui
## rendait k fois la valeur -- une multiplication par un compte de lignes, pas
## une agregation. k est imprime ici pour que le facteur de gonflement soit un
## CHIFFRE du paquet et non une deduction de lecture.
{
  kk <- BACI[, .N, by = .(year, iso3_o, iso3_d, product_code_ISICRev3)]$N
  cat(sprintf("[k] lignes par groupe : k=1 %.2f %% | k>=2 %.2f %% | median %g | max %d\n",
              100 * mean(kk == 1L), 100 * mean(kk > 1L), median(kk), max(kk)))
  cat(sprintf("    denominateur gonfle d'un facteur moyen %.4f sous sum() au lieu de first()\n",
              mean(kk)))
  print(table(pmin(kk, 6L), dnn = "k (6 = 6 et plus)"))
  rm(kk); invisible(gc(verbose = FALSE))
}

BACI <- BACI[, .(theta_k_weighted = mean(theta_k),
                 total = sum(value, na.rm = TRUE),
                 total_gloria = first(value_gloria),
                 gdp_orig = mean(gdp_orig), gdp_dest = mean(gdp_dest)),
             by = .(year, iso3_o, iso3_d, product_code_ISICRev3)]
setnames(BACI, c("iso3_o", "iso3_d"), c("orig", "dest"))

BACI <- merge(BACI, gravity, by = c("orig", "dest", "year"), all.x = TRUE)
BACI <- merge(BACI, EMDAT,   by = c("orig", "year"),         all.x = TRUE)
BACI[is.na(nb_events), (DISASTER_CONTROLS) := 0]

BACI <- BACI[!is.na(total_gloria) & total_gloria > 0]
BACI[, explained := total / total_gloria]
BACI <- BACI[is.finite(explained)]

## theta_k * log(y_orig / y_dest): exporter over importer.
BACI[, GDP_diff := theta_k_weighted * log(gdp_orig / gdp_dest)]
BACI <- BACI[is.finite(GDP_diff)]

## Balanced panel: climate identifies BETWEEN pairs, and an unbalanced panel
## would inject composition into that comparison.
BACI[, nb_years := uniqueN(year), by = .(orig, dest)]
BACI <- BACI[nb_years == TRADE_N_YEARS][, nb_years := NULL]

BACI <- BACI[!product_code_ISICRev3 %in% c(ISIC_EXTRACTIVE, ISIC_SERVICES)]

## Ratio cap: removes entrepot flows and accounting mismatches.
thresh <- min(quantile(BACI$explained, DV_QUANTILE_CAP, na.rm = TRUE), DV_HARD_CAP)
BACI <- BACI[explained <= thresh]
cat(sprintf("[BACI] N = %s | %d secteurs | plafond DV = %.2f\n",
            format(nrow(BACI), big.mark = ","),
            uniqueN(BACI$product_code_ISICRev3), thresh))

saveRDS(BACI, file.path(DIR_TMP_OUT, "BACI_raw_pairs_annual.rds"))

##################################################################################
# 5. CLIMATE -- daily distribution aggregated to the country, no recoding
##################################################################################

climate <- fread(file.path(DIR_TMP, "ERA_tas_annual.csv"))
setnames(climate, c("GID_0", "year", "value", "count", "share"))
climate <- climate[year >= TRADE_YEAR_MIN & year <= TRADE_YEAR_MAX]
climate <- climate[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
saveRDS(climate, file.path(DIR_TMP_OUT, "climate_long_annual.rds"))

##################################################################################
# 6. STAGE 1 -- GRAVITY WITHOUT TEMPERATURE
# Fixed effects: origin x year x product and destination x year x product.
##################################################################################

cat("=== Stage 1 ===\n")
{
  dt <- as.data.frame(BACI)
  dt$fe_o <- paste(dt$orig, dt$year, dt$product_code_ISICRev3, sep = "_")
  dt$fe_d <- paste(dt$dest, dt$year, dt$product_code_ISICRev3, sep = "_")
  s1 <- feglm(as.formula(paste0("explained ~ ",
                paste(TRADE_STAGE3_CONTROLS, collapse = " + "))),
              data = dt, family = poisson(), cluster = c("orig", "dest"),
              fixef = c("fe_o", "fe_d"))
}

emit_stage1_table <- function(mod, data, output_path) {
  sq <- round(fitstat(mod, "sq.cor")[[1]], 3)
  ll <- round(as.numeric(logLik(mod)), 0); n <- nobs(mod)
  w  <- wald(mod, GRAVITY_CONTROLS[GRAVITY_CONTROLS %in% names(coef(mod))])
  disp <- round(sum((mod$y - fitted(mod))^2 / pmax(fitted(mod), 1e-10)) /
                  (length(mod$y) - length(coef(mod))), 3)
  zeros <- round(mean(data$explained == 0, na.rm = TRUE), 3)

  ## CSV keyed by a STABLE IDENTIFIER, never by a display label. Parsing the
  ## rendered table by regex breaks as soon as the table changes shape.
  ## n_sectors ajoute le 2026-09-07. La macro \\EMPtradeNsectors n'avait plus
  ## aucune source : l'emetteur la cherchait par regex sur une ligne
  ## "N. sectors (ISIC Rev.3)" que cette table n'ecrit plus, et la cle manquait
  ## ici. Elle decrit l'echantillon au manuscrit (l.478) ; mesure du 07/09 : 20,
  ## soit le litteral publie. L'emetteur prefere le CSV : poser la cle suffit.
  n_sec <- if ("product_code_ISICRev3" %in% names(data))
             length(unique(data$product_code_ISICRev3)) else NA_integer_
  if (is.na(n_sec))
    warning("stage1_trade_diagnostics.csv : product_code_ISICRev3 absent, ",
            "n_sectors non calcule -- \\EMPtradeNsectors restera non emise.")
  fwrite(data.frame(
    key = c("sq_cor", "loglik", "nobs", "wald_f", "wald_p", "dispersion",
            "share_zeros", "n_sectors"),
    value = c(sq, ll, n, w$stat, w$p, disp, zeros, n_sec)),
    file.path(DIR_OUT, "stage1_trade_diagnostics.csv"))

  pa <- data.frame(
    s = c("Pseudo-$R^2$ (sq. corr.)", "Log-likelihood", "Observations",
          "Wald $F$ (gravity controls)", "$p$-value (Wald)",
          "PPML dispersion", "Share zero flows"),
    v = c(sprintf("%.3f", sq), formatC(ll, format = "d", big.mark = ","),
          formatC(n, format = "d", big.mark = ","), sprintf("%.2f", w$stat),
          sprintf("%.4f", w$p), sprintf("%.3f", disp), sprintf("%.3f", zeros)))
  nm <- c(distw_harmonic = "Log distance (harmonic)",
          distcap = "Log distance (capitals)", comlang_off = "Common language",
          col_dep_ever = "Colonial link", fta_wto = "FTA (WTO)",
          GDP_diff = "$\\theta_k\\cdot\\log(y_o/y_d)$")
  cf <- coef(mod); se <- sqrt(diag(vcov(mod)))
  keep <- names(cf)[names(cf) %in% names(nm)]
  pb <- data.frame(v = nm[keep], c = sprintf("%.4f", cf[keep]),
                   s = sprintf("(%.4f)", se[keep]))

  writeLines(c(
    "%% GENERATED by 01a_productivity_drf.R -- DO NOT EDIT BY HAND.",
    "\\begin{table}[H]", "\\centering\\small",
    "\\caption{Stage 1 gravity regression diagnostics --- Trade}",
    "\\label{tab:stage1_trade}", "\\bigskip",
    "\\textbf{Panel A: Diagnostics}\\\\[4pt]",
    "\\begin{tabular}{lc}", "\\hline\\hline", "Statistic & Value \\\\", "\\hline",
    apply(pa, 1, function(r) paste0(r[1], " & ", r[2], " \\\\")), "\\hline\\hline",
    "\\multicolumn{2}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    "PPML with origin$\\times$year$\\times$product and",
    "destination$\\times$year$\\times$product fixed effects.",
    "Dispersion $=\\sum(y_i-\\hat y_i)^2/\\hat y_i/(N-K)$.}",
    "\\end{tabular}", "\\par\\bigskip",
    "\\textbf{Panel B: Gravity control coefficients}\\\\[4pt]",
    "\\begin{tabular}{lcc}", "\\hline\\hline",
    "Variable & Coefficient & Std. Error \\\\", "\\hline",
    apply(pb, 1, function(r) paste0(r[1], " & ", r[2], " & ", r[3], " \\\\")),
    "\\hline\\hline",
    "\\multicolumn{3}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
    "Standard errors clustered by origin and destination.}",
    "\\end{tabular}", "\\end{table}"), output_path)
  cat(sprintf("[write] %s\n", basename(output_path)))
}
emit_stage1_table(s1, BACI, file.path(DIR_OUT, "tab_stage1_trade.tex"))

##################################################################################
# 7. IDENTIFICATION APPENDIX -- three measured arguments
#
# 1. pair x product fixed effects absorb almost all the climate variation;
# 2. gravity controls absorb part of the between-pair confounders;
# 3. the curve is stable with and without those controls (Oster bound).
##################################################################################

data_ref <- prepare_diffs(BACI, climate, REF_SPEC, REF_DEG, REF_WLOWER, REF_WUPPER,
                          scale = SCALE_PROD)
data_ref <- stage3_climate_filter(data_ref, REF_SPEC, REF_DEG)   # lot 41 : ecrit, plus herite
data_ref[, fe_orig_dest_prod := paste(orig, dest, product_code_ISICRev3, sep = "_")]

r2_on <- function(col, fe) {
  y <- as.numeric(data_ref[[col]])
  r <- as.numeric(residuals(felm(as.formula(paste0(col, " ~ 1 | ", fe)),
                                 data = as.data.frame(data_ref))))
  list(r2 = 1 - var(r) / var(y, na.rm = TRUE), sd_ratio = sqrt(var(r) / var(y, na.rm = TRUE)))
}

arg1 <- rbindlist(lapply(seq_len(REF_DEG), function(k) {
  col <- paste0("diff_", REF_SPEC, "_", k)
  a <- r2_on(col, "fe_orig_dest_prod"); b <- r2_on(col, "fe_prod_year")
  c_ <- r2_on(col, "fe_orig_dest_prod + fe_prod_year")
  ## sd_within_ratio is computed HERE on the unrounded residuals. Re-deriving it
  ## from an r2 rounded to four decimals propagates a few hundredths of a point
  ## through the square root and makes the first decimal unrecoverable.
  data.table(term = k, r2_pair = round(a$r2, 4), r2_pyr = round(b$r2, 6),
             r2_both = round(c_$r2, 4),
             r2_pair_raw = a$r2, r2_pyr_raw = b$r2, sd_within_ratio = a$sd_ratio)
}))
fwrite(arg1, file.path(DIR_OUT, "fe_absorption.csv"))

arg2 <- rbindlist(lapply(seq_len(REF_DEG), function(k) {
  col <- paste0("diff_", REF_SPEC, "_", k)
  r2_no <- r2_on(col, "fe_prod_year")$r2
  r2_ct <- summary(felm(as.formula(paste0(col, " ~ ",
             paste(TRADE_STAGE3_CONTROLS, collapse = " + "), " | fe_prod_year")),
             data = as.data.frame(data_ref)))$r2
  data.table(term = k, r2_fe_only = round(r2_no, 6), r2_fe_ctrl = round(r2_ct, 4),
             r2_absorbed = round(r2_ct - r2_no, 4), r2_absorbed_raw = r2_ct - r2_no)
}))
fwrite(arg2, file.path(DIR_OUT, "ctrl_absorption.csv"))
cat("[write] fe_absorption.csv + ctrl_absorption.csv\n")

## Les deux membres de la comparaison portent leurs controles EXPLICITEMENT.
## ⛔ `m_ct` s'en remettait au defaut de run_stage3(), qui vaut character(0) :
## les deux appels etaient donc identiques, la figure tracait deux fois la meme
## courbe sous deux etiquettes, et le log rendait « correlation avec/sans
## controles = 1.0000 » — un vert qui ne verifiait rien. Corrige le 2026-09-04.
## Cause : le refactor qui a unifie les quatre variantes en run_stage3() a
## remplace le defaut `controls = gravity_controls` de drf_prod_final.r l.577
## par character(0), et les sites d'appel de 01a s'en remettaient au defaut.
## Un defaut qui change de valeur ne previent pas ses appelants : les controles
## s'ecrivent desormais a chaque site, comme 01b l.269 et 01c le font deja.
m_no <- run_stage3(data_ref, REF_SPEC, REF_DEG, controls = character(0))
m_ct <- run_stage3(data_ref, REF_SPEC, REF_DEG, controls = TRADE_STAGE3_CONTROLS)
stage3_nobs_check(m_no, data_ref, "reference sans controles")
stage3_nobs_check(m_ct, data_ref, "reference avec controles")
for (m in c("m_no", "m_ct")) {
  assign(m, `[[<-`(get(m), "wval_lower", REF_WLOWER))
  assign(m, `[[<-`(get(m), "wval_upper", REF_WUPPER))
}
c_no <- get_effect_curve(m_no, climate)[, controls := "Without gravity controls"]
c_ct <- get_effect_curve(m_ct, climate)[, controls := "With gravity controls"]
corr_stability <- cor(c_no$y, c_ct$y, use = "complete.obs")
cat(sprintf("[stab] corrélation avec/sans contrôles = %.4f\n", corr_stability))

p_stab <- ggplot(rbind(c_no, c_ct), aes(x, y, colour = controls, linetype = controls)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = controls), alpha = .12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c("Without gravity controls" = "#d73027",
                                 "With gravity controls" = "#2166ac")) +
  scale_fill_manual(values = c("Without gravity controls" = "#d73027",
                               "With gravity controls" = "#2166ac")) +
  scale_linetype_manual(values = c("Without gravity controls" = "dashed",
                                   "With gravity controls" = "solid")) +
  scale_x_continuous(breaks = seq(-5, 35, 5), labels = function(x) paste0(x, "°C")) +
  labs(x = "Daily mean temperature", y = FAN_YLAB[["prod"]], colour = NULL,
       fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 14) + theme(legend.position = "bottom")
ggsave(file.path(DIR_OUT, "stability_controls_prod.pdf"), p_stab,
       width = 9, height = 6, dpi = 300)

## ---- the two tables and the appendix paragraph -------------------------------
## Every value in these notes is INTERPOLATED from the CSV files, never typed. A
## hand-written note had quoted the most favourable term as if it were the bound,
## and had read a fraction as a percentage.
writeLines(c(
  "%% GENERATED by 01a_productivity_drf.R -- DO NOT EDIT BY HAND.",
  "\\begin{table}[H]", "\\centering\\small",
  "\\caption{Variance absorption by fixed effects}", "\\label{tab:fe_absorption}",
  "\\begin{tabular}{lccc}", "\\hline\\hline",
  "Climate term & Pair$\\times$prod FE & Prod$\\times$year FE & Both FE \\\\", "\\hline",
  sprintf("$\\Delta_%d$ & %.4f & %.6f & %.4f \\\\",
          arg1$term, arg1$r2_pair, arg1$r2_pyr, arg1$r2_both),
  "\\hline\\hline",
  "\\multicolumn{4}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
  "$R^2$ from projecting the bilateral climate differential",
  "$\\Delta_k = g_k(T(s)) - g_k(T(r))$ on fixed effects alone (OLS).",
  "Origin$\\times$destination$\\times$product fixed effects absorb",
  sprintf("at least %.1f\\%% of the variance in bilateral climate differentials",
          100 * min(arg1$r2_pair_raw)),
  sprintf("(%.1f\\%%--%.1f\\%% across the three terms), leaving insufficient",
          100 * min(arg1$r2_pair_raw), 100 * max(arg1$r2_pair_raw)),
  "residual variation for PPML identification.",
  sprintf("Product$\\times$year fixed effects absorb at most %.2f\\%%,",
          100 * max(arg1$r2_pyr_raw)),
  "preserving the between-pair variation that identifies $g(\\cdot)$.}",
  "\\end{tabular}", "\\end{table}"),
  file.path(DIR_OUT, "tab_fe_absorption.tex"))

writeLines(c(
  "%% GENERATED by 01a_productivity_drf.R -- DO NOT EDIT BY HAND.",
  "\\begin{table}[H]", "\\centering\\small",
  "\\caption{Absorption of between-pair confounders by gravity controls}",
  "\\label{tab:ctrl_absorption}",
  "\\begin{tabular}{lccc}", "\\hline\\hline",
  "Climate term & $R^2$ (FE only) & $R^2$ (FE + controls) & Absorbed \\\\", "\\hline",
  sprintf("$\\Delta_%d$ & %.6f & %.4f & %.4f \\\\",
          arg2$term, arg2$r2_fe_only, arg2$r2_fe_ctrl, arg2$r2_absorbed),
  "\\hline\\hline",
  "\\multicolumn{4}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
  "OLS projection of $\\Delta_k$ on gravity controls conditional on",
  "product$\\times$year fixed effects. Gravity controls (log bilateral distance,",
  ## 🚨 23/09/2026 (lot 25) : `$\theta^k$` RETABLI. Le retrait du 22/09
  ## reposait sur une confusion : « libre » qualifiait le COEFFICIENT de
  ## `GDP_diff`, pas sa CONSTRUCTION. Le regresseur est theta_k * log(y_o/y_d)
  ## (01a l.411), theta_k variant par division (17 valeurs sur 20 secteurs) ;
  ## un coefficient scalaire n'en absorbe que R2 = 0,38 (RECAP §3.197). La
  ## legende nomme donc le regresseur TEL QU'IL EST CONSTRUIT ; son
  ## coefficient reste libre, ce que dit l'eq. (5) (zeta_y).
  ## ⛔ Le commentaire ci-dessous est la justification RETIREE, gardee comme trace :
  ## 22/09/2026 : `$\theta^k\cdot\log(y_s/y_r)$` -> `$\log(y_s/y_r)$`.
  ## DEUX erreurs dans le meme symbole, et c'est le troisieme site de la famille
  ## (eq. 5 du manuscrit en portait deux autres) :
  ##   1. l'indice k. Le secteur n'entre dans l'etage 3 QUE par `fe_prod_year`
  ##      (lib_drf.R l.157) ; les coefficients climatiques sont un vecteur
  ##      POOLE, et aucun theta^k n'est impose terme a terme.
  ##   2. le coefficient contraint. `GDP_diff` est un membre de
  ##      TRADE_STAGE3_CONTROLS (00_parameters.R l.417), donc un regresseur
  ##      LIBRE passe a feglm : rien n'egalise son coefficient a theta. Ecrire
  ##      `theta^k * log(y_s/y_r)` dans une liste de controles affirmait une
  ##      restriction que l'estimateur n'impose pas.
  ## La legende nomme donc le REGRESSEUR, sans prejuger de son coefficient.
  "common language, colonial link, WTO FTA, $\\theta^k\\cdot\\log(y_s/y_r)$)",
  sprintf("absorb %.1f--%.1f\\%% of the between-pair climate variation,",
          100 * min(arg2$r2_absorbed_raw), 100 * max(arg2$r2_absorbed_raw)),
  "partially controlling for the main observable confounders.}",
  "\\end{tabular}", "\\end{table}"),
  file.path(DIR_OUT, "tab_ctrl_absorption.tex"))

## No \label and no \subsection* here: the manuscript defines both just before
## its \input. Carrying them in the included file too would define them twice.
## The anchor and the title belong to the including document.
writeLines(c(
  "%% GENERATED by 01a_productivity_drf.R -- DO NOT EDIT BY HAND.", "",
  "The productivity specification replaces",
  "origin$\\times$destination$\\times$product fixed effects with",
  "product$\\times$year fixed effects and gravity controls, exploiting",
  "between-pair variation in bilateral temperature differentials",
  "$[g(T_t(s)) - g(T_t(r))]$ to identify the dose-response $g(\\cdot)$.",
  "We provide three pieces of evidence supporting this choice.", "",
  "\\paragraph{Fixed effects absorb identifying variation.}",
  "Table~\\ref{tab:fe_absorption} shows that",
  sprintf("origin$\\times$destination$\\times$product fixed effects absorb at least %.1f\\%%",
          100 * min(arg1$r2_pair_raw)),
  sprintf("of the variance in bilateral climate differentials $\\Delta_k$ (%.1f\\%%--%.1f\\%%",
          100 * min(arg1$r2_pair_raw), 100 * max(arg1$r2_pair_raw)),
  "across the three terms). This near-collinearity arises because bilateral",
  "temperature differentials are quasi-permanent: the within-pair standard",
  sprintf("deviation of $\\Delta_k$ is at most %.1f\\%% of its total standard deviation",
          100 * max(arg1$sd_within_ratio)),
  sprintf("(%.1f\\%% for $\\Delta_1$, rising to %.1f\\%% for $\\Delta_%d$),",
          100 * arg1$sd_within_ratio[1],
          100 * arg1$sd_within_ratio[nrow(arg1)], nrow(arg1)),
  "reflecting the geographic stability of cross-country temperature gaps.",
  sprintf("Product$\\times$year fixed effects, by contrast, absorb at most %.2f\\%%",
          100 * max(arg1$r2_pyr_raw)),
  "of $\\Delta_k$ variance, preserving the between-pair variation.", "",
  "\\paragraph{Gravity controls absorb observable confounders.}",
  "The main threat to cross-sectional identification is that country pairs with",
  "large temperature differentials differ on dimensions beyond climate.",
  "Table~\\ref{tab:ctrl_absorption} shows that the gravity controls",
  sprintf("absorb %.1f--%.1f\\%% of the between-pair climate variation",
          100 * min(arg2$r2_absorbed_raw), 100 * max(arg2$r2_absorbed_raw)),
  "conditional on product$\\times$year fixed effects. The normalisation of",
  "bilateral flows by domestic absorption $X_t^k(rr)$ eliminates importer-side",
  "multilateral resistance \\citep{head2014gravity}.", "",
  "\\paragraph{Coefficient stability.}",
  "Figure~\\ref{fig:stability_controls} compares the estimated dose-response",
  "curves with and without gravity controls.",
  sprintf("The two curves are highly correlated ($\\rho = %.2f$),", corr_stability),
  "confirming that the shape of $g(\\cdot)$ is not driven by observable",
  "geographic confounders. Following \\citet{oster2019unobservable}, stability",
  "to observable controls suggests that unobservable confounders of similar",
  "importance would not materially alter the estimated dose-response.", "",
  "\\begin{figure}[H]", "\\centering",
  "\\includegraphics[width=0.85\\textwidth]{stability_controls_prod.pdf}",
  "\\caption{Dose-response stability: with vs.\\ without gravity controls.",
  "Polynomial degree 3, winsorization [1\\%, 1\\%].}",
  "\\label{fig:stability_controls}", "\\end{figure}", "",
  "\\input{output/tab_fe_absorption.tex}",
  "\\input{output/tab_ctrl_absorption.tex}"),
  file.path(DIR_OUT, "appendix_identification_trade.tex"))
cat("[write] appendix_identification_trade.tex + les deux tables\n")

##################################################################################
# 8. STAGE 3 -- the grid of 12 specifications, two arms
##################################################################################

estimate_grid <- function(is_annual) {
  res <- list()
  for (cfg in SPEC_CONFIGS) for (w in WINSOR_CONFIGS) {
    key <- sprintf("%s%d_wl%swu%s", cfg$spec, cfg$degree, w$lower, w$upper)
    cat(sprintf("  Stage 3 [%s] %s\n", if (is_annual) "annual" else "daily", key))
    dp <- prepare_diffs(BACI, climate, cfg$spec, cfg$degree, w$lower, w$upper, is_annual,
                        scale = SCALE_PROD)
    dp <- stage3_climate_filter(dp, cfg$spec, cfg$degree, is_annual)   # lot 41
    ## controls EXPLICITE. ⛔ Cet appel s'en remettait au defaut de run_stage3(),
    ## qui vaut character(0), et estimait donc la grille SANS controles
    ## gravitaires -- alors que le pipeline racine l'estimait avec (defaut
    ## `controls = gravity_controls`, drf_prod_final.r l.577, applique aux deux
    ## grilles l.1606 et l.2068), et que 01c les passe explicitement a ses huit
    ## sites, dont l.672 qui est CETTE grille. Le meme objet etait estime de deux
    ## facons dans le meme paquet, et le manuscrit met les deux tables face a
    ## face (l.555-565).
    ## Mesure du 2026-09-04 sur spline4_wl0.5wu0.5, bras daily : Delta_north
    ## +22,9 % sans les controles, ecart d'AIC 137 502 -- ils ne sont pas
    ## absorbes par les effets fixes produit x annee. Corrige le 2026-09-04.
    r <- run_stage3(dp, cfg$spec, cfg$degree, is_annual,
                    controls = TRADE_STAGE3_CONTROLS)
    stage3_nobs_check(r, dp, key)
    r$wval_lower <- w$lower; r$wval_upper <- w$upper
    res[[key]] <- r
  }
  res
}
cat("=== Stage 3 ===\n")
results_daily  <- estimate_grid(FALSE)
results_annual <- estimate_grid(TRUE)

raw_table <- function(res) {
  rbindlist(lapply(names(res), function(nm) {
    r <- res[[nm]]
    data.table(model_name = nm, spec = r$spec, degree = as.integer(r$degree),
               wval_lower = r$wval_lower, wval_upper = r$wval_upper,
               p_joint = as.numeric(r$wald_joint$p),      # RAW, not rounded
               p_nonlin = as.numeric(r$wald_nonlin$p),    # RAW, not rounded
               stat_joint = as.numeric(r$wald_joint$stat),
               stat_nonlin = as.numeric(r$wald_nonlin$stat),
               aic = as.numeric(AIC(r$model)),
               nobs = tryCatch(as.integer(nobs(r$model)), error = function(e) NA_integer_))
  }))
}
raw_daily  <- raw_table(results_daily)
raw_annual <- raw_table(results_annual)
fwrite(raw_daily,  file.path(DIR_OUT, "prod_joint_tests_daily.csv"))
fwrite(raw_annual, file.path(DIR_OUT, "prod_joint_tests_annual.csv"))

## The AIC of the two arms is comparable only if N is identical: same data, same
## fixed effects, same controls, with only the exposure metric changing. If N
## differs the gap stops being interpretable, and that is said rather than
## rounded away.
cat(sprintf("[AIC-CMP] N daily = %s | N annual = %s | %s\n",
            paste(unique(raw_daily$nobs), collapse = "/"),
            paste(unique(raw_annual$nobs), collapse = "/"),
            if (isTRUE(all.equal(unique(raw_daily$nobs), unique(raw_annual$nobs))))
              "comparables" else "⚠ NON COMPARABLES"))

##################################################################################
# tab_robustness_prod_joint.tex N'EST PLUS EMISE ICI (18/09/2026, lot 10).
#
# CE QUI A CHANGE. La table dit « across specifications ». Ce script estime
# SPEC_CONFIGS, c'est-a-dire la grille du CORPS -- quatre formes x trois
# winsorisations, douze specifications. La grille du CONTRASTE, celle que le
# manuscrit appelle la grille, en compte QUINZE : GRID_SPECS y ajoute les trois
# polynomes de degre 4. La table publiait donc douze lignes sous un titre qui en
# annoncait quinze, et aucun controle ne pouvait le voir : les douze lignes
# etaient justes, et leur nombre n'etait ecrit nulle part.
#
# OU ELLE EST EMISE DESORMAIS : 04_contrast/ct4_emit_tab_contribA_grid.R, qui
# tourne APRES ct1 et dispose donc des quinze ajustements. Rien n'y est
# re-estime : ct1 appelle le meme prepare_diffs() puis le meme run_stage3() que
# ci-dessus, aux memes controles gravitaires et au meme clustering
# bidimensionnel, et ct4 verifie VALEUR PAR VALEUR que les douze specifications
# communes y reproduisent celles emises ici.
#
# CE QUI RESTE ICI, ET POURQUOI. Les deux CSV prod_joint_tests_{daily,annual}.csv
# continuent d'etre ecrits juste au-dessus. Ils ne sont plus la source de la
# table : ils en sont l'ANCRE. ct4 refuse d'emettre si l'un des deux manque, et
# refuse si un seul chiffre des douze lignes communes s'en ecarte. Les supprimer
# rendrait le controle de non-regression insatisfaisable -- c'est-a-dire qu'il
# passerait au orange, jamais au vert.
##################################################################################
cat("[note] tab_robustness_prod_joint.tex : emise par 04_contrast/",
    "ct4_emit_tab_contribA_grid.R sur les 15 specifications de GRID_SPECS.\n",
    sep = "")

##################################################################################
# 9. MONTE CARLO -- two layers of uncertainty
#   layer 1: the specification is drawn according to the Akaike weights;
#   couche 2 : beta* ~ N(beta_hat, V_hat).
##################################################################################

mc_daily  <- mc_run(results_daily,  climate, MC_N_SIM, MC_X_GRID,
                    MC_SEED[["prod_daily"]],  FALSE, "Productivity-Daily")
mc_annual <- mc_run(results_annual, climate, MC_N_SIM, MC_X_GRID,
                    MC_SEED[["prod_annual"]], TRUE,  "Productivity-Annual")

s_daily   <- mc_summarise(mc_daily)
s_annual  <- mc_summarise(mc_annual)
mc_export(s_daily,  file.path(DIR_OUT, "MC_damage_prod_daily_quantiles.csv"))
mc_export(s_annual, file.path(DIR_OUT, "MC_damage_prod_annual_quantiles.csv"))

## Central damage functions, weighted by the Akaike weights. They are inputs to
## the GE and to 01e_diagnostics.R, which applies them to the ERA5 distributions
## to compute eta and the non-stationarity.
## BOTH arms are emitted. The annual arm is not a duplicate: it carries the
## metric the paper compares the daily one against, and without it the gap in
## eta between the two cannot be measured.
for (arm in c("daily", "annual")) {
  res <- if (arm == "daily") results_daily else results_annual
  avg <- aic_curve(res, climate, is_annual = (arm == "annual"))
  setnames(avg, c("x", "y_avg"), c("temperature_C", "f_z_raw"))
  f <- file.path(DIR_OUT, sprintf("damage_prod%s_bin.csv",
                                  if (arm == "daily") "" else "_annual"))
  fwrite(avg, f)
  cat(sprintf("[write] %s (%d points, [%.1f, %.1f] °C)\n",
              basename(f), nrow(avg), min(avg$temperature_C),
              max(avg$temperature_C)))
}

##################################################################################
# 9bis. LES DOUZE COURBES, UNE PAR SPECIFICATION
#
# POURQUOI CE BLOC EXISTE (09/09/2026).
# `damage_prod_bin.csv` ci-dessus est la MOYENNE ponderee des douze. Deux
# consommateurs ont besoin des douze SEPAREMENT, et ni l'un ni l'autre ne peut
# les reconstruire depuis la moyenne :
#   - proxy_pe_signtest.r integre CHAQUE courbe contre une distribution de pays,
#     puis repondere : integrer la moyenne donnerait un autre objet des que la
#     distribution n'est pas la meme d'un pays a l'autre ;
#   - le controle A2 de ct1_delta_north_grid.R (l.653) confronte sa propre
#     re-estimation aux douze dumps. Sans eux il s'inscrivait « NON EVALUEE —
#     aucune etape du paquet ne produit output/drf_prod_specs/ », c'est-a-dire
#     que la non-regression du canal productivite n'etait pas mesuree.
#
# AUCUNE ESTIMATION NEUVE. `results_daily` et `aic_weights()` sont deja en
# memoire ; ce bloc les ecrit.
#
# LES NOMS SONT CEUX QUE LES LECTEURS ATTENDENT. `temperature_C` / `f_z_raw`
# sont les colonnes de damage_prod_bin.csv, et `model_name` la cle de
# aic_weights_Productivity_Daily.csv. Le nom de fichier suit le motif
# `^drf_prod_spec_\\d+\\.csv$` que ct1 l.664 filtre.
#
# CHAQUE COURBE EST NORMALISEE A SON PROPRE MAXIMUM -- c'est la convention de
# get_effect_curve, et ct1 l.676 la recale explicitement. Elle est SANS EFFET
# sur proxy_pe_signtest.r, qui ne lit que des DIFFERENCES (projection moins
# base) sur une meme courbe : une constante additive s'y annule.
##################################################################################

{
  dumpdir <- file.path(DIR_OUT, "drf_prod_specs")
  dir.create(dumpdir, showWarnings = FALSE, recursive = TRUE)
  aw_d <- aic_weights(results_daily)
  n_dump <- 0L
  for (i in seq_along(results_daily)) {
    nm <- names(results_daily)[i]
    cc <- get_effect_curve(results_daily[[nm]], climate)
    if (is.null(cc)) {
      ## Un mecanisme qui ne fait rien doit dire qu'il n'a rien fait.
      cat(sprintf("  [!] %s : aucune courbe (coefficients absents) — NON dumpee\n", nm))
      next
    }
    out <- cc[, .(temperature_C = round(x, 1), f_z_raw = y)]
    out[, `:=`(model_name = nm,
               aic_weight = aw_d$weight[aw_d$model_name == nm],
               spec       = results_daily[[nm]]$spec,
               degree     = results_daily[[nm]]$degree,
               wval_lower = results_daily[[nm]]$wval_lower,
               wval_upper = results_daily[[nm]]$wval_upper)]
    fwrite(out, file.path(dumpdir, sprintf("drf_prod_spec_%02d.csv", i)))
    n_dump <- n_dump + 1L
  }
  ## Le compte est un RESULTAT, pas une cible : il vaut 12 quand les douze
  ## specifications rendent une courbe, et moins si l'une d'elles n'a pas de
  ## coefficients. Il est imprime pour que le cas degrade se lise au log.
  cat(sprintf("[write] drf_prod_specs/ : %d courbes sur %d specifications\n",
              n_dump, length(results_daily)))
  ## Les poids sont ecrits AVEC les courbes : un consommateur qui reponderait
  ## depuis un autre fichier pourrait apparier sur un millesime different.
  cat(sprintf("        poids AIC : somme %.6f | max %.3f\n",
              sum(aw_d$weight), max(aw_d$weight)))
  rm(dumpdir, aw_d, n_dump)
}

##################################################################################
# 10. PAPER FIGURES -- the two fan charts
#
# Redrawn from the quantile CSV files, with no re-estimation:
#   - the support is DERIVED from the artefact by plateau detection, never coded;
#   - the two panels share a window and a scale;
#   - the shift is SCALAR: the median is zero at the reference temperature, but
#     the band is not pinched there, because each draw was normalised to its own
#     maximum. The axis says "MC median = 0", never "normalised".
##################################################################################

## LA FIGURE DU PAPIER A CHANGE DE PRODUCTEUR (2026-09-17, lot 8).
## Elle porte desormais, en fond, la courbe estimee de chaque specification de
## la grille -- un objet archive par ct1, donc par une etape POSTERIEURE a
## celle-ci. La produire ici la produirait sans son fond, et sa legende
## annoncerait vingt-quatre courbes qu'elle ne montrerait pas.
## Producteur unique : 04_contrast/ct6_fig_drf_paper.R, qui relit
## climate_long_annual.rds ecrit quelques centaines de lignes plus haut.
## Le canal AMENITE n'est pas concerne : 01b appelle toujours plot_fans.
cat("\n[note] MC_prod_{daily,annual}_fan_paper.pdf : emis par",
    "04_contrast/ct6_fig_drf_paper.R (apres ct1).\n")
cat("\n=== 03a_productivity_drf.R terminé ===\n")
