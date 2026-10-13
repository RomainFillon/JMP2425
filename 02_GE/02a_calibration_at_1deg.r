####################################################################################################
# 02a_calibration_at_1deg.r -- builds the one-degree gridded inputs of the GE.
#
# Produces the initial state and the exogenous paths the solver needs:
#   initial population, wage and productivity
#   the indices of the cells kept, those positive in all three variables
#   lambda0, the initial trade shares
#   exogenous population and productivity paths
#   descriptive statistics on initial productivity and migration
#
# Reads:  data_inputs/ (GADM, GHSL, Kummu, ILO, SSP2, Wittgenstein, World Bank),
#         the Cruz and Rossi-Hansberg archive via lib_eggw.R
# Writes: intermediate_outputs/identifier_xy_final.csv, indices_drop.csv,
#         population0.csv, wage0.csv, productivity0.csv, lambda.h5,
#         data_final_{population,productivity}.csv,
#         output/myplot_migration.pdf, output/no_activity_map_smooth_light.pdf
#
# THIS SCRIPT RENUMBERS THE CELL IDS. identifier_xy_final.csv assigns a new id to
# every cell, and the other outputs here are keyed by that numbering. Any file
# keyed by the OLD ids becomes silently wrong: the id vectors still match, so no
# check based on ids alone reports it.
# In particular, 00_temperature_projections_JMP.R must be re-run after this
# script, or every cell will receive the temperature of an unrelated cell.

rm(list=ls())
libraries <- c("rhdf5","akima","hdf5r","smoothr","spData","tidyterra","DescTools",
               "rnaturalearth","fuzzyjoin","patchwork","ggplot2","data.table",
               "R.matlab","readxl","sf","raster","Matrix","dplyr","tidyr",
               "terra","countrycode","exactextractr","bigmemory","RcppCNPy")
lapply(libraries, require, character.only = TRUE)

winsor_choice_main <- 0.05    # population: lower tail only; productivity and wage: both tails
winsor_choice_plot <- 0.005   # plots only

# Paths
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
## PROVENANCE check, when the CRH package is present. Read from INSIDE the zip:
## nothing is unpacked, and the package is not copied into the repository.
identifier <- try({
  con <- unz(file.path(DIR_IN, "EGGW_Replication_Package.zip"), file.path("EGGW_Replication_Package", "Data/HDI_GDPpc.csv"))
  g <- data.table::fread(text = readLines(con, warn = FALSE))
  close(con)
  unique(g[, .(y = lat_17048, x = lon_17048)])
}, silent = TRUE)

identifier <- copy(identifier)[, id := seq_len(.N)]
write.csv(identifier, file.path(DIR_TMP, "identifier_xy_final.csv"), row.names=FALSE)

####################################################################################################
# 1. LOAD BASE GRID
####################################################################################################
identifier <- read.csv(file.path(DIR_TMP, "identifier_xy_final.csv"))
points     <- vect(identifier, geom=c("x","y"), crs="EPSG:4326")
r_template <- rast(ext(-180,180,-90,90), resolution=1, crs="EPSG:4326")
r          <- rasterize(points, r_template, field="id", fun=mean)

####################################################################################################
# 2. LOAD COUNTRIES AND RASTERIZE
####################################################################################################
countries      <- st_read(file.path(DIR_IN, "gadm36_shp/gadm_simplified.shp"))
countries      <- st_transform(countries, crs(r))
country_raster <- rasterize(vect(countries), r, field="GID_0")
country_list   <- unique(countries$GID_0)

####################################################################################################
# 3. LOAD RAW RASTERS (no winsorization yet)
####################################################################################################
# GDP, gridded at one degree
GDP_gridded <- rast(file.path(DIR_IN, "Kummu/rast_gdpTot_1990_2022_30arcmin.tif"))
GDP_gridded <- aggregate(GDP_gridded, fact=2, fun=sum, na.rm=TRUE)
GDP_gridded <- resample(GDP_gridded, country_raster, method="near")
GDP_gridded <- GDP_gridded[["gdp_tot_2015"]]

# Population, gridded at one degree
population_gridded <- rast(file.path(DIR_IN, "GHSL_POP/tif/GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif"))
population_gridded <- aggregate(population_gridded, fact=120, fun=sum, na.rm=TRUE)
population_gridded <- resample(population_gridded, country_raster, method="near")
population_all=as.data.frame(population_gridded, xy=TRUE)
colnames(population_all)=c("x","y","value")
write.csv(population_all, file.path(DIR_TMP, "population_all.csv"), row.names=FALSE)
#population_all <- read.csv(file.path(mydirection_data, "population_all.csv"))

# Land area (for density filter)
land     <- rast(file.path(DIR_IN, "GPW/gpw_v4_land_water_area_rev11_landareakm_1_deg.tif"))
area_r   <- cellSize(r, unit="km")
area_r   <- as.data.frame(area_r, xy=TRUE)
colnames(area_r) <- c("x","y","cell_area_km2")

####################################################################################################
# 4. POPULATION WEIGHTS PER COUNTRY PER ONE-DEGREE CELL
####################################################################################################
## ⛔ REMPLACE LE 28/09/2026 (lot 51). Cette section faisait `exact_extract` de la
## population DEJA AGREGEE AU DEGRE sur les 339 127 polygones de GADM 3.6, sans
## employer `coverage_area` : chaque polygone touchant une cellule recevait la
## population ENTIERE de la cellule, et la « part » d'un pays etait la part de ses
## POLYGONES (coverage_id / population = entier dans 15 402 cellules sur 15 402 ;
## RECAP 3.230). Les poids viennent desormais de 00d_poids_population.py :
## population de l'INTERSECTION cellule x pays (GHSL 30'' dans gadm_410 non
## simplifie, periode 2015), normalisee PAR CELLULE.
source(file.path(.root, "02_GE", "lib_poids.R"))
lengths <- fread(file.path(DIR_TMP, "poids_par_cellule.csv"), na.strings = NULL)
lengths <- lengths[GID_0 != ""]
temoin_part_de_la_cellule(lengths)
lengths[, `:=`(x = round(x, 1), y = round(y, 1))]
lengths <- merge(as.data.table(identifier)[, .(id, x = round(x, 1), y = round(y, 1))],
                 lengths[, .(x, y, GID_0, pop, part_de_la_cellule)], by = c("x", "y"))
lengths <- as.data.frame(lengths[, .(id, GID_0, pop, part_de_la_cellule)])
cat(sprintf("Poids : %d couples (cellule, pays), %d cellules, %d pays\n",
            nrow(lengths), length(unique(lengths$id)), length(unique(lengths$GID_0))))

####################################################################################################
# 5. LABOR SHARE FROM ILO (with regional imputation)
####################################################################################################
labor_share_raw <- read.csv(file.path(DIR_IN, "ILO_labor_share_income/LAP_2GDP_NOC_RT_A-20240922T2147.csv"))
labor_share_raw <- labor_share_raw[, c("ref_area.label","obs_value","time")]
colnames(labor_share_raw) <- c("country","share_income","time")
labor_share_raw <- labor_share_raw[labor_share_raw$time == 2015, ]

# Manual name harmonization (ILO -> WB naming)
labor_share_raw$country[labor_share_raw$country == "Côte d'Ivoire"]                                <- "Cote d'Ivoire"
labor_share_raw$country[labor_share_raw$country == "Congo, Democratic Republic of the"]            <- "Democratic Republic of the Congo"
labor_share_raw$country[labor_share_raw$country == "Hong Kong, China"]                             <- "China, Hong Kong SAR"
labor_share_raw$country[labor_share_raw$country == "Macao, China"]                                 <- "China, Macao SAR"
labor_share_raw$country[labor_share_raw$country == "Korea (the Democratic People's Republic of)"]  <- "Dem. People's Rep. of Korea"
labor_share_raw$country[labor_share_raw$country == "Occupied Palestinian Territory"]               <- "Palestine"
labor_share_raw$country[labor_share_raw$country == "Taiwan, China"]                                <- "China, Taiwan Province of China"
labor_share_raw$country[labor_share_raw$country == "Tanzania, United Republic of"]                 <- "United Republic of Tanzania"
labor_share_raw$country[labor_share_raw$country == "Channel Islands"]                              <- NA

# Remove regional aggregates
regional_aggregates <- c("World","Africa","Americas","Caribbean","ASEAN","BRICS",
                         "G20","G7","MENA","APEC","CARICOM","Arab League",
                         "Arab States","European Union 28","European Union 27",
                         "Northern America","South America","Central America",
                         "Eastern Asia","Southern Asia","Western Asia","Central Asia",
                         "Northern Europe","Southern Europe","Western Europe","Eastern Europe",
                         "Southern Africa","Eastern Africa","Western Africa",
                         "Central Africa","Northern Africa","Sub-Saharan Africa",
                         "Pacific Islands","South-Eastern Asia",
                         "Latin America and the Caribbean",
                         "Europe and Central Asia","Asia and the Pacific",
                         "Central and Western Asia",
                         "Northern, Southern and Western Europe",
                         "South-Eastern Asia and the Pacific")

labor_share_raw <- labor_share_raw %>%
  filter(!country %in% regional_aggregates,
         !grepl("income|excluding|BRICS|World", country, ignore.case=TRUE),
         !is.na(country))

cat("Rows after cleaning:", nrow(labor_share_raw), "\n")

# Join to ISO3
## CORRIGE (lot 60). Le fichier OIT ne porte AUCUN code pays, seulement
## ref_area.label. L'appariement se faisait par egalite stricte avec les noms de
## total_UN_code.csv : 7 pays dont la part est OBSERVEE (CIV, HKG, MAC, NLD, PRK,
## PSE, TWN) sortaient avec GID_0 = NA, puis recevaient une moyenne regionale.
## Le code ISO3 est desormais construit par le dictionnaire de countrycode
## (country.name -> iso3c) : les 181 pays deja apparies gardent le meme code, les
## 7 sont gagnes, aucun ne reste sans code, aucun code en double.
labor_share_raw$GID_0 <- countrycode(labor_share_raw$country, "country.name", "iso3c")
temoin_code_source(labor_share_raw, "part du revenu du travail (OIT, 2015)")

# Regional averages (for imputation)
labor_share_raw <- labor_share_raw %>%
  mutate(region = countrycode(GID_0, origin="iso3c", destination="region"))

labor_share_region <- labor_share_raw %>%
  filter(!is.na(share_income), !is.na(region)) %>%
  group_by(region) %>%
  summarise(share_income_region = mean(share_income, na.rm=TRUE), .groups="drop")

cat("Regional averages:\n")
print(as.data.frame(labor_share_region))

# Impute missing share_income by region
labor_share_raw <- labor_share_raw %>%
  left_join(labor_share_region, by="region") %>%
  mutate(share_income = ifelse(is.na(share_income), share_income_region, share_income))

# Handle GID_0 in lengths not covered by ILO at all
global_mean <- mean(labor_share_region$share_income_region, na.rm=TRUE)

missing_gid <- lengths %>%
  distinct(GID_0) %>%
  anti_join(labor_share_raw %>% filter(!is.na(GID_0)), by="GID_0") %>%
  mutate(region = countrycode(GID_0, origin="iso3c", destination="region")) %>%
  left_join(labor_share_region, by="region") %>%
  mutate(share_income = case_when(
    !is.na(share_income_region)      ~ share_income_region,
    GID_0 == "XKO"                   ~ labor_share_region$share_income_region[
                                          labor_share_region$region == "Europe & Central Asia"],
    GID_0 == "MYT"                   ~ labor_share_region$share_income_region[
                                          labor_share_region$region == "Sub-Saharan Africa"],
    GID_0 %in% c("XCA","XNC","XAD") ~ global_mean,
    TRUE                             ~ global_mean
  )) %>%
  dplyr::select(GID_0, share_income)

cat("Additional GID_0 imputed:", nrow(missing_gid), "\n")
temoin_imputation_evitable(missing_gid$GID_0, labor_share_raw$GID_0,
                           "part du revenu du travail (OIT, 2015)")

# Combine and aggregate to id level
labor_share_pays <- bind_rows(
    labor_share_raw %>% dplyr::select(GID_0, share_income),
    missing_gid
  )
labor_share <- labor_share_pays %>%
  merge(lengths, by="GID_0") %>%
  group_by(id) %>%
  summarise(share_income = sum(share_income * part_de_la_cellule, na.rm=TRUE) / 100,
            .groups="drop")
## Garde (lot 51) : la part du travail d'une cellule est la repartition
## proportionnelle a la population des parts nationales.
identite_aller_retour("part du revenu du travail (OIT, 2015)", lengths,
                      transform(labor_share_pays, r = share_income / 100),
                      transform(labor_share, v = share_income))

# Diagnostics
cat("Remaining NA:", sum(is.na(labor_share$share_income)), "\n")
cat("Remaining zeros (should all have zero population):",
    sum(labor_share$share_income == 0, na.rm=TRUE), "\n")

####################################################################################################
# 6. BUILD PRODUCTIVITY AND WAGE (raw, before winsorization)
####################################################################################################
# Density filter: zero out low-density cells
data_prod <- as.data.frame(GDP_gridded, xy=TRUE)
data_prod <- left_join(identifier, data_prod, by=c("x","y"))
data_pop  <- as.data.frame(population_gridded, xy=TRUE)
data_pop  <- left_join(identifier, data_pop, by=c("x","y"))
colnames(data_pop)[colnames(data_pop) == names(data_pop)[ncol(data_pop)]] <- "Population"
colnames(data_prod)[colnames(data_prod) == "gdp_tot_2015"] <- "GDP"

data3 <- left_join(data_prod, data_pop[, c("id","Population")], by="id")
data3 <- left_join(data3, area_r, by=c("x","y"))
data3$density    <- data3$Population / data3$cell_area_km2

####################################################################################################
# 7. WINSORIZE EACH VARIABLE ON ITS OWN DISTRIBUTION (once, on final variable)
####################################################################################################
# Population (post density filter)
# winsorized on non-null 
pop_final        <- data3[, c("id","Population","GDP")]
pop_final$Population[pop_final$Population<1] = 0

pop_bounds       <- quantile(pop_final$Population[pop_final$Population > 0],
                             probs=c(winsor_choice_main), na.rm=TRUE)
pop_final$value <- pop_final$Population
mask_nonzero <- pop_final$Population > 0

pop_final$prod       <- pop_final$GDP / pop_final$value  
prod_final = pop_final[c("id","prod")]
prod_final$prod = ifelse(is.finite(prod_final$prod), prod_final$prod, 0)
prod_bounds       <- quantile(prod_final$prod[prod_final$prod > 0 & is.finite(prod_final$prod)],
                             probs=c(winsor_choice_main, 1-winsor_choice_main), na.rm=TRUE)
prod_final$value <- prod_final$prod
mask_nonzero <- prod_final$value > 0
prod_final$value[mask_nonzero] <- pmin(pmax(prod_final$value[mask_nonzero], prod_bounds[1]), prod_bounds[2])

# Raw wage = productivity * labor share
wage_final <- left_join(prod_final[, c("id","prod")], labor_share, by="id")
wage_final$wage <- wage_final$prod * wage_final$share_income
wage_final$wage <- ifelse(is.na(wage_final$wage), 0, wage_final$wage)

wage_final = wage_final[c("id","wage")]
wage_final$wage = ifelse(is.finite(wage_final$wage), wage_final$wage, 0)
wage_bounds       <- quantile(wage_final$wage[wage_final$wage > 0 & is.finite(wage_final$wage)],
                             probs=c(winsor_choice_main, 1-winsor_choice_main), na.rm=TRUE)
wage_final$value <- wage_final$wage
mask_nonzero <- wage_final$value > 0
#no winsoring for low wages
wage_final$value[mask_nonzero] <- pmin(wage_final$value[mask_nonzero], wage_bounds[2])

####################################################################################################
# 8. INDICES TO KEEP (non-zero in all three variables)
####################################################################################################
active_ids <- Reduce(intersect, list(
  prod_final$id[prod_final$value > 0 & !is.na(prod_final$value) & is.finite(prod_final$value)],
  wage_final$id[wage_final$value > 0 & !is.na(wage_final$value) & is.finite(wage_final$value)],
  pop_final$id[pop_final$value > 1 & !is.na(pop_final$value) & is.finite(pop_final$value)]))

idrop <- data.frame(x = sort(active_ids))
write.csv(idrop, file.path(DIR_TMP, "indices_drop.csv"), row.names=FALSE)
cat("Number of active cells:", nrow(idrop), "out of 17048\n")

####################################################################################################
# 9. SAVE WAGE, POPULATION, PRODUCTIVITY (same format)
####################################################################################################
# wage0.csv
wage_mean <- mean(wage_final$value[wage_final$id %in% active_ids & wage_final$value > 0], na.rm=TRUE)
wage_out  <- wage_final[, c("id","value")]
wage_out$value <- wage_out$value / wage_mean
colnames(wage_out) <- c("id","value")
write.csv(wage_out, file.path(DIR_TMP, "wage0.csv"), row.names=FALSE)
cat("wage_mean used for normalization:", wage_mean, "\n")

# population0.csv
pop_out <- pop_final[, c("id","value")]
colnames(pop_out) <- c("id","value")
write.csv(pop_out, file.path(DIR_TMP, "population0.csv"), row.names=FALSE)

# productivity0.csv
prod_out <- prod_final[, c("id","value")]
colnames(prod_out) <- c("id","value")
write.csv(prod_out, file.path(DIR_TMP, "productivity0.csv"), row.names=FALSE)

####################################################################################################
# 10. TRADE SHARES -> HDF5
####################################################################################################
EGGW_ZIP    <- file.path(DIR_IN, "EGGW_Replication_Package.zip")
EGGW_PREFIX <- "EGGW_Replication_Package"

## The .mat readers need a REAL file: hdf5r and R.matlab open a path and
## cannot read from a connection. Only that one file is extracted, once, and
## cached. The CSV files of the same package are read directly in the archive.
eggw_mat <- function(nm) {
  dst <- file.path(DIR_TMP, "eggw_mat", nm)
  if (!file.exists(dst)) {
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  extraction de %s...\n", nm))
    unzip(EGGW_ZIP, files = file.path(EGGW_PREFIX, "Data", nm),
          exdir = dirname(dst), junkpaths = TRUE)
    if (!file.exists(dst))
      stop("⛔ ", nm, " introuvable dans l'archive — arborescence modifiée ?")
  }
  dst
}

#trmult_reduced.mat is already ^6.5
file_tau <- eggw_mat("trmult_reduced.mat")   # already raised to the power theta

h5_tau    <- H5File$new(file_tau, mode="r")
dset      <- h5_tau[["trmult_reduced"]]
n         <- dset$dims[1]

theta <- 6.5
mu    <- 0.6
# Normalisation before lambda is computed. Productivity is divided by its mean
# over the active locations, so Z has mean one. wage0.csv is already normalised
# by the wage mean of the active locations, so it is not divided again here.
# Normalising the wage a second time gives a different lambda0, so the choice
# is not neutral.

Z_mean <- mean(prod_out$value[prod_out$value > 0], na.rm=TRUE)
Z      <- prod_out$value / Z_mean   # normalised, mean one

x_vec <- wage_out$value^mu

lambda_h5_file <- file.path(DIR_TMP, "lambda.h5")
if (file.exists(lambda_h5_file)) file.remove(lambda_h5_file)
h5_out <- H5File$new(lambda_h5_file, mode="w")
h5_out$create_dataset("lambda",
                       dims    = c(n, n),
                       dtype   = h5types$H5T_NATIVE_FLOAT,
                       chunk_dims = c(100, n))   # chunked by row, for fast reads

for (row_r in 1:n) {
  tau_neg_theta <- dset[row_r, ]
  tau_neg_theta[is.na(tau_neg_theta) | is.infinite(tau_neg_theta)] <- 0.0
  
  ## 🚨 CORRIGE LE 23/09/2026 (lot 25) : Z ENTRE A LA PUISSANCE theta.
  ## Le solveur traite la productivite comme z, avec l'echelle de Frechet
  ## T = z^theta (forward.py l.97 : (z/x)**param_theta), et l'eq. (5) du
  ## papier identifie theta*f_z. Ecrire Z lineaire ici faisait de ce proxy
  ## l'echelle T elle-meme : exposant net de Z * x^-theta sur le PIB/tete
  ## mesure a -3,36 (theorique 1 - 0,6*6,5 = -2,9), donc un gradient INVERSE --
  ## les cellules riches vendaient 3 % de leur revenu (RECAP §3.203).
  ## Sous la forme corrigee, la grille calibree est largement domestique
  ## (part importee moyenne 0,43 %) ; c'est une propriete du cadre de DNRH avec
  ## cette matrice de couts, verifiee sur LEURS fondamentaux (§3.208).
  ## Reconstruction independante : output_lambda0_probe/rebuild_lambda0.py.
  raw <- Z^theta * (x_vec)^(-theta) * tau_neg_theta
  raw[is.na(raw) | is.infinite(raw)] <- 0.0
  
  denom <- sum(raw)
  
  if (denom == 0) {
    lambda_row        <- rep(0, n)
    lambda_row[row_r] <- 1.0
  } else {
    lambda_row <- raw / denom
  }
  
  h5_out[["lambda"]][row_r, ] <- lambda_row
  if (row_r %% 500 == 0) cat("Done row", row_r, "/", n, "\n")
}

h5_tau$close_all()
h5_out$close_all()
cat("lambda.h5 saved.\n")

h5_read <- H5File$new(lambda_h5_file, mode="r")
dset_lambda <- h5_read[["lambda"]]

active_indices <- sort(active_ids)
wage_active <- wage_out$value[match(active_indices, wage_out$id)]
pop_active  <- pop_out$value[match(active_indices, pop_out$id)]
income_vec  <- wage_active * pop_active

cat("wage normalisé: min=", min(wage_active[wage_active>0]),
    "median=", median(wage_active[wage_active>0]),
    "mean=", mean(wage_active[wage_active>0]), "\n")

# Test on five locations only
set.seed(42)
sample_k <- sample(seq_along(active_indices), 5)

w_implied_sample <- numeric(5)
for (k in seq_along(sample_k)) {
  r_global <- active_indices[sample_k[k]]
  
  # Read the whole row r at once, which is much faster
  lambda_row_r <- dset_lambda[r_global, ]   # lambda[r, s] for every s
  
  # Keep only the active locations
  lambda_col_r <- lambda_row_r[active_indices]  # lambda[r, s] for active s
  
  # Note: here lambda[r,s] is the purchases of r from s
  # For clearing: receipts of r are sum_s lambda[s,r] * w(s) * L(s).
  # Approximated here by lambda[r,s] ~ lambda[s,r]: coarse, but fast.
  w_implied_sample[k] <- sum(lambda_col_r * income_vec, na.rm=TRUE) /
                          max(pop_active[sample_k[k]], 1e-8)
}

w_implied_norm <- w_implied_sample / mean(w_implied_sample[w_implied_sample>0])
w_obs_norm     <- wage_active[sample_k] / mean(wage_active[wage_active>0])
ratio          <- w_implied_norm / pmax(w_obs_norm, 1e-8)

cat("ratio w_implied/w_obs: median=", median(ratio),
    "min=", min(ratio), "max=", max(ratio), "\n")

h5_read$close_all()

####################################################################################################
# Population and Productivity Growth Rates (SSP2)
####################################################################################################

####################################################################################################
# 1. POPULATION GROWTH RATES (Wittgenstein SSP2)
####################################################################################################

# Load raw data
pop_raw <- read.csv(file.path(DIR_IN, "Wittgenstein_population_SSP2/PROJresult_AGE_SSP2mig0_V13.csv"))
pop_raw <- pop_raw %>%
  group_by(region, Time) %>%
  summarise(population = sum(pop, na.rm=TRUE), .groups="drop")
colnames(pop_raw) <- c("varval","year","population")

# Region codes to names
data_id <- read.csv(file.path(DIR_IN, "Wittgenstein_population_SSP2/recode_dictionary.csv"))
pop_raw <- merge(pop_raw, data_id, by="varval")
pop_raw <- pop_raw[, c("year","population","varvaldesc")]

# Names to ISO3
iso_map <- read.csv(file.path(DIR_IN, "Wittgenstein_population_SSP2/UNSD_Methodology.csv"), sep=";")
iso_map <- iso_map[, c("ISO.alpha3.Code","Country.or.Area")]
colnames(iso_map) <- c("GID_0","varvaldesc")

# Manual name harmonization
iso_map$varvaldesc[iso_map$varvaldesc == "China, Hong Kong Special Administrative Region"] <- "China, Hong Kong SAR"
iso_map$varvaldesc[iso_map$varvaldesc == "Democratic Republic of the Congo"]               <- "Dem. Republic of the Congo"
iso_map$varvaldesc[iso_map$varvaldesc == "China"]                                           <- "China, Taiwan Province of China"
iso_map$varvaldesc[iso_map$varvaldesc == "Côte d'Ivoire"]                                  <- "Cote d'Ivoire"
iso_map$varvaldesc[iso_map$varvaldesc == "Lao People's Democratic Republic"]               <- "Lao People's Dem. Republic"
iso_map$varvaldesc[iso_map$varvaldesc == "China, Macao Special Administrative Region"]     <- "China, Macao SAR"
iso_map$varvaldesc[iso_map$varvaldesc == "Democratic People's Republic of Korea"]          <- "Dem. People's Rep. of Korea"
iso_map$varvaldesc[iso_map$varvaldesc == "Netherlands (Kingdom of the)"]                   <- "Netherlands"
iso_map$varvaldesc[iso_map$varvaldesc == "Türkiye"]                                        <- "Turkey"
iso_map$varvaldesc[iso_map$varvaldesc == "United Kingdom of Great Britain and Northern Ireland"] <- "United Kingdom"
iso_map$varvaldesc[iso_map$varvaldesc == "Réunion"]                                        <- "Reunion"
iso_map$varvaldesc[iso_map$varvaldesc == "Western Sahara"]                                 <- "Morocco"

pop_raw <- merge(pop_raw, iso_map, by="varvaldesc")
colnames(pop_raw) <- c("NAME_0","year","population","GID_0")
pop_raw <- pop_raw %>%
  group_by(GID_0, year) %>%
  summarise(population = sum(population, na.rm=TRUE), .groups="drop")

# --- Disaggregate to the one-degree grid using population weights ---
## CORRIGE (lot 51, author's decision (a)). La cle etait (id, year, GID_0) : chaque
## groupe n'avait qu'une ligne, share_norm valait 1 partout, et une cellule frontaliere
## recevait la SOMME des populations nationales entieres des pays qui la touchent. La
## cle est desormais (id, year), comme la croissance de productivite plus bas.
pop_parts <- merge(lengths, pop_raw, by="GID_0") %>%
  group_by(id, year) %>%
  mutate(total_share = sum(part_de_la_cellule, na.rm=TRUE),
         share_norm  = part_de_la_cellule / total_share)
temoin_part_non_triviale(pop_parts$share_norm, "share_norm, population Wittgenstein")
pop_grid <- pop_parts %>%
  group_by(id, year) %>%
  summarise(population_weighted = sum(population * share_norm, na.rm=TRUE),
            .groups="drop")
identite_aller_retour("population Wittgenstein SSP2 (niveau)", lengths,
                      transform(pop_raw, r = population),
                      transform(pop_grid, v = population_weighted), par = "year")

# --- Growth rates (5-year steps from Wittgenstein: 2020, 2025, ...) ---
pop_growth <- pop_grid %>%
  filter(year != 2015) %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(growth_rate = 1 + (population_weighted - lag(population_weighted)) / lag(population_weighted)) %>%
  ungroup()

pop_growth$growth_rate <- ifelse(is.na(pop_growth$growth_rate), 1, pop_growth$growth_rate)

# Fill missing ids with growth = 1 (no change)
max_id <- max(identifier$id, na.rm=TRUE)
years_witt <- unique(pop_growth$year)
pop_growth <- left_join(
  expand.grid(id=1:max_id, year=years_witt, stringsAsFactors=FALSE),
  pop_growth[, c("id","year","growth_rate")],
  by=c("id","year")
) %>% mutate(growth_rate = ifelse(is.na(growth_rate), 1, growth_rate))

# --- First period 2015-2020: World Bank observed data ---
wb_pop <- read.csv(file.path(DIR_IN, "WB_pop_country/population_data.csv"))
wb_pop <- wb_pop[, c("Country.Code","X2015","X2020")]
colnames(wb_pop) <- c("GID_0","pop_2015","pop_2020")

## CORRIGE (lot 51, decision (a)) : cle (id, GID_0) -> (id), meme defaut que ci-dessus.
wb_parts <- merge(wb_pop, lengths, by="GID_0") %>%
  group_by(id) %>%
  mutate(total_share = sum(part_de_la_cellule, na.rm=TRUE),
         share_norm  = part_de_la_cellule / total_share)
temoin_part_non_triviale(wb_parts$share_norm, "share_norm, population Banque mondiale")
wb_grid <- wb_parts %>%
  group_by(id) %>%
  summarise(
    pop_2015_w = sum(pop_2015 * share_norm, na.rm=TRUE),
    pop_2020_w = sum(pop_2020 * share_norm, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(growth_rate = 1 + (pop_2020_w - pop_2015_w) / pop_2015_w,
         growth_rate = ifelse(is.na(growth_rate) | is.infinite(growth_rate), 1, growth_rate))
identite_aller_retour("population Banque mondiale 2015 (niveau)", lengths,
                      transform(wb_pop, r = pop_2015),
                      transform(wb_grid, v = pop_2015_w))

wb_out <- left_join(
  data.frame(id=1:max_id),
  wb_grid[, c("id","growth_rate")],
  by="id"
) %>% mutate(growth_rate = ifelse(is.na(growth_rate), 1, growth_rate))
colnames(wb_out) <- c("id","2020")

# --- Combine WB 2020 + Wittgenstein 2025+ and pivot wide ---
pop_growth_wide <- pop_growth %>%
  filter(year > 2020) %>%
  pivot_wider(names_from=year, id_cols=id, values_from=growth_rate)

pop_final_wide <- left_join(wb_out, pop_growth_wide, by="id")
pop_final_wide <- pop_final_wide[, -1]   # drop id, keep only growth columns

write.csv(pop_final_wide,
          file.path(DIR_TMP, "data_final_population.csv"),
          row.names=FALSE)
cat("Population growth rates saved.\n")

####################################################################################################
# 2. PRODUCTIVITY GROWTH RATES (SSP2 IAMC) # nolint
####################################################################################################

prod_raw <- read_xlsx(file.path(DIR_IN, "SSP2_productivity_growth/iamc_db.xlsx"))
prod_raw <- as.data.frame(prod_raw)

## PLUSIEURS MODELES PAR PAYS (lot 60, author's decision : on garde la moyenne,
## et on le dit). iamc_db.xlsx porte GDP|PPP en SSP2 de DEUX modeles, OECD Env-Growth
## et IIASA GDP, pour la plupart des pays ; un seul pour quelques-uns ; et un
## troisieme, PIK GDP-32, pour USA seulement. Plus bas, share_norm se divise sur
## toutes les lignes d'un (pays, annee) : la valeur d'une cellule est donc la
## MOYENNE des modeles disponibles pour son pays. Le compte par pays est ecrit
## pour l'emission.
ssp2_modeles <- as.data.table(prod_raw)[!is.na(Scenario) & !is.na(Region),
  .(n_modeles = .N, modeles = paste(sort(Model), collapse = " + ")), by = .(GID_0 = Region)]
ssp2_modeles[, dans_les_poids := GID_0 %in% lengths$GID_0]
write.csv(ssp2_modeles[order(n_modeles, GID_0)], file.path(DIR_OUT, "ssp2_modeles_par_pays.csv"),
          row.names = FALSE)
.tab <- ssp2_modeles[dans_les_poids == TRUE, .N, by = n_modeles][order(n_modeles)]
cat("[SSP2 IAMC] pays des poids par nombre de modeles :",
    paste0(.tab$n_modeles, " modele(s) : ", .tab$N, " pays", collapse = " ; "), "\n")

# Pivot to long
prod_raw <- prod_raw %>%
  pivot_longer(cols=6:ncol(prod_raw), names_to="year", values_to="value") %>%
  mutate(year = as.numeric(substr(year, 1, 4))) %>%
  filter(year != 2010) %>%
  select(GID_0=Region, year, value)

# --- Disaggregate to the one-degree grid using population weights ---
prod_parts <- merge(prod_raw, lengths, by="GID_0") %>%
  group_by(id, year) %>%
  mutate(total_share = sum(part_de_la_cellule, na.rm=TRUE),
         share_norm  = part_de_la_cellule / total_share)
temoin_part_non_triviale(prod_parts$share_norm, "share_norm, productivite SSP2")
prod_grid <- prod_parts %>%
  summarise(value = sum(value * share_norm, na.rm=TRUE), .groups="drop")
## Garde (lot 51) : repartition proportionnelle, renormalisee sur les pays
## apparies d'une cellule (SSP2 IAMC ne couvre pas tous les pays).
identite_aller_retour("productivite SSP2 IAMC", lengths,
                      transform(prod_raw, r = value),
                      transform(prod_grid, v = value), par = "year")

# --- Growth rates ---
prod_growth <- prod_grid %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(growth_rate = 1 + (value - lag(value)) / lag(value)) %>%
  filter(!is.na(growth_rate)) %>%
  ungroup()

# Fill missing ids with growth = 1
years_prod <- unique(prod_growth$year)
prod_growth <- left_join(
  expand.grid(id=1:max_id, year=years_prod, stringsAsFactors=FALSE),
  prod_growth[, c("id","year","growth_rate")],
  by=c("id","year")
) %>% mutate(growth_rate = ifelse(is.na(growth_rate), 1, growth_rate))

# Pivot wide (drop id column)
prod_final_wide <- prod_growth %>%
  pivot_wider(names_from=year, id_cols=id, values_from=growth_rate) %>%
  select(-id)

write.csv(prod_final_wide,
          file.path(DIR_TMP, "data_final_productivity.csv"),
          row.names=FALSE)
cat("Productivity growth rates saved.\n")

####################################################################################################
# Initial productivities, gridded at one degree, used to calibrate the sectoral
####################################################################################################

library(Matrix)

#load data
#load grid from Cruz & Rossi-Hansberg
identifier=read.csv(file.path(DIR_TMP,"identifier_xy_final.csv"))
points <- vect(identifier, geom=c("x", "y"), crs="EPSG:4326")
new_extent <- ext(-180, 180, -90, 90)
new_resolution <- 1
r_template <- rast(ext=new_extent, resolution=new_resolution, crs="EPSG:4326")
r <- rasterize(points, r_template, field="id", fun=mean) 

#match countries to grid 
countries <- st_read(file.path(DIR_IN, "gadm36_shp/gadm_simplified.shp"))
countries <- st_transform(countries, crs(r))
country_raster <- rasterize(vect(countries), r, field="GID_0")  
country_list <- unique(countries$GID_0)

# GDP, gridded at one degree
GDP_gridded = rast(file.path(DIR_IN, "Kummu/rast_gdpTot_1990_2022_30arcmin.tif"))
GDP_gridded = aggregate(GDP_gridded, fact=2, fun=sum, na.rm=TRUE)
GDP_gridded <- resample(GDP_gridded, country_raster, method = "near")
GDP_gridded <- GDP_gridded[["gdp_tot_2015"]]

#labour workforce country level #when 2015 not available, take closest
shareL_per_country = read.csv(file.path(DIR_IN, "ILO_labor_share_income/LAP_2GDP_NOC_RT_A-20240922T2147.csv"), sep=",")
classification = read.csv(file.path(DIR_IN, "ILO_labor_share_income/classification_country_groupings_en.csv"), sep=",")
classification=classification[,c("ISO3.Code","Country")]
colnames(classification)= c("Country.Code","Country.Name")
shareL_per_country = shareL_per_country[,c("time","obs_value","ref_area.label")]
shareL_per_country <- shareL_per_country %>%
  group_by(ref_area.label) %>%
  summarize(
    selected_value = {  
      years <- time
      values <- obs_value
      if (2015 %in% years) {
        values[years == 2015]
      } else {
        idx <- which.min(abs(years - 2015))
        values[idx]
      }
    },
    selected_year = {
      years <- time
      values <- obs_value
      if (2015 %in% years) {
        2015
      } else {
        years[which.min(abs(years - 2015))]
      }
    },
    .groups = "drop")
colnames(shareL_per_country) = c("Country.Name", "selected_value", "selected_year")
shareL_per_country=left_join(shareL_per_country,classification,by="Country.Name")

shareL_per_country <- shareL_per_country %>%
  filter(Country.Code %in% country_list)

#show areas where population and gdp = 0 
# Population, gridded at one degree
population_gridded = rast(file.path(DIR_IN, "GHSL_POP/tif/GHS_POP_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif"))
population_gridded = aggregate(population_gridded, fact=120, fun=sum, na.rm=TRUE)

#country shapes
#land area
land = rast(file.path(DIR_IN, "GPW/gpw_v4_land_water_area_rev11_landareakm_1_deg.tif"))
land_mask <- land > 0

world <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world)
world_vect <- project(world_vect, crs(land))

#identifier
identifier=read.csv(file.path(DIR_TMP,"identifier_xy_final.csv"))
#id_vect <- vect(identifier, geom = c("x", "y"), crs = crs(no_activity))
id_vect <- vect(identifier, geom = c("x", "y"), crs = crs(land_mask))

pop_aligned <- resample(population_gridded, land, method = "near")
gdp_aligned <- resample(GDP_gridded, land, method = "near")
no_activity <- (pop_aligned == 0) & land_mask
id_raster <- rasterize(id_vect, no_activity, field = 1)
no_act_in_id <- no_activity & !is.na(id_raster)
no_act_not_in_id <- no_activity & is.na(id_raster)

no_activity_bin <- ifel(no_activity, 1, NA)

# Convert to polygons (terra), then to sf
no_activity_poly <- as.polygons(no_activity_bin, dissolve = TRUE)
no_activity_sf <- st_as_sf(no_activity_poly)

# Convert world_vect to sf if it is not already
world_sf <- st_as_sf(world_vect)

library(sf)
library(ggplot2)

# Projection Robinson (plus publication)
world_sf  <- st_transform(st_as_sf(world_vect), "+proj=robin")
no_activity_sf <- st_transform(no_activity_sf, "+proj=robin")

# Combine without fully merging
no_activity_combined <- no_activity_sf |>
  st_make_valid() |>
  st_combine()

# Light smoothing
no_activity_smooth <- smooth(
  no_activity_combined,
  method = "ksmooth",
  smoothness = 2   # beaucoup plus subtil
) |>
  st_simplify(dTolerance = 10000)  # light simplification

# Projection Robinson
world_sf  <- st_transform(st_as_sf(world_vect), "+proj=robin")
no_activity_smooth <- st_transform(no_activity_smooth, "+proj=robin")

p <- ggplot() +
  geom_sf(data = world_sf, fill = "grey97", color = "grey70", linewidth = 0.15) +
  geom_sf(data = no_activity_smooth,
          fill = "#2C7FB8",
          color = NA,
          alpha = 0.9) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Land Areas with No Population or Economic Output (2015)",
    caption = "Shaded areas have been smoothed from 1° cell areas for visual clarity"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    plot.caption = element_text(
      size = 10,
      hjust = 1,          # right-aligned
      color = "grey40",
      face = "italic"
    ),
    plot.margin = margin(10, 15, 10, 15)
  )
ggsave(
  file.path(DIR_OUT,"no_activity_map_smooth_light.pdf"),
  p,
  width = 10,
  height = 5.5,
  device = cairo_pdf
)

####################################################################################################
#now try to plot migration flows wrt flows in Cruz & Rossi Hansberg
# Migration stock per one-degree band of latitude, over time
#data from Kummu is per 1000 people
bilat = read.csv(file.path(DIR_IN, "abel_cohen_migration/bilat_mig.csv"))
bilat <- bilat %>%
  group_by(year0,orig, dest)%>%
  summarise(da_pb_closed=sum(da_pb_closed, na.rm=T))

net_mig_full <- bilat %>%
  group_by(year0, country = orig) %>%
  summarise(outflow = sum(da_pb_closed, na.rm = TRUE), .groups = "drop") %>%
  right_join(
    bilat %>%
      group_by(year0, country = dest) %>%
      summarise(inflow = sum(da_pb_closed, na.rm = TRUE), .groups = "drop"),
    by = c("year0", "country")
  ) %>%
  mutate(net_migration = inflow - outflow)

world <- ne_countries(scale = "medium", returnclass = "sf")

world_iso <- world %>%
  mutate(lat_mean = st_coordinates(st_centroid(geometry))[,2]) %>%
  dplyr::select(adm0_a3, iso_a3, gu_a3, lat_mean)

world_iso_long <- world_iso %>%
  pivot_longer(cols = c(adm0_a3, iso_a3, gu_a3),
               names_to = "type_code", values_to = "country") %>%
  filter(!is.na(country)) %>%
  distinct(country, .keep_all = TRUE) %>%
  dplyr::select(country, lat_mean)
missing_codes <- setdiff(unique(net_mig_full$country), world_iso_long$country)
missing_codes
manual_map <- tibble::tibble(
  country = c("GIB","GLP","GUF","MTQ","MYT","REU"),
  lat_mean = c(36.14,16.25,3.93,14.63,-12.83,-21.12))
historical_map <- tibble::tibble(
  country = c("SCG","SUD"),
  lat_mean = c(44.0, 12.5))
tokelau_map <- tibble::tibble(
  country = "TKL",
  lat_mean = -9.2)
manual_lat <- bind_rows(manual_map, historical_map, tokelau_map)
net_mig_full <- net_mig_full %>%
  left_join(as.data.frame(world_iso_long[,c("country","lat_mean")])[,c("country","lat_mean")], by = "country") %>%  # join standard
  rows_update(manual_lat, by = "country")       # replace NA by the manual values
net_mig_full$country[is.na(net_mig_full$lat_mean)]

#merge with population
population = read.csv(file.path(DIR_IN, "WB_pop_country/population_data.csv"),    # skip the first line
  sep = ",",   # ou ";" si besoin
  header = TRUE)

population <- population %>%
  pivot_longer(
    cols = starts_with("X"),          # toutes les colonnes "X1960", "X1961", ...
    names_to = "year",                # new column for the year
    values_to = "value"               # new column for the values
  ) %>%
  mutate(
    year = as.integer(sub("X", "", year))  # enlever le "X" et convertir en entier
  ) %>%
  dplyr::select(Country.Code, year, value)  

population <- population %>%
  mutate(period = case_when(
    year >= 1990 & year < 1995 ~ "1990-1995",
    year >= 1995 & year < 2000 ~ "1995-2000",
    year >= 2000 & year < 2005 ~ "2000-2005",
    year >= 2005 & year < 2010 ~ "2005-2010",
    year >= 2010 & year < 2015 ~ "2010-2015",
    year >= 2015 & year < 2020 ~ "2015-2020",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(Country.Code, period) %>%
  summarise(mean_population = mean(value, na.rm = TRUE), .groups = "drop")

colnames(population) = c("country","period","mean_population")
population <- population %>%
  left_join(as.data.frame(world_iso_long[,c("country","lat_mean")])[,c("country","lat_mean")], by = "country") %>%  # join standard
  rows_update(manual_lat, by = "country", unmatched = "ignore")       # replace NA by the manual values

net_mig_full <- net_mig_full %>%
  mutate(period = cut(year0, breaks = seq(1990, 2020, by = 5),
                        right = FALSE,
                      labels = c("1990-1995","1995-2000", "2000-2005","2005-2010","2010-2015", "2015-2020")))
## ---- Une seule figure : stock net de migration par bande de latitude et periode ----
## Cinq bandes fixes ; pays affecte par la latitude de son centroide (lat_mean).
band_breaks <- c(-90, -23.5, 23.5, 45, 60, 90)
band_labels <- c("South of 23.5\u00b0S", "Tropics (23.5\u00b0S\u201323.5\u00b0N)",
                 "23.5\u201345\u00b0N", "45\u201360\u00b0N", "North of 60\u00b0N")

mig_bands <- net_mig_full %>%
  filter(!is.na(lat_mean), !is.na(period)) %>%
  mutate(band = cut(lat_mean, breaks = band_breaks, labels = band_labels,
                    include.lowest = TRUE, right = FALSE)) %>%
  group_by(period, band) %>%
  summarise(net_stock_m = sum(net_migration, na.rm = TRUE) / 1e6, .groups = "drop")

## Controle : les flux internationaux se compensent, la somme des bandes doit
## etre nulle a chaque periode. Un ecart signale des pays sans latitude ou des
## flux vers des entites hors carte ; il est imprime, pas masque.
check_sum <- mig_bands %>% group_by(period) %>%
  summarise(sum_bands = sum(net_stock_m), .groups = "drop")
print(check_sum)
dropped <- net_mig_full %>% filter(is.na(lat_mean)) %>% distinct(country) %>% pull(country)
if (length(dropped) > 0) message("Pays sans latitude, exclus : ", paste(dropped, collapse = ", "))

p_bands <- ggplot(mig_bands, aes(x = period, y = net_stock_m,
                                  group = band, colour = band, shape = band)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("#7f7f7f", "#c8302b", "#e0a030", "#1f5fbf", "#0b2f6b")) +
  scale_shape_manual(values = c(16, 15, 17, 18, 25)) +
  labs(x = NULL, y = "Net migration stock\n(millions per five-year period)",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_blank(), legend.key.height = unit(0.4, "cm"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(DIR_OUT, "myplot_migration.pdf"), plot = p_bands,
       width = 12, height = 8, units = "cm", device = cairo_pdf)
