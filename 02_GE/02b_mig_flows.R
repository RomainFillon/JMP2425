#####################################################################################
# 02b_mig_flows.R -- bilateral migration flow matrix at one-degree resolution.
#
# Combines the Abel and Cohen international bilateral flows with the Kummu et al.
# net internal stocks, allocates both to one-degree cells, caps the result by
# population, and writes the transition matrix the GE solver reads.
#
# Reads:  data_inputs/abel_cohen_migration/bilat_mig.csv (year0 = 2010, column da_pb_closed),
#         data_inputs/Kummu/raster_netMgr_2010_2014_5yrSum.tif,
#         intermediate_outputs/{identifier_xy_final,poids_par_cellule,
#         population_all,population0}.csv
# Writes: intermediate_outputs/mig_flows.h5, five-year transition probabilities
#           with rows summing to one and the diagonal holding the stayers
#         intermediate_outputs/mig_matrix_final_raw.h5, the raw five-year counts
#         plus the intermediate .rds objects of each step
# Needs:  02a must have run, since the cell grid and the population come from it.
#####################################################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(hdf5r)
  library(tidyr)
  library(terra)
  library(data.table)
})

#####################################################################################
# CONFIGURATION
#####################################################################################

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

PARAMS <- list(
  # Time horizon of flows (Abel-Cohen 2010-2015 and Kummu 2010-2014 are 5-year sums)
  period_years      = 5,

  # How to allocate international inflows to cells within a destination country:
  #   "pop"       : proportional to population (neutral, default)
  #   "stock_pos" : proportional to positive net stock
  # "stock_pos" mixes internal and international signals - prefer "pop" unless
  # you have a specific reason.
  inflow_allocation = "pop",

  # LOT 68 (30/09/2026) -- the internal block no longer crosses borders.
  # A country allocates its internal (Kummu residual) flows only among the cells
  # where it holds at least this share of the cell population. Measured on the
  # 29/09 inputs: with the whole-cell population and no restriction, 39.5 M of
  # the 114 M internal flows crossed a border once cells were attributed to
  # countries by population share (Syria +7.2 M "internal" inflow); with the
  # country's own population but no restriction, 16.5 M (Germany's whole internal
  # inflow landed in the Prague cell, 13 % German); with this restriction, 5.7 M.
  internal_min_part = 0.5,

  # Declared limit on the witness (STEP 15b): a country's net in mu_0 may differ
  # from its Abel-Cohen net by the error of the international block alone (caps,
  # domain, cell projection -- measured in the same run) plus this fraction of
  # the country's population, which bounds the residual leakage of the internal
  # block through mixed cells (measured maximum 2.7 %, Syria, on the 29/09 inputs).
  internal_leak_frac = 0.05,

  verbose           = TRUE
)

log_msg <- function(...) if (isTRUE(PARAMS$verbose)) cat(...)

# Period multiplier: flows are 5-year cumulated, population is a stock.
# We compare flows to pop * period_years ("person-years" available to migrate).
PERIOD_SCALE <- PARAMS$period_years

#####################################################################################
# STEP 1 - Load data
#####################################################################################
log_msg("=== STEP 1: Loading data ===\n")

# Identifier: (id, x, y) where id is the final 1:17048 index used for output
identifier <- read.csv(file.path(DIR_TMP, "identifier_xy_final.csv"))

# Population share of each country GID_0 within each (x, y) cell
## Lot 51 : population_weight_all.csv (parts = comptes de polygones GADM 3.6, RECAP 3.230)
## est remplace par poids_par_cellule.csv (00d_poids_population.py) : population de
## l'intersection cellule x pays, normalisee PAR CELLULE.
source(file.path(.root, "02_GE", "lib_poids.R"))
population_weight <- read.csv(file.path(DIR_TMP, "poids_par_cellule.csv"), na.strings = character(0))
population_weight <- population_weight[population_weight$GID_0 != "",
                                       c("GID_0", "x", "y", "part_de_la_cellule")]
temoin_part_de_la_cellule(population_weight)

# Abel-Cohen international bilateral flows 2010-2015
bilat_mig_1015 <- read.csv(file.path(DIR_IN, "/abel_cohen_migration/bilat_mig.csv"))
bilat_mig_1015 = bilat_mig_1015[bilat_mig_1015$year0==2010,]
bilat_mig_1015 <- bilat_mig_1015 %>%
  group_by(orig, dest) %>%
  summarise(da_pb_closed=sum(da_pb_closed,na.rm=T))

# Kummu net migration stocks at 1 degree (aggregated from 5' to 1 degree)
M_intra_raster <- rast(file.path(DIR_IN, "Kummu/raster_netMgr_2010_2014_5yrSum.tif"))
M_intra_raster <- aggregate(M_intra_raster, fact = 12, fun = sum)
M_intra_z      <- as.data.frame(M_intra_raster, xy = TRUE)
colnames(M_intra_z) <- c("x", "y", "value_z")
M_intra_z <- M_intra_z[!is.na(M_intra_z$value_z), ]

# Population per cell (independent of country breakdown)
population_all <- read.csv(file.path(DIR_TMP, "population_all.csv"))
population_all <- population_all[, c("x", "y", "value")]
colnames(population_all)[3] <- "pop_cell"

log_msg("  identifier:        ", nrow(identifier), " rows\n")
log_msg("  population_weight: ", nrow(population_weight), " rows\n")
log_msg("  bilat_mig_1015:    ", nrow(bilat_mig_1015), " rows\n")
log_msg("  M_intra_z:         ", nrow(M_intra_z), " rows\n")

#####################################################################################
# STEP 2 - Allocate net migration stocks to (zone x country) pairs
# value_zc = net stock in cell z attributable to country c
#          = net stock in z * population share of c in z
#####################################################################################
log_msg("\n=== STEP 2: Allocating net stocks to (zone, country) ===\n")

M_intra_zc <- left_join(population_weight, M_intra_z, by = c("x", "y"))
M_intra_zc$value_zc <- M_intra_zc$value_z * M_intra_zc$part_de_la_cellule

M_intra_zc <- M_intra_zc %>%
  group_by(GID_0, x, y) %>%
  summarise(value_zc = sum(value_zc, na.rm = TRUE),
            value_z  = sum(value_z,  na.rm = TRUE),
            .groups  = "drop") %>%
  mutate(zone1_id = as.integer(factor(paste(x, y, sep = "_"))))

log_msg("  rows M_intra_zc:        ", nrow(M_intra_zc), "\n")
log_msg("  unique zone1_id:        ", length(unique(M_intra_zc$zone1_id)), "\n")
log_msg("  unique GID_0:           ", length(unique(M_intra_zc$GID_0)), "\n")

#####################################################################################
# STEP 3 - Build cell-level population lookup (zone1_id -> pop)
# Built from population_all (cell-level) to avoid double-counting border cells
# that appear multiple times in M_intra_zc (one row per country sharing the cell).
#####################################################################################
log_msg("\n=== STEP 3: Building population lookup (zone-level) ===\n")

xy_to_zone <- distinct(M_intra_zc, x, y, zone1_id)

pop_by_cell <- left_join(population_all, xy_to_zone, by = c("x", "y"))
pop_by_cell <- pop_by_cell[!is.na(pop_by_cell$zone1_id), ]

pop_lookup <- pop_by_cell %>%
  group_by(zone1_id) %>%
  summarise(pop = sum(pop_cell, na.rm = TRUE), .groups = "drop")

log_msg("  pop_lookup rows:        ", nrow(pop_lookup), "\n")
log_msg("  total population:       ", round(sum(pop_lookup$pop) / 1e9, 3), " B\n")

#####################################################################################
# STEP 4 - Dedup (GID_0, zone1_id) and compute inflow/outflow probabilities
# Outflow: proportional to population of (country, zone) pair
# Inflow:  proportional to population (PARAMS$inflow_allocation = "pop"), default.
#####################################################################################
log_msg("\n=== STEP 4: Computing inflow/outflow probabilities ===\n")

# Country-weighted cell population, for outflow probas (double-counting is legitimate
# here: each country "claims" its share of a border cell's population).
pop_zc <- population_weight %>%
  left_join(population_all, by = c("x", "y")) %>%
  mutate(pop_zc = part_de_la_cellule * pop_cell) %>%
  left_join(xy_to_zone, by = c("x", "y")) %>%
  filter(!is.na(zone1_id)) %>%
  group_by(GID_0, zone1_id) %>%
  summarise(pop_zc = sum(pop_zc, na.rm = TRUE), .groups = "drop")

M_intra_zc <- M_intra_zc %>%
  group_by(GID_0, zone1_id) %>%
  summarise(
    value_zc = sum(value_zc, na.rm = TRUE),
    value_z  = sum(value_z,  na.rm = TRUE),
    x        = first(x),
    y        = first(y),
    .groups  = "drop"
  ) %>%
  left_join(pop_zc, by = c("GID_0", "zone1_id"))

M_intra_zc$pop_zc[is.na(M_intra_zc$pop_zc)] <- 0

log_msg("  unique (GID_0, zone1_id): ", nrow(distinct(M_intra_zc, GID_0, zone1_id)), "\n")

M_intra_zc <- M_intra_zc %>%
  group_by(GID_0) %>%
  mutate(
    total_pop_c     = sum(pop_zc, na.rm = TRUE),
    P_out           = ifelse(total_pop_c > 0, pop_zc / total_pop_c, 0),
    stock_pos       = pmax(value_zc, 0),
    total_stock_pos = sum(stock_pos, na.rm = TRUE),
    P_in = if (PARAMS$inflow_allocation == "stock_pos") {
      ifelse(total_stock_pos > 0,
             stock_pos / total_stock_pos,
             ifelse(total_pop_c > 0, pop_zc / total_pop_c, 0))
    } else {
      ifelse(total_pop_c > 0, pop_zc / total_pop_c, 0)
    }
  ) %>%
  select(-total_pop_c, -total_stock_pos) %>%
  ungroup()

# Check: probas sum to 1 per country
check_probs <- M_intra_zc %>%
  group_by(GID_0) %>%
  summarise(sum_out = sum(P_out, na.rm = TRUE),
            sum_in  = sum(P_in,  na.rm = TRUE),
            .groups = "drop")

n_bad_out <- sum(abs(check_probs$sum_out - 1) > 0.01 & check_probs$sum_out > 0)
n_bad_in  <- sum(abs(check_probs$sum_in  - 1) > 0.01 & check_probs$sum_in  > 0)
log_msg("  countries with sum(P_out) != 1: ", n_bad_out, "\n")
log_msg("  countries with sum(P_in)  != 1: ", n_bad_in, "\n")


#####################################################################################
# STEP 5 - Merge with identifier, handle cells not in Abel-Cohen country list
#####################################################################################
log_msg("\n=== STEP 5: Merging with identifier ===\n")

M_intra_zc <- left_join(M_intra_zc, identifier, by = c("x", "y"))

country_list <- unique(c(bilat_mig_1015$orig, bilat_mig_1015$dest))

df_all               <- M_intra_zc
M_intra_zc           <- M_intra_zc[M_intra_zc$GID_0 %in% country_list, ]
country_list_updated <- unique(M_intra_zc$GID_0)

# Cells with an identifier id but not matched to any Abel-Cohen country -> ZZZ_missing
missing_ids <- setdiff(identifier$id, M_intra_zc$id)
df_missing  <- df_all[df_all$id %in% missing_ids, ]
if (nrow(df_missing) > 0) {
  df_missing$GID_0     <- "ZZZ_missing"
  df_missing$P_out     <- 0
  df_missing$P_in      <- 0
  df_missing$pop_zc    <- 0
  df_missing$stock_pos <- 0
  M_intra_zc <- rbind(M_intra_zc, df_missing)
}

translation <- distinct(M_intra_zc[, c("id", "zone1_id")])

log_msg("  rows after merge:       ", nrow(M_intra_zc), "\n")
log_msg("  cells ZZZ_missing:      ", sum(M_intra_zc$GID_0 == "ZZZ_missing"), "\n")


#####################################################################################
# DIAGNOSTICS ALIGNEMENT id / zone1_id / GID_0
#####################################################################################
cat("\n=== DIAGNOSTICS ALIGNEMENT ===\n")

# 1. how many unique (id, zone1_id) pairs?
cat("N (id, zone1_id) uniques dans M_intra_zc:", 
    nrow(distinct(M_intra_zc[!is.na(M_intra_zc$id), ], id, zone1_id)), "\n")

# 2. ids present in identifier but absent from M_intra_zc
cat("IDs dans identifier mais pas dans M_intra_zc:", 
    length(setdiff(identifier$id, M_intra_zc$id[!is.na(M_intra_zc$id)])), "\n")

# 3. ids mapping to several zone1_id, that is, border cells
multi_zone <- M_intra_zc %>%
  filter(!is.na(id)) %>%
  group_by(id) %>%
  summarise(n_zone = n_distinct(zone1_id), .groups="drop") %>%
  filter(n_zone > 1)
cat("IDs avec plusieurs zone1_id:", nrow(multi_zone), "\n")
if (nrow(multi_zone) > 0) {
  cat("Exemples:\n")
  print(head(multi_zone, 5))
}

# 4. zone1_id mapping to several ids, the reverse alignment problem
multi_id <- M_intra_zc %>%
  filter(!is.na(id)) %>%
  group_by(zone1_id) %>%
  summarise(n_id = n_distinct(id), .groups="drop") %>%
  filter(n_id > 1)
cat("zone1_id avec plusieurs id:", nrow(multi_id), "\n")
if (nrow(multi_id) > 0) {
  cat("Exemples:\n")
  print(head(multi_id, 5))
}

# 5. translation table: check its coherence
cat("Lignes dans translation:", nrow(translation), "\n")
cat("IDs uniques dans translation:", n_distinct(translation$id), "\n")
cat("zone1_id uniques dans translation:", n_distinct(translation$zone1_id), "\n")
cat("IDs NA dans translation:", sum(is.na(translation$id)), "\n")

# 6. check that list_zone1_id_vec covers every active id
# (done after step 12)

#####################################################################################
# STEP 6 - Build sparse probability matrices (zone x country)
# Rows: zone1_id sorted, Columns: country GID_0 sorted alphabetically
#####################################################################################
log_msg("\n=== STEP 6: Building sparse probability matrices ===\n")

to_sparse_zone_country <- function(df, value_col) {
  wide <- df %>%
    select(zone1_id, GID_0, all_of(value_col)) %>%
    group_by(zone1_id, GID_0) %>%
    summarise(v = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = GID_0, values_from = v, values_fill = 0) %>%
    arrange(zone1_id)

  zone1_order <- wide$zone1_id
  wide <- wide[, colnames(wide) != "zone1_id", drop = FALSE]
  wide <- wide[, sort(colnames(wide)), drop = FALSE]
  mat  <- Matrix(data.matrix(wide), sparse = TRUE)
  list(matrix = mat, zone1_id = zone1_order)
}

out_obj <- to_sparse_zone_country(M_intra_zc, "P_out")
in_obj  <- to_sparse_zone_country(M_intra_zc, "P_in")

Matrix_probaout_zc <- out_obj$matrix
Matrix_probain_zc  <- in_obj$matrix
list_zone1_id_vec  <- out_obj$zone1_id

stopifnot(identical(out_obj$zone1_id, in_obj$zone1_id))
stopifnot(identical(colnames(Matrix_probaout_zc), colnames(Matrix_probain_zc)))

log_msg("  Matrix_probaout_zc dim: ", paste(dim(Matrix_probaout_zc), collapse = " x "), "\n")
log_msg("  Matrix_probain_zc  dim: ", paste(dim(Matrix_probain_zc),  collapse = " x "), "\n")

#####################################################################################
# STEP 7 - Build country x country matrix from Abel-Cohen
# Rows = origin (sorted alphabetically), Cols = destination (sorted alphabetically)
#####################################################################################
log_msg("\n=== STEP 7: Building country-country matrix ===\n")

bilat_mig_1015 <- bilat_mig_1015 %>%
  filter(orig %in% country_list_updated, dest %in% country_list_updated)

# ZZZ_missing as placeholder country so matrix dimensions align with P_in / P_out.
# Flows to/from ZZZ_missing are zero.
bilat_mig_1015 <- rbind(
  bilat_mig_1015,
  data.frame(orig = c("ZZZ_missing", country_list_updated[1]),
             dest = c(country_list_updated[1], "ZZZ_missing"),
             da_pb_closed = 0)
)

Matrix_cc <- bilat_mig_1015 %>%
  pivot_wider(names_from = dest, values_from = da_pb_closed,
              values_fill = 0, values_fn = sum) %>%
  arrange(orig)

orig_names <- Matrix_cc$orig
Matrix_cc  <- as.data.frame(Matrix_cc[, colnames(Matrix_cc) != "orig"])
rownames(Matrix_cc) <- orig_names

Matrix_cc <- Matrix_cc[order(rownames(Matrix_cc)), order(colnames(Matrix_cc)), drop = FALSE]
Matrix_cc <- Matrix(data.matrix(Matrix_cc), sparse = TRUE)

# Safety: names must align before matrix products
stopifnot(identical(colnames(Matrix_probaout_zc), rownames(Matrix_cc)))
stopifnot(identical(colnames(Matrix_cc),          colnames(Matrix_probain_zc)))

log_msg("  Matrix_cc dim:          ", paste(dim(Matrix_cc), collapse = " x "), "\n")
log_msg("  total Abel-Cohen flows: ", round(sum(Matrix_cc) / 1e6, 2), " M\n")

#####################################################################################
# STEP 8 - Build international zone-zone flow matrix
# migration_matrix[z1, z2] = sum_{c1, c2} P_out[z1, c1] * M_cc[c1, c2] * P_in[z2, c2]
#####################################################################################
log_msg("\n=== STEP 8: Building international zone-zone matrix ===\n")

step1 <- Matrix_probaout_zc %*% Matrix_cc
log_msg("  after step1 (~ Abel-Cohen): ", round(sum(step1) / 1e6, 2), " M\n")

migration_matrix <- step1 %*% t(Matrix_probain_zc)
log_msg("  migration_matrix total:     ", round(sum(migration_matrix) / 1e6, 2), " M\n")

# Diagonal first (before capping): zeroes out "country X -> country X" ghost flows
# that leaked through probability allocation so they do not inflate row sums.
diag(migration_matrix) <- 0

#####################################################################################
# STEP 9 - Cap international flows by zone population
# Flows are 5-year cumulated, population is a stock -> cap against pop * period_years.
#####################################################################################
log_msg("\n=== STEP 9: Capping international flows by population ===\n")

pop_migration <- pop_lookup$pop[match(list_zone1_id_vec, pop_lookup$zone1_id)]
pop_migration[is.na(pop_migration)] <- 0

stopifnot(length(pop_migration) == nrow(migration_matrix))

pop_available <- pop_migration * PERIOD_SCALE
row_sums_mig  <- rowSums(migration_matrix)
scale_mig     <- ifelse(row_sums_mig > pop_available & row_sums_mig > 0,
                        pop_available / row_sums_mig,
                        1)
scale_mig[is.na(scale_mig)]   <- 0
scale_mig[pop_migration == 0] <- 0

migration_matrix <- Diagonal(x = scale_mig) %*% migration_matrix

log_msg("  zones capped:           ", sum(scale_mig < 1 & scale_mig > 0), "\n")
log_msg("  zones zeroed (pop=0):   ", sum(pop_migration == 0), "\n")
log_msg("  total after cap:        ", round(sum(migration_matrix) / 1e6, 2), " M\n")

# Persist intermediate artefacts
saveRDS(migration_matrix,
        file.path(DIR_TMP, "migration_matrix.rds"), compress = "xz")
saveRDS(M_intra_zc,
        file.path(DIR_TMP, "M_intra_zc.rds"), compress = "xz")
saveRDS(list_zone1_id_vec,
        file.path(DIR_TMP, "list_zone1_id_vec.rds"), compress = "xz")
saveRDS(pop_lookup,
        file.path(DIR_TMP, "pop_lookup.rds"), compress = "xz")

#####################################################################################
# STEP 10 - Compute residual stocks after removing international flows
# stock[z, c] = value_zc[z, c] - net_international_flow[z, c]
#####################################################################################
log_msg("\n=== STEP 10: Computing residual stocks ===\n")

inflows  <- colSums(migration_matrix)
outflows <- rowSums(migration_matrix)
net_flow <- data.frame(
  zone1_id = list_zone1_id_vec,
  net_flow = as.numeric(inflows - outflows)
)

M_intra_zc <- left_join(M_intra_zc, net_flow, by = "zone1_id")
M_intra_zc$net_flow[is.na(M_intra_zc$net_flow)] <- 0
## LOT 68 : the international net flow of a zone is SHARED among the countries of
## the zone in proportion to pop_zc. It was subtracted in full from each of them,
## so that a border cell carried its whole net flow once per country touching it.
M_intra_zc <- M_intra_zc %>%
  group_by(zone1_id) %>%
  mutate(part_z = { s <- sum(pop_zc, na.rm = TRUE); if (s > 0) pop_zc / s else 0 }) %>%
  ungroup()
M_intra_zc$stock <- M_intra_zc$value_zc - M_intra_zc$net_flow * M_intra_zc$part_z
M_intra_zc <- M_intra_zc[, c("zone1_id", "GID_0", "stock", "pop_zc", "part_z")]

log_msg("  rows with stock != 0:   ",
        sum(M_intra_zc$stock != 0, na.rm = TRUE), "\n")

#####################################################################################
# STEP 11 - Allocate internal (within-country) flows
# For each country: distribute flows between zones proportional to
#   outflow: population share
#   inflow:  positive residual stock share (fallback to population if total_pos = 0)
# Total internal flow mass = max(total_pos, |total_neg|), capped at total_pop.
# Uses outer() instead of expand.grid() for speed; diagonal set to zero.
#####################################################################################
log_msg("\n=== STEP 11: Allocating internal flows ===\n")

allocate_internal_flows <- function(df_country, pop_lookup_local) {
  ## LOT 68 : the country's OWN population in the cell (pop_zc), never the whole
  ## cell (pop_lookup), and only the cells where the country holds at least
  ## PARAMS$internal_min_part of the population. pop_lookup_local is kept in the
  ## signature and no longer read.
  df_country <- df_country[df_country$part_z >= PARAMS$internal_min_part, ]
  if (nrow(df_country) == 0) return(NULL)

  df_country$pop <- df_country$pop_zc
  df_country$pop[is.na(df_country$pop)] <- 0

  total_pos <- sum(df_country$stock[df_country$stock > 0], na.rm = TRUE)
  total_neg <- sum(df_country$stock[df_country$stock < 0], na.rm = TRUE)
  total_pop <- sum(df_country$pop, na.rm = TRUE)

  # Cap total outflow volume by total population (flows <= population)
  maximum_stock <- min(max(total_pos, abs(total_neg)), total_pop)
  if (maximum_stock == 0) return(NULL)

  total_pop_c <- sum(df_country$pop, na.rm = TRUE)
  P_outflow <- if (total_pop_c > 0) df_country$pop / total_pop_c else rep(0, nrow(df_country))

  stock_pos       <- pmax(df_country$stock, 0)
  total_stock_pos <- sum(stock_pos, na.rm = TRUE)
  P_inflow <- if (total_stock_pos > 0) stock_pos / total_stock_pos else P_outflow

  # Outer product builds the flow matrix for this country directly
  flow_mat <- maximum_stock * (P_outflow %o% P_inflow)
  diag(flow_mat) <- 0

  nz <- which(flow_mat != 0, arr.ind = TRUE)
  if (nrow(nz) == 0) return(NULL)

  data.frame(
    orig          = df_country$zone1_id[nz[, 1]],
    dest          = df_country$zone1_id[nz[, 2]],
    flow_internal = flow_mat[nz]
  )
}

split_data    <- split(M_intra_zc, M_intra_zc$GID_0)
results       <- lapply(split_data, allocate_internal_flows, pop_lookup_local = pop_lookup)
dataset_final <- bind_rows(results)

log_msg("  internal flow rows:     ", nrow(dataset_final), "\n")
log_msg("  total internal flows:   ", round(sum(dataset_final$flow_internal) / 1e6, 2), " M\n")

#####################################################################################
# STEP 12 - Convert internal flows to sparse matrix indexed by zone1_id
#####################################################################################
log_msg("\n=== STEP 12: Building sparse internal-flow matrix ===\n")

dt <- as.data.table(dataset_final)
dt <- dt[flow_internal != 0 & !is.na(orig) & !is.na(dest) & !is.nan(flow_internal)]
dt_sum <- dt[, .(flow_internal = sum(flow_internal, na.rm = TRUE)), by = .(orig, dest)]
rm(dt)

# list_zone1_id: translation table zone1_id <-> row-position <-> id
list_zone1_id <- data.frame(
  zone1_id = list_zone1_id_vec,
  row      = seq_along(list_zone1_id_vec)
)
list_zone1_id <- left_join(list_zone1_id, translation, by = "zone1_id")
list_zone1_id <- distinct(list_zone1_id, zone1_id, row, id)

saveRDS(list_zone1_id, file.path(DIR_TMP, "list_zone1_id.rds"))

zones_valides <- sort(unique(list_zone1_id$zone1_id))
orig_idx      <- match(dt_sum$orig, zones_valides)
dest_idx      <- match(dt_sum$dest, zones_valides)

keep     <- !is.na(orig_idx) & !is.na(dest_idx)
orig_idx <- orig_idx[keep]
dest_idx <- dest_idx[keep]
x_vals   <- dt_sum$flow_internal[keep]

dataset_final_sparse <- sparseMatrix(
  i        = orig_idx,
  j        = dest_idx,
  x        = x_vals,
  dims     = c(length(zones_valides), length(zones_valides)),
  dimnames = list(as.character(zones_valides), as.character(zones_valides))
)

log_msg("  dataset_final_sparse dim: ", paste(dim(dataset_final_sparse), collapse = " x "), "\n")
log_msg("  total internal flows:     ", round(sum(dataset_final_sparse) / 1e6, 2), " M\n")


# 7. Verifier list_zone1_id
cat("\nlist_zone1_id:\n")
cat("  N lignes:", nrow(list_zone1_id), "\n")
cat("  IDs uniques:", n_distinct(list_zone1_id$id[!is.na(list_zone1_id$id)]), "\n")
cat("  IDs NA:", sum(is.na(list_zone1_id$id)), "\n")
cat("  zone1_id uniques:", n_distinct(list_zone1_id$zone1_id), "\n")

# ids present in identifier but absent from list_zone1_id
missing_from_list <- setdiff(
  identifier$id, 
  list_zone1_id$id[!is.na(list_zone1_id$id)]
)
cat("  IDs dans identifier mais pas dans list_zone1_id:", length(missing_from_list), "\n")

# 8. check that the indices of mig_flows_final match the ids
# After step 14

#####################################################################################
# STEP 13 - Subset both matrices to valid zones, align and sum
#####################################################################################
log_msg("\n=== STEP 13: Aligning and summing international + internal ===\n")

zones_to_keep <- sort(unique(list_zone1_id$row[!is.na(list_zone1_id$id)]))

migration_matrix_kept <- migration_matrix[zones_to_keep, zones_to_keep]

zone1_ids_to_keep <- list_zone1_id$zone1_id[match(zones_to_keep, list_zone1_id$row)]
zone1_ids_to_keep <- zone1_ids_to_keep[!is.na(zone1_ids_to_keep)]
zone1_ids_to_keep <- sort(zone1_ids_to_keep)

dataset_final_kept <- dataset_final_sparse[
  as.character(zone1_ids_to_keep),
  as.character(zone1_ids_to_keep)
]

# Re-label migration_matrix rows/cols with zone1_id so both matrices share the same index
dimnames(migration_matrix_kept) <- list(
  as.character(zone1_ids_to_keep),
  as.character(zone1_ids_to_keep)
)

final_matrix <- dataset_final_kept + migration_matrix_kept

log_msg("  final_matrix dim:       ", paste(dim(final_matrix), collapse = " x "), "\n")
log_msg("  total final flows:      ", round(sum(final_matrix) / 1e6, 2), " M\n")
log_msg("    international:        ", round(sum(migration_matrix_kept) / 1e6, 2), " M\n")
log_msg("    internal:             ", round(sum(dataset_final_kept)    / 1e6, 2), " M\n")

#####################################################################################
# STEP 14 - Re-index from zone1_id to identifier id (1:17048)
#####################################################################################
log_msg("\n=== STEP 14: Re-indexing to identifier id ===\n")

z1_to_id <- setNames(list_zone1_id$id, as.character(list_zone1_id$zone1_id))

nz          <- summary(final_matrix)
z1_rownames <- rownames(final_matrix)

new_i <- as.integer(z1_to_id[z1_rownames[nz$i]])
new_j <- as.integer(z1_to_id[z1_rownames[nz$j]])
keep  <- !is.na(new_i) & !is.na(new_j)

log_msg("  non-zero entries kept:  ", sum(keep), " / ", nrow(nz), "\n")

max_id <- max(identifier$id)

mig_flows_final <- sparseMatrix(
  i        = new_i[keep],
  j        = new_j[keep],
  x        = nz$x[keep],
  dims     = c(max_id, max_id),
  dimnames = list(as.character(seq_len(max_id)), as.character(seq_len(max_id)))
)

#####################################################################################
# STEP 15 - Final cap by identifier-level population
# Zones with pop=0 cannot emit or receive flows; outflows capped at pop * period.
#####################################################################################
log_msg("\n=== STEP 15: Final cap by identifier-level population ===\n")

pop_id <- read.csv(file.path(DIR_TMP, "population0.csv"))
pop_id$value[is.na(pop_id$value)] <- 0

pop_vec <- numeric(max_id)
pop_vec[pop_id$id] <- pop_id$value

zones_no_pop <- which(pop_vec == 0)
mig_flows_final[zones_no_pop, ] <- 0
mig_flows_final[, zones_no_pop] <- 0
log_msg("  zones pop=0 zeroed:     ", length(zones_no_pop), "\n")

pop_available_final <- pop_vec * PERIOD_SCALE
row_sums_final      <- rowSums(mig_flows_final)
scale_final <- ifelse(
  row_sums_final > pop_available_final & row_sums_final > 0,
  pop_available_final / row_sums_final,
  1
)
scale_final[pop_vec == 0] <- 0

mig_flows_final <- Diagonal(x = scale_final) %*% mig_flows_final

log_msg("  zones capped:           ", sum(scale_final < 1 & scale_final > 0), "\n")
log_msg("  total after final cap:  ", round(sum(mig_flows_final) / 1e6, 2), " M\n")

# Sanity diagnostic
row_sums_final <- rowSums(mig_flows_final)
ratio_final    <- row_sums_final / pmax(pop_available_final, 1)
log_msg("  zones flux > pop*T:     ", sum(ratio_final > 1,    na.rm = TRUE), "\n")
log_msg("  zones flux > 5*pop*T:   ", sum(ratio_final > 5,    na.rm = TRUE), "\n")
log_msg("  max ratio:              ", round(max(ratio_final, na.rm = TRUE), 4), "\n")

#####################################################################################
# STEP 15b - WITNESS (lot 68): the net of each country in mu_0 against Abel-Cohen
#
# Countries are read off cells as ge_calibration_diagnostic.py reads them: each
# cell's in- and outflows are attributed to its countries by part_de_la_cellule.
#   src_c    net international migration of c in the source, as read in STEP 1
#   intl_c   net of c in the INTERNATIONAL block of mu_0 alone, after STEP 13
#            (domain) and STEP 15 (zero-population cells, cap) -- same scale
#   total_c  net of c in mu_0 (international + internal)
# The international block reproduces the source exactly at the (zone, country)
# level; E_c = intl_c - src_c is therefore the error that the caps, the domain
# and the projection of countries onto cells cause, measured in this run.
# The internal block moves no one across a border by construction, but a mixed
# cell is attributed to several countries; its residual leakage is bounded by the
# declared limit PARAMS$internal_leak_frac x pop_c.
# REFUSES if, for any country, |total_c - src_c| > |E_c| + internal_leak_frac x pop_c.
# KNOWN-BAD, run every time: the same test on a copy where Syria's net is the
# opposite of its source net (the lot 67 defect: -5.24 M in the source, +2.58 M
# in the old mu_0) must refuse.
#####################################################################################
log_msg("\n=== STEP 15b: Witness -- net by country, mu_0 against the source ===\n")

# international block alone, on the id grid, with the same zeroing and the same scale
z_ids   <- as.integer(z1_to_id[rownames(migration_matrix_kept)])
ok_z    <- !is.na(z_ids)
Mk      <- migration_matrix_kept[ok_z, ok_z]
z_ids   <- z_ids[ok_z]
colmask <- as.numeric(pop_vec[z_ids] > 0)
Mk      <- Mk %*% Diagonal(x = colmask)
s_z     <- scale_final[z_ids]
out_int <- numeric(max_id); in_int <- numeric(max_id)
out_int[z_ids] <- s_z * as.numeric(rowSums(Mk))
in_int[z_ids]  <- as.numeric(crossprod(Mk, s_z))
rm(Mk)

out_tot <- as.numeric(rowSums(mig_flows_final)) - as.numeric(diag(mig_flows_final))
in_tot  <- as.numeric(colSums(mig_flows_final)) - as.numeric(diag(mig_flows_final))

# cell -> country, by population share
pw_id <- merge(population_weight, identifier[, c("x", "y", "id")], by = c("x", "y"))
w_cty <- sort(unique(c(pw_id$GID_0, rownames(Matrix_cc), colnames(Matrix_cc))))
w_cty <- setdiff(w_cty, "ZZZ_missing")
pw_id <- pw_id[pw_id$GID_0 %in% w_cty, ]
P_idc <- sparseMatrix(i = pw_id$id, j = match(pw_id$GID_0, w_cty), x = pw_id$part_de_la_cellule,
                      dims = c(max_id, length(w_cty)), dimnames = list(NULL, w_cty))
by_cty <- function(v) setNames(as.numeric(crossprod(P_idc, v)), w_cty)

cc_names <- intersect(rownames(Matrix_cc), w_cty)
src <- setNames(numeric(length(w_cty)), w_cty)
Mcc_off <- Matrix_cc; diag(Mcc_off) <- 0
src[cc_names] <- as.numeric(colSums(Mcc_off)[cc_names] - rowSums(Mcc_off)[cc_names])

intl  <- by_cty(in_int - out_int)
total <- by_cty(in_tot - out_tot)
pop_c <- by_cty(pop_vec)
E_c   <- intl - src
tol   <- abs(E_c) + PARAMS$internal_leak_frac * pop_c

witness_fails <- function(tot) names(which(abs(tot - src) > tol))

wit <- data.frame(GID_0 = w_cty, pop = pop_c, src = src, intl = intl,
                  internal = total - intl, total = total, E = E_c, tol = tol,
                  ok = abs(total - src) <= tol)
write.csv(wit, file.path(DIR_TMP, "mig_witness_by_country.csv"), row.names = FALSE)

o <- order(total)
fmt <- function(k) sprintf("    %s  mu0 %+6.2f M | intl %+6.2f | interne %+6.2f | source %+6.2f | pop %8.2f M",
                           w_cty[k], total[k]/1e6, intl[k]/1e6, (total[k]-intl[k])/1e6, src[k]/1e6, pop_c[k]/1e6)
log_msg("  5 emetteurs nets de mu_0 :\n"); for (k in head(o, 5)) log_msg(fmt(k), "\n")
log_msg("  5 recepteurs nets de mu_0 :\n"); for (k in rev(tail(o, 5))) log_msg(fmt(k), "\n")
log_msg(sprintf("  |E_c| max %.2f M ; |interne_c| max %.2f M ; |interne_c|/pop_c max %.2f %% (pays > 0,1 M hab)\n",
                max(abs(E_c))/1e6, max(abs(total - intl))/1e6,
                100 * max(abs((total - intl)[pop_c > 1e5] / pop_c[pop_c > 1e5]))))
log_msg(sprintf("  masse transfrontaliere : sum |interne_c| / 2 = %.2f M\n", sum(abs(total - intl))/2e6))

stopifnot("SYR" %in% w_cty, src[["SYR"]] < 0)
bad <- total; bad[["SYR"]] <- -src[["SYR"]]
kb <- witness_fails(bad)
if (!("SYR" %in% kb))
  stop("[STOP] STEP 15b : le temoin ne refuse pas le connu-faux (Syrie inversee). La garde est morte.")
log_msg("  [known-bad] Syrie inversee (", round(bad[["SYR"]]/1e6, 2), " M) : refusee, comme attendu\n")

fails <- witness_fails(total)
if (length(fails)) {
  print(wit[wit$GID_0 %in% fails, ])
  stop("[STOP] STEP 15b : solde net de mu_0 hors tolerance pour ", paste(fails, collapse = ", "),
       ". mig_flows.h5 n'est PAS ecrit. Detail : intermediate_outputs/mig_witness_by_country.csv")
}
log_msg("  [ok] temoin : ", length(w_cty), " pays dans la tolerance (|E_c| + ",
        PARAMS$internal_leak_frac, " x pop_c)\n")

#####################################################################################
# STEP 16 - Write raw flows and probability matrix to HDF5
# mig_matrix_final_raw.h5 : raw 5-year flows (counts)
# mig_flows.h5            : 5-year transition probabilities (rows sum to 1, diag = stay)
##########################################################################################################################################################################
# STEP 16 - Write raw flows and probability matrix to HDF5
#####################################################################################
log_msg("\n=== STEP 16: Writing HDF5 outputs ===\n")

write_h5_matrix <- function(mat, filepath, dataset_name = "mig_flows",
                            rownames_vec = NULL, chunk_rows = 100,
                            block_rows = 200, transform = NULL) {
  n <- nrow(mat)
  if (file.exists(filepath)) file.remove(filepath)
  h5 <- H5File$new(filepath, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  h5$create_dataset(name = dataset_name, dims = c(n, n),
                    dtype = h5types$H5T_NATIVE_FLOAT,
                    chunk_dims = c(chunk_rows, n))
  if (!is.null(rownames_vec)) {
    h5$create_dataset(name = "rownames", robj = rownames_vec)
  }
  for (row_start in seq(1, n, by = block_rows)) {
    row_end <- min(row_start + block_rows - 1, n)
    bloc    <- as.matrix(mat[row_start:row_end, ])
    if (!is.null(transform)) bloc <- transform(bloc, row_start, row_end)
    h5[[dataset_name]][row_start:row_end, ] <- bloc
    if (row_end %% 2000 == 0 || row_end == n)
      log_msg("    rows ", row_start, " -> ", row_end, " / ", n, "\n")
  }
}

# ---- Raw flows ----
log_msg("  writing mig_matrix_final_raw.h5\n")
write_h5_matrix(
  mig_flows_final,
  filepath     = file.path(DIR_TMP, "mig_matrix_final_raw.h5"),
  rownames_vec = dimnames(mig_flows_final)[[1]]
)

# ---- Threshold ----
mig_flows_final@x[mig_flows_final@x < 1] <- 0
mig_flows_final <- drop0(mig_flows_final)
cat("Sparsity after threshold:", round(nnzero(mig_flows_final) / prod(dim(mig_flows_final)), 4), "\n")
cat("nnz after threshold:", nnzero(mig_flows_final), "\n")

# ---- Diagnostic before writing ----
#pop_denom <- pmax(pop_available_final, 1)
pop_denom <- pmax(pop_vec, 1)
cat("pop_denom min:", min(pop_denom), "max:", max(pop_denom), "\n")
cat("pop_denom zeros:", sum(pop_denom == 1 & pop_available_final == 0), "\n")

# Test on the first block
bloc_test  <- as.matrix(mig_flows_final[1:5, ])
rs_test    <- rowSums(bloc_test)
stay_test  <- pmax(pop_denom[1:5] - rs_test, 0)
cat("Test bloc 1-5:\n")
cat("  row_sums:", round(rs_test, 2), "\n")
cat("  pop_denom:", round(pop_denom[1:5], 2), "\n")
cat("  stay:", round(stay_test, 2), "\n")
cat("  stay/pop:", round(stay_test / pop_denom[1:5], 4), "\n")

# ---- probability_transform ----

probability_transform <- function(bloc, row_start, row_end) {
  pop_bloc      <- pop_denom[row_start:row_end]
  row_sums_bloc <- rowSums(bloc, na.rm = TRUE)

  over <- row_sums_bloc > pop_bloc
  if (any(over)) {
    scale <- ifelse(over, pop_bloc / row_sums_bloc, 1)
    bloc  <- bloc * scale
    row_sums_bloc <- rowSums(bloc, na.rm = TRUE)
  }

  stay      <- pmax(pop_bloc - row_sums_bloc, 0)
  bloc_prob <- bloc / pop_bloc

  # global column index = row_start + i - 1
  n_rows <- nrow(bloc_prob)
  for (i in seq_len(n_rows)) {
    col_global <- row_start + i - 1
    bloc_prob[i, col_global] <- bloc_prob[i, col_global] + (stay[i] / pop_bloc[i])
  }

  rs <- rowSums(bloc_prob)
  bloc_prob / pmax(rs, 1e-12)
}


# ? AJOUTER ICI
bloc_test  <- as.matrix(mig_flows_final[1:200, ])
result     <- probability_transform(bloc_test, 1, 200)
cat("Test direct probability_transform bloc 1-200:\n")
cat("  diag min:", min(diag(result)), "\n")
cat("  diag median:", median(diag(result)), "\n")
cat("  diag zeros:", sum(diag(result) == 0), "\n")
cat("  row sums min:", min(rowSums(result)), "max:", max(rowSums(result)), "\n")


# ---- Probabilites ----
log_msg("  writing mig_flows.h5 (probabilities)\n")
write_h5_matrix(
  mig_flows_final,
  filepath     = file.path(DIR_TMP, "mig_flows.h5"),
  rownames_vec = dimnames(mig_flows_final)[[1]],
  transform    = probability_transform
)

log_msg("=== DONE ===\n")

# ---- Verification diagonale ----
cat("\n--- Verification diagonale mig_flows.h5 ---\n")
h5_check <- H5File$new(file.path(DIR_TMP, "mig_flows.h5"), mode = "r")
n_check  <- h5_check[["mig_flows"]]$dims[1]
cat("Dimension:", n_check, "x", n_check, "\n")

diag_mig <- numeric(n_check)
chunk_diag <- 200
for (row_start in seq(1, n_check, by = chunk_diag)) {
  row_end  <- min(row_start + chunk_diag - 1, n_check)
  bloc     <- h5_check[["mig_flows"]][row_start:row_end, row_start:row_end]
  idx      <- seq_len(row_end - row_start + 1)
  diag_mig[row_start:row_end] <- bloc[cbind(idx, idx)]
}
h5_check$close_all()

cat("N diagonales = 0:", sum(diag_mig == 0), "\n")
cat("N diagonales < 0.5:", sum(diag_mig < 0.5), "\n")
cat("Diagonale: min=", round(min(diag_mig), 4),
    "median=", round(median(diag_mig), 4),
    "max=", round(max(diag_mig), 4), "\n")

pop_id_check <- read.csv(file.path(DIR_TMP, "population0.csv"))
pop_vec_check <- numeric(n_check)
pop_vec_check[pop_id_check$id] <- pop_id_check$value

pop_zero_diag <- pop_vec_check[diag_mig == 0]
cat("Population des locations diag=0: min=", min(pop_zero_diag),
    "median=", median(pop_zero_diag),
    "max=", max(pop_zero_diag), "\n")
cat("N locations diag=0 avec pop < 20:", sum(pop_zero_diag < 20), "\n")
cat("N locations diag=0 avec pop > 1000:", sum(pop_zero_diag > 1000), "\n")
