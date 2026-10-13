#!/usr/bin/env Rscript
# =============================================================================
# 00_temperature_projections_JMP.R
#
# The single producer of every FUTURE temperature used downstream, for the GE and
# for the diagnostics. It is the only script in the package that reads netCDF.
#
# Reads:  ISIMIP3b bias-adjusted (W5E5) netCDF, SSP2-4.5, daily tas, 0.5 degree,
#         four models over 2015-2100. The models are discovered by listing the
#         directory, then checked against CMIP6_MODELS (00_parameters.R): a
#         model or a decade that was never downloaded stops the run.
#         intermediate_outputs/identifier_xy_final.csv, for the cell grid.
# Writes: tas_ssp2_daily.csv, tas_ssp2_annual.csv, CMIP6_tas.csv
#
# THE OUTPUTS ARE KEYED BY THE CELL ID OF identifier_xy_final.csv, and
# 02a_calibration_at_1deg.r RENUMBERS THOSE IDS when it runs. This script must
# therefore be re-run after 02a. Otherwise every cell receives the temperature of
# an unrelated cell, the id vectors still match, and nothing reports it.
#
# PROCESSING
#   netCDF -> aggregate 0.5 to 1 degree (mean) -> K to degC -> integer bins
#   -> mask onto the cell grid -> count days per (id, temperature bin, 5-year
#      period) PER MODEL -> average ACROSS MODELS on the share of days.
#   2015 is included: the bias-adjusted W5E5 input is already on the observed
#   scale, so no splice with ERA5 is needed.
#
# OUTPUT SHAPES
#   tas_ssp2_daily.csv   daily distribution, 5-year FORWARD blocks, multi-model
#                        wide: id, value, id_value, <2015, 2020, ..., 2095>
#                        counts are about 365 days per cell-year; the labels are
#                        block starts, and the last block is 2095-2099
#   tas_ssp2_annual.csv  annual mean temperature
#                        long: id, year, temp_mean = sum(value * share)
#   CMIP6_tas.csv        distribution at the cell (x, y), for the diagnostics
#                        long: x, y, year, value, count (= share)
#
# Runs on the server: the netCDF inputs are not shipped with the package.
# =============================================================================
suppressMessages({ library(terra); library(data.table) })


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
#source(file.path(Sys.getenv("JMP_ROOT"), "01_estimation", "lib_drf.R"))


# ------- PATHS (editable, or passed as input_dir identifier_csv out_dir) ----
#a <- commandArgs(TRUE)
input_dir = file.path(DIR_IN, "CMIP6")
identifier_csv = file.path(DIR_TMP, "identifier_xy_final.csv")

bin_choice  <- 0L          # bin rounding (0 = 1 degC)
PERIOD_LEN  <- 5L          # model period = 5 years
YEAR_MIN    <- 2015L       # first period label = the base block {2015-2019}
YEAR_MAX    <- 2095L       # DERNIERE etiquette = dernier bloc PLEIN {2095-2099}.
                           # 2100 ouvrirait un bloc degenere d'1 an ({2100}) : EXCLU
                           # The grid is the GE grid, not the raster grid.

cat("=== temperature_projections_JMP.R ===\n")
cat("input_dir      :", input_dir, "\n")
cat("identifier_csv :", identifier_csv, "\n")
stopifnot(dir.exists(input_dir), file.exists(identifier_csv), dir.exists(DIR_OUT))

# ------- 1. dynamic discovery of the files and models -----------------------
nc_files <- list.files(input_dir, pattern="_w5e5_ssp245_tas_.*\\.nc$", full.names=TRUE)
stopifnot(length(nc_files) > 0)
model_of <- function(f) sub("_r1i1p1f1.*", "", basename(f))
## The set must be CMIP6_MODELS (00_parameters.R), with as many files per model:
## discovery alone cannot see a model, or a decade, that was never downloaded.
models   <- cmip6_garde_modeles(nc_files)
cat(sprintf("Fichiers .nc : %d\n", length(nc_files)))
cat(sprintf("Modeles decouverts (%d) : %s\n\n", length(models), paste(models, collapse=", ")))

# ------- 2. one-degree grid and mask, from the identifier -------------------
identifier <- fread(identifier_csv)
setnames(identifier, tolower(names(identifier)))
stopifnot(all(c("x","y","id") %in% names(identifier)))
r <- rast(ncols=360, nrows=180, xmin=-180, xmax=180, ymin=-90, ymax=90)
values(r) <- NA
r[cellFromXY(r, as.matrix(identifier[, .(x, y)]))] <- 1
cat(sprintf("Masque : %d cellules id\n\n", nrow(identifier)))

# Five-year FORWARD convention, aligned on the GE time steps, blocks starting
# label L = debut du bloc -> fenetre [L .. L+4]. Etiquettes 2015, 2020, ..., 2095.
# base = {2015-2019} (label 2015) ; dernier point PLEIN = {2095-2099} (label 2095).
# 2100 would be a degenerate one-year block and is EXCLUDED by the filter.
# This is the grid the GE consumes, symmetric between baseline and future.
period_of <- function(y) (y %/% PERIOD_LEN) * PERIOD_LEN

# ------- 3. loop over files -> counts (model, id, value, period) ------------
#   Each file is reduced immediately and collected, which keeps memory bounded.
acc <- vector("list", length(nc_files))
for (k in seq_along(nc_files)) {
  f  <- nc_files[k]; mdl <- model_of(f)
  ra <- rast(f)
  yr <- as.integer(format(time(ra), "%Y"))
  # agrege 0.5->1 deg si besoin, K->degC, masque
  if (res(ra)[1] < 0.99) ra <- aggregate(ra, fact=round(1/res(ra)[1]), fun=mean, na.rm=TRUE)
  ra <- ra - 273.15
  ra <- ra * r
  df <- as.data.table(as.data.frame(ra, xy=TRUE, na.rm=TRUE))   # x, y + one column per day
  daycols <- setdiff(names(df), c("x","y"))
  setnames(df, daycols, paste0("d", seq_along(daycols)))
  per <- period_of(yr)                                          # period of each day
  long <- melt(df, id.vars=c("x","y"), variable.name="d", value.name="val")
  rm(df)
  long[, per := per[as.integer(sub("d","",d))]]
  long[, val := round(val, bin_choice)]
  long <- long[per >= YEAR_MIN & per <= YEAR_MAX & is.finite(val)]
  # count days per cell (x, y) x bin x period, for THIS file
  cnt <- long[, .(count = .N), by=.(x, y, val, per)]
  rm(long)
  cnt[, model := mdl]
  acc[[k]] <- cnt
  cat(sprintf("  [%2d/%2d] %-16s %s  (%d lignes)\n", k, length(nc_files), mdl,
              paste(range(yr), collapse="-"), nrow(cnt)))
  rm(ra, cnt); gc()
}
acc <- rbindlist(acc)

# ------- 4. map (x,y) to id, aggregate per model, then average across models -
# Both x,y and id are kept, one to one through the mask, so that the cell-level
# diagnostics file and the id-keyed GE outputs can be emitted from one pass.
acc <- merge(acc, identifier[, .(x, y, id)], by=c("x","y"))
acc <- acc[, .(count = sum(count)), by=.(model, x, y, id, val, per)]  # one model can cover a period through two files
# share of days per (model, cell, period)
acc[, share := count / sum(count), by=.(model, x, y, per)]
# Average ACROSS models: a bin absent for one model counts as zero, so the
avg <- acc[, .(s = sum(share)), by=.(x, y, id, val, per)]
avg[, share := s / sum(s), by=.(x, y, per)]                     # sums to 1 per (cell, period)
avg[, s := NULL]
n_models <- uniqueN(acc$model)
cat(sprintf("\nMoyenne sur %d modeles. Periodes : %s\n",
            n_models, paste(sort(unique(avg$per)), collapse=" ")))

# ------- 4bis. OUTPUT CMIP6_tas.csv, for the diagnostics --------------------
# Schema expected by the diagnostics: x, y, year, value, count (share, summing
# One unified convention: the same periods, 2015 to 2100, as the GE. Any
# diagnostic that reads this file must be re-run and re-validated when the
# convention changes, even though the model ensemble is unchanged.
cmip6_tas <- avg[, .(x, y, year = per, value = val, count = share)]
setorder(cmip6_tas, x, y, year, value)
fwrite(cmip6_tas, file.path(DIR_OUT, "CMIP6_tas.csv"))
cat(sprintf("[OK] CMIP6_tas.csv : %d lignes (diagnostics, x,y,year,value,count)\n", nrow(cmip6_tas)))

# ------- 5. DAILY OUTPUT (sc2): wide id, value, id_value, <periods> ---------
daily <- avg[, .(id, value = val, id_value = paste(id, val, sep="_"),
                 per, count = 365 * share)]                     # sums to ~365 days per cell-year
daily_wide <- dcast(daily, id + value + id_value ~ per, value.var="count")
ycols <- sort(setdiff(names(daily_wide), c("id","value","id_value")))
setcolorder(daily_wide, c("id","value","id_value", ycols))
fwrite(daily_wide, file.path(DIR_OUT, "tas_ssp2_daily.csv"))
cat(sprintf("[OK] tas_ssp2_daily.csv  : %d lignes, periodes = %s\n",
            nrow(daily_wide), paste(ycols, collapse=" ")))

# ------- 6. ANNUAL OUTPUT (sc3): long id, year, temp_mean ------------------
annual <- avg[, .(temp_mean = sum(val * share)), by=.(id, per)]
setnames(annual, "per", "year"); setorder(annual, id, year)
fwrite(annual, file.path(DIR_OUT, "tas_ssp2_annual.csv"))
cat(sprintf("[OK] tas_ssp2_annual.csv : %d lignes (long id,year,temp_mean)\n", nrow(annual)))

# ------- 7. VALIDATION -------------------------------------------------------
cat("\n=== VALIDATION ===\n")
cat(sprintf("ids        : %d\n", uniqueN(avg$id)))
cat(sprintf("temp bins  : %d a %d degC\n", min(avg$val), max(avg$val)))
sc <- daily[, .(s = sum(count)), by=.(id, per)]
cat(sprintf("count/cellule-an (attendu ~365) : min=%.1f med=%.1f max=%.1f\n",
            min(sc$s), median(sc$s), max(sc$s)))
g <- annual[, .(t = mean(temp_mean)), by=year][order(year)]
print(g)
if (all(c(2015,2095) %in% g$year))
  cat(sprintf("Rechauffement moyen {2015-2019}->{2095-2099} : %+.2f degC (SSP2-4.5 attendu ~ +2.3 a +2.8)\n",
              g[year==2095, t] - g[year==2015, t]))
cat("\nDONE.\n")
