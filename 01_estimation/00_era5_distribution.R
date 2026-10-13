##################################################################################
# 00_era5_distribution.R -- daily ERA5 temperature distributions.
#
# ONE PASS OVER THE DATA, TWO OUTPUTS:
#
#   ERA_tas_distribution_grid_annual.csv   x, y, year, value, count
#     count = NUMBER OF DAYS the cell spends in the bin that year.
#     Consumed by the Cruz and Rossi-Hansberg triangulation, whose fundamentals
#     are at one degree.
#
#   ERA_tas_annual.csv                     GID_0, year, value, count, share
#     count = SUM OF POPULATION SHARES of the country in the bin that year.
#     The single climate input of both arms of both damage functions.
#
# WHY ONE PASS. The expensive work -- reading netCDF, rotating, aggregating to
# one degree, converting units, rounding -- is IDENTICAL for the two outputs.
# Only the last step differs. Two passes would reread about 30 GB for nothing,
# and, more importantly, they could DIVERGE if one is re-run without the other.
#
# THE WEIGHTING IS NOT A DETAIL. At country level, `count` is not a number of
# days: it is a sum of population shares. An empty Siberian cell and a Moscow
# cell do not weigh the same in a country exposure. On the grid there is nothing
# to weight; the weighting belongs to the cell-to-country step, not to this one.
#
# THE PROCESSING PLAN IS MEASURED HERE, from the files themselves: rotation,
# aggregation factor, unit conversion and the source of the year. A parameter
# derived from the data cannot diverge from it, whereas a parameter read from a
# side file can, and four instrument traps all returned a plausible value where
# the information was in fact missing.
#
# Reads:  ERA5 daily netCDF at 0.25 degree, one file per year,
#         data_inputs/gadm36_shp/gadm_simplified.shp, the land mask,
#         intermediate_outputs/poids_par_pays.csv, the population weights (lot 51)
# Writes: intermediate_outputs/ERA_tas_distribution_grid_annual.csv,
#         intermediate_outputs/ERA_tas_annual.csv
# Needs:  00_population_weights.R must have run.
#
# Runs on the server: the netCDF files are not shipped with the package.
##################################################################################

suppressMessages({
  library(sf); library(terra); library(data.table)
  library(parallel); library(doParallel); library(foreach)
})

t_start <- proc.time()

##################################################################################
# 0. CHEMINS ET PLAN
##################################################################################

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

# On the server, the raw netCDF files live under the directory given by ERA5_DIR.
# Set ERA5_DIR to point there.

## DIR_IN and DIR_TMP are R objects, defined by 00_parameters.R above. JMP_ROOT
## is the name of the ENVIRONMENT VARIABLE: used bare, R looks for a symbol of
## that name, does not find one, and raises "object 'JMP_ROOT' not found", an
## error that says nothing about paths.
DIR_NC   <- Sys.getenv("ERA5_DIR", file.path(DIR_IN, "ERA"))

## The two ERA CSV files are INTERMEDIATES, not paper outputs: 01a, 01b, 01d and
## 01e all read them back from intermediate_outputs/. Leaving DIR_OUT on output/
## would write them where nothing reads them, and nothing would report it.
DIR_OUT  <- Sys.getenv("ERA5_OUT", DIR_TMP)
DIR_YEAR <- file.path(DIR_OUT, "by_year")
for (d in c(DIR_OUT, DIR_YEAR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

OUT_GRID  <- file.path(DIR_OUT, "ERA_tas_distribution_grid_annual.csv")
OUT_CTRY  <- file.path(DIR_OUT, "ERA_tas_annual.csv")


## =============================================================================
## PROCESSING PLAN -- MEASURED HERE, NOT READ FROM A FILE
##
## This block replaces era5_plan.csv, which carried four parameters -- rotation,
## aggregation factor, unit conversion, cadence -- measured by a separate
## inspection script. Two defects made that arrangement worse than useless:
##
##   - the producer of the plan was not in the repository, so the plan could
##     not be regenerated and the whole chain stopped on a missing file;
##   - the plan also carried a YEAR COVERAGE (1985-2008) that this script never
##     reads: it scans whatever .nc files sit in DIR_NC. A stale plan could
##     therefore describe a different series than the one being processed, and
##     nothing would say so.
##
## The four parameters are now measured on the first of the files that will
## actually be read (see the warning below: only that one). A parameter derived
## from the data cannot go stale.
## =============================================================================

files <- list.files(DIR_NC, pattern = "\\.nc$", full.names = TRUE)
if (!length(files)) stop("[STOP] no .nc files under ", DIR_NC)

measure_plan <- function(f) {
  r <- terra::rast(f)
  e <- as.vector(terra::ext(r))
  u <- tryCatch(terra::units(r)[1], error = function(e) NA_character_)

  ## The DECLARED unit wins over the values. minmax() often returns NA on a
  ## compressed netCDF, and relying on it once led to concluding "already in
  ## celsius" for files in kelvin -- a silent 273-degree shift.
  conv <- if (!is.na(u) && nzchar(u) && toupper(u) %in% c("K", "KELVIN")) TRUE
          else if (!is.na(u) && nzchar(u)) FALSE
          else { v <- suppressWarnings(terra::minmax(r[[1]]))
                 if (all(is.finite(v))) v[1] > 150 else NA }
  if (is.na(conv))
    stop("[STOP] unit undetermined for ", basename(f),
         ": the conversion would be a guess, not a measurement.")

  ## Cadence is read off the LAYER COUNT, which depends on no attribute:
  ## 365 layers over one year are daily, whatever the timestamps claim.
  ## An hourly series preserves the MEANS -- so any guard on the mean would
  ## pass -- while inflating the SHAPE of the distribution through the diurnal
  ## cycle. That is fatal to a specification whose object IS the shape.
  n <- terra::nlyr(r)
  if (!(n %in% c(365L, 366L)))
    stop(sprintf("[STOP] %s has %d layers; 365 or 366 expected.\n",
                 basename(f), n),
         "  Check that the download used the Daily mean statistic.")

  list(rotation = e[2] > 180.5,          # 0:360 extent -> rotation required
       conv_K   = conv,
       agg_fac  = round(1 / terra::res(r)[1]),
       res      = terra::res(r)[1],
       nlyr     = n,
       unit     = u)
}

P        <- measure_plan(files[1])
ROTATION <- P$rotation
CONV_K   <- P$conv_K
AGG_FAC  <- P$agg_fac

cat(sprintf("  plan measured on %s\n", basename(files[1])))
cat(sprintf("    resolution %.2f deg -> aggregation factor %d\n", P$res, AGG_FAC))
cat(sprintf("    longitudes %s -> rotation %s\n",
            if (ROTATION) "0:360" else "-180:180", ROTATION))
cat(sprintf("    unit '%s' -> conversion %s\n",
            ifelse(is.na(P$unit), "(undeclared)", P$unit),
            if (CONV_K) "K -> C" else "none"))
cat(sprintf("    %d layers -> daily\n", P$nlyr))

## ⚠️ The plan is measured on the FIRST file only, and NO other file is checked
## against it: process_year() applies ROTATION, CONV_K and AGG_FAC as measured
## here. A file of another resolution, unit or longitude convention would be
## processed with this plan without a word. (Corrected 27/09/2026: this comment
## used to promise a per-file check that the code never made. The check itself
## is on the post-submission list: it is a server step.)

## Rounded to the whole degree: the same bins as the CMIP6 file, so that the
## projections read on the same temperature grid as the estimation.
BIN_DIGITS <- 0L
## Layers processed at once. Beyond this, the memory footprint of a global
## raster at 0.25 degree becomes the binding constraint, not the computation.
BATCH_SIZE <- 10L
N_CORES    <- max(1L, min(detectCores() %/% 2L, 8L))

cat(sprintf("=== ERA5 | %d coeurs | rotation=%s | conversion=%s | agg=%d ===\n",
            N_CORES, ROTATION, CONV_K, AGG_FAC))
cat(sprintf("    -> %s\n    -> %s\n", basename(OUT_GRID), basename(OUT_CTRY)))

##################################################################################
# 1. MASQUE TERRESTRE
#
# terra objects do not serialise across processes, so the mask is written to disk
# and each process reloads it. That is the reason for the intermediate .tif.
##################################################################################

MASK_FILE <- file.path(DIR_OUT, "land_mask.tif")
if (!file.exists(MASK_FILE)) {
  cat("Construction du masque terrestre...\n")
  shp <- st_read(file.path(DIR_IN, "gadm36_shp/gadm_simplified.shp"), quiet = TRUE)
  tpl <- rast(ext = ext(-180, 180, -90, 90), resolution = 1, crs = "EPSG:4326")
  values(tpl) <- NA
  writeRaster(rasterize(vect(shp), tpl, field = 1, background = NA),
              MASK_FILE, overwrite = TRUE)
  rm(shp); gc()
}

##################################################################################
# 2. POIDS DE POPULATION
#
# Pre-filtered by quinquennial period, which turns the match into a plain join on
# (x, y) instead of a rolling join on an interval. That is what makes the match
# tractable over 35 years of daily data.
#
# The period retained for a year is the LAST one starting before it: the weights
# of 1995-1999 apply to 1995 through 1999 and no further.
##################################################################################

## Lot 51 (28/09/2026) : poids_par_pays.csv (00d_poids_population.py) remplace
## lengths_country.csv, qui comptait la population ENTIERE de chaque cellule pour
## chaque pays qui la touche (RECAP 3.231). Colonne `part_du_pays` : population de
## l'intersection cellule x pays sur la population du pays, somme 1 PAR PAYS.
## Le chemin LOCAL equivalent (sans les .nc) est 00e_era5_pays_local.py.
LC_PATH <- file.path(DIR_TMP, "poids_par_pays.csv")
if (!file.exists(LC_PATH))
  stop("[STOP] poids_par_pays.csv absent : ", LC_PATH, "\n",
       "  Il vient de 00d_poids_population.py, a lancer en local, puis a\n",
       "  remonter sur le serveur.")

lc <- fread(LC_PATH, na.strings = NULL)[GID_0 != ""]
lc[, year := as.integer(year)]
LC_PERIODS <- sort(unique(lc$year))
LC_BY_PERIOD <- lapply(LC_PERIODS, function(q) {
  d <- lc[year == q, .(x, y, GID_0, part_du_pays)]
  setkey(d, x, y); d
})
names(LC_BY_PERIOD) <- as.character(LC_PERIODS)
cat(sprintf("  poids : %d pays | %d pixels | periodes %s\n",
            uniqueN(lc$GID_0), uniqueN(lc[, paste(x, y)]),
            paste(LC_PERIODS, collapse = ", ")))
rm(lc); gc()

##################################################################################
# 3. PROCESSING ONE YEAR
##################################################################################

## Rotation is applied only when the extent calls for it. A raster already
## centred on zero, rotated again, would be shifted by half a turn, silently. The
## expectation is one thing; the extent is the fact.
safe_rotate <- function(r) {
  e <- ext(r)
  if (e[2] > 180.5) rotate(r) else r
}

## These files carry NO usable timestamp: terra::time() returns the Unix epoch,
## and the layer names carry layer indices rather than seconds. The year is
## therefore taken from the FILE NAME, which is a dependency on a naming
## convention and is declared as such.
year_of <- function(path) {
  m <- regmatches(basename(path), regexpr("(19|20)\\d{2}", basename(path)))
  if (!length(m)) stop("Annee indeterminable pour : ", basename(path))
  as.integer(m[1])
}

process_year <- function(f) {
  library(terra); library(data.table)
  t0 <- proc.time()[["elapsed"]]
  yr <- year_of(f)

  ## Resume after an interruption: a year already done is skipped. Over 35 years
  ## in parallel, a cut no longer loses everything.
  fg <- file.path(DIR_YEAR, sprintf("grid_%d.csv", yr))
  fc <- file.path(DIR_YEAR, sprintf("country_%d.csv", yr))
  if (file.exists(fg) && file.exists(fc)) return(c(fg, fc))

  mask_r <- rast(MASK_FILE)
  tpl <- rast(ext = ext(-180, 180, -90, 90), resolution = 1, crs = "EPSG:4326")
  values(tpl) <- NA

  ## Weights: the last period starting before this year.
  qs <- as.integer(names(LC_BY_PERIOD))
  q  <- max(qs[qs <= yr])
  lc_yr <- LC_BY_PERIOD[[as.character(q)]]

  res <- tryCatch({
    tas <- rast(f)
    ## Aggregated to one degree by SPATIAL MEAN, so that the mean of the cell
    ## daily distribution reproduces the one-degree mean. That is what makes the
    ## "only the metric changes" contrast clean on the CRH side.
    if (AGG_FAC > 1L)
      tas <- terra::aggregate(tas, fact = AGG_FAC, fun = mean, na.rm = TRUE)
    tas <- mask(resample(safe_rotate(tas), tpl, method = "near"), mask_r)

    acc_g <- acc_c <- NULL
    for (ix in split(seq_len(nlyr(tas)), ceiling(seq_len(nlyr(tas)) / BATCH_SIZE))) {
      df <- as.data.frame(tas[[ix]], xy = TRUE, na.rm = TRUE)
      if (!nrow(df)) next
      setDT(df)
      dc <- paste0("d", seq_along(ix))
      setnames(df, c("x", "y", dc))
      df[, (dc) := lapply(.SD, function(v)
        round(if (CONV_K) v - 273.15 else v, BIN_DIGITS)), .SDcols = dc]
      lg <- melt(df, id.vars = c("x", "y"), measure.vars = dc,
                 variable.name = "day", value.name = "value")
      rm(df)

      ## --- OUTPUT 1: grid. count = number of DAYS. -------------------------
      g <- lg[, .(count = .N), by = .(x, y, value)]
      acc_g <- rbindlist(list(acc_g, g))[, .(count = sum(count)), by = .(x, y, value)]

      ## --- OUTPUT 2: country. count = sum of POPULATION SHARES. ------------
      ## A pixel belongs to several countries on borders, so the join
      ## duplicates it, deliberately: it contributes to each in proportion.
      lg[, `:=`(x = as.numeric(x), y = as.numeric(y))]
      j <- lc_yr[lg, on = .(x, y), nomatch = 0L]
      if (nrow(j)) {
        cc <- j[, .(count = sum(part_du_pays, na.rm = TRUE)), by = .(GID_0, value)]
        acc_c <- rbindlist(list(acc_c, cc))[, .(count = sum(count)), by = .(GID_0, value)]
      }
      rm(lg, g); gc()
    }
    list(grid = acc_g, country = acc_c)
  }, error = function(e) {
    message(sprintf("[%d] ERREUR : %s", yr, conditionMessage(e))); NULL
  })

  if (is.null(res) || is.null(res$grid) || !nrow(res$grid)) {
    message(sprintf("[%d] aucune donnee produite", yr)); return(NULL)
  }
  res$grid[, year := yr]; fwrite(res$grid, fg)
  if (!is.null(res$country) && nrow(res$country)) {
    res$country[, year := yr]; fwrite(res$country, fc)
  } else {
    ## A mechanism that does nothing must say so. An empty country output means
    ## the weight join matched no pixel at all.
    message(sprintf("[%d] ⚠ SORTIE PAYS VIDE : aucun pixel apparie aux poids", yr))
    return(NULL)
  }
  message(sprintf("[%d] %.0f s | grille %d lignes | pays %d lignes", yr,
                  proc.time()[["elapsed"]] - t0, nrow(res$grid), nrow(res$country)))
  c(fg, fc)
}

##################################################################################
# 4. EXECUTION PARALLELE
##################################################################################

years <- vapply(files, year_of, 1L)
cat(sprintf("\n%d fichiers : %d-%d\n", length(files), min(years), max(years)))

cl <- makeCluster(N_CORES); registerDoParallel(cl)
clusterExport(cl, c("DIR_YEAR", "MASK_FILE", "LC_BY_PERIOD", "AGG_FAC", "CONV_K",
                    "BIN_DIGITS", "BATCH_SIZE", "safe_rotate", "year_of"),
              envir = environment())
res <- foreach(f = files, .errorhandling = "pass") %dopar% process_year(f)
stopCluster(cl)

done <- vapply(res, function(x) is.character(x) && length(x) == 2L &&
                                all(file.exists(x)), TRUE)
cat(sprintf("\n%d/%d annees produites\n", sum(done), length(files)))
if (!all(done)) {
  cat(sprintf("  ⚠ echouees : %s\n", paste(sort(years[!done]), collapse = ", ")))
  cat("    Relancer : les annees deja faites seront sautees.\n")
}
if (!any(done)) stop("[STOP] aucune annee produite -- rien a fusionner.")

##################################################################################
# 5. FUSION
##################################################################################

fg <- list.files(DIR_YEAR, pattern = "^grid_\\d{4}\\.csv$",    full.names = TRUE)
fc <- list.files(DIR_YEAR, pattern = "^country_\\d{4}\\.csv$", full.names = TRUE)

grid <- rbindlist(lapply(fg, fread))
setcolorder(grid, c("x", "y", "year", "value", "count"))
setorder(grid, x, y, year, value)
fwrite(grid, OUT_GRID)

ctry <- rbindlist(lapply(fc, fread))
## The share is computed HERE, not downstream: two consumers would otherwise
## recompute it, and one of the two would eventually diverge.
ctry[, share := count / sum(count, na.rm = TRUE), by = .(GID_0, year)]
setcolorder(ctry, c("GID_0", "year", "value", "count", "share"))
setorder(ctry, GID_0, year, value)
fwrite(ctry, OUT_CTRY)

##################################################################################
# 6. CONTROLES
##################################################################################

cat("\n=== CONTROLES ===\n")

cat(sprintf("\n--- %s ---\n", basename(OUT_GRID)))
cat(sprintf("  %d lignes | %d cellules | %d annees (%d-%d) | %d bins\n",
            nrow(grid), uniqueN(grid[, paste(x, y)]), uniqueN(grid$year),
            min(grid$year), max(grid$year), uniqueN(grid$value)))
nd <- grid[, .(n = sum(count)), by = .(x, y, year)]
cat(sprintf("  jours par cellule-annee : mediane %.0f | min %.0f | max %.0f\n",
            median(nd$n), min(nd$n), max(nd$n)))
## A median far from 365 would signal HOURLY data taken for daily data, the
## defect that preserves the means and inflates the SHAPE.
if (median(nd$n) < 360 || median(nd$n) > 370)
  cat("  🚨 mediane hors [360, 370] : verifier la cadence des .nc.\n")

cat(sprintf("\n--- %s ---\n", basename(OUT_CTRY)))
cat(sprintf("  %d lignes | %d pays | %d annees (%d-%d) | %d bins\n",
            nrow(ctry), uniqueN(ctry$GID_0), uniqueN(ctry$year),
            min(ctry$year), max(ctry$year), uniqueN(ctry$value)))
s <- ctry[, .(tot = sum(share)), by = .(GID_0, year)]
cat(sprintf("  parts sommant a 1 : %d ecart(s) > 1e-6 sur %d pays-annees\n",
            sum(abs(s$tot - 1) > 1e-6), nrow(s)))

## Plausibility check: the weighted means must be ordered as geography orders
## them. An inversion would signal a failed rotation, the error that never raises
## an exception.
m <- ctry[, .(T = weighted.mean(value, count)), by = GID_0][order(-T)]
cat("  moyennes ponderees, extremes :\n")
print(rbind(head(m, 3), tail(m, 3)))

cat(sprintf("\n[write] %s\n[write] %s\n", OUT_GRID, OUT_CTRY))
cat(sprintf("=== termine en %.1f min ===\n",
            (proc.time() - t_start)[["elapsed"]] / 60))
