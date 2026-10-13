##################################################################################
# 00_fetch/oh_extract_hdi.R — country x year Human Development Index from the
# gridded HDI of Kummu et al. (2018).
#
# One-time extraction, not a step of run_all.sh; its output
# output/oh_hdi_country_year.csv is shipped with the package.
#
# Reads  : data_inputs/Kummu/HDI_1990_2015_v2.nc (Kummu et al. 2018, Dryad
#          doi:10.5061/dryad.dk1j0; 5 arc-minutes, layers 1990..2015)
#          data_inputs/Kummu/polyg_adm0_dataNetMgr.gpkg (national polygons,
#          field `iso3`)
#          Neither file is listed in DATA.md: they are needed only to re-run
#          this extraction, not the chain.
# Method : coverage-fraction-weighted (area) country mean via exactextractr.
#          The HDI enters only as a residualized development control, so area
#          weighting rather than population weighting is a simplification of a
#          robustness control.
# Writes : output/oh_hdi_country_year.csv (columns iso3, year, HDI), the
#          observed years 1990-2015 only.
##################################################################################

suppressMessages({ library(terra); library(sf); library(exactextractr); library(data.table) })

ROOT    <- Sys.getenv("JMP_ROOT", ".")
HDI_NC  <- file.path(ROOT, "data_inputs", "Kummu", "HDI_1990_2015_v2.nc")
ADM0    <- file.path(ROOT, "data_inputs", "Kummu", "polyg_adm0_dataNetMgr.gpkg")
OUT     <- file.path(ROOT, "output")

cat("Loading HDI raster ...\n")
hdi <- terra::rast(HDI_NC, subds = "HDI")
if (is.na(terra::crs(hdi)) || terra::crs(hdi) == "")
  terra::crs(hdi) <- "EPSG:4326"
yrs <- 1990:2015
if (terra::nlyr(hdi) != length(yrs))
  stop(sprintf("HDI layers=%d != 26 expected", terra::nlyr(hdi)))
names(hdi) <- paste0("HDI_", yrs)
cat(sprintf("  HDI: %d layers, %d x %d cells, ext=[%s]\n",
            terra::nlyr(hdi), terra::nrow(hdi), terra::ncol(hdi),
            paste(round(as.vector(terra::ext(hdi)), 1), collapse = ", ")))

cat("Loading adm0 polygons ...\n")
adm0 <- sf::st_read(ADM0, quiet = TRUE)
adm0 <- sf::st_make_valid(adm0)
stopifnot("iso3" %in% names(adm0))
cat(sprintf("  adm0: %d features, %d with non-NA iso3\n",
            nrow(adm0), sum(!is.na(adm0$iso3))))

cat("Zonal extraction (coverage-weighted mean) ...\n")
m <- exactextractr::exact_extract(hdi, adm0, "mean", progress = FALSE)
# m: data.frame [n_features x 26], cols mean.HDI_1990 ...
m <- as.data.table(m)
setnames(m, names(m), sub("^mean\\.", "", names(m)))
m[, iso3 := adm0$iso3]
m <- m[!is.na(iso3)]

# collapse to one row per iso3 (a few iso3 appear in >1 polygon row): area-mean
long <- melt(m, id.vars = "iso3", variable.name = "var", value.name = "HDI")
long[, year := as.integer(sub("HDI_", "", var))]
long <- long[is.finite(HDI), .(HDI = mean(HDI, na.rm = TRUE)), by = .(iso3, year)]
setorder(long, iso3, year)

outf <- file.path(OUT, "oh_hdi_country_year.csv")
fwrite(long[, .(iso3, year, HDI)], outf)
cat(sprintf("Wrote %s : %d rows, %d distinct iso3, years %d-%d\n",
            outf, nrow(long), length(unique(long$iso3)),
            min(long$year), max(long$year)))
cat(sprintf("  HDI range: [%.3f, %.3f], mean %.3f\n",
            min(long$HDI), max(long$HDI), mean(long$HDI)))
