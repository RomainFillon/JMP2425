# 00_fetch/extract_netmgr_2010_2014.R — one-time extraction, not a step of run_all.sh.
# Reads the net migration raster of Niva et al. (2023), Zenodo doi:10.5281/zenodo.7997134,
# file raster_netMgr_2000_2019_5yrSum.tif (four 5-year sums), keeps the layer 2010_2014
# and writes data_inputs/Kummu/raster_netMgr_2010_2014_5yrSum.tif, read by 02_GE/02b_mig_flows.R.
# Usage: Rscript 00_fetch/extract_netmgr_2010_2014.R [input.tif] [output.tif]
suppressMessages(library(terra))
ROOT <- Sys.getenv("JMP_ROOT", ".")
a    <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(a) >= 1) a[1] else file.path(ROOT, "data_inputs/Kummu/raster_netMgr_2000_2019_5yrSum.tif")
out  <- if (length(a) >= 2) a[2] else file.path(ROOT, "data_inputs/Kummu/raster_netMgr_2010_2014_5yrSum.tif")
r    <- rast(inp)
stopifnot("2010_2014" %in% names(r))
writeRaster(r[["2010_2014"]], out, filetype = "GTiff", datatype = "FLT4S",
            NAflag = -3.4e38, gdal = "COMPRESS=LZW", overwrite = TRUE)
cat("Wrote", out, "\n")
