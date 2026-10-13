# _dependencies.R -- declares to renv the packages it cannot detect (lot 84, 03/10/2026).
# This file is never run. renv::dependencies() reads it.
#
# 02_GE/02a_calibration_at_1deg.r loads its packages through a vector:
#   lapply(libraries, require, character.only = TRUE)
# renv does not see that form, and under renv a require() on a missing package returns
# FALSE without an error, so 02a would carry on without it. Keep this list identical to
# the vector `libraries` of 02a (l.26-29).
library(rhdf5)
library(akima)
library(hdf5r)
library(smoothr)
library(spData)
library(tidyterra)
library(DescTools)
library(rnaturalearth)
library(fuzzyjoin)
library(patchwork)
library(ggplot2)
library(data.table)
library(R.matlab)
library(readxl)
library(sf)
library(raster)
library(Matrix)
library(dplyr)
library(tidyr)
library(terra)
library(countrycode)
library(exactextractr)
library(bigmemory)
library(RcppCNPy)

# terra::writeCDF() needs ncdf4, which terra only suggests (renv does not lock suggestions).
# Used by tools/test_isimip_annee_utc.R; found when that test failed under renv (lot 84).
library(ncdf4)

# rnaturalearth::ne_countries(scale = "medium") reads its polygons from rnaturalearthdata,
# which rnaturalearth only suggests and loads at run time (requireNamespace). Called by 01c,
# 01d, 01e, 02a and 03_inference/step3, step5, step7; found when 01c failed on a fresh clone
# (lot 85d). No call uses scale = "large", so rnaturalearthhires is not needed.
library(rnaturalearthdata)
