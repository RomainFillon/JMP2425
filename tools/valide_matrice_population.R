## valide_matrice_population.R -- applique le temoin des populations a une matrice
## (x, y, GID_0, year, pop), periode 2015 (lot 46). Code 0 si acceptee, 1 sinon.
## Usage : Rscript tools/valide_matrice_population.R <matrice.csv>
suppressMessages({ library(data.table) })
a <- commandArgs(trailingOnly = TRUE)
## Racine derivee de l'emplacement du script, jamais du repertoire courant (lot 51).
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
.root <- if (nzchar(Sys.getenv("JMP_ROOT"))) Sys.getenv("JMP_ROOT") else
  dirname(dirname(normalizePath(.self, winslash = "/")))
suppressMessages(source(file.path(.root, "01_estimation", "00_parameters.R")))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "tools", "temoin_population_pays.R"))
m <- fread(a[1], na.strings = NULL)[year == 2015 & GID_0 != ""]
wb <- load_wb(file.path(DIR_IN, "WB_pop_country/population_data.csv"), "population",
              period = "annual", year_min = 2015L, year_max = 2015L,
              verbose = FALSE)[, .(GID_0, population)]
tout  <- fread(a[1], na.strings = NULL)[year == 2015]
monde <- sum(tout$pop)                                   # GHSL 2015, hors pays compris
hors  <- tout[GID_0 == "", sum(pop)]                     # pixels hors de tout polygone
r <- tryCatch({ temoin_population_pays(m[, .(pop = sum(pop)), by = GID_0], monde, wb,
                                       basename(a[1]), hors = hors); 0L },
              error = function(e) { message(conditionMessage(e)); 1L })
quit(status = r)
