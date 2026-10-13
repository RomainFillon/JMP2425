##################################################################################
# lib_eggw.R -- access to the Cruz & Rossi-Hansberg replication package
#
# Reads data files from the published CRH replication archive. The archive is
# public and citable; it is not redistributed with this package.
#   Cruz, J. L. & Rossi-Hansberg, E., "The Economic Geography of Global
#   Warming", replication package, Zenodo, DOI 10.5281/zenodo.7681424.
#
# Reads:  data_inputs/EGGW_Replication_Package.zip
# Writes: intermediate_outputs/eggw_mat/ (cache of extracted .mat files)
# Needs from 00_parameters.R: EGGW_ZIP, EGGW_PREFIX, DIR_IN, DIR_TMP.
#
# CSV files are read from inside the archive without unpacking it. Only .mat
# files are extracted, because hdf5r and R.matlab open a path and cannot read
# from a connection.
#
# The Zenodo download is a zip inside a zip: "EGGW Replication Package
# Final.zip" (2.3 GB) holds a single entry, EGGW_Replication_Package.zip. R
# cannot read through two levels, so EGGW_ZIP points at the INNER archive,
# which must have been extracted from the outer one first.
##################################################################################

##################################################################################
# 1. ARCHIVE AVAILABILITY
##################################################################################

## No automatic download. Pulling 2.3 GB from a running script is a surprise
## rather than a service, and an interrupted transfer leaves a truncated file
## that unzip fails to read without saying why. The user is told what to do.
require_eggw <- function() {
  if (file.exists(EGGW_ZIP)) return(invisible(EGGW_ZIP))
  stop("⛔ Archive CRH absente : ", EGGW_ZIP, "\n",
       "   1. Télécharger « EGGW Replication Package Final.zip » (2,3 Go) sur\n",
       "      https://zenodo.org/records/7681424\n",
       "   2. En extraire EGGW_Replication_Package.zip vers ", DIR_IN, "\n",
       "   (l'archive Zenodo est un zip dans un zip ; c'est l'INTERNE qu'il faut)")
}

##################################################################################
# 2. READING FROM INSIDE THE ARCHIVE
##################################################################################

## Path inside the archive. The __MACOSX/ entries are compression noise: they
## carry the same names prefixed with ._ and must never be read.
eggw_path <- function(rel) file.path(EGGW_PREFIX, rel)

## `rel` is the path relative to the package root, e.g. "Data/C.csv".
##
## The CSV files in this archive have no final newline. readLines() then drops
## the last line, so the file is read as raw bytes and split here instead.
## Nothing is lost whatever the line terminator.
eggw_csv <- function(rel, ...) {
  require_eggw()
  entry <- eggw_path(rel)

  z <- utils::unzip(EGGW_ZIP, list = TRUE)
  if (!entry %in% z$Name)
    stop("⛔ ", rel, " absent de l'archive (entree attendue : ", entry, ")")

  con <- unz(EGGW_ZIP, entry, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  raw <- readBin(con, "raw", n = z$Length[z$Name == entry] + 1L)
  close(con)
  if (!length(raw)) stop("⛔ ", rel, " est vide dans l'archive.")

  txt <- strsplit(rawToChar(raw), "\r?\n")[[1]]

  ## Completeness check. A file with no final newline holds one more line than
  ## it has newline bytes, and strsplit returns no trailing empty element. Both
  ## terminations therefore give the same count, so any mismatch means the read
  ## was truncated. Without this check a truncated read is silent: the caller
  ## gets a shorter matrix and any downstream indexing shifts.
  n_nl   <- sum(raw == as.raw(10L))
  n_reel <- n_nl + as.integer(utils::tail(raw, 1) != as.raw(10L))
  if (length(txt) != n_reel)
    stop(sprintf("⛔ %s : %d lignes decoupees, %d dans l'archive.", rel,
                 length(txt), n_reel))

  data.table::fread(text = txt, ...)
}

## The .mat readers need a real file, so those are extracted one at a time and
## cached. The cache belongs to intermediate_outputs/ and should be gitignored:
## it can always be rebuilt from the archive.
eggw_mat <- function(nm) {
  require_eggw()
  dst <- file.path(DIR_TMP, "eggw_mat", nm)
  if (!file.exists(dst)) {
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  extraction de %s...\n", nm))
    ## junkpaths = TRUE flattens the output; without it the file would land
    ## under eggw_mat/EGGW_Replication_Package/Data/.
    unzip(EGGW_ZIP, files = eggw_path(file.path("Data", nm)),
          exdir = dirname(dst), junkpaths = TRUE)
    ## unzip(files=) raises no error when the name matches nothing: it extracts
    ## zero files silently. Hence this test.
    if (!file.exists(dst))
      stop("⛔ ", nm, " introuvable dans l'archive — arborescence modifiée ?")
  }
  dst
}

##################################################################################
# 3. READING AN HDF5 (v7.3) .mat FILE
#
# The .mat files in this archive are not all in the same format. H0_areal,
# pop_Gecon, wage_Gecon and trmult_reduced are v7.3, hence HDF5 and readable by
# hdf5r; CO2_EDGAR and map_grid are v7 and need R.matlab::readMat. The format
# cannot be guessed from the name, and using the wrong reader returns
# "file signature not found".
##################################################################################

read_mat_h5 <- function(path, var, dims) {
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  if (length(dims) == 2L) f[[var]][1:dims[1], 1:dims[2]]
  else                    f[[var]][1:dims[1], 1:dims[2], 1:dims[3]]
}
