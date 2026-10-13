#####################################################################################
# mig_internal_cross.R -- T1 de tools/REEMISSION_FINALE.md (lots 70-72, 01/10/2026).
#
# La masse BRUTE transfrontaliere du bloc interne de mu_0, et sa part du bloc
# interne, pour \GEmigInternalCrossM et \GEmigInternalCrossShare (E4).
#
# Objet (rapport_lot70 section 4) : somme des flux internes dont l'origine et la
# destination, attribuees aux pays par parts, sont deux pays differents, soit la
# masse hors diagonale de W' F W, F = bloc interne zone x zone (STEP 11-12 de 02b),
# W = part de chaque pays dans la population de la zone (part_z de STEP 10).
# Niveau ZONE, AVANT STEP 13 (domaine) et AVANT STEP 15 (plafond).
#
# Lecture seule. Lit les RDS que 02b ecrit en STEP 9 (migration_matrix,
# M_intra_zc, list_zone1_id_vec, pop_lookup) et en STEP 12 (list_zone1_id), plus
# les CSV que 02b lit. Recalcule STEP 10-11 avec les regles de 02b (lot 68) : ce
# recalcul REPRODUIT l'allocation de 02b au lieu de la lire, d'ou deux assertions :
#   A1  la masse interne recalculee = celle que 02b a journalisee en STEP 11
#       (argument --masse-step11, en M, a 0,005 M pres : le journal arrondit a 2 decimales) ;
#   A2  le solde interne net par pays, porte par le meme chemin que 02b (STEP 13
#       domaine, STEP 15 zeros et plafond, attribution par part_de_la_cellule),
#       = colonne `internal` de mig_witness_by_country.csv, a 1 personne pres.
#   Connu-faux de A2, a chaque passage : la Syrie au solde interne inverse -> refus.
# Ecrit intermediate_outputs/mig_internal_cross.csv. N'emet aucune macro.
#
# Usage : Rscript tools/mig_internal_cross.R --masse-step11 89.65
#####################################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
})

.a <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
ROOT <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
TMP  <- file.path(ROOT, "intermediate_outputs")
args <- commandArgs(trailingOnly = TRUE)
i_m  <- match("--masse-step11", args)
if (is.na(i_m)) stop("[STOP] --masse-step11 <M> requis (STEP 11 du journal de 02b).")
MASSE_STEP11 <- as.numeric(args[i_m + 1])

## Les deux parametres de 02b (lot 68), lus dans le fichier, pas recopies.
src02b <- readLines(file.path(ROOT, "02_GE", "02b_mig_flows.R"))
lire_param <- function(nom) {
  l <- grep(paste0("^\\s*", nom, "\\s*="), src02b, value = TRUE)
  if (length(l) != 1) stop("[STOP] parametre ", nom, " introuvable ou ambigu dans 02b")
  as.numeric(sub(",.*$", "", sub("^[^=]*=\\s*", "", l)))
}
MIN_PART <- lire_param("internal_min_part")
PERIOD   <- lire_param("period_years")
cat(sprintf("[02b] internal_min_part = %g ; period_years = %g (lus dans 02b_mig_flows.R)\n",
            MIN_PART, PERIOD))

md5 <- function(p) unname(tools::md5sum(p))
for (f in c("migration_matrix.rds", "M_intra_zc.rds", "list_zone1_id_vec.rds",
            "list_zone1_id.rds", "pop_lookup.rds", "mig_witness_by_country.csv",
            "poids_par_cellule.csv", "identifier_xy_final.csv", "population0.csv"))
  cat(sprintf("[md5] %-28s %s  %s\n", f, md5(file.path(TMP, f)),
              format(file.mtime(file.path(TMP, f)), "%d/%m %H:%M")))

t0 <- Sys.time()
mm         <- readRDS(file.path(TMP, "migration_matrix.rds"))
M_intra_zc <- readRDS(file.path(TMP, "M_intra_zc.rds"))
zvec       <- readRDS(file.path(TMP, "list_zone1_id_vec.rds"))
lz         <- readRDS(file.path(TMP, "list_zone1_id.rds"))
cat(sprintf("[load] migration_matrix %d x %d, nnz %d (%.0f s)\n",
            nrow(mm), ncol(mm), length(mm@x), as.numeric(Sys.time() - t0, units = "secs")))
stopifnot(nrow(mm) == length(zvec))

#####################################################################################
# STEP 10 (copie de 02b l.441-458)
#####################################################################################
net_flow <- data.frame(zone1_id = zvec, net_flow = as.numeric(colSums(mm) - rowSums(mm)))
M <- left_join(M_intra_zc, net_flow, by = "zone1_id")
M$net_flow[is.na(M$net_flow)] <- 0
M <- M %>% group_by(zone1_id) %>%
  mutate(part_z = { s <- sum(pop_zc, na.rm = TRUE); if (s > 0) pop_zc / s else 0 }) %>%
  ungroup()
M$stock <- M$value_zc - M$net_flow * M$part_z
M <- M[, c("zone1_id", "GID_0", "stock", "pop_zc", "part_z")]

#####################################################################################
# STEP 11 (copie de 02b l.473-515)
#####################################################################################
alloc <- function(d) {
  d <- d[d$part_z >= MIN_PART, ]
  if (nrow(d) == 0) return(NULL)
  d$pop <- d$pop_zc; d$pop[is.na(d$pop)] <- 0
  tp <- sum(d$stock[d$stock > 0], na.rm = TRUE)
  tn <- sum(d$stock[d$stock < 0], na.rm = TRUE)
  tpop <- sum(d$pop, na.rm = TRUE)
  ms <- min(max(tp, abs(tn)), tpop)
  if (ms == 0) return(NULL)
  Po <- if (tpop > 0) d$pop / tpop else rep(0, nrow(d))
  sp <- pmax(d$stock, 0); tsp <- sum(sp, na.rm = TRUE)
  Pi <- if (tsp > 0) sp / tsp else Po
  fm <- ms * (Po %o% Pi); diag(fm) <- 0
  nz <- which(fm != 0, arr.ind = TRUE)
  if (nrow(nz) == 0) return(NULL)
  data.frame(orig = d$zone1_id[nz[, 1]], dest = d$zone1_id[nz[, 2]], flow = fm[nz])
}
Fl <- bind_rows(lapply(split(M, M$GID_0), alloc))
Fl <- Fl[Fl$flow != 0 & !is.na(Fl$orig) & !is.na(Fl$dest) & !is.nan(Fl$flow), ]
Fl <- Fl %>% group_by(orig, dest) %>% summarise(flow = sum(flow), .groups = "drop")
masse_interne <- sum(Fl$flow)
cat(sprintf("[STEP 11] masse interne recalculee : %.4f M (journal 02b : %.2f M)\n",
            masse_interne / 1e6, MASSE_STEP11))
if (abs(masse_interne / 1e6 - MASSE_STEP11) > 0.005)
  stop("[STOP] A1 : la masse interne recalculee ne reproduit pas STEP 11 de 02b.")
cat("[ok] A1 : STEP 11 reproduit.\n")

#####################################################################################
# LA MESURE : masse hors diagonale de W' F W (niveau zone, avant STEP 13 et 15)
#####################################################################################
zones <- sort(unique(c(M$zone1_id, Fl$orig, Fl$dest)))
ctys  <- sort(unique(M$GID_0[M$part_z > 0]))
Wz <- sparseMatrix(i = match(M$zone1_id[M$part_z > 0], zones),
                   j = match(M$GID_0[M$part_z > 0], ctys),
                   x = M$part_z[M$part_z > 0], dims = c(length(zones), length(ctys)),
                   dimnames = list(NULL, ctys))
Fz <- sparseMatrix(i = match(Fl$orig, zones), j = match(Fl$dest, zones), x = Fl$flow,
                   dims = c(length(zones), length(zones)))
## Une zone qui emet ou recoit doit etre attribuee en entier (somme des parts = 1).
zs <- rowSums(Wz)
act <- unique(c(match(Fl$orig, zones), match(Fl$dest, zones)))
if (any(abs(zs[act] - 1) > 1e-9))
  stop("[STOP] des zones portant des flux internes n'ont pas des parts sommant a 1.")
CC <- as.matrix(t(Wz) %*% Fz %*% Wz)
cross <- sum(CC) - sum(diag(CC))
net_z <- colSums(CC) - rowSums(CC)
cat(sprintf("[mesure] masse brute transfrontaliere : %.4f M sur %.4f M = %.3f %%\n",
            cross / 1e6, masse_interne / 1e6, 100 * cross / masse_interne))
cat(sprintf("[mesure] desequilibre net niveau zone, sum |solde| / 2 : %.4f M\n",
            sum(abs(net_z)) / 2e6))

#####################################################################################
# A2 : le meme bloc interne porte par le chemin de 02b (STEP 12-15), contre le temoin
#####################################################################################
zones_valides <- sort(unique(lz$zone1_id))
zones_to_keep <- sort(unique(lz$row[!is.na(lz$id)]))
z1_keep <- lz$zone1_id[match(zones_to_keep, lz$row)]
z1_keep <- sort(z1_keep[!is.na(z1_keep)])
stopifnot(length(z1_keep) == length(zones_to_keep))
cat(sprintf("[STEP 13] list_zone1_id_vec trie : %s (02b etiquette les lignes gardees par le zone1_id TRIE)\n",
            !is.unsorted(zvec)))
z1_to_id <- setNames(lz$id, as.character(lz$zone1_id))

identifier <- read.csv(file.path(TMP, "identifier_xy_final.csv"))
max_id <- max(identifier$id)
pop_id <- read.csv(file.path(TMP, "population0.csv")); pop_id$value[is.na(pop_id$value)] <- 0
pop_vec <- numeric(max_id); pop_vec[pop_id$id] <- pop_id$value
mask <- as.numeric(pop_vec > 0)

## bloc international garde : sommes de lignes sur les colonnes non masquees, sans copie
id_row <- as.integer(z1_to_id[as.character(z1_keep)])       # id de la k-ieme ligne gardee
v <- numeric(nrow(mm)); v[zones_to_keep] <- mask[id_row]
Mrow <- numeric(max_id); Mrow[id_row] <- as.numeric(mm %*% v)[zones_to_keep]
rm(mm); invisible(gc())

## bloc interne garde, re-indexe sur id (STEP 12-14)
fk <- Fl[Fl$orig %in% zones_valides & Fl$dest %in% zones_valides &
         Fl$orig %in% z1_keep & Fl$dest %in% z1_keep, ]
fk$i <- as.integer(z1_to_id[as.character(fk$orig)])
fk$j <- as.integer(z1_to_id[as.character(fk$dest)])
fk <- fk[!is.na(fk$i) & !is.na(fk$j), ]
fk$x <- fk$flow * mask[fk$i] * mask[fk$j]                   # STEP 15 : zeros
Frow <- numeric(max_id); r <- tapply(fk$x, fk$i, sum); Frow[as.integer(names(r))] <- r
R <- (Frow + Mrow) * mask
cap <- pop_vec * PERIOD
s <- ifelse(R > cap & R > 0, cap / R, 1); s[pop_vec == 0] <- 0  # STEP 15 : plafond
out_i <- numeric(max_id); in_i <- numeric(max_id)
o <- tapply(s[fk$i] * fk$x, fk$i, sum); out_i[as.integer(names(o))] <- o
n <- tapply(s[fk$i] * fk$x, fk$j, sum); in_i[as.integer(names(n))] <- n
cat(sprintf("[STEP 13-15] bloc interne apres domaine et plafond : %.4f M (journal 02b STEP 13 : interne 85,24 M avant plafond)\n",
            sum(out_i) / 1e6))

pw <- read.csv(file.path(TMP, "poids_par_cellule.csv"), na.strings = character(0))
pw <- pw[pw$GID_0 != "", c("GID_0", "x", "y", "part_de_la_cellule")]
pw <- merge(pw, identifier[, c("x", "y", "id")], by = c("x", "y"))
wit <- read.csv(file.path(TMP, "mig_witness_by_country.csv"), na.strings = character(0))
pw <- pw[pw$GID_0 %in% wit$GID_0, ]
P <- sparseMatrix(i = pw$id, j = match(pw$GID_0, wit$GID_0), x = pw$part_de_la_cellule,
                  dims = c(max_id, nrow(wit)))
int_c <- as.numeric(crossprod(P, in_i - out_i))
ecart <- abs(int_c - wit$internal)
cat(sprintf("[A2] %d pays ; |recalcule - temoin 02b| max %.3g personnes (%s)\n",
            nrow(wit), max(ecart), wit$GID_0[which.max(ecart)]))
a2 <- function(x) wit$GID_0[abs(x - wit$internal) > 1]
bad <- int_c; k <- match("SYR", wit$GID_0); bad[k] <- -bad[k]
if (!("SYR" %in% a2(bad))) stop("[STOP] A2 ne refuse pas le connu-faux (Syrie inversee). Garde morte.")
cat(sprintf("[known-bad] Syrie au solde interne inverse (%+.3f M) : refusee, comme attendu\n", bad[k] / 1e6))
if (length(a2(int_c))) stop("[STOP] A2 : le recalcul diverge de 02b pour ", paste(a2(int_c), collapse = ", "))
cat("[ok] A2 : soldes internes par pays = mig_witness_by_country.csv a 1 personne pres.\n")

res <- data.frame(
  quantity = c("internal_mass", "internal_cross_gross", "internal_cross_share_pct",
               "internal_net_zone_half_abs", "internal_net_cell_half_abs_witness"),
  value = c(masse_interne, cross, 100 * cross / masse_interne, sum(abs(net_z)) / 2,
            sum(abs(wit$internal)) / 2),
  unit = c("persons, 5 years", "persons, 5 years", "% of internal_mass", "persons", "persons"),
  definition = c("bloc interne de mu_0 (STEP 11 de 02b), niveau zone, avant STEP 13 et 15",
                 "masse hors diagonale de W'FW, W = part_z (pop_zc), niveau zone, avant STEP 13 et 15",
                 "internal_cross_gross / internal_mass",
                 "sum_c |solde interne net de c| / 2, niveau zone, avant STEP 13 et 15",
                 "sum_c |internal| / 2 de mig_witness_by_country.csv (cellule, apres STEP 13 et 15)"),
  source_md5 = md5(file.path(TMP, "migration_matrix.rds")))
write.csv(res, file.path(TMP, "mig_internal_cross.csv"), row.names = FALSE)
cat("[write] intermediate_outputs/mig_internal_cross.csv\n")
print(res[, 1:2], row.names = FALSE, digits = 10)
