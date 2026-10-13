##################################################################################
# ct6_fig_drf_paper.R -- FIGURE 1 du manuscrit : les deux fans de productivite,
# avec la grille des specifications en fond.
#
# Reads:  output/MC_damage_prod_{daily,annual}_quantiles.csv   <- 01a
#         output/contribA/spec_curves_{daily,annual}.csv       <- ct1
#         output/contribA_grid.csv                             <- ct4
#         intermediate_outputs/climate_long_annual.rds         <- 01a
# Writes: output/MC_prod_{daily,annual}_fan_paper.pdf
#
# POURQUOI CE SCRIPT EXISTE, ET POURQUOI IL N'EST PAS DANS 01a
# ------------------------------------------------------------
# La figure porte desormais, derriere la courbe retenue, la courbe estimee de
# CHAQUE specification de la grille. Ces courbes sont archivees par ct1, qui
# tourne apres 01a et dont les intrants viennent de 01a. La dependance s'est
# donc inversee : la figure du papier ne peut plus etre emise a l'etape 01a.
# Elle est emise ici, une seule fois, par un seul producteur -- l'appel
# correspondant a ete retire de 01a, qui ne peut pas ecrire le meme fichier.
#
# Le canal AMENITE n'a pas de grille archivee (ct1 n'archive que prod) : ses
# fans restent produits par 01b, sans fond, et leur legende ne promet rien
# d'autre.
#
# CE QUI EST TRACE, ET SON NOM EXACT
# ----------------------------------
# Les courbes grises sont les courbes de POINT (get_effect_curve sur beta
# chapeau), pas des medianes Monte Carlo : il n'existe aucun tirage par
# specification. La legende du manuscrit ecrit « the estimated curve of each of
# the \EMPgridProdNspec specifications », et c'est cet objet-la.
#
# LE COMPTE EST ASSERTE, PAS SUPPOSE
# ----------------------------------
# Le nombre de courbes grises doit valoir le nombre de lignes de
# contribA_grid.csv -- la meme passage qui remplit \EMPgridProdNspec. Un
# panneau qui en montrerait douze sous une legende qui en annonce quinze
# serait un mensonge silencieux : le script ARRETE. Depuis le 18/09/2026 les
# deux comptes sont en outre confrontes a GRID_SPECS de 00_parameters.R, qui
# porte la decision d'auteur et ne descend d'aucun run.
##################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(splines); library(MASS)
})

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f))
    stop("Impossible de deduire la racine. Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))

CHANNEL <- "prod"

cat("=== ct6_fig_drf_paper.R : Figure 1, fans + grille ===\n")

##################################################################################
# 1. LE CLIMAT -- L'OBJET DE 01a, PAS UNE COPIE
##################################################################################
## climate_long_annual.rds est ecrit par 01a (l.414) apres les trois memes
## operations d'agregation. Le relire est la seule facon de garantir que
## l'histogramme d'exposition de la figure est celui du panel d'estimation :
## refaire les trois lignes ici creerait une SECONDE source de climat, et ce
## depot en a deja paye le prix.
f_clim <- file.path(DIR_TMP, "climate_long_annual.rds")
if (!file.exists(f_clim))
  stop("Absent : ", f_clim, "\n  Produit par 01a_productivity_drf.R. Lancer 01a d'abord.")
climate <- readRDS(f_clim); setDT(climate)
cat(sprintf("[lu] climate_long_annual.rds : %d lignes, %d pays\n",
            nrow(climate), uniqueN(climate$GID_0)))

##################################################################################
# 2. LA GRILLE -- ARCHIVE DE ct1, COMPTE ASSERTE CONTRE ct4
##################################################################################
f_grid <- file.path(DIR_OUT, "contribA_grid.csv")
if (!file.exists(f_grid))
  stop("Absent : ", f_grid, "\n  Produit par ct4_emit_tab_contribA_grid.R.")
.grid <- fread(f_grid)
N_SPEC <- nrow(.grid)
## LE COMPTE EST CONFRONTE A LA GRILLE DECIDEE, PAS SEULEMENT A LUI-MEME
## (18/09/2026). Ce script comparait le nombre de courbes archivees au nombre de
## lignes du CSV : deux objets produits par la MEME chaine. Les deux pouvaient
## donc valoir douze sous une legende qui en annonce quinze, et l'assertion
## serait passee. GRID_SPECS est le tiers : il porte la decision d'auteur, il ne
## descend d'aucun run.
if (!setequal(.grid$spec_label, GRID_LABELS) || N_SPEC != GRID_NSPEC)
  stop("[STOP] contribA_grid.csv ne porte pas GRID_SPECS : ", N_SPEC,
       " lignes contre ", GRID_NSPEC, " specifications decidees.\n",
       "  ecart : ", paste(sort(union(setdiff(.grid$spec_label, GRID_LABELS),
                                      setdiff(GRID_LABELS, .grid$spec_label))),
                           collapse = ", "), "\n",
       "  Rejouer ct1, ct2 puis ct4.")
cat(sprintf("[lu] contribA_grid.csv : %d specifications attendues par panneau (= GRID_NSPEC)\n",
            N_SPEC))

GRID <- list()
for (fr in c("daily", "annual")) {
  f <- file.path(DIR_OUT, "contribA", sprintf("spec_curves_%s.csv", fr))
  if (!file.exists(f))
    stop("Absent : ", f, "\n",
         "  Produit par ct1_delta_north_grid.R (archivage des courbes des deux\n",
         "  bras). La legende de la Figure 1 annonce la grille : sans elle, la\n",
         "  figure et sa legende diraient deux choses differentes.")
  g <- fread(f); setDT(g)
  n <- uniqueN(g$model_name)
  if (n != N_SPEC)
    stop(sprintf("bras %s : %d courbes archivees pour %d specifications dans %s.",
                 fr, n, N_SPEC, basename(f_grid)))
  ## GARDE ANTI-DEGENERESCENCE -- ECRITE APRES L'AVOIR VUE ARRIVER.
  ## Quinze specifications ne peuvent pas rendre quinze fois la meme courbe.
  ## Le premier archivage du bras annuel l'a fait : il recopiait une variable
  ## de boucle perimee, et le panneau montrait UNE courbe sous une legende qui
  ## en annonce quinze. Le compte etait juste, et faux.
  ## Compter ne suffit donc pas : il faut de la dispersion.
  disp <- g[, .(s = stats::sd(y)), by = temperature_C][, max(s, na.rm = TRUE)]
  if (!is.finite(disp) || disp < 1e-8)
    stop(sprintf(paste0("bras %s : les %d courbes archivees sont IDENTIQUES ",
                        "(dispersion maximale %.2e sur le support). Une grille ",
                        "degeneree tracee en fond ferait lire un accord qui ",
                        "n'existe pas."), fr, n, disp))
  cat(sprintf("[lu] spec_curves_%s.csv : %d courbes, support [%g, %g], dispersion max %.3f\n",
              fr, n, min(g$temperature_C), max(g$temperature_C), disp))
  GRID[[fr]] <- g
}

##################################################################################
# 3. LA FIGURE
##################################################################################
plot_fans(CHANNEL, climate, grid_curves = GRID)

cat("=== ct6_fig_drf_paper.R termine ===\n")
