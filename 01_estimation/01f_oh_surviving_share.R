##################################################################################
# 01f_oh_surviving_share.R — la part du differentiel identifiant qui SURVIT au
# demeaning exportateur x annee.
#
# SORTIE   output/oh_surviving_share.csv  (une ligne par terme + la moyenne)
# LIT      intermediate_outputs/BACI_raw_pairs_annual.rds   (panneau d'analyse)
#          intermediate_outputs/climate_long_annual.rds     (la distribution)
# N'ECRIT  rien d'autre. AUCUNE ESTIMATION : que des projections OLS sur des
#          effets fixes, c'est-a-dire des demeanings.
#
# ---------------------------------------------------------------------------
# LA DEFINITION, ET POURQUOI ELLE DOIT ETRE ECRITE
# ---------------------------------------------------------------------------
# Le manuscrit ecrit « X\% of the identifying differential survives the
# demeaning ». Trois objets pouvaient porter ce nom -- une colonne de base, le
# vecteur des termes climatiques, ou l'indice ajuste -- et ils ne donnent pas le
# meme nombre. L'arbitrage du 16/09/2026 retient :
#
#     part survivante = moyenne sur k de (1 - R2_k),
#     R2_k = R2 de la projection de Delta_k sur les EF exportateur x annee,
#
# avec Delta_k = g_k(T(s)) - g_k(T(r)), les memes colonnes que tab_fe_absorption
# et la meme fonction r2_on. C'est donc la MEME construction, avec un autre jeu
# d'effets fixes -- ce que le brief demande.
#
# ---------------------------------------------------------------------------
# LE CONTROLE QUI REND LE CHIFFRE CROYABLE
# ---------------------------------------------------------------------------
# Ce script reconstruit la matrice de design hors de 01a. Une reconstruction qui
# differe d'un filtre ou d'un merge donnerait un nombre plausible et faux. On
# recalcule donc AUSSI les R2 paire x produit, que 01a ecrit au meme passage
# dans fe_absorption.csv (r2_pair_raw), et on refuse d'ecrire si on ne les
# retrouve pas. Le chiffre neuf n'est publiable que si le chiffre connu revient.
# Lot 52 : ces R2 etaient ecrits en dur, photographie de l'archive ; ils se
# lisent desormais (garde_ancre_absorption(), lib_drf.R).
##################################################################################

t0 <- proc.time(); cat("=== 01f_oh_surviving_share.R ===\n")

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Racine indeterminable. Poser JMP_ROOT.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
suppressPackageStartupMessages({ library(data.table); library(lfe) })

OUT <- file.path(.root, "output")
INT <- file.path(.root, "intermediate_outputs")

## L'ancre se LIT (lot 52) : les R2 paire x produit que 01a vient d'ecrire au
## meme passage, a pleine precision. Voir garde_ancre_absorption() (lib_drf.R).
FE_ABS <- file.path(OUT, "fe_absorption.csv")

cat("  lecture du panneau d'analyse...\n")
BACI <- readRDS(file.path(INT, "BACI_raw_pairs_annual.rds")); setDT(BACI)
climate <- readRDS(file.path(INT, "climate_long_annual.rds")); setDT(climate)
cat(sprintf("  panneau : %s lignes\n", format(nrow(BACI), big.mark = " ")))

cat(sprintf("  base de reference : %s degre %d, winsorisation [%g, %g]\n",
            REF_SPEC, REF_DEG, REF_WLOWER, REF_WUPPER))
## `scale` EST OBLIGATOIRE DEPUIS LE 20/09/2026, ET CE SITE LA RECLAMAIT.
## Ce script etait le seul appel de prepare_diffs() qui vivait du defaut
## `scale = 1` -- quatorzieme site, alors que le releve du 19/09 en comptait
## treize et concluait « zero site sans echelle ».
## AUCUN CHIFFRE PUBLIE N'EN DESCEND : ce script ne rend que des R2 de
## projection sur effets fixes, et le R2 de Delta_k sur des EF ne change pas
## quand Delta_k est divise par une constante. L'ancre bloquante ci-dessous
## (fe_absorption.csv du meme passage) le verifie a chaque passage.
## On passe SCALE_PROD quand meme : le canal est celui de la productivite, et
## un site qui ne declare pas son echelle est un site qui la perdra un jour.
d <- prepare_diffs(BACI, climate, REF_SPEC, REF_DEG, REF_WLOWER, REF_WUPPER,
                   scale = SCALE_PROD)
setDT(d)
cat(sprintf("  matrice de design : %s lignes\n", format(nrow(d), big.mark = " ")))

## Les deux jeux d'effets fixes. `fe_prod_year` existe deja dans le design ;
## l'exportateur x annee est construit ici, et il ne l'est nulle part ailleurs.
d[, fe_orig_dest_prod := paste(orig, dest, product_code_ISICRev3, sep = "_")]
d[, fe_orig_year      := paste(orig, year, sep = "_")]

r2_on <- function(col, fe) {
  y <- as.numeric(d[[col]])
  r <- as.numeric(residuals(felm(as.formula(paste0(col, " ~ 1 | ", fe)),
                                 data = as.data.frame(d))))
  1 - var(r) / var(y, na.rm = TRUE)
}

res <- rbindlist(lapply(seq_len(REF_DEG), function(k) {
  col <- sprintf("diff_%s_%d", REF_SPEC, k)
  if (!col %in% names(d)) stop("[STOP] colonne absente du design : ", col)
  a <- r2_on(col, "fe_orig_dest_prod")
  b <- r2_on(col, "fe_orig_year")
  cat(sprintf("  terme %d : R2 paire x produit = %.4f | ", k, a))
  cat(sprintf("R2 exportateur x annee = %.4f -> survit %.1f %%\n", b, 100 * (1 - b)))
  data.table(term = k, r2_pair_prod = a, r2_exporter_year = b,
             surviving_share_pct = 100 * (1 - b))
}))

## ---- LE CONTROLE, BLOQUANT --------------------------------------------------
res[, anchor_pair_prod := garde_ancre_absorption(r2_pair_prod, FE_ABS)]

part <- mean(res$surviving_share_pct)
cat(sprintf("\n  >>> part survivante (moyenne sur les %d termes) = %.1f %%\n",
            REF_DEG, part))

res <- rbind(res, data.table(
  term = NA_integer_, r2_pair_prod = NA_real_,
  r2_exporter_year = mean(res$r2_exporter_year),
  surviving_share_pct = part, anchor_pair_prod = NA_real_), fill = TRUE)
res[, definition := paste("mean over terms of 1 - R2 of the projection of the",
                          "bilateral climate differential on exporter x year FE")]
fwrite(res, file.path(OUT, "oh_surviving_share.csv"))
cat(sprintf("  wrote oh_surviving_share.csv\n"))

cat(sprintf("\nTemps total: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))
cat("=== Done ===\n")
