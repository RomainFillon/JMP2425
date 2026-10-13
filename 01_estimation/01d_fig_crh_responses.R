##################################################################################
# 01d_fig_crh_responses.R — les quatre reponses re-estimees dans le pipeline de
# Cruz et Rossi-Hansberg : deux canaux x deux metriques d'exposition.
#
# Reads:  intermediate_outputs/crh/daily_{prod,amen}_bins.csv     <- 01d, C4(a)
#         intermediate_outputs/crh/daily_{prod,amen}_curve.csv    <- 01d, C4(a)
#         intermediate_outputs/crh/mean_{prod,amen}_curve.csv     <- 01d, C4(b)
#         intermediate_outputs/crh/mean_{prod,amen}_coeff.csv     <- 01d, C4(b)
#         intermediate_outputs/crh/crh_support_bounds.csv         <- 01d, C5
# Writes: output/fig_crh_responses.pdf
#         output/crh_daily_response_peaks.csv
#         output/crh_mean_response_optima.csv
# Lu par: JMP_V11.tex (fig:crh_response) ;
#         02_GE/emit_values_empirics.py (\EMPcrhPeakProd, \EMPcrhPeakAmen,
#         \EMPcrhMeanOpt, \EMPcrhMeanOptAmen, \EMPcrhDailyRefBin)
#
# ------------------------------------------------------------------------------
# POURQUOI CETTE FIGURE A SON PROPRE PRODUCTEUR (18/09/2026, lot 10)
# ------------------------------------------------------------------------------
# Elle remplace fig_crh_daily_response.pdf, qui etait dessinee dans la section C7
# de 01d_crh_triangulation.R. Elle n'a jamais eu besoin d'y etre : les quatre
# panneaux sont des LECTURES de CSV que 01d a deja ecrits, et pas un seul
# coefficient n'y est re-estime. Les laisser dans 01d obligeait a rejouer une
# heure de calcul -- inversion du modele spatial comprise, avec un .mat de
# 2,1 Go -- pour deplacer un trait.
# ⚠️ ET LE REJEU N'AURAIT PAS ETE NEUTRE : DIR_CRH ne suit PAS JMP_OUTPUT_DIR,
# de sorte qu'un rejeu de 01d reecrit intermediate_outputs/crh/ meme sous
# redirection. Sortir la figure evite ce rejeu, donc evite ce risque.
#
# ------------------------------------------------------------------------------
# CE QUE LES DEUX COLONNES MONTRENT
# ------------------------------------------------------------------------------
# Gauche  : la metrique JOURNALIERE (01d, C4(a)). delta(T) est le coefficient du
#           bin T sur la PART des jours de la cellule-periode qui y tombent. Son
#           unite est donc le point de log PAR UNITE DE PART DE JOURS, et non
#           « par jour » : les parts somment a un par cellule-periode, et le
#           nombre de jours va de 365 a 3653 selon la longueur du bloc. Elle se
#           lit RELATIVEMENT AU BIN DE REFERENCE, dont le coefficient est nul par
#           construction (il est retire du membre de droite, les parts etant
#           colineaires a la constante des effets fixes).
# Droite  : la metrique MOYENNE ANNUELLE (01d, C4(b)), c'est-a-dire le bras que
#           le manuscrit compare au journalier et que la projection C6 consomme.
#           La courbe est la pente estimee -- effet marginal de +1 degC de
#           moyenne annuelle de la distribution journaliere sur le log du
#           fondamental -- lissee par spline naturelle ponderee. Son unite est le
#           point de log PAR DEGRE.
#
# ⚠️ CE N'EST PAS LA COLONNE QUE CE SCRIPT TRACAIT LE 18/09 AU MATIN, et le
# changement est une CORRECTION, pas un choix de presentation. Elle tracait
# alors run_crh_baseline (01d, C2), la specification d'ORIGINE de Cruz et
# Rossi-Hansberg, dont le regresseur est `temp_10y_1h` : la moyenne mobile
# decennale de JANVIER au nord et de JUILLET au sud -- une temperature de mois
# d'hiver, que leur papier nomme « winter temperature » (p. 20) et leurs axes
# « January-July Temperature ». Mesure du 18/09 sur 68 192 cellules-periodes :
# mediane 0,5 degC contre 10,4 degC pour la moyenne annuelle decennale, ecart
# moyen -10,6 degC. L'etiqueter « decadal mean » et la mettre en regard d'une
# abscisse journaliere melangeait deux temperatures distantes d'une dizaine de
# degres. La colonne C2 reste l'objet des points de felicite
# (\EMPcrhOptProd, \EMPcrhOptAmen), qui sont cites dans la PROSE de la section
# et n'ont pas besoin d'un panneau.
#
# 🔑 L'ABSCISSE EST DESORMAIS PARTAGEE SANS ABUS. Les deux colonnes portent la
# meme grandeur : la colonne de droite indexe la MOYENNE de la distribution que
# la colonne de gauche decompose. C'est la propriete que 01d C4 met en avant --
# « the scalar is the annual mean OF THE DAILY DISTRIBUTION, hence in THE SAME
# UNIT as (a) ».
#
# ⚠️ LES ORDONNEES RESTENT SEPAREES, ET C'EST MESURE. Les unites different (par
# part de jours a gauche, par degre a droite) et les magnitudes avec elles :
# productivite, [-7,06 ; +4,56] contre [-0,057 ; +0,364] ; amenite,
# [-1,15 ; +0,33] contre [-0,002 ; +0,028]. Une ordonnee commune ecraserait la
# colonne de droite sur la ligne de zero.
#
# ⛔ ET CE N'EST PAS Lambda. mean_{tag}_curve.csv porte DEUX colonnes : `delta`,
# la pente, qui est tracee ici, et `Lambda`, son integrale trapezoidale, qui est
# une fonction de NIVEAU (0 -> 8,83 en productivite) et que la projection C6
# utilise. Tracer `Lambda` donnerait une courbe de niveau sous une legende qui
# annonce un effet marginal.
##################################################################################

suppressMessages(library(data.table))

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Impossible de deduire la racine. Poser JMP_ROOT.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))

DIR_CRH <- file.path(DIR_TMP, "crh")
OUT     <- DIR_OUT
CHANNELS <- c("prod", "amen")
CH_TITLE <- c(prod = "Productivity", amen = "Amenity")

cat("=== 01d_fig_crh_responses.R ===\n")

.intrants <- c(
  file.path(DIR_CRH, "crh_support_bounds.csv"),
  file.path(DIR_CRH, sprintf("daily_%s_bins.csv", CHANNELS)),
  file.path(DIR_CRH, sprintf("daily_%s_curve.csv", CHANNELS)),
  file.path(DIR_CRH, sprintf("mean_%s_curve.csv", CHANNELS)),
  file.path(DIR_CRH, sprintf("mean_%s_coeff.csv", CHANNELS)))
.abs <- .intrants[!file.exists(.intrants)]
if (length(.abs))
  stop("[STOP] intrant(s) absent(s) :\n  ", paste(.abs, collapse = "\n  "),
       "\n  Tous sont produits par 01_estimation/01d_crh_triangulation.R.")

sb <- fread(file.path(DIR_CRH, "crh_support_bounds.csv"))

## LE ZERO DE LA COURBE DE DROITE, RAFFINE.
## mean_{tag}_coeff.csv porte `Topt_mean`, mais c'est le POINT DE GRILLE qui
## precede le changement de signe, pas la racine : sur une grille de 400 points
## au pas de 0,11 degC l'ecart atteint 0,025 degC. La figure trace la courbe, la
## macro doit donc porter le zero DE LA COURBE. On interpole lineairement sur
## l'intervalle de changement de signe, ce qui est exact pour une fonction
## echantillonnee et reste dans le meme fichier que la courbe.
zero_of <- function(mc) {
  cr <- which(diff(sign(mc$delta)) != 0)
  if (length(cr) != 1L)
    stop("[STOP] la pente de la metrique annuelle change de signe ", length(cr),
         " fois : « l'optimum » ne serait pas defini, et un premier passage a\n",
         "  zero se lirait comme un optimum unique.")
  i <- cr[1]
  approx(mc$delta[i:(i + 1)], mc$mu[i:(i + 1)], xout = 0)$y
}

pdf(file.path(OUT, "fig_crh_responses.pdf"), width = 9, height = 7.4)
op <- par(mfrow = c(2, 2), mar = c(4.0, 4.4, 2.4, 0.8), mgp = c(2.5, 0.8, 0),
          cex.axis = 0.9)
peaks <- list(); optima <- list()

for (tag in CHANNELS) {
  b  <- fread(file.path(DIR_CRH, sprintf("daily_%s_bins.csv", tag)))
  cv <- fread(file.path(DIR_CRH, sprintf("daily_%s_curve.csv", tag)))
  mc <- fread(file.path(DIR_CRH, sprintf("mean_%s_curve.csv", tag)))
  pts <- b[is.finite(t_center) & !is_ref]
  ## LA MEME ABSCISSE POUR LES DEUX PANNEAUX DE LA LIGNE, et elle est legitime :
  ## la droite indexe la moyenne de la distribution que la gauche decompose.
  xr <- range(c(cv$t, pts$t_center, mc$mu), finite = TRUE)
  ref_T <- b[is_ref == TRUE, t_center]

  ## ---------- gauche : la metrique journaliere ------------------------------
  ## L'ORDONNEE EST CELLE DE LA COURBE, pas celle des points. Deux coefficients
  ## de bins extremes -- les bins les moins peuples -- valent plusieurs fois
  ## l'amplitude de la fonction ; les laisser fixer l'echelle reduirait la
  ## fonction a un trait horizontal.
  ins <- cv$t >= sb$k_min & cv$t <= sb$k_max
  plot(cv$t, cv$delta, type = "n", xlim = xr,
       xlab = "Daily mean temperature (°C)",
       ## LE TIRET DOIT ETRE DANS UNE CHAINE. En plotmath, day-share est lu
       ## comme une SOUSTRACTION et rend « day - share » ; verifie a l image.
       ylab = expression(hat(delta)(T) * " (log points per unit of day-share)"),
       main = sprintf("%s, daily distribution", CH_TITLE[[tag]]),
       font.main = 1, cex.main = 1.0)
  abline(h = 0, col = "grey60", lty = 2)
  points(pts$t_center, pts$delta, pch = 16, cex = 0.6, col = "grey30")
  lines(cv$t[ins],  cv$delta[ins],  lwd = 2)
  lines(cv$t[!ins], cv$delta[!ins], lwd = 2, lty = 3)
  peaks[[length(peaks) + 1L]] <-
    data.table(channel = tag, peak_T = cv$t[which.max(cv$delta)])

  ## ---------- droite : la metrique moyenne annuelle -------------------------
  ## LA COURBE S'ARRETE A SON SUPPORT. L'abscisse est celle de la ligne, plus
  ## large que [min(mu), max(mu)] ; prolonger le trait jusqu'au bord ferait
  ## passer une extrapolation pour une estimation.
  z <- zero_of(mc)
  optima[[length(optima) + 1L]] <-
    data.table(channel = tag, T_opt = z, grid_step = diff(mc$mu)[1],
               support_lo = min(mc$mu), support_hi = max(mc$mu),
               daily_ref_bin_T = ref_T)
  plot(mc$mu, mc$delta, type = "n", xlim = xr,
       xlab = "Annual mean of the daily distribution (°C)",
       ylab = "Marginal effect of +1°C (log points)",
       main = sprintf("%s, annual mean", CH_TITLE[[tag]]),
       font.main = 1, cex.main = 1.0)
  abline(h = 0, col = "grey60", lty = 2)
  lines(mc$mu, mc$delta, lwd = 2)
  abline(v = z, col = "#d73027", lwd = 1.2, lty = 2)
  ## L'etiquette se pose du cote ou il reste de la place, et le cote est CHOISI
  ## sur la position du zero dans l'abscisse, pas fixe.
  yr <- range(mc$delta, finite = TRUE)
  .side <- if (z < mean(xr)) 4 else 2
  text(z, yr[2] - 0.06 * diff(yr), sprintf(" %.1f°C ", z),
       pos = .side, col = "#d73027", cex = 0.85)
}
par(op)
dev.off()

pk <- rbindlist(peaks)
fwrite(pk, file.path(OUT, "crh_daily_response_peaks.csv"))
op_dt <- rbindlist(optima)
fwrite(op_dt, file.path(OUT, "crh_mean_response_optima.csv"))
cat(sprintf("  fig_crh_responses.pdf  (support journalier [%.1f, %.1f], tirete = extrapolation)\n",
            sb$k_min, sb$k_max))
for (i in seq_len(nrow(pk)))
  cat(sprintf("    %-5s pic journalier %.1f °C | zero de la pente annuelle %.3f °C (pas %.4f)\n",
              pk$channel[i], pk$peak_T[i], op_dt$T_opt[i], op_dt$grid_step[i]))
cat(sprintf("    bin de reference journalier : %.2f °C (delta = 0 par construction)\n",
            op_dt$daily_ref_bin_T[1]))
cat("  wrote crh_daily_response_peaks.csv, crh_mean_response_optima.csv\n")
cat("=== Done ===\n")
