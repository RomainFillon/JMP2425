##################################################################################
# 01d_crh_clamp_sensitivity.R — la sensibilite de la difference nordique au
# point de clamp de delta(T).
#
# Reads:  intermediate_outputs/CMIP6_tas.csv                  [SERVEUR]
#         intermediate_outputs/crh/daily_{prod,amen}_curve.csv   <- 01d, C4
#         intermediate_outputs/crh/mean_{prod,amen}_curve.csv    <- 01d, C4
#         intermediate_outputs/crh/mean_{prod,amen}_coeff.csv    <- 01d, C4
#         intermediate_outputs/crh/amen_prod.csv                 <- 01d, C1
#         intermediate_outputs/crh/tab_crh_northern_reversal_daily.csv  <- 01d, C6
# Writes: output/crh_clamp_sensitivity.csv
#         output/crh_north_daymass.csv
# Lu par: 02_GE/emit_values_empirics.py
#         (\EMPcrhClampMove, \EMPcrhClampThirtyMove,
#          \EMPcrhNorthDayMassAboveThirty)
#
# ------------------------------------------------------------------------------
# CE QUE CE SCRIPT FAIT, ET CE QU'IL NE FAIT PAS
# ------------------------------------------------------------------------------
# Il n'estime RIEN. delta(T) et Lambda(mu) sont lus tels que 01d les a ecrits ;
# ce script ne fait que les EVALUER sur les distributions CMIP6 archivees, une
# fois telles quelles et une fois tenues plates au-dela d'un seuil. C'est
# exactement le calcul de project_channel() (01d, section C6) avec son argument
# `reliable_hi`, sorti de la boucle pour pouvoir le refaire a plusieurs seuils
# sans relancer l'inversion du modele spatial, qui coute une heure et charge
# 2,1 Go de .mat.
#
# ------------------------------------------------------------------------------
# LE CONTROLE QUI AUTORISE CETTE SORTIE DE BOUCLE
# ------------------------------------------------------------------------------
# « Ce script refait le meme calcul que 01d » est une affirmation. On la PROUVE
# a chaque passage : le bras non clampe est recalcule ici et compare a
# tab_crh_northern_reversal_daily.csv, l'artefact que 01d a ecrit. La tolerance
# est 1e-12, c'est-a-dire l'egalite a la precision machine et non « proche ».
# Mesure du 18/09/2026 : ecart 3,3e-16 sur le bras moyen, 4,6e-16 sur le bras
# journalier, dans les deux canaux. Si 01d bouge et pas ce script, la sortie ne
# se produit pas.
#
# ------------------------------------------------------------------------------
# LE CLAMP N'EST PAS DANS LA PROJECTION PUBLIEE
# ------------------------------------------------------------------------------
# project_channel() a `reliable_hi = Inf` par defaut : la projection du
# manuscrit utilise la fonction TELLE QU'ESTIMEE, sans clamp. Les deux seuils
# de ce script sont donc deux SENSIBILITES, et le « deplacement » qu'ils
# rapportent se lit contre ce bras non clampe.
#   CRH_CLAMP_SUPPORT (40 degC) = la borne du support bien peuple, celle que le
#     manuscrit nomme \EMPcrhSupportHi ;
#   CRH_CLAMP_SEVERE  (30 degC) = la borne plus severe, demandee pour savoir ce
#     que couterait un support declare beaucoup plus prudent.
#
# ⚠️ LE SIGNE EST L'INFORMATION, et les deux macros ne le portent pas de la
# meme facon, deliberement :
#   \EMPcrhClampMove       est une MAGNITUDE. La phrase du supplement dit
#     « clamping moves the Northern aggregate BY x log points » : un signe y
#     serait grammaticalement faux.
#   \EMPcrhClampThirtyMove est SIGNE. Elle porte le deplacement de la DIFFERENCE
#     nordique, et tout son contenu est que ce deplacement est POSITIF : couper
#     la queue chaude ne rapproche pas les deux metriques, il les ecarte. Une
#     magnitude laisserait croire que la severite du support joue contre le
#     resultat, alors qu'elle joue pour lui.
##################################################################################

suppressMessages(library(data.table))

## Racine : posee par .Rprofile, par l'orchestrateur, ou deduite du chemin du
## script. Meme amorce que ct4, pour la meme raison -- aucun chemin absolu.
.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .a <- commandArgs(trailingOnly = FALSE)
  .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
  if (!length(.f)) stop("Impossible de deduire la racine. Poser JMP_ROOT.")
  .root <- dirname(dirname(normalizePath(.f[1], winslash = "/")))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))

DIR_CRH <- file.path(DIR_TMP, "crh")
OUT     <- file.path(.root, "output")
CHANNELS <- c("prod", "amen")
CLAMPS   <- c(CRH_CLAMP_SUPPORT, CRH_CLAMP_SEVERE)

cat("=== 01d_crh_clamp_sensitivity.R ===\n")

.intrants <- c(file.path(DIR_TMP, "CMIP6_tas.csv"),
               file.path(DIR_CRH, "amen_prod.csv"),
               file.path(DIR_CRH, "tab_crh_northern_reversal_daily.csv"),
               file.path(DIR_CRH, sprintf("daily_%s_curve.csv", CHANNELS)),
               file.path(DIR_CRH, sprintf("mean_%s_curve.csv", CHANNELS)),
               file.path(DIR_CRH, sprintf("mean_%s_coeff.csv", CHANNELS)))
.abs <- .intrants[!file.exists(.intrants)]
if (length(.abs))
  stop("[STOP] intrant(s) absent(s) :\n  ", paste(.abs, collapse = "\n  "),
       "\n  Les fichiers de intermediate_outputs/crh/ sont produits par\n",
       "  01_estimation/01d_crh_triangulation.R ; CMIP6_tas.csv vient du serveur.")

## ---- les distributions CMIP6, en PARTS par cellule et par annee --------------
cmip <- fread(file.path(DIR_TMP, "CMIP6_tas.csv"))[
  year %in% c(CRH_BASE_YEAR, CRH_FUT_YEAR)]
cmip[, c := count / sum(count), by = .(x, y, year)]

## Poids de population, ceux de l'annexe CRH : la coupe 2005 de amen_prod.
pop05 <- fread(file.path(DIR_CRH, "amen_prod.csv"))[
  year == 2005, .(lat, lon, pop = as.numeric(pop_x))]

## ---- les deux bras, cellule par cellule --------------------------------------
## Portage litteral de project_channel() (01d, C6). La seule difference est que
## `hi` y est un argument de boucle et non un defaut.
daily_arm <- function(tag, hi = Inf) {
  cv  <- fread(file.path(DIR_CRH, sprintf("daily_%s_curve.csv", tag)))
  cvt <- cv[t <= hi]
  d <- copy(cmip)[, dlt := approx(cvt$t, cvt$delta, xout = value, rule = 2)$y]
  gw <- dcast(d[, .(g = sum(dlt * c)), by = .(x, y, year)],
              x + y ~ year, value.var = "g")
  setnames(gw, as.character(c(CRH_BASE_YEAR, CRH_FUT_YEAR)), c("g_base", "g_fut"))
  gw[, .(x, y, value = g_fut - g_base)]
}

mean_arm <- function(tag) {
  mc <- fread(file.path(DIR_CRH, sprintf("mean_%s_curve.csv", tag)))
  co <- fread(file.path(DIR_CRH, sprintf("mean_%s_coeff.csv", tag)))
  smin <- co[param == "mu_support_min", value]
  smax <- co[param == "mu_support_max", value]
  Lam <- function(m) approx(mc$mu, mc$Lambda,
                            xout = pmin(pmax(m, smin), smax), rule = 2)$y
  mu <- cmip[, .(mu = sum(value * c)), by = .(x, y, year)]
  mw <- dcast(mu, x + y ~ year, value.var = "mu")
  setnames(mw, as.character(c(CRH_BASE_YEAR, CRH_FUT_YEAR)), c("mu_base", "mu_fut"))
  mw <- mw[is.finite(mu_base) & is.finite(mu_fut)]
  mw[, .(x, y, value = Lam(mu_fut) - Lam(mu_base))]
}

## Agregat nordique, pondere population. Meme masque que aggregate_north() :
## jointure INTERNE avec les poids, puis y >= CRH_NORTH_LAT.
north_agg <- function(tab) {
  r <- merge(tab, pop05, by.x = c("x", "y"), by.y = c("lon", "lat"))
  r <- r[y >= CRH_NORTH_LAT]
  list(n = nrow(r), v = sum(r$value * r$pop) / sum(r$pop))
}

## ---- le controle de non-regression contre 01d --------------------------------
ref <- fread(file.path(DIR_CRH, "tab_crh_northern_reversal_daily.csv"))
TOL <- 1e-12
rows <- list()
for (tag in CHANNELS) {
  m0 <- north_agg(mean_arm(tag))
  d0 <- north_agg(daily_arm(tag))
  r  <- ref[channel == tag & iso3 == "NORTH>=45N"]
  if (nrow(r) != 1L)
    stop("[STOP] canal ", tag, " : ligne NORTH>=45N introuvable dans\n",
         "  tab_crh_northern_reversal_daily.csv.")
  .em <- abs(m0$v - r$mean); .ed <- abs(d0$v - r$daily)
  if (.em > TOL || .ed > TOL)
    stop(sprintf(paste0("[STOP] canal %s : ce script ne reproduit pas la\n",
                        "  projection de 01d. ecart mean %.3e | daily %.3e\n",
                        "  (tolerance %.0e). L'evaluation sur archives n'est\n",
                        "  legitime que tant qu'elle redonne le bras publie."),
                 tag, .em, .ed, TOL))
  if (d0$n != r$n_cells)
    stop(sprintf("[STOP] canal %s : %d cellules contre %d au manifeste de 01d.",
                 tag, d0$n, r$n_cells))
  cat(sprintf("  [non-regression %-4s] mean %+.8f | daily %+.8f | %d cellules",
              tag, m0$v, d0$v, d0$n),
      sprintf(" | ecart max %.1e\n", max(.em, .ed)))

  rows[[length(rows) + 1L]] <- data.frame(
    channel = tag, clamp_hi = Inf, n_cells = d0$n,
    mean_north = m0$v, daily_north = d0$v, diff_north = m0$v - d0$v,
    move_daily = 0, move_diff = 0, stringsAsFactors = FALSE)

  for (hi in CLAMPS) {
    dh <- north_agg(daily_arm(tag, hi))
    rows[[length(rows) + 1L]] <- data.frame(
      channel = tag, clamp_hi = hi, n_cells = dh$n,
      mean_north = m0$v, daily_north = dh$v, diff_north = m0$v - dh$v,
      move_daily = dh$v - d0$v,
      move_diff  = (m0$v - dh$v) - (m0$v - d0$v), stringsAsFactors = FALSE)
    cat(sprintf(paste0("  [clamp %2.0f degC   %-4s] daily %+.8f | deplacement ",
                       "daily %+.4e | de la difference %+.4e\n"),
                hi, tag, dh$v, dh$v - d0$v, (m0$v - dh$v) - (m0$v - d0$v)))
  }
}
sens <- do.call(rbind, rows)
write.csv(sens, file.path(OUT, "crh_clamp_sensitivity.csv"), row.names = FALSE)
cat("  wrote crh_clamp_sensitivity.csv\n")

## ---- la masse-jour nordique au-dessus de chaque seuil -------------------------
## C'est ce qui rend le deplacement interpretable : un deplacement faible parce
## que la masse est faible n'est pas le meme enonce qu'un deplacement faible
## malgre une masse importante. La colonne d'abord, le verdict ensuite.
cf <- merge(cmip[year == CRH_FUT_YEAR], pop05,
            by.x = c("x", "y"), by.y = c("lon", "lat"))
nn <- cf[y >= CRH_NORTH_LAT]
dm <- do.call(rbind, lapply(CLAMPS, function(hi) {
  per_cell <- nn[, .(above = sum(c[value > hi])), by = .(x, y, pop)]
  data.frame(threshold_C = hi, year = CRH_FUT_YEAR,
             n_cells = nrow(per_cell),
             share_pct = 100 * sum(per_cell$above * per_cell$pop) /
                         sum(per_cell$pop),
             stringsAsFactors = FALSE)
}))
write.csv(dm, file.path(OUT, "crh_north_daymass.csv"), row.names = FALSE)
for (i in seq_len(nrow(dm)))
  cat(sprintf("  [masse-jour >=%.0fN] au-dessus de %2.0f degC en %d : %.4f %% (%d cellules)\n",
              CRH_NORTH_LAT, dm$threshold_C[i], dm$year[i], dm$share_pct[i],
              dm$n_cells[i]))
cat("  wrote crh_north_daymass.csv\n")
cat("=== Done ===\n")
