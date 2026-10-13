## temoin_population_pays.R -- la population implicite des pays confrontee au MONDE (lot 45).
##
## CE QU'IL FERME. Les poids de population de la chaine (lengths_country.csv,
## population_weight_all.csv) comptaient la population ENTIERE d'une cellule pour
## chaque pays ou polygone qui la touche. Tous les controles existants verifiaient
## une normalisation (« les parts somment a 1 ») que ce defaut PRESERVE. Ce temoin
## confronte le calcul a des grandeurs du monde, pas a un autre fichier du depot :
##
##   1. MONDE : somme des populations implicites / population GHSL de la meme
##      periode. Tolerance 1 % : une repartition correcte PARTITIONNE la population
##      GHSL entre pays disjoints ; elle ne peut perdre que les pixels hors des
##      polygones (franges de la simplification), et ne peut rien compter deux fois.
##   2. PAYS : ecart relatif a la population Banque mondiale 2015 (load_wb, meme
##      fichier que 01b). Tolerance : moyenne PONDEREE PAR LA POPULATION des |ecarts|
##      <= 5 % (GHSL R2023A est cale sur les totaux ONU, qui suivent la Banque
##      mondiale a quelques pourcents).
##   Il NOMME TOUJOURS les cinq pays aux plus grands ecarts : c'est la que se
##   logeront les defauts suivants.
##   Tolerances posees le 27/09/2026 AVANT la construction de la matrice corrigee.
##
## Usage : Rscript tools/temoin_population_pays.R            (test connu-faux)
##         source() puis temoin_population_pays(tab, ref_monde, wb)

TOL_MONDE <- 0.01
TOL_PAYS  <- 0.05
## POPULATION HORS DE TOUT POLYGONE (lot 46, addendum). La grandeur globale qui dit
## l'etat de la couverture : pixels GHSL dont le centre ne tombe dans aucun pays.
## Seuil 1 % : a 30'' (~1 km) sur des polygones NON simplifies, seuls les pixels
## cotiers dont le centre est en mer restent dehors ; le fichier simplifie a 0,01°
## en perdait 2,5 % (184 M), le non simplifie 0,4 % (27/09/2026).
TOL_HORS <- 0.01
## LES EXCEPTIONS, NOMMEES ET NON EXEMPTEES EN SILENCE. Tout pays dont l'ecart a la
## Banque mondiale depasse SEUIL_ISOLE doit figurer ici avec sa cause ; un nouveau venu
## est un defaut, pas la geographie : la liste ne doit pas s'allonger.
SEUIL_ISOLE <- 0.25
## VIDE depuis le lot 51 (28/09/2026) : la regle de rattachement a un pas de trame a
## ramene les six anciennes exceptions (MHL, KIR, TUV, MDV, NRU, FSM) sous le seuil.
## Une regle qui rend ses propres exceptions inutiles est le signe qu'elle etait la
## bonne. MINIMALITE : une exception dont l'ecart est revenu sous le seuil est un
## defaut (elle exempterait en silence une regression future) : le temoin refuse.
EXCEPTIONS <- character(0)
## TABLE DE CORRESPONDANCE GADM -> Banque mondiale, pour les codes qui ne s'alignent pas
## un pour un. On ne touche ni aux donnees ni a la definition des pays : on compare la
## somme des codes GADM a l'agregat que publie la Banque mondiale. Chaque entree porte
## sa raison ; une entree qui ne change plus rien est un defaut (table minimale).
CORRESPONDANCE <- list(
  CYP = list(gadm = c("CYP", "ZNC", "XAD"),
             raison = paste("GADM 4.1 decoupe l'ile en CYP, ZNC (nord) et XAD (bases",
                            "souveraines) ; la Banque mondiale publie un seul CYP pour l'ile.")))

temoin_population_pays <- function(tab, ref_monde, wb, etiquette = "", hors = NULL) {
  ## tab  : data.table(GID_0, pop) -- population implicite par pays
  ## wb   : data.table(GID_0, population) -- Banque mondiale 2015
  ## hors : population des pixels hors de tout polygone (NULL = non mesuree)
  tot <- sum(tab$pop, na.rm = TRUE)
  ecart_monde <- tot / ref_monde - 1
  message("  table de correspondance GADM -> Banque mondiale :")
  tab <- copy(tab)
  for (cle in names(CORRESPONDANCE)) {
    e <- CORRESPONDANCE[[cle]]
    autres <- setdiff(e$gadm, cle)
    apport <- tab[GID_0 %in% autres, sum(pop)]
    message(sprintf("    %s <- %s (+%.0f) : %s", cle, paste(e$gadm, collapse = " + "),
                    apport, e$raison))
    if (!(apport > 0))
      stop("[STOP] table de correspondance : l'entree ", cle, " ne change plus rien ",
           "(aucune population sous ", paste(autres, collapse = ", "), ") -- a retirer.")
    tab[GID_0 == cle, pop := pop + apport]
    tab <- tab[!GID_0 %in% autres]
  }
  if (!is.null(hors)) {
    message(sprintf("  population hors de tout polygone : %.4g = %.2f %% du monde (seuil %.0f %%)",
                    hors, 100 * hors / ref_monde, 100 * TOL_HORS))
    if (hors / ref_monde > TOL_HORS)
      stop("[STOP] population hors de tout polygone : ", sprintf("%.2f %%", 100 * hors / ref_monde),
           " du monde, au-dela de ", 100 * TOL_HORS, " %.")
  }
  m <- merge(tab, wb, by = "GID_0")
  m <- m[is.finite(population) & population > 0 & is.finite(pop)]
  m[, ecart := pop / population - 1]
  pond <- sum(abs(m$ecart) * m$population) / sum(m$population)
  top <- m[order(-abs(ecart))][1:min(5L, .N)]
  msg <- sprintf(paste0(
    "[temoin population %s] MONDE : %.4g implicites pour %.4g (GHSL), ecart %+.1f %% ",
    "(tolerance %.0f %%) | PAYS (%d apparies a la Banque mondiale 2015) : ecart ",
    "absolu moyen pondere %.1f %% (tolerance %.0f %%) | cinq plus grands ecarts : %s"),
    etiquette, tot, ref_monde, 100 * ecart_monde, 100 * TOL_MONDE, nrow(m),
    100 * pond, 100 * TOL_PAYS,
    paste(sprintf("%s %+.0f %%", top$GID_0, 100 * top$ecart), collapse = ", "))
  message(msg)
  if (abs(ecart_monde) > TOL_MONDE)
    stop("[STOP] ", msg, "\n  => la somme des pays n'est pas la population du monde : ",
         "une cellule est comptee dans plusieurs pays (ou perdue).")
  if (pond > TOL_PAYS)
    stop("[STOP] ", msg, "\n  => les populations des pays s'ecartent de la Banque mondiale.")
  isoles <- m[abs(ecart) > SEUIL_ISOLE]
  message(sprintf("  ecarts > %.0f %% : %s", 100 * SEUIL_ISOLE,
                  paste(sprintf("%s %+.0f %%", isoles$GID_0, 100 * isoles$ecart), collapse = ", ")))
  for (c in names(EXCEPTIONS))
    message(sprintf("    exception %s : %s", c, EXCEPTIONS[[c]]))
  inutiles <- setdiff(names(EXCEPTIONS), isoles$GID_0)
  if (length(inutiles))
    stop("[STOP] ", msg, "\n  => exception(s) devenue(s) inutile(s) (ecart <= ",
         100 * SEUIL_ISOLE, " %) : ", paste(inutiles, collapse = ", "),
         " -- a retirer de EXCEPTIONS (liste minimale).")
  nouveaux <- setdiff(isoles$GID_0, names(EXCEPTIONS))
  if (length(nouveaux))
    stop("[STOP] ", msg, "\n  => ecart > ", 100 * SEUIL_ISOLE, " % hors des exceptions nommees : ",
         paste(nouveaux, collapse = ", "), ". La liste s'allonge : c'est un defaut, pas la geographie.")
  invisible(list(ecart_monde = ecart_monde, ecart_pays = pond, top = top))
}

## ---- mode script : test connu-faux sur lengths_country.csv d'aujourd'hui -----
if (sys.nframe() == 0L) {
  suppressMessages({ library(data.table) })
  ## Racine derivee de l'emplacement du script, jamais du repertoire courant (lot 51).
  .self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  .rt <- if (nzchar(Sys.getenv("JMP_ROOT"))) Sys.getenv("JMP_ROOT") else
    dirname(dirname(normalizePath(.self, winslash = "/")))
  suppressMessages(source(file.path(.rt, "01_estimation", "00_parameters.R")))
  source(file.path(.rt, "01_estimation", "lib_drf.R"))
  wb <- load_wb(file.path(DIR_IN, "WB_pop_country/population_data.csv"), "population",
                period = "annual", year_min = 2015L, year_max = 2015L, verbose = FALSE)
  wb <- wb[, .(GID_0, population)]
  cell <- fread(file.path(.rt, "tools", "connus_faux", "population_all_2015.csv"))  # GHSL 2015 au degre, fige
  cell[, `:=`(x = round(x, 1), y = round(y, 1))]
  lc <- fread(file.path(.rt, "tools", "connus_faux", "lengths_country_FAUTIF_2015.csv"))  # fige (lot 51)
  lc[, `:=`(x = round(x, 1), y = round(y, 1))]
  j <- merge(lc, cell[value > 0], by = c("x", "y"))
  implic <- j[, .(pop = median(value / share_population)), by = GID_0]
  r <- tryCatch({ temoin_population_pays(implic, sum(cell$value, na.rm = TRUE), wb,
                                         "lengths_country.csv 2015"); "PASSE" },
                error = function(e) conditionMessage(e))
  ok <- !identical(r, "PASSE") && grepl("+23", r, fixed = TRUE)
  cat(sprintf("\n%s  connu-faux : lengths_country.csv doit etre REFUSE en nommant +23 %% -> %s\n",
              if (ok) "OK   " else "ECHEC", if (identical(r, "PASSE")) "passe (!)" else "refuse"))
  quit(status = if (ok) 0L else 1L)
}
