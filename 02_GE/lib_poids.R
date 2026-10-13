## lib_poids.R -- les poids de population, cote R (lot 51, 28/09/2026).
##
## temoin_part_de_la_cellule() : jumeau R de tools/temoins_poids.py -- la part doit
##   sommer a 1 PAR CELLULE, et le message NOMME cette normalisation.
##
## identite_aller_retour() : la garde qui aurait arrete le projet le premier jour.
##   Un taux national r_c reparti entre les cellules en proportion de la population
##   donne, pour la cellule z, v_z = somme_c pop(z,c) r_c / somme_c pop(z,c), la somme
##   portant sur les pays APPARIES (ceux qui ont un taux). Alors, EXACTEMENT :
##       somme_z  W_z v_z  =  somme_(z,c) pop(z,c) r_c,   W_z = somme_(c apparies) pop(z,c).
##   C'est une identite, pas une tolerance : on exige l'egalite a la precision machine
##   (ecart relatif <= 1e-9). Avec des parts qui ne sont pas des parts de population
##   (comptes de polygones, ou une normalisation qui vaut 1 partout), elle ne se
##   ferme pas, et la garde REFUSE en nommant le taux et l'ecart.

temoin_part_de_la_cellule <- function(dt, col = "part_de_la_cellule", tol = 1e-9) {
  s <- data.table::as.data.table(dt)[, .(s = sum(get(col))), by = .(x, y)]
  bad <- s[abs(s - 1) > tol]
  if (nrow(bad))
    stop(sprintf(paste0("[STOP] temoin PAR CELLULE (%s) : %d cellule(s) sur %d ne somment ",
                        "pas a 1 (min %.4g, max %.4g). Une part normalisee PAR PAYS fournie a ",
                        "la place d'une part PAR CELLULE produit exactement ceci."),
                 col, nrow(bad), nrow(s), min(s$s), max(s$s)))
  message(sprintf("[temoin PAR CELLULE] %s : %d cellules, chacune somme a 1", col, nrow(s)))
  temoin_part_non_triviale(data.table::as.data.table(dt)[[col]], col)
  invisible(nrow(s))
}

## temoin_part_non_triviale() (lot 51) : une part qui vaut 1 PARTOUT signe une cle de
## groupement qui contient la dimension sur laquelle la normalisation somme -- chaque
## groupe n'a qu'une ligne. Portee : ce defaut-la (les deux taux demographiques de 02a),
## PAS les comptes de polygones ni la population entiere comptee par pays, dont les
## parts ne valaient pas 1 (RECAP 3.238).
temoin_part_non_triviale <- function(x, nom, tol = 1e-12) {
  x <- x[is.finite(x)]
  if (length(x) && all(abs(x - 1) <= tol))
    stop(sprintf(paste0("[STOP] part non triviale (%s) : les %d valeurs valent 1. La cle ",
                        "de groupement contient la dimension sur laquelle la part se ",
                        "normalise ; chaque groupe n'a qu'une ligne."), nom, length(x)))
  message(sprintf("[temoin part non triviale] %s : %d valeurs, %d differentes de 1",
                  nom, length(x), sum(abs(x - 1) > tol)))
  invisible(TRUE)
}

## temoin_code_source() / temoin_imputation_evitable() (lot 60) : un taux national
## OBSERVE ne doit jamais etre remplace par une imputation parce qu'une jointure a
## echoue. Le premier refuse toute ligne de la source restee sans code pays (c'est
## la qu'une observation se perd) ou un code en double ; le second refuse toute
## imputation d'un code que la source contient.
temoin_code_source <- function(src, nom) {
  src <- data.table::as.data.table(src)
  perdus <- src[is.na(GID_0), unique(country)]
  if (length(perdus))
    stop(sprintf("[STOP] code source %s : %d pays de la source sans code ISO3 (%s). Leur valeur observee serait remplacee par une imputation.",
                 nom, length(perdus), paste(head(perdus, 10), collapse = " ; ")))
  dbl <- src[duplicated(GID_0), unique(GID_0)]
  if (length(dbl))
    stop(sprintf("[STOP] code source %s : codes ISO3 en double (%s).", nom, paste(dbl, collapse = " ")))
  message(sprintf("[temoin code source] %s : %d pays, tous codes, aucun double", nom, nrow(src)))
  invisible(TRUE)
}

temoin_imputation_evitable <- function(imputes, source, nom) {
  evitables <- intersect(imputes, source[!is.na(source)])
  if (length(evitables))
    stop(sprintf("[STOP] imputation evitable %s : %d pays imputes alors que la source les porte (%s).",
                 nom, length(evitables), paste(evitables, collapse = " ")))
  message(sprintf("[temoin imputation evitable] %s : %d imputes, aucun present dans la source", nom, length(imputes)))
  invisible(TRUE)
}

identite_aller_retour <- function(nom, poids, taux, valeurs, par = character(0), tol = 1e-9) {
  ## poids   : (id, GID_0, pop)          population de l'intersection cellule x pays
  ## taux    : (GID_0, r [, par])        le taux national reparti
  ## valeurs : (id, v [, par])           la valeur que 02a a attribuee a la cellule
  p <- data.table::as.data.table(poids)[, c("id", "GID_0", "pop"), with = FALSE]
  t <- data.table::as.data.table(taux)[, c("GID_0", "r", par), with = FALSE]
  v <- data.table::as.data.table(valeurs)[, c("id", "v", par), with = FALSE]
  ## Cle complete (lot 59). Un taux qui porte plusieurs annees passe SANS par = "year",
  ## ou des valeurs en double par (cellule, par) : on refuse ici, en le nommant.
  ## Les doublons du TAUX sur (pays, par) ne sont pas refuses : SSP2 IAMC porte deux
  ## modeles par pays, que 02a moyenne ; l'identite verifie alors la repartition des
  ## lignes fournies, et reste exacte.
  if (!"year" %in% par && "year" %in% names(taux) && data.table::uniqueN(taux$year) > 1)
    stop(sprintf(paste0("[STOP] identite aller-retour %s : cle incomplete, les taux portent ",
                        "%d annees et par ne contient pas \"year\"."),
                 nom, data.table::uniqueN(taux$year)))
  if (anyDuplicated(v, by = c("id", par)))
    stop(sprintf(paste0("[STOP] identite aller-retour %s : cle incomplete, les valeurs ont des ",
                        "doublons sur (%s)."), nom, paste(c("id", par), collapse = ", ")))
  t <- t[is.finite(r)]
  ## Plusieurs-a-plusieurs VOULU quand le taux porte par (les poids n'ont pas d'annee) :
  ## chaque (cellule, pays) se joint a chaque annee de son pays. allow.cartesian leve le
  ## plafond nrow(x)+nrow(i) de data.table, qui refusait des 3 annees (lot 59, 02a).
  a <- merge(p, t, by = "GID_0", allow.cartesian = TRUE)         # couples apparies
  droite <- a[, .(droite = sum(pop * r)), by = par]
  w <- a[, .(W = sum(pop)), by = c("id", par)]
  gauche <- merge(w, v, by = c("id", par))[, .(gauche = sum(W * v)), by = par]
  k <- if (length(par)) merge(gauche, droite, by = par) else cbind(gauche, droite)
  k[, ecart := abs(gauche - droite) / pmax(abs(droite), 1e-300)]
  message(sprintf("[identite aller-retour] %s : ecart relatif max %.2e (tolerance %.0e) sur %d %s",
                  nom, max(k$ecart), tol, nrow(k), if (length(par)) par else "agregat"))
  if (any(!is.finite(k$ecart)) || max(k$ecart) > tol)
    stop(sprintf(paste0("[STOP] identite aller-retour %s : somme_z W_z v_z = %.10g contre ",
                        "somme pop x taux = %.10g (ecart relatif %.2e). La valeur des cellules ",
                        "n'est pas la repartition PROPORTIONNELLE A LA POPULATION du taux national."),
                 nom, k$gauche[which.max(k$ecart)], k$droite[which.max(k$ecart)], max(k$ecart)))
  invisible(k)
}
