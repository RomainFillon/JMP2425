##################################################################################
# emit_tab_contribA_grid.R — table de la contribution A : les deux bras et leur
# ecart, une ligne par specification.
#
# SORTIES
#   output/contribA_grid.csv        CSV producteur (source machine de la table)
#   output/contribA_grid_counts.csv les comptes, pour emit_values_empirics.py
#   output/tab_contrast_grid.tex    l'exhibit  (⚠️ PAS tab_contribA_grid.tex : ce
#                                   nom-la designe une table RETIREE le 11/09/2026,
#                                   et l'en-tete l'a porte jusqu'au 18/09)
#   output/fig_contrast_grid.pdf    la figure, horodatage neutralise
#
# ENTREES -- PRODUITES PAR LE PAQUET depuis la migration du 2026-09-04
#   output/contribA/test1_delta_north.csv   <- ct1_delta_north_grid.R
#   output/contribA/tache1_diff_bras.csv    <- ct2_between_arm_test.R
#
# ⚠ LA LARGEUR DE LA TABLE A ETE MESUREE SUR JMP_V9.tex, PAS SUR JMP_V10.tex.
# Les 398,2 pt retenus ci-dessous valent pour la geometrie de V9 (letterpaper,
# left/right 1.2in). A re-mesurer sur V10 avant de conclure que la table ne
# deborde pas : une mesure faite dans le mauvais cadre a deja donne un faux
# vert de 29 pt sur ce meme exhibit.
#
# ---------------------------------------------------------------------------
# CE QUE LA TABLE PORTE, ET POURQUOI
# ---------------------------------------------------------------------------
# Trois colonnes de POINTS (annuel, journalier, ecart) + l'INTERVALLE DE L'ECART.
# ⚠️ Les intervalles PAR BRAS sont volontairement absents : ils alourdiraient sans
#    porter d'argument. L'intervalle de l'ecart, lui, est indispensable -- sans
#    lui la colonne d'ecart ne montre pas ce qu'elle affirme, et le lecteur
#    devrait croire la note sur parole pour le 15/15 (mesure du 18/09/2026).
# 🔑 C'est aussi ce qui justifie l'existence de cet emetteur : les POINTS sont
#    deja tous dans tache1_diff_bras.csv. Sans les intervalles, la jointure
#    n'aurait rien ajoute.
#
# ---------------------------------------------------------------------------
# LARGEUR — MESUREE, PAS ESTIMEE (27/08)
# ---------------------------------------------------------------------------
# Justification du manuscrit : 440,85 pt = 15,49 cm (JMP_V9.tex l. 10 :
# letterpaper, left/right 1.2in). ⚠️ Une premiere mesure faite en margin=1in
# donnait 469,8 pt : mesurer dans le mauvais cadre aurait donne un faux vert de
# 29 pt. Largeurs relevees par \settowidth dans la classe et la police reelles :
#   libelles longs + en-tetes longs   467,3 pt   DEBORDE de 26,5 pt
#   libelles courts (« NS7 [2.5%] »)  374,6 pt   marge +66,3
#   libelles LONGS + en-tetes COURTS  398,2 pt   marge +42,6  <- RETENU
# ⇒ On raccourcit les EN-TETES, pas les libelles : « natural spline, 7 df, 2.5%
#   winsor. » se lit sans legende, « NS7 [2.5%] » ne se lit pas sans elle.
# ⚠️ \dimexpr n'est PAS une reserve de largeur : il impose une largeur a UNE
#    colonne (ici la note en multicolumn). Il ne dispense pas de mesurer les
#    cinq colonnes de donnees.
#
# ---------------------------------------------------------------------------
# LA NOTE — sens MESURE, et il est l'inverse de l'intuition
# ---------------------------------------------------------------------------
# Mesure du 2026-09-18, grille a QUINZE specifications (GRID_SPECS) :
# rho entre bras ∈ [0,49 ; 0,74], TOUTES POSITIVES. Or var(A-B) = var(A)+var(B)
# -2cov : une covariance positive REDUIT la variance de la difference.
#   se_diff / se_diff_naive ∈ [54 % ; 77 %] sur les quinze specifications.
# ⇒ Traiter les bras comme independants est CONSERVATEUR, pas permissif.
# ⛔ Ecrire « le calcul naif surestime la precision » inverserait l'argument :
#    cela ferait croire que la covariance SAUVE le 15/15, donc qu'il en depend.
# ⛔ Et « le signe est determine sous l'un comme sous l'autre calcul » est FAUX :
#    15/15 avec covariance, 12/15 sans -- les trois exceptions sont les
#    polynomes de degre 2.
#    ⚠️ AUCUN DE CES NOMBRES N'EST ECRIT DANS LA SORTIE : la note de la table
#    les DERIVE, et la liste des exceptions avec eux. Ce commentaire les repete
#    pour le lecteur du code, et il n'a donc aucune autorite sur l'exhibit.
#    🚨 C'EST EXACTEMENT POUR CELA QU'IL EST DERIVE. Ce commentaire a porte
#    « 21/24 ... les trois exceptions sont les polynomes de degre 2 » jusqu'au
#    2026-09-16, et les deux moities etaient fausses sur la grille d'alors : la
#    mesure donnait 18/24 et SIX exceptions -- les trois polynomes de degre 2,
#    ET NS 6 [0.5%], NS 7 [0.5%], NS 7 [1.0%], ces trois dernieres de justesse
#    (bornes basses -0,0011 / -0,0006 / -0,0024). La macro publiee
#    \EMPdiffExclNaive valait deja 18 : c'est la PROSE qui avait derive, pas le
#    chiffre, et c'est pourquoi rien ne l'avait signale.
#    ⚠️ Et la phrase « les trois polynomes de degre 2 » redevient VRAIE le
#    18/09 -- non parce qu'elle l'avait toujours ete, mais parce que la grille a
#    change : les trois splines qui faisaient exception en sont sorties. Une
#    phrase qui redevient vraie pour une autre raison reste une phrase qu'on ne
#    recopie pas.
#    Les deux comptes sont EMIS separement -- « sous l'un comme sous l'autre »
#    serait une unanimite deguisee : elle affirme deux comptes en n'en montrant
#    aucun.
##################################################################################

t0 <- proc.time(); cat("=== ct4_emit_tab_contribA_grid.R ===\n")

## Racine deduite de l'emplacement de CE fichier. Meme idiome que 01a l.46-59.
.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .self <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                    error = function(e) NA_character_)
  if (is.na(.self)) {
    .a <- commandArgs(trailingOnly = FALSE)
    .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
    if (length(.f)) .self <- normalizePath(.f[1], winslash = "/")
  }
  if (is.na(.self))
    stop("Impossible de deduire la racine. Poser JMP_ROOT, ou lancer par Rscript.")
  .root <- dirname(dirname(.self))
}
## LA GRILLE VIENT DE 00_parameters.R, ET CE SCRIPT LA VERIFIE (18/09/2026).
## Il ne l'a jamais portee : il comptait ses lignes sur le CSV recu, ce qui
## rendait \EMPgridProdNspec vrai de ce que ct1 avait tourne, jamais de ce que
## l'auteur avait decide. Une grille reduite en amont serait passee dans la
## table, dans la figure et dans la macro sans qu'aucune des trois ne bronche.
## ⚠️ source() de 00_parameters.R cree intermediate_outputs/ et output/ s'ils
## manquent et pose DIR_OUT ; OUT reste defini ici, comme avant, car cet
## emetteur ecrit dans output/ et non dans une redirection de rejeu.
source(file.path(.root, "01_estimation", "00_parameters.R"))
OUT <- file.path(.root, "output")
CA  <- file.path(OUT, "contribA")
CHANNEL <- "prod"
Z <- 1.96

## Les deux entrees, et leur producteur. Un fichier absent n'est pas un
## probleme de chemin : c'est une etape qui n'a pas tourne.
.entrees <- c("test1_delta_north.csv" = "ct1_delta_north_grid.R",
              "tache1_diff_bras.csv"  = "ct2_between_arm_test.R")
.abs <- names(.entrees)[!file.exists(file.path(CA, names(.entrees)))]
if (length(.abs))
  stop("[STOP] entree(s) absente(s) dans output/contribA/ :\n",
       paste(sprintf("  %s  <- produit par %s", .abs, .entrees[.abs]),
             collapse = "\n"))

t1 <- read.csv(file.path(CA, "test1_delta_north.csv"), stringsAsFactors = FALSE)
db <- read.csv(file.path(CA, "tache1_diff_bras.csv"),  stringsAsFactors = FALSE)

## ---- REFUS DE PROVENANCE ----------------------------------------------------
## Cette table JOINT les deux artefacts : les points viennent de ct1, les
## intervalles de l'ecart viennent de ct2. Si les deux n'ont pas tourne sur le
## meme etat d'entree, la jointure reste bien formee -- memes etiquettes, memes
## colonnes -- et la table publie un intervalle qui n'appartient pas a son point.
## Le garde-fou de jointure ci-dessous compare les VALEURS et l'attraperait par
## chance ; l'empreinte l'attrape par construction.
if (!all(c("inputs_fingerprint") %in% names(t1)) ||
    !all(c("inputs_fingerprint") %in% names(db)))
  stop("[STOP] un des deux CSV ne porte pas inputs_fingerprint : artefact\n",
       "  anterieur a la migration du 2026-09-04. Relancer ct1 puis ct2.")
.fp1 <- unique(as.character(t1$inputs_fingerprint))
.fp2 <- unique(as.character(db$inputs_fingerprint))
if (length(.fp1) != 1L || length(.fp2) != 1L || !identical(.fp1, .fp2))
  stop("[STOP] APPARIEMENT CROISE entre ct1 et ct2.\n",
       "  ct1 : ", paste(.fp1, collapse = ", "), "\n",
       "  ct2 : ", paste(.fp2, collapse = ", "), "\n",
       "  Les points et les intervalles ne porteraient pas sur le meme etat\n",
       "  d'entree. Relancer ct1 puis ct2 sur l'etat courant.")
cat(sprintf("  [provenance] empreinte commune : %s\n", .fp1))
t1 <- t1[t1$channel == CHANNEL, ]; db <- db[db$channel == CHANNEL, ]
ta <- t1[t1$arm == "annual", ]; td <- t1[t1$arm == "daily", ]

# ---- garde-fou de jointure : les ETIQUETTES *et* les VALEURS -------------------
# Le lot M9 a montre qu'une verification qui compare etiquettes et format sans
# les VALEURS laisse passer un port incomplet. On controle les deux.
stopifnot(setequal(ta$spec_label, td$spec_label), setequal(ta$spec_label, db$spec_label),
          nrow(ta) == nrow(td), nrow(ta) == nrow(db))
ta <- ta[order(ta$spec_label), ]; td <- td[order(td$spec_label), ]; db <- db[order(db$spec_label), ]

## ---- LA GRILLE RECUE EST-ELLE LA GRILLE DECIDEE ? -----------------------------
## Le controle porte sur les ETIQUETTES, pas sur leur nombre. Quinze lignes dont
## une serait NS 5 au lieu de Poly 4 donneraient le bon compte, la bonne table et
## la bonne figure -- et une grille qui n'est pas celle du manuscrit.
if (!setequal(db$spec_label, GRID_LABELS) || nrow(db) != GRID_NSPEC) {
  .trop <- setdiff(db$spec_label, GRID_LABELS)
  .manq <- setdiff(GRID_LABELS, db$spec_label)
  stop("[STOP] la grille recue de ct1/ct2 n'est pas GRID_SPECS.\n",
       "  recues     : ", nrow(db), " | GRID_SPECS : ", GRID_NSPEC, "\n",
       if (length(.trop)) paste0("  en trop    : ", paste(.trop, collapse = ", "), "\n") else "",
       if (length(.manq)) paste0("  manquantes : ", paste(.manq, collapse = ", "), "\n") else "",
       "  Rejouer ct1 puis ct2 apres le changement de grille.")
}
cat(sprintf("[grille] %d specifications, identiques a GRID_SPECS de 00_parameters.R\n",
            GRID_NSPEC))
.da <- max(abs(ta$delta_north - db$delta_annual)); .dd <- max(abs(td$delta_north - db$delta_daily))
cat(sprintf("[jointure] %d specs | ecart max annual %.3e | daily %.3e\n", nrow(ta), .da, .dd))
if (.da > 1e-10 || .dd > 1e-10) stop("[FATAL] les deux sources ne portent pas les memes points.")

# ---- libelle DERIVE de form/df/winsor, jamais recopie --------------------------
# FORME ABREGEE : « Poly 2 [0.5%] » / « NS 4 [0.5%] ».
# 🔑 C'est une ABREVIATION de la convention des autres tables de specification,
#    pas une convention nouvelle -- et la legende doit le dire ainsi.
# ⚠️ Trois mesures ont conduit ici, dans cet ordre :
#    (i)  « NS » n'apparait NULLE PART ailleurs dans le flux vivant (0 occurrence
#         de NS[0-9]) : la forme abregee doit donc etre EXPLICITEE, sans quoi elle
#         est illisible ;
#    (ii) le mot du papier pour le chiffre est « Deg. » (tab_spec_invariance,
#         tab_iv_latitude, tab_robustness_prod_joint) ou « Degree »
#         (prod_robustness) -- JAMAIS « df ». Une version anterieure de CETTE
#         table portait les 15 seules occurrences de « df » du papier entier ;
#    (iii) la winsorisation est BILATERALE partout ailleurs ([0.5%,0.5%]). La
#         forme abregee n'en garde qu'une borne : la legende doit donc dire
#         qu'elle est SYMETRIQUE, sinon l'abreviation perd une information.
# ⚠️ « Deg. » reste impropre pour une spline naturelle (ns(df=7) n'a pas de degre
#    7, elle est cubique par morceaux) -- mais l'imprecision est celle du papier,
#    et cette table l'abrege sans l'aggraver ni la corriger unilateralement.
if (!all(db$form %in% c("spline", "polynomial")))
  stop(sprintf("[FATAL] forme non prevue : %s. Le libelle serait faux en silence.",
               paste(setdiff(unique(db$form), c("spline", "polynomial")), collapse = ", ")))
.w <- format(db$winsor_lo, trim = TRUE)
spec <- sprintf("%s %s [%s\\%%]", ifelse(db$form == "spline", "NS", "Poly"), db$df, .w)
# la winsorisation abregee ne garde qu'une borne : verifier qu'elle est SYMETRIQUE
if (!all(abs(db$winsor_lo - db$winsor_hi) < 1e-12))
  stop("[FATAL] winsorisation asymetrique : « [0.5%] » perdrait la borne haute en silence.")

# ---- les comptes, EMIS ---------------------------------------------------------
exc <- function(d, se) (d - Z * se) * (d + Z * se) > 0
n_spec  <- nrow(db)
n_ann   <- sum(ta$excludes_zero %in% c(TRUE, "TRUE"))
n_day   <- sum(td$excludes_zero %in% c(TRUE, "TRUE"))
n_diff  <- sum(db$excludes_zero %in% c(TRUE, "TRUE"))
n_naive <- sum(exc(db$diff, db$se_diff_naive))
rho_lo <- min(db$rho); rho_hi <- max(db$rho); gain <- db$se_diff / db$se_diff_naive
cat(sprintf("[comptes] annual %d/%d | daily %d/%d | diff %d/%d | diff naif %d/%d\n",
            n_ann, n_spec, n_day, n_spec, n_diff, n_spec, n_naive, n_spec))

## ---- f_z aux deux temperatures de reference, TOUTE LA GRILLE ----------------
## Source : contribA/spec_curves_daily.csv (ct1, get_effect_curve). Les courbes
## y sont deja normalisees a leur propre maximum -- c'est la convention de
## get_effect_curve, et c'est elle qui rend les lignes comparables entre elles :
## une normalisation commune a une temperature fixe ferait dependre le classement
## du point d'ancrage.
.fcurve <- file.path(OUT, "contribA", "spec_curves_daily.csv")
if (!file.exists(.fcurve))
  stop("[STOP] ", .fcurve, " absent.\n",
       "  <- ct1_delta_north_grid.R (archivage des courbes de la grille)\n",
       "  Sans lui, les colonnes f_z a 25 et 30 degres n'ont pas de source, et on\n",
       "  ne les remplace PAS par celles des douze de tab_spec_invariance : ce\n",
       "  serait melanger deux sous-ensembles dans une meme colonne.")
.cv <- read.csv(.fcurve, stringsAsFactors = FALSE)
.ev <- do.call(rbind, lapply(split(.cv, .cv$model_name), function(d) {
  d <- d[order(d$temperature_C), ]
  v <- approx(d$temperature_C, d$y, xout = c(30, 25), rule = 1)$y
  data.frame(spec_label = d$model_name[1], fz30 = v[1], fz25 = v[2],
             stringsAsFactors = FALSE)
}))
.man <- setdiff(db$spec_label, .ev$spec_label)
if (length(.man))
  stop(sprintf("[STOP] %d specification(s) de la grille absente(s) des courbes : %s",
               length(.man), paste(.man, collapse = ", ")))
cat(sprintf("  [courbes] %d specifications evaluees a 25 et 30 degC\n", nrow(.ev)))

prod <- data.frame(spec_label = db$spec_label, spec_display = spec,
  form = db$form, df = db$df, winsor = db$winsor_lo,
  annual = ta$delta_north, annual_lo = ta$ci_lo, annual_hi = ta$ci_hi,
  annual_excl = ta$excludes_zero,
  daily = td$delta_north, daily_lo = td$ci_lo, daily_hi = td$ci_hi,
  daily_excl = td$excludes_zero,
  diff = db$diff, diff_lo = db$ci_lo, diff_hi = db$ci_hi, diff_excl = db$excludes_zero,
  se_diff = db$se_diff, se_diff_naive = db$se_diff_naive, se_ratio = gain, rho = db$rho,
  diff_excl_naive = exc(db$diff, db$se_diff_naive), stringsAsFactors = FALSE)
## Les BORNES naives, et pas seulement leur verdict. diff_excl_naive dit si
## l'intervalle exclut zero ; il ne dit pas ou il passe. La figure trace le
## segment, donc elle a besoin des bornes -- et les recalculer chez elle ferait
## exactement la recopie entre scripts que la deuxieme regle interdit.
prod$diff_lo_naive <- prod$diff - Z * prod$se_diff_naive
prod$diff_hi_naive <- prod$diff + Z * prod$se_diff_naive
## Controle de coherence : le verdict deja emis doit se lire sur les bornes.
if (!all((prod$diff_lo_naive * prod$diff_hi_naive > 0) == prod$diff_excl_naive))
  stop("[FATAL] diff_excl_naive ne s'accorde pas avec les bornes naives : ",
       "le verdict et le segment ne porteraient pas sur le meme intervalle.")
prod <- merge(prod, .ev, by = "spec_label", all.x = TRUE, sort = FALSE)
prod <- prod[match(db$spec_label, prod$spec_label), ]
write.csv(prod, file.path(OUT, "contribA_grid.csv"), row.names = FALSE)
cat(sprintf("  wrote contribA_grid.csv (%d lignes)\n", nrow(prod)))

# ---- .tex : 5 colonnes, en-tetes COURTS -----------------------------------------
# ---- UN ZERO SIGNE EST UN SIGNE SANS MAGNITUDE ---------------------------------
# sprintf("%+.3f", -0.0004) rend « -0.000 » : le SIGNE survit a un arrondi qui a
# detruit la magnitude, et se lit comme une valeur negative. Sur la colonne
# journaliere c'est le pire endroit possible -- le signe y est la seule
# information, et « -0.000 » l'affirme la ou la valeur ne l'etablit pas.
# ⇒ En deca du dernier chiffre affichable, on rend un zero NON signe.
f3 <- function(v) ifelse(abs(v) < 5e-4, "0.000", sprintf("%+.3f", v))
## LA COLONNE DAILY A QUATRE DECIMALES (lot 34, 25/09/2026). Toutes ses valeurs
## sont a moins de 0,003 de zero : a trois decimales, f3() imprimait 0.000 pour
## NS 3 [2.5%] = -0,000341 et la colonne ne montrait qu'UN negatif la ou
## \EMPgridDailyNneg en compte deux. Valeur ponctuelle : zero sans signe sous le
## plancher d'affichage. ⚠️ f3() sert aussi aux BORNES de la colonne 95 % CI :
## un second suppresseur de zero signe, comme s3() avant le lot 32 -- inscrit a
## la liste post-soumission (aucune borne sous 0,0024 aujourd'hui).
f4 <- function(v) ifelse(abs(v) < 5e-5, "0.0000", sprintf("%+.4f", v))

## Enumeration anglaise d'une liste de libelles, pour une note qui DERIVE ses
## exceptions au lieu de les nommer en dur. Renvoie « A », « A and B », ou
## « A, B and C ». Une liste vide serait une phrase sans objet : on le dit.
.enum_specs <- function(v) {
  v <- as.character(v)
  if (!length(v)) stop("[FATAL] note : aucune exception a enumerer, mais la phrase en annonce.")
  if (length(v) == 1L) return(v)
  paste(paste(utils::head(v, -1), collapse = ", "), "and", utils::tail(v, 1))
}

## ---- la PHRASE des exceptions, derivee ---------------------------------------
## Le manuscrit decrit les exceptions en toutes lettres, pas par une liste de
## libelles. Une phrase en dur est un litteral de resultat : elle a deja menti
## pendant des semaines (« the three degree-two polynomials », alors qu'il y en
## a six). On la DERIVE.
## Deux regimes : quand les exceptions se decrivent par famille -- tous les
## polynomes de degre d, plus des splines -- on rend la description ; sinon on
## rend la liste courte, qui est toujours vraie meme si elle est moins lisible.
.nombre_en_lettres <- function(k) {
  mots <- c("one", "two", "three", "four", "five", "six", "seven", "eight",
            "nine", "ten", "eleven", "twelve")
  if (k >= 1L && k <= length(mots)) mots[k] else as.character(k)
}
.phrase_exceptions <- function(d, grille) {
  ## d      : les lignes qui N'excluent PAS zero sous le calcul naif.
  ## grille : la grille ENTIERE, necessaire pour savoir si « tous les polynomes
  ##          de degre deux » est vrai ou si la description sur-generalise.
  if (!nrow(d)) stop("[FATAL] aucune exception : la phrase serait sans objet.")
  pol <- d[d$form == "polynomial", ]; spl <- d[d$form == "spline", ]
  morceaux <- character(0)
  ## Les polynomes : « the three degree-two polynomials » seulement si TOUS les
  ## polynomes exceptes partagent le meme degre ET que TOUS ceux de ce degre
  ## dans la grille sont exceptes. Sinon la description sur-generalise.
  ok_pol <- nrow(pol) == 0L ||
    (length(unique(pol$df)) == 1L &&
     nrow(pol) == sum(grille$form == "polynomial" & grille$df == pol$df[1]))
  ok_spl <- nrow(spl) == 0L || length(unique(spl$df)) <= 2L
  if (!ok_pol || !ok_spl) return(.enum_specs(d$spec_display))
  if (nrow(pol))
    morceaux <- c(morceaux, sprintf("the %s degree-%s polynomial%s",
                                    .nombre_en_lettres(nrow(pol)),
                                    .nombre_en_lettres(as.integer(pol$df[1])),
                                    if (nrow(pol) > 1L) "s" else ""))
  if (nrow(spl)) {
    dfs <- sort(unique(as.integer(spl$df)))
    morceaux <- c(morceaux, sprintf(
      "%s natural spline%s with %s degree%s of freedom",
      .nombre_en_lettres(nrow(spl)), if (nrow(spl) > 1L) "s" else "",
      paste(vapply(dfs, .nombre_en_lettres, ""), collapse = " or "),
      if (length(dfs) > 1L || dfs[1] > 1L) "s" else ""))
  }
  paste(morceaux, collapse = " and ")
}
L <- c("%% GENERATED by 04_contrast/ct4_emit_tab_contribA_grid.R -- DO NOT EDIT BY HAND.",
  "\\begin{table}[!ht]\\centering\\small",
  # L'OBJET et la GRILLE, dans la legende : deux tables de l'annexe montrent
  # « des specifications x une quantite qui exclut zero », et le lecteur ne peut
  # pas deviner qu'elles ne portent pas sur le meme objet.
  # ⚠️ La grille de tab_spec_invariance est INCLUSE dans celle-ci (12 sur 12,
  #    verifie ligne a ligne) : ce sont la meme famille de specifications, l'une
  #    etendue de df<=4 a df<=7. Ecrire « deux grilles distinctes » serait une
  #    affirmation FAUSSE et verifiable par tout lecteur du replication package.
  paste0("\\caption{Projected change in the northern log productivity fundamental ",
         "under the annual-mean and daily-distribution exposure metrics, their ",
         "difference, and the estimated dose-response behind them. The first three ",
         "columns are the projected CHANGE; the last two are $\\hat f_z$ itself at ",
         "two reference temperatures, each relative to its own specification's ",
         "maximum, in log points. No specification is preferred: the table reports ",
         "the grid, not a selection.}"),
  "\\label{tab:contrast_grid}",
  "\\begin{tabular}{lccccrr}", "\\toprule",
  paste0("Specification & Annual & Daily & Difference & 95\\% CI ",
         "& $\\hat f_z(25^\\circ)$ & $\\hat f_z(30^\\circ)$ \\\\"),
  "\\midrule")
## « --- » et jamais une valeur de bord : hors du support estime, f_z n'existe
## pas. approx(rule = 1) rend NA, on l'imprime tel quel.
.fz <- function(v) if (is.na(v)) "---" else sprintf("%.1f", 100 * v)
for (i in seq_len(n_spec))
  L <- c(L, sprintf("%s & %s & %s & %s & {\\footnotesize $[%s,\\,%s]$} & %s & %s \\\\",
                    spec[i], f3(prod$annual[i]), f4(prod$daily[i]), f3(prod$diff[i]),
                    f3(prod$diff_lo[i]), f3(prod$diff_hi[i]),
                    .fz(prod$fz25[i]), .fz(prod$fz30[i])))
L <- c(L, "\\bottomrule",
  # ⚠️ p{\dimexpr...} et JAMAIS p{Xcm} : le motif en centimetres absolus est celui
  #    qui a tronque la Table 12 (ANCRAGE §4).
  "\\multicolumn{7}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
  # La convention EN TETE de note : sans elle « NS 4 [0.5\%] » n'est pas lisible.
  # ⚠️ Formulee comme une ABREVIATION de la convention des autres tables, pas
  #    comme une definition independante -- deux definitions du meme objet
  #    seraient exactement l'homonymie qu'on cherche a ne pas ouvrir.
  "  \\textit{Note.} Specification labels abbreviate the convention used in the other",
  "  specification tables of this appendix: \\emph{Poly}~$k$ and \\emph{NS}~$k$ are the",
  "  polynomial and natural-spline entries with \\emph{Deg.}~$k$, and the bracketed",
  "  percentage is the winsorization, applied symmetrically to both tails.}",
  ## 22/09/2026 : LE PARAGRAPHE DE COVARIANCE EST RETIRE DE CETTE NOTE.
  ## Il etait ecrit TROIS FOIS -- §2.5, B.4, et ici -- avec les memes chiffres,
  ## et les deux dernieres sur des pages voisines. La note garde ce qui lui est
  ## PROPRE, la convention d'etiquetage ; l'argument de covariance vit dans le
  ## texte de B.4, immediatement au-dessus du tableau.
  ## ⛔ Les grandeurs restent CALCULEES (rho_lo, rho_hi, gain, n_naive) : elles
  ##    alimentent contribA_grid_counts.csv et de la les macros du corps. On ne
  ##    retire qu'un rendu en double, jamais un calcul.
  ## 25/09/2026 (lot 33) : la ligne « being <exceptions>. » qui suivait ici etait
  ## la QUEUE de ce paragraphe (« ...the exceptions being ... »), restee seule
  ## apres son retrait et imprimee en orphelin sous « both tails. » (PDF p.37).
  ## Retiree. La phrase des exceptions vit dans le corps (B.4) par
  ## \EMPdiffExclNaiveExceptions, toujours derivee plus bas (.exc_txt).
  "\\end{tabular}", "\\end{table}")
##################################################################################
# GARDE-FOU DE CONCEPTION — une table d'invariance ne DESIGNE aucune ligne.
# Il s'execute a chaque emission, pas seulement dans une verification humaine :
# si quelqu'un modifie cet emetteur plus tard pour mettre une ligne en gras ou
# ajouter une colonne de poids, l'emission ECHOUE.
# ⚠️ Il tourne AVANT writeLines : une premiere version verifiait apres, et laissait
#    donc le fichier fautif sur disque en s'arretant. Un garde-fou qui s'execute
#    apres l'ecriture ne protege pas l'artefact, il commente sa presence.
##################################################################################
.tex <- paste(L, collapse = "\n")
.corps <- paste(L[(which(L == "\\midrule")[1] + 1):(which(L == "\\bottomrule")[1] - 1)],
                collapse = "\n")
.interdits <- c("\\\\textbf", "\\\\mathbf", "\\\\bf\\b", "\\\\underline", "\\\\rowcolor",
                "\\\\cellcolor", "\\\\dagger", "\\\\star", "\\\\ast", "\\*\\}")
.vus <- .interdits[vapply(.interdits, function(p) grepl(p, .corps), logical(1))]
if (length(.vus))
  stop(sprintf("[FATAL] mise en evidence dans le corps de la table : %s. Une table d'invariance ne designe aucune ligne.",
               paste(.vus, collapse = ", ")))
.colonnes <- c("aic", "AIC", "weight", "poids", "rank", "preferred", "selected")
.entete <- paste(L[grep("Specification &", L)], collapse = " ")   # any(), pas L[grep()] : plusieurs lignes casseraient vapply
.vus2 <- .colonnes[vapply(.colonnes, function(p) any(grepl(p, .entete)), logical(1))]
if (length(.vus2))
  stop(sprintf("[FATAL] colonne de selection dans l'en-tete : %s.", paste(.vus2, collapse = ", ")))

# ---- test de MUTATION du garde-fou, DANS LES DEUX SENS -------------------------
# ⚠️ Deux versions anterieures ont echoue ICI, et chacune sur un defaut different :
#   (1) la premiere verifiait que la MUTATION mute (grepl textbf sur .mut) : elle
#       ne pouvait pas echouer sur un detecteur casse -- le defaut meme que nomme
#       la quatorzieme regle. On teste le DETECTEUR, pas la mutation.
#   (2) la seconde mutait le mot « polynomial », disparu du corps au passage aux
#       libelles abreges (« Poly »). Elle a REFUSE d'emettre, ce qui est le bon
#       comportement, mais son critere etait ancre sur une forme d'USAGE.
# ⇒ La mutation porte desormais sur la PREMIERE LIGNE du corps, quel que soit son
#   contenu, et par concatenation -- pas par sub() sur un motif.
.detecte <- function(txt) any(vapply(.interdits, function(p) grepl(p, txt), logical(1)))
.lignes <- strsplit(.corps, "\n", fixed = TRUE)[[1]]
if (length(.lignes) < 1L || !nzchar(.lignes[1])) stop("[FATAL] corps de table vide : rien a muter.")
.mut_gras <- paste(c(paste0("\\textbf{", .lignes[1], "}"), .lignes[-1]), collapse = "\n")
.mut_marq <- paste(c(paste0(.lignes[1], "$^{\\dagger}$"), .lignes[-1]), collapse = "\n")
if (identical(.mut_gras, .corps) || identical(.mut_marq, .corps))
  stop("[FATAL] la mutation ne mute rien : le test est vide (14e regle).")
if (!.detecte(.mut_gras)) stop("[FATAL] le garde-fou ne voit pas une ligne mise en gras.")
if (!.detecte(.mut_marq)) stop("[FATAL] le garde-fou ne voit pas une ligne marquee d'une dague.")
if (.detecte(.corps))     stop("[FATAL] le garde-fou se declenche sur le corps sain (faux positif).")

# ---- GARDE-FOU : aucun zero signe ne doit subsister -----------------------------
.zs <- "-0\\.0+([^0-9]|$)"
if (grepl(.zs, .corps)) {
  .l <- .lignes[grepl(.zs, .lignes)]
  stop(sprintf("[FATAL] zero signe dans la table (%d ligne(s)) : %s. Un zero signe est un signe sans magnitude.",
               length(.l), paste(substr(.l, 1, 40), collapse = " | ")))
}
# ...et le detecteur doit savoir echouer : mutation dans les deux sens
if (!grepl(.zs, "X & -0.000 & Y")) stop("[FATAL] le detecteur de zero signe ne voit pas « -0.000 ».")
if (!grepl(.zs, "X & -0.00 & Y"))  stop("[FATAL] le detecteur de zero signe ne voit pas « -0.00 ».")
if (grepl(.zs, "X & -0.001 & Y"))  stop("[FATAL] le detecteur de zero signe se declenche sur « -0.001 ».")
.marge <- min(abs(c(prod$annual, prod$daily, prod$diff)))
cat(sprintf("  [garde-fou] aucun zero signe | valeur la plus proche de zero : %.6f (%.1fx le seuil d'arrondi)\n",
            .marge, .marge / 5e-4))
cat(sprintf("  [garde-fou] aucune ligne designee | %d motifs | mutation gras VUE, dague VUE, sain NON declenche\n",
            length(.interdits) + length(.colonnes)))

writeLines(L, file.path(OUT, "tab_contrast_grid.tex"))
cat("  wrote tab_contrast_grid.tex\n")

##################################################################################
# FIGURE — output/fig_contrast_grid.pdf
#
# Le contraste passe de l'annexe au corps. La table reste : elle porte les
# valeurs, la figure porte la FORME de l'argument -- vingt-quatre fois le meme
# ordre entre trois intervalles.
#
# ---------------------------------------------------------------------------
# UNE SEULE SOURCE, ET ELLE EST RELUE SUR DISQUE
# ---------------------------------------------------------------------------
# La figure ne recalcule rien : elle relit contribA_grid.csv, le fichier que la
# table vient d'ecrire. Le faire depuis `prod` en memoire serait plus court et
# strictement moins bon -- une erreur d'ECRITURE du CSV resterait invisible,
# alors qu'ici la figure et la table divergeraient et l'assertion le dirait.
#
# ---------------------------------------------------------------------------
# CE QUE LES ASSERTIONS PROTEGENT
# ---------------------------------------------------------------------------
# Les comptes publies (\EMPgridAnnualExcl, \EMPgridDailyExcl, \EMPdiffExcl) sont
# des enonces SUR LA FIGURE : « tant d'intervalles excluent zero ». On compte
# donc les segments TRACES, pas les lignes du tableau, et on exige l'egalite.
# Un segment qui ne serait pas trace -- une borne NA, une ligne perdue par un
# ordre -- ferait diverger les deux comptes et arreterait l'emission.
##################################################################################

## Les etiquettes d'horizon viennent de 00_parameters.R, jamais epelees ici.
## Source dans un environnement SEPARE : 00_parameters.R definit une centaine de
## noms (SPEC_CONFIGS entre autres) et ct4 tourne apres ct1, qui ECRASE
## deliberement cette grille. Un source() nu les ferait rentrer par la fenetre.
.penv <- new.env()
sys.source(file.path(.root, "01_estimation", "00_parameters.R"), envir = .penv)
HORIZON_LABEL <- get("HORIZON_LABEL", envir = .penv)

g <- read.csv(file.path(OUT, "contribA_grid.csv"), stringsAsFactors = FALSE)
if (nrow(g) != n_spec)
  stop(sprintf("[FATAL] contribA_grid.csv relu : %d lignes, %d attendues.",
               nrow(g), n_spec))

## ---- ordre : famille, puis degre/df, puis winsorisation -----------------------
## Les polynomes d'abord, les splines naturelles ensuite. match() sur un vecteur
## explicite, et non factor() : un niveau inattendu donnerait NA et l'ordre
## silencieusement faux plutot qu'une erreur.
.fam <- match(g$form, c("polynomial", "spline"))
if (anyNA(.fam)) stop("[FATAL] forme hors {polynomial, spline} : ordre indefini.")
g <- g[order(.fam, g$df, g$winsor), ]
n <- nrow(g)
## HAS_NAIVE : l intervalle d independance n existe que si ct2 a ecrit ses
## deux colonnes. Son absence retire une serie, elle n invente rien.
HAS_NAIVE <- all(c("diff_lo_naive", "diff_hi_naive") %in% names(g))

## ---- etiquettes courtes d'abscisse : « P2 .5 », ... , « NS7 2.5 » -----------
## La winsorisation s'ecrit sans zero de tete ni decimale inutile : c'est la
## seule concession de place, et elle ne cree pas d'ambiguite (trois seuils).
.wtxt <- function(w) sub("^0\\.", ".", sub("\\.0$", "", format(w, trim = TRUE)))
xlab_spec <- sprintf("%s%s %s",
                     ifelse(g$form == "spline", "NS", "P"),
                     g$df, .wtxt(g$winsor))

xpos <- seq_len(n)                       # la premiere specification a GAUCHE
DX   <- 0.15                             # decalage horizontal des deux bras
.kpol <- sum(g$form == "polynomial")     # derniere colonne polynomiale

COL <- c(annual = "#1f5fbf", daily = "#c8302b",
         diff   = "#000000", naive = gray(0.72))

## ---- bornes verticales ------------------------------------------------------
## Chaque panneau prend SES intervalles, plus la marge que sa legende exige.
## Une legende calee dans un ylim dimensionne sur les seules donnees mordrait
## sur un segment des qu'une borne bougerait ; la marge est donc explicite.
.yt <- c(g$annual_lo, g$annual_hi, g$daily_lo, g$daily_hi)
.yb <- c(g$diff_lo, g$diff_hi, g$diff_lo_naive, g$diff_hi_naive)
if (anyNA(.yt) || anyNA(.yb))
  stop("[FATAL] borne d'intervalle manquante : un segment ne serait pas trace.")
## Haut : la legende est en BAS a droite -> marge sous le minimum.
.rt <- range(.yt); ylim_top <- c(.rt[1] - 0.30 * diff(.rt), .rt[2] + 0.06 * diff(.rt))
## Bas : la legende est en HAUT a gauche -> marge au-dessus du maximum.
.rb <- range(.yb); ylim_bot <- c(.rb[1] - 0.06 * diff(.rb), .rb[2] + 0.34 * diff(.rb))

xlim <- c(0.4, n + 0.6)

fig <- file.path(OUT, "fig_contrast_grid.pdf")
## 16 x 10 cm, police 9 pt, PDF vectoriel.
pdf(fig, width = 16 / 2.54, height = 10 / 2.54, pointsize = 9)
## Deux panneaux empiles, ABSCISSE COMMUNE : seul celui du bas porte les
## etiquettes. Le haut garde 1.4 ligne de marge pour « Polynomials » /
## « Natural splines », le bas 3.2 pour les etiquettes a 45 degres.
layout(matrix(1:2, nrow = 2), heights = c(1, 1))
par(oma = c(0, 0, 0, 0), mgp = c(2.1, 0.6, 0), tcl = -0.25, xaxs = "i")

## ============================ PANNEAU DU HAUT ==============================
par(mar = c(0.3, 4.6, 1.5, 0.6))
plot(NA, xlim = xlim, ylim = ylim_top, axes = FALSE, xlab = "", ylab = "")
abline(h = 0, col = "black", lwd = 0.6)
if (.kpol > 0 && .kpol < n)
  abline(v = .kpol + 0.5, col = gray(0.72), lwd = 0.6, lty = 3)

## lend = 1 : embouts francs. Le defaut arrondi ajoute un demi-trait de chaque
## cote, donc allonge visuellement chaque intervalle.
segments(xpos - DX, g$annual_lo, xpos - DX, g$annual_hi,
         col = COL["annual"], lwd = 1, lend = 1)
segments(xpos + DX, g$daily_lo,  xpos + DX, g$daily_hi,
         col = COL["daily"],  lwd = 1, lend = 1)
points(xpos - DX, g$annual, pch = 16, cex = 0.55, col = COL["annual"])
points(xpos + DX, g$daily,  pch = 16, cex = 0.55, col = COL["daily"])

axis(2, las = 1, lwd = 0, lwd.ticks = 0.5, cex.axis = 0.8)
mtext("Log points", side = 2, line = 2.9, cex = 0.8)

## Les deux familles, nommees AU-DESSUS du panneau du haut, centrees sur leur
## propre plage de colonnes -- donc jamais a re-caler si la grille s'allonge.
if (.kpol > 0 && .kpol < n) {
  mtext("Polynomials",     side = 3, line = 0.15, cex = 0.8,
        at = (1 + .kpol) / 2)
  mtext("Natural splines", side = 3, line = 0.15, cex = 0.8,
        at = (.kpol + 1 + n) / 2)
}

legend("bottomright", bty = "n", cex = 0.72,
       legend = c("Annual mean", "Daily distribution"),
       col = c(COL["annual"], COL["daily"]),
       pch = 16, pt.cex = 0.55, lty = 1, lwd = 1, seg.len = 1.6,
       y.intersp = 1.15)

## ============================= PANNEAU DU BAS ==============================
par(mar = c(3.2, 4.6, 0.3, 0.6))
plot(NA, xlim = xlim, ylim = ylim_bot, axes = FALSE, xlab = "", ylab = "")
abline(h = 0, col = "black", lwd = 0.6)
if (.kpol > 0 && .kpol < n)
  abline(v = .kpol + 0.5, col = gray(0.72), lwd = 0.6, lty = 3)

## L'intervalle d'INDEPENDANCE d'abord : il passe DERRIERE le noir, et il est
## plus epais, de sorte que le noir se lise comme un retrecissement de lui.
if (HAS_NAIVE)
  segments(xpos, g$diff_lo_naive, xpos, g$diff_hi_naive,
           col = COL["naive"], lwd = 2.2, lend = 1)
segments(xpos, g$diff_lo, xpos, g$diff_hi,
         col = COL["diff"], lwd = 1, lend = 1)
points(xpos, g$diff, pch = 18, cex = 0.95, col = COL["diff"])

axis(2, las = 1, lwd = 0, lwd.ticks = 0.5, cex.axis = 0.8)
mtext("Log points", side = 2, line = 2.9, cex = 0.8)

## Etiquettes a 45 degres, ancrees a droite : l'ancre est le haut-droit du
## texte, donc l'etiquette pend sous sa colonne sans jamais la deborder.
text(x = xpos, y = ylim_bot[1] - 0.035 * diff(ylim_bot),
     labels = xlab_spec, srt = 45, adj = c(1, 1), cex = 0.7, xpd = NA)

.leg_b <- c("Difference, covariance carried",
            "Difference, arms treated as independent")
legend("topleft", bty = "n", cex = 0.72,
       legend = .leg_b,
       col = c(COL["diff"], COL["naive"]),
       pch = c(18, NA), pt.cex = 0.95,
       lty = 1, lwd = c(1, 2.2), seg.len = 1.6, y.intersp = 1.15)

## La legende ne doit toucher AUCUN segment : on le verifie au lieu de le
## supposer. Sa boite est mesuree dans le device, puis confrontee aux colonnes
## qu'elle recouvre. Une marge posee en dur redeviendrait fausse a la premiere
## retouche de police.
.lb <- legend("topleft", plot = FALSE, bty = "n", cex = 0.72,
              legend = .leg_b, col = c(COL["diff"], COL["naive"]),
              pch = c(18, NA), lty = 1, lwd = c(1, 2.2), seg.len = 1.6,
              y.intersp = 1.15)
.cov <- which(xpos >= .lb$rect$left & xpos <= .lb$rect$left + .lb$rect$w)
.top <- if (length(.cov))
          max(c(g$diff_hi[.cov], if (HAS_NAIVE) g$diff_hi_naive[.cov])) else -Inf
if (length(.cov) && .top >= .lb$rect$top - .lb$rect$h)
  stop(sprintf(paste("[FATAL] la legende du panneau du bas mord sur un segment :",
                     "sommet %.4f sous %.4f, colonnes %s."),
               .lb$rect$top - .lb$rect$h, .top, paste(.cov, collapse = ",")))

invisible(dev.off())

##################################################################################
# ASSERTIONS — les comptes publies sont des enonces SUR LES SEGMENTS TRACES
##################################################################################
.excl_seg <- function(lo, hi) sum(lo * hi > 0)
n_ann_fig   <- .excl_seg(g$annual_lo, g$annual_hi)
n_day_fig   <- .excl_seg(g$daily_lo,  g$daily_hi)
n_diff_fig  <- .excl_seg(g$diff_lo,   g$diff_hi)
n_naive_fig <- if (HAS_NAIVE) .excl_seg(g$diff_lo_naive, g$diff_hi_naive) else NA_integer_
n_dayneg    <- sum(g$daily < 0)

.check <- function(nom, vu, attendu) {
  if (!identical(as.integer(vu), as.integer(attendu)))
    stop(sprintf("[FATAL] %s : %d segments traces, %d annonces par la table.",
                 nom, vu, attendu))
}
.check("segments bleus (annual)", n_ann_fig,   n_ann)
.check("segments rouges (daily)", n_day_fig,   n_day)
.check("segments noirs (diff)",   n_diff_fig,  n_diff)
if (HAS_NAIVE) .check("segments gris (diff naif)", n_naive_fig, n_naive)
cat(sprintf("  [figure] segments excluant zero : annual %d/%d | daily %d/%d | diff %d/%d | diff naif %d/%d\n",
            n_ann_fig, n, n_day_fig, n, n_diff_fig, n, n_naive_fig, n))

## Le detecteur doit savoir echouer : un intervalle qui CONTIENT zero ne compte pas.
if (.excl_seg(c(-1, 1), c(1, 2)) != 1L)
  stop("[FATAL] le compteur d'exclusion ne distingue pas [-1,1] de [1,2].")

##################################################################################
# REPRODUCTIBILITE OCTET A OCTET
# pdf() ecrit /CreationDate et /ModDate. Deux emissions identiques donneraient
# deux md5 differents, et « la figure a change » deviendrait ininterpretable.
# Les deux dates font QUATORZE caracteres : on les remplace par un horodatage
# constant de meme longueur, donc sans deplacer un seul offset de la table xref.
##################################################################################
.stamp_pdf <- function(path, stamp = "19700101000000") {
  if (nchar(stamp) != 14L) stop("[FATAL] horodatage de longueur != 14 : xref deplacee.")
  n  <- file.info(path)$size
  rw <- readBin(path, "raw", n)
  key <- charToRaw("(D:")
  lim <- min(n - 16L, 4000L)          # l'Info dict est en tete de fichier
  hits <- 0L
  for (i in seq_len(lim))
    if (identical(rw[i:(i + 2L)], key)) {
      rw[(i + 3L):(i + 16L)] <- charToRaw(stamp); hits <- hits + 1L
    }
  if (hits == 0L)
    stop("[FATAL] aucune date trouvee dans le PDF : le PDF ne serait pas reproductible.")
  writeBin(rw, path)
  hits
}
.nd <- .stamp_pdf(fig)
cat(sprintf("  wrote fig_contrast_grid.pdf (%d horodatage(s) neutralise(s))\n", .nd))

##################################################################################
# LES COMPTES, EMIS — depuis la structure qui porte la figure ET la table
#
# Ces six nombres etaient calcules par emit_values_empirics.py sur les CSV
# AMONT (test1_delta_north.csv, tache1_diff_bras.csv). Ils le sont desormais
# ici, sur `g`, c'est-a-dire sur l'objet meme dont la figure trace les segments
# et dont la table imprime les lignes -- ce que la deuxieme regle demande.
# emit_values_empirics.py lit ce fichier et ne recompte plus.
##################################################################################
## La PHRASE des exceptions est emise comme les comptes le sont : derivee de la
## meme grille, au meme passage. Le manuscrit la porte en toutes lettres et ne
## doit plus la redire lui-meme.
.exc_txt <- .phrase_exceptions(g[!g$diff_excl_naive, ], g)
cat(sprintf("  [exceptions] %s\n", .exc_txt))
## Garde-fou : la phrase doit compter ce que le compte compte. Le nombre ecrit
## en lettres dans la phrase est reconstruit et compare a n - n_naive.
.n_exc <- sum(!g$diff_excl_naive)
if (.n_exc != n - n_naive_fig)
  stop(sprintf("[FATAL] %d exceptions listees contre %d annoncees par le compte.",
               .n_exc, n - n_naive_fig))

.counts <- data.frame(
  macro = c("EMPgridProdNspec", "EMPgridAnnualExcl", "EMPgridDailyExcl",
            "EMPgridDailyNneg", "EMPdiffExcl", "EMPdiffExclNaive",
            "EMPdiffExclNaiveExceptions"),
  value = c(as.character(c(n, n_ann_fig, n_day_fig, n_dayneg, n_diff_fig,
                           n_naive_fig)), .exc_txt),
  source = "contribA_grid.csv (ct4_emit_tab_contribA_grid.R, meme passage que tab_contrast_grid.tex et fig_contrast_grid.pdf)",
  stringsAsFactors = FALSE)
write.csv(.counts, file.path(OUT, "contribA_grid_counts.csv"), row.names = FALSE)
cat("  wrote contribA_grid_counts.csv\n")


##################################################################################
# LA TABLE DE ROBUSTESSE SUR LES QUINZE SPECIFICATIONS
#   output/prod_joint_tests_grid_{daily,annual}.csv   (source machine)
#   output/tab_robustness_prod_joint.tex              (l'exhibit)
#
# POURQUOI LE PRODUCTEUR EST ICI, ET PAS DANS 01a (18/09/2026, lot 10)
# --------------------------------------------------------------------
# La table porte douze lignes depuis toujours parce que 01a estime
# SPEC_CONFIGS -- la grille du CORPS, quatre formes x trois winsorisations. La
# grille du CONTRASTE en compte quinze : les trois polynomes de degre 4 lui
# manquaient. Elle disait donc « across specifications » d'un sous-ensemble.
#
# Les trois ajustements manquants existent deja, et il n'y a rien a re-estimer :
# ct1 parcourt GRID_SPECS et appelle le MEME chemin que 01a -- prepare_diffs()
# puis run_stage3() avec les controles gravitaires du commerce -- de sorte que
# p_joint et p_nonlin y sont les memes tests, au meme clustering bidimensionnel
# origine x destination. Le seul obstacle etait l'ORDRE : 01a tourne avant ct1
# et ne peut pas lire une sortie de ct1. On deplace donc l'emetteur en aval, la
# ou les quinze lignes sont disponibles, plutot que de faire estimer trois
# modeles de plus a un script qui lit douze gigaoctets de BACI pour cela.
#
# ⚠️ LE CONTROLE QUI RENDAIT CE DEPLACEMENT PUBLIABLE EST PERMANENT CI-DESSOUS.
# « ct1 refait le meme calcul que 01a » est une affirmation, pas une preuve, et
# une re-implementation voisine passerait toutes les verifications de forme. On
# compare donc, VALEUR PAR VALEUR, les douze specifications communes aux deux
# sources. Mesure du 18/09/2026 : ecart maximal nul sur p_joint, p_nonlin, AIC
# et N, dans les deux bras. Le jour ou les deux chemins divergeront, la table
# ne sortira pas.
##################################################################################

.pj_arm <- function(d) {
  d <- d[order(d$spec_label), ]
  data.frame(model_name = d$spec_label, spec = d$form, degree = as.integer(d$df),
             wval_lower = d$winsor_lo, wval_upper = d$winsor_hi,
             p_joint = as.numeric(d$p_joint), p_nonlin = as.numeric(d$p_nonlin),
             nobs = as.integer(d$n_obs), aic = as.numeric(d$aic),
             stringsAsFactors = FALSE)
}
gd <- .pj_arm(td)          # bras journalier, quinze lignes
ga <- .pj_arm(ta)          # bras annuel, quinze lignes

## ---- CONTROLE DE NON-REGRESSION CONTRE 01a ---------------------------------
## Polarite : ARRET. Ce n'est pas une propriete du resultat, c'est la condition
## qui autorise a melanger douze lignes d'une source et trois d'une autre dans
## un meme tableau. Les DOUZE communes doivent coincider AU DERNIER CHIFFRE.
.legacy <- c(daily = "prod_joint_tests_daily.csv",
             annual = "prod_joint_tests_annual.csv")
for (.arm in names(.legacy)) {
  .f <- file.path(OUT, .legacy[[.arm]])
  if (!file.exists(.f)) {
    cat(sprintf("  [01a] %s absent : controle de non-regression NON CONCLUANT.\n",
                .legacy[[.arm]]))
    stop("[STOP] la table des quinze ne s'emet pas sans l'ancre des douze.\n",
         "  ", .legacy[[.arm]], " est produit par 01a_productivity_drf.R.\n",
         "  Un controle qui ne peut pas conclure ne passe pas au vert.")
  }
  .old <- read.csv(.f, stringsAsFactors = FALSE)
  .new <- if (.arm == "daily") gd else ga
  .com <- intersect(.old$model_name, .new$model_name)
  if (!setequal(.old$model_name, .com))
    stop("[STOP] ", .legacy[[.arm]], " porte une specification absente de\n",
         "  GRID_SPECS : ", paste(setdiff(.old$model_name, .com), collapse = ", "),
         "\n  Les deux grilles ne sont pas emboitees, la table melangerait\n",
         "  deux echantillons de specifications.")
  .o <- .old[match(.com, .old$model_name), ]
  .n <- .new[match(.com, .new$model_name), ]
  .ec <- sapply(c("p_joint", "p_nonlin", "aic"),
                function(k) max(abs(.o[[k]] - .n[[k]])))
  .en <- max(abs(as.numeric(.o$nobs) - as.numeric(.n$nobs)))
  if (any(.ec > 0) || .en > 0)
    stop(sprintf(paste0("[STOP] bras %s : ct1 et 01a ne rendent pas les memes\n",
                        "  tests sur les %d specifications communes.\n",
                        "  ecart max  p_joint %.3e | p_nonlin %.3e | AIC %.3e | N %.0f\n",
                        "  Les quinze lignes ne peuvent pas etre presentees comme\n",
                        "  une seule grille tant que les douze ne coincident pas."),
                 .arm, length(.com), .ec[["p_joint"]], .ec[["p_nonlin"]],
                 .ec[["aic"]], .en))
  cat(sprintf(paste0("  [non-regression %-6s] %2d specifications communes ",
                     "reproduites au bit (p_joint, p_nonlin, AIC, N)\n"),
              .arm, length(.com)))
}

write.csv(gd, file.path(OUT, "prod_joint_tests_grid_daily.csv"), row.names = FALSE)
write.csv(ga, file.path(OUT, "prod_joint_tests_grid_annual.csv"), row.names = FALSE)
cat("  wrote prod_joint_tests_grid_{daily,annual}.csv\n")

## ---- l'exhibit --------------------------------------------------------------
.m <- merge(gd[, c("model_name", "spec", "degree", "wval_lower", "wval_upper",
                   "p_joint", "p_nonlin", "nobs")],
            ga[, c("model_name", "p_joint", "p_nonlin", "nobs")],
            by = "model_name", suffixes = c("_d", "_a"))
.m <- .m[order(.m$spec, .m$degree, .m$wval_lower), ]
## Le meme garde-fou de recopie que portait 01a : deux colonnes identiques a la
## precision machine sont la signature d'un bras copie sur l'autre, pas d'un
## resultat.
if (isTRUE(all.equal(.m$p_nonlin_d, .m$p_nonlin_a)))
  stop("[FATAL] p_nonlin IDENTIQUE entre les deux bras -- signature de la recopie.")

## Quatre decimales et pas trois : a trois, une cellule a 0.0325 serait indecise.
## format() est VECTORISE et aligne ses sorties sur une largeur commune :
## format(c(0.5, 1, 2.5)) rend "1.0" la ou la table portait "1" depuis 01a.
## On formate donc element par element, sinon le seul passage a ct4 deplacerait
## six libelles de winsorisation sans qu'aucun chiffre ait bouge.
.fmtw <- function(v) vapply(v, function(z) format(z, trim = TRUE), "")
.fmtp <- function(x) ifelse(is.na(x), "---",
                     ifelse(x < 0.001, "$<0.001$", sprintf("%.4f", x)))
.nd <- unique(.m$nobs_d); .na <- unique(.m$nobs_a)
.nsent <- (if (length(.nd) == 1L && length(.na) == 1L && isTRUE(.nd == .na))
  sprintf("$N = %s$ in every cell, identical across the two metrics.",
          formatC(.nd, format = "d", big.mark = "{,}"))
else "Sample sizes differ across cells; see the replication package.")

.nrow <- nrow(.m); .ncell <- 2L * .nrow
.allnl <- c(.m$p_nonlin_d, .m$p_nonlin_a)
.alljt <- c(.m$p_joint_d,  .m$p_joint_a)

## LA NOTE EST DERIVEE, PHRASE PAR PHRASE. Le compte de rejets a 5 % etait
## « in all cells » tant que la grille s'arretait au degre 3 ; il ne l'est plus
## avec les polynomes de degre 4, dont la non-linearite est la plus faible de la
## grille sous la metrique journaliere. La note DIT le compte qu'elle mesure, et
## c'est exactement pourquoi elle est construite ici et pas ecrite a la main.
.n05 <- sum(.allnl < 0.05, na.rm = TRUE)
.sent05 <- (if (.n05 == .ncell)
  sprintf(paste0("Non-linearity is significant at $p<0.05$ in all %d ",
                 "specification--metric cells ($%d$ specifications $\\times$ ",
                 "two exposure metrics)"), .ncell, .nrow)
else
  sprintf(paste0("Non-linearity is significant at $p<0.05$ in %d of the %d ",
                 "specification--metric cells ($%d$ specifications $\\times$ ",
                 "two exposure metrics)"), .n05, .ncell, .nrow))

writeLines(c(
  "%% GENERATED by 04_contrast/ct4_emit_tab_contribA_grid.R -- DO NOT EDIT BY HAND.",
  "%% Source: prod_joint_tests_grid_{daily,annual}.csv, derived from",
  "%% output/contribA/test1_delta_north.csv (ct1_delta_north_grid.R).",
  "\\begin{table}[H]", "\\centering\\small",
  paste0("\\caption{Robustness: productivity dose-response across specifications,",
         " daily distribution versus annual mean}"),
  "\\label{tab:robustness_prod_joint}",
  "\\begin{tabular}{llccc}", "\\toprule",
  "& & & \\multicolumn{2}{c}{$p$ (non-linearity)} \\\\", "\\cmidrule(lr){4-5}",
  "Functional form & Deg. & Winsorization & Daily distribution & Annual mean \\\\",
  "\\midrule",
  sprintf("%s & %d & $[%s\\%%, %s\\%%]$ & %s & %s \\\\",
          ifelse(.m$spec == "polynomial", "Polynomial", "Natural spline"),
          .m$degree, .fmtw(.m$wval_lower), .fmtw(.m$wval_upper),
          .fmtp(.m$p_nonlin_d), .fmtp(.m$p_nonlin_a)),
  "\\bottomrule",
  "\\multicolumn{5}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
  "\\textit{Note.} Stage~3 PPML estimates with product$\\times$year fixed effects",
  "and bilateral gravity controls. Standard errors clustered by origin and",
  paste("destination.", .nsent),
  ## Lot 85x (11/10/2026) : casse et ponctuation de la note portees ici depuis la
  ## retouche manuelle de l'auteur sur la table emise (09/10, lot 85u).
  "Entries are $p$-values for the non-linearity test: the joint Wald test on",
  "the degree-$\\geq 2$ climate terms alone.",
  ## 22/09/2026 : seuil secondaire aligne sur celui de la PROSE de B.3, qui
  ## porte 1 % en macro (EMPprodNnonlinRejOnePct). Deux seuils secondaires pour
  ## un meme tableau faisaient lire deux resultats la ou il y en a un.
  sprintf(paste0("%s, and at $p<0.01$ in %d of the %d (%d of %d under the ",
                 "daily distribution, %d of %d under the annual mean)."),
          .sent05,
          sum(.allnl < 0.01, na.rm = TRUE), .ncell,
          sum(.m$p_nonlin_d < 0.01, na.rm = TRUE), .nrow,
          sum(.m$p_nonlin_a < 0.01, na.rm = TRUE), .nrow),
  sprintf(paste0("A separate Wald test on all climate terms jointly (linear ",
                 "and higher-order, which is not the quantity tabulated here) ",
                 "rejects at $p<0.001$ in %d of the %d cells."),
          sum(.alljt < 0.001, na.rm = TRUE), .ncell),
  "The two metrics are estimated on the same sample, with the same dependent",
  "variable, the same fixed effects and the same controls: only the",
  "right-hand-side exposure measure differs.}",
  "\\end{tabular}", "\\end{table}"),
  file.path(OUT, "tab_robustness_prod_joint.tex"), useBytes = TRUE)
##################################################################################
# LES COMPTES DE LA TABLE, EMIS -- ET LA PHRASE DES EXCEPTIONS AVEC EUX
#
# « linearity is rejected in ALL cells » etait vrai de la grille de douze et
# faux de celle de quinze. Une affirmation d'unanimite ne peut etre vraie que
# par un COMPTE, et une liste d'exceptions ne peut etre vraie que si elle est
# DERIVEE : c'est exactement le defaut ferme le 16/09 sur
# \EMPdiffExclNaiveExceptions, ou le manuscrit disait « the three degree-two
# polynomials » pour six exceptions. On reutilise le meme constructeur, avec
# une seule addition -- le BRAS -- puisque les cellules exceptees sont ici des
# couples (specification, metrique) et non des specifications.
##################################################################################
.exc_idx <- which(.allnl >= 0.05)
.exc_arm <- rep(c("daily", "annual"), each = .nrow)[.exc_idx]
.exc_lab <- rep(.m$model_name, 2L)[.exc_idx]
.exc_frm <- rep(.m$spec, 2L)[.exc_idx]
.exc_df  <- rep(.m$degree, 2L)[.exc_idx]
.exc_dsp <- sprintf("%s %s [%s\\%%]",
                    ifelse(.exc_frm == "spline", "NS", "Poly"), .exc_df,
                    .fmtw(rep(.m$wval_lower, 2L)[.exc_idx]))
.ARM_WORD <- c(daily = "the daily distribution", annual = "the annual mean")

if (!length(.exc_idx)) {
  ## Pas d'exception : la phrase doit le DIRE, pas rester vide. Une chaine vide
  ## dans une macro rendrait « the exceptions at 5 % are . » dans le manuscrit.
  .nl_exc_txt <- "there are no exceptions"
  .nl_exc_lo <- NA_real_; .nl_exc_hi <- NA_real_
} else {
  ## La grille de reference du constructeur est la grille du CANAL et du BRAS
  ## considere : « the three degree-four polynomials » n'est exact que si les
  ## trois Poly 4 de CE bras sont toutes exceptees.
  .g1 <- data.frame(form = .exc_frm, df = .exc_df, spec_display = .exc_dsp,
                    stringsAsFactors = FALSE)
  .gall <- data.frame(form = .m$spec, df = .m$degree,
                      spec_display = sprintf("%s %s [%s\\%%]",
                                             ifelse(.m$spec == "spline", "NS", "Poly"),
                                             .m$degree, .fmtw(.m$wval_lower)),
                      stringsAsFactors = FALSE)
  if (length(unique(.exc_arm)) == 1L) {
    .nl_exc_txt <- sprintf("%s under %s",
                           .phrase_exceptions(.g1, .gall),
                           .ARM_WORD[[.exc_arm[1]]])
  } else {
    ## Exceptions reparties sur les deux bras : on enumere les CELLULES, sans
    ## quoi la description perdrait la metrique et sur-generaliserait.
    .nl_exc_txt <- .enum_specs(sprintf("%s under %s", .exc_dsp,
                                       .ARM_WORD[.exc_arm]))
  }
  .nl_exc_lo <- min(.allnl[.exc_idx]); .nl_exc_hi <- max(.allnl[.exc_idx])
}
cat(sprintf("  [exceptions non-lin.] %s\n", .nl_exc_txt))

.pjcounts <- data.frame(
  macro = c("EMPprodNnonlinRejFivePct", "EMPprodNnonlinExceptions",
            "EMPprodPnonlinExcLo", "EMPprodPnonlinExcHi"),
  value = c(as.character(.n05), .nl_exc_txt,
            if (is.na(.nl_exc_lo)) "---" else sprintf("%.3f", .nl_exc_lo),
            if (is.na(.nl_exc_hi)) "---" else sprintf("%.3f", .nl_exc_hi)),
  source = paste0("prod_joint_tests_grid_{daily,annual}.csv ",
                  "(ct4_emit_tab_contribA_grid.R, meme passage que ",
                  "tab_robustness_prod_joint.tex)"),
  stringsAsFactors = FALSE)
write.csv(.pjcounts, file.path(OUT, "prod_joint_tests_grid_counts.csv"),
          row.names = FALSE)
cat("  wrote prod_joint_tests_grid_counts.csv\n")

## Garde-fou : le compte de la NOTE et le compte de la MACRO sont le meme
## nombre, donc ils sortent du meme objet. On verifie qu'ils s'accordent avec
## la liste d'exceptions -- trois enonces, une seule mesure.
if (.n05 + length(.exc_idx) != .ncell)
  stop(sprintf("[FATAL] %d rejets + %d exceptions != %d cellules.",
               .n05, length(.exc_idx), .ncell))

cat(sprintf(paste0("  wrote tab_robustness_prod_joint.tex : %d lignes, %d cellules",
                   " | non-lin. <0.05 : %d/%d | <0.01 : %d/%d\n"),
            .nrow, .ncell, .n05, .ncell,
            sum(.allnl < 0.01, na.rm = TRUE), .ncell))
if (.n05 < .ncell)
  cat(sprintf(paste0("  [!] la non-linearite n'est PAS rejetee a 5 %% dans %d ",
                     "cellule(s) : %s\n"),
              .ncell - .n05,
              paste(sprintf("%s (%s, p=%.4f)",
                            c(.m$model_name, .m$model_name)[.allnl >= 0.05],
                            rep(c("daily", "annual"), each = .nrow)[.allnl >= 0.05],
                            .allnl[.allnl >= 0.05]), collapse = " ; ")))

cat(sprintf("\nTemps total: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))
cat("=== Done ===\n")

