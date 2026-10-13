##################################################################################
# ct1_delta_north_grid.R — CONTRASTE ANNUEL/JOURNALIER, ETAPE 1 : la grille
# etendue, avec intervalles.
#
# Reads:  intermediate_outputs/BACI_raw_pairs_annual.rds, step2_migration.rds,
#         climate_long_annual.rds, climate_long.rds, CMIP6_tas.csv,
#         intermediate_outputs/crh/amen_prod.csv
# Writes: output/contribA/test1_delta_north.csv
# Sourced by: nothing. ct2 le lit ; ct3 en derive les bornes de support.
#
# CE QUE CE SCRIPT CALCULE
# ------------------------
# Pour chaque (canal, bras, specification) : l'INCREMENT NORDIQUE Delta_north,
# son erreur-type analytique, son intervalle a 95 %, et la masse-jour nordique
# hors support qui conditionne l'interpretation de cet intervalle.
#
#   canal : prod (productivite via commerce)  |  amen (amenite via migration)
#   bras  : daily (distribution journaliere)  |  annual (moyenne annuelle)
#   spec  : 5 formes x 3 winsorisations = 15 modeles par (canal, bras)
#           (GRID_SPECS de 00_parameters.R -- decision d'auteur du 18/09/2026 ;
#            ce script ne porte plus de liste de specifications)
#
# L'OBJET
# -------
# Bras DAILY   : Delta_north = sum_c omega_c [ sum_b (n_b^fut(c) - n_b^base(c)) . b(T_b) ]' beta
# Bras ANNUAL  : Delta_north = sum_c omega_c [ b_ann(Tbar_c^fut) - b_ann(Tbar_c^base) ]' beta
#
# ou c parcourt les cellules 1° de latitude y >= 45N, omega_c est la population 2005
# de amen_prod.csv (les poids de l'annexe CRH), b(.) est la base d'estimation de la
# specification consideree, et beta le vecteur de coefficients climatiques du modele.
# Base 2015 et futur 2095 proviennent TOUS DEUX de CMIP6_tas.csv : aucune difference
# ne melange ERA et CMIP6.
#
# Les deux bras vivent dans des bases DIFFERENTES (noeuds sur la distribution
# journaliere vs sur la distribution des moyennes annuelles). Ce n'est pas une gene a
# contourner : c'est le canal FORME rendu concret. Les deux metriques n'estiment pas
# la meme fonction, et Delta_north est le meme fonctionnel applique a chacune.
#
# Les deux formes sont LINEAIRES en beta, donc Var = w' V w est analytique et aucun
# bootstrap n'est requis pour ce test.
#
# INVARIANCE A LA NORMALISATION — MESUREE, PAS SUPPOSEE
# ----------------------------------------------------
# Le pipeline normalise les courbes a max = 0. Delta_north est invariant a toute
# constante additive k des lors que sum_b Delta n_b = 0 par cellule. Cette identite a
# ete MESUREE sur CMIP6_tas.csv (_work_contribA/mesure_A_jourcount.R, 2026-07-31) :
# `count` y est deja une part, sum_b count_b = 1 exactement pour les 17 048 cellules
# x 2 annees ; l'ecart residuel max est 2,11e-15 sur les 7 349 cellules nordiques.
# Borne du biais pour une constante k : |k| * 2,11e-15. L'amplitude de la DRF etant
# d'ordre 1, le residu est a 1e-15. La propriete s'ecrit telle quelle.
# Le design dit la meme chose par un autre chemin : la base sommant a 1 en tout T, la
# direction constante s'annule identiquement dans la difference bilaterale — une
# colonne est structurellement colineaire. Le design NE PEUT PAS voir une constante.
#
# LA COLONNE HORS SUPPORT EST OBLIGATOIRE, ET ELLE COMPTE DANS LES DEUX SENS
# --------------------------------------------------------------------------
# La DRF est ajustee sur [lo, hi] issus de la winsorisation, et clampee plate au-dela
# (convention rule=2 du pipeline aval, damage_projection.r l.123). Les cellules
# nordiques passent une part de l'hiver SOUS lo. Cette masse-jour clampee attenue
# Delta_north^daily MECANIQUEMENT, sans que ce soit une preuve d'indetermination.
# Un intervalle qui couvre zero est ININTERPRETABLE sans cette colonne : un refere
# ecrirait "votre indetermination est un artefact de clamping", et il aurait raison
# de le demander. Masse faible => la revendication est renforcee. Masse importante =>
# on l'apprend AVANT d'ecrire. La colonne d'abord, le verdict ensuite.
#
# PRE-ENGAGEMENT SUR LA NON-DEFINIE-POSITIVITE — FIXE AVANT DE VOIR LES RESULTATS
# -------------------------------------------------------------------------------
# La vcov clusteree bidimensionnelle est de forme Cameron-Gelbach-Miller V1 + V2 - V12.
# Avec 147 x 138 clusters, une valeur propre negative est un cas courant, et w' V w
# peut sortir NEGATIF — auquel cas se = sqrt(.) n'existe pas. Le projet a deja
# rencontre exactement cela sur la vcov de Conley pendant l'inference CRH.
# TRAITEMENT FIXE ICI, AVANT TOUT RESULTAT :
#   (1) on calcule w' V w sur le bloc climatique BRUT de V ; s'il est >= 0, on le publie
#       tel quel et la colonne psd_fix vaut "none" ;
#   (2) s'il est < 0, on tronque a zero les valeurs propres negatives du bloc
#       climatique (projection sur le cone PSD, correction CGM standard), on recalcule,
#       et psd_fix vaut "eigen_trunc" ;
#   (3) la variance brute et la variance corrigee sont TOUTES DEUX reportees
#       (var_raw, var_psd), et le nombre de valeurs propres tronquees aussi.
# Aucune ligne n'est retiree pour cette raison. La decision ne se prend pas devant le
# resultat : elle est ecrite ici.
#
# VALIDATION CROISEE — L'ANALOGUE DU TEST A LAMBDA = 0
# -----------------------------------------------------
# Ce script re-estime, entre autres, les 12 specifications D'ORIGINE (deg <= 4) du
# canal prod. Elles doivent reproduire les dumps publies output/drf_prod_specs/
# drf_prod_spec_*.csv. C'est ce qui distingue une re-implementation d'une
# reconstruction ad hoc : elle se referme sur un artefact independant, produit par un
# autre run, avant cette session. Deux assertions permanentes :
#   A1  ma base reproduit les TERMES D'ESTIMATION de make_climate_terms_wide,
#       c'est-a-dire la base ou vit beta_chapeau — sinon w et beta ne vivent pas
#       dans la meme base et Delta_north = w'beta n'a pas de sens ;
#   A2  mes 12 courbes d'origine reproduisent les 12 dumps publies sur [-25, 40] degC,
#       la comparaison etant faite dans la base d'AFFICHAGE des deux cotes (voir §7).
# A1 est BLOQUANTE. A2 est rapportee et bloquante au-dela de 0.02 log-point.
#
# 🚨 CE SEUIL RENVOYAIT A UN FICHIER QUI N'EXISTE PAS (corrige le 2026-09-19).
# La ligne disait « la tolerance de check_nonregression_specs.r (0.02
# log-point) ». Ce script n'est nulle part dans le depot : l'assertion
# deleguait son seuil a une autorite inexistante, c'est-a-dire qu'elle
# n'avait pas de seuil. La valeur 0.02 est celle qui circulait dans le
# renvoi ; on l'ecrit ici, ou elle est lue, plutot que de pointer vers le
# vide. ⚠ L'autorite reste A ECRIRE : tant qu'elle ne l'est pas, 0.02 est
# une convention heritee et non une tolerance derivee.
#
# ✅ A2 EST REDEVENUE EVALUABLE, ET LE BLOC CI-DESSOUS EST PERIME (2026-09-19).
# Mesure du 19/09 : un rejeu propre de 01a ECRIT output/drf_prod_specs/ --
# « [write] drf_prod_specs/ : 12 courbes sur 12 specifications » -- et les
# douze reproduisent les archives. Le dernier run de ct1 porte d'ailleurs
# a2_status = « evaluee » et a2_max_diff = 5.107e-15 : l'assertion tourne et
# passe. Le paragraphe qui suit decrit l'etat du 2026-09-04, quand 01a
# n'emettait pas encore les douze courbes ; il est conserve pour la trace du
# raisonnement, PAS comme description de l'etat courant.
#
# [PERIME 2026-09-04] A2 n'etait pas evaluable dans le paquet :
# Elle se compare a output/drf_prod_specs/drf_prod_spec_*.csv, les douze dumps
# par specification. AUCUNE etape du paquet ne les produit : 01a emet la courbe
# AIC-ponderee (damage_prod_bin.csv) et les quantiles Monte-Carlo, jamais les
# douze courbes separement. Le §7 le detecte et n'evalue pas ; la colonne
# a2_status de la sortie le dit, et a2_max_diff vaut NA.
# ⚠ Ne PAS verser les dumps comme intrant fourni pour faire passer A2 : ce serait
# comparer la chaine a un artefact sans producteur, c'est-a-dire exactement la
# faute que la migration corrige. Le retablir demande d'emettre les douze courbes
# depuis 01a — une addition a 01a, pas une donnee a verser.
#
# 🚨 LES DEUX ONT ETE RE-ANCREES LE 2026-07-31, ET C'EST LE DEFAUT LE PLUS VICIEUX
# RENCONTRE SUR CE PROJET. Dans leur version d'origine, elles etaient toutes deux
# ancrees sur get_effect_curve — l'AFFICHAGE — alors que beta_chapeau est estime sur
# make_climate_terms_wide — l'ESTIMATION. Elles PASSAIENT, a la precision machine,
# sur le mauvais objet. Une assertion qui passe est un signal de surete, donc
# personne ne la rouvre. Le defaut n'a ete trouve que parce qu'une assertion
# INDEPENDANTE (C2 de tache2_pspline.R) comparait a un artefact tiers et a refuse de
# tourner. Regle qui en sort : une assertion doit s'ancrer sur l'objet qui gouverne
# la revendication, pas sur l'artefact le plus commode.
#
# LES FONCTIONS D'ESTIMATION VIENNENT DE lib_drf.R  (migration 2026-09-04)
# ------------------------------------------------------------------------
# HORS DU PAQUET, ce script extrayait par parsing seize definitions de
# drf_prod_final.r et drf_amenity_final.r, et epinglait le md5 des deux sources.
# Le motif etait exact : ces deux fichiers sont des scripts monolithiques a code
# de premier niveau, qu'on ne peut pas source() sans tout executer. Recopier
# leurs fonctions ici aurait cree une seconde source de verite.
#
# CE MOTIF DISPARAIT DANS LE PAQUET. lib_drf.R est une bibliotheque : elle ne
# porte que des definitions, elle est sourcee par 01a, 01b et 01c, donc exercee a
# chaque run. On l'appelle. Le paquet a par ailleurs UNIFIE les quatre variantes
# trade/migration x daily/annual en deux fonctions generiques, prepare_diffs() et
# run_stage3(), parametrees par is_annual / sign_dest_minus_orig / scale /
# fe_cols / controls. Les huit noms extraits ici (prepare_climate_diffs,
# run_stage3_annual_trade, prepare_climate_diffs_mig, ...) n'existent plus.
#
# 🚨 LE PIEGE DE LA SUBSTITUTION, ET IL EST SILENCIEUX. L'ancienne API portait
# les controles EN VALEUR PAR DEFAUT (drf_prod_final.r l.577 controls =
# gravity_controls ; drf_amenity_final.r l.858 controls = c("GDP_diff",
# "populationdiff")). La generique du paquet a controls = character(0). Les
# omettre estimerait sans controles gravitaires SANS RIEN SIGNALER. Les quatre
# arguments a risque sont donc ecrits explicitement ci-dessous, et RECOPIES DES
# SITES D'APPEL VIVANTS 01a l.541-542 et 01b l.264-270 — pas deduits des
# anciennes signatures.
#
# VALIDATION A BLANC — 2026-09-04, avant migration. Une specification
# (spline4_wl0.5wu0.5), les quatre canal x bras, sur les ANCIENS intrants, sous
# les DEUX ecritures dans la meme session. Ecart neuve - ancienne sur
# Delta_north : -2,3e-15 a +1,1e-15 ; sur se : -2,5e-15 a +8,3e-17 ; bornes de
# support identiques au bit ; nombre d'observations ecartees pour NA identique
# (224 792 prod, 46 818 amen), ce qui atteste que le membre de DROITE est le
# meme et pas seulement le resultat. Le bras temoin reproduisait par ailleurs le
# CSV archive du 31/07 au dernier chiffre, alors que drf_prod_final.r a change
# deux fois depuis : sans lui, la reproduction n'aurait pas ete attribuable.
#
# ⚠ PORTEE DE CETTE VALIDATION : une specification sur 15, et une SPLINE. Les
# NEUF specifications polynomiales empruntent une autre branche (poly() dans
# make_climate_terms_wide, et dcw = NULL dans daily_basis §4). A verifier au
# premier run complet — l'assertion A1 les couvre toutes.
#
# POURQUOI IL N'Y A PLUS DE md5 EN DUR  (migration 2026-09-04)
# ------------------------------------------------------------
# L'epinglage protegeait d'un DANGER PRECIS, qu'il faut nommer avant de le
# remplacer : l'APPARIEMENT CROISE. Ce script tire beta d'une estimation et w
# d'une base reconstruite ; si les deux ne viennent pas du meme etat d'entree,
# Delta_north = w'beta reste arithmetiquement bien forme et ne veut plus rien
# dire. Hors du paquet, les entrees vivaient dans code_jmp/ et a la racine,
# regenerees a la main entre deux sessions : rien ne les tenait ensemble.
#
# DANS LE PAQUET, SIX DES HUIT ENTREES ONT UN PRODUCTEUR VERSIONNE, execute par
# une etape numerotee de run_all.sh, dans un ordre fixe : 01a ecrit
# BACI_raw_pairs_annual.rds et climate_long_annual.rds, 01b step2_migration.rds
# et climate_long.rds, 01d amen_prod.csv ; CMIP6_tas.csv est une sortie serveur
# fournie et verifiee en tete de run_all.sh. La mutabilite contre laquelle le
# md5 protegeait n'existe plus.
#
# ⛔ ET DES md5 GELES SERAIENT PIRES QU'INUTILES ICI. Le script le disait
# lui-meme hors du paquet : ses empreintes portaient sur les copies RACINE de
# drf_*.r, qui ne sont plus executees depuis le 24/08 — « ce ne sont pas des
# garde-fous VIVANTS, c'est un SCEAU D'ARCHIVE ». Un sceau d'archive dans un
# paquet de replication est un garde-fou qu'un refere apprend a ignorer : tout
# rejeu legitime le fait echouer.
#
# CE QUI REMPLACE : une EMPREINTE DE PROVENANCE (§0). Elle est calculee sur les
# six entrees, inscrite en colonne de la sortie, et re-verifiee par ct2 et ct3.
# Le refus ne porte plus sur « les md5 valent ceci » mais sur « ces artefacts ne
# viennent pas du meme etat d'entree » — c'est-a-dire sur le danger vise. Meme
# idiome qu'INPUTS_MD5 dans la grille GE de run_all.sh.
#
# 🔑 ET A1 RESTE BLOQUANTE. C'est elle qui attrape reellement l'appariement
# base / beta, elle est interne, elle ne depend d'aucun md5, et elle couvre les
# les specifications de GRID_SPECS a chaque run.
#
# SORTIE : output/contribA/test1_delta_north.csv
#
# Romain Fillon / Contribution A, Test 1 — 2026-07-31
# Migre au paquet de replication — 2026-09-04
##################################################################################

suppressMessages({
  library(data.table); library(splines); library(fixest); library(tools)
})
options(width = 200)
t0 <- proc.time(); cat("=== ct1_delta_north_grid.R ===\n")

## Racine deduite de l'emplacement de CE fichier, jamais du repertoire courant.
## Meme idiome que 01a l.46-59 ; JMP_ROOT l'emporte quand il est pose.
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
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "04_contrast", "lib_contrast.R"))

DIR_CONTRIB <- file.path(DIR_OUT, "contribA")
if (!dir.exists(DIR_CONTRIB)) dir.create(DIR_CONTRIB, recursive = TRUE)

BASE_YEAR <- 2015L; FUT_YEAR <- 2095L; NORTH_LAT <- 45
POP_YEAR  <- 2005L                     # poids population de l'annexe CRH

##################################################################################
# 0. LES ENTREES, ET LEUR EMPREINTE DE PROVENANCE
#
# Six entrees, six producteurs, tous dans le paquet sauf la sortie serveur.
# Voir l'en-tete pour pourquoi ce n'est plus un epinglage a des md5 geles.
#
# L'EMPREINTE agrege les six md5 en un jeton unique, inscrit en colonne de la
# sortie. ct2 et ct3 la recalculent sur LEURS entrees et refusent de tourner si
# elle differe. C'est ce qui interdit l'appariement croise que l'ancien
# epinglage visait : une decomposition assemblee de millesimes melanges
# telescope exactement, et aucun controle de somme ne peut donc la voir.
#
# ⚠ CE QUE L'EMPREINTE NE VOIT PAS : deux artefacts du MEME run mais estimes
# sous des SPECIFICATIONS differentes. C'est un danger reel, distinct, et il se
# traite en amont -- en ecrivant les controles a chaque site d'appel plutot
# qu'en s'en remettant a un defaut. Voir §1.
##################################################################################

PATHS <- c(
  "BACI_raw_pairs_annual.rds" = file.path(DIR_TMP, "BACI_raw_pairs_annual.rds"),
  "climate_long_annual.rds"   = file.path(DIR_TMP, "climate_long_annual.rds"),
  "step2_migration.rds"       = file.path(DIR_TMP, "step2_migration.rds"),
  "climate_long.rds"          = file.path(DIR_TMP, "climate_long.rds"),
  "CMIP6_tas.csv"             = file.path(DIR_TMP, "CMIP6_tas.csv"),
  "amen_prod.csv"             = file.path(DIR_CRH, "amen_prod.csv")
)
PRODUCTEUR <- c(
  "BACI_raw_pairs_annual.rds" = "01a_productivity_drf.R  l.274",
  "climate_long_annual.rds"   = "01a_productivity_drf.R  l.284",
  "step2_migration.rds"       = "01b_amenity_drf.R       l.141",
  "climate_long.rds"          = "01b_amenity_drf.R       l.149",
  "CMIP6_tas.csv"             = "00_temperature_projections_JMP.R  [SERVEUR, fourni]",
  "amen_prod.csv"             = "01d_crh_triangulation.R l.332"
)


cat("\n--- entrees et provenance ---\n")
manque <- names(PATHS)[!file.exists(PATHS)]
if (length(manque)) {
  for (n in manque)
    cat(sprintf("  ABSENT   %-28s produit par : %s\n", n, PRODUCTEUR[n]))
  stop("[STOP] ", length(manque), " entree(s) absente(s).\n",
       "  Lancer d'abord l'etape qui les produit — voir run_all.sh.\n",
       "  Une entree absente n'est pas un incident de chemin : c'est une etape\n",
       "  qui n'a pas tourne, et la grille n'a rien a estimer.")
}
OBSERVED_MD5 <- setNames(as.character(tools::md5sum(PATHS)), names(PATHS))
for (n in names(PATHS))
  cat(sprintf("  ok  %-28s %s\n      <- %s\n", n, OBSERVED_MD5[n], PRODUCTEUR[n]))
INPUTS_FP <- empreinte(OBSERVED_MD5)
cat(sprintf("\n  empreinte de provenance : %s\n", INPUTS_FP))

##################################################################################
# 1. LES FONCTIONS D'ESTIMATION — lib_drf.R, ET LES ARGUMENTS QUI LES SPECIALISENT
#
# prepare_diffs() et run_stage3() sont generiques. Ce qui distingue les quatre
# combinaisons canal x bras tient dans CINQ arguments, ecrits ici explicitement.
#
# ⛔ AUCUN DEFAUT N'EST ACCEPTE, ET C'EST LA REGLE DU BLOC. run_stage3() a
# controls = character(0) ; l'ancienne API avait controls = gravity_controls
# (drf_prod_final.r l.577). Un site qui s'en remet au defaut estime sans
# controles gravitaires sans rien signaler -- mesure le 2026-09-04 :
# Delta_north +22,9 %, ecart d'AIC 137 502. C'est exactement le defaut qui
# dormait dans 01a l.542 et l.402, corrige le meme jour.
#
# D'OU VIENNENT CES VALEURS, LIGNE A LIGNE :
#   amen : les quatre sont recopies du site vivant 01b l.264-270.
#   prod : fe / sign sont les defauts de prepare_diffs() ; `scale` n'en a
#          plus depuis le 20/09/2026 et est passe explicitement, qui sont ce
#          que 01a l.541 utilise. controls vient du defaut de l'ANCIENNE
#          run_stage3 (drf_prod_final.r l.577, applique aux deux grilles
#          l.1606 et l.2068) -- et non de 01a l.542, qui ne passait rien
#          jusqu'au correctif du 2026-09-04. C'est la valeur que 01c passe
#          explicitement a ses huit sites, dont l.672 qui est la meme grille,
#          et c'est celle qui reproduit le CSV archive du 31/07 au dernier
#          chiffre (validation a blanc du 2026-09-04).
##################################################################################

ARGS <- list(
  prod = list(
    controls = TRADE_STAGE3_CONTROLS,            # = c(GRAVITY_CONTROLS, "GDP_diff")
    fe       = "fe_prod_year",                   # defaut de prepare_diffs / run_stage3
    sign_dmo = FALSE,                            # commerce : orig (exportateur) - dest
    ## 19/09/2026 : etait `scale = 1`, commente « aucune mise a l'echelle ».
    ## Le PPML identifie theta * f_z (eq. 5) et le GE reapplique theta dans
    ## trade_weights : sans cette division le choc etait amplifie deux fois.
    ## Le canal amenite faisait deja le geste, avec Omega, deux lignes plus bas.
    scale    = SCALE_PROD),                               # aucune mise a l'echelle
  amen = list(
    controls = AMENITY_STAGE3_CONTROLS,          # 01b l.269
    fe       = c("fe_orig_dest", "fe_year"),     # 01b l.266 et l.270
    sign_dmo = TRUE,                             # migration : dest - orig
    scale    = Omega)                            # 01b l.267
)

## Garde-fou d'ecriture : les cinq arguments doivent etre non vides la ou ils
## portent quelque chose. Un `controls` vide cote prod serait precisement le
## defaut qu'on vient de bannir, et il repasserait inapercu.
stopifnot(length(ARGS$prod$controls) == 6L, length(ARGS$amen$controls) == 2L,
          isFALSE(ARGS$prod$sign_dmo), isTRUE(ARGS$amen$sign_dmo),
          ARGS$prod$scale == SCALE_PROD, ARGS$amen$scale == SCALE_AMEN,
          ## aucune liste ne peut omettre son echelle : un canal ajoute
          ## demain heriterait de NULL et scalerait en silence.
          all(vapply(ARGS, function(z) is.numeric(z$scale) &&
                       length(z$scale) == 1L && is.finite(z$scale), TRUE)))
cat("\n--- arguments d'estimation (explicites, aucun defaut) ---\n")
for (ch in names(ARGS))
  cat(sprintf("  %-5s controls=%s | fe=%s | sign_dest_minus_orig=%s | scale=%g\n",
              ch, paste(ARGS[[ch]]$controls, collapse = "+"),
              paste(ARGS[[ch]]$fe, collapse = "+"),
              ARGS[[ch]]$sign_dmo, ARGS[[ch]]$scale))

##################################################################################
# 2. LA GRILLE — LUE DANS 00_parameters.R, ARRETEE AVANT TOUT RESULTAT
##################################################################################
# Aucune specification n'est retiree apres avoir vu son resultat, aucune ajoutee
# pour en rattraper un. La grille elle-meme est une DECISION D'AUTEUR, prise hors
# de ce script et datee : le 2026-09-18 elle passe de vingt-quatre a QUINZE
# specifications -- polynomes de degre 2 a 4, splines naturelles a 3 et 4 df,
# trois winsorisations. Les splines a 5, 6 et 7 df sont retirees parce qu'elles
# oscillent dans la partie DENSE du support (20-28 degC).
#
# 🚨 CE BLOC PORTAIT SA PROPRE LISTE JUSQU'AU 2026-09-18, ET C'ETAIT LA FAUTE.
# Il ecrasait deliberement SPEC_CONFIGS de 00_parameters.R, avec un commentaire
# expliquant pourquoi l'ecrasement etait obligatoire. Le motif etait juste -- la
# grille du CORPS (12) n'est pas celle du CONTRASTE -- mais la reponse ne l'etait
# pas : ct2 portait la MEME liste, recopiee, et ct4 comme ct6 ne comptaient les
# specifications que sur le CSV recu. Changer la grille demandait donc deux
# editions coordonnees dans deux fichiers, et rien n'aurait signale la seconde
# oubliee : les deux scripts auraient tourne, chacun sur sa grille, et la
# jointure de ct4 se serait arretee sur un message d'etiquettes -- au mieux.
# ⇒ La liste vit desormais dans 00_parameters.R sous le nom GRID_SPECS, et ce
#   script la LIT. tools/check_macros.py (controle 9) echoue si un producteur
#   en reintroduit une en dur.
#
# ⚠ NEUF des 15 specifications sont POLYNOMIALES (poly 2, 3, 4 x 3 winsorisations).
# Elles empruntent une autre branche que les splines, dans make_climate_terms_wide
# comme dans daily_basis (§4, ou dcw vaut NULL). La validation a blanc du 2026-09-04
# n'a exerce qu'une spline : c'est A1 qui couvre les neuf autres, a chaque run.

SPEC_CONFIGS   <- GRID_FORMS
WINSOR_CONFIGS <- GRID_WINSOR

## L'IDENTITE AVEC GRID_SPECS EST VERIFIEE, PAS SUPPOSEE. Les deux boucles
## imbriquees ci-dessous reconstruisent le produit croise ; GRID_SPECS le porte
## deja, a plat. Si les deux divergeaient -- une forme ajoutee ici, un ordre
## change la-bas -- ce script estimerait une grille et les comptes en
## decriraient une autre. On compare les ETIQUETTES, c'est-a-dire l'objet meme
## qui joint ct1, ct2 et ct4.
.labels_locaux <- unlist(lapply(SPEC_CONFIGS, function(cfg)
  vapply(WINSOR_CONFIGS, function(wc)
    sprintf("%s%d_wl%gwu%g", cfg$spec, cfg$degree, wc$lower, wc$upper),
    character(1))))
if (!setequal(.labels_locaux, GRID_LABELS) ||
    length(.labels_locaux) != GRID_NSPEC)
  stop("[FATAL] la grille parcourue par ce script n'est pas GRID_SPECS.\n",
       "  parcourue : ", paste(sort(.labels_locaux), collapse = ", "), "\n",
       "  GRID_SPECS: ", paste(sort(GRID_LABELS),   collapse = ", "))
stopifnot(length(WINSOR_CONFIGS) == 3L)

# fige les 3 winsorisations du papier AVANT toute reduction (mode smoke) : A2 en a
# besoin pour reconstruire la grille commune des dumps, qui ne depend que d'elles.
WINSOR_CONFIGS_ORIG <- WINSOR_CONFIGS
# Mode smoke : valide la plomberie (extraction, geographie, A1, w, PSD) sur une seule
# specification par canal x bras. La sortie va dans un AUTRE fichier : une table smoke
# ne doit jamais pouvoir etre prise pour la table du test.
SMOKE <- nzchar(Sys.getenv("T1_SMOKE"))
if (SMOKE) {
  ## 🚨 CHOIX PAR NOM, ET JAMAIS PAR INDICE. Cette ligne a porte
  ## `SPEC_CONFIGS[3]` avec le commentaire « spline df=3 » : c'etait vrai de
  ## l'ancienne liste, ou les quatre formes d'origine venaient en tete. Dans
  ## GRID_FORMS, le rang 3 est le POLYNOME de degre 4. Un indice ne dit pas ce
  ## qu'il designe, et le smoke aurait valide la plomberie sur une autre forme
  ## sans que rien ne l'indique -- le compte « 1 x 1 » etant reste juste.
  .i <- which(vapply(SPEC_CONFIGS,
                     function(z) z$spec == "spline" && z$degree == 3L, TRUE))
  if (length(.i) != 1L)
    stop("[FATAL] mode smoke : la spline df=3 n'est pas dans GRID_FORMS, ou y ",
         "figure ", length(.i), " fois.")
  SPEC_CONFIGS   <- SPEC_CONFIGS[.i]    # spline df=3, specification d'origine
  WINSOR_CONFIGS <- WINSOR_CONFIGS[2]   # [1%, 1%]
  cat("\n*** MODE SMOKE : grille reduite a 1 x 1. Sortie NON publiable. ***\n")
}
OUT_CSV <- file.path(DIR_CONTRIB,
                     if (SMOKE) "test1_SMOKE.csv" else "test1_delta_north.csv")

cat(sprintf("\nGrille : %d formes x %d winsorisations = %d modeles par (canal, bras)\n",
            length(SPEC_CONFIGS), length(WINSOR_CONFIGS),
            length(SPEC_CONFIGS) * length(WINSOR_CONFIGS)))

##################################################################################
# 3. GEOGRAPHIE NORDIQUE : cellules, poids population, distributions base/futur
##################################################################################

cat("\n--- geographie nordique ---\n")
cm <- fread(PATHS["CMIP6_tas.csv"])[year %in% c(BASE_YEAR, FUT_YEAR) & y >= NORTH_LAT]
cm[, c := count / sum(count), by = .(x, y, year)]
cat(sprintf("  cellules y >= %g : %d\n", NORTH_LAT, uniqueN(cm[, .(x, y)])))

pop <- fread(PATHS["amen_prod.csv"])[year == POP_YEAR, .(lat, lon, pop = as.numeric(pop_x))]
cm <- merge(cm, pop, by.x = c("x", "y"), by.y = c("lon", "lat"), all = FALSE)
ncell <- uniqueN(cm[, .(x, y)])
cat(sprintf("  cellules appariees aux poids population %d : %d\n", POP_YEAR, ncell))
# une cellule ne compte que si elle a SES DEUX annees
ok2 <- cm[, .(ny = uniqueN(year)), by = .(x, y)][ny == 2L, .(x, y)]
cm <- merge(cm, ok2, by = c("x", "y"))
ncell <- nrow(ok2)
cat(sprintf("  cellules avec base ET futur : %d\n", ncell))

# poids population normalises (constants dans le temps)
wcell <- unique(cm[, .(x, y, pop)])
wcell[, omega := pop / sum(pop)]
stopifnot(abs(sum(wcell$omega) - 1) < 1e-12)

# les poids entrent une fois pour toutes
cm <- merge(cm, wcell[, .(x, y, omega)], by = c("x", "y"))

# DAILY : masse-jour signee par (cellule, temperature) : Delta c_b(cellule),
# puis agregee sur le Nord. C'est le seul objet dont Delta_north^daily depend
# (linearite : Delta = (sum_cellules omega_c Delta c_b)' B beta).
cm[, sgn := fifelse(year == FUT_YEAR, 1, -1)]
dmass <- cm[, .(m = sum(omega * sgn * c)), by = value][order(value)]
cat(sprintf("  temperatures entieres presentes : [%d, %d]\n", min(dmass$value), max(dmass$value)))
cat(sprintf("  somme des masses signees (doit etre ~0) : %.3e\n", sum(dmass$m)))

# part de masse-jour nordique par temperature, base et futur separement :
# c'est la source de la colonne hors support
share_by_T <- cm[, .(s = sum(omega * c)), by = .(year, value)]
stopifnot(all(abs(share_by_T[, .(t = sum(s)), by = year]$t - 1) < 1e-9))

# ANNUAL : moyenne annuelle par cellule, base et futur
mu <- cm[, .(mu = sum(value * c)), by = .(x, y, year)]
muw <- dcast(mu, x + y ~ year, value.var = "mu")
setnames(muw, c(as.character(BASE_YEAR), as.character(FUT_YEAR)), c("mu_base", "mu_fut"))
muw <- merge(muw, wcell[, .(x, y, omega)], by = c("x", "y"))
muw <- muw[is.finite(mu_base) & is.finite(mu_fut)]
cat(sprintf("  Tbar nordique : base [%.1f, %.1f]  futur [%.1f, %.1f] degC\n",
            min(muw$mu_base), max(muw$mu_base), min(muw$mu_fut), max(muw$mu_fut)))
rm(cm); gc()

##################################################################################
# 4. OUTILS : base d'estimation, Delta_north, variance PSD
#
# daily_basis / annual_basis / delta_var vivent dans lib_contrast.R, sourcee en
# tete : ct2 les emploie aussi, et une definition dupliquee divergerait en
# silence. C'est la meme raison qui met prepare_diffs dans lib_drf.R.
##################################################################################

##################################################################################
# 5. BOUCLE PRINCIPALE
##################################################################################

## Le dictionnaire ne porte plus de NOMS DE FONCTIONS : les quatre variantes ont
## disparu au profit d'une generique. Il ne reste que les donnees du canal et les
## motifs de colonnes ; les arguments d'estimation sont dans ARGS (§1).
CHANNELS <- list(
  prod = list(
    data_path = PATHS["BACI_raw_pairs_annual.rds"],
    clim_path = PATHS["climate_long_annual.rds"], year_min = NA_integer_,
    pat_daily = "^diff_%s_[1-9]", pat_annual = "^diff_annual_%s_[1-9]"),
  amen = list(
    data_path = PATHS["step2_migration.rds"],
    clim_path = PATHS["climate_long.rds"], year_min = YEAR_MIN_STAGE3,
    pat_daily = "^diff_%s_[1-9]", pat_annual = "^diff_annual_%s_[1-9]")
)

RES <- list(); CURVES <- list(); k <- 0L
## CURVES : les douze specs d'origine, pour l'assertion A2 (ne pas elargir,
## A2 itere sur ses noms). CURVES_ALL : TOUTE la grille du bras prod/daily,
## pour l'archivage. Deux usages, deux accumulateurs.
CURVES_ALL <- list()
## CURVES_ANN : toute la grille du bras prod/ANNUEL. Meme geste, autre bras.
CURVES_ANN <- list()

for (chan in names(CHANNELS)) {
  CH <- CHANNELS[[chan]]; AG <- ARGS[[chan]]
  cat(sprintf("\n################ CANAL %s ################\n", toupper(chan)))
  dat  <- readRDS(CH$data_path); setDT(dat)
  if (!is.na(CH$year_min)) dat <- dat[year >= CH$year_min]
  clim <- readRDS(CH$clim_path); setDT(clim)
  cat(sprintf("  donnees : %d lignes | climat : %d lignes\n", nrow(dat), nrow(clim)))

  for (arm in c("daily", "annual")) {
    for (cfg in SPEC_CONFIGS) {
      spec <- cfg$spec; deg <- cfg$degree
      for (wc in WINSOR_CONFIGS) {
        wl <- wc$lower; wu <- wc$upper
        key <- sprintf("%s%d_wl%gwu%g", spec, deg, wl, wu)
        k <- k + 1L
        cat(sprintf("\n[%02d] %s / %s / %s ... ", k, chan, arm, key)); flush.console()
        tk <- proc.time()

        ## Un seul chemin d'appel pour les quatre combinaisons : `is_annual` porte
        ## le bras, ARGS porte le canal. Les cinq arguments sont passes
        ## explicitement — voir §1 pour pourquoi aucun defaut n'est accepte ici.
        fit <- tryCatch({
          is_ann <- (arm == "annual")
          dp <- prepare_diffs(dat, clim, spec, deg, wl, wu, is_ann,
                              fe_cols              = AG$fe,
                              sign_dest_minus_orig = AG$sign_dmo,
                              scale                = AG$scale)
          m  <- run_stage3(dp, spec, deg, is_ann,
                           controls = AG$controls,
                           fe_str   = AG$fe)
          rm(dp); gc(verbose = FALSE); m
        }, error = function(e) { cat("ECHEC:", conditionMessage(e), " "); NULL })

        if (is.null(fit)) {
          RES[[length(RES)+1L]] <- data.table(
            channel = chan, arm = arm, spec_label = key, form = spec, df = deg,
            winsor_lo = wl, winsor_hi = wu, origine = cfg$origine,
            status = "FIT_FAILED")
          next
        }

        ## `gec` EST REMIS A NULL A CHAQUE ITERATION (2026-09-17).
        ## Il n'etait affecte que dans la branche journaliere : sur une
        ## iteration annuelle il gardait silencieusement la courbe de la
        ## DERNIERE specification journaliere. Rien ne le lisait jusqu'ici --
        ## la branche annuelle pose a1_disp <- NA -- de sorte qu'aucune sortie
        ## publiee n'en descend ; mais le premier lecteur ajoute en a herite
        ## autant de copies d'une courbe journaliere etiquetees « annuel ».
        ## Une variable de boucle qui survit a son bloc est un piege, pas une
        ## optimisation : elle est desormais effacee avant chaque usage.
        gec <- NULL

        mod <- fit$model
        pat <- sprintf(if (arm == "daily") CH$pat_daily else CH$pat_annual, spec)
        cc  <- grep(pat, names(coef(mod)), value = TRUE)
        beta <- coef(mod)[cc]
        V    <- vcov(mod)[cc, cc, drop = FALSE]

        # ---- base d'estimation + poids w ----
        if (arm == "daily") {
          B  <- daily_basis(clim, spec, deg, wl, wu)
          # ---- A1, RE-ANCREE LE 2026-07-31 ------------------------------------
          # Version anterieure : « ma base doit reproduire get_effect_curve ».
          # Elle n'etait PAS fausse — elle passait a la precision machine — mais
          # elle etait ancree sur le MAUVAIS OBJET. get_effect_curve est
          # l'AFFICHAGE ; les coefficients, eux, vivent dans la base de
          # make_climate_terms_wide, l'ESTIMATION. Une assertion qui PASSE est un
          # signal de surete, donc personne ne la rouvre : c'est le mode de
          # defaillance le plus vicieux rencontre sur ce projet, et il n'a ete
          # trouve que parce qu'une assertion INDEPENDANTE (C2 de la Tache 2) a
          # compare a un artefact tiers et a refuse de tourner.
          # A1 s'ancre desormais sur l'objet qui gouverne la revendication : les
          # termes d'estimation eux-memes.
          twref <- make_climate_terms_wide(clim, spec, deg, wl, wu)
          dcw <- B$dcw
          if (is.null(dcw)) {          # polynomes : memes poids, base insensible
            dcw <- copy(as.data.table(clim))
            dcw[value < B$lo, value := B$lo]; dcw[value > B$hi, value := B$hi]
            dcw <- dcw[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
            dcw[, count := 365 * count / sum(count), by = .(GID_0, year)]
          }
          tmine <- dcw[, {
            bv <- predict(B$basis, value); wt <- count / 365
            as.list(setNames(colSums(bv * wt, na.rm = TRUE), paste0("t", seq_len(ncol(bv)))))
          }, by = .(GID_0, year)]
          setkey(twref, GID_0, year); setkey(tmine, GID_0, year)
          pref <- if (spec == "polynomial") "polynomial_" else "spline_"
          cr <- paste0(pref, seq_len(deg)); cmn <- paste0("t", seq_len(deg))
          stopifnot(identical(twref$GID_0, tmine$GID_0), identical(twref$year, tmine$year))
          a1 <- max(abs(as.matrix(twref[, ..cr]) - as.matrix(tmine[, ..cmn])))
          if (!is.finite(a1) || a1 > 1e-12)
            stop(sprintf(paste0("A1 VIOLEE (%s/%s/%s) : ma base ne reproduit pas les ",
                                "TERMES D'ESTIMATION (%.3e). w et beta ne vivraient pas ",
                                "dans la meme base."), chan, arm, key, a1))
          # --- diagnostic NON bloquant : l'ecart affichage / estimation ---------
          # Conserve plutot que supprime : c'est la mesure du §5, et elle est
          # informative en elle-meme (elle chiffre de combien la courbe PUBLIEE
          # s'ecarte de la fonction ESTIMEE).
          fit$wval_lower <- wl; fit$wval_upper <- wu
          gec <- tryCatch(get_effect_curve(fit, clim), error = function(e) NULL)
          a1_disp <- NA_real_
          if (!is.null(gec)) {
            yg <- as.numeric(predict(B$basis, gec$x) %*% beta)
            yg <- yg - max(yg, na.rm = TRUE)
            a1_disp <- max(abs(yg - gec$y), na.rm = TRUE)
          }
          Tcl <- pmin(pmax(dmass$value, B$lo), B$hi)          # clamp plat, rule=2
          X   <- predict(B$basis, Tcl)
          w   <- as.numeric(t(X) %*% dmass$m)
          oos <- sum(share_by_T[year == FUT_YEAR & (value < B$lo | value > B$hi), s])
          oos_base <- sum(share_by_T[year == BASE_YEAR & (value < B$lo | value > B$hi), s])
        } else {
          B  <- annual_basis(clim, spec, deg, wl, wu, mean_climate)
          # bras annuel : annual_basis est DEJA identique a make_annual_terms_trade
          # (memes bornes min/max, memes quantiles, meme x_ref) — verifie par lecture
          # le 31/07. C'est le bras journalier qui etait desaligne, pas celui-ci.
          a1 <- NA_real_; a1_disp <- NA_real_

          ## LA COURBE DU BRAS ANNUEL (2026-09-17).
          ## `is_annual = TRUE` est OBLIGATOIRE : sans lui get_effect_curve
          ## cherche les colonnes ^diff_<spec>_[1-9] et ne trouve rien sur un
          ## modele annuel, dont les coefficients s'appellent
          ## diff_annual_<spec>_[1-9] ; il rendrait NULL.
          ## a1_disp reste NA : c'est un diagnostic du bras journalier, et le
          ## remplir ici changerait une colonne publiee. L'ecart affichage /
          ## estimation du bras annuel est MESURE quand meme, et ecrit au log.
          fit$wval_lower <- wl; fit$wval_upper <- wu
          gec <- tryCatch(get_effect_curve(fit, clim, is_annual = TRUE),
                          error = function(e) NULL)
          if (!is.null(gec)) {
            yg <- as.numeric(predict(B$basis, gec$x) %*% beta)
            yg <- yg - max(yg, na.rm = TRUE)
            cat(sprintf("[gec-ann ecart affichage/estimation = %.2e] ",
                        max(abs(yg - gec$y), na.rm = TRUE)))
          }

          Xf <- predict(B$basis, pmin(pmax(muw$mu_fut,  B$lo), B$hi))
          Xb <- predict(B$basis, pmin(pmax(muw$mu_base, B$lo), B$hi))
          w  <- as.numeric(t(Xf - Xb) %*% muw$omega)
          oos      <- sum(muw$omega[muw$mu_fut  < B$lo | muw$mu_fut  > B$hi])
          oos_base <- sum(muw$omega[muw$mu_base < B$lo | muw$mu_base > B$hi])
        }
        stopifnot(length(w) == length(beta))

        dn <- as.numeric(crossprod(w, beta))
        vv <- delta_var(w, V)
        se <- sqrt(vv$var_psd)

        RES[[length(RES)+1L]] <- data.table(
          channel = chan, arm = arm, spec_label = key, form = spec, df = deg,
          winsor_lo = wl, winsor_hi = wu, origine = cfg$origine, status = "OK",
          delta_north = dn, se = se,
          ci_lo = dn - 1.96 * se, ci_hi = dn + 1.96 * se,
          excludes_zero = (dn - 1.96 * se) * (dn + 1.96 * se) > 0,
          var_raw = vv$var_raw, var_psd = vv$var_psd,
          psd_fix = vv$fix, n_eigen_trunc = vv$ntrunc,
          oos_mass_fut = oos, oos_mass_base = oos_base,
          support_lo = B$lo, support_hi = B$hi,
          aic = AIC(mod), n_obs = nobs(mod), n_coef_clim = length(beta),
          n_cells_north = ncell,
          cluster_unit = "orig x dest (two-way CGM)",
          p_joint = fit$wald_joint$p, p_nonlin = fit$wald_nonlin$p,
          a1_maxdiff = a1, a1_display_gap = a1_disp,
          basis_convention = "estimation (make_climate_terms_wide)",
          secs = as.numeric((proc.time()-tk)["elapsed"]))

        # courbes, pour l'assertion A2 (canal prod, bras daily, specs d'origine)
        if (chan == "prod" && arm == "daily" && cfg$origine) {
          # ---- A2, RE-ANCREE LE 2026-07-31, meme defaut qu'A1 ----------------
          # Les dumps publies sont produits par save_damage_csv a partir de
          # get_effect_curve, donc dans la base d'AFFICHAGE. Depuis que
          # daily_basis suit la base d'ESTIMATION, comparer MA courbe a ces dumps
          # ne teste plus une non-regression : ca mesure la difference de
          # conventions — deja mesuree ailleurs (tache2_diag_courbe.csv), et sur
          # une grille plus fine.
          # Ce qu'A2 doit tester est que MON beta_chapeau reproduit celui du run
          # publie. On compare donc la courbe d'AFFICHAGE construite sur MON beta
          # — c'est exactement `gec` — aux dumps. Si les coefficients coincident,
          # cette comparaison-la passe ; si elle echoue, c'est l'estimation qui a
          # derive, ce qu'A2 est faite pour attraper.
          if (!is.null(gec))
            CURVES[[key]] <- data.table(temperature_C = round(gec$x, 1),
                                        y = gec$y, model_name = key)
        }

        ## ARCHIVAGE : TOUTE LA GRILLE, pas seulement les specs d'origine.
        ## `gec` est calcule pour toutes (a1_display_gap est rempli partout) ;
        ## il etait jete pour les douze autres, qui sont justement celles
        ## qu'aucune table ne pouvait citer.
        if (chan == "prod" && arm == "daily" && !is.null(gec))
          CURVES_ALL[[key]] <- data.table(temperature_C = round(gec$x, 1),
                                          y = gec$y, model_name = key)

        ## LE MEME GESTE POUR LE BRAS ANNUEL (2026-09-17).
        ## `gec` est calcule pour les DEUX bras, quelques lignes plus haut, et
        ## n'etait garde que pour daily. Le panneau (b) de la Figure 1 n'avait
        ## donc aucune source de courbes par specification, alors que le panneau
        ## (a) portait toute la grille : une asymetrie d'ARCHIVAGE, pas de calcul.
        ## Rien n'est reestime ici non plus ; on cesse de jeter, des deux cotes.
        if (chan == "prod" && arm == "annual" && !is.null(gec))
          CURVES_ANN[[key]] <- data.table(temperature_C = round(gec$x, 1),
                                          y = gec$y, model_name = key)

        cat(sprintf("Dn=%+.4f se=%.4f oos=%.3f%% %s (%.0fs)",
                    dn, se, 100*oos, vv$fix, (proc.time()-tk)["elapsed"]))
        rm(fit, mod); gc(verbose = FALSE)
      }
    }
  }
  rm(dat, clim); gc(verbose = FALSE)
}

##################################################################################
# 6. POIDS AIC SUR LA GRILLE ETENDUE (par canal x bras)
##################################################################################

R <- rbindlist(RES, fill = TRUE)
R[status == "OK", delta_aic := aic - min(aic, na.rm = TRUE), by = .(channel, arm)]
R[status == "OK", aic_weight := exp(-delta_aic/2) / sum(exp(-delta_aic/2), na.rm = TRUE),
  by = .(channel, arm)]

##################################################################################
# 7. ASSERTION A2 — non-regression contre les dumps publies
##################################################################################

# PORTEE DE LA COMPARAISON — corrige apres le smoke du 2026-07-31.
# Le dump N'EST PAS la courbe brute de get_effect_curve : save_damage_csv l.2420-2423
# interpole chaque spec sur la GRILLE COMMUNE = intersection des supports des 12 specs,
# puis (l.2456-2458) ECRETE PLAT au-dela. La grille commune est donc pilotee par la
# winsorisation la plus agressive (2,5 %). Comparer sur le support propre d'une spec a
# 1 % revient a confronter une courbe qui varie a une courbe constante sur les ailes :
# le premier smoke a ainsi renvoye 0,72 sur du code pourtant exact (AIC et p-values
# reproduits au chiffre pres). A2 ne vaut donc que SUR LA GRILLE COMMUNE, qu'on
# recalcule ici a partir des seuils de winsorisation eux-memes — pas devinee.
cat("\n\n--- A2 : non-regression contre output/drf_prod_specs/ ---\n")
a2_max <- NA_real_
dumpdir <- file.path(DIR_OUT, "drf_prod_specs")
if (dir.exists(dumpdir) && length(CURVES)) {
  clim_p <- readRDS(PATHS["climate_long_annual.rds"]); setDT(clim_p)
  ## LA FORME ET LE DEGRE NE SERVENT A RIEN ICI, ET C'EST MESURE (18/09/2026).
  ## Seules B$lo et B$hi sont lues, et winsor_bounds ne depend que des deux
  ## seuils. La ligne portait pourtant « "spline", 3L » en clair : une
  ## specification ecrite en dur dans un producteur, c'est-a-dire ce que la
  ## regle du 18/09 interdit -- et, si les bornes en avaient dependu, une grille
  ## commune silencieusement fausse. On prend donc la PREMIERE forme de la
  ## grille, et on VERIFIE l'independance sur la DERNIERE : sans ce controle,
  ## « ne depend pas de la forme » resterait une affirmation.
  .f1 <- GRID_FORMS[[1]]; .f2 <- GRID_FORMS[[length(GRID_FORMS)]]
  sup <- vapply(WINSOR_CONFIGS_ORIG, function(wc) {
    B1 <- daily_basis(clim_p, .f1$spec, .f1$degree, wc$lower, wc$upper)
    B2 <- daily_basis(clim_p, .f2$spec, .f2$degree, wc$lower, wc$upper)
    if (!isTRUE(all.equal(c(B1$lo, B1$hi), c(B2$lo, B2$hi))))
      stop("[FATAL] les bornes de support dependent de la forme ",
           "(", .f1$spec, .f1$degree, " vs ", .f2$spec, .f2$degree, ") : ",
           "la grille commune de A2 ne serait pas definie.")
    c(B1$lo, B1$hi)
  }, numeric(2))
  common_lo <- max(sup[1, ]); common_hi <- min(sup[2, ])
  rm(clim_p)
  cat(sprintf("  grille commune des 12 dumps : [%.1f, %.1f] degC\n", common_lo, common_hi))
  pub <- rbindlist(lapply(list.files(dumpdir, pattern = "^drf_prod_spec_\\d+\\.csv$",
                                     full.names = TRUE), fread))
  diffs <- c()
  for (nm in names(CURVES)) {
    p <- pub[model_name == nm, .(temperature_C, f_z_raw)]
    if (!nrow(p)) { cat(sprintf("  %-24s absent des dumps\n", nm)); next }
    mine <- CURVES[[nm]]
    cmp <- merge(p, mine[, .(temperature_C = round(temperature_C, 1), y)],
                 by = "temperature_C")
    cmp <- cmp[temperature_C >= max(-25, ceiling(common_lo)) &
               temperature_C <= min(40, floor(common_hi))]
    if (!nrow(cmp)) { cat(sprintf("  %-24s pas de recouvrement\n", nm)); next }
    # les dumps publies sont normalises max=0 sur leur propre grille : on recale
    d <- max(abs((cmp$f_z_raw - max(cmp$f_z_raw)) - (cmp$y - max(cmp$y))))
    diffs <- c(diffs, d)
    cat(sprintf("  %-24s max|diff| = %.5f  %s\n", nm, d,
                if (d <= 0.02) "OK" else "*** HORS TOLERANCE ***"))
  }
  a2_max <- if (length(diffs)) max(diffs) else NA_real_
  cat(sprintf("\n  A2 max global = %.5f (tolerance 0.02)\n", a2_max))
  if (is.finite(a2_max) && a2_max > 0.02)
    cat("  *** A2 VIOLEE : la re-estimation ne reproduit pas les dumps publies. ***\n")
} else cat("  dumps absents — A2 non evaluee\n")

##################################################################################
# 8. SORTIE
##################################################################################

## LES COURBES DE TOUTE LA GRILLE SONT ARCHIVEES (11/09/2026).
## CURVES[[key]] porte la sortie de get_effect_curve() -- la MEME fonction qui
## produit les douze valeurs de pe_spec_ranking.csv. Elle servait au seul
## diagnostic A2, puis disparaissait, et les douze specifications absentes de
## tab_spec_invariance n'avaient donc aucune source. Rien n'est reestime : on
## cesse de jeter.
if (length(CURVES_ALL)) {
  .cv <- rbindlist(CURVES_ALL, fill = TRUE)
  ## SMOKE ecrit sous un autre nom : une grille reduite a 1x1 ne doit pas
  ## atterrir sous le nom que l'emission lira.
  .cvf <- if (SMOKE) "spec_curves_daily_SMOKE.csv" else "spec_curves_daily.csv"
  fwrite(.cv, file.path(DIR_CONTRIB, .cvf))
  cat(sprintf("\nEcrit: contribA/%s (%d specifications, %d points)\n",
              .cvf, uniqueN(.cv$model_name), nrow(.cv)))
} else {
  ## Un mecanisme qui ne fait rien doit le dire.
  cat("\n[!] CURVES_ALL vide : spec_curves_daily.csv NON ecrit.\n")
}

if (length(CURVES_ANN)) {
  .cva <- rbindlist(CURVES_ANN, fill = TRUE)
  .cvaf <- if (SMOKE) "spec_curves_annual_SMOKE.csv" else "spec_curves_annual.csv"
  fwrite(.cva, file.path(DIR_CONTRIB, .cvaf))
  cat(sprintf("Ecrit: contribA/%s (%d specifications, %d points)\n",
              .cvaf, uniqueN(.cva$model_name), nrow(.cva)))
} else {
  cat("\n[!] CURVES_ANN vide : spec_curves_annual.csv NON ecrit.\n")
}

md5tab <- paste(sprintf("%s=%s", names(OBSERVED_MD5), OBSERVED_MD5), collapse = "; ")
R[, md5_inputs := md5tab]
## L'EMPREINTE est la colonne que ct2 et ct3 relisent. md5_inputs reste, mais
## comme trace lisible : c'est le jeton, pas la liste, qui gouverne le refus.
R[, inputs_fingerprint := INPUTS_FP]
R[, a2_max_diff := a2_max]
R[, a2_status := if (is.na(a2_max))
     "NON EVALUEE — aucune etape du paquet ne produit output/drf_prod_specs/"
   else "evaluee"]
R[, generated := "ct1_delta_north_grid.R"]

setcolorder(R, c("channel","arm","spec_label","form","df","winsor_lo","winsor_hi",
                 "origine","status","delta_north","se","ci_lo","ci_hi","excludes_zero",
                 "aic_weight","oos_mass_fut","oos_mass_base","support_lo","support_hi",
                 "psd_fix","n_eigen_trunc","var_raw","var_psd","n_obs","cluster_unit"))
fwrite(R, OUT_CSV)
cat("\nEcrit:", OUT_CSV, "\n")

cat("\n=============== TABLE ===============\n")
print(R[status == "OK", .(channel, arm, spec_label, delta_north = round(delta_north, 4),
                          se = round(se, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4),
                          excl0 = excludes_zero, aic_w = round(aic_weight, 3),
                          oos = round(100*oos_mass_fut, 2), psd = psd_fix)], nrows = 200)

cat("\n--- synthese : combien de specs excluent zero, par canal x bras ---\n")
print(R[status == "OK", .(n = .N, n_excl0 = sum(excludes_zero),
                          n_pos = sum(excludes_zero & delta_north > 0),
                          n_neg = sum(excludes_zero & delta_north < 0),
                          oos_med = round(100*median(oos_mass_fut), 2)),
        by = .(channel, arm)])

cat("\nTemps total:", round((proc.time()-t0)["elapsed"]/60, 1), "min\n=== Done ===\n")
