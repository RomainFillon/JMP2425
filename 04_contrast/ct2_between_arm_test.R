##################################################################################
# tache1_diff_bras.R — CONTRIBUTION A, TACHE 1 : LE TEST DE DIFFERENCE ENTRE BRAS
#
# CE QUE CE SCRIPT CALCULE
# ------------------------
# Pour chaque (canal, specification) : l'ECART entre l'increment nordique du bras
# annuel et celui du bras journalier,
#
#       D = Delta_north^annual - Delta_north^daily,
#
# son erreur-type, son intervalle a 95 %, et la CORRELATION rho entre les deux
# increments. Le Test 1 a etabli que le bras annuel est positif et significatif
# (Delta = +0,133) la ou le journalier est nul (+0,000) ; il n'a PAS etabli que les
# deux different, parce qu'il n'a jamais forme la covariance entre eux.
#
# POURQUOI CE TEST N'EST PAS UNE COQUETTERIE
# ------------------------------------------
# « significatif ici, pas la » est precisement la forme d'argument que Gelman et
# Stern (2006) montrent invalide : la difference entre « significatif » et « non
# significatif » n'est pas elle-meme significative. Un refere le soulevera. Si D
# exclut zero, la revendication change de nature — ce n'est plus une asymetrie de
# significativite, c'est l'estimation journaliere qui REJETTE l'estimation annuelle.
#
# CE QUI REND LE CALCUL NON TRIVIAL
# ---------------------------------
# Les deux bras sont deux regressions sur les MEMES observations, ajustees avec les
# MEMES effets fixes et les MEMES clusters. Leurs estimateurs sont donc fortement
# correles, et cette correlation n'apparait dans AUCUNE des deux vcov marginales.
# Var(D) = Var_a + Var_d - 2 Cov(a,d) : ignorer Cov revient a publier une variance
# fausse, presque surement trop grande (le precedent maison mesure rho a 0,67-0,89).
#
# LA ROUTE, ET POURQUOI ELLE EST CELLE-CI
# ---------------------------------------
# Pour un PPML, la fonction d'influence est psi_i = (X'WX)^-1 x_i (y_i - mu_i).
# fixest expose les deux morceaux : `scores` (les x_i (y_i - mu_i), sur les
# regresseurs demoyennes des FE absorbes) et `hessian` (le X'WX correspondant).
# On agrege les scores des DEUX bras sur les MEMES groupes de clusters et on forme
# le bloc croise exactement comme fixest forme les blocs propres.
#
# ⚠️ Cela suppose de savoir EXACTEMENT ce que fixest met dans un bloc propre —
# bread, viande, et corrections petit echantillon. Une erreur de ssc contaminerait
# le bloc croise du meme facteur SANS RIEN DECLENCHER. La formule a donc ete
# etablie par sonde, pas supposee : _work_contribA/tache1_probe_ssc.R la reproduit
# a 5,5e-16 sur donnees synthetiques. Elle est :
#
#   V_CGM = [(N-1)/(N-nparams)] . [G_min/(G_min-1)] . B (M_or + M_de - M_orXde) B
#
# avec B = solve(hessian), M_g = crossprod(rowsum(scores, g)), nparams COMPRENANT
# les effets fixes absorbes, et G_min = min(G_orig, G_dest).
#
# La meme sonde (section D) valide le BLOC CROISE lui-meme contre un systeme
# empile — la route du precedent maison `_work_crh_inference/step3_vcov.R` — et
# les deux coincident a 1,3e-4, du meme ordre que l'ecart des blocs PROPRES entre
# les deux routes (7,7e-5), c'est-a-dire au plancher de convergence de l'empilement
# et non a un ecart de structure.
#
#################################################################################
# PRE-ENGAGEMENT — ECRIT AVANT LA PREMIERE LIGNE DE CALCUL, ET AVANT TOUT RESULTAT
#################################################################################
#
# (P1) L'ECART SE CALCULE PAR PAIRE (forme, df, winsorisation) : autant de
#      differences par canal que GRID_SPECS porte de specifications (quinze
#      depuis le 18/09/2026). Aucune paire ne sera retiree apres avoir vu son
#      resultat, et le compte n'est ecrit nulle part en dur.
#
# (P2) L'APPARIEMENT EST NOMINAL, ET DOIT ETRE DIT COMME TEL. Les deux bras ne
#      vivent pas dans le meme espace de fonctions : la base daily a ses noeuds sur
#      la distribution journaliere, la base annual sur la distribution des moyennes
#      annuelles. « ns df=4 » ne designe donc pas la meme famille de fonctions dans
#      les deux bras. Apparier par etiquette est une convention de presentation, pas
#      une equivalence mathematique. C'est precisement pourquoi (P3).
#
# (P3) LA PRESENTATION EST L'ETENDUE SUR LA GRILLE, JAMAIS UNE VALEUR PONDEREE PAR
#      L'AIC. Le Test 1 a etabli qu'a N = 3,58 M les ecarts d'AIC entre formes se
#      comptent en dizaines de milliers de points : les poids sont degeneres PAR
#      CONSTRUCTION dans le canal productivite. Une moyenne ponderee par l'AIC y
#      serait la valeur d'une seule specification deguisee en synthese. L'AIC est
#      calcule et reporte en colonne, il n'entre dans aucune agregation.
#
# (P4) SPECIFICATION DE REFERENCE, DESIGNEE AVANT TOUT RESULTAT :
#         ns df = 4, winsorisation 0,5 %   (spec_label "spline4_wl0.5wu0.5")
#      Motif, entierement anterieur aux chiffres : c'est la specification du papier
#      publie, et 0,5 % est le reglage le MOINS ecretant des trois. Elle n'est pas
#      choisie parce qu'un resultat l'aurait recommandee. Si son resultat est
#      defavorable, il est rapporte comme le resultat de reference.
#
# (P5) NON-DEFINIE-POSITIVITE — TRAITEMENT FIXE AVANT DE SAVOIR SI LE CAS SE PRODUIT.
#      La vcov jointe est de forme CGM V1 + V2 - V12 ; une valeur propre negative
#      est un cas courant, et Var(D) = u' V u peut sortir NEGATIVE.
#        (1) on calcule Var(D) sur la matrice jointe 2k x 2k BRUTE ; si >= 0, on la
#            publie et psd_fix vaut "none" ;
#        (2) si < 0, on tronque a zero les valeurs propres negatives de la matrice
#            JOINTE — et non des blocs separement, car tronquer bloc par bloc ne
#            garantit pas la semi-definie-positivite de la jointe — on recalcule, et
#            psd_fix vaut "eigen_trunc" ;
#        (3) var_raw ET var_psd sont toutes deux reportees, avec le nombre de valeurs
#            propres tronquees et la plus negative d'entre elles ;
#        (4) rho est rapporte sur la matrice brute ET sur la matrice publiee.
#      Aucune ligne n'est retiree pour cette raison.
#
# (P6) SSC DU BLOC CROISE. Il n'existe pas de convention canonique pour la correction
#      petit echantillon d'un bloc CROISE entre deux modeles. Pre-engagement : on
#      applique la MOYENNE GEOMETRIQUE des facteurs scalaires propres aux deux bras,
#      sqrt(c_daily . c_annual). Trois raisons, toutes anterieures aux chiffres :
#        (i)   elle est symetrique dans les deux bras, ce qu'une covariance doit etre ;
#        (ii)  elle se reduit EXACTEMENT au facteur commun quand les deux bras ont
#              meme N et meme nparams — le cas attendu, et qui sera RAPPORTE ;
#        (iii) la jointe corrigee s'ecrit D V_brut D avec D = diag(sqrt(c_d) I,
#              sqrt(c_a) I) : c'est un rescaling qui PRESERVE la semi-definie-
#              positivite, donc la correction ne peut pas a elle seule fabriquer une
#              variance negative.
#      Var(D) est en outre reportee SANS aucune correction (var_diff_nossc) pour que
#      la sensibilite a ce choix soit lisible et non a reconstituer.
#
# (P7) ECHANTILLONS. Les blocs PROPRES sont calcules sur l'echantillon propre de
#      chaque bras (c'est la seule facon de reproduire vcov() de fixest, cf. B1).
#      Le bloc CROISE est calcule sur l'INTERSECTION des deux echantillons — ce qui
#      est exact et non une approximation : une observation absente d'un bras y a un
#      score nul, donc elle contribue zero au produit croise. Le recouvrement est
#      rapporte en colonne (n_obs_daily, n_obs_annual, n_obs_common).
#
#################################################################################
# ASSERTIONS PERMANENTES (pas des verifications de passage)
#################################################################################
#
# B1 — BLOQUANTE, et c'est l'analogue gratuit du test a lambda = 0.
#      Ma reconstruction manuelle du bloc PROPRE de chaque bras doit reproduire
#      vcov(mod) de fixest sur le meme objet ajuste, A UN SCALAIRE PRES et a 1e-8 en
#      norme relative. Le scalaire est la correction petit echantillon, une
#      convention ; ce qui doit etre exact est la STRUCTURE — viande, bread, groupes
#      de clusters — parce que c'est elle, et non le scalaire, que le bloc croise
#      reutilise. Si les blocs propres sont structurellement exacts, le bloc croise,
#      bati par le MEME code sur les MEMES scores, est credible. Sinon le script
#      s'arrete : un bloc croise faux ne se signalerait nulle part ailleurs.
#      Le scalaire recupere est confronte a sa forme analytique et l'ecart est
#      rapporte plutot que corrige (il est reel, et documente : cf. ssc_recover).
#
# B2 — BLOQUANTE. Les Delta_north que ce script recalcule doivent reproduire
#      _work_contribA/test1_delta_north.csv a 1e-9. Ce script re-derive la geographie
#      nordique et les poids w ; B2 est la preuve que cette re-derivation n'a pas
#      diverge de celle du Test 1. Sans elle, une difference pourrait venir d'un
#      ecart de geographie plutot que d'un fait sur les donnees.
#
# B3 — RAPPORTEE. Les deux bras tournent-ils sur le meme echantillon ? Attendu oui.
#      Le calcul est exact dans les deux cas (P7), mais l'ecart doit etre VU.
#
# Les fonctions d'estimation ne sont pas recopiees : elles sont EXTRAITES verbatim
# de drf_prod_final.r / drf_amenity_final.r, et les constructeurs de base
# daily_basis / annual_basis / delta_var sont EXTRAITS verbatim de
# test1_delta_north.R. Une seule source de verite, et son md5 est bloquant.
#
# SORTIE : _work_contribA/tache1_diff_bras.csv    (aucune ecriture dans output/)
#          _work_contribA/tache1_scores.rds       (viandes k x k + scores par cluster)
#
# Romain Fillon / Contribution A, Tache 1 — 2026-07-31
##################################################################################

suppressMessages({
  library(data.table); library(splines); library(fixest); library(tools)
})
options(width = 220)
t0 <- proc.time(); cat("=== ct2_between_arm_test.R ===\n")
cat("fixest", as.character(packageVersion("fixest")), "|", R.version.string, "\n")

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
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "04_contrast", "lib_contrast.R"))

DIR_CONTRIB <- file.path(DIR_OUT, "contribA")
if (!dir.exists(DIR_CONTRIB)) dir.create(DIR_CONTRIB, recursive = TRUE)
DIR_TMP_CONTRIB <- file.path(DIR_TMP, "contribA")
if (!dir.exists(DIR_TMP_CONTRIB)) dir.create(DIR_TMP_CONTRIB, recursive = TRUE)

BASE_YEAR <- 2015L; FUT_YEAR <- 2095L; NORTH_LAT <- 45
POP_YEAR  <- 2005L
## (P4), designee avant tout resultat -- et LUE, non retapee, depuis le
## 18/09/2026. Deux raisons : ct3 (Python) en portait une seconde copie, et le
## nom `REF_SPEC` ECRASAIT ici celui de 00_parameters.R, qui vaut "polynomial"
## et designe la reference des diagnostics d'identification. Deux objets sous un
## nom, dans le meme espace de noms, apres un source() : c'est le genre de
## collision qui ne se voit qu'au moment ou elle a deja menti.
## => la designation vit dans GRID_REF_SPEC (00_parameters.R). Aucune
## affectation ici : ce script la LIT.

##################################################################################
# 0. LES ENTREES, ET LE REFUS DE PROVENANCE
#
# Sept entrees : les six de ct1, plus la sortie de ct1 elle-meme.
#
# ⛔ CE SCRIPT REFUSE DE TOURNER SI SON EMPREINTE DIFFERE DE CELLE DE ct1. C'est
# le remplacant de l'epinglage md5, et il vise le danger que l'epinglage visait :
# l'APPARIEMENT CROISE. Ce test compose les increments des DEUX bras et leur
# covariance croisee ; si les deux ne viennent pas du meme etat d'entree, la
# difference reste arithmetiquement bien formee et ne veut plus rien dire.
#
# ⚠ CE QUE L'EMPREINTE NE VOIT PAS : deux artefacts du meme run estimes sous des
# SPECIFICATIONS differentes. Cela se traite en ecrivant les controles a chaque
# site d'appel — voir §1, et le correctif de 01a du 2026-09-04.
##################################################################################

PATHS <- c(
  "BACI_raw_pairs_annual.rds" = file.path(DIR_TMP, "BACI_raw_pairs_annual.rds"),
  "climate_long_annual.rds"   = file.path(DIR_TMP, "climate_long_annual.rds"),
  "step2_migration.rds"       = file.path(DIR_TMP, "step2_migration.rds"),
  "climate_long.rds"          = file.path(DIR_TMP, "climate_long.rds"),
  "CMIP6_tas.csv"             = file.path(DIR_TMP, "CMIP6_tas.csv"),
  "amen_prod.csv"             = file.path(DIR_CRH, "amen_prod.csv")
)
T1_CSV <- file.path(DIR_CONTRIB, "test1_delta_north.csv")

cat("\n--- entrees et provenance ---\n")
manque <- names(PATHS)[!file.exists(PATHS)]
if (length(manque))
  stop("[STOP] entree(s) absente(s) : ", paste(manque, collapse = ", "),
       "\n  Lancer d'abord l'etape qui les produit — voir run_all.sh.")
if (!file.exists(T1_CSV))
  stop("[STOP] ", basename(T1_CSV), " absent.\n",
       "  Il est produit par ct1_delta_north_grid.R, qui doit tourner avant ce\n",
       "  script : la covariance croisee se calcule ENTRE les deux bras que ct1\n",
       "  a estimes, pas a partir de rien.")

OBSERVED_MD5 <- setNames(as.character(tools::md5sum(PATHS)), names(PATHS))
for (n in names(PATHS)) cat(sprintf("  ok  %-28s %s\n", n, OBSERVED_MD5[n]))
INPUTS_FP <- empreinte(OBSERVED_MD5)
cat(sprintf("  empreinte de provenance : %s\n", INPUTS_FP))

## LE REFUS. ct1 inscrit son empreinte en colonne ; on la relit et on compare.
.t1_head <- data.table::fread(T1_CSV, nrows = 1L)
if (!"inputs_fingerprint" %in% names(.t1_head))
  stop("[STOP] ", basename(T1_CSV), " ne porte pas de colonne inputs_fingerprint.\n",
       "  C'est un artefact anterieur a la migration du 2026-09-04. Le relancer\n",
       "  par ct1_delta_north_grid.R plutot que de desarmer ce controle.")
.t1_fp <- as.character(.t1_head$inputs_fingerprint[1])
if (!identical(.t1_fp, INPUTS_FP))
  stop("[STOP] APPARIEMENT CROISE.\n",
       "  ct1 a tourne sur l'empreinte ", .t1_fp, "\n",
       "  ct2 lit aujourd'hui l'empreinte ", INPUTS_FP, "\n",
       "  Les increments par bras et la covariance croisee ne porteraient pas sur\n",
       "  le meme etat d'entree. La difference resterait bien formee et n'aurait\n",
       "  aucun sens. Relancer ct1 sur l'etat courant.")
cat("  ✓ empreinte identique a celle de ct1\n")

##################################################################################
# 1. LES FONCTIONS D'ESTIMATION — lib_drf.R, ET LES ARGUMENTS QUI LES SPECIALISENT
#
# Voir l'en-tete de ct1 pour le dossier complet. Le resume tient en une phrase :
# les quatre variantes trade/migration x daily/annual ont ete unifiees en
# prepare_diffs() / run_stage3(), et les cinq arguments qui les specialisent
# s'ecrivent ici plutot que de venir d'un defaut.
#
# ⛔ controls SURTOUT. run_stage3() a controls = character(0) ; l'ancienne API
# avait controls = gravity_controls. Mesure du 2026-09-04 : Delta_north +22,9 %
# sans les controles cote prod, ecart d'AIC 137 502.
#
# daily_basis / annual_basis viennent de lib_contrast.R et non plus d'un parsing
# de ct1. delta_var de ct1 n'est PAS employee ici : elle traite un bloc k x k, la
# ou (P5) traite la matrice JOINTE 2k x 2k. Importer une fonction inutilisee
# laisserait croire qu'elle agit.
##################################################################################

ARGS <- list(
  prod = list(
    controls = TRADE_STAGE3_CONTROLS,
    fe       = "fe_prod_year",
    sign_dmo = FALSE,
    ## 19/09/2026 : etait `scale = 1`, commente « aucune mise a l'echelle ».
    ## Le PPML identifie theta * f_z (eq. 5) et le GE reapplique theta dans
    ## trade_weights : sans cette division le choc etait amplifie deux fois.
    ## Le canal amenite faisait deja le geste, avec Omega, deux lignes plus bas.
    scale    = SCALE_PROD),
  amen = list(
    controls = AMENITY_STAGE3_CONTROLS,          # 01b l.269
    fe       = c("fe_orig_dest", "fe_year"),     # 01b l.266 et l.270
    sign_dmo = TRUE,                             # migration : dest - orig
    scale    = Omega)                            # 01b l.267
)
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
# 2. LA GRILLE — identique au Test 1 PARCE QU'ELLE VIENT DU MEME ENDROIT
##################################################################################
# 🚨 « identique au Test 1 » ETAIT UNE PROMESSE, PAS UNE PROPRIETE, JUSQU'AU
# 2026-09-18 : ce bloc portait sa propre copie de la liste des huit formes, et
# ct1 la sienne. Deux copies d'une meme decision sont deux endroits ou la
# changer. ct4 joint les deux sorties par etiquette et se serait arrete sur un
# `setequal` en cas de divergence -- au mieux ; au pire deux grilles de meme
# taille auraient joint proprement des objets differents.
# ⇒ Les deux scripts LISENT desormais GRID_SPECS (00_parameters.R). La grille
#   compte QUINZE specifications depuis la decision d'auteur du 18/09/2026 :
#   polynomes de degre 2 a 4, splines naturelles a 3 et 4 df, trois
#   winsorisations. Aucune paire n'est retiree apres avoir vu son resultat.

SPEC_CONFIGS   <- GRID_FORMS
WINSOR_CONFIGS <- GRID_WINSOR

## Meme controle qu'en ct1 : les etiquettes que CE script va produire doivent
## etre exactement celles de GRID_SPECS. C'est par elles que ct4 joint.
.labels_locaux <- unlist(lapply(SPEC_CONFIGS, function(cfg)
  vapply(WINSOR_CONFIGS, function(wc)
    sprintf("%s%d_wl%gwu%g", cfg$spec, cfg$degree, wc$lower, wc$upper),
    character(1))))
if (!setequal(.labels_locaux, GRID_LABELS) ||
    length(.labels_locaux) != GRID_NSPEC)
  stop("[FATAL] la grille parcourue par ce script n'est pas GRID_SPECS.\n",
       "  parcourue : ", paste(sort(.labels_locaux), collapse = ", "), "\n",
       "  GRID_SPECS: ", paste(sort(GRID_LABELS),   collapse = ", "))

SMOKE <- nzchar(Sys.getenv("T1D_SMOKE"))
if (SMOKE) {
  ## PAR NOM, JAMAIS PAR INDICE : cette ligne a porte `SPEC_CONFIGS[4]` sous le
  ## commentaire « spline df=4 ». Le rang 4 de GRID_FORMS est la spline df=3.
  ## La reference du pre-engagement (P4) est GRID_REF_SPEC : on la retrouve
  ## par sa forme et son degre, pas par sa place dans une liste.
  .i <- which(vapply(SPEC_CONFIGS,
                     function(z) z$spec == "spline" && z$degree == 4L, TRUE))
  if (length(.i) != 1L)
    stop("[FATAL] mode smoke : la spline df=4 n'est pas dans GRID_FORMS, ou y ",
         "figure ", length(.i), " fois.")
  SPEC_CONFIGS   <- SPEC_CONFIGS[.i]     # spline df=4 : la specification de reference
  WINSOR_CONFIGS <- WINSOR_CONFIGS[1]    # 0,5 %
  cat("\n*** MODE SMOKE : grille reduite a 1 x 1. Sortie NON publiable. ***\n")
}
## Le CSV est un LIVRABLE (emit_values_empirics.py et ct4 le lisent) : il va
## dans output/contribA/. Le .rds est un objet de travail que personne ne lit :
## il va dans intermediate_outputs/contribA/. Meme partage que la chaine CRH
## entre output/crh_inference/ et intermediate_outputs/crh/.
OUT_CSV <- file.path(DIR_CONTRIB,
                     if (SMOKE) "tache1_SMOKE.csv" else "tache1_diff_bras.csv")
OUT_RDS <- file.path(DIR_TMP_CONTRIB,
                     if (SMOKE) "tache1_SMOKE_scores.rds" else "tache1_scores.rds")
cat(sprintf("\nGrille : %d formes x %d winsorisations = %d paires par canal\n",
            length(SPEC_CONFIGS), length(WINSOR_CONFIGS),
            length(SPEC_CONFIGS) * length(WINSOR_CONFIGS)))
cat(sprintf("Specification de reference (P4) : %s\n", GRID_REF_SPEC))

##################################################################################
# 3. GEOGRAPHIE NORDIQUE — re-derivee ici, VALIDEE par B2 contre le Test 1
##################################################################################

cat("\n--- geographie nordique ---\n")
cm <- fread(PATHS["CMIP6_tas.csv"])[year %in% c(BASE_YEAR, FUT_YEAR) & y >= NORTH_LAT]
cm[, c := count / sum(count), by = .(x, y, year)]
pop <- fread(PATHS["amen_prod.csv"])[year == POP_YEAR, .(lat, lon, pop = as.numeric(pop_x))]
cm  <- merge(cm, pop, by.x = c("x", "y"), by.y = c("lon", "lat"), all = FALSE)
ok2 <- cm[, .(ny = uniqueN(year)), by = .(x, y)][ny == 2L, .(x, y)]
cm  <- merge(cm, ok2, by = c("x", "y"))
ncell <- nrow(ok2)
cat(sprintf("  cellules y >= %g avec base ET futur : %d\n", NORTH_LAT, ncell))

wcell <- unique(cm[, .(x, y, pop)]); wcell[, omega := pop / sum(pop)]
stopifnot(abs(sum(wcell$omega) - 1) < 1e-12)
cm <- merge(cm, wcell[, .(x, y, omega)], by = c("x", "y"))

cm[, sgn := fifelse(year == FUT_YEAR, 1, -1)]
dmass <- cm[, .(m = sum(omega * sgn * c)), by = value][order(value)]
share_by_T <- cm[, .(s = sum(omega * c)), by = .(year, value)]
stopifnot(all(abs(share_by_T[, .(t = sum(s)), by = year]$t - 1) < 1e-9))

mu  <- cm[, .(mu = sum(value * c)), by = .(x, y, year)]
muw <- dcast(mu, x + y ~ year, value.var = "mu")
setnames(muw, c(as.character(BASE_YEAR), as.character(FUT_YEAR)), c("mu_base", "mu_fut"))
muw <- merge(muw, wcell[, .(x, y, omega)], by = c("x", "y"))
muw <- muw[is.finite(mu_base) & is.finite(mu_fut)]
rm(cm); gc(verbose = FALSE)

##################################################################################
# 4. MACHINERIE CGM — la formule etablie par tache1_probe_ssc.R
##################################################################################

# indices des lignes CONSERVEES par fixest, dans l'ordre des donnees fournies
rowsel <- function(mod, n_full) {
  r <- seq_len(n_full)
  if (!is.null(mod$obs_selection) && length(mod$obs_selection))
    for (s in mod$obs_selection) r <- r[s]
  r
}

# viande clusteree : somme des produits exterieurs des scores agreges par groupe.
# Deux matrices de scores => bloc CROISE ; une seule => bloc PROPRE.
meat_cross <- function(S1, S2, g_or, g_de, g_pair) {
  ag <- function(S, g) rowsum(S, group = g, reorder = TRUE)
  crossprod(ag(S1, g_or),   ag(S2, g_or)) +
  crossprod(ag(S1, g_de),   ag(S2, g_de)) -
  crossprod(ag(S1, g_pair), ag(S2, g_pair))
}

# --- facteur scalaire de correction petit echantillon ---------------------------
# Forme ANALYTIQUE, etablie par tache1_probe_ssc.R sur un modele a UN effet fixe :
# ⚠️ ssc_analytic N'ENTRE DANS AUCUN CALCUL. Sa valeur est affichee et ecrite en
# colonne (ssc_*_analytic) a cote du scalaire reellement utilise (ssc_*), et les
# deux diffèrent de 20 % sur le canal amenite. Deux colonnes voisines dont l'une
# porte « analytic » invitent a croire qu'elles se controlent : elles ne le font
# pas, et B1 ne teste PAS le scalaire (verifie par mutation, cf.
# tache1_mutation_B1.R : l'injection d'une erreur scalaire de 20 % sort VERTE).
# Ce n'est pas une lacune : ce qui est publie est c_RECUPERE . V_brut, ou c_recupere
# est extrait de la vcov de fixest elle-meme ; joint a B1 (norme relative < 1e-8),
# cela garantit que le bloc propre publie EST celui de fixest. Il n'y a donc rien a
# controler sur le scalaire — mais il faut que ce soit ecrit, pas deduit.
ssc_analytic <- function(N, nparams, G_or, G_de) {
  Gmin <- min(G_or, G_de)
  ((N - 1) / (N - nparams)) * (Gmin / (Gmin - 1))
}
# ⚠️ MAIS elle ne vaut PAS telle quelle partout, et le smoke du 2026-07-31 l'a
# montre : le canal amenite l'a violee de 20 %. Cause etablie, non supposee —
# fixest applique ssc(fixef.K = "nested") : un effet fixe IMBRIQUE dans une variable
# de clustering ne compte pas dans K. `fe_orig_dest` est exactement orig x dest,
# donc imbrique dans les deux clusters ; ses ~14 000 parametres sont exclus de K par
# fixest et comptes par mod$nparams. Avec N = 86 172, l'ecart sur (N-1)/(N-K) vaut
# 86171/72172 = 1,194 — les 20 % observes, au chiffre pres. Le canal productivite
# passait parce que fe_prod_year n'est imbrique dans aucun cluster.
#
# On ne replique donc PAS la comptabilite de K de fixest — ce serait reconstruire
# une convention mouvante et se tromper ailleurs. On RECUPERE le facteur, et on
# asserte ce qui doit vraiment etre vrai pour que le bloc croise soit valide :
# que ma viande et ma bread reproduisent celles de fixest A UN SCALAIRE PRES.
# C'est une assertion PLUS forte que « le scalaire est le bon », car elle porte sur
# la STRUCTURE de la matrice, qui est ce que le bloc croise reutilise.
#
# ⚠️ CHOIX DE LA METRIQUE — corrige le 31/07 apres que B1 a bloque a df=6.
# La premiere version mesurait l'ecart ENTREE PAR ENTREE, max|r_ij - c|/c sur le
# ratio r_ij = V_fixest_ij / V_brut_ij. Cette metrique est dominee par les plus
# PETITES entrees retenues, dont l'erreur relative est la plus grande : elle mesure
# le conditionnement du TEST autant que celui de la matrice. A df=6 elle a rendu
# 2,23e-9 alors que le scalaire recupere egalait l'analytique au chiffre pres
# (1,007443) — ce qu'un desaccord de structure ne produit pas.
#
# La metrique retenue est la NORME RELATIVE, ||V_fixest - c V_brut||_max / ||V_fixest||_max,
# et ce n'est pas un assouplissement : c'est la quantite qui gouverne reellement
# l'erreur sur Var(D). Var(D) est une forme quadratique u' V u ; une erreur de eps
# en norme relative sur V y produit au plus ~eps d'erreur relative. Le seuil est
# fixe a 1e-8, ce qui borne l'erreur sur l'erreur-type a ~5e-9 en relatif, soit six
# ordres de grandeur sous la precision rapportee (4 decimales) — tout en restant
# des ordres de grandeur en dessous de ce que produirait une erreur de structure,
# qui serait O(1) ou au mieux O(1e-2).
# L'ecart entree-par-entree reste CALCULE ET REPORTE en colonne : rien n'est masque.
#
# 🔑 B1 EST DEMONTRABLEMENT CAPABLE D'ECHOUER : elle a bloque le run du 31/07 a la
# paire 19/48 (prod/spline6, ecart 2,23e-9 contre seuil 1e-9), et elle avait deja
# bloque le smoke sur le canal amenite (ecart 2,0e-1, FE imbrique dans un cluster).
# « Un test qui ne peut pas echouer ne vaut rien » — celui-ci a echoue deux fois,
# pour deux causes distinctes, dont une VRAIE (l'amenite).
ssc_recover <- function(V_raw, V_fixest) {
  # c par moindres carres sur la matrice entiere : minimise ||V_fixest - c V_brut||_F
  c_hat <- sum(V_fixest * V_raw) / sum(V_raw * V_raw)
  rel_norm <- max(abs(V_fixest - c_hat * V_raw)) / max(abs(V_fixest))
  sel <- abs(V_raw) > max(abs(V_raw)) * 1e-8
  r   <- V_fixest[sel] / V_raw[sel]
  list(c = c_hat, rel_norm = rel_norm,
       spread = max(abs(r - median(r))) / abs(median(r)), n_used = sum(sel))
}

##################################################################################
# 5. BOUCLE PRINCIPALE — les DEUX bras dans la MEME iteration
##################################################################################
# Les objets modeles ne sont persistes nulle part (constat de l'inventaire de
# Phase 1 : le Test 1 a estime 96 modeles et les a jetes). Le bloc croise exige les
# scores des deux bras simultanement, sur les memes groupes : il ne peut donc PAS
# etre reconstruit a partir des sorties du Test 1. On ajuste les deux bras dans la
# meme iteration et on forme le bloc croise immediatement.

## Plus de noms de fonctions dans le dictionnaire : les quatre variantes ont
## disparu au profit d'une generique. Il reste les donnees et les motifs de
## colonnes ; les arguments d'estimation sont dans ARGS (§1).
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

RES <- list(); STORE <- list(); k_it <- 0L

for (chan in names(CHANNELS)) {
  CH <- CHANNELS[[chan]]; AG <- ARGS[[chan]]
  cat(sprintf("\n################ CANAL %s ################\n", toupper(chan)))
  dat <- readRDS(CH$data_path); setDT(dat)
  if (!is.na(CH$year_min)) dat <- dat[year >= CH$year_min]
  clim <- readRDS(CH$clim_path); setDT(clim)
  cat(sprintf("  donnees : %d lignes | climat : %d lignes\n", nrow(dat), nrow(clim)))

  for (cfg in SPEC_CONFIGS) {
    spec <- cfg$spec; deg <- cfg$degree
    for (wc in WINSOR_CONFIGS) {
      wl <- wc$lower; wu <- wc$upper
      key <- sprintf("%s%d_wl%gwu%g", spec, deg, wl, wu)
      k_it <- k_it + 1L
      cat(sprintf("\n[%02d] %s / %s\n", k_it, chan, key)); flush.console()
      tk <- proc.time()

      # ---------- les deux bras, dans cette iteration ----------
      ARM <- list()
      failed <- NULL
      for (arm in c("daily", "annual")) {
        cat(sprintf("     %-6s ajustement ... ", arm)); flush.console()
        # wald() de run_stage3 imprime son resultat sur stdout : on l'absorbe pour
        # garder le log lisible. Les erreurs passent par les conditions, pas par
        # stdout, donc tryCatch les voit toujours.
        got <- tryCatch({
          # wald() de run_stage3 imprime sur stdout : capture.output l'absorbe tout
          # en evaluant dans CET environnement (donc les affectations persistent).
          # Les erreurs passent par les conditions et non par stdout : tryCatch les
          # voit toujours.
          invisible(capture.output({
            if (arm == "daily") {
              dp <- prepare_diffs(dat, clim, spec, deg, wl, wu, FALSE,
                                  fe_cols              = AG$fe,
                                  sign_dest_minus_orig = AG$sign_dmo,
                                  scale                = AG$scale)
              m  <- run_stage3(dp, spec, deg, FALSE,
                               controls = AG$controls, fe_str = AG$fe)
            } else {
              dp <- prepare_diffs(dat, clim, spec, deg, wl, wu, TRUE,
                                  fe_cols              = AG$fe,
                                  sign_dest_minus_orig = AG$sign_dmo,
                                  scale                = AG$scale)
              m  <- run_stage3(dp, spec, deg, TRUE,
                               controls = AG$controls, fe_str = AG$fe)
            }
          }))
          mod <- m$model
          rws <- rowsel(mod, nrow(dp))
          stopifnot(length(rws) == nrow(mod$scores))
          # cles de cluster alignees sur les lignes de scores, PUIS on libere dp
          keys <- list(or = as.character(dp$orig[rws]), de = as.character(dp$dest[rws]))
          rm(dp); gc(verbose = FALSE)
          list(mod = mod, rows = rws, keys = keys, wald = m$wald_joint)
        }, error = function(e) { cat("ECHEC:", conditionMessage(e), "\n"); NULL })
        if (is.null(got)) { failed <- arm; break }
        cat(sprintf("N = %d, %d coefs\n", nobs(got$mod), length(coef(got$mod))))
        ARM[[arm]] <- got
      }

      if (!is.null(failed)) {
        RES[[length(RES)+1L]] <- data.table(
          channel = chan, spec_label = key, form = spec, df = deg,
          winsor_lo = wl, winsor_hi = wu, origine = cfg$origine,
          status = paste0("FIT_FAILED_", failed))
        rm(ARM); gc(verbose = FALSE); next
      }

      # ---------- B3 : les deux bras tournent-ils sur le meme echantillon ? ----------
      rd <- ARM$daily$rows; ra <- ARM$annual$rows
      common <- intersect(rd, ra)
      same_sample <- identical(rd, ra)
      cat(sprintf("     echantillons : daily = %d | annual = %d | commun = %d | identiques = %s\n",
                  length(rd), length(ra), length(common), same_sample))

      # ---------- coefficients climatiques, poids w, Delta ----------
      PACK <- list()
      for (arm in c("daily", "annual")) {
        mod <- ARM[[arm]]$mod
        pat <- sprintf(if (arm == "daily") CH$pat_daily else CH$pat_annual, spec)
        cc  <- grep(pat, names(coef(mod)), value = TRUE)
        if (arm == "daily") {
          B   <- daily_basis(clim, spec, deg, wl, wu)
          Tcl <- pmin(pmax(dmass$value, B$lo), B$hi)
          w   <- as.numeric(t(predict(B$basis, Tcl)) %*% dmass$m)
          oos <- sum(share_by_T[year == FUT_YEAR & (value < B$lo | value > B$hi), s])
        } else {
          B  <- annual_basis(clim, spec, deg, wl, wu, mean_climate)
          Xf <- predict(B$basis, pmin(pmax(muw$mu_fut,  B$lo), B$hi))
          Xb <- predict(B$basis, pmin(pmax(muw$mu_base, B$lo), B$hi))
          w  <- as.numeric(t(Xf - Xb) %*% muw$omega)
          oos <- sum(muw$omega[muw$mu_fut < B$lo | muw$mu_fut > B$hi])
        }
        stopifnot(length(w) == length(cc))
        # w est etendu par des zeros sur les controles : la bread est pleine (k x k)
        wfull <- setNames(numeric(length(coef(mod))), names(coef(mod)))
        wfull[cc] <- w
        PACK[[arm]] <- list(cc = cc, beta = coef(mod)[cc], w = w, wfull = wfull,
                            delta = as.numeric(crossprod(w, coef(mod)[cc])),
                            oos = oos, lo = B$lo, hi = B$hi,
                            aic = AIC(mod), nobs = nobs(mod), nparams = mod$nparams)
      }

      # ---------- blocs propres, sur l'echantillon PROPRE de chaque bras (P7) ----------
      OWN <- list(); B1max <- 0; B1spr <- 0
      for (arm in c("daily", "annual")) {
        mod <- ARM[[arm]]$mod
        S   <- as.matrix(mod$scores)
        kk  <- ARM[[arm]]$keys
        gp  <- paste(kk$or, kk$de, sep = "\r")
        Br  <- solve(mod$hessian)
        M   <- meat_cross(S, S, kk$or, kk$de, gp)
        Vraw <- Br %*% M %*% Br
        Vfx  <- vcov(mod)
        # --- B1 : BLOQUANTE. Ma structure doit egaler celle de fixest a un SCALAIRE
        #     pres, et ce scalaire doit etre UNIFORME sur toutes les entrees. ---
        rc  <- ssc_recover(Vraw, Vfx)
        B1max <- max(B1max, rc$rel_norm)
        B1spr <- max(B1spr, rc$spread)
        cA  <- ssc_analytic(nobs(mod), mod$nparams, uniqueN(kk$or), uniqueN(kk$de))
        cat(sprintf("     B1 %-6s : norme rel = %.3e | entree-par-entree = %.3e | ssc = %.6f | analytique = %.6f%s\n",
                    arm, rc$rel_norm, rc$spread, rc$c, cA,
                    if (abs(rc$c - cA)/cA > 1e-6) "  (FE imbrique dans un cluster)" else ""))
        if (!is.finite(rc$rel_norm) || rc$rel_norm >= 1e-8)
          stop(sprintf(paste0("B1 VIOLEE (%s/%s/%s) : ||V_fixest - c V_brut|| / ||V_fixest|| ",
                              "= %.3e >= 1e-8. Ma viande et ma bread ne reproduisent pas ",
                              "celles de fixest a un scalaire pres : la STRUCTURE differe, ",
                              "et le bloc croise n'a aucune garantie."), chan, arm, key, rc$rel_norm))
        OWN[[arm]] <- list(V = rc$c * Vraw, bread = Br, ssc = rc$c, ssc_ana = cA,
                           meat = M, G_or = uniqueN(kk$or), G_de = uniqueN(kk$de))
      }

      # ---------- bloc croise, sur l'INTERSECTION (P7) ----------
      pos_d <- match(common, rd); pos_a <- match(common, ra)
      # quand les echantillons coincident (cas attendu), on evite deux copies de
      # ~400 Mo en indexant a l'identique
      Sd <- if (same_sample) as.matrix(ARM$daily$mod$scores)
            else as.matrix(ARM$daily$mod$scores)[pos_d, , drop = FALSE]
      Sa <- if (same_sample) as.matrix(ARM$annual$mod$scores)
            else as.matrix(ARM$annual$mod$scores)[pos_a, , drop = FALSE]
      g_or <- ARM$daily$keys$or[pos_d]; g_de <- ARM$daily$keys$de[pos_d]
      # controle : les cles de cluster de l'intersection doivent coincider entre bras
      stopifnot(identical(g_or, ARM$annual$keys$or[pos_a]),
                identical(g_de, ARM$annual$keys$de[pos_a]))
      g_pair <- paste(g_or, g_de, sep = "\r")
      M_cross <- meat_cross(Sd, Sa, g_or, g_de, g_pair)
      ssc_x   <- sqrt(OWN$daily$ssc * OWN$annual$ssc)                      # (P6)
      V_cross <- ssc_x * (OWN$daily$bread %*% M_cross %*% OWN$annual$bread)
      rm(Sd, Sa); gc(verbose = FALSE)

      # ---------- matrice jointe et variance de la difference ----------
      kd <- length(coef(ARM$daily$mod)); ka <- length(coef(ARM$annual$mod))
      Vj <- matrix(0, kd + ka, kd + ka)
      Vj[1:kd, 1:kd]                     <- OWN$daily$V
      Vj[(kd+1):(kd+ka), (kd+1):(kd+ka)] <- OWN$annual$V
      Vj[1:kd, (kd+1):(kd+ka)]           <- V_cross
      Vj[(kd+1):(kd+ka), 1:kd]           <- t(V_cross)
      # D = annual - daily  =>  u = (-w_daily, +w_annual)
      u <- c(-PACK$daily$wfull, PACK$annual$wfull)

      var_raw <- as.numeric(t(u) %*% Vj %*% u)
      # (P5) : traitement fixe avant de savoir si le cas se produit
      if (is.finite(var_raw) && var_raw >= 0) {
        var_pub <- var_raw; psd_fix <- "none"; ntr <- 0L; minev <- NA_real_; Vpub <- Vj
      } else {
        ei <- eigen((Vj + t(Vj))/2, symmetric = TRUE)
        ntr <- sum(ei$values < 0); minev <- min(ei$values)
        Vpub <- ei$vectors %*% diag(pmax(ei$values, 0), nrow(Vj)) %*% t(ei$vectors)
        var_pub <- as.numeric(t(u) %*% Vpub %*% u); psd_fix <- "eigen_trunc"
      }
      # variance SANS aucune correction petit echantillon (P6, sensibilite)
      Vj0 <- Vj
      Vj0[1:kd, 1:kd]                     <- OWN$daily$V  / OWN$daily$ssc
      Vj0[(kd+1):(kd+ka), (kd+1):(kd+ka)] <- OWN$annual$V / OWN$annual$ssc
      Vj0[1:kd, (kd+1):(kd+ka)]           <- V_cross / ssc_x
      Vj0[(kd+1):(kd+ka), 1:kd]           <- t(V_cross / ssc_x)
      var_nossc <- as.numeric(t(u) %*% Vj0 %*% u)

      # variances marginales et covariance de (Delta_daily, Delta_annual)
      ud <- c(PACK$daily$wfull, numeric(ka)); ua <- c(numeric(kd), PACK$annual$wfull)
      v_d   <- as.numeric(t(ud) %*% Vpub %*% ud)
      v_a   <- as.numeric(t(ua) %*% Vpub %*% ua)
      cov_da<- as.numeric(t(ud) %*% Vpub %*% ua)
      rho   <- cov_da / sqrt(v_d * v_a)
      # rho sur la matrice BRUTE aussi (P5-4)
      v_d0  <- as.numeric(t(ud) %*% Vj %*% ud); v_a0 <- as.numeric(t(ua) %*% Vj %*% ua)
      rho0  <- as.numeric(t(ud) %*% Vj %*% ua) / sqrt(v_d0 * v_a0)

      D  <- PACK$annual$delta - PACK$daily$delta
      se <- sqrt(var_pub)

      RES[[length(RES)+1L]] <- data.table(
        channel = chan, spec_label = key, form = spec, df = deg,
        winsor_lo = wl, winsor_hi = wu, origine = cfg$origine,
        is_ref = (key == GRID_REF_SPEC), status = "OK",
        delta_daily = PACK$daily$delta, delta_annual = PACK$annual$delta,
        diff = D, se_diff = se, ci_lo = D - 1.96*se, ci_hi = D + 1.96*se,
        excludes_zero = (D - 1.96*se) * (D + 1.96*se) > 0,
        rho = rho, rho_raw = rho0, cov_da = cov_da,
        se_daily = sqrt(v_d), se_annual = sqrt(v_a),
        se_diff_naive = sqrt(v_d + v_a),          # ce qu'on aurait ecrit sans le bloc croise
        var_diff_raw = var_raw, var_diff_psd = var_pub, var_diff_nossc = var_nossc,
        psd_fix = psd_fix, n_eigen_trunc = ntr, min_eigen = minev,
        ssc_daily = OWN$daily$ssc, ssc_annual = OWN$annual$ssc, ssc_cross = ssc_x,
        ssc_daily_analytic = OWN$daily$ssc_ana, ssc_annual_analytic = OWN$annual$ssc_ana,
        n_obs_daily = PACK$daily$nobs, n_obs_annual = PACK$annual$nobs,
        n_obs_common = length(common), same_sample = same_sample,
        nparams_daily = PACK$daily$nparams, nparams_annual = PACK$annual$nparams,
        G_orig = OWN$daily$G_or, G_dest = OWN$daily$G_de,
        oos_daily = PACK$daily$oos, oos_annual = PACK$annual$oos,
        support_lo_daily = PACK$daily$lo, support_hi_daily = PACK$daily$hi,
        support_lo_annual = PACK$annual$lo, support_hi_annual = PACK$annual$hi,
        aic_daily = PACK$daily$aic, aic_annual = PACK$annual$aic,
        b1_max_relnorm = B1max, b1_max_spread = B1spr, n_cells_north = ncell,
        cluster_unit = "orig x dest (two-way CGM)",
        secs = as.numeric((proc.time()-tk)["elapsed"]))

      # ---------- persistance des scores agreges (demande explicite) ----------
      # Les viandes k x k suffisent a reconstruire n'importe quelle vcov de ces deux
      # bras. Les scores agreges par PAIRE, plus lourds, permettent en outre de
      # re-clusteriser autrement ; on ne les garde que pour la reference (P4).
      st <- list(channel = chan, spec_label = key, cc_daily = PACK$daily$cc,
                 cc_annual = PACK$annual$cc,
                 meat_daily = OWN$daily$meat, meat_annual = OWN$annual$meat,
                 meat_cross = M_cross,
                 bread_daily = OWN$daily$bread, bread_annual = OWN$annual$bread,
                 w_daily = PACK$daily$wfull, w_annual = PACK$annual$wfull,
                 ssc = c(daily = OWN$daily$ssc, annual = OWN$annual$ssc, cross = ssc_x))
      if (key == GRID_REF_SPEC) {
        ag <- function(S, g) rowsum(S, group = g, reorder = TRUE)
        st$scores_pair_daily  <- ag(as.matrix(ARM$daily$mod$scores)[pos_d, , drop=FALSE], g_pair)
        st$scores_pair_annual <- ag(as.matrix(ARM$annual$mod$scores)[pos_a, , drop=FALSE], g_pair)
      }
      STORE[[paste(chan, key, sep = "|")]] <- st

      cat(sprintf("     D = %+.4f  se = %.4f  IC [%+.4f, %+.4f]  rho = %.3f  %s  (%.0fs)\n",
                  D, se, D - 1.96*se, D + 1.96*se, rho,
                  if ((D-1.96*se)*(D+1.96*se) > 0) "EXCLUT ZERO" else "couvre zero",
                  (proc.time()-tk)["elapsed"]))
      rm(ARM, OWN, PACK, Vj, Vj0, Vpub, V_cross, M_cross); gc(verbose = FALSE)
    }
  }
  rm(dat, clim); gc(verbose = FALSE)
}

R <- rbindlist(RES, fill = TRUE)

##################################################################################
# 6. ASSERTION B2 — les Delta doivent reproduire le Test 1
##################################################################################

cat("\n\n--- B2 : reproduction des Delta_north du Test 1 ---\n")
t1 <- fread(T1_CSV)[status == "OK",
        .(channel, arm, spec_label, delta_ref = delta_north)]
mine <- melt(R[status == "OK", .(channel, spec_label, daily = delta_daily,
                                 annual = delta_annual)],
             id.vars = c("channel", "spec_label"), variable.name = "arm",
             value.name = "delta_mine")
mine[, arm := as.character(arm)]
cmp <- merge(mine, t1, by = c("channel", "arm", "spec_label"))
cat(sprintf("  paires comparees : %d\n", nrow(cmp)))
b2 <- if (nrow(cmp)) max(abs(cmp$delta_mine - cmp$delta_ref)) else NA_real_
cat(sprintf("  max|ecart| = %.3e (tolerance 1e-9)\n", b2))
if (nrow(cmp) && (!is.finite(b2) || b2 > 1e-9)) {
  print(cmp[abs(delta_mine - delta_ref) > 1e-9])
  stop(sprintf("B2 VIOLEE : ecart max %.3e. La geographie re-derivee a diverge du Test 1.", b2))
}
if (!nrow(cmp)) {
  cat("  ⚠ aucune paire comparable — B2 NON EVALUEE (le garde-fou ne s'est pas execute)\n")
} else {
  cat(sprintf("  ✓ B2 : la re-derivation reproduit le Test 1 sur %d paires\n", nrow(cmp)))
}
R[, b2_max_diff := b2]

##################################################################################
# 7. SORTIE
##################################################################################

R[, md5_inputs := paste(sprintf("%s=%s", names(OBSERVED_MD5), OBSERVED_MD5), collapse = "; ")]
R[, inputs_fingerprint := INPUTS_FP]
R[, generated := "ct2_between_arm_test.R"]
R[, ref_spec := GRID_REF_SPEC]
setcolorder(R, c("channel","spec_label","form","df","winsor_lo","winsor_hi","origine",
                 "is_ref","status","delta_daily","delta_annual","diff","se_diff",
                 "ci_lo","ci_hi","excludes_zero","rho","se_daily","se_annual",
                 "se_diff_naive","psd_fix","same_sample","n_obs_common"))
fwrite(R, OUT_CSV);       cat("\nEcrit:", OUT_CSV, "\n")
saveRDS(STORE, OUT_RDS, compress = "xz")
cat(sprintf("Ecrit: %s (%.1f Mo)\n", OUT_RDS, file.size(OUT_RDS)/1024^2))

cat("\n=============== TABLE ===============\n")
print(R[status == "OK", .(channel, spec_label, ref = is_ref,
                          d_daily = round(delta_daily, 4), d_annual = round(delta_annual, 4),
                          diff = round(diff, 4), se = round(se_diff, 4),
                          ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4),
                          excl0 = excludes_zero, rho = round(rho, 3),
                          se_naif = round(se_diff_naive, 4), psd = psd_fix)], nrows = 100)

cat("\n--- ETENDUE SUR LA GRILLE (P3 : jamais de ponderation AIC) ---\n")
print(R[status == "OK", .(
  n = .N, n_excl0 = sum(excludes_zero),
  diff_min = round(min(diff), 4), diff_max = round(max(diff), 4),
  se_min = round(min(se_diff), 4), se_max = round(max(se_diff), 4),
  rho_min = round(min(rho), 3), rho_max = round(max(rho), 3),
  gain_vs_naif = round(mean(1 - se_diff/se_diff_naive), 3)), by = channel])

cat("\n--- SPECIFICATION DE REFERENCE (P4) ---\n")
print(R[is_ref == TRUE & status == "OK",
        .(channel, spec_label, diff = round(diff, 4), se = round(se_diff, 4),
          ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4),
          excl0 = excludes_zero, rho = round(rho, 3))])

cat("\nTemps total:", round((proc.time()-t0)["elapsed"]/60, 1), "min\n=== Done ===\n")
