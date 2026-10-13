##################################################################################
# 00_parameters.R -- exogenous parameters of the dose-response estimations
#
# Every parameter that the paper does NOT estimate lives here, with its source.
# No other script in the package should hard-code one.
#
# Convention: a parameter borrowed from the literature carries its citation; a
# construction parameter carries its reason. A parameter with neither is a
# defect, not a choice.
#
# Reads:  nothing (definitions only)
# Writes: creates intermediate_outputs/ and output/ if absent
# Sourced by every estimation and GE script in the package.
##################################################################################

##################################################################################
# PACKAGE ROOT -- DERIVED FROM THE LOCATION OF THIS FILE
#
# This file lives in <root>/01_estimation/, so the root is the parent of its own
# directory. It is measured, not declared: a hard-coded path, the working
# directory and a required environment variable can all be wrong, whereas the
# location of the file cannot.
#
# JMP_ROOT takes precedence when set. The shell orchestrator sets it explicitly,
# and when it does, the deduction is not needed.
##################################################################################

.jmp_find_self <- function() {
  ## 1) source(). The innermost source() frame carries the path of THIS file.
  ##    The search runs backwards to catch it before the calling script's frame,
  ##    which would otherwise return 01a_productivity_drf.R and make the derived
  ##    root correct only by accident.
  fr <- sys.frames()
  if (length(fr))
    for (i in rev(seq_along(fr))) {
      f <- fr[[i]]$ofile
      if (!is.null(f) && nzchar(f))
        return(normalizePath(f, winslash = "/", mustWork = FALSE))
    }
  ## 2) Rscript <file>. The path is in --file=.
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m))
    return(normalizePath(sub("^--file=", "", m[1]), winslash = "/", mustWork = FALSE))
  NA_character_
}

.jmp_root   <- Sys.getenv("JMP_ROOT", unset = "")
.jmp_origine <- "JMP_ROOT (environment variable)"

if (!nzchar(.jmp_root)) {
  .jmp_self <- .jmp_find_self()
  if (is.na(.jmp_self))
    stop("00_parameters.R: cannot determine the root.\n",
         "  Neither JMP_ROOT set, nor a readable file location.\n",
         "  Known case: code pasted into the console rather than source().\n",
         "  Set JMP_ROOT, or source this file by its path.")
  .jmp_root    <- dirname(dirname(.jmp_self))   # <root>/01_estimation/00_parameters.R
  .jmp_origine <- paste0("derived from ", .jmp_self)
}

## Sanity check on the derived root. An empty root makes file.path("",
## "01_estimation", ...) return "/01_estimation/...", an absolute path from the
## disk root. source() fails loudly there, but fread() or dir.create() would
## silently look in, or create, C:/01_estimation/. This check runs BEFORE any
## dir.create: a wrong root must not create anything.
.jmp_attendus <- c("01_estimation", "data_inputs")
.jmp_manque   <- .jmp_attendus[!dir.exists(file.path(.jmp_root, .jmp_attendus))]
if (length(.jmp_manque))
  stop("00_parameters.R: the root does not look like the repository.\n",
       "  root      : ", .jmp_root, "\n",
       "  origin    : ", .jmp_origine, "\n",
       "  missing   : ", paste(.jmp_manque, collapse = ", "), "\n",
       "  No directory was created. Fix the root before running again.")

message("[00_parameters] root = ", .jmp_root, "\n",
        "                origin = ", .jmp_origine)

DIR_IN  <- file.path(.jmp_root, "data_inputs")
DIR_TMP <- file.path(.jmp_root, "intermediate_outputs")

## ---------------------------------------------------------------------------
## DIR_OUT, ET SA REDIRECTION POUR UN REJEU EN PARALLELE
## ---------------------------------------------------------------------------
## Par defaut : <racine>/output. C'est le comportement de toujours, et une
## variable non posee ne doit RIEN changer -- douze scripts sourcent ce fichier.
##
## JMP_OUTPUT_DIR redirige les ECRITURES vers un autre dossier. Elle sert a
## rejouer un script d'estimation cote a cote avec l'etat publie, pour comparer
## table par table AVANT de remplacer quoi que ce soit. Sans elle, verifier la
## reproductibilite d'un script de deux heures obligerait a ecraser d'abord et
## a comparer ensuite -- c'est-a-dire a ne plus avoir de point de comparaison.
##
## Un chemin RELATIF se lit depuis la racine du depot, jamais depuis le
## repertoire courant : R peut etre lance de n'importe ou, et un rejeu qui
## ecrirait dans un « output_01c_rerun » cree a cote du terminal serait
## introuvable au matin.
.jmp_outdir <- Sys.getenv("JMP_OUTPUT_DIR")
if (nzchar(.jmp_outdir)) {
  DIR_OUT <- if (grepl("^([A-Za-z]:|/|\\\\)", .jmp_outdir)) .jmp_outdir
             else file.path(.jmp_root, .jmp_outdir)
  message("[00_parameters] DIR_OUT REDIRECTED (JMP_OUTPUT_DIR) = ", DIR_OUT,
          "\n                output/ is NOT touched by this run.")
} else {
  DIR_OUT <- file.path(.jmp_root, "output")
}
DIR_CRH  <- file.path(DIR_TMP, "crh")
DIR_EGGW <- file.path(DIR_IN, "EGGW_Replication_Package")

## ---------------------------------------------------------------------------
## THETA_GE — L'ELASTICITE COMMERCIALE DU MODELE, LUE A SA SOURCE UNIQUE.
## Posee le 19/09/2026.
## ---------------------------------------------------------------------------
## CE QU'ELLE CONVERTIT. Les flux commerciaux sont une loupe : une baisse de 1 %
## de la productivite fondamentale deplace les exportations d'environ theta %.
## Ce que le PPML de l'etage 3 identifie est donc theta fois la reponse de
## productivite -- ce que l'eq. (5) du papier ecrit correctement. Le GE, lui,
## reapplique theta dans trade_weights = (Z/x)^theta. Sans division, le choc
## etait amplifie DEUX FOIS et toutes les magnitudes GE valaient ~6,5 fois trop.
##
## POURQUOI CETTE VALEUR ET PAS UNE AUTRE. Les theta_k sectoriels de
## Caliendo-Parro s'etalent de 0,37 a 51,08 (moyenne 6,95, ponderee 7,08,
## harmonique 2,67) : aucun scalaire tire de cette distribution ne serait
## justifiable. Le diviseur est EPINGLE PAR LE MODELE -- trade_weights eleve a
## param_theta, et la chaine n'est coherente que pour cette valeur exacte.
## (La ponderation sectorielle scale = 1/theta_k a ete testee et rejetee :
##  +105 325 d'AIC et perte de la significativite jointe. Voir RECAP 3.149.)
##
## ⛔ LUE, JAMAIS REECRITE EN DUR. La valeur vit dans le fichier de configuration
## du solveur. La recopier ici creerait une seconde source de verite pour la
## constante dont le mauvais usage EST le defaut qu'on corrige.
## ⛔ AUCUN REPLI. Fichier introuvable ou cle absente => ARRET. Une constante de
## conversion qui se degrade en silence est exactement le mode de defaillance
## qu'on ferme.
.theta_file <- file.path(.jmp_root, "02_GE", "02c_model", "specification_study.csv")
if (!file.exists(.theta_file))
  stop("[FATAL] specification_study.csv not found: ", .theta_file,
       "\n  THETA_GE cannot be read, and there is no default value.")
.theta_raw <- readLines(.theta_file, warn = FALSE)
.theta_hit <- grep("^\\s*param_theta\\s*;", .theta_raw, value = TRUE)
if (length(.theta_hit) != 1L)
  stop("[FATAL] key 'param_theta' missing or repeated in ", .theta_file,
       " (", length(.theta_hit), " occurrence(s)). No fallback.")
THETA_GE <- suppressWarnings(as.numeric(strsplit(.theta_hit, ";")[[1]][2]))
if (!is.finite(THETA_GE) || THETA_GE <= 0)
  stop("[FATAL] param_theta unreadable or not positive in ", .theta_file,
       ": '", .theta_hit, "'")
rm(.theta_file, .theta_raw, .theta_hit)

## SCALE_PROD et SCALE_AMEN sont poses PLUS BAS, juste apres Omega : Omega est
## defini a la section des parametres de migration et n'existe pas encore ici.

## ---------------------------------------------------------------------------
## DIR_TMP_OUT — LES ECRITURES INTERMEDIAIRES SUIVENT LA REDIRECTION, PAS LES
## LECTURES. Pose le 19/09/2026.
## ---------------------------------------------------------------------------
## 🚨 LE DEFAUT QUE CECI CORRIGE. JMP_OUTPUT_DIR promet qu'un run de sonde ne
## touche pas l'etat vivant, et le message ci-dessus l'affirme : « output/ n'est
## PAS touche par ce run ». C'etait vrai de output/ et FAUX de
## intermediate_outputs/ : DIR_TMP ne suivait pas la redirection, donc un 01a
## redirige REECRIVAIT BACI_raw_pairs_annual.rds et climate_long_annual.rds --
## deux fichiers dont le md5 entre dans l'empreinte de provenance de ct1, ct2 et
## ct3. Une sonde pouvait donc invalider la chaine du contraste en silence.
## Mesure du 19/09 : les deux ont ete reecrits a l'identique au bit, donc
## l'exposition etait reelle et sans degat CETTE FOIS. C'est la seule fois ou on
## peut le dire, et ce n'est pas une garantie.
##
## ⛔ POURQUOI PAS « DIR_TMP SUIT JMP_OUTPUT_DIR », QUI TIENDRAIT EN UNE LIGNE.
## Recensement du 19/09 : intermediate_outputs est LU par sept scripts (01a,
## 01b, 01c, 01d_clamp, 01e, step2_maps, step7) et ECRIT par deux seulement
## (01a, 01b). Faire suivre la redirection a DIR_TMP casserait les sept
## lectures : sous redirection, 01a chercherait ERA_tas_annual.csv dans un
## dossier vide et s'arreterait, et 01c perdrait ses quatre intrants. La
## redirection naive serait une regression immediate, pas un correctif.
##
## ⇒ PARTAGE LECTURE / ECRITURE. DIR_TMP reste le chemin CANONIQUE en lecture.
##   DIR_TMP_OUT suit la redirection et ne sert qu'aux ECRITURES. Sans
##   JMP_OUTPUT_DIR les deux coincident, donc rien ne change par defaut.
if (nzchar(.jmp_outdir)) {
  DIR_TMP_OUT <- file.path(DIR_OUT, "intermediate_outputs")
  message("[00_parameters] DIR_TMP_OUT REDIRECTED = ", DIR_TMP_OUT,
          "\n                intermediate_outputs/ is read, never written, ",
          "by this run.")
} else {
  DIR_TMP_OUT <- DIR_TMP
}

for (d in c(DIR_TMP, DIR_OUT, DIR_TMP_OUT))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## Cruz, J. L. & Rossi-Hansberg, E., "The Economic Geography of Global
## Warming", replication package, Zenodo, DOI 10.5281/zenodo.7681424.
## Supplies bilateral trade costs, the (x, y) cell identifier, and the
## framework replicated in the robustness appendix.
CRH_ZENODO_ID    <- "7681424"
CRH_OUTER_ZIP    <- file.path(DIR_IN, "EGGW Replication Package Final.zip")
CRH_OUTER_MD5    <- "343ad34579940a46c4b71b5eda3bb422"
CRH_INNER_NAME   <- "EGGW_Replication_Package.zip"
EGGW_ZIP <- file.path(DIR_IN, CRH_INNER_NAME)
EGGW_PREFIX      <- "EGGW_Replication_Package"

## LA TRIANGULATION CRH -- constantes partagees (remontees ici le 18/09/2026).
## Elles etaient declarees dans 01d_crh_triangulation.R, son unique lecteur a
## l'epoque. Le lot 10 leur en donne deux autres -- 01d_crh_clamp_sensitivity.R
## et 01d_fig_crh_responses.R -- et une latitude de coupure recopiee dans trois
## scripts est la forme habituelle des divergences silencieuses de ce paquet.
##   CRH_FUT_YEAR reste 2095 : c'est le dernier bloc PLEIN {2095-2099}, que le
##   manuscrit etiquette 2100 par HORIZON_LABEL. Ce n'est pas une etiquette,
##   c'est l'indice que le modele evalue -- voir tools/check_horizon.py.
CRH_BASE_YEAR <- 2015L
CRH_FUT_YEAR  <- 2095L
CRH_COLD      <- c("RUS", "CAN", "FIN", "DEU", "FRA")
CRH_NORTH_LAT <- 45
## Point de clamp de robustesse sur delta(T) : la fonction est tenue plate
## au-dela. Il n'y a PAS de clamp par defaut dans la projection publiee
## (reliable_hi = Inf) ; ces deux valeurs sont celles des deux sensibilites.
##   40 degC = la borne du support bien peuple (regle du 1 %, k_max) ;
##   30 degC = la borne plus severe demandee au lot 10.
CRH_CLAMP_SUPPORT <- 40
CRH_CLAMP_SEVERE  <- 30

##################################################################################
# AMENITY CHANNEL -- bilateral migration
##################################################################################

## Congestion elasticity and dispersion parameter, both taken from
## Desmet & Rossi-Hansberg (2018).
Lambda <- 0.32   # congestion elasticity
Omega  <- 0.50   # dispersion parameter

## ---------------------------------------------------------------------------
## L'ECHELLE DE BASE PAR CANAL — LES DEUX COTE A COTE, ET C'EST VOULU.
## ---------------------------------------------------------------------------
## prepare_diffs() divise le differentiel climatique par `scale`, donc le
## coefficient estime sort MULTIPLIE par `scale`. Les deux canaux convertissent
## vers leurs unites FONDAMENTALES par le MEME mecanisme, au MEME endroit :
##   amenite      : l'equation porte f_a'/Omega   => scale = Omega        (01b l.267)
##   productivite : l'equation porte theta * f_z  => scale = 1/THETA_GE   (19/09/2026)
##
## ⚠️ C'EST L'ASYMETRIE ENTRE CES DEUX LIGNES QUI A LAISSE VIVRE LE DEFAUT.
## Le canal amenite faisait le geste depuis toujours ; le canal productivite
## posait `scale = 1`, commente « aucune mise a l'echelle », et le GE
## reappliquait theta. Les deux se lisent desormais l'une sous l'autre : un
## lecteur qui trouve l'une trouve l'autre, et une divergence future se voit.
##
## ⛔ NE PAS CONVERTIR AILLEURS. La division se fait ICI et nulle part ailleurs :
## ni a l'export de la courbe, ni a l'emission d'une macro. Une seconde
## conversion s'appliquerait en cascade et serait silencieuse. Voir RECAP 3.149.
SCALE_PROD <- 1 / THETA_GE
SCALE_AMEN <- Omega
message("[00_parameters] THETA_GE = ", THETA_GE,
        " (read from 02c_model/specification_study.csv)",
        " | SCALE_PROD = 1/", THETA_GE, " = ", signif(SCALE_PROD, 8),
        " | SCALE_AMEN = Omega = ", Omega)

## Abel & Cohen quinquennial panel: Stage 3 is estimated from 1990 onward.
YEAR_MIN_STAGE3  <- 1990L
PERIODS_EXPECTED <- c(1985L, 1990L, 1995L, 2000L, 2005L, 2010L, 2015L)

##################################################################################
# PRODUCTIVITY CHANNEL -- bilateral trade
##################################################################################

## Sectoral theta_k values come from Caliendo & Parro (2015) and enter through a
## file (data_inputs/caliendo_parro_2015_thetak/theta_k.csv). They are not
## constants and must not be copied here.

## Annual panel: GLORIA and ERA5 both cover 1995-2019.
TRADE_YEAR_MIN <- 1995L
TRADE_YEAR_MAX <- 2019L
TRADE_N_YEARS  <- TRADE_YEAR_MAX - TRADE_YEAR_MIN + 1L   # complete pairs only

## Excluded sectors.
##   extractive: no theta_k in Caliendo-Parro, and no TFP channel;
##   services:   X_rs/X_ss undefined in the GLORIA final-demand matrix.
ISIC_EXTRACTIVE <- c("02", "05", "10", "11", "12", "13", "14")
ISIC_SERVICES   <- c("40", "41", "45", "50", "51", "52", "55",
                     "60", "61", "62", "63", "64", "65", "66", "67",
                     "70", "71", "72", "73", "74", "75",
                     "80", "85", "90", "91", "92", "93", "95")

## Cap on the X_rs/X_ss ratio. GLORIA final demand understates domestic
## absorption in manufacturing, where intermediate use dominates final demand,
## which produces ratios above one. The series is cut at the 99th percentile
## with a hard ceiling at 20.
DV_QUANTILE_CAP <- 0.99
DV_HARD_CAP     <- 20

## Reference specification of the identification diagnostics (appendix C).
REF_SPEC   <- "polynomial"
REF_DEG    <- 3L
REF_WLOWER <- 1
REF_WUPPER <- 1

## TACHE B DE 01c (controles IDH, section R4) : SAUTEE PAR DECLARATION.
## Author's decision, 28/09/2026 (lot 56, option A). La tache B n'a JAMAIS
## tourne dans le paquet : son lecteur (01c l.~627) attend `GID_0` et `hdi` alors
## que oh_hdi_country_year.csv porte `iso3,year,HDI` ; `hdi` y resout vers la table
## elle-meme ; et il ne reporte pas 2015 sur 2016-2019 comme robustness_oh.R. L'etat
## publie (\EMPohSchemes*) vient des taches A et C. L'executer pour la premiere fois
## au milieu de la reprise des poids melangerait deux changements. Elle ne tourne
## pas parce que CE PARAMETRE le dit, et non plus parce qu'un fichier manque.
## A rouvrir comme changement ISOLE (LISTE_POST_SOUMISSION n° 32).
OH_TACHE_B_HDI <- FALSE

##################################################################################
# COMMON TO BOTH CHANNELS
##################################################################################

## 4 functional forms x 3 winsorisations = 12 specifications per arm.
SPEC_CONFIGS <- list(
  list(spec = "polynomial", degree = 2),
  list(spec = "polynomial", degree = 3),
  list(spec = "spline",     degree = 3),
  list(spec = "spline",     degree = 4)
)

WINSOR_CONFIGS <- list(
  list(lower = 0.5, upper = 0.5),
  list(lower = 1.0, upper = 1.0),
  list(lower = 2.5, upper = 2.5)
)

##################################################################################
# GRID_SPECS -- LA GRILLE DE CONTRASTE (contribution A), DEFINIE ICI ET NULLE PART
# AILLEURS
#
# DECISION D'AUTEUR DU 2026-09-18 : la grille passe de VINGT-QUATRE a QUINZE
# specifications. Les splines naturelles a 5, 6 et 7 degres de liberte sont
# retirees -- elles oscillent dans la partie DENSE du support (20-28 degC), ou la
# masse-jour est la plus forte et ou la revendication du papier se joue.
# Ce qui reste : polynomes de degre 2 a 4, splines naturelles a 3 et 4 df,
# soit 5 formes x 3 winsorisations = 15.
#
# ⚠️ DEUX GRILLES COEXISTENT DANS CE PAQUET, ET ELLES NE SE CONFONDENT PAS.
#   SPEC_CONFIGS (ci-dessus, 4 formes x 3 winsor. = 12) : la grille du CORPS,
#     celle qu'estiment 01a, 01b et 01c, celle des poids d'Akaike, celle dont
#     descendent le GE, la decomposition et l'annexe CRH. ELLE N'EST PAS
#     TOUCHEE par la decision du 18/09 : \EMPgridNspec vaut toujours 12.
#   GRID_SPECS (ici, 5 formes x 3 winsor. = 15) : la grille du CONTRASTE
#     annuel/journalier, celle de ct1, ct2, ct4 et ct6, celle que compte
#     \EMPgridProdNspec et que la Table du supplement imprime ligne a ligne.
# La premiere est INCLUSE dans la seconde (12 sur 15 : poly 4 est l'ajout).
# Ecrire « la grille » sans dire laquelle est la faute que ce commentaire evite.
#
# 🔑 POURQUOI ICI PLUTOT QUE DANS ct1. Jusqu'au 18/09, ct1 et ct2 portaient
# CHACUN leur copie de la liste des huit formes, avec un commentaire disant que
# l'ecrasement de SPEC_CONFIGS etait delibere. Deux copies d'une meme decision
# sont deux endroits ou la changer, et l'une des deux aurait fini par ne pas
# l'etre : ct4 et ct6, eux, ne comptaient les specifications que sur le CSV
# recu, donc une grille reduite d'un cote leur serait passee inapercue. Depuis
# le 18/09 la liste vit ici, les quatre scripts la lisent, et tools/check_macros.py
# (controle 9) echoue si un producteur en reintroduit une en dur.
##################################################################################

## Les FORMES. `origine` marque les quatre formes du papier publie (df <= 4 hors
## poly 4) : ct1 s'en sert pour savoir quelles courbes entrent dans l'assertion
## A2, et cette colonne est ecrite telle quelle dans test1_delta_north.csv.
GRID_FORMS <- list(
  list(spec = "polynomial", degree = 2L, origine = TRUE),
  list(spec = "polynomial", degree = 3L, origine = TRUE),
  list(spec = "polynomial", degree = 4L, origine = FALSE),
  list(spec = "spline",     degree = 3L, origine = TRUE),
  list(spec = "spline",     degree = 4L, origine = TRUE)
)

## Les WINSORISATIONS. Les memes trois seuils que la grille du corps : la
## decision du 18/09 ne porte que sur les formes. On les DERIVE de
## WINSOR_CONFIGS au lieu de les recopier -- une seconde liste de seuils serait
## exactement la duplication que ce bloc supprime.
GRID_WINSOR <- WINSOR_CONFIGS

## L'ETIQUETTE, ET C'EST LA MEME FONCTION QUE CELLE DES PRODUCTEURS.
## Format historique de ct1 : sprintf("%s%d_wl%gwu%g", ...) -> "spline4_wl0.5wu0.5".
## `%g` rend 0.5, 1 et 2.5 : ne pas passer a %s, qui rendrait "1.0" et casserait
## la jointure avec tous les artefacts deja ecrits.
grid_spec_label <- function(spec, degree, w_lo, w_hi)
  sprintf("%s%d_wl%gwu%g", spec, as.integer(degree), w_lo, w_hi)

## LA GRILLE, A PLAT : une entree par specification, etiquette comprise.
## Produit croise DERIVE, jamais reecrit a la main -- 15 lignes recopiees sont
## 15 occasions d'une faute de frappe qu'aucune assertion ne verrait.
GRID_SPECS <- unlist(lapply(GRID_FORMS, function(f)
  lapply(GRID_WINSOR, function(w)
    list(spec = f$spec, degree = f$degree, origine = f$origine,
         winsor_lo = w$lower, winsor_hi = w$upper,
         spec_label = grid_spec_label(f$spec, f$degree, w$lower, w$upper)))),
  recursive = FALSE)

GRID_NSPEC  <- length(GRID_SPECS)
GRID_LABELS <- vapply(GRID_SPECS, function(z) z$spec_label, character(1))

## LA SPECIFICATION DE REFERENCE DU CONTRASTE, DESIGNEE AVANT TOUT RESULTAT.
## C'est le pre-engagement (P4) de ct2_between_arm_test.R : ns df = 4,
## winsorisation 0,5 %. Motif entierement anterieur aux chiffres -- c'est la
## specification du papier publie, et 0,5 % est le reglage le moins ecretant des
## trois. Elle n'est PAS choisie parce qu'un resultat l'aurait recommandee.
## Elle vit ici depuis le 18/09/2026 parce que ct2 (R) et ct3 (Python) en
## portaient chacun une copie en dur : deux copies d'une designation qui doit
## etre unique, et dans deux langages, donc introuvables l'une par l'autre.
## ⚠️ HOMONYMIE ASSUMEE ET SIGNALEE : REF_SPEC, plus haut, vaut "polynomial" et
## designe la specification de reference des DIAGNOSTICS D'IDENTIFICATION
## (annexe C), un autre objet. Les deux noms different pour cette raison.
GRID_REF_SPEC <- grid_spec_label("spline", 4L, 0.5, 0.5)

## Les quatre proprietes que les consommateurs supposent, verifiees ICI, une
## fois, au chargement -- pas dans chacun des quatre scripts.
##   (1) le compte est bien le produit croise ;
##   (2) les etiquettes sont DISTINCTES : deux identiques feraient une jointure
##       ct1 x ct2 silencieusement fausse, avec le bon nombre de lignes ;
##   (3) les formes sont dans {polynomial, spline} : ct4 derive « Poly »/« NS »
##       de cette colonne et rendrait un libelle faux pour toute autre ;
##   (4) la winsorisation est SYMETRIQUE : la table n'imprime qu'une borne.
stopifnot(GRID_NSPEC == length(GRID_FORMS) * length(GRID_WINSOR),
          !anyDuplicated(GRID_LABELS),
          GRID_REF_SPEC %in% GRID_LABELS,
          all(vapply(GRID_SPECS, function(z) z$spec %in% c("polynomial", "spline"), TRUE)),
          all(vapply(GRID_SPECS, function(z) abs(z$winsor_lo - z$winsor_hi) < 1e-12, TRUE)))

message("[00_parameters] GRID_SPECS = ", GRID_NSPEC, " specifications (",
        length(GRID_FORMS), " forms x ", length(GRID_WINSOR), " winsorisations)")

## ---- Controls ------------------------------------------------------------------
GRAVITY_CONTROLS  <- c("distw_harmonic", "distcap", "comlang_off",
                       "col_dep_ever", "fta_wto")
DISASTER_CONTROLS <- c("nb_events", "total_duration_days", "total_affected")
## Productivity: the Stage 3 controls include the structural income term.
TRADE_STAGE3_CONTROLS   <- c(GRAVITY_CONTROLS, "GDP_diff")
## Amenity: the two structural terms only; the rest is absorbed by the pair.
AMENITY_STAGE3_CONTROLS <- c("GDP_diff", "populationdiff")

## ---- Monte Carlo ---------------------------------------------------------------
MC_N_SIM  <- 300L                    # draws per arm

## MC_N_SIM_OH — LE BLOC OH SEUL, ET SA VALEUR EST MESUREE, PAS CHOISIE.
##
## Mesure du 2026-09-07 (logs/nsim_*/exp_nsim.csv) : les memes ajustements,
## deux graines de Monte-Carlo, deux valeurs de N.
##   N =  300 : ecart median entre graines sur q95 = 1,043 ; les deux graines
##              ne s'accordent sur le verdict d'identification que pour 2 des
##              5 schemas. Zero tombe DANS l'ecart : « identifie » n'est alors
##              pas une propriete des donnees, c'est le tirage qui decide.
##   N = 3000 : ecart median 0,142 (facteur ~7), et les deux graines
##              s'accordent sur les 5 schemas.
## Le verdict lui-meme change : 3 schemas identifies sur 5 a N = 3000, contre
## 0 sur 5 a N = 300. Le « zero identifie » etait un artefact du nombre de
## tirages. Aucun critere, aucun seuil, aucune graine n'a ete ajuste.
##
## ⛔ POURQUOI PAS GLOBALEMENT. MC_N_SIM gouverne six sites : 01a x2, 01b x2,
## 01c x2. Le monter partout reecrirait MC_damage_{prod,amenity}_*_quantiles.csv,
## que le GE CONSOMME — il faudrait refaire les onze runs (1 h 38 mesurees) — et
## deplacerait les eventails publies de 01a et 01b, qui ne sont pas en cause.
## ⚠️ Ce qui N'A PAS ete mesure : si N = 300 suffit pour 01a et 01b. L'instabilite
## constatee porte sur une statistique de QUEUE au voisinage d'un seuil ; une
## enveloppe de courbe est un autre objet. La question se teste de la meme facon
## et reste ouverte.
##
## COUT : 5 appels mc_run dans 01c (3 en R3, 2 en R4) a ~6 min contre ~1 min,
## soit 01c de ~25 min a ~50 min. mc_run est SOUS-lineaire en N : mesure
## x5,07 pour x10 tirages (70,9 s -> 359,2 s).
MC_N_SIM_OH <- 3000L                 # bloc OH de 01c (R3, R4) UNIQUEMENT
MC_X_GRID <- seq(-5, 35, by = 0.5)   # evaluation grid (deg C)
## Seeds are fixed and DISTINCT across arms: the Monte Carlo is reproducible to
## the bit, and two arms of the same channel do not share a draw sequence.
MC_SEED <- c(prod_daily    = 44L, prod_annual    = 45L,
             amenity_daily = 42L, amenity_annual = 43L)
## LOT 68 (30/09/2026) -- la graine du bras productivite journalier est
## SURCHARGEABLE par l'environnement, comme JMP_OUTPUT_DIR. Seule l'etape
## 01c_graine2 de run_all.sh la pose, et toujours avec JMP_OUTPUT_DIR : un rejeu
## a une autre graine n'ecrit jamais dans output/. On refuse la surcharge sans
## redirection plutot que de laisser une seconde graine ecraser la premiere.
.jmp_seed_env <- Sys.getenv("JMP_MC_SEED_PROD_DAILY")
if (nzchar(.jmp_seed_env)) {
  if (!grepl("^[0-9]+$", .jmp_seed_env))
    stop("[STOP] JMP_MC_SEED_PROD_DAILY must be an integer: ", .jmp_seed_env)
  if (!nzchar(Sys.getenv("JMP_OUTPUT_DIR")))
    stop("[STOP] JMP_MC_SEED_PROD_DAILY set without JMP_OUTPUT_DIR: the second seed ",
         "would write into output/.")
  MC_SEED[["prod_daily"]] <- as.integer(.jmp_seed_env)
  message("[00_parameters] MC_SEED prod_daily OVERRIDDEN (JMP_MC_SEED_PROD_DAILY) = ",
          MC_SEED[["prod_daily"]])
}
MC_EXPORT_GRID <- -100L:100L         # grid of the CSV handed to the GE (constant extension)

## ---- Bande GE « arm B » : LE NOMBRE DE TIRAGES, DECLARE UNE SEULE FOIS ---------
## Pose le 20/09/2026, apres la nuit perdue du 19.
##
## CE QUE CE BLOC REMPLACE, ET POURQUOI IL EXISTE.
## ct5 lisait son nombre de tirages dans une VARIABLE D'ENVIRONNEMENT AVEC
## DEFAUT (`ARMB_EXTRA`, defaut FALSE). Le pilote de la nuit du 19 l'a appelee
## nue : ct5 a rendu trente colonnes rb, le pilote en a demande soixante, et la
## passe est morte au deuxieme run sur une colonne absente -- neuf heures plus
## tard. Pire que la panne : avec trente colonnes, TOUS les gardes en aval
## passaient, parce que l'emetteur lit la taille attendue DANS le meme artefact
## que ct5 vient de retrecir.
##
## N N'EST PLUS UNE VARIABLE D'ENVIRONNEMENT. Il est ici, il n'a pas de repli,
## et les deux consommateurs -- ct5 qui ECRIT la grille, emit_values.py qui la
## LIT -- le comparent a ce qu'ils voient. Le pilote ne le declare pas non plus :
## il le lit dans ce fichier et le passe en poignee de main a ct5.
##
## POURQUOI SOIXANTE ET PAS TRENTE. C'est un choix methodologique publie : la
## granularite de la part negative par unite est de 1/N. Passer a trente
## deplacerait le chiffre lui-meme, ce qui serait changer un resultat pour
## eviter un cout.
ARMB_BASE_N <- 30L    # la grille de BASE : percentiles 5..95 equirepartis
ARMB_N      <- 60L    # la grille PUBLIEE : la base, plus des INTERCALAIRES
## La grille etendue n'est pas une grille reguliere de soixante points : c'est
## la base de trente, PUIS 2,5 et 97,5 (les deux bouts, pris en premier), PUIS
## les milieux. Les trente premiers percentiles ne bougent donc jamais quand on
## etend, et rb01..rb30 gardent leur sens. Voir ct5, qui construit le vecteur.
stopifnot(ARMB_N %in% c(ARMB_BASE_N, 2L * ARMB_BASE_N))

## ---- Paper figures -------------------------------------------------------------
FAN_T_REF <- 0        # reference point shared by both arms (deg C)
FAN_COLOR <- c(prod = "#d73027", amenity = "#2166ac")
FAN_YLAB  <- c(prod    = "Relative effect on productivity",
               amenity = "Relative effect on amenity value")

## ---- Horizon labels ------------------------------------------------------------
## THE LABEL IS NOT THE COMPUTED YEAR, AND THE TWO ARE DELIBERATELY DIFFERENT.
##
## The quinquennial grid labels each block by its FIRST year, so the last FULL
## block is {2095-2099} and carries the label 2095 inside the computation.
## A block labelled 2100 would hold a single year and be degenerate, so the
## model never evaluates one: every FUT_YEAR, YEAR_MAX and CRH_FUT_YEAR stays
## at 2095L and none of them is touched by what follows.
##
## What the READER is shown is the END of that terminal block, 2100. These
## three constants are the ONLY place that decision is taken. Every producer
## that writes a caption, a table header, an axis title or a file name reads
## them; no producer spells a horizon out itself.
HORIZON_LABEL      <- "2100"        # a level AT the horizon: "by 2100"
HORIZON_SPAN       <- "2015--2100"  # a trajectory: en dash, LaTeX source form
HORIZON_SPAN_PLAIN <- "2015-2100"   # same, for plot labels and file names
HORIZON_BLOCK      <- "2095-2099"   # the terminal block itself, when named

## ---- Climate projection ensemble (pose le 27/09/2026, lot 37) -------------------
## THE FOUR ISIMIP3b GCMs, AS THE FILE NAMES SPELL THEM. Measured on the server
## on 04/07/2026 (MANIFESTE_data_definitif.md, section 1: nine decadal files
## each, all dated <= 2024-10-13), and they are the four the manuscript names.
##
## 🚨 WHY A LIST, WHEN THE PRODUCER DISCOVERS THE MODELS BY LISTING THE DIRECTORY.
## Discovery cannot OMIT a file that is present, and that is all it guarantees.
## A download that brought three models, or one model with a missing decade,
## would be averaged without a word, and the text would still say "four". The
## guard below makes the directory and this list agree, or stops.
CMIP6_MODELS <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0")

## nc_files: paths or basenames. Returns the sorted models, or stops.
cmip6_garde_modeles <- function(nc_files) {
  m <- sub("_r1i1p1f1.*", "", basename(nc_files))
  trouves <- sort(unique(m))
  attendus <- sort(CMIP6_MODELS)
  if (!identical(trouves, attendus))
    stop("[STOP] CMIP6: models present != CMIP6_MODELS.\n",
         "  present  : ", paste(trouves, collapse = ", "), "\n",
         "  expected : ", paste(attendus, collapse = ", "), "\n",
         "  The text names the expected ones: complete the download, or ",
         "change the list AND the manuscript's sentence together.")
  n <- table(m)
  if (length(unique(as.integer(n))) != 1L)
    stop("[STOP] CMIP6: unequal number of files across models (",
         paste(sprintf("%s=%d", names(n), as.integer(n)), collapse = ", "),
         ") -- a partial download would be averaged without a word.")
  ## LES NEUF TRANCHES, PAR MODELE (lot 41). Des comptes egaux ne suffisent pas :
  ## quatre modeles a neuf fichiers dont un porte une tranche fausse passeraient.
  ## Tranches mesurees sur le serveur le 27/09/2026 : 2015_2020, puis 2021_2030
  ## ... 2091_2100, soit 2015-2100 sans trou.
  tr <- sub(".*_(\\d{4}_\\d{4})\\.nc$", "\\1", basename(nc_files))
  for (mo in attendus) {
    t_mo <- sort(tr[m == mo])
    if (!identical(t_mo, sort(CMIP6_SLICES)))
      stop("[STOP] CMIP6: slices of ", mo, " != CMIP6_SLICES.\n",
           "  present  : ", paste(t_mo, collapse = ", "), "\n",
           "  expected : ", paste(CMIP6_SLICES, collapse = ", "))
  }
  trouves
}
CMIP6_SLICES <- c("2015_2020", sprintf("%d_%d", seq(2021L, 2091L, 10L), seq(2030L, 2100L, 10L)))

## ---- Guardrail tolerances ------------------------------------------------------
TOL_BASIS_MATCH <- 1e-12   # MC basis == estimation basis
