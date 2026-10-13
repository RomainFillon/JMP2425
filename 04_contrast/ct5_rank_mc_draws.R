##################################################################################
# ct5_rank_mc_draws.R — la bande GE sur des TIRAGES CLASSES, bras par bras.
#
# SORTIES
#   output/mc_draws_ranked_daily.csv    draw_id, rank, dN, is_q05/50/95
#   output/mc_draws_ranked_annual.csv   idem
#   output/mc_draws_selected.csv        les coefficients (perturbations) retenus
#   output/MC_damage_prod_{daily,annual}_quantiles.csv   + colonnes rq05/rq50/rq95
#   output/MC_damage_amenity_{daily,annual}_quantiles_NULL.csv  + memes colonnes
#
# ENTREES
#   intermediate_outputs/contribA/tache1_scores.rds   <- ct2 (pains et viandes)
#   output/contribA/test1_delta_north.csv             <- ct1 (Delta_N chapeau)
#   output/aic_weights_Productivity_{Daily,Annual}.csv <- 01a (la spec tiree)
#   intermediate_outputs/climate_long_annual.rds      <- la base, reconstruite
#
# ---------------------------------------------------------------------------
# CE QUE CE SCRIPT NE FAIT PAS, ET POURQUOI
# ---------------------------------------------------------------------------
# ⛔ AUCUNE COVARIANCE CROISEE ENTRE SPECIFICATIONS. Les deux bras du GE ne
#    tournent pas sur la meme specification -- journalier spline4, annuel
#    spline3, poids d'Akaike 1,000 chacun. ct2 ne construit sa jointe 2k x 2k
#    qu'A L'INTERIEUR d'une cle de specification, et n'archive les scores que
#    pour REF_SPEC. La covariance croisee du couple REEL n'existe nulle part et
#    la fabriquer demanderait de re-ajuster. On tire donc CHAQUE BRAS DANS SA
#    PROPRE SPECIFICATION, et l'encadrement de la marge « extra » est fait de
#    facon CONSERVATRICE (bras traites comme independants), ce qui est le bon
#    sens de l'erreur : leur correlation mesuree est POSITIVE, donc la vraie
#    variance de la difference est PLUS PETITE que celle qu'on suppose.
#
# ⛔ AUCUNE RE-ESTIMATION. Les pains et les viandes sont lus tels qu'archives.
#
# ---------------------------------------------------------------------------
# LA GRILLE SUR LAQUELLE ON COMPARE
# ---------------------------------------------------------------------------
# ⚠️ Les fichiers MC_damage_* sont sur MC_EXPORT_GRID = -100:100, A EXTENSION
#    CONSTANTE hors de MC_X_GRID = [-5, 35] (lib_drf.R::mc_export l.397-406).
#    Toute comparaison de courbes doit se faire sur [-5, 35] : sur [-21, 36]
#    elle oppose une courbe plate a une courbe qui descend, et le classement
#    des specifications s'inverse. Mesure du 16/09/2026 : sur [-21, 36] la
#    mediane du bras journalier « colle » a spline5 ; sur [-5, 35] elle colle a
#    spline4, qui est la specification que les poids d'Akaike designent.
##################################################################################

t0 <- proc.time(); cat("=== ct5_rank_mc_draws.R ===\n")

.root <- Sys.getenv("JMP_ROOT")
if (!nzchar(.root)) {
  .self <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                    error = function(e) NA_character_)
  if (is.na(.self)) {
    .a <- commandArgs(trailingOnly = FALSE)
    .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
    if (length(.f)) .self <- normalizePath(.f[1], winslash = "/")
  }
  if (is.na(.self)) stop("Racine indeterminable. Poser JMP_ROOT.")
  .root <- dirname(dirname(.self))
}
OUT  <- file.path(.root, "output")
INT  <- file.path(.root, "intermediate_outputs")

suppressPackageStartupMessages({
  library(data.table); library(splines)
})
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))
source(file.path(.root, "04_contrast", "lib_contrast.R"))

N_DRAWS   <- 300L
SEED_RANK <- 20260916L      # graine de CE classement, consignee au manifeste
CHANNEL   <- "prod"
Z_SUPPORT <- range(MC_X_GRID)   # [-5, 35] : la seule grille de comparaison

cat(sprintf("  graine = %d | tirages = %d | support de comparaison [%g, %g] degC\n",
            SEED_RANK, N_DRAWS, Z_SUPPORT[1], Z_SUPPORT[2]))

## ---- 1. quelle specification chaque bras tire-t-il ? -------------------------
## Lue dans les poids que 01a archive, jamais supposee.
spec_of_arm <- function(arm) {
  f <- file.path(OUT, sprintf("aic_weights_Productivity_%s.csv",
                              if (arm == "daily") "Daily" else "Annual"))
  if (!file.exists(f)) stop("[STOP] ", basename(f), " absent <- 01a.")
  a <- read.csv(f, stringsAsFactors = FALSE)
  a <- a[order(-a$weight), ]
  if (a$weight[1] < 0.999)
    stop(sprintf("[STOP] bras %s : poids maximal %.4f < 0,999. Le Monte Carlo ",
                 "melange plusieurs specifications ; « la specification du ",
                 "bras » n'a plus de referent.", arm, a$weight[1]))
  a$model_name[1]
}
SPEC <- vapply(c("daily", "annual"), spec_of_arm, "")
cat(sprintf("  bras journalier -> %s\n  bras annuel     -> %s\n",
            SPEC[["daily"]], SPEC[["annual"]]))

## ---- 2. pains, viandes, poids d'exposition ----------------------------------
ST <- readRDS(file.path(INT, "contribA", "tache1_scores.rds"))
T1 <- read.csv(file.path(OUT, "contribA", "test1_delta_north.csv"),
               stringsAsFactors = FALSE)
T1 <- T1[T1$channel == CHANNEL & T1$status == "OK", ]

## La vcov sandwich du bras, telle que ct2 la construit : B M B', corrigee du
## facteur d'echantillon fini. On ne garde que le BLOC CLIMATIQUE : les
## controles ont un poids d'exposition nul, donc ils ne portent ni Delta_N ni
## la courbe, et les inclure ferait tirer du bruit qui ne va nulle part.
arm_pieces <- function(arm) {
  key <- paste(CHANNEL, SPEC[[arm]], sep = "|")
  if (!key %in% names(ST))
    stop("[STOP] ", key, " absent de tache1_scores.rds <- ct2.")
  st <- ST[[key]]
  B <- if (arm == "daily") st$bread_daily else st$bread_annual
  M <- if (arm == "daily") st$meat_daily  else st$meat_annual
  w <- if (arm == "daily") st$w_daily     else st$w_annual
  cc <- if (arm == "daily") st$cc_daily   else st$cc_annual
  V <- st$ssc[[arm]] * (B %*% M %*% t(B))
  idx <- match(cc, names(w))
  if (anyNA(idx)) stop("[STOP] ", arm, " : termes climatiques introuvables dans w.")
  if (any(abs(w[-idx]) > 1e-12))
    stop("[STOP] ", arm, " : un controle porte un poids d'exposition non nul ; ",
         "le bloc climatique ne suffirait plus.")
  dn <- T1$delta_north[T1$arm == arm & T1$spec_label == SPEC[[arm]]]
  se <- T1$se[T1$arm == arm & T1$spec_label == SPEC[[arm]]]
  if (length(dn) != 1L) stop("[STOP] ", arm, " : Delta_N introuvable pour ", SPEC[[arm]])
  list(V = V[idx, idx, drop = FALSE], w = as.numeric(w[idx]), cc = cc,
       dn = dn, se = se, k = length(idx))
}

## ---- 3. tirages, et classement ----------------------------------------------
## La racine de V est le facteur de Cholesky R (R'R = V) : E = Z R a pour
## covariance V. Lot 85l (06/10/2026) : elle remplace la racine spectrale
## tronquee (eigen, pmax(lambda, 0)), discontinue en V -- le signe et l'ordre des
## vecteurs propres peuvent basculer sous une perturbation de 1e-13 et deplacer
## les tirages d'un coup (lot 85k, test_mc_cholesky.R). Memes normales, dans le
## meme ordre. Une V non definie positive ARRETE : aucun jitter, aucune
## troncature (tools/test_ct5_cholesky.R).
ct5_chol_root <- function(V, arm) {
  V <- (V + t(V)) / 2
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R))
    stop(sprintf(paste0("[STOP] ct5 %s : V non definie positive (plus petite ",
                        "valeur propre %.3e) ; chol refuse, aucun jitter."),
                 arm, min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)))
  R
}

draw_arm <- function(arm) {
  P <- arm_pieces(arm)
  R <- ct5_chol_root(P$V, arm)
  ## Graine DERIVEE du nom du bras : les deux bras ne partagent pas une
  ## sequence, et la graine reste reproductible a partir de SEED_RANK seul.
  set.seed(SEED_RANK + if (arm == "daily") 0L else 1L)
  E <- matrix(rnorm(N_DRAWS * P$k), nrow = N_DRAWS) %*% R          # (B, k)
  dN <- as.numeric(P$dn + E %*% P$w)
  ord <- order(dN)
  rk <- integer(N_DRAWS); rk[ord] <- seq_len(N_DRAWS)
  ## Les rangs des percentiles : ceil(p * B), la convention la plus simple et
  ## la seule qui rende un TIRAGE EXISTANT plutot qu'une interpolation.
  sel <- c(q05 = ord[ceiling(0.05 * N_DRAWS)],
           q50 = ord[ceiling(0.50 * N_DRAWS)],
           q95 = ord[ceiling(0.95 * N_DRAWS)])
  cat(sprintf("  [%s] Delta_N chapeau %+.4f (se %.4f) | tirages : q05 %+.4f | q50 %+.4f | q95 %+.4f\n",
              arm, P$dn, P$se, dN[sel["q05"]], dN[sel["q50"]], dN[sel["q95"]]))
  list(P = P, E = E, dN = dN, rank = rk, sel = sel)
}

D <- list(daily = draw_arm("daily"), annual = draw_arm("annual"))

## ---- 4. GARDE : le tirage median reproduit-il la Table 18 ? ------------------
## Le q50 du classement doit tomber au centieme sur Delta_N chapeau : c'est ce
## qui dit que le classement porte bien sur la meme fonctionnelle que la table.
for (arm in c("daily", "annual")) {
  g <- abs(D[[arm]]$dN[D[[arm]]$sel["q50"]] - D[[arm]]$P$dn)
  cat(sprintf("  [garde] %s : |q50 - Table 18| = %.5f\n", arm, g))
  if (g > 0.005)
    stop(sprintf("[HOLD] bras %s : le tirage median s'ecarte de %.4f de la ",
                 "Table 18, au-dela du centieme. Le classement ne porte pas ",
                 "sur la meme fonctionnelle.", arm, g))
}

## ---- 5. ecriture des classements --------------------------------------------
for (arm in c("daily", "annual")) {
  d <- D[[arm]]
  out <- data.frame(
    draw_id = seq_len(N_DRAWS), rank = d$rank, dN = d$dN,
    is_q05 = seq_len(N_DRAWS) == d$sel["q05"],
    is_q50 = seq_len(N_DRAWS) == d$sel["q50"],
    is_q95 = seq_len(N_DRAWS) == d$sel["q95"],
    spec_key = SPEC[[arm]], mc_seed = SEED_RANK,
    stringsAsFactors = FALSE)
  write.csv(out, file.path(OUT, sprintf("mc_draws_ranked_%s.csv", arm)),
            row.names = FALSE)
  cat(sprintf("  wrote mc_draws_ranked_%s.csv\n", arm))
}

sel_rows <- do.call(rbind, lapply(c("daily", "annual"), function(arm) {
  d <- D[[arm]]
  do.call(rbind, lapply(names(d$sel), function(q) {
    b <- d$sel[[q]]
    data.frame(arm = arm, spec_key = SPEC[[arm]], percentile = q, draw_id = b,
               rank = d$rank[b], dN = d$dN[b], mc_seed = SEED_RANK,
               term = d$P$cc, perturbation = as.numeric(d$E[b, ]),
               exposure_weight = d$P$w, stringsAsFactors = FALSE)
  }))
}))
write.csv(sel_rows, file.path(OUT, "mc_draws_selected.csv"), row.names = FALSE)
cat("  wrote mc_draws_selected.csv\n")

##################################################################################
# LE SIDECAR DE MANIFESTE — arm, spec_key, draw_id, draw_rank, mc_seed
#
# ⚠️ CES COLONNES NE VONT PAS DANS ge_runs_manifest.csv, ET C'EST DELIBERE.
# run_one_ge.sh REECRIT la ligne d'un run a chaque execution, avec son propre
# en-tete fige. Des colonnes ajoutees a la main y survivraient jusqu'au prochain
# run de ce fichier, puis dispararaitraient sans bruit -- exactement le genre de
# dispositif que ce projet passe son temps a demonter. Le sidecar est ecrit par
# le producteur des tirages, se joint au manifeste par `file`, et ne peut pas
# etre efface par un run.
#
# La baseline sc1 n'est PAS rejouee : elle est invariante au tirage et a
# l'amenite (meme md5 entre welfare_sc1_q50.csv et son jumeau _NULLamen), donc
# les six runs partagent une seule baseline, nommee ici une fois.
##################################################################################
BASELINE_FILE <- sprintf("welfare_sc1_q50%s.csv", "_NULLamen")
side <- do.call(rbind, lapply(c("daily", "annual"), function(arm) {
  sc <- if (arm == "daily") 2L else 3L      # sc2 = journalier, sc3 = annuel
  d <- D[[arm]]
  do.call(rbind, lapply(names(d$sel), function(q) {
    b <- d$sel[[q]]
    data.frame(
      file = sprintf("welfare_sc%d_r%s_NULLamen.csv", sc, q),
      arm = arm, scenario = sc, spec_key = SPEC[[arm]],
      draw_id = b, draw_rank = d$rank[b], draw_percentile = q,
      dN_partial_equilibrium = d$dN[b], mc_seed = SEED_RANK,
      baseline_file = BASELINE_FILE,
      stringsAsFactors = FALSE)
  }))
}))
write.csv(side, file.path(OUT, "ge_runs_band_ranked.csv"), row.names = FALSE)
cat(sprintf("  wrote ge_runs_band_ranked.csv (%d lignes)\n", nrow(side)))

##################################################################################
# 6. LES COURBES DES TIRAGES RETENUS, DANS LES FICHIERS QUE LE SOLVEUR LIT
#
# main.py choisit une COLONNE du fichier de dose-reponse (`quantile=`). On
# ajoute donc trois colonnes rq05 / rq50 / rq95 plutot que d'ecrire un fichier
# a part : le solveur n'est pas touche, et QUANTILE=rq05 suffit a lancer le run.
#
# LA COURBE D'UN TIRAGE. La normalisation du Monte Carlo est LINEAIRE (chaque
# tirage est ramene a zero en 0 degC), donc la perturbation se propage sans
# terme croise :
#       y_b(T) = y_q50(T) + [X(T) - X(0)] eps_b
# On ANCRE donc la bande sur la colonne q50 PUBLIEE -- celle que les runs
# _NULLamen ont consommee -- au lieu de reconstruire beta chapeau. C'est ce qui
# garantit que le centre de la bande est le run publie et non une re-derivation
# qui en differerait au troisieme chiffre.
#
# ⚠️ La base est CELLE DE L'ESTIMATION (lib_contrast.R), pas une base
#    d'affichage : ct1 assertait deja son accord a 1e-12 sur les termes
#    d'estimation. Melanger les deux est le bug de juillet (noeuds journaliers
#    appliques aux beta annuels) et il ne rend aucune erreur.
##################################################################################
cat("\n--- courbes des tirages retenus ---\n")
clim <- readRDS(file.path(INT, "climate_long_annual.rds")); setDT(clim)

basis_of <- function(arm) {
  sp <- SPEC[[arm]]
  m <- regmatches(sp, regexec("^(polynomial|spline)([0-9]+)_wl([0-9.]+)wu([0-9.]+)$", sp))[[1]]
  if (!length(m)) stop("[STOP] libelle de specification non reconnu : ", sp)
  form <- m[2]; deg <- as.integer(m[3]); wl <- as.numeric(m[4]); wu <- as.numeric(m[5])
  if (arm == "daily") daily_basis(clim, form, deg, wl, wu)
  else                annual_basis(clim, form, deg, wl, wu, mean_climate)
}

add_cols <- function(path, arm, d, B) {
  if (!file.exists(path)) stop("[STOP] ", basename(path), " absent.")
  g <- read.csv(path, stringsAsFactors = FALSE)
  Tgrid <- g$temperature_C
  ## Hors du support d'ESTIMATION la base n'existe pas : on ecrete, exactement
  ## comme ct1 (`pmin(pmax(., lo), hi)`). Hors du support du MONTE CARLO la
  ## colonne publiee est deja plate ; la perturbation doit l'etre aussi, sinon
  ## la bande s'ouvrirait la ou la mediane ne bouge plus.
  Tmc <- pmin(pmax(Tgrid, Z_SUPPORT[1]), Z_SUPPORT[2])
  Tcl <- pmin(pmax(Tmc, B$lo), B$hi)
  X  <- predict(B$basis, Tcl)
  X0 <- predict(B$basis, pmin(pmax(0, B$lo), B$hi))
  Xc <- sweep(X, 2, as.numeric(X0), "-")
  for (q in c("q05", "q50", "q95")) {
    eps <- as.numeric(d$E[d$sel[[q]], ])
    g[[paste0("r", q)]] <- g$q50 + as.numeric(Xc %*% eps)
  }
  write.csv(g, path, row.names = FALSE)
  cat(sprintf("  %-46s + rq05/rq50/rq95\n", basename(path)))
  invisible(g)
}

for (arm in c("daily", "annual")) {
  B <- basis_of(arm)
  cat(sprintf("  [base] %s : support estime [%.2f, %.2f] degC, %d colonnes\n",
              arm, B$lo, B$hi, ncol(as.matrix(B$basis))))
  gp <- add_cols(file.path(OUT, sprintf("MC_damage_prod_%s_quantiles.csv", arm)),
                 arm, D[[arm]], B)
  ## Le bras d'amenite tourne a la borne basse : sa courbe est nulle et ne
  ## depend d'aucun tirage. On ajoute les MEMES noms de colonnes pour que
  ## QUANTILE_PROD et QUANTILE_AMENITY coincident -- sans quoi main.py
  ## fabrique un suffixe de sortie `_ppX_paY` et les fichiers changent de nom.
  fa <- file.path(OUT, sprintf("MC_damage_amenity_%s_quantiles_NULL.csv", arm))
  if (file.exists(fa)) {
    ga <- read.csv(fa, stringsAsFactors = FALSE)
    for (q in c("q05", "q50", "q95")) ga[[paste0("r", q)]] <- ga$q50
    if (max(abs(ga$q50)) > 1e-12)
      stop("[STOP] ", basename(fa), " : la borne basse d'amenite n'est pas nulle ",
           "(max |q50| = ", format(max(abs(ga$q50))), "). Copier q50 serait faux.")
    write.csv(ga, fa, row.names = FALSE)
    cat(sprintf("  %-46s + rq05/rq50/rq95 (copie de q50, nulle)\n", basename(fa)))
  } else {
    stop("[STOP] ", basename(fa), " absent : le bras null ne pourrait pas tourner.")
  }
  ## Controle : en 0 degC la courbe d'un tirage vaut celle de la mediane.
  z <- gp[gp$temperature_C == 0, ]
  for (q in c("rq05", "rq50", "rq95"))
    if (abs(z[[q]] - z$q50) > 1e-10)
      stop(sprintf("[STOP] %s : %s(0) - q50(0) = %.3e, la normalisation ",
                   "en 0 degC n'est pas preservee.", arm, q, z[[q]] - z$q50))
  cat(sprintf("  [garde] %s : les trois courbes coincident avec q50 en 0 degC\n", arm))
}

##################################################################################
# 7. ARM B — LES TIRAGES DU BRAS JOURNALIER
#
# Le total journalier nord et la part gagnante ne dependent QUE du bras
# journalier (sc2 contre la baseline). Des tirages de ce seul bras donnent donc
# une VRAIE distribution de ces deux quantites -- pas un encadrement.
#
# ⛔ LE NOMBRE DE TIRAGES N'EST PAS ECRIT ICI. Il vient de ARMB_N, declare dans
# 01_estimation/00_parameters.R et nulle part ailleurs, et repete cote Python
# dans 00_parameters.py pour l'emetteur. Aucun compte en dur dans ce fichier :
# c'est la lecon du 19/09/2026, ou « trente » etait ecrit dans quatre endroits
# qui ne disaient pas tous la meme chose.
#
# Les colonnes s'appellent rb01, rb02, ... : `run_one_ge.sh` en fait un suffixe
# de nom de fichier, donc elles doivent rester courtes et sans point.
##################################################################################
## ---- la grille de l'arm B : UNE BASE, PUIS DES INTERCALAIRES -----------------
## N VIENT DE 00_parameters.R. Il n'y a plus de variable d'environnement, plus
## de defaut, et donc plus de facon d'obtenir une demi-grille par omission.
##
## CE QUI S'EST PASSE LE 19/09/2026, ET QUE CE BLOC FERME.
## La taille de la grille etait `ARMB_EXTRA <- as.logical(Sys.getenv(
## "ARMB_EXTRA", "FALSE"))`. Le pilote de la nuit a appele ct5 nu tout en
## demandant soixante runs. ct5 a ecrit trente colonnes -- sans erreur, sans
## avertissement, en annoncant trente dans son propre message -- et la passe est
## morte neuf heures plus tard, au deuxieme run, sur `assert quantile in
## damage_consolidated.columns`. Une incompatibilite entre le pilote et ct5 doit
## tomber DANS ct5, en quatre secondes, et pas trente runs plus loin.
##
## LA POIGNEE DE MAIN. Le pilote lit ARMB_N dans 00_parameters.R et le repasse
## dans ARMB_N_PILOTE. Si les deux different, on ARRETE : cela veut dire qu'un
## pilote a ete ecrit contre une autre grille que celle qui est declaree.
## Absente, la variable ne vaut pas accord tacite -- elle vaut « lance a la
## main », et ct5 dit alors bien fort quel N il applique.
.armb_pilote <- Sys.getenv("ARMB_N_PILOTE", "")
if (nzchar(.armb_pilote)) {
  .armb_pilote <- suppressWarnings(as.integer(.armb_pilote))
  if (!is.finite(.armb_pilote))
    stop("[STOP] ARMB_N_PILOTE illisible : '", Sys.getenv("ARMB_N_PILOTE"), "'.")
  if (.armb_pilote != ARMB_N)
    stop("[STOP] le pilote demande ", .armb_pilote, " tirages de bande et ",
         "00_parameters.R en declare ", ARMB_N, ". Un pilote ecrit contre une ",
         "autre grille produirait des runs sans colonne. Aligner les deux ",
         "AVANT de relancer : N se change dans 00_parameters.R, nulle part ",
         "ailleurs.")
  cat(sprintf("  [accord] pilote et 00_parameters.R s'accordent sur N = %d\n",
              ARMB_N))
} else {
  cat(sprintf("  [accord] pas de pilote (ARMB_N_PILOTE absente) : N = %d, lu ",
              ARMB_N))
  cat("dans 00_parameters.R\n")
}

armb_base <- seq(5, 95, length.out = ARMB_BASE_N)
if (ARMB_N == ARMB_BASE_N) {
  armb_p <- armb_base
} else {
  ## Les DEUX bouts d'abord, puis les milieux. Une premiere version mettait
  ## 97,5 en dernier et tronquait : elle perdait exactement la borne haute,
  ## c'est-a-dire la couverture que le garde exige. Ce qu'on sacrifie quand il
  ## faut sacrifier, c'est un milieu, jamais un bout.
  mids <- (head(armb_base, -1) + tail(armb_base, -1)) / 2   # BASE_N - 1 milieux
  armb_p <- c(armb_base, 2.5, 97.5, mids)                   # BASE_N + 2 + (BASE_N-1)
  armb_p <- armb_p[seq_len(ARMB_N)]                         # un milieu en moins
}
if (length(armb_p) != ARMB_N)
  stop("[STOP] la grille construite fait ", length(armb_p), " percentiles pour ",
       "N = ", ARMB_N, ".")
armb_rank <- ceiling(armb_p / 100 * N_DRAWS)
if (anyDuplicated(armb_rank))
  stop("[STOP] deux percentiles tombent sur le meme rang : la grille de ",
       "l'arm B serait plus courte qu'annoncee.")
## Les ARMB_BASE_N premieres colonnes DOIVENT rester ce qu'elles etaient : un
## decalage invaliderait trente runs sans le dire.
if (!isTRUE(all.equal(armb_p[seq_len(ARMB_BASE_N)], armb_base)))
  stop("[STOP] l'extension a deplace les percentiles de base : les runs ",
       "rb01..rb", ARMB_BASE_N, " ne correspondraient plus a leurs colonnes.")

Bd <- basis_of("daily")
gp <- read.csv(file.path(OUT, "MC_damage_prod_daily_quantiles.csv"),
               stringsAsFactors = FALSE)
ga <- read.csv(file.path(OUT, "MC_damage_amenity_daily_quantiles_NULL.csv"),
               stringsAsFactors = FALSE)
Tmc <- pmin(pmax(gp$temperature_C, Z_SUPPORT[1]), Z_SUPPORT[2])
Xb  <- predict(Bd$basis, pmin(pmax(Tmc, Bd$lo), Bd$hi))
X0b <- predict(Bd$basis, pmin(pmax(0, Bd$lo), Bd$hi))
Xcb <- sweep(Xb, 2, as.numeric(X0b), "-")

ordd <- order(D$daily$dN)
armb <- data.frame()
## seq_along(armb_p), PAS seq_len(ARMB_BASE_N) : la grille peut etre etendue, et
## une boucle calee sur la taille de BASE ecrirait trente colonnes en annoncant
## soixante -- le message et l'artefact se contrediraient sans rien casser.
for (i in seq_along(armb_p)) {
  b <- ordd[armb_rank[i]]
  nm <- sprintf("rb%02d", i)
  gp[[nm]] <- gp$q50 + as.numeric(Xcb %*% as.numeric(D$daily$E[b, ]))
  ga[[nm]] <- ga$q50
  armb <- rbind(armb, data.frame(
    file = sprintf("welfare_sc2_%s_NULLamen.csv", nm),
    column = nm, arm = "daily", scenario = 2L, spec_key = SPEC[["daily"]],
    percentile = armb_p[i], draw_id = b, draw_rank = D$daily$rank[b],
    dN_partial_equilibrium = D$daily$dN[b], mc_seed = SEED_RANK,
    baseline_file = BASELINE_FILE, stringsAsFactors = FALSE))
}
write.csv(gp, file.path(OUT, "MC_damage_prod_daily_quantiles.csv"), row.names = FALSE)
write.csv(ga, file.path(OUT, "MC_damage_amenity_daily_quantiles_NULL.csv"), row.names = FALSE)
write.csv(armb, file.path(OUT, "ge_runs_band_armB.csv"), row.names = FALSE)

## ---- RELECTURE SUR DISQUE, ET NON SUR L'OBJET EN MEMOIRE --------------------
## La regle du paquet : tout script qui ecrit dans un fichier partage le RELIT
## avant d'annoncer un succes. Ici elle a un mordant particulier -- c'est
## exactement le controle qui aurait arrete la nuit du 19 a 00 h 42 au lieu de
## 01 h 16, et il porte sur les TROIS fichiers que le solveur va ouvrir.
for (.f in c("MC_damage_prod_daily_quantiles.csv",
             "MC_damage_amenity_daily_quantiles_NULL.csv")) {
  .cols <- names(read.csv(file.path(OUT, .f), nrows = 1, stringsAsFactors = FALSE))
  .rb <- grep("^rb[0-9]+$", .cols, value = TRUE)
  if (length(.rb) != ARMB_N)
    stop("[STOP] ", .f, " porte ", length(.rb), " colonnes rb sur disque, ",
         "ARMB_N = ", ARMB_N, ". Le solveur echouerait sur la premiere ",
         "colonne manquante, et seulement quand il y arriverait.")
  .attendu <- sprintf("rb%02d", seq_len(ARMB_N))
  if (!identical(sort(.rb), sort(.attendu)))
    stop("[STOP] ", .f, " : les noms de colonnes rb ne sont pas rb01..rb",
         sprintf("%02d", ARMB_N), " exactement.")
}
.side <- read.csv(file.path(OUT, "ge_runs_band_armB.csv"), stringsAsFactors = FALSE)
if (nrow(.side) != ARMB_N)
  stop("[STOP] ge_runs_band_armB.csv porte ", nrow(.side), " lignes sur disque ",
       "pour ARMB_N = ", ARMB_N, ".")
cat(sprintf("  [relu] les deux tables de dose-reponse et le sidecar portent %d colonnes rb\n",
            ARMB_N))

cat(sprintf(paste0("  arm B : %d colonnes rb01..rb%02d | percentiles %.1f a ",
                   "%.1f | Delta_N de %+.4f a %+.4f\n"),
            length(armb_p), length(armb_p), min(armb_p), max(armb_p),
            min(armb$dN_partial_equilibrium),
            max(armb$dN_partial_equilibrium)))
cat("  wrote ge_runs_band_armB.csv\n")

cat(sprintf("\nTemps total: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))
cat("=== Done ===\n")
