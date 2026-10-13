##################################################################################
# lib_contrast.R -- outils communs a la chaine du contraste annuel/journalier.
#
# Reads:  nothing (definitions only)
# Writes: nothing
# Sourced by: ct1_delta_north_grid.R, ct2_between_arm_test.R
#
# POURQUOI CETTE BIBLIOTHEQUE EXISTE
# ----------------------------------
# Hors du paquet, ct2 EXTRAYAIT daily_basis et annual_basis de ct1 par parsing,
# et epinglait le md5 de ct1. Le motif etait le meme que pour lib_drf : ne pas
# creer une seconde source de verite en recopiant. Le remede est le meme :
# une bibliotheque que les deux scripts sourcent, exercee a chaque run.
#
# ⚠ CES TROIS CONSTRUCTEURS NE SONT PAS DES DOUBLONS DE lib_drf.R. Ils
# reconstruisent la base d'ESTIMATION pour en tirer le vecteur de poids w, la
# ou lib_drf construit les TERMES qui entrent dans la regression. Les deux
# doivent coincider, et c'est precisement ce que l'assertion A1 de ct1 mesure a
# chaque specification, a 1e-12. Fusionner les deux effacerait l'assertion.
##################################################################################

suppressMessages({ library(data.table); library(splines); library(tools) })

## L'empreinte se calcule sur les md5 tries PAR NOM, pas dans l'ordre de
## declaration : reordonner PATHS ne doit pas changer le jeton. Le fichier
## intermediaire est ecrit en binaire pour que la valeur ne depende pas de la
## convention de fin de ligne de la plateforme -- ct3 la recalcule en Python et
## les deux doivent tomber sur le meme jeton.
empreinte <- function(md5vec) {
  ## 🚨 CORRIGE LE 2026-09-07 — method = "radix", ET POURQUOI C'EST BLOQUANT.
  ## Version anterieure : order(names(md5vec)), sans methode. order() suit la
  ## COLLATION DE LA LOCALE ; sorted() de Python suit les POINTS DE CODE. Sur les
  ## six cles, mesure sous French_France.utf8 :
  ##   R      : amen_prod.csv | BACI_... | climate_long.rds |
  ##            climate_long_annual.rds | CMIP6_tas.csv | step2_migration.rds
  ##   Python : BACI_... | CMIP6_tas.csv | amen_prod.csv | climate_long.rds |
  ##            climate_long_annual.rds | step2_migration.rds
  ## Les six md5 etaient pourtant IDENTIQUES des deux cotes : seul l'assemblage
  ## differait. Le jeton R (eae054fa...) et le jeton Python (5cd8256c...) ne
  ## pouvaient donc jamais coincider, et le refus de ct3 tirait a CHAQUE run.
  ##
  ## ⛔ Un garde-fou qui refuse toujours n'est pas un garde-fou : c'est celui
  ## qu'on finit par desarmer. Il n'avait jamais laisse passer ct3 -- le smoke
  ## n'exercait que ct1, et ct3 n'avait jamais tourne derriere un ct1 du paquet.
  ##
  ## method = "radix" trie par octets (locale C), ce que fait sorted() en ASCII.
  ## Effet second et voulu : le jeton cesse de dependre de la locale de la
  ## machine, donc deux replicateurs sous des locales differentes obtiennent le
  ## meme. Verifie le 07/09 : R et Python rendent 5cd8256c00... sur les six.
  txt <- paste(sprintf("%s=%s", names(md5vec), md5vec)[order(names(md5vec), method = "radix")],
               collapse = "\n")
  tf <- tempfile(); on.exit(unlink(tf))
  con <- file(tf, "wb"); writeBin(charToRaw(txt), con); close(con)
  unname(tools::md5sum(tf))
}

##################################################################################
# 4. OUTILS : base d'estimation, Delta_north, variance PSD
##################################################################################

# --- reconstruction de la base DAILY, miroir exact de get_effect_curve l.638-651 ---
daily_basis <- function(climate_data, spec, deg, w_lo, w_up) {
  dc <- as.data.table(climate_data)[order(value)]
  tot <- sum(dc$count); cc <- cumsum(dc$count)
  lo <- dc$value[which.max(cc >= w_lo / 100 * tot)]
  hi <- dc$value[which.max(cc >= (1 - w_up / 100) * tot)]
  x_ref <- seq(lo, hi, by = 1)
  if (spec == "polynomial") {
    br <- poly(x_ref, degree = deg, raw = FALSE)
    dcw <- NULL
  } else {
    # 🚨 CORRIGE LE 2026-07-31 — LES NOEUDS SUIVENT L'ESTIMATION, PAS L'AFFICHAGE.
    # Version anterieure (convention d'AFFICHAGE, celle de get_effect_curve l.644-647) :
    #     dcw <- dc[value >= lo & value <= hi]                       # FILTRE le hors-support
    #     tv  <- rep(dcw$value, times = pmax(1L, round(dcw$count)))  # comptages BRUTS
    # Elle est fausse ici, et ce n'est pas une preference : beta_chapeau est estime sur
    # la base de make_climate_terms_wide (l.488-491), qui ECRETE le hors-support sur
    # [lo, hi] et RENORMALISE les comptages a 365 par pays-annee. Ecreter deplace la
    # masse SUR les bornes, filtrer la SUPPRIME : les quantiles different, donc les
    # noeuds, donc la base. Delta_north = w'beta melangeait un vecteur de poids d'une
    # base a des coefficients d'une autre.
    # Mesure de l'ecart (tache2_diag_bases.csv) : 15 des 24 specifications touchees
    # (les 9 polynomiales ne le sont pas, poly() ne dependant que de x_ref) ; ecart sur
    # Delta_north dans [-0,00327 ; +0,00376] ; CONCLUSION INVARIANTE, 0/24 excluent zero
    # sous les deux conventions.
    dcw <- copy(dc)
    dcw[value < lo, value := lo]; dcw[value > hi, value := hi]
    dcw <- dcw[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
    dcw[, count := 365 * count / sum(count), by = .(GID_0, year)]
    # 🚨 CORRIGE LE 2026-09-07 — LE pmax MANQUANT, ET CE QU'IL COUTAIT.
    # Version anterieure :  tv <- rep(dcw$value, times = round(dcw$count))
    # make_climate_terms_wide (lib_drf.R l.~88), qui produit les termes sur
    # lesquels beta est REELLEMENT estime, ecrit  pmax(1L, round(dc$count)).
    # Sans le pmax, toute ligne dont le comptage renormalise arrondit a 0 est
    # SUPPRIMEE du vecteur des quantiles au lieu d'y peser 1. Ces lignes sont
    # 15,2 % du total cote prod et 14,9 % cote amenite : les quantiles bougent,
    # donc les noeuds, donc la base. C'etait le dernier des quatre ecarts entre
    # les deux constructeurs ; le portage du 31/07 en avait ferme trois
    # (ecretage, agregation, renormalisation a 365) et laisse celui-ci.
    #
    # ⚠ CE QUE LE COMPTE DES SPECS QUI PASSENT NE DIT PAS. Mesure du 07/09
    # (logs/.../diag_a1.csv) : 7 des 48 (canal x spec x winsorisation) violaient
    # A1 — prod spline6/0,5 % ; amen spline5 et spline6 aux trois
    # winsorisations — et les 7 tombent a 0,000e+00 avec le pmax.
    # Les 41 autres n'etaient PAS correctes pour autant : le mauvais vecteur
    # servait dans les 48. Les temperatures vivent sur une grille ENTIERE, donc
    # le decalage de quantile ne devient visible que lorsqu'il franchit une
    # borne entiere — auquel cas le noeud saute d'exactement 1,0 (mesure :
    # 7->6 a df=6, 9->8 a df=5). Les 41 passaient par chance de grille, pas par
    # construction. C'est aussi pourquoi l'ecart etait IDENTIQUE entre canaux et
    # winsorisations (1,934e-02 a df=6 ; 2,637e-02 a df=5) : le meme saut d'un
    # degre sur le meme noeud produit le meme ecart de base.
    #
    # ⛔ NE PAS relacher la tolerance de A1 : elle a fait exactement son travail.
    # Le bras ANNUEL, lui, ne porte pas d'assertion (ct1 l.572, a1 <- NA_real_) ;
    # son alignement a ete MESURE le 07/09, 48/48 identiques a 0,000e+00.
    tv  <- rep(dcw$value, times = pmax(1L, round(dcw$count)))
    kp  <- seq(0, 1, length.out = deg + 1)[-c(1, deg + 1)]
    ck  <- quantile(tv, probs = kp, na.rm = TRUE)
    br  <- ns(x_ref, knots = ck, Boundary.knots = c(lo, hi))
  }
  list(basis = br, lo = lo, hi = hi, dcw = dcw)
}

# --- reconstruction de la base ANNUAL, miroir exact de make_annual_terms_trade ---
annual_basis <- function(climate_data, spec, deg, w_lo, w_up, mean_climate_fun) {
  ann <- mean_climate_fun(climate_data, w_lo, w_up)[, .(GID_0, year, annual_temperature)]
  lo <- min(ann$annual_temperature, na.rm = TRUE)
  hi <- max(ann$annual_temperature, na.rm = TRUE)
  x_ref <- seq(floor(lo), ceiling(hi), by = 0.5)
  if (spec == "polynomial") {
    br <- poly(x_ref, degree = deg, raw = FALSE)
  } else {
    kp <- seq(0, 1, length.out = deg + 1)[-c(1, deg + 1)]
    ck <- quantile(ann$annual_temperature, probs = kp, na.rm = TRUE)
    br <- ns(x_ref, knots = ck, Boundary.knots = c(lo, hi))
  }
  list(basis = br, lo = lo, hi = hi)
}

# --- variance de w'beta avec pre-engagement PSD ---
delta_var <- function(w, V) {
  vr <- as.numeric(t(w) %*% V %*% w)
  if (is.finite(vr) && vr >= 0) return(list(var_raw = vr, var_psd = vr, fix = "none", ntrunc = 0L))
   e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  nt <- sum(e$values < 0)
  Vp <- e$vectors %*% diag(pmax(e$values, 0), nrow(V)) %*% t(e$vectors)
  list(var_raw = vr, var_psd = as.numeric(t(w) %*% Vp %*% w), fix = "eigen_trunc", ntrunc = nt)
}

