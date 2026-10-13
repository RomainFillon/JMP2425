##################################################################################
# lib_drf.R -- shared machinery for the dose-response functions
#
# Holds the climate bases, the bilateral differencing, the PPML estimation, the
# Monte Carlo, and the LaTeX and figure helpers used by both channels.
#
# Reads:  nothing directly; every function takes its data as an argument
# Writes: output/aic_weights_*.csv, output/MC_*_fan_paper.pdf
# Needs:  00_parameters.R (DIR_OUT, MC_EXPORT_GRID, FAN_* constants)
# Sourced by 01a_productivity_drf.R, 01b_amenity_drf.R and 01e_diagnostics.R.
#
# The climate basis is defined here once and only here. The Monte Carlo and the
# estimation must build the basis the same way, otherwise beta_hat is applied to
# a basis it was not estimated on. A single definition makes that impossible.
# Nothing channel-specific belongs in this file.
##################################################################################

suppressPackageStartupMessages({
  library(data.table); library(fixest); library(splines); library(MASS)
})


##################################################################################
# 1. CLIMATE BASES
##################################################################################

## Winsorisation bounds, taken from the count-weighted distribution.
winsor_bounds <- function(dc, w_lo, w_hi) {
  dc <- as.data.table(dc)[order(value)]
  cum <- cumsum(dc$count); total <- sum(dc$count)
  list(lower = dc$value[which.max(cum >= w_lo / 100 * total)],
       upper = dc$value[which.max(cum >= (1 - w_hi / 100) * total)])
}

## DAILY basis: integrated over the distribution of days.
## Values are clipped to the bounds and the counts are then renormalised to 365.
## This is the ESTIMATION convention. The alternative -- drop the tails and keep
## the raw counts -- is a DISPLAY convention: it yields a different basis, and
## beta_hat does not live on it.
make_climate_terms_wide <- function(dc, specification = c("polynomial", "spline"),
                                    degree = 3, winsor_lower = 1, winsor_upper = 1) {
  .brut <- dc   # lot 54 : le climat tel que recu, pour l'empreinte de la base ecrite
  specification <- match.arg(specification)
  b <- winsor_bounds(dc, winsor_lower, winsor_upper)
  dc <- copy(as.data.table(dc))
  dc[value < b$lower, value := b$lower][value > b$upper, value := b$upper]
  dc <- dc[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
  dc[, count := 365 * count / sum(count), by = .(GID_0, year)]

  x_ref <- seq(b$lower, b$upper, by = 1)
  if (specification == "polynomial") {
    basis_ref <- poly(x_ref, degree = degree, raw = FALSE); prefix <- "polynomial_"
  } else {
    knots <- quantile(rep(dc$value, times = pmax(1L, round(dc$count))),
                      probs = seq(0, 1, length.out = degree + 1)[-c(1, degree + 1)],
                      na.rm = TRUE)
    basis_ref <- ns(x_ref, knots = knots, Boundary.knots = c(b$lower, b$upper))
    prefix <- "spline_"
    ## Lot 54 : la base d'ESTIMATION est ecrite telle quelle ; le Monte-Carlo et les
    ## courbes la RELISENT (mc_basis_ref) au lieu de la reconstruire.
    enregistre_base("spline", degree, b$lower, b$upper, .brut, knots, basis_ref)
  }
  res <- dc[, as.list(setNames(colSums(predict(basis_ref, value) * (count / 365),
                                       na.rm = TRUE),
                               paste0(prefix, seq_len(degree)))), by = .(GID_0, year)]
  res[, `:=`(lower_threshold = b$lower, upper_threshold = b$upper)][]
}

## Annual mean temperature, winsorised at the same bounds.
mean_climate <- function(dc, winsor_lower = 1, winsor_upper = 1) {
  b <- winsor_bounds(dc, winsor_lower, winsor_upper)
  dc <- copy(as.data.table(dc))
  dc[value < b$lower, value := b$lower][value > b$upper, value := b$upper]
  dc <- dc[, .(count = sum(count, na.rm = TRUE)), by = .(GID_0, year, value)]
  dc[, .(annual_temperature = sum(value * count, na.rm = TRUE) / sum(count, na.rm = TRUE)),
     by = .(GID_0, year)]
}

## ANNUAL basis: evaluated at T_bar, with knots at the quantiles of the ANNUAL
## distribution and boundaries at the min and max of T_bar. This reference is
## deliberately distinct from the daily one; the two are not interchangeable.
make_annual_terms <- function(dc, spec, degree, w_lo, w_hi) {
  annual <- mean_climate(dc, w_lo, w_hi)
  lo <- min(annual$annual_temperature, na.rm = TRUE)
  hi <- max(annual$annual_temperature, na.rm = TRUE)
  x_ref <- seq(floor(lo), ceiling(hi), by = 0.5)
  if (spec == "polynomial") {
    basis_ref <- poly(x_ref, degree = degree, raw = FALSE); prefix <- "annual_poly_"
  } else {
    knots <- quantile(annual$annual_temperature,
                      probs = seq(0, 1, length.out = degree + 1)[-c(1, degree + 1)],
                      na.rm = TRUE)
    basis_ref <- ns(x_ref, knots = knots, Boundary.knots = c(lo, hi))
    prefix <- "annual_spline_"
  }
  annual[, annual_temperature := pmax(pmin(annual_temperature, hi), lo)]
  bv <- as.matrix(predict(basis_ref, annual$annual_temperature))
  colnames(bv) <- paste0(prefix, seq_len(ncol(bv)))
  res <- as.data.table(cbind(annual[, .(GID_0, year)], as.data.frame(bv)))
  res[, `:=`(lower_threshold = lo, upper_threshold = hi)][]
}

##################################################################################
# 2. BILATERAL DIFFERENCES
#
# The direction of the difference is not the same in the two channels, and the
# asymmetry is structural:
#   trade:     exporter = orig, importer = dest, difference orig - dest;
#   migration: origin = r, destination = s, difference (dest - orig) / Omega.
# The arguments `sign_dest_minus_orig` and `scale` carry that choice explicitly
# rather than hiding it in the caller.
##################################################################################

prepare_diffs <- function(data, climate_data, spec, deg, w_lo, w_hi,
                          is_annual = FALSE, fe_cols = "fe_prod_year",
                          sign_dest_minus_orig = FALSE, scale,
                          fe_builder = NULL) {
  ##############################################################################
  # `scale` N'A PAS DE VALEUR PAR DEFAUT, ET C'EST LE CORRECTIF DU 20/09/2026.
  #
  # Elle valait 1. Un appelant qui l'omettait n'obtenait donc pas une erreur,
  # il obtenait la convention « aucune conversion » -- silencieusement, et sous
  # le meme nom de fonction que celui qui convertit. C'est exactement par la que
  # le canal productivite a vecu deux ans avec theta applique deux fois pendant
  # que le canal amenite faisait le bon geste trois lignes plus loin.
  #
  # Desormais un site ajoute plus tard ARRETE au lieu d'heriter. C'est la regle
  # « aucun defaut » deja en vigueur pour run_stage3() et pour THETA_GE.
  ##############################################################################
  if (missing(scale))
    stop("[STOP] prepare_diffs() : l'argument `scale` n'a pas de valeur par ",
         "defaut. Passer SCALE_PROD pour le canal productivite, SCALE_AMEN ",
         "pour le canal amenite (00_parameters.R). Un appel sans echelle ",
         "rendrait des coefficients dans les unites de la regression et non ",
         "dans les unites fondamentales du modele.")
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale == 0)
    stop("[STOP] prepare_diffs() : `scale` doit etre un scalaire fini non nul, ",
         "recu : ", paste(format(scale), collapse = ", "))
  if (is_annual) {
    prefix <- if (spec == "polynomial") "annual_poly_" else "annual_spline_"
    tw <- make_annual_terms(climate_data, spec, deg, w_lo, w_hi)
    dp <- paste0("diff_annual_", spec, "_")
  } else {
    prefix <- paste0(spec, "_")
    tw <- make_climate_terms_wide(climate_data, spec, deg, w_lo, w_hi)
    dp <- paste0("diff_", spec, "_")
  }
  tc <- paste0(prefix, seq_len(deg))
  o <- copy(tw); setnames(o, c("GID_0", tc), c("orig", paste0(tc, "_orig")))
  d <- copy(tw); setnames(d, c("GID_0", tc), c("dest", paste0(tc, "_dest")))
  out <- merge(data, o[, c("orig", "year", paste0(tc, "_orig")), with = FALSE],
               by = c("orig", "year"), all.x = TRUE)
  out <- merge(out,  d[, c("dest", "year", paste0(tc, "_dest")), with = FALSE],
               by = c("dest", "year"), all.x = TRUE)
  for (k in seq_len(deg)) {
    a <- out[[paste0(tc[k], "_orig")]]; b <- out[[paste0(tc[k], "_dest")]]
    out[[paste0(dp, k)]] <- (if (sign_dest_minus_orig) b - a else a - b) / scale
  }
  if (is.null(fe_builder)) {
    if ("fe_prod_year" %in% fe_cols)
      out[, fe_prod_year := paste(product_code_ISICRev3, year, sep = "_")]
    if ("fe_orig_dest" %in% fe_cols) out[, fe_orig_dest := paste(orig, dest, sep = "_")]
    if ("fe_year" %in% fe_cols)      out[, fe_year := as.character(year)]
  } else fe_builder(out)
  out[]
}

##################################################################################
# 2 bis. THE STAGE-3 SAMPLE, WRITTEN RATHER THAN INHERITED (27/09/2026, lot 41)
#
# prepare_diffs() joins the climate terms with all.x = TRUE, so a flow whose
# exporter or importer has no ERA climate carries NA climate differences. Until
# this lot, feglm() dropped those rows by its own default, and the effective
# sample of equation (5) -- 215,092 flows fewer than stage 1, from 20 countries
# with no climate record -- was decided by an estimator's NA handling, not by a
# written filter. stage3_climate_filter() writes the filter; stage3_nobs_check()
# asserts that feglm() then received every row and dropped none. Used by 01a
# only: the other callers of run_stage3() may legitimately lose rows to their
# fixed effects, and the assertion would stop them.
##################################################################################

stage3_climate_filter <- function(data, spec, deg, is_annual = FALSE) {
  dp   <- if (is_annual) paste0("diff_annual_", spec, "_") else paste0("diff_", spec, "_")
  cols <- paste0(dp, seq_len(deg))
  cols <- cols[cols %in% names(data)]
  if (!length(cols))
    stop("[STOP] stage3_climate_filter() : aucune colonne ", dp, "* dans les donnees.")
  keep <- stats::complete.cases(as.data.frame(data)[, cols, drop = FALSE])
  message(sprintf("[stage3] filtre climat : %d lignes gardees, %d sans climat a une extremite",
                  sum(keep), sum(!keep)))
  data[keep]
}

stage3_nobs_check <- function(r, data, label = "") {
  n_in <- nrow(data); n_est <- as.integer(stats::nobs(r$model))
  if (!identical(n_in, n_est))
    stop(sprintf(paste0("[STOP] stage 3 %s : feglm() a estime sur %d lignes, %d lui ",
                        "ont ete fournies apres le filtre climat. Le filtre ecrit et ",
                        "l'echantillon estime divergent."), label, n_est, n_in))
  invisible(n_est)
}

##################################################################################
# 3. PPML ESTIMATION AND JOINT TESTS
##################################################################################

run_stage3 <- function(data, spec, deg, is_annual = FALSE,
                       controls = character(0), fe_str = "fe_prod_year",
                       y = "explained", cluster = c("orig", "dest")) {
  dp <- if (is_annual) paste0("diff_annual_", spec, "_") else paste0("diff_", spec, "_")
  cols <- paste0(dp, seq_len(deg)); cols <- cols[cols %in% names(data)]
  mod <- feglm(as.formula(paste0(y, " ~ ", paste(c(cols, controls), collapse = " + "))),
               data = as.data.frame(data), family = poisson(),
               cluster = cluster, fixef = fe_str)
  in_mod <- cols[cols %in% names(coef(mod))]
  nonlin <- in_mod[grepl("_[2-9]+$", in_mod)]
  list(model = mod, spec = spec, degree = deg,
       wald_joint  = if (length(in_mod)) wald(mod, in_mod) else list(stat = NA, p = NA),
       wald_nonlin = if (length(nonlin)) wald(mod, nonlin) else list(stat = NA, p = NA))
}

##################################################################################
# 4. DOSE-RESPONSE CURVES AND MONTE CARLO
##################################################################################

## Normalisation: each curve is shifted to its OWN maximum, y - max(y).
## Three different normalisations coexist in this work and they are NOT
## reconciled: g(15 deg C) = 0 is the identifying restriction stated in the
## paper, y - max(y) is what the code applies here, and the figures shift the
## MC median to zero at 0 deg C. Any comparison of levels across the three is
## meaningless until one convention is adopted throughout.
get_effect_curve <- function(m, climate_data, is_annual = FALSE) {
  spec <- m$spec; deg <- m$degree
  pre <- if (is_annual) paste0("^diff_annual_", spec, "_[1-9]")
         else           paste0("^diff_", spec, "_[1-9]")
  cols <- grep(pre, names(coef(m$model)), value = TRUE)
  if (!length(cols)) return(NULL)
  if (is_annual) {
    annual <- mean_climate(climate_data, m$wval_lower, m$wval_upper)
    lo <- min(annual$annual_temperature); hi <- max(annual$annual_temperature)
    x_ref <- seq(floor(lo), ceiling(hi), by = 0.5)
    basis_ref <- if (spec == "polynomial") poly(x_ref, degree = deg, raw = FALSE)
      else ns(x_ref, knots = quantile(annual$annual_temperature,
              probs = seq(0, 1, length.out = deg + 1)[-c(1, deg + 1)], na.rm = TRUE),
              Boundary.knots = c(lo, hi))
  } else {
    b <- winsor_bounds(climate_data, m$wval_lower, m$wval_upper)
    lo <- b$lower; hi <- b$upper
    basis_ref <- mc_basis_ref(spec, deg, lo, hi, climate_data)
  }
  x <- seq(lo, hi, by = 0.1); X <- predict(basis_ref, x)
  y <- as.numeric(X %*% coef(m$model)[cols]); y <- y - max(y, na.rm = TRUE)
  se <- sqrt(pmax(rowSums((X %*% vcov(m$model)[cols, cols, drop = FALSE]) * X), 0))
  data.table(x = x, y = y, y_lower = y - 1.96 * se, y_upper = y + 1.96 * se,
             fval_joint = m$wald_joint$stat, pval_joint = m$wald_joint$p,
             fval_nonlin = m$wald_nonlin$stat, pval_nonlin = m$wald_nonlin$p,
             spec = spec, degree = deg)
}

## The Monte Carlo basis must reproduce the estimation basis, never a display
## variant. Falling back to equidistant knots would apply beta_hat to a basis it
## was not estimated on, and nothing in the output would show it.
##################################################################################
# LA BASE SE CALCULE UNE FOIS (lot 54, 28/09/2026).
#
# mc_basis_ref() RECONSTRUISAIT la base de make_climate_terms_wide() avec
# `round(count)` la ou l'estimation ecrit `pmax(1L, round(count))`. Les deux ne
# coincidaient que par accident des poids ; les poids corriges du lot 51 ont
# deplace les noeuds (14/24 contre 15/24) et assert_mc_basis a refuse (RECAP 3.240).
# Desormais l'estimation ECRIT sa base (enregistre_base) et le Monte-Carlo la RELIT
# (lis_base). La garde ne verifie plus une coincidence mais une PROVENANCE :
#   - le fichier relu a l'empreinte enregistree a son ecriture (REGISTRE.tsv) ;
#   - il a ete ecrit pour CE climat (md5 du climat recu) ;
#   - il a ete ecrit dans CE passage : par ce processus (une etape = un Rscript),
#     apres le debut de l'etape.
# ⚠️ Les noeuds sont les quantiles d'un echantillon replique a comptes arrondis,
# avec un plancher pmax(1L, ...) : ce n'est pas un quantile pondere exact. C'est la
# specification d'estimation, inchangee ici (LISTE_POST_SOUMISSION).
##################################################################################
if (!exists(".JMP_T0")) .JMP_T0 <- Sys.time()
.BASES_ECRITES <- new.env(); .EMPREINTES_CLIMAT <- new.env()

dir_bases <- function() {
  racine <- if (exists("DIR_TMP_OUT")) DIR_TMP_OUT else
    file.path(Sys.getenv("JMP_ROOT"), "intermediate_outputs")
  getOption("jmp.dir_bases", file.path(racine, "bases_estimation"))
}

## md5 du climat recu (GID_0, year, value, count, trie). Calcule une fois par
## signature ; la signature n'est qu'un index de cache, la comparaison porte sur le md5.
empreinte_climat <- function(dc) {
  d <- as.data.table(dc)
  sig <- paste(nrow(d), format(sum(as.numeric(d$count)), digits = 17),
               format(sum(as.numeric(d$value)), digits = 17),
               format(sum(as.numeric(d$year)), digits = 17), sep = "|")
  if (!is.null(.EMPREINTES_CLIMAT[[sig]])) return(.EMPREINTES_CLIMAT[[sig]])
  d <- d[, .(GID_0 = as.character(GID_0), year = as.integer(year),
             value = as.numeric(value), count = as.numeric(count))]
  setorder(d, GID_0, year, value)
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  saveRDS(lapply(d, identity), f, compress = FALSE)
  fp <- unname(tools::md5sum(f))
  assign(sig, fp, envir = .EMPREINTES_CLIMAT)
  fp
}

cle_base <- function(spec, deg, lo, hi, fp)
  sprintf("%s%d_lo%s_hi%s_%s", spec, as.integer(deg), format(lo), format(hi), substr(fp, 1, 16))

enregistre_base <- function(spec, deg, lo, hi, dc, knots, basis_ref) {
  fp <- empreinte_climat(dc); cle <- cle_base(spec, deg, lo, hi, fp)
  if (!is.null(.BASES_ECRITES[[cle]])) return(invisible(cle))   # deja ecrite ce passage
  d <- dir_bases(); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  p <- file.path(d, paste0(cle, ".rds"))
  obj <- list(spec = spec, degree = as.integer(deg), lo = lo, hi = hi, knots = knots,
              basis_ref = basis_ref, climat_md5 = fp, ecrit_le = Sys.time(),
              pid = Sys.getpid())
  tmp <- paste0(p, ".tmp"); saveRDS(obj, tmp)
  if (!file.rename(tmp, p)) { unlink(p); if (!file.rename(tmp, p))
    stop("[STOP] enregistre_base : ecriture impossible de ", p) }
  m <- unname(tools::md5sum(p))
  reg <- file.path(d, "REGISTRE.tsv")
  if (!file.exists(reg)) cat("cle\tmd5\tclimat_md5\tecrit_le\tpid\n", file = reg)
  cat(sprintf("%s\t%s\t%s\t%s\t%d\n", cle, m, fp,
              format(obj$ecrit_le, "%Y-%m-%dT%H:%M:%OS3"), obj$pid), file = reg, append = TRUE)
  assign(cle, m, envir = .BASES_ECRITES)
  invisible(cle)
}

## La garde de PROVENANCE. Refuse en nommant l'ecart.
lis_base <- function(spec, deg, lo, hi, climate_data) {
  fp <- empreinte_climat(climate_data); cle <- cle_base(spec, deg, lo, hi, fp)
  d <- dir_bases(); p <- file.path(d, paste0(cle, ".rds"))
  pre <- sprintf("[STOP] provenance de la base MC (%s) : ", cle)
  if (!file.exists(p))
    stop(pre, "aucune base d'estimation ecrite pour ce climat et ces bornes. ",
         "L'estimation n'a pas tourne sur ce climat dans ce passage.")
  m <- unname(tools::md5sum(p))
  reg <- if (file.exists(file.path(d, "REGISTRE.tsv")))
    data.table::fread(file.path(d, "REGISTRE.tsv"), colClasses = "character") else NULL
  enr <- if (is.null(reg)) character(0) else reg$md5[reg$cle == cle]
  if (!length(enr) || !identical(m, enr[length(enr)]))
    stop(pre, sprintf("empreinte du fichier %s, empreinte enregistree a l'ecriture %s. ",
                      substr(m, 1, 8), if (length(enr)) substr(enr[length(enr)], 1, 8) else "<aucune>"),
         "Le fichier relu n'est pas celui qui a ete ecrit.")
  obj <- readRDS(p)
  if (!identical(obj$climat_md5, fp))
    stop(pre, sprintf("base ecrite pour le climat %s, le Monte-Carlo recoit %s.",
                      substr(obj$climat_md5, 1, 8), substr(fp, 1, 8)))
  if (obj$pid != Sys.getpid() || obj$ecrit_le < .JMP_T0)
    stop(pre, sprintf(paste0("base ecrite le %s par le processus %d ; cette etape a ",
                             "commence le %s (processus %d). Base d'un autre passage."),
                      format(obj$ecrit_le, "%Y-%m-%d %H:%M:%S"), obj$pid,
                      format(.JMP_T0, "%Y-%m-%d %H:%M:%S"), Sys.getpid()))
  if (!identical(obj$spec, spec) || obj$degree != as.integer(deg) || obj$lo != lo || obj$hi != hi)
    stop(pre, "l'objet relu ne porte pas la specification demandee.")
  obj
}

mc_basis_ref <- function(spec, deg, lo, hi, climate_data) {
  x_ref <- seq(lo, hi, by = 1)
  if (spec == "polynomial") return(poly(x_ref, degree = deg, raw = FALSE))
  if (is.null(climate_data)) stop("mc_basis_ref : climate_data requis pour une spline.")
  lis_base(spec, deg, lo, hi, climate_data)$basis_ref
}

aic_weights <- function(res) {
  d <- rbindlist(lapply(names(res), function(nm)
    data.table(model_name = nm, aic = AIC(res[[nm]]$model))))
  d[, delta_aic := aic - min(aic, na.rm = TRUE)]
  d[, weight := exp(-delta_aic / 2) / sum(exp(-delta_aic / 2), na.rm = TRUE)][]
}

##################################################################################
# 5. GUARDRAIL: the Monte Carlo basis must reproduce the estimation basis
#
# The check exercises the shipped code on both sides and compares the AGGREGATED
# terms, which is the object the coefficients live on. It is blocking: without
# it, a basis mismatch produces plausible curves and no error.
# Declared scope: the daily arm, splines. Polynomials depend only on x_ref, which
# is shared, and the annual arm has its own reference.
# LOT 54 : le bras JOURNALIER ne compare plus deux constructions -- il n'y en a
# plus qu'une. Il verifie la PROVENANCE de la base relue (lis_base). Le bras
# annuel garde sa comparaison.
##################################################################################
assert_mc_basis <- function(res, climate_data, tol = 1e-12, label = "",
                            is_annual = FALSE) {
  n_ok <- 0L; n_poly <- 0L; worst <- 0; worst_nm <- NA_character_
  for (nm in names(res)) {
    r <- res[[nm]]
    ## Polynomials depend only on x_ref, which both constructions share, so
    ## there is nothing to check.
    if (identical(r$spec, "polynomial")) { n_poly <- n_poly + 1L; next }

    if (is_annual) {
      ## The annual arm is anchored on make_annual_terms: knots at the quantiles
      ## of the ANNUAL distribution, boundaries at the min and max of T_bar. A
      ## different reference from the daily one, but just as checkable.
      ref <- make_annual_terms(climate_data, r$spec, r$degree,
                               r$wval_lower, r$wval_upper)
      annual <- mean_climate(climate_data, r$wval_lower, r$wval_upper)
      lo <- min(annual$annual_temperature); hi <- max(annual$annual_temperature)
      x_ref <- seq(floor(lo), ceiling(hi), by = 0.5)
      kn <- quantile(annual$annual_temperature,
                     probs = seq(0, 1, length.out = r$degree + 1)[-c(1, r$degree + 1)],
                     na.rm = TRUE)
      tb <- pmax(pmin(annual$annual_temperature, hi), lo)
      A <- as.matrix(predict(ns(x_ref, knots = kn, Boundary.knots = c(lo, hi)), tb))
      B <- as.matrix(ref[, paste0("annual_spline_", seq_len(r$degree)), with = FALSE])
      d <- max(abs(A - B))
    } else {
      ## Lot 54 : bras journalier = garde de PROVENANCE. La base n'est plus
      ## reconstruite : lis_base() refuse si l'objet relu n'est pas celui que
      ## l'estimation a ecrit dans ce passage, pour ce climat.
      b <- winsor_bounds(climate_data, r$wval_lower, r$wval_upper)
      invisible(lis_base(r$spec, r$degree, b$lower, b$upper, climate_data))
      d <- 0
    }

    if (!is.finite(d) || d > tol)
      stop(sprintf("BASE MC != BASE D'ESTIMATION (%s, bras %s) : ecart max %.3e.",
                   nm, if (is_annual) "annuel" else "journalier", d))
    if (d > worst) { worst <- d; worst_nm <- nm }
    n_ok <- n_ok + 1L
  }
  ## A check that verified nothing must say so. Silence would read as success.
  if (n_ok == 0L)
    cat("  AUCUNE spline verifiee -- garde-fou INERTE, ce n'est pas un succes.\n")
  else if (is_annual) cat(sprintf("  [%s] %d spline(s) verifiee(s) -- ecart max %.3e (%s)\n",
                   label, n_ok, worst, worst_nm))
  else cat(sprintf(paste0("  [%s] %d spline(s) : base relue = base ecrite par l'estimation ",
                          "de ce passage (empreinte du fichier, climat, processus, date)\n"),
                   label, n_ok))
  cat(sprintf("  %d polynome(s) saute(s) : base insensible a la convention.\n", n_poly))
}

## One draw beta* ~ N(beta_hat, V_hat): beta_hat + z' R, with R'R = V_hat
## (Cholesky) and z = rnorm(p) -- the same p normals, in the same order, that
## MASS::mvrnorm consumed. Until lot 85k (05/10/2026) the root was mvrnorm's
## eigen decomposition, whose eigenvector signs are not continuous in V_hat: a
## 1e-13 run-to-run noise in the two-FE schemes of 01c gave another realisation
## of every draw (rapport_lot85i, rapport_lot85j). The Cholesky factor is
## continuous in V_hat (tools/test_mc_cholesky.R). No jitter: a V_hat that is not
## positive definite stops the run and says by how much, where mvrnorm used to
## fall back silently on beta_hat.
mc_draw_beta <- function(mu, V, what = "") {
  R <- tryCatch(chol(V), error = function(e) {
    ev <- tryCatch(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
                   error = function(e2) NA_real_)
    stop(sprintf(paste0("[STOP] mc_run [%s] : V_hat non definie positive ",
                        "(plus petite valeur propre %.3e) ; chol refuse, aucun jitter."),
                 what, ev), call. = FALSE)
  })
  mu + drop(rnorm(length(mu)) %*% R)
}

## Two layers of uncertainty: the specification is drawn according to the Akaike
## weights, then beta* ~ N(beta_hat, V_hat).
mc_run <- function(res, climate_data, n_sim, x_grid, seed,
                   is_annual = FALSE, label = "") {
  set.seed(seed)
  ## The blocking guardrail is called here rather than by the caller, so that it
  ## runs on every Monte Carlo. Left to the calling scripts it eventually goes
  ## missing from one of them, and a basis mismatch passes unnoticed.
  assert_mc_basis(res, climate_data, label = label, is_annual = is_annual)
  aw <- aic_weights(res)
  out <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    nm <- aw$model_name[sample(nrow(aw), 1, prob = aw$weight)]
    r <- res[[nm]]
    pre <- if (is_annual) paste0("^diff_annual_", r$spec, "_[1-9]")
           else           paste0("^diff_", r$spec, "_[1-9]")
    cols <- grep(pre, names(coef(r$model)), value = TRUE)
    if (!length(cols)) next
    if (is_annual) {
      annual <- mean_climate(climate_data, r$wval_lower, r$wval_upper)
      lo <- min(annual$annual_temperature); hi <- max(annual$annual_temperature)
      x_ref <- seq(floor(lo), ceiling(hi), by = 0.5)
      basis_ref <- if (r$spec == "polynomial") poly(x_ref, degree = r$degree, raw = FALSE)
        else ns(x_ref, knots = quantile(annual$annual_temperature,
                probs = seq(0, 1, length.out = r$degree + 1)[-c(1, r$degree + 1)],
                na.rm = TRUE), Boundary.knots = c(lo, hi))
    } else {
      b <- winsor_bounds(climate_data, r$wval_lower, r$wval_upper)
      lo <- b$lower; hi <- b$upper
      basis_ref <- mc_basis_ref(r$spec, r$degree, lo, hi, climate_data)
    }
    x_in <- x_grid[x_grid >= lo & x_grid <= hi]
    if (length(x_in) < 5) next
    X <- predict(basis_ref, x_in)
    beta <- mc_draw_beta(coef(r$model)[cols], vcov(r$model)[cols, cols, drop = FALSE],
                         sprintf("%s / %s", label, nm))
    y <- as.numeric(X %*% beta); y <- y - max(y, na.rm = TRUE)
    yo <- rep(NA_real_, length(x_grid)); ins <- x_grid >= lo & x_grid <= hi
    yo[ins] <- approx(x_in, y, xout = x_grid[ins], rule = 2)$y
    yo[x_grid < lo] <- y[which.min(abs(x_in - lo))]
    yo[x_grid > hi] <- y[which.min(abs(x_in - hi))]
    out[[i]] <- data.table(x = x_grid, y = yo, sim_id = i)
  }
  ## The Akaike weights decide which specification is drawn, so they are written
  ## out with the draws: publishing a weighted average without its weights would
  ## not be reproducible.
  if (nzchar(label))
    fwrite(aw, file.path(DIR_OUT, sprintf("aic_weights_%s.csv",
                                          gsub("[^A-Za-z0-9]", "_", label))))
  mc <- rbindlist(out)
  cat(sprintf("  MC [%s] : %d/%d tirages valides\n", label, uniqueN(mc$sim_id), n_sim))
  mc
}

## Central curve, weighted by the Akaike weights, on the grid common to all
## specifications. The grid is derived from the supports, never hard-coded.
aic_curve <- function(res, climate_data, is_annual = FALSE) {
  aw <- aic_weights(res)
  cur <- rbindlist(lapply(names(res), function(nm) {
    cc <- get_effect_curve(res[[nm]], climate_data, is_annual)
    if (is.null(cc)) return(NULL)
    cc[, `:=`(model_name = nm, weight = aw$weight[aw$model_name == nm])]
  }))
  g <- seq(max(cur[, min(x), by = model_name]$V1),
           min(cur[, max(x), by = model_name]$V1), by = 0.1)
  itp <- cur[, .(x = g, y = approx(x, y, xout = g, rule = 2)$y, weight = weight[1]),
             by = model_name]
  itp[, .(y_avg = weighted.mean(y, weight, na.rm = TRUE)), by = x][order(x)]
}

##################################################################################
# 6. LATEX EMISSION HELPERS
##################################################################################

## Table note of adaptive width. A hard-coded width such as p{15cm} truncates
## the note as soon as the text layout changes; \dimexpr tracks \linewidth.
tex_note <- function(ncol, lines) {
  c(sprintf("\\multicolumn{%d}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
            ncol), lines, "}")
}

## Provenance header. A generated file must say that it is generated: otherwise
## a hand edit survives until the next run silently overwrites it.
tex_header <- function(producer, source = NULL) {
  c(sprintf("%%%% GENERATED by %s -- DO NOT EDIT BY HAND.", producer),
    if (!is.null(source)) sprintf("%%%% Source: %s", source))
}


##################################################################################
# 7. FAN FIGURES -- shared by both channels
#
# The two channels differ only in the axis label, which FAN_YLAB carries in
# 00_parameters.R.
##################################################################################

mc_summarise <- function(mc) {
  mc[, .(y_mean = mean(y, na.rm = TRUE),
         y_q05 = quantile(y, .05, na.rm = TRUE), y_q25 = quantile(y, .25, na.rm = TRUE),
         y_q50 = quantile(y, .50, na.rm = TRUE), y_q75 = quantile(y, .75, na.rm = TRUE),
         y_q95 = quantile(y, .95, na.rm = TRUE), n_sims = .N), by = x][order(x)]
}

## Export on the wide grid the GE reads. Outside the estimated support the curve
## is extended by a CONSTANT, never extrapolated. The resulting plateau is what
## lets a reader recover the real support from the file alone.
mc_export <- function(s, output_path) {
  x_min <- min(s$x); x_max <- max(s$x)
  r_lo <- s[which.min(abs(x - x_min))]; r_hi <- s[which.min(abs(x - x_max))]
  g <- data.table(temperature_C = MC_EXPORT_GRID)
  for (q in c("q05", "q25", "q50", "q75", "q95", "y_mean")) {
    src <- if (q == "y_mean") "y_mean" else paste0("y_", q)
    g[[q]] <- fifelse(g$temperature_C < x_min, r_lo[[src]],
              fifelse(g$temperature_C > x_max, r_hi[[src]],
                      approx(s$x, s[[src]], xout = g$temperature_C, rule = 2)$y))
  }
  fwrite(g, output_path)
  cat(sprintf("[write] %s | support [%.1f, %.1f]°C\n",
              basename(output_path), x_min, x_max))
}

## The support is derived from the artefact by detecting the plateau, never
## hard-coded, so it cannot drift away from the file it describes.
real_support <- function(q) {
  setorder(q, temperature_C)
  moved <- which(abs(diff(q$q50)) > 1e-12)
  if (!length(moved)) stop("q50 constante sur toute la grille : artefact dégénéré.")
  c(q$temperature_C[min(moved)], q$temperature_C[max(moved) + 1])
}

expo_daily <- function(cl) {
  td <- cl[, .(count = sum(count)), by = value]
  td[, density := count / sum(count)][]
}
expo_annual <- function(cl) {
  td <- cl[, .(count = sum(count)), by = .(value, GID_0)]
  ta <- td[, .(t = sum(value * count) / sum(count)), by = GID_0][
          , .(density = .N), by = .(value = round(t, 0))]
  ta[, density := density / sum(density)][]
}

## The two panels share a window and a scale, so they can be compared by moving
## the eye. The shift applied is SCALAR: the median is zero at T_REF, but the
## band is not pinched there, because each draw was normalised to its own
## maximum. The axis label therefore says "MC median = 0 at 0 deg C" and never
## "normalised", which would suggest the band collapses at that point.
## LA GRILLE EN FOND (2026-09-17, lot 8).
## `grid_curves` : NULL, ou une liste nommee `daily`/`annual` de data.table
## (temperature_C, y, model_name) -- les courbes de POINT de chaque
## specification de la grille, archivees par ct1 dans contribA/spec_curves_*.csv.
## Ce ne sont PAS des medianes Monte Carlo : il n'y a pas de tirage par
## specification. La legende du manuscrit dit « the estimated curve », et c'est
## cette phrase-la qui decrit l'objet trace.
## ⚠ ORDRE. L'archive vient de ct1, qui tourne APRES 01a. La figure du papier
## n'est donc plus emise par 01a mais par 04_contrast/ct6_fig_drf_paper.R :
## une figure dont la legende annonce les vingt-quatre courbes ne peut pas
## etre produite avant que les vingt-quatre existent.
plot_fans <- function(channel, climate_data, grid_curves = NULL) {
  colour    <- FAN_COLOR[[channel]]
  ylab_stub <- FAN_YLAB[[channel]]
  Q <- SUPP <- list()
  for (fr in c("daily", "annual")) {
    Q[[fr]] <- fread(file.path(DIR_OUT,
                 sprintf("MC_damage_%s_%s_quantiles.csv", channel, fr)))
    SUPP[[fr]] <- real_support(Q[[fr]])
  }
  win <- c(max(SUPP$daily[1], SUPP$annual[1]), min(SUPP$daily[2], SUPP$annual[2]))
  stopifnot(FAN_T_REF >= win[1], FAN_T_REF <= win[2])
  for (fr in names(SUPP)) if (win[1] < SUPP[[fr]][1] || win[2] > SUPP[[fr]][2])
    stop("la fenêtre commune déborde un support réel : du padding serait tracé.")

  yr <- range(unlist(lapply(names(Q), function(fr) {
    d <- Q[[fr]][temperature_C >= win[1] & temperature_C <= win[2]]
    sh <- d[temperature_C == FAN_T_REF, q50]; stopifnot(length(sh) == 1L)
    c(d$q05 - sh, d$q95 - sh)
  })), na.rm = TRUE)
  yr <- yr + c(-0.30, 0.06) * diff(yr)

  for (fr in c("daily", "annual")) {
    d <- Q[[fr]][temperature_C >= win[1] & temperature_C <= win[2]]
    shift <- d[temperature_C == FAN_T_REF, q50]
    for (cc in c("q05", "q25", "q50", "q75", "q95", "y_mean"))
      if (cc %in% names(d)) d[[cc]] <- d[[cc]] - shift
    td <- if (fr == "daily") expo_daily(climate_data) else expo_annual(climate_data)
    sec_lab <- if (fr == "daily") "Exposure: share of country-year days"
               else              "Exposure: share of countries"
    x_lab <- if (fr == "daily") "Daily mean temperature" else "Annual mean temperature"
    td <- td[value >= win[1] & value <= win[2]]; stopifnot(nrow(td) > 0)
    y0 <- yr[1]; sf <- diff(yr) * 0.28 / max(td$density)

    ## Les courbes de la grille, recalees COMME la selectionnee : chacune est
    ## translatee pour passer par zero a FAN_T_REF, puis coupee a la fenetre
    ## commune. Une courbe qui ne couvre pas FAN_T_REF ne peut pas etre recalee
    ## et ARRETE : la tracer sur une autre origine la rendrait incomparable.
    gg <- NULL
    if (!is.null(grid_curves)) {
      gg <- data.table::copy(grid_curves[[fr]])
      stopifnot(all(c("temperature_C", "y", "model_name") %in% names(gg)))
      gg[, temperature_C := round(temperature_C, 1)]
      ref <- gg[temperature_C == FAN_T_REF, .(s = y), by = model_name]
      manquantes <- setdiff(unique(gg$model_name), ref$model_name)
      if (length(manquantes))
        stop(sprintf("grille %s : %d courbe(s) sans valeur a %g degC (%s)",
                     fr, length(manquantes), FAN_T_REF,
                     paste(manquantes, collapse = ", ")))
      gg <- merge(gg, ref, by = "model_name")
      gg[, y := y - s]
      gg <- gg[temperature_C >= win[1] & temperature_C <= win[2]]
      stopifnot(nrow(gg) > 0)
    }

    ## geom_rect and not geom_col: geom_col anchors bars at y = 0 and would flip
    ## them whenever the panel floor is negative.
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = td, ggplot2::aes(xmin = value - .45, xmax = value + .45,
                         ymin = y0, ymax = y0 + density * sf),
                         fill = "grey85", alpha = .75) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = .3) +
      ggplot2::geom_vline(xintercept = FAN_T_REF, linetype = "dotted", colour = "grey30", linewidth = .4)

    ## La grille passe AVANT le ruban et la courbe coloree : elle est derriere.
    ## Pas de bande, pas de legende -- l'objet trace est une courbe de point, et
    ## la legende du manuscrit le dit en toutes lettres.
    if (!is.null(gg))
      p <- p + ggplot2::geom_line(data = gg,
                 ggplot2::aes(temperature_C, y, group = model_name),
                 colour = "grey75", linewidth = .5)

    p <- p +
      ggplot2::geom_ribbon(data = d, ggplot2::aes(temperature_C, ymin = q25, ymax = q75),
                           fill = colour, alpha = .18) +
      ggplot2::geom_line(data = d, ggplot2::aes(temperature_C, q25, linetype = "Q25 / Q75"),
                         colour = colour, linewidth = .6) +
      ggplot2::geom_line(data = d, ggplot2::aes(temperature_C, q75, linetype = "Q25 / Q75"),
                         colour = colour, linewidth = .6) +
      ggplot2::geom_line(data = d, ggplot2::aes(temperature_C, q05, linetype = "Q05 / Q95"),
                         colour = colour, linewidth = .55) +
      ggplot2::geom_line(data = d, ggplot2::aes(temperature_C, q95, linetype = "Q05 / Q95"),
                         colour = colour, linewidth = .55) +
      ggplot2::geom_line(data = d, ggplot2::aes(temperature_C, q50, linetype = "MC median"),
                         colour = colour, linewidth = 1.4) +
      ggplot2::scale_linetype_manual(name = NULL,
        values = c("MC median" = "solid", "Q25 / Q75" = "dashed", "Q05 / Q95" = "dotted"),
        breaks = c("MC median", "Q25 / Q75", "Q05 / Q95")) +
      ggplot2::scale_x_continuous(labels = function(x) paste0(x, "°C"),
                                  breaks = seq(-5, 30, 5)) +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ (. - y0) / sf, name = sec_lab)) +
      ## coord_cartesian and not `limits`: `limits` DROPS observations outside
      ## the range, which removes the bars at both ends of the exposure histogram.
      ggplot2::coord_cartesian(xlim = win, ylim = yr, expand = TRUE) +
      ggplot2::labs(x = x_lab,
        y = sprintf("%s\n(MC median = 0 at %g°C)", ylab_stub, FAN_T_REF)) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(legend.position = "bottom",
            legend.key.width = grid::unit(2.2, "lines"),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_blank(),
            axis.title.y.right = ggplot2::element_text(colour = "grey45", size = 12),
            axis.text.y.right  = ggplot2::element_blank())

    f_out <- file.path(DIR_OUT, sprintf("MC_%s_%s_fan_paper.pdf", channel, fr))
    ggplot2::ggsave(f_out, plot = p, width = 11, height = 7, dpi = 300)
    cat(sprintf("[write] %s | fenêtre [%g, %g] | translation %+.5f | grille %s\n",
                basename(f_out), win[1], win[2], -shift,
                if (is.null(gg)) "aucune"
                else sprintf("%d courbes, excursion hors cadre %.3f",
                             length(unique(gg$model_name)),
                             max(0, max(gg$y) - yr[2], yr[1] - min(gg$y)))))
  }
}

##################################################################################
# WORLD BANK LOADER -- wide to long, two temporal aggregations
#
# The two channels do not use the same panel, and the difference is structural:
#   trade:     annual 1995-2019, gaps interpolated, panel not constrained;
#   migration: quinquennial 1985-2015, BALANCED panel required. Stage 3 is
#              identified within pair, so an unbalanced panel would inject
#              composition into the estimate.
# A single function imposing either rule would distort the other channel. The
# `period` argument carries that difference instead of hiding it.
##################################################################################

load_wb <- function(path, value_name, period = c("annual", "quinquennial"),
                    year_min, year_max, periods_expected = NULL,
                    verbose = TRUE) {
  period <- match.arg(period)
  if (!file.exists(path)) stop("[STOP] fichier absent : ", path)

  ## ---- 1. Header: detected, not assumed ------------------------------------------
  ## Two export formats coexist: column names either on the first DATA row (raw
  ## World Bank export) or in a proper file header. The test is whether the first
  ## row holds a recognisable country code, not whether it holds years, since
  ## both formats contain years.
  probe <- data.table::fread(path, header = FALSE, nrows = 1,
                             na.strings = c("", "NA", ".."), colClasses = "character")
  hdr_in_row1 <- any(grepl("Country|iso|GID_0|Code", as.character(unlist(probe)),
                           ignore.case = TRUE))
  if (hdr_in_row1) {
    ## colClasses = "character" is required here. Without it fread types the
    ## header row as numeric, and as.character then returns denormalised floats
    ## ("9.68e-321" instead of "1960").
    raw <- data.table::fread(path, header = FALSE, na.strings = c("", "NA", ".."),
                             colClasses = "character")
    hdr <- trimws(as.character(unlist(raw[1])))
    raw <- raw[-1]
    hdr <- sub("\\.0+$", "", sub("\\s*\\[.*\\]$", "", hdr))
    keep <- which(nzchar(hdr) & !is.na(hdr))
    raw  <- raw[, ..keep]
    data.table::setnames(raw, hdr[keep])
  } else {
    raw <- data.table::fread(path, header = TRUE, na.strings = c("", "NA", ".."))
    data.table::setnames(raw, sub("\\.0+$", "",
                                  sub("\\s*\\[.*\\]$", "", trimws(names(raw)))))
  }
  ## Some exports name the columns "1985 [YR1985]" or "1985.0".

  ## ---- 2. Year and country columns -----------------------------------------------
  yc <- names(raw)[grepl("^[0-9]{4}$", names(raw))]
  if (!length(yc))
    stop("[STOP] ", basename(path), " : aucune colonne d'année sur 4 chiffres.\n",
         "  Colonnes lues : ", paste(head(names(raw), 8), collapse = ", "), " ...")
  cc <- grep("^(Country Code|iso3|iso3c|GID_0|code)$", names(raw),
             ignore.case = TRUE, value = TRUE)[1]
  if (is.na(cc))
    stop("[STOP] ", basename(path), " : colonne de code pays introuvable.\n",
         "  Colonnes lues : ", paste(head(names(raw), 8), collapse = ", "), " ...")

  dt <- raw[, c(cc, yc), with = FALSE]
  data.table::setnames(dt, cc, "GID_0")
  dt[, (yc) := lapply(.SD, as.numeric), .SDcols = yc]
  dt <- data.table::melt(dt, id.vars = "GID_0", variable.name = "year",
                         value.name = value_name)
  dt[, year := as.integer(as.character(year))]

  ## ---- 3. Time window --------------------------------------------------------------
  span <- range(dt$year, na.rm = TRUE)
  dt <- dt[year >= year_min & year <= year_max]
  if (!nrow(dt))
    stop(sprintf("[STOP] %s : le fichier couvre %d-%d, hors de la fenêtre %d-%d.",
                 basename(path), span[1], span[2], year_min, year_max))
  if (verbose)
    cat(sprintf("  [%s] fichier %d-%d | fenêtre %d-%d | %d pays\n",
                value_name, span[1], span[2], year_min, year_max, uniqueN(dt$GID_0)))

  ## ---- 4a. Annual ------------------------------------------------------------------
  if (period == "annual") {
    ## Linear interpolation INSIDE gaps, with no extrapolation at the edges
    ## (rule = 1). Extending a series beyond what it covers would invent a value
    ## where there is no observation.
    data.table::setorder(dt, GID_0, year)
    dt[, (value_name) := {
      v <- get(value_name)
      if (sum(!is.na(v)) < 2L) v
      else stats::approx(seq_along(v), v, xout = seq_along(v), rule = 1)$y
    }, by = GID_0]
    out <- dt[!is.na(get(value_name)) & get(value_name) > 0]
    if (!nrow(out))
      stop(sprintf("[STOP] %s : aucune valeur positive après interpolation.",
                   basename(path)))
    if (verbose) cat(sprintf("  [%s] %d lignes | %d pays | %d-%d\n", value_name,
                             nrow(out), uniqueN(out$GID_0),
                             min(out$year), max(out$year)))
    return(out[])
  }

  ## ---- 4b. Quinquennial -------------------------------------------------------------
  ## A period is labelled by its first year.
  dt[, year := recode_year5(year)]
  dt <- dt[!is.na(get(value_name)) & get(value_name) > 0]
  dt <- dt[, .(v = mean(get(value_name), na.rm = TRUE)), by = .(GID_0, year)]
  data.table::setnames(dt, "v", value_name)
  if (!nrow(dt))
    stop(sprintf("[STOP] %s : aucune valeur positive dans %d-%d.",
                 basename(path), year_min, year_max))

  if (is.null(periods_expected)) {
    if (verbose) cat(sprintf("  [%s] %d lignes | %d pays | périodes %s\n", value_name,
                             nrow(dt), uniqueN(dt$GID_0),
                             paste(sort(unique(dt$year)), collapse = ", ")))
    return(dt[])
  }

  ## ---- 5. Balanced panel -------------------------------------------------------------
  ## Keep only the countries observed in EVERY expected period.
  present <- sort(unique(dt$year))
  manquantes <- setdiff(periods_expected, present)
  ## An expected period missing from the FILE makes a balanced panel impossible
  ## for everyone. That is insufficient coverage, not an over-strict filter. The
  ## two cases look alike and call for different fixes, so they are reported
  ## separately.
  if (length(manquantes))
    stop(sprintf("[STOP] %s : période(s) attendue(s) ABSENTE(S) du fichier : %s\n",
                 basename(path), paste(manquantes, collapse = ", ")),
         sprintf("  Présentes : %s\n", paste(present, collapse = ", ")),
         "  Le panel équilibré est impossible : élargir la fenêtre, ou retirer\n",
         "  ces périodes de PERIODS_EXPECTED.")

  n_par_pays <- dt[year %in% periods_expected, .(n = uniqueN(year)), by = GID_0]
  keep  <- n_par_pays[n == length(periods_expected), GID_0]
  n_out <- uniqueN(dt$GID_0) - length(keep)
  if (!length(keep))
    stop(sprintf("[STOP] %s : AUCUN pays observé sur les %d périodes attendues.\n",
                 basename(path), length(periods_expected)),
         sprintf("  Couverture maximale atteinte : %d période(s), par %d pays.\n",
                 max(n_par_pays$n), sum(n_par_pays$n == max(n_par_pays$n))),
         "  Le filtre de panel équilibré est trop strict pour ce fichier.")
  if (verbose)
    cat(sprintf("  [%s] panel équilibré : %d pays retenus, %d écartés\n",
                value_name, length(keep), n_out))
  dt[GID_0 %in% keep][]
}

##################################################################################
# Quinquennial recoding -- a period is labelled by its first year.
#
# fcase() requires a `default` of length 1, which y is not. The result is
# therefore initialised to y and only the intermediate years are overwritten.
##################################################################################

recode_year5 <- function(y) {
  y <- as.integer(y); out <- y
  map <- list(`1985` = 1986:1989, `1990` = 1991:1994, `1995` = 1996:1999,
              `2000` = 2001:2004, `2005` = 2006:2009, `2010` = 2011:2014,
              `2015` = 2016:2019)
  for (p in names(map)) out[y %in% map[[p]]] <- as.integer(p)
  out
}

##################################################################################
# GLORIA : L'ANNEE SE LIT DANS LE NOM DU FICHIER (lot 34, 25/09/2026)
#
# 01a donnait a la i-eme table du listage l'annee TRADE_YEAR_MIN + i - 1, sur un
# list.files() SANS motif. Juste tant que le repertoire contient exactement
# 1995..2019 et qu'ils se trient dans l'ordre ; un millesime de plus, une table
# d'un autre type (T-Results, V-Results...) ou une annee manquante decalait
# TOUTES les annees, sans un message. Defaut trouve par tools/audit_graphe.py.
#
# Renvoie les chemins TRIES PAR ANNEE et leurs annees, lues dans le nom.
# S'ARRETE si le repertoire porte un fichier qui n'est pas une table Y-Results de
# la release attendue, si une annee est en double, ou si l'ensemble des annees
# n'est pas exactement `years`. Test : tools/test_gloria_annees.R.
##################################################################################
gloria_year_files <- function(dir, years, release = "057") {
  motif <- sprintf("_Y-Results_([0-9]{4})_%s_Markup001\\(full\\)\\.csv$", release)
  tous <- list.files(dir)
  if (!length(tous))
    stop("⛔ aucun fichier GLORIA dans ", dir)
  etrangers <- tous[!grepl(motif, tous)]
  if (length(etrangers))
    stop(sprintf("⛔ %d fichier(s) de %s ne sont pas des tables Y-Results de la release %s : %s",
                 length(etrangers), dir, release, paste(etrangers, collapse = ", ")))
  an <- as.integer(sub(paste0(".*", motif), "\\1", tous))
  if (anyDuplicated(an))
    stop("⛔ annee(s) GLORIA en double : ", paste(sort(unique(an[duplicated(an)])), collapse = ", "))
  manque <- setdiff(years, an); trop <- setdiff(an, years)
  if (length(manque) || length(trop))
    stop(sprintf("⛔ millesimes GLORIA : attendus %d-%d ; manquants [%s] ; en trop [%s]",
                 min(years), max(years), paste(manque, collapse = ", "), paste(trop, collapse = ", ")))
  o <- order(an)
  list(files = file.path(dir, tous[o]), years = an[o])
}

##################################################################################
# garde_ancre_absorption() -- lot 52, 28/09/2026. L'ancre de 01f se LIT.
#
# Elle etait ecrite en dur (0.9970 / 0.9938 / 0.9778) : une photographie de
# l'archive, qui refusait tout changement legitime des donnees. Elle compare
# desormais deux calculs de la meme quantite dans le meme passage :
#   01a l.534-557 : R2 paire x produit sur les objets EN MEMOIRE de 01a
#                   -> fe_absorption.csv, colonne r2_pair_raw ;
#   01f           : le meme R2 recalcule depuis les .rds que 01a a serialises.
# 01f ne lit fe_absorption.csv qu'ici, jamais pour calculer : la garde n'est pas
# circulaire. Elle ne voit pas un defaut de prepare_diffs() (commun aux deux).
# Tolerance : r2_pair_raw est ecrit par fwrite a 15 chiffres significatifs.
##################################################################################
garde_ancre_absorption <- function(r2, chemin, tol = 1e-12) {
  if (!file.exists(chemin))
    stop("[FATAL] ancre de 01f : ", chemin, " absent (produit par 01a au meme passage).")
  a <- data.table::fread(chemin)
  if (!all(c("term", "r2_pair_raw") %in% names(a)))
    stop("[FATAL] ancre de 01f : colonnes term / r2_pair_raw absentes de ", chemin)
  ref <- a$r2_pair_raw[match(seq_along(r2), a$term)]
  if (anyNA(ref))
    stop(sprintf("[FATAL] ancre de 01f : terme(s) %s absent(s) de %s",
                 paste(which(is.na(ref)), collapse = ", "), chemin))
  ecart <- abs(r2 - ref)
  cat(sprintf("  [ancre] ecart max a fe_absorption.csv (r2_pair_raw, 01a meme passage) : %.2e (tolerance %.0e)\n",
              max(ecart), tol))
  bad <- which(ecart > tol)
  if (length(bad))
    stop(sprintf(paste0("[FATAL] la reconstruction de 01f ne retrouve PAS le R2 paire x produit ",
                        "de 01a pour le(s) terme(s) %s : %s. La matrice de design differe de ",
                        "celle de 01a : le chiffre neuf ne serait pas comparable."),
                 paste(bad, collapse = ", "),
                 paste(sprintf("terme %d : 01f %.15g contre fe_absorption.csv %.15g",
                               bad, r2[bad], ref[bad]), collapse = " ; ")))
  invisible(ref)
}
