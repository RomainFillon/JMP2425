##################################################################################
# 02d_emit_latex.R -- emits the LaTeX tables and macros of the manuscript.
#
# Reads:  output/table_eta.csv, output/table_eta_quinquennial.csv,
#         output/pe_spec_ranking.csv, output/pe_spec_invariance_share.csv
# Writes: output/tab_eta.tex
#         (tab_spec_invariance.tex retiree le 11/09/2026 : absorbee par
#          tab_contrast_grid, emise par ct4_emit_tab_contribA_grid.R)
#         then calls the three Python emitters, which write values.tex,
#         values_empirics.tex, the tab_ge_*.tex tables and the two GE figures
# Needs:  every upstream step, since the emitters read their outputs.
#
# The rule of this step: no computed value is typed into a .tex by hand. Every
# number comes from a CSV produced upstream. That is what makes a stale figure
# impossible: if the run changes, the rendered text changes with it.
##################################################################################

suppressMessages(library(data.table))

## No root is set here. 00_parameters.R derives it from its own location and
## sets it for every script; redefining it here would restore a hard-coded path.
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
    stop("Cannot determine the root. Set JMP_ROOT, or run through Rscript.")
  .root <- dirname(dirname(.self))
}
source(file.path(.root, "01_estimation", "00_parameters.R"))
source(file.path(.root, "01_estimation", "lib_drf.R"))


## Table note of adaptive width. A hard-coded width truncates the note as soon
## as the text layout changes; \dimexpr tracks \linewidth.
tex_note <- function(ncol, lines)
  c(sprintf("\\multicolumn{%d}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}{\\footnotesize",
            ncol), lines, "}")

##################################################################################
# 1. tab_eta.tex -- shared share of the temporal variance
##################################################################################

## ---------------------------------------------------------------------
## PARAMETRES DECLARES -> CSV, pour que les emetteurs les LISENT
## ---------------------------------------------------------------------
## Certaines macros portent un parametre du dispositif et non un resultat : le
## nombre de tirages Monte Carlo, la graine, le nombre de pays. Aucun artefact
## ne les portait -- mc_export laisse tomber la colonne n_sims -- et ils
## finissaient donc saisis a la main dans le .tex, ou ils redeviennent faux au
## premier changement, et silencieusement (MC_N_SIM_OH est passe de 300 a 3000
## le 07/09).
## 02d source deja 00_parameters.R : il est l'endroit ou les deposer.
{
  .pv <- list(MC_N_SIM = MC_N_SIM, MC_N_SIM_OH = MC_N_SIM_OH,
              N_COUNTRIES = if (exists("N_COUNTRIES")) N_COUNTRIES else NA,
              MC_X_GRID_LO = min(MC_X_GRID), MC_X_GRID_HI = max(MC_X_GRID))
  .pp <- data.frame(param = names(.pv),
                    value = vapply(.pv, function(x) as.character(x[1]), character(1)),
                    stringsAsFactors = FALSE)
  write.csv(.pp, file.path(DIR_OUT, "config_params.csv"), row.names = FALSE)
  cat(sprintf("  [write] config_params.csv (%d declared parameters)\n", nrow(.pp)))
}

cat("=== tab_eta.tex ===\n")

f_eta <- file.path(DIR_OUT, "table_eta.csv")
if (!file.exists(f_eta)) {
  cat("  [OMITTED] table_eta.csv missing -> tab_eta.tex NOT produced. Run 06.\n")
} else {
  d <- fread(f_eta)

  ## The indexing is POSITIONAL: row 1 is the benchmark on the raw climate
  ## driver and rows 2 to 5 are the damage indices. 01e sets that order and
  ## declares it structural. Changing it on one side only would move the note
  ## silently onto the wrong row.
  if (nrow(d) != 5L)
    stop("[STOP] table_eta.csv: ", nrow(d), " rows instead of 5. ",
         "The order of ETA_SPECS in 06 has changed; the note would be wrong.")
  if (!grepl("^Annual mean", d$metric[1]))
    stop("[STOP] the first row is not the mean-T benchmark: ", d$metric[1])

  pct <- function(m, ch) {
    s <- d[metric == m & channel == ch, eta_pct]
    if (!length(s)) stop("[STOP] not found: ", m, " / ", ch)
    sprintf("%.1f", s[1])
  }

  ## Read when present, omitted when absent, and the omission is announced. A
  ## mechanism that can do nothing must say whether it acted.
  f_q <- file.path(DIR_OUT, "table_eta_quinquennial.csv")
  if (file.exists(f_q)) {
    q <- sprintf("%.3f", fread(f_q)$eta[1])
    quinq <- paste0(", rising to ", q, " on quinquennial aggregates")
    cat(sprintf("  quinquennial eta read: %s\n", q))
  } else {
    quinq <- ""
    cat("  [OMITTED] table_eta_quinquennial.csv missing -> the clause\n",
        "        \"rising to X on quinquennial aggregates\" is NOT written.\n",
        "        The former 0.790 was a literal without an artefact: not carried over.\n")
  }

  ## The four values are named CHANNEL BY CHANNEL, never as a range. Written
  ## "18.2--24.3\%" they are cross-channel bounds, which a reader takes for a
  ## confidence interval; eta has none.
  note <- c(
    "\\textit{Note.} The first row is a benchmark on the raw climate driver:",
    "$\\eta$ computed on annual mean temperature itself, before any damage",
    "transformation (hence no channel). Even this exogenous input is only",
    sprintf("%s\\%% common across countries at annual frequency%s.",
            sprintf("%.1f", d$eta_pct[1]), quinq),
    "Rows 2--5 apply the same decomposition to the estimated damage indices",
    "$D^{j}_{t}(r)$; the gap below the first-row benchmark measures how much",
    "smaller a share of the common climate signal each exposure metric retains",
    "once passed through the estimated dose-response.",
    "Channel by channel, the daily-distribution indices preserve",
    sprintf("%s\\%% of the common component on productivity and",
            pct("Damage index (daily bins)", "Productivity")),
    sprintf("%s\\%% on amenities, against",
            pct("Damage index (daily bins)", "Amenity")),
    sprintf("%s\\%% and %s\\%% for the annual-mean indices.",
            pct("Damage index (annual mean)", "Productivity"),
            pct("Damage index (annual mean)", "Amenity")),
    "These are point diagnostics computed on a single historical sample; they",
    "carry no uncertainty band, and the four values are not interval estimates.")

  rows <- character(0)
  for (i in seq_len(nrow(d))) {
    metric <- sub("Annual mean T", "Annual mean $T$", d$metric[i])
    ch <- if (trimws(as.character(d$channel[i])) %in% c("---", "-", "NA", "")) "---"
          else d$channel[i]
    rows <- c(rows, sprintf("%s & %s & %.3f & %.1f \\\\",
                            metric, ch, d$eta[i], d$eta_pct[i]))
    ## Spacing between the mean-temperature, amenity and productivity blocks.
    if (i %in% c(1L, 3L)) rows <- c(rows, "\\addlinespace")
  }

  writeLines(c(
    "% GENERATED by 02d_emit_latex.R -- DO NOT EDIT BY HAND.",
    "\\begin{table}[H]\\centering",
    paste("\\caption{Share of within-country temporal damage variance absorbed by",
          "year fixed effects ($\\eta^j$), by climate metric and damage channel.",
          "Equal country weights; ERA5; all entries on annual data.}"),
    "\\label{tab:eta}",
    "\\begin{tabular}{llcc}", "\\toprule",
    "Climate metric & Channel & $\\eta^j$ & $\\eta^j$ (\\%)\\\\", "\\midrule",
    rows, "\\bottomrule", tex_note(4, note),
    "\\end{tabular}", "\\end{table}"),
    file.path(DIR_OUT, "tab_eta.tex"))
  cat(sprintf("  [write] tab_eta.tex (%d rows)\n", nrow(d)))
}

##################################################################################
# 2. tab_spec_invariance.tex -- the integrated quantity does not depend on the
#    specification
#
# This table replaces a justification by AIC selection. At this sample size an
# Akaike weight of one is a property of N, not evidence about functional form.
#
# What it deliberately does not do:
#   - it marks no row as selected. An invariance exhibit that names a winner
#     puts back, through the layout, the argument the text removes;
#   - it hard-codes no unanimity ("all twelve", "every specification"). The
#     counts are emitted, so if unanimity breaks the rendered text follows;
#   - it types no percentage where a macro carries the value.
##################################################################################

## ---------------------------------------------------------------------
## tab_spec_invariance.tex N'EST PLUS EMISE  (11/09/2026)
## ---------------------------------------------------------------------
## Absorbee par tab_contrast_grid (ct4_emit_tab_contribA_grid.R) : ses
## douze specifications sont les douze premieres des vingt-quatre, et ses
## deux colonnes -- f_z a 25 et 30 degC, relatives au maximum de chaque
## specification -- y figurent. Les douze valeurs communes sont identiques
## au chiffre pres, et c'est mesure : tools/test_reproduit_warmdrop.R rend
## un ecart de 0,000e+00.
##
## pe_spec_ranking.csv CONTINUE d'etre produit par 01c (R6) : il est la
## REFERENCE de ce controle de reproduction. Le retirer avec la table
## supprimerait le controle en silence.
cat("=== tab_spec_invariance.tex: NOT EMITTED (absorbed by tab_contrast_grid) ===\n")

##################################################################################
# 3. THE PYTHON EMITTERS
#
# They carry the LaTeX macros, the six GE tables and the two GE figures.
#
# A failure here propagates. Printing the failure and moving on would let the
# orchestrator record 02d as a success while no macro and no GE exhibit had been
# written, and the manuscript would then fail to compile with no step blamed.
#
# The emitters are WRITE-LAST: they build everything, run their checks, and only
# then write. A failure therefore leaves no partial output, which is why the
# chain can stop here without cleaning up.

cat("=== Python emitters ===\n")

PY_BIN <- Sys.getenv("PYTHON", unset = "python")

EMETTEURS <- list(
  ## Order matters: the empirical macros depend on no GE run, the other two do.
  ## A failure of the first is therefore an estimation problem, a failure of the
  ## others a GE problem.
  list(f = "emit_values_empirics.py", quoi = "empirical macros (values_empirics.tex)"),
  list(f = "emit_values.py",          quoi = "82 GE macros + 6 tab_ge_*.tex"),
  list(f = "emit_figures.py",         quoi = "ge_welfare_trajectories.pdf, map_welfare_2100_country.pdf")
)

echecs <- character(0)
for (e in EMETTEURS) {
  p <- file.path(.jmp_root, "02_GE", e$f)
  if (!file.exists(p)) {
    ## A missing emitter is not a benign omission: without it the manuscript
    ## loses targets and stops compiling.
    cat(sprintf("  [ABSENT] %-28s %s\n", e$f, e$quoi))
    echecs <- c(echecs, sprintf("%s (missing)", e$f))
    next
  }
  cat(sprintf("  -- %s: %s\n", e$f, e$quoi))
  rc <- system2(PY_BIN, shQuote(p), stdout = TRUE, stderr = TRUE)
  st <- attr(rc, "status")
  ## Le code 3 est un VERDICT rendu (au moins une macro en non-concordance),
  ## pas un plantage. La distinction n'existait pas : emit_values_empirics.py
  ## sortait en 1 pour son verdict comme pour ses neuf gardes dures et pour
  ## toute exception. Elle existe depuis le 2026-09-07.
  ## L'etape s'arrete dans les DEUX cas -- un mismatch laisse la macro NON
  ## emise (elle sort en commentaire dans values_empirics.tex) et le
  ## manuscrit ne compile pas -- mais le motif annonce est desormais exact.
  if (is.null(st) || st == 0L) {
    cat(sprintf("     ok\n"))
    ## LOT 85o : un avertissement non bloquant de l'emetteur (table classee) doit
    ## se lire dans le journal de la course, pas seulement en dry-run.
    av <- grep("[WARNING]", rc, fixed = TRUE, value = TRUE, useBytes = TRUE)
    if (length(av)) cat(paste0("     ", av), sep = "\n")
  } else if (identical(st, 3L)) {
    cat("     VERDICT: at least one macro does not match.\n")
    cat("       The macros concerned are NOT emitted. See the log.\n")
    cat(paste0("       ", utils::tail(rc, 8)), sep = "\n")
    echecs <- c(echecs, sprintf("%s (macro mismatch)", e$f))
  } else {
    cat(sprintf("     FAILURE (code %d):\n", st))
    cat(paste0("       ", utils::tail(rc, 8)), sep = "\n")
    echecs <- c(echecs, sprintf("%s (crash, code %d)", e$f, st))
  }
}

if (length(echecs)) {
  stop(sprintf("[STOP] %d emitter(s) failed: %s\n",
               length(echecs), paste(echecs, collapse = ", ")),
       "  The exhibits they produce are NOT in output/, and the\n",
       "  manuscript will not compile. The failure is passed up to the orchestrator:\n",
       "  do not mask it by declaring the step successful.")
}

## ====================================================================
## TESTS D'ACCEPTATION DE LA COUCHE D'EMISSION
## ====================================================================
## Ils tournent ICI, apres les emetteurs, a chaque passage. Un test qu'on ne
## lance que sur demande ne protege que le jour ou on y pense -- et les deux
## defauts qu'ils gardent sont SILENCIEUX : l'emetteur sort en erreur, ce qui a
## l'air correct, pendant que values.tex est detruit ; ou le PDF compile en
## rendant « ?? » sans que rien ne le signale.
.run_tool <- function(script, label, bloquant) {
  p <- file.path(.root, "tools", script)
  if (!file.exists(p)) {
    cat(sprintf("  [!] %s missing -- check NOT RUN (a gap, not a success)\n",
                script))
    return(invisible(NULL))
  }
  cat(sprintf("\n--- %s ---\n", label))
  out <- system2(PY_BIN, shQuote(p), stdout = TRUE, stderr = TRUE)
  st  <- attr(out, "status"); if (is.null(st)) st <- 0L
  cat(paste0("  ", out), sep = "\n")
  if (st != 0L) {
    if (bloquant)
      stop(sprintf("[STOP] %s failed (code %d). ", label, st),
           "The property it guards is violated SILENTLY: do not declare ",
           "the step successful.")
    cat(sprintf("  [!] %s: %d macro(s) still to be produced.\n", label, st))
  }
  invisible(st)
}

.run_tool("known_bad_emission.py",
          "known-bad: a refusal leaves the artefacts intact", TRUE)
.run_tool("check_macros.py",
          "macros called / defined", FALSE)

cat("\n=== 02d_emit_latex.R done ===\n")
