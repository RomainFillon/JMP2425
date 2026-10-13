#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
emit_values.py -- single source of the general-equilibrium magnitudes.

Reads one arm of the GE run grid and emits:
  - values.tex           the 
ewcommand macros, namespaced \GE and \CALIB
  - output/tab_ge_*.tex  the \input-able tables

Reads:  output/welfare_sc{1,2,3}_*.csv and the matching hat-path files,
        intermediate_outputs/{indices_drop,identifier_xy_final,population0,
        poids_par_cellule}.csv, and the calibration diagnostics when present
Writes: values.tex and up to six tab_ge_*.tex
Needs:  the solver grid must have run.

Idempotent. No GE value is ever written by hand into the manuscript.

Usage:
    python emit_values.py
    RUNS_DIR=<dir> CAP_TAG=<tag> python emit_values.py --dry-run

--------------------------------------------------------------------------
HORIZON CONVENTION -- fixed HERE, once. No table reinterprets it.
--------------------------------------------------------------------------
The model runs over T = 17 five-year periods, from 2015 to 2095.
welfare_sc*.csv has SIXTEEN rows: those are the TRANSITIONS between periods,
not the levels. Row i carries the growth of period i+1, that is
                          year(i) = 2015 + 5*(i+1)
so i=2 is 2030, i=6 is 2050, and i=15 is 2095, the last row.

The cumulative effect at a horizon is cumsum(log(W_A / W_B), axis=0)[i, :].

The LAST row, the terminal 2095-2099 period, is the horizon the paper LABELS
"2100" for readability. The computation never evaluates a one-year 2100 window,
which the quinquennial structure would make degenerate.

--------------------------------------------------------------------------
RUN ARMS
--------------------------------------------------------------------------
The solver leaves the output suffix empty for the headline value of the
column-mass cap, so in a run directory the files WITHOUT a suffix are the
headline arm. The control arms carry their cap in the suffix.
"""

import os
import sys
import numpy as np
import pandas as pd
import rasterio

HERE = os.path.dirname(os.path.abspath(__file__))          # .../02_GE/
PKG  = os.path.dirname(HERE)                               # .../replication/
REPO = PKG

## Horizon labels. Single source, shared with 00_parameters.R; see that file for
## why the label and the computed year differ. The module name starts with a
## digit, so it is imported by name rather than with a plain import statement.
sys.path.insert(0, os.path.join(PKG, '01_estimation'))
from importlib import import_module as _import_module
_PARAMS = _import_module('00_parameters')
HORIZON_LABEL = _PARAMS.HORIZON_LABEL
HORIZON_SPAN  = _PARAMS.HORIZON_SPAN
HORIZON_BLOCK = _PARAMS.HORIZON_BLOCK
INPUTS     = os.path.join(PKG, 'intermediate_outputs')
_RUNS_ENV  = os.environ.get('RUNS_DIR')
## RUNS_DIR arrives in three forms: an absolute path, "output" (the directory
## itself, which is what the orchestrator exports), or the name of a subdirectory
## of output/. Joining it naively to PKG/output would give output/output.
_rd = os.environ.get('RUNS_DIR', '')
if os.path.isabs(_rd):
    RUNS_DIR = _rd
elif _rd in ('', 'output'):
    RUNS_DIR = os.path.join(PKG, 'output')
else:
    RUNS_DIR = os.path.join(PKG, 'output', _rd)
RUNS_DIR = os.path.normpath(RUNS_DIR)
VALUES_TEX = os.path.join(PKG, 'output', 'values.tex')
OUTPUT_DIR = os.path.join(PKG, 'output')
## LOT 85o : BAC A SABLE du controle known-bad. EMIT_SANDBOX_DIR, s'il est pose,
## remplace output/ en lecture ET en ecriture (values.tex, tables, CSV compagnons,
## manifeste). Seul tools/known_bad_emission.py le pose, sur une COPIE d'output/
## dans un repertoire temporaire : le controle ne touche plus l'arbre vivant.
## Refuse s'il designe output/ lui-meme ou un repertoire absent.
_SANDBOX = os.environ.get('EMIT_SANDBOX_DIR', '')
if _SANDBOX:
    _SANDBOX = os.path.normpath(os.path.abspath(_SANDBOX))
    if (os.path.normcase(_SANDBOX) == os.path.normcase(os.path.normpath(OUTPUT_DIR))
            or not os.path.isdir(_SANDBOX)):
        raise SystemExit("[STOP] EMIT_SANDBOX_DIR=%s: neither output/ itself nor a "
                         "missing directory." % _SANDBOX)
    OUTPUT_DIR = _SANDBOX
    VALUES_TEX = os.path.join(_SANDBOX, 'values.tex')
    print("[SANDBOX] EMIT_SANDBOX_DIR=%s: reading and writing outside output/."
          % _SANDBOX)
QUANTILE = os.environ.get('QUANTILE', 'q50')
CAP_TAG = os.environ.get('CAP_TAG', '')

# ----------------------------------------------------------------------
# CE QUE LE BRAS EST, LU AU MANIFESTE -- PAS SUPPOSE
# ----------------------------------------------------------------------
# CAP_TAG est un SUFFIXE DE NOM DE FICHIER. Il ne dit rien de ce que le solveur
# a fait. L'emetteur affichait `CAP_TAG or "C=1.5 (headline)"` : un defaut
# d'f-string qui AFFIRMAIT C=1.5 des que le suffixe etait vide.
#
# Mesure du 11/09/2026 sur output/ge_runs_manifest.csv : les 11 runs portent
# cap_colmass=none. L'etiquette publiee etait donc fausse, et elle l'etait
# partout -- bandeau de values.tex et legende de tab_ge_latitude.
#
# Desormais l'etiquette est LUE. Quand le manifeste manque, elle dit qu'elle ne
# sait pas, ce qui est la seule autre reponse honnete.
def _arm_label_from_manifest():
    p = os.path.join(OUTPUT_DIR, 'ge_runs_manifest.csv')
    if not os.path.exists(p):
        return "cap INCONNU (ge_runs_manifest.csv absent)", None
    try:
        man = pd.read_csv(p)
        caps = sorted({str(v) for v in man['cap_colmass'].fillna('(vide)')})
    except Exception as exc:
        return "cap ILLISIBLE au manifeste (%s)" % exc, None
    if len(caps) != 1:
        return "cap HETEROGENE au manifeste : %s" % ', '.join(caps), None
    cap = caps[0]
    return ("C=%s" % cap if cap not in ('none', 'None', '(vide)')
            else "SANS cap de column-mass (cap_colmass=none)"), cap


ARM_LABEL, ARM_CAP = _arm_label_from_manifest()


# ----------------------------------------------------------------------
# L'ASSERTION QUI MANQUAIT : UNE TABLE, UN BRAS
# ----------------------------------------------------------------------
# Toute table dont les lignes viennent de plusieurs runs doit lire des runs qui
# partagent le cap ET la specification d'amenite. Sans cela, deux colonnes de la
# meme page decrivent deux modeles, et rien ne le dit.
#
# Elle porte sur les fichiers REELLEMENT lus, passes par l'appelant. Une liste
# gelee resterait verte apres un realignement de chemin -- c'est le defaut que ce
# paquet a deja publie une fois.
#
# Polarite : ARRET. Ce n'est pas une propriete scientifique, c'est la coherence
# du montage ; une divergence est une panne.
def assert_one_arm(files, label):
    """files : noms de base des welfare lus par une meme table."""
    p = os.path.join(OUTPUT_DIR, 'ge_runs_manifest.csv')
    if not os.path.exists(p):
        print("[bras] %s: ge_runs_manifest.csv missing -- check NOT RUN, "
              "which is a gap, not a success." % label)
        return None
    man = pd.read_csv(p).set_index('file')
    noms = [os.path.basename(f) for f in files]
    absents = [n for n in noms if n not in man.index]
    if absents:
        sys.exit("STOP [%s]: missing from the manifest: %s. A run that is not "
                 "listed there has no known arm." % (label, ', '.join(absents)))
    sub = man.loc[noms]
    caps = sorted({str(v) for v in sub['cap_colmass'].fillna('(vide)')})
    amen = sorted({str(v) for v in sub['amenity_spec'].fillna('(vide=null)')})
    runs = sorted({str(v) for v in sub['run_id'].fillna('(vide)')})
    ##########################################################################
    # LE MILLESIME, AU MEME TITRE QUE LE BRAS. Ajoute le 20/09/2026.
    #
    # Cette fonction verifiait que les runs d'une table partagent un BRAS
    # (cap_colmass, amenity_spec). Elle ne verifiait pas qu'ils partagent une
    # CONVENTION. Le 20/09 le depot portait les deux a la fois : sc2/sc3 du
    # 19/09 refaits apres la division par theta, et sc2_q05/q95, sc3_q05/q95,
    # sc4_dailycurve, les six rq*, tous d'avant. La colonne q50 d'une bande
    # aurait ete 6,5 fois plus petite que ses colonnes q05 et q95, sur la meme
    # ligne de la meme table, sans qu'aucun controle puisse le voir : les six
    # runs partagent bien un bras.
    #
    # `inputs_md5` hache desormais SCALE_PROD (run_one_ge.sh), donc deux
    # conventions ne peuvent plus partager une empreinte. C'est la meme cle que
    # 02e_north_contrast_band.py et decompose_metric_channels.py emploient deja
    # pour dire « un seul millesime » ; elle est ici etendue a TOUTE table.
    #
    # Polarite : ARRET. Une table faite de deux conventions n'est pas une table
    # imprecise, c'est une table qui compare deux modeles.
    ##########################################################################
    fps = sorted({str(v) for v in sub['inputs_md5'].fillna('(vide)')})
    if len(fps) != 1:
        det = '\n'.join('     %-46s %s  %s' % (n, str(man.loc[n, 'inputs_md5'])[:8],
                                               str(man.loc[n, 'timestamp']))
                        for n in noms)
        sys.exit("STOP [%s]: the runs of this table carry %d input "
                 "fingerprints, hence %d conventions.\n%s\n"
                 "   inputs_md5 hashes SCALE_PROD: a different "
                 "fingerprint is a different scale, not a different "
                 "timestamp. Re-run the runs that carry "
                 "the old one." % (label, len(fps), len(fps), det))
    if len(caps) != 1 or len(amen) != 1:
        sys.exit("STOP [%s]: the runs of this table do not share one arm.\n"
                 "   cap_colmass  : %s\n   amenity_spec : %s\n"
                 "   Two columns of the same page would describe two models."
                 % (label, ', '.join(caps), ', '.join(amen)))
    print("[bras] %-22s %d run(s) | cap=%s | amenity=%s | md5=%s | run_id=%s"
          % (label, len(noms), caps[0], amen[0], fps[0][:8],
             runs[0] if len(runs) == 1 else 'HETEROGENEOUS:' + ','.join(runs)))
    return {'cap': caps[0], 'amenity': amen[0], 'run_ids': runs,
            'inputs_md5': fps[0]}

# ----------------------------------------------------------------------
# PROVENANCE OF THE RUN GRID
# ----------------------------------------------------------------------
# All three scenarios are read from ONE directory, the one the orchestrator
# writes. An earlier layout composed the headline arm from three separate
# directories, because three vintages of welfare coexisted on disk. A composed
# path is dangerous in a way a single one is not: someone realigns it in good
# faith onto one directory and silently brings back pre-rerun values. The guard
# below exists for that case.
_HYBRIDE = False

# ----------------------------------------------------------------------
# GUARD ON AN UNDECLARED RUN DIRECTORY
# ----------------------------------------------------------------------
# The dangerous case is emitting without saying where the runs come from, while
# more than one vintage sits on disk. The emitter would then regenerate every
# macro from an old grid and exit green, which is how months of stale figures
# reach a manuscript with no signal at all.
#
# The guard is CONDITIONAL: it fires only on that case. Once RUNS_DIR is set the
# behaviour is unchanged, and the mutation test checks both directions.
# Reproducing an earlier state stays possible by declaring it explicitly. What is
# forbidden is the silence, not the composition.
if _HYBRIDE and os.environ.get('RUNS_HYBRIDE') != '1':
    sys.exit(
        "\n[STOP] emit_values.py: HYBRID layout active and RUNS_DIR not set.\n"
        "\n"
        "  Without RUNS_DIR, the 82 \\GE* macros would be regenerated from:\n"
        "      sc1  runs_cap_sweep/welfare_sc1_q50.csv   (05/07)\n"
        "      sc2  runs_2208/welfare_sc2_q50.csv        (05/08)\n"
        "      sc3  runs_cap_sweep/welfare_sc3_q50.csv   (19/07)\n"
        "  ...and the emitter would exit green. After a re-run, the\n"
        "  OLD values would land in values.tex.\n"
        "\n"
        "  WHAT TO SET -- choose ONE of the two:\n"
        "\n"
        "  1) AFTER A RE-RUN, point at the fresh pull (normal case):\n"
        "         RUNS_DIR=runs_2808 python ge_analysis/emit_values.py\n"
        "     The directory must contain the SIX welfare_sc{1,2,3}_q50.csv\n"
        "     and welfare_sc{1,2,3}_q50_cap1p2.csv. Check: the six md5s of\n"
        "     data_processed/BASELINE_values_GE.tsv must ALL have changed.\n"
        "\n"
        "  2) REPRODUCE THE PUBLISHED STATE (July-August), declaring it:\n"
        "         RUNS_HYBRIDE=1 python ge_analysis/emit_values.py\n"
        "\n"
        "  ⛔ RUNS_DIR=runs_2208 fails: runs_2208 is a DELTA (sc2 only),\n"
        "     no sc1, no cap arm, no null arm. This is not an anomaly.\n"
        "\n"
        "  Full procedure: replication_package/PROCEDURE_post_rejeu.md\n")

## One directory for all three scenarios: the one the orchestrator writes.
SCENARIO_RUNS = {sc: RUNS_DIR for sc in (1, 2, 3)}
# A mechanism that can do nothing must say whether it acted.
print("[runs] arm %s: sc1=%s sc2=%s sc3=%s"
      % ("COMPOSED BY SCENARIO (headline)" if _HYBRIDE else "homogeneous (override)",
         *(os.path.basename(SCENARIO_RUNS[sc]) for sc in (1, 2, 3))))

# ----------------------------------------------------------------------
# THE RUN DIRECTORY OF EACH MACRO FAMILY, AND A PER-RUN GUARDRAIL
# ----------------------------------------------------------------------
# What catches a stale figure is not a column name: it is the RUN. Each family of
# macros therefore declares below the files it descends from, with their
# fingerprints. check_run_fingerprints() STOPS the emitter when one of those
# files changes without a deliberate update.
#
# Some families can legitimately sit on a different vintage from the headline,
# because the control arms are not re-run as often. The table SAYS so rather than
# hiding it, and the banner of values.tex reports it to whoever reads the emitted
# file.
_R = lambda d, f: os.path.join(HERE, d, f)

RUN_FAMILIES = [
    # (family, macros or exhibits, files read, vintage, note)
    ('headline',
     'GEnorth* GEagg* GEcty* GEtier* GEpopLoss* GEpopGain* GEnorthLat*',
     # DERIVED from SCENARIO_RUNS, never copied. If the emitter is realigned onto
     # another directory these paths move with it and the manifest below stops
     # recognising them. A frozen list would instead stay green, attesting to
     # files that are no longer read.
     [os.path.join(SCENARIO_RUNS[sc], f'welfare_sc{sc}_{QUANTILE}{CAP_TAG}.csv')
      for sc in (1, 2, 3)],
     'headline',
     'configuration composed by scenario (see the block above)'),

    ('endogenous channel',
     'GEchanIncomeTrop GEchanIncomeMid GEchanIncomeSub GEchanIncomeArc',
     [os.path.join(RUNS_DIR, f'amenity_hat_path_sc2_{QUANTILE}.csv'),
      os.path.join(RUNS_DIR, f'amenity_hat_path_sc3_{QUANTILE}.csv')],
     'read from the manifest',
     'MIX: chan_income = extra(sc2 re-run) - cum_amen(sc2 before re-run). '
     'runs_2208/ CONTAINS the re-run hat-paths -- this family is the only one '
     'of the three that can be repaired without a GE re-run. Not repaired here: the '
     'repair moves published values and awaits a decision.'),


    ('null arm g=0',
     'GEsensNorthLo GEsensPreserved* + columns g^LB=0 of tab_ge_north_sens '
     'and row " Lower bound " of tab_ge_sensitivity',
     [os.path.join(RUNS_DIR, 'welfare_sc1_q50_NULLamen.csv'),
      os.path.join(RUNS_DIR, 'welfare_sc2_q50_NULLamen.csv'),
      os.path.join(RUNS_DIR, 'welfare_sc3_q50_NULLamen.csv')],
     'read from the manifest',
     'MIX: in tab_ge_north_sens the q50 columns are post re-run and the '
     'g^LB=0 columns are from before the re-run, side by side. runs_2208/ contains '
     'NO _NULLamen: cannot be repaired without a GE re-run.'),

    ('mask decomposition',
     'GEmask*',
     [],   # read through output/exposure_mask_decomposition.csv (02i)
     'read from the manifest',
     'four masked runs, own manifest: intermediate_outputs/'
     'exposure_masks/runs_manifest.csv, checked by 02i (base fingerprint = '
     'headline runs). The emitter requires the full column to reproduce the body '
     'to 1e-9.'),

    ('sc4 channels',
     'GEnorthChanEval GEnorthChanShape',
     [],   # read by the channel decomposition, which pins its own md5
     'read from the manifest',
     'MIX: decompose_metric_channels.py pins SC2=runs_cap_sweep by md5 '
     'and REFUSES to run if the input changes -- locked on the pre-re-run state. '
     'Its cross-check compares to a hard-coded PUBLISHED_DAILY_EXTRA_NORTH=-0.105 '
     '(l. 98) and prints OK while the emitted value is -0.131.'),
]

# Fingerprints, measured on the files themselves.
RUN_FINGERPRINTS = {
    'runs_cap_sweep/welfare_sc1_q50.csv':          '3b2b238fdf256b11f50392f445ba7c92',
    'runs_2208/welfare_sc2_q50.csv':               '4a60ee2fb06c796492b1dbc68c50da14',
    'runs_cap_sweep/welfare_sc3_q50.csv':          '3ac954e5b091d8389bc9c4515ee20563',
    'runs_cap_sweep/welfare_sc1_q50_cap1p2.csv':   '3cf2041efce4d756ffa8c4f2939a3655',
    'runs_cap_sweep/welfare_sc2_q50_cap1p2.csv':   '85abf1a6d2fd8afa52348b1b82bc5a5a',
    'runs_cap_sweep/welfare_sc3_q50_cap1p2.csv':   '6bc44082695d03db65d1223503b9d932',
    'runs_cap_sweep/welfare_sc1_q50_NULLamen.csv': '3b2b238fdf256b11f50392f445ba7c92',
    'runs_cap_sweep/welfare_sc2_q50_NULLamen.csv': 'f36ae800ffef5efa2ace80e0d557764c',
    'runs_cap_sweep/welfare_sc3_q50_NULLamen.csv': '0e05d08918c3958b6098213a3fc4fcfa',
    'runs_cap_sweep/amenity_hat_path_sc2_q50.csv': 'e21bab676aa48054abfd0432e44c7f37',
    'runs_cap_sweep/amenity_hat_path_sc3_q50.csv': 'bbade436021c00047a32ed084f911019',
}


def check_run_fingerprints():
    """Garde-fou PAR RUN. Arrete si un fichier de run lu par une famille de
    macros n'a plus l'empreinte gelee, OU si l'emetteur lit un fichier qui
    n'est pas au manifeste.

    Les deux moities comptent, et la seconde est celle qui manquait : le defaut
    publie pendant quatre mois n'etait pas un fichier modifie, c'etait un
    fichier REMPLACE PAR UN AUTRE. Un controle qui verifie une liste gelee ne
    voit que la premiere moitie.

    (!) Ce garde-fou n'assert AUCUNE conclusion : il compare des empreintes de
        fichiers d'entree. Contrat SS0, regle jumelle. Les CSV de welfare n'ont
        pas d'horodatage interne, le md5 y est un controle valide (a la
        difference des PDF et des tables xtable)."""
    ## The cap and null arms are not produced by the current run grid. The
    ## families that read them are skipped when their files are missing, and the
    ## omission is ANNOUNCED: a published arm that disappears in silence is worse
    ## than one that is absent.
    if not _HYBRIDE:
        # Say when nothing is done, and why.
        print("[fingerprint] NOT RUN: homogeneous arm (RUNS_DIR or CAP_TAG "
              "override). The frozen fingerprints describe the headline arm.")
        return
    import hashlib
    lus = []
    for fam, _macros, files, _mill, _note in RUN_FAMILIES:
        for f in files:
            rel = os.path.relpath(f, HERE).replace(os.sep, '/')
            if (fam, rel) not in lus:
                lus.append((fam, rel))

    hors = [(fam, rel) for fam, rel in lus if rel not in RUN_FINGERPRINTS]
    hors_rel = {rel for _fam, rel in hors}
    bad, missing, seen = [], [], 0
    for fam, rel in lus:
        if rel in hors_rel:
            continue
        p = os.path.join(HERE, rel)
        if not os.path.exists(p):
            missing.append(rel)
            continue
        h = hashlib.md5()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        seen += 1
        if h.hexdigest() != RUN_FINGERPRINTS[rel]:
            bad.append((rel, RUN_FINGERPRINTS[rel], h.hexdigest()))

    print("[fingerprint] %d files read, %d checked by md5%s%s"
          % (len(lus), seen,
             (", %d OUTSIDE THE MANIFEST" % len(hors)) if hors else "",
             (", %d missing" % len(missing)) if missing else ""))
    for fam, rel in hors:
        print("   OUTSIDE MANIFEST [%s]  %s" % (fam, rel))
    for rel in missing:
        print("   MISSING         %s" % rel)
    for rel, want, got in bad:
        print("   DRIFT           %s\n       frozen=%s\n       read  =%s"
              % (rel, want, got))

    if hors:
        sys.exit(
            "STOP: the emitter reads a run file that is not in the "
            "manifest. That is the signature of a path REALIGNMENT -- the "
            "move by which pre-re-run values would come back without "
            "anything signalling it. If the realignment is intended: update "
            "RUN_FINGERPRINTS AND the RUN_FAMILIES line that carries it, "
            "noting the run, then re-run BOTH mutation tests.")
    if bad:
        sys.exit(
            "STOP: a run file changed without a deliberate update of "
            "RUN_FINGERPRINTS. That is exactly the defect this emitter "
            "published for four months -- a replaced sc2, macros emitted "
            "without anything signalling it. If the change is intended: update "
            "the fingerprint AND the RUN_FAMILIES line that carries it, "
            "noting the run, then re-run BOTH mutation tests.")
    if missing:
        print("   (the missing files belong to an incomplete arm, not to a drift)")


check_run_fingerprints()
DRY_RUN  = '--dry-run' in sys.argv

# Horizons visible in the tables. The key is the DISPLAYED label, so the
# terminal key is HORIZON_LABEL and not the year the model indexes on: the last
# row is the transition 2095-2099 and the model never evaluates a 2100 window,
# but the reader is shown the END of that block. The decision lives in
# 00_parameters.R / .py; nothing here restates it.
# TERMINAL_ROW is a ROW INDEX and is unaffected: no computation moves.
# Les NOMS DE MACRO sont inchanges (HORIZON_WORD garde 'Hundred') : renommer les
# macros deplacerait le .tex sans rien corriger.
HORIZON_ROW = {'2030': 2, '2050': 6, HORIZON_LABEL: 15}
TERMINAL_ROW = 15

# LaTeX forbids digits in a command name, so \GEaggDaily2030 is illegal. The
# horizons are therefore spelled out in the macro names.
HORIZON_WORD = {'2030': 'Thirty', '2050': 'Fifty', HORIZON_LABEL: 'Hundred'}

NORTH_LAT = 45.0

# Country categories are a FRAMING DECISION, not an output of this generator.
# The generator checks that they remain true; see the assertions below.
REVERSAL     = [('RUS', 'Russia'), ('FIN', 'Finland'),
                ('NOR', 'Norway'), ('SWE', 'Sweden')]
AMPLIFICATION = [('CAN', 'Canada'), ('DEU', 'Germany')]
FRANCE        = [('FRA', 'France')]
MACRO_NAME = {'RUS': 'Rus', 'FIN': 'Fin', 'NOR': 'Nor', 'SWE': 'Swe',
              'CAN': 'Can', 'DEU': 'Deu', 'FRA': 'Fra'}
# L'ordre du tableau des economies froides. DEFINI ICI, et non plus pres de son
# premier usage d'affichage : aggregates_for() en a besoin pour compter les
# economies sur CHAQUE bras, et il tourne bien avant.
COLD_ORDER_ISO = REVERSAL + AMPLIFICATION + FRANCE

# ----------------------------------------------------------------------
# CALIB_FROM_NOTE WAS REMOVED ON 25/09/2026, with its five macros:
# CALIBconcNoCap, CALIBconcIncome, CALIBallocCorr, CALIBclearWstar,
# CALIBinnerResid. They were result literals transcribed from a working note and
# never recomputed in this tree. They described the lambda0 that 02a produced
# before its correction of 23/09 and, for three of them, a C=1.5 capped arm that
# no run has. Zero calls in the manuscript. Their values and what they claimed
# to describe are recorded in RECAP.md, entry 3.215.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Column-mass diagnostic of lambda0 and mu0. Its producer pins the md5 of its
# inputs and refuses to run if lambda0 does not reproduce the reference. Its
# absence is traced rather than filled in; see the emission block below.
# ----------------------------------------------------------------------
_CAP_DIAG_CSV = os.path.join(REPO, 'output', 'cap_colmass_diagnostic.csv')
CAP_COLMASS_DIAG = (pd.read_csv(_CAP_DIAG_CSV)
                    if os.path.exists(_CAP_DIAG_CSV) else None)
print(f"[load] cap_colmass_diagnostic.csv : "
      + ("OK" if CAP_COLMASS_DIAG is not None else "ABSENT (4 macros non emises)"))


# ======================================================================
# 1. Loading
# ======================================================================
def load_welfare(sc):
    f = os.path.join(SCENARIO_RUNS[sc], f'welfare_sc{sc}_{QUANTILE}{CAP_TAG}.csv')
    df = pd.read_csv(f)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return df.values


def load_hatpath(name, sc):
    f = os.path.join(SCENARIO_RUNS[sc], f'{name}_hat_path_sc{sc}_{QUANTILE}{CAP_TAG}.csv')
    df = pd.read_csv(f)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return df.values.reshape(-1, N) if df.shape[1] == 1 else df.values


HEADLINE_ARM = assert_one_arm(
    [f'welfare_sc{sc}_{QUANTILE}{CAP_TAG}.csv' for sc in (1, 2, 3)],
    'headline (sc1,sc2,sc3)')
W = {sc: load_welfare(sc) for sc in (1, 2, 3)}
T, N = W[1].shape
print(f"[load] welfare: T={T} rows (transitions), N={N} cells | "
      f"RUNS_DIR={os.path.basename(RUNS_DIR)} | CAP_TAG suffix='{CAP_TAG}' | "
      f"ACTUAL arm (manifest): {ARM_LABEL}")
assert T == 16, f"16 transitions expected (T=17 periods), got {T}"

# The three margins, named as in the paper. Note that the paper scenario numbers
# and the file scenario numbers are not the same:
#   annual = Scenario 1 vs baseline, the annual-mean response  = sc3 / sc1
#   daily  = Scenario 2 vs baseline, the daily response        = sc2 / sc1
#   extra  = Sc.2 vs Sc.1      (increment daily)                = sc2 vs sc3
MARGINS = {'daily': (2, 1), 'annual': (3, 1), 'extra': (2, 3)}

##################################################################################
# LES DEUX NUMEROTATIONS, ET LA REGLE QUI LES SEPARE
#
# Le PAPIER numerote Sc.1 = annuel, Sc.2 = journalier. Le CODE numerote sc2 =
# journalier, sc3 = annuel, avec sc1 pour la baseline. Les deux sont figees :
# renumeroter le code invaliderait onze runs archives, renumeroter le papier
# deplacerait une convention publiee.
#
# ⇒ REGLE : aucune legende, aucune note, aucun en-tete GENERE ne porte un
#   numero de scenario. Ils disent « daily arm » et « annual arm », qui sont
#   vrais dans les deux numerotations. Trois briefs de suite ont apparie un
#   numero au mauvais bras ; le numero est le probleme, pas le brief.
#   La table ci-dessous existe pour que le code puisse TRADUIRE quand il doit
#   lire un manifeste ou un nom de fichier, jamais pour composer du texte.
##################################################################################
CODE_TO_PAPER = {2: 'daily', 3: 'annual'}
PAPER_TO_CODE = {v: k for k, v in CODE_TO_PAPER.items()}
ARM_PHRASE = {'daily': 'the daily arm', 'annual': 'the annual arm'}
CUM = {m: np.cumsum(np.log(W[a] / W[b]), axis=0) for m, (a, b) in MARGINS.items()}


# ======================================================================
# 2. Pixel metadata: id, lat/lon, population, 1990 GDP per capita, country
# ======================================================================
kept = pd.read_csv(os.path.join(INPUTS, 'indices_drop.csv')).iloc[:, 0].astype(int).values
assert len(kept) == N, f"indices_drop ({len(kept)}) != N welfare ({N})"

xy = pd.read_csv(os.path.join(INPUTS, 'identifier_xy_final.csv'))
xy.columns = [c.strip().strip('"') for c in xy.columns]
id_xy = dict(zip(xy['id'].astype(int), zip(xy['y'], xy['x'])))

pop_df = pd.read_csv(os.path.join(INPUTS, 'population0.csv'))
id_pop = dict(zip(pop_df['id'].astype(int), pop_df['value']))

pix = pd.DataFrame({'id': kept})
pix['y'] = [round(id_xy[i][0], 1) for i in kept]
pix['x'] = [round(id_xy[i][1], 1) for i in kept]
pix['pop'] = [id_pop[i] for i in kept]

# Provenance guardrail: population0 of the run grid has a minimum of 1.0052 over
# the retained pixels. Another copy in circulation is floored much higher, so
# this check identifies the file rather than trusting its name.
# Lot 53 : 1.0052 est un RESULTAT (le domaine a change, 13 990 -> 13 986). En
# reprise deliberee la garde RAPPORTE (tools/reprise.py) ; hors reprise elle bloque.
# Son remplacement par une comparaison de deux calculs du meme passage est soumis.
sys.path.insert(0, os.path.join(PKG, 'tools'))
from reprise import reprise_declaree, rapporte   # noqa: E402
import hashlib   # noqa: E402  (E4, lot 74 : md5 de migration_matrix.rds)
import gardes_lot74 as _g74   # noqa: E402


def _garde_phrase(ok, garde, attendu, obtenu, message):
    """Lot 74 : garde de PHRASE, meme semantique que garde_phrase() de
    emit_values_empirics.py. Hors reprise : elle ARRETE. En reprise declaree :
    elle rapporte (ecarts_reprise_gardes.tsv) et l'emission continue."""
    if ok:
        return
    if reprise_declaree():
        rapporte('emit_values.py', garde, attendu, obtenu)
        return
    raise SystemExit(message)


if abs(pix['pop'].min() - 1.0052) >= 1e-3:
    if reprise_declaree():
        rapporte('emit_values.py', 'min population0 over the retained cells',
                 '1.0052', '%.4f' % pix['pop'].min())
    else:
        raise AssertionError(
            f"unexpected population0.csv (min={pix['pop'].min():.4f}, 1.0052 expected). "
            "Three vintages of this file exist; this check identifies the one the "
            "run grid used.")

for m in MARGINS:
    pix[m] = CUM[m][TERMINAL_ROW, :]

# 1990 GDP per pixel, aggregated from a half-degree raster to one degree
with rasterio.open(os.path.join(PKG, 'data_inputs', 'Kummu',
                                'rast_gdpTot_1990_2022_30arcmin.tif')) as r:
    band = np.nan_to_num(r.read(1), nan=0.0)
gdp_1deg = band.reshape(180, 2, 360, 2).sum(axis=(1, 3))


def sample_gdp(y, x):
    row = min(max(int(np.floor(90 - y)), 0), 179)
    col = min(max(int(np.floor(x + 180)), 0), 359)
    return gdp_1deg[row, col]


pix['gdp_1990'] = [sample_gdp(a, b) for a, b in zip(pix['y'], pix['x'])]

##################################################################################
# PAYS D'UNE CELLULE : AGREGATION PROPORTIONNELLE (lot 51, 28/09/2026)
#
# ⛔ Ce bloc prenait l'argmax, par cellule, de `share_population` de
# lengths_country.csv (annee 1995) sous le commentaire « dominant country per
# pixel ». Or cette colonne etait normalisee PAR PAYS : l'argmax designait le pays
# dont la cellule porte la plus grande part de SA propre population -- en
# frontiere, le plus petit (RECAP 3.229). Et la population comptee etait celle de
# la cellule ENTIERE pour chaque pays (3.231).
#
# Desormais une cellule contribue a CHAQUE pays en proportion de sa part de
# population (poids_par_cellule.csv, `part_de_la_cellule`, somme 1 par cellule,
# periode 2015). `pix` reste UNE ligne par cellule -- plusieurs blocs le
# remplissent par POSITION -- et n'a plus de colonne GID_0 : tout usage oublie
# echoue au lieu de reprendre l'argmax. par_pays() deplie une table de cellules en
# (cellule x pays) et multiplie `pop` et `gdp_1990` par la part, juste avant un
# regroupement par pays. Les agregats mondiaux, nordiques et par latitude sont
# inchanges par construction : les parts d'une cellule somment a 1.
##################################################################################
sys.path.insert(0, os.path.join(PKG, 'tools'))
from temoins_poids import temoin_part_de_la_cellule   # noqa: E402
PARTS = pd.read_csv(os.path.join(INPUTS, 'poids_par_cellule.csv'), keep_default_na=False)
PARTS = PARTS[PARTS['GID_0'] != ''].copy()
temoin_part_de_la_cellule(PARTS)
PARTS['x'] = PARTS['x'].round(1)
PARTS['y'] = PARTS['y'].round(1)
PARTS = pix[['id', 'x', 'y']].merge(PARTS[['x', 'y', 'GID_0', 'part_de_la_cellule']],
                                    on=['x', 'y'])[['id', 'GID_0', 'part_de_la_cellule']]
print(f"[meta] pixels without a country: {N - PARTS['id'].nunique()} / {N}; "
      f"(cell, country) pairs: {len(PARTS)}")


def par_pays(df):
    """Deplie une table de CELLULES (colonne `id`) en lignes (cellule x pays) :
    `pop` et `gdp_1990` sont multiplies par part_de_la_cellule ; les autres
    colonnes (marges, latitude) sont celles de la cellule."""
    d = df.merge(PARTS, on='id', how='inner')
    for col in ('pop', 'gdp_1990'):
        if col in d.columns:
            d[col] = d[col] * d['part_de_la_cellule']
    return d


def _pays(d, iso):
    """Selecteur d'unite PAYS : les lignes (cellule x pays) du pays `iso`, poids
    = population de la cellule x part_de_la_cellule. Les unites monde, nord et
    latitude restent sur les cellules et ne passent pas par ici."""
    e = par_pays(d)
    return e[e['GID_0'] == iso]


# ======================================================================
# 3. Aggregations
# ======================================================================
def wmean(df, col, wcol='pop'):
    if df[wcol].sum() <= 0:
        return np.nan
    return float(np.average(df[col], weights=df[wcol]))


def margins_of(df):
    return {m: wmean(df, m) for m in MARGINS}


def _band_pop(unit):
    """Population de l'unite, sur la meme grille `pix` que les marges. Sert au
    seul controle de structure de la bande."""
    if unit == 'world':
        return float(pix['pop'].sum())
    if unit == 'north':
        return float(pix.loc[pix['y'] >= NORTH_LAT, 'pop'].sum())
    if unit in ('latFF', 'latA'):
        return float(pix.loc[pix['y'] >= BAND_LAT[unit], 'pop'].sum())
    pp = par_pays(pix[['id', 'pop']])
    return float(pp.loc[pp['GID_0'] == UNIT_ISO[unit], 'pop'].sum())


world = margins_of(pix)
north = margins_of(pix[pix['y'] >= NORTH_LAT])

nrt = pix[pix['y'] >= NORTH_LAT]
north_pop_loss = 100.0 * nrt.loc[nrt['extra'] < 0, 'pop'].sum() / nrt['pop'].sum()
north_pop_M = nrt['pop'].sum() / 1e6

# world population share in loss, by margin, at the terminal horizon
world_pop_loss = {m: 100.0 * pix.loc[pix[m] < 0, 'pop'].sum() / pix['pop'].sum()
                  for m in MARGINS}

# World aggregate at the three horizons
agg_horizon = {}
for label, row in HORIZON_ROW.items():
    agg_horizon[label] = {m: float(np.average(CUM[m][row, :], weights=pix['pop']))
                          for m in MARGINS}
# Amplification ratio: the daily extra as a share of the annual effect, at the
# terminal horizon.
# Amplification ratio: the daily extra as a share of the annual effect, at the
amp_ratio = 100.0 * abs(world['extra']) / abs(world['annual'])

# ---- countries ----
cty = (par_pays(pix)
          .groupby('GID_0')
          .apply(lambda g: pd.Series({
              'pop': g['pop'].sum(), 'gdp': g['gdp_1990'].sum(),
              **margins_of(g)}), include_groups=False)
          .reset_index())
cty['gdppc'] = cty['gdp'] / cty['pop'].replace(0, np.nan)
cty = cty[(cty['pop'] > 0) & (cty['gdppc'] > 0) & np.isfinite(cty['gdppc'])]

# ---- quartiles of 1990 GDP per capita, population-weighted ----
cty = cty.sort_values('gdppc').reset_index(drop=True)
pcum = np.cumsum(cty['pop'].values)
breaks = [cty['gdppc'].values[np.searchsorted(pcum, q * pcum[-1])] for q in (.25, .5, .75)]


def tier(g):
    return ('Q1' if g <= breaks[0] else 'Q2' if g <= breaks[1]
            else 'Q3' if g <= breaks[2] else 'Q4')


cty['tier'] = cty['gdppc'].map(tier)
tiers = {q: margins_of(cty[cty['tier'] == q]) for q in ('Q1', 'Q2', 'Q3', 'Q4')}

# TEMOIN DES QUARTILES (lot 51). Les populations des quatre quartiles somment a
# celle des pays classes, et chacun porte environ le quart. Tolerance, posee avant
# la mesure : la frontiere d'un quartile tombe sur un pays ENTIER, donc l'ecart
# d'un quartile a 25 % ne peut depasser la part du plus grand pays. Au-dela, la
# ponderation des quartiles n'est pas celle qu'on croit.
_tp = cty.groupby('tier')['pop'].sum()
_tot = cty['pop'].sum()
_pmax = cty['pop'].max() / _tot
_parts = {t: _tp.get(t, 0.0) / _tot for t in ('Q1', 'Q2', 'Q3', 'Q4')}
print('[temoin quartiles] population shares: %s | largest country %.1f %% | '
      'classified countries %.4g out of %.4g pixel population'
      % (', '.join('%s %.1f %%' % (t, 100 * v) for t, v in _parts.items()),
         100 * _pmax, _tot, pix['pop'].sum()))
if abs(_tp.sum() / _tot - 1) > 1e-12 or any(abs(v - .25) > _pmax for v in _parts.values()):
    raise SystemExit('[STOP] quartile witness: the quartiles do not split the '
                     'population of the classified countries as expected.')

# ---- endogenous income channel, by latitude band ----
# welfare_extra = amen + income (income = residu GE endogene, wage x labor).
## ORIENTATION. The hat-paths are written by the solver in (T, N), one row per
## transition and one column per cell, exactly like welfare. The whole matrix is
## read and transposed here.
## Reading them flattened and reshaping to (N, T) would interleave cells and
## time, and would raise no error whenever the dimensions happen to agree.
def load_hatpath_NT(name, sc):
    arr = pd.read_csv(os.path.join(RUNS_DIR,
          f'{name}_hat_path_sc{sc}_{QUANTILE}{CAP_TAG}.csv')).values   # (T, N)
    if arr.shape != (T, N):
        raise SystemExit(f"[STOP] {name}_hat_path_sc{sc}: shape {arr.shape}, "
                         f"({T}, {N}) expected.")
    return arr.T                          # (N, T) : orientation attendue en aval

    
amen_h = {sc: load_hatpath_NT('amenity', sc) for sc in (2, 3)}
cum_amen = np.cumsum(np.log(amen_h[2] / amen_h[3]), axis=1)[:, -1]   # per-cell (N,)
pix['chan_income'] = pix['extra'] - cum_amen

# Latitude bands: robust to orientation, and the only publishable channel object
LAT_BANDS = [('Trop', -90, 30), ('Mid', 30, 45), ('Sub', 45, 60), ('Arc', 60, 90)]
chan_lat = {tag: wmean(pix[(pix['y'] >= lo) & (pix['y'] < hi)], 'chan_income')
            for tag, lo, hi in LAT_BANDS}

# ---- bornes de latitude ----
lat_rows = []
for thr in (40, 45, 50, 55, 60, 66.5):
    s = pix[pix['y'] >= thr]
    lat_rows.append((thr, s['pop'].sum() / 1e6, *[wmean(s, m) for m in ('annual', 'extra', 'daily')]))

# CAP INVARIANCE -- RETIREE le 11/09/2026. `cap_moves` et les deux macros
# qu'il alimentait decrivaient une comparaison entre un bras C=1,5 et un bras
# C=1,2. output/ge_runs_manifest.csv porte cap_colmass=none sur les onze runs :
# il n'y a qu'un bras, et il est sans cap. Rien a restaurer ici -- il faudrait
# d'abord des runs cappes sur la grille.

# ---- latitude bound of the reversal ---------------------------------------
# EMITTED from lat_rows, that is from tab:ge_latitude itself, which the
# manuscript already designates as the source.
# lat_rows = (threshold, pop_M, annual, extra, daily)
_neg = [r for r in lat_rows if r[4] < 0]
_pos = [r for r in lat_rows if r[4] >= 0]
lat_last_neg = max(r[0] for r in _neg) if _neg else float('nan')
lat_first_pos = min(r[0] for r in _pos) if _pos else float('nan')
print(f"[lat] daily total negative up to {lat_last_neg}N inclusive; "
      f"positive from {lat_first_pos}N")


# ======================================================================
# 4. THE SIGN OBSERVATION on the Northern aggregate
#



print("[obs] Sign of the Northern reversal (contribution B) — OBSERVATION, NON-blocking")
# A guardrail NEVER encodes a conclusion. That the North reverses is a scientific
# result which can legitimately change; letting it kill the emitter would
# guarantee that the pipeline gets "repaired" until it comes back, that is, that
# the false result is reinstated.
#
# Hence the polarity below: the hard stop belongs to NUMERICAL drift, which is a
# bug, and never to a qualitative property of the result. The sign of the
# Northern aggregate is therefore reported as an OBSERVATION and blocks nothing,
# while the magnitudes are asserted.
#
# Country-level REVERSAL/AMPLIFICATION labels are deliberately not invariants.
# Classifying a country on the sign of a fourth decimal is fragile, and some
# countries sit at the noise level. The latitude gradient replaced that
# classification.
print(f"  north.extra = {north['extra']:+.4f} — "
      f"{'negative sign, as in the published result' if north['extra'] < 0 else 'POSITIVE SIGN: the result has changed. This is not a pipeline failure; it is a fact to examine.'}")

# THE NUMERICAL NON-REGRESSION BLOCK HAS BEEN REMOVED, and nothing replaces
# it. It compared six aggregates against a reference arm; that comparison is
# no longer made.
#
# WHAT THAT COSTS, stated because a silent removal is worse than an absent
# check: a run of this emitter that ends without error no longer certifies
# its magnitudes. It certifies that the inputs were found, that the shapes
# match, and that the sign of the Northern aggregate is what the printed
# observation says. Numerical drift in the welfare matrices passes through.
#
# The sign of the Northern aggregate is reported above as an OBSERVATION and
# blocks nothing, by design: a hard stop on a qualitative property of the
# result would guarantee that the pipeline gets "repaired" until the property
# comes back, that is, that a false result is reinstated. That polarity is
# unchanged. What has gone is the numerical side of the pair.


# ======================================================================
# 4bis. NULL AMENITY SPECIFICATION (lower bound g = 0), loaded ONLY in the
#       headline arm. It serves the Manski identification range: the amenity
#       channel is not point-identified from bilateral migration, so the effect
#       is bounded between
#         - lower bound: g^LB = 0, welfare flowing through productivity alone
#         - upper bound: the interior estimate
#
#   THE QUESTION THIS BLOCK ANSWERS: does the Northern reversal survive g = 0?
#   If it does, the reversal does NOT depend on the unidentified amenity
#   channel; only the MAGNITUDE moves, not the SIGN.
# ======================================================================
NULL_TAG  = '_NULLamen'
HAVE_NULL = (CAP_TAG == '' and
             os.path.exists(os.path.join(RUNS_DIR, f'welfare_sc1_q50{NULL_TAG}.csv')))


# Latitude bands carrying a macro of their own, beyond the >=45N aggregate.
# Same thresholds as lat_rows, hence the same object as tab:ge_latitude.
BAND_LAT = {'latFF': 55.0, 'latA': 66.5}


def aggregates_for(tag, quant=None):
    """Recharge le welfare avec un suffixe et renvoie world / north / lat / cty
    (marges terminales) plus la part de population gagnante, en REUTILISANT la
    meta-grille `pix` deja construite. Aucune re-derivation de pixel : meme
    masque, meme pop, meme lat, meme pays, meme cumul, meme periode terminale.

    `quant` : quantile PAR SCENARIO, {1: 'q50', 2: 'q05', 3: 'q05'}. Absent, les
    trois scenarios prennent le QUANTILE global -- le comportement d'avant, que
    le bras null continue d'utiliser tel quel.

    (!) SC1 RESTE EN q50 PARTOUT. Le scenario 1 est le contrefactuel sans
    rechauffement : main.py y force les dommages a 1 (l. 349-350) et QUANTILE
    n'a aucun effet sur lui. Il n'existe qu'en q50 sur la grille, et ce n'est
    pas une omission. Sans consequence pour la marge `extra` (sc2/sc3), qui ne
    fait pas intervenir sc1."""
    if quant is None:
        quant = {sc: QUANTILE for sc in (1, 2, 3)}
    Wt = {}
    for sc in (1, 2, 3):
        d = pd.read_csv(os.path.join(RUNS_DIR, f'welfare_sc{sc}_{quant[sc]}{tag}.csv'))
        d = d.loc[:, ~d.columns.str.contains('^Unnamed')]
        if d.values.shape != (T, N):
            raise SystemExit(f"[STOP] welfare_sc{sc}_{quant[sc]}{tag}: shape "
                             f"{d.values.shape}, ({T}, {N}) expected.")
        Wt[sc] = d.values
    cum = {m: np.cumsum(np.log(Wt[a] / Wt[b]), axis=0) for m, (a, b) in MARGINS.items()}
    q = pix[['id', 'y', 'x', 'pop']].copy()
    for m in MARGINS:
        q[m] = cum[m][TERMINAL_ROW, :]
    w = {m: wmean(q, m) for m in MARGINS}
    n = {m: wmean(q[q['y'] >= NORTH_LAT], m) for m in MARGINS}
    lat = {tg: {m: wmean(q[q['y'] >= thr], m) for m in MARGINS}
           for tg, thr in BAND_LAT.items()}
    c = (par_pays(q).groupby('GID_0')
           .apply(lambda g: pd.Series({m: wmean(g, m) for m in MARGINS}),
                  include_groups=False))
    # Part de la population MONDIALE qui GAGNE, par marge : complement exact de
    # world_pop_loss, calcule sur les memes colonnes.
    gain = {m: 100.0 - 100.0 * q.loc[q[m] < 0, 'pop'].sum() / q['pop'].sum()
            for m in MARGINS}

    ##################################################################
    # LES OBJETS DU CORPS, SUR CE BRAS
    #
    # Jusqu'au 16/09/2026 cette fonction ne rendait que les marges, et tout le
    # reste -- quartiles, horizons, parts de population, canal de revenu,
    # bornes de latitude -- n'existait que sur le bras de tete. Les macros du
    # corps n'avaient donc PAS LE CHOIX : elles sortaient du bras interieur
    # meme quand le corps tournait a f_a = 0. On construit ici les memes
    # objets, sur le meme masque `pix`, pour n'importe quel bras.
    ##################################################################
    horizon = {lab: {m: float(np.average(cum[m][row, :], weights=q['pop']))
                     for m in MARGINS}
               for lab, row in HORIZON_ROW.items()}
    pop_loss = {m: 100.0 * q.loc[q[m] < 0, 'pop'].sum() / q['pop'].sum()
                for m in MARGINS}
    nq = q[q['y'] >= NORTH_LAT]
    npop_loss = 100.0 * nq.loc[nq['extra'] < 0, 'pop'].sum() / nq['pop'].sum()
    npop_M = nq['pop'].sum() / 1e6
    amp = 100.0 * abs(w['extra']) / abs(w['annual'])

    # Quartiles : l'APPARTENANCE au quartile depend du PIB/tete de 1990, qui
    # ne depend pas du bras. On reprend donc l'affectation du module plutot
    # que de la recalculer, ce qui garantit que les deux bras comparent les
    # memes groupes de pays et non deux partitions voisines.
    tier_map = dict(zip(cty['GID_0'], cty['tier']))
    qc = par_pays(q)
    qc['tier'] = qc['GID_0'].map(tier_map)
    tiers_arm = {t: {m: wmean(qc[qc['tier'] == t], m) for m in MARGINS}
                 for t in ('Q1', 'Q2', 'Q3', 'Q4')}

    # Canal de revenu endogene : extra moins la part d'amenite, par cellule.
    ah = {}
    for sc in (2, 3):
        f = os.path.join(RUNS_DIR,
                         f'amenity_hat_path_sc{sc}_{quant[sc]}{tag}.csv')
        arr = pd.read_csv(f).values
        if arr.shape != (T, N):
            raise SystemExit(f"[STOP] {os.path.basename(f)}: shape "
                             f"{arr.shape}, ({T}, {N}) expected.")
        ah[sc] = arr.T
    chan = q['extra'].values - np.cumsum(np.log(ah[2] / ah[3]), axis=1)[:, -1]
    qch = q.copy()
    qch['chan_income'] = chan
    chan_arm = {tg: wmean(qch[(qch['y'] >= lo) & (qch['y'] < hi)], 'chan_income')
                for tg, lo, hi in LAT_BANDS}

    lrows = []
    for thr in (40, 45, 50, 55, 60, 66.5):
        sl = q[q['y'] >= thr]
        lrows.append((thr, sl['pop'].sum() / 1e6,
                      *[wmean(sl, m) for m in ('annual', 'extra', 'daily')]))

    rows_c = [(nm, c.loc[iso]) for iso, nm in COLD_ORDER_ISO]
    counts = {'n': len(rows_c),
              'gain': [nm for nm, r in rows_c if r['annual'] > 0],
              'loser': [nm for nm, r in rows_c if r['annual'] <= 0],
              'extra_nonneg': [nm for nm, r in rows_c if r['extra'] >= 0]}

    return {'world': w, 'north': n, 'lat': lat, 'cty': c, 'gain': gain,
            'horizon': horizon, 'pop_loss': pop_loss, 'tiers': tiers_arm,
            'north_pop_loss': npop_loss, 'north_pop_M': npop_M,
            'amp_ratio': amp, 'chan_lat': chan_arm, 'lat_rows': lrows,
            'counts': counts}


if HAVE_NULL:
    NULL_ARM = assert_one_arm(
        [f'welfare_sc{sc}_{QUANTILE}{NULL_TAG}.csv' for sc in (1, 2, 3)],
        'null arm g=0')
    NULLagg = aggregates_for(NULL_TAG)
    ##################################################################
    # BODY : LE BRAS DU CORPS. Le corps du papier tourne a la borne basse
    # d'amenite, donc TOUTE macro qu'il consomme sort d'ici. Les objets du
    # bras interieur gardent leurs noms de module (world, north, cty, tiers,
    # agg_horizon...) et ne sont plus lus que pour les macros *Qfifty et pour
    # les tables de sensibilite de S5.
    ##################################################################
    BODY = NULLagg
    print("[null ] spec g=0 loaded (amenity lower bound) -- ARM OF THE BODY")
    # --- VERDICT ---
    nz = NULLagg['north']
    print(f"  VERDICT Northern reversal under g=0:")
    print(f"    north.extra  base={north['extra']:+.4f}  null={nz['extra']:+.4f}"
          f"  -> {'NEGATIVE (reversal HOLDS)' if nz['extra'] < 0 else 'POSITIVE (BROKEN)'}")
    print(f"    north.net    base={north['daily']:+.4f}  null={nz['daily']:+.4f}"
          f"  -> {'net loss kept' if nz['daily'] < 0 else 'NO MORE NET LOSS'}")
    if not (nz['extra'] < 0 and nz['daily'] < 0):
        print("  (!) WARNING: the Northern reversal does NOT hold under g=0. "
              "Major information — to raise before any regeneration.")
    # --- Northern range: net aggregate bounded by the two specifications ---
    _nb, _nn = north['daily'], nz['daily']
    sens_north_lo, sens_north_hi = min(_nb, _nn), max(_nb, _nn)
    # --- share of magnitude preserved, Northern aggregate, three margins ---
    _pres = [100.0 * nz[m] / north[m] for m in ('annual', 'daily', 'extra')
             if north[m] != 0]
    sens_pres_lo, sens_pres_hi = min(_pres), max(_pres)
    # --- share preserved AT COUNTRY LEVEL (extra margin, the object of the
    #     the cold economies, France excluded as erratic. LOWER than the
    #     aggregate: the amenity share by country exceeds the aggregate share.
    _cty_pres = []
    for iso, _ in REVERSAL + AMPLIFICATION:
        b_ext = cty[cty['GID_0'] == iso].iloc[0]['extra']
        z_ext = NULLagg['cty'].loc[iso]['extra']
        if b_ext != 0:
            _cty_pres.append(100.0 * z_ext / b_ext)
    sens_pres_cty_lo, sens_pres_cty_hi = min(_cty_pres), max(_cty_pres)
    ##########################################################################
    # DE COMBIEN LE BRAS INTERIEUR EST PLUS PROFOND QUE LA BORNE. 20/09/2026.
    #
    # Avant la division par theta, le net nordique passait de -0,043 a -0,091
    # entre les deux bras -- un facteur 2,1 -- et le supplement pouvait ecrire
    # que le contraste nordique ne dependait pas du canal. Apres division il
    # passe de -0,0067 a -0,055, un facteur 8,3 : l'amenite etait en unites
    # fondamentales pendant que la productivite portait theta, donc le bras
    # interieur SOUS-PONDERAIT l'amenite d'un facteur 6,5. La phrase tombe, et
    # ces deux macros disent DE COMBIEN plutot que de la remplacer par un
    # adjectif.
    #
    # RAPPORT DE MAGNITUDES, DONC SIGNE DES DEUX VERIFIE : un rapport de deux
    # quantites de signes opposes ne serait pas une profondeur relative, et la
    # phrase qui l'emploie serait fausse.
    ##########################################################################
    def _depth_ratio(nom, base, interieur):
        ## Lot 57 : garde de PHRASE (« profondeur relative »), pas de calcul. En
        ## reprise declaree : rapportee avec ses deux valeurs, le rapport est emis.
        if base * interieur <= 0 and reprise_declaree():
            rapporte('emit_values.py', '%s: "relative depth" (arms of the same sign)' % nom,
                     'same sign', 'f_a=0 %.5f, interior %.5f' % (base, interieur))
            return abs(interieur) / abs(base) if base != 0 else float('nan')
        if base * interieur <= 0:
            raise SystemExit(
                '[STOP] %s: the two arms have opposite signs '
                '(%.5f at f_a=0, %.5f at the interior estimate). A ratio '
                'of magnitudes would no longer describe a relative depth.'
                % (nom, base, interieur))
        return abs(interieur) / abs(base)
    sens_net_ratio = _depth_ratio('GEsensNetRatio', nz['daily'], north['daily'])
    sens_extra_ratio = _depth_ratio('GEsensExtraRatio', nz['extra'], north['extra'])
    print('  interior arm / lower bound: net %.2fx, extra %.2fx'
          % (sens_net_ratio, sens_extra_ratio))
    print(f"  Northern net range: [{sens_north_lo:+.3f}, {sens_north_hi:+.3f}]"
          f"  | aggregate magnitude preserved: {sens_pres_lo:.0f}-{sens_pres_hi:.0f}%"
          f"  | countries: {sens_pres_cty_lo:.0f}-{sens_pres_cty_hi:.0f}%")
    # THE WORLD daily-extra AGGREGATE CHANGES SIGN under g=0. The robust object
    # is the North, not the world. It is reported as it stands in the sensitivity
    # table and has to be stated in prose.
    print(f"  [note] world.extra base={world['extra']:+.4f}  null={NULLagg['world']['extra']:+.4f}"
          f"  ({'SIGN CHANGE' if world['extra'] * NULLagg['world']['extra'] < 0 else 'same sign'})")
else:
    NULLagg = None
    NULL_ARM = None
    BODY = None
    print("[null ] spec g=0 missing from this arm -> macros/range left PROVISOIRE")


# ======================================================================
# 4ter. QUANTILES OF THE MC BAND ON THE DOSE-RESPONSES  (q05 / q50 / q95)
#
#   The grid runs six solves at three quantiles of the Monte Carlo band on the
#   dose-responses. Until this block only q50 was read, so four of the six runs
#   carried no exhibit. They are read HERE, in ONE pass: three separate runs of
#   the emitter under QUANTILE=q05 then q95 would give three values.tex, not
#   three columns of one table.
#
#   ---------------------------------------------------------------------------
#   WHAT A "q05 CURVE" IS, MEASURED IN THE PRODUCER, NOT ASSUMED
#   ---------------------------------------------------------------------------
#   lib_drf.R::mc_summarise computes `quantile(y, .05) ... by = x`: the column is
#   a POINTWISE quantile across the MC draws, taken independently at each
#   temperature. NO SINGLE DRAW TRACES THAT CURVE. The q05 column is an envelope,
#   not a feasible dose-response.
#
#   Consequence, and it governs the wording of the caption: a functional of an
#   envelope is NOT a quantile of the functional's distribution. Welfare computed
#   on the q05 curve is not the 5th percentile of welfare. The three solves
#   bracket a range of assumed dose-response severity; they do not bound a
#   confidence interval, and MONOTONICITY OF THE RESULT IS NOT A PROPERTY TO BE
#   EXPECTED. Where it holds it is a fact about these runs, not a guarantee.
#
#   The four DRF tables moreover come from FOUR INDEPENDENT DRAW SEQUENCES --
#   00_parameters.R l.207, MC_SEED = {prod_daily 44, prod_annual 45,
#   amenity_daily 42, amenity_annual 43}, whose own comment states that "two arms
#   of the same channel do not share a draw sequence". main.py passes the SAME
#   COLUMN LABEL to all of them (l. 365-366 for sc2, l. 408-409 for sc3, via
#   QUANTILE_PROD / QUANTILE_AMENITY which both default to QUANTILE). What the
#   runs share is therefore a LEVEL OF SEVERITY, not a draw index: there is no
#   draw-level covariance between the curves to preserve, by construction.
#
#   ---------------------------------------------------------------------------
#   WHY THE BAND COLUMN CARRIES THE DAILY EXTRA
#   ---------------------------------------------------------------------------
#   `extra` (sc2/sc3) is the contrast object of the paper, and it is the only
#   margin whose numerator and denominator are evaluated at the SAME severity
#   level, so the severity assumption differences out of the comparison between
#   metrics. `daily` (sc2/sc1) has sc1 pinned at q50 by construction, so its
#   three columns mix a varying numerator with a fixed denominator.
#
#   The Northern daily TOTAL is not dropped for that: the question the six runs
#   were run to answer -- does the North still lose on the daily total across the
#   three solves? -- is answered in the table note, as a RANGE, from the macros
#   emitted below.
#
#   ---------------------------------------------------------------------------
#   THE CONSTRAINT THAT MATTERS
#   ---------------------------------------------------------------------------
#   Every column goes through aggregates_for(), the function the null arm already
#   uses: same cell mask, same population weights, same cumulation, same terminal
#   period. Nothing is re-implemented. The closure check below PROVES it rather
#   than asserting it -- the q50 column is recomputed through the band path and
#   must land on the module-level aggregates to floating-point identity.
#
#   (!) The magnitude assertions of this emitter pin the MEDIAN arm and are NOT
#   re-run on q05/q95: they would fail for the wrong reason.
# ======================================================================
BAND_Q = ('q05', 'q50', 'q95')

# The margin the band column carries. See the block header.
BAND_MARGIN = 'extra'

# The units, in table order. Each maps to the macro that already carries its q50
# on BAND_MARGIN, so a band macro is unambiguously attached to its own median.
# The convention is the package's own -- GEagg*Hundred for the world,
# GEnorthDailyExtra, GEcty<Nm>DailyExtra -- measured against values.tex.
# LES NOMS SONT CEUX QUE LE MANUSCRIT APPELLE DEJA -- \GEworldExtraQlo,
# \GEnorthExtraQfifty, ... Les sites d'appel font foi ; un nom invente ici
# obligerait le texte a bouger pour suivre l'emetteur, ce qui est l'inverse de
# la regle. Mesure : check_macros.py listait ces noms comme appeles et non
# definis.
BAND_UNITS = [
    ('world',  'GEworldExtra'),         # world, terminal horizon
    ('north',  'GEnorthExtra'),         # >= 45N
    ('latFF',  'GElatFFExtra'),         # >= 55N
    ('latA',   'GElatAExtra'),          # >= 66.5N
    ('ctyRus', 'GEctyRusExtra'),
    ('ctyFin', 'GEctyFinExtra'),
    ('ctyCan', 'GEctyCanExtra'),
    ('ctyDeu', 'GEctyDeuExtra'),
]
# Le q50 est emis pour TOUTES les unites sous le suffixe Qfifty : le manuscrit
# cite \GEworldExtraQfifty et \GEnorthExtraQfifty a cote de leurs bornes, et
# une bande dont le median porte un autre nom se relit mal.
BAND_NEW_Q50 = tuple(u for u, _m in BAND_UNITS)
UNIT_ISO = {'ctyRus': 'RUS', 'ctyFin': 'FIN', 'ctyNor': 'NOR',
            'ctySwe': 'SWE', 'ctyCan': 'CAN', 'ctyDeu': 'DEU'}


def _band_val(agg, unit, margin):
    """Marge d'une unite, depuis un agregat rendu par aggregates_for(). Une
    seule fonction pour les trois colonnes et pour les deux marges."""
    if unit == 'world':
        return agg['world'][margin]
    if unit == 'north':
        return agg['north'][margin]
    if unit in agg['lat']:
        return agg['lat'][unit][margin]
    return float(agg['cty'].loc[UNIT_ISO[unit]][margin])


HAVE_BAND = (CAP_TAG == '' and QUANTILE == 'q50' and all(
    os.path.exists(os.path.join(RUNS_DIR, f'welfare_sc{sc}_{qq}.csv'))
    for qq in BAND_Q for sc in (2, 3))
    and os.path.exists(os.path.join(RUNS_DIR, 'welfare_sc1_q50.csv')))

BAND = {}
BAND_OK = False
BAND_WHY = ''
BAND_NONMONO = []          # unites dont le q50 sort de [q05, q95] : RAPPORTE
if not HAVE_BAND:
    BAND_WHY = ("not the headline arm, or the six welfare_sc{2,3}_q{05,50,95}.csv "
                "are not all present")
    print("[bande] q05/q95 NOT loaded: %s -> macros left PROVISOIRE"
          % BAND_WHY)
else:
    # La bande lit six runs plus sc1 : un seul bras pour les sept, sinon la
    # colonne q05 et la colonne q95 ne sont pas comparables entre elles.
    _BAND_ARM = assert_one_arm(
        ['welfare_sc1_q50.csv']
        + [f'welfare_sc{sc}_{qq}.csv' for qq in BAND_Q for sc in (2, 3)],
        'band q05/q50/q95')
    BAND_AMENITY = _BAND_ARM['amenity'] if _BAND_ARM else None
    for qq in BAND_Q:
        BAND[qq] = aggregates_for('', {1: 'q50', 2: qq, 3: qq})
    print("[bande] q05/q50/q95 loaded | band column on margin '%s' "
          "(sc1 stays at q50, no effect on this margin)" % BAND_MARGIN)

    # ---- CONTROLE DE FERMETURE : le chemin de bande redonne le median -------
    # Si cette egalite tombe, la colonne q50 de la bande n'est pas le q50 du
    # papier, et les colonnes q05/q95 ne sont comparables a rien. C'est le seul
    # controle de ce bloc qui ARRETE : c'est une panne, pas un resultat.
    _clos = []
    for unit, _mac in BAND_UNITS:
        if unit in BAND_NEW_Q50:
            continue                      # aucune reference module-level
        ref = (agg_horizon[HORIZON_LABEL][BAND_MARGIN] if unit == 'world'
               else north[BAND_MARGIN] if unit == 'north'
               else float(cty[cty['GID_0'] == UNIT_ISO[unit]].iloc[0][BAND_MARGIN]))
        got = _band_val(BAND['q50'], unit, BAND_MARGIN)
        _clos.append((unit, ref, got, abs(ref - got)))
    _clos.append(('popGain', 100.0 - world_pop_loss[BAND_MARGIN],
                  BAND['q50']['gain'][BAND_MARGIN],
                  abs((100.0 - world_pop_loss[BAND_MARGIN])
                      - BAND['q50']['gain'][BAND_MARGIN])))
    _clos.append(('northDaily', north['daily'],
                  _band_val(BAND['q50'], 'north', 'daily'),
                  abs(north['daily'] - _band_val(BAND['q50'], 'north', 'daily'))))
    _worst = max(d for _u, _r, _g, d in _clos)
    print("[bande] closure q50 vs module aggregates: max gap %.3e over %d units"
          % (_worst, len(_clos)))
    if _worst > 1e-10:
        for u, r, g, d in _clos:
            if d > 1e-10:
                print("   DIVERGE  %-12s module=%+.6f  band=%+.6f  gap=%.3e"
                      % (u, r, g, d))
        sys.exit("STOP: the band path does not reproduce the median. The "
                 "q50 column of the table would not be the paper's.")
    BAND_OK = True

    # ---- STRUCTURE DE LA BANDE : RAPPORTEE, NON BLOQUANTE -------------------
    # La monotonie n'est PAS une propriete attendue : la colonne q05 est une
    # enveloppe ponctuelle, pas un tirage (cf. l'entete de ce bloc). Une unite
    # dont le median sort de [q05, q95] est une information sur la structure du
    # probleme, pas une panne -- elle est IMPRIMEE, et l'emission continue.
    print("[bande] structure, margin '%s' (Qlo = column q05, Qhi = column q95):"
          % BAND_MARGIN)
    for unit, _mac in BAND_UNITS + [('popGain', 'GEpopGainExtra')]:
        if unit == 'popGain':
            lo, md, hi = (BAND[qq]['gain'][BAND_MARGIN] for qq in BAND_Q)
        else:
            lo, md, hi = (_band_val(BAND[qq], unit, BAND_MARGIN) for qq in BAND_Q)
        inside = min(lo, hi) <= md <= max(lo, hi)
        if not inside:
            BAND_NONMONO.append((unit, lo, md, hi))
        print("   %-8s q05=%+.4f  q50=%+.4f  q95=%+.4f  |q95-q05|=%.4f  %s"
              % (unit, lo, md, hi, abs(hi - lo),
                 "median inside the interval" if inside
                 else "*** MEDIAN OUTSIDE [q05,q95] ***"))
    if BAND_NONMONO:
        print("[bande] (!) %d unit(s) where the median leaves [q05, q95]. "
              "Expected: a functional of a pointwise envelope has no "
              "reason to be ordered. To STATE in the note, not to hide."
              % len(BAND_NONMONO))
    else:
        print("[bande] median inside [q05, q95] for every unit of this draw "
              "-- a fact about these runs, not a guarantee.")

    # ---- LA QUESTION QUI A MOTIVE LES SIX RUNS ------------------------------
    # Le total journalier nord reste-t-il negatif sur les trois solves ?
    _nd = [_band_val(BAND[qq], 'north', 'daily') for qq in BAND_Q]
    NORTH_DAILY_LO, NORTH_DAILY_HI = min(_nd), max(_nd)
    NORTH_DAILY_ALLNEG = all(v < 0 for v in _nd)
    print("[bande] NORTHERN DAILY TOTAL over the three solves: "
          "q05=%+.4f  q50=%+.4f  q95=%+.4f  -> range [%+.3f, %+.3f], %s"
          % (*_nd, NORTH_DAILY_LO, NORTH_DAILY_HI,
             "NEGATIVE AT ALL THREE" if NORTH_DAILY_ALLNEG
             else "NOT negative everywhere -- a fact to examine before publication"))

    # ---- STRUCTURE : largeur contre taille de l'unite -----------------------
    # Rapporte, jamais bloquant. Le quantile est un choc COMMUN (une severite
    # assumee), pas du bruit d'echantillonnage : l'agregation ne moyenne rien, et
    # la largeur n'a aucune raison de decroitre avec la taille de l'unite.
    _w = sorted(((u, abs(_band_val(BAND['q95'], u, BAND_MARGIN)
                         - _band_val(BAND['q05'], u, BAND_MARGIN)), _band_pop(u))
                 for u, _m in BAND_UNITS), key=lambda t: t[2])
    print("[bande] width |q95-q05| by increasing population:")
    for u, wd, pp in _w:
        print("   %-8s pop=%12.0f  width=%.4f" % (u, pp, wd))


# ======================================================================
# 5. LaTeX emission
# ======================================================================
def _amenity_clause(arm=None):
    """Comment une table nomme la specification d'amenite DU BRAS QU'ELLE LIT.

    Lue dans output/ge_runs_manifest.csv via assert_one_arm(). Le manifeste
    laisse la colonne VIDE pour le bras null -- ce n'est pas une absence
    d'information, c'est ainsi que run_one_ge.sh l'enregistre -- donc le vide
    se lit f_a = 0 et rien d'autre n'est suppose.

    Le parametre `arm` est ce qui manquait jusqu'au 17/09/2026 : la clause
    etait calculee UNE FOIS sur HEADLINE_ARM puis collee sur TOUTES les
    tables, y compris les quatre que S5.1-S5.3 accompagne d'un texte dont
    les macros sortent, elles, du bras f_a = 0. Une table nomme desormais
    le bras dont ses propres cellules sortent.
    """
    a = ((HEADLINE_ARM if arm is None else arm) or {}).get('amenity')
    # assert_one_arm() remplit les vides par le SENTINEL '(vide=null)' avant de
    # dedupliquer : ce n'est pas une specification, c'est la facon dont
    # run_one_ge.sh enregistre le bras null. Le manquer faisait imprimer
    # « amenity specification 	exttt{(vide=null)} » dans quatre legendes.
    if a in (None, '', 'nan', '(vide=null)', '(vide)'):
        return "amenity lower bound $f_a = 0$"
    if str(a) == 'baseline':
        return "interior amenity estimate $" + chr(92) + "hat f_a^{" + chr(92) + ",q50}$"
    return "amenity specification \\texttt{%s}" % a


def _amenity_tag(arm):
    """L'etiquette BRUTE du manifeste. C'est elle que le commentaire en tete
    de chaque table emise porte, et que tools/check_macros.py (controle 8)
    relit pour verifier qu une table \\input dans un bloc a f_a = 0 sort
    bien d un run a f_a = 0. Le manifeste laisse la colonne VIDE pour le bras
    null : on ecrit `none`, le mot du contrat, et non la chaine vide, qui ne
    se relit pas."""
    a = (arm or {}).get('amenity')
    return ('none' if a in (None, '', 'nan', '(vide=null)', '(vide)')
            else str(a))


AMEN_CLAUSE = _amenity_clause()
# ----------------------------------------------------------------------
# LE BRAS DU CORPS, ET LES QUATRE TABLES QUI LE SUIVENT
# ----------------------------------------------------------------------
# Le corps du papier tourne a la borne basse f_a = 0. Les quatre tables de
# S5.1-S5.3 (aggregate, tiers, north, latitude) sont lues EN REGARD d un texte
# qui cite des macros emises depuis ce bras : elles doivent en sortir aussi.
# Jusqu au 17/09/2026 elles sortaient du bras interieur, et le lecteur avait
# sous les yeux une phrase et une cellule qui ne parlaient pas du meme modele.
# Les tables de S5.4 (sensitivity, north_sens) COMPARENT les deux bras : elles
# gardent AMEN_CLAUSE, qui nomme le bras interieur de leur colonne de gauche.
BODY_ARM    = NULL_ARM if HAVE_NULL else HEADLINE_ARM
BODY_CLAUSE = _amenity_clause(BODY_ARM)
BODY_TAG    = _amenity_tag(BODY_ARM)
HEAD_TAG    = _amenity_tag(HEADLINE_ARM)
print("[amenite] headline arm (S5.4, sensitivity): %s  [tag=%s]"
      % (AMEN_CLAUSE, HEAD_TAG))
print("[amenite] arm of the BODY (S5.1-S5.3)     : %s  [tag=%s]"
      % (BODY_CLAUSE, BODY_TAG))


# NO SIGN ON A ZERO -- FOR POINT VALUES ONLY (25/09/2026). A point value that
# rounds to zero prints 0.000: there the sign is noise. NEVER for an interval
# bound, whose sign says on which side of zero the interval lies: bounds go
# through fb(), which keeps the sign and adds decimals instead. The first
# version applied s3 to bounds too and erased the '+' of a positive q95.
def s3(v):
    t = f"{v:+.3f}"
    return t[1:] if float(t) == 0 else t


def f3(v):
    """Delta log U, VALEUR PONCTUELLE : 3 decimales, signe sauf sur zero."""
    return f"${s3(v)}$"


def fb(v):
    """Delta log U, BORNE d'intervalle : signe toujours ; si la valeur
    s'arrondit a zero a 3 decimales, on en ajoute jusqu'a ce que le chiffre
    imprime ne soit plus nul (6 au plus)."""
    for d in range(3, 7):
        t = f"{v:+.{d}f}"
        if float(t) != 0:
            return f"${t}$"
    return f"${v:+.6f}$"


def n3(v):
    """Trois decimales, signe toujours. LOT 74 : si la valeur s'arrondit a zero,
    on ajoute des decimales jusqu'a un chiffre non nul (6 au plus), comme fb()
    pour les bornes. Le 01/10, quatre macros s'imprimaient « -0.000 » / « +0.000 »
    (GEdailyMedNorthArctic -0.00014, GEdailyHiNorthFiftyFive +0.00038,
    GEtierQtwoExtra, GEchanIncomeMid) : un zero signe imprime. Toute autre valeur
    sort IDENTIQUE a l'ancien format. macro() refuse en outre tout zero signe."""
    return _g74.decimales_non_nulles(v, 3)


# FOUR DECIMALS, FOR THE DIFFERENCES THE READER ADDS UP (25/09/2026).
# Rule: a quantity the paper states to be the difference of two others prints
# at a precision where the difference closes. Family: GEnorthAnnual,
# GEnorthNet, GEnorthDailyExtra, GEnorthChanEval, GEnorthChanShape, at four
# decimals everywhere they appear. At three, +0.003 - 0.019 = -0.016
# while the total printed -0.017 (unrounded +0.002797 - 0.019421 =
# -0.016624): each term was rounded correctly and the displayed sum was
# false. At four it closes. Rounding a term the other way to force the sum
# is excluded: a published figure is the rounding of what was computed.
# Where one of them stands in a sentence next to other three-decimal values,
# the manuscript calls GEnorthExtraLB or GEnorthNetLB, the same numbers at
# three decimals.
def n4(v):
    return f"{v:+.4f}"


def n1(v):
    return f"{v:.1f}"


def _provenance_banner():
    """Bandeau emis EN TETE de values.tex : le runs_dir de CHAQUE famille de
    macros, et le millesime qu'il porte.

    Le bandeau precedent annoncait un unique " bras : RUNS_DIR=runs_cap_sweep "
    alors que le bras headline est compose par scenario et que trois familles
    lisent un autre millesime. Un fichier genere qui declare le mauvais bras est
    pire qu'un fichier genere qui n'en declare aucun : il donne au lecteur une
    provenance a laquelle se fier."""
    L = ["% " + "=" * 68,
         "% GENERATED by 02_GE/emit_values.py",
         "% DO NOT EDIT BY HAND. Every GE value in the paper comes from here.",
         "% terminal horizon = the last transition (" + HORIZON_BLOCK + "), labelled "
         + HORIZON_LABEL + ".",
         "% " + "-" * 68,
         "%% PROVENANCE BY MACRO FAMILY",
         "%%   CAP_TAG (filename suffix) : '%s'" % CAP_TAG,
         "%%   ARM AS RECORDED IN output/ge_runs_manifest.csv : %s" % ARM_LABEL,
         "%%   (the suffix names a file; only the manifest says what the solver did)",
         "%"]
    ## The cap and null arms are not produced by the current run grid. The
    ## families that read them are skipped when their files are missing, and the
    ## omission is ANNOUNCED: a published arm that disappears in silence is worse
    ## than one that is absent.
    if not _HYBRIDE:
        L += ["%% (!) HOMOGENEOUS OVERRIDE ARM: everything is read from %s."
              % os.path.basename(RUNS_DIR),
              "%   The vintages below describe the headline arm, not this one.",
              "%"]
    # LE MILLESIME EST LU, PAS GELE. Pour chaque famille, on va chercher au
    # manifeste le run_id des fichiers qu'elle lit REELLEMENT. Une chaine gelee
    # survit a un rejeu et atteste alors d'un run qui n'existe plus.
    _man = None
    _mp = os.path.join(OUTPUT_DIR, 'ge_runs_manifest.csv')
    if os.path.exists(_mp):
        try:
            _man = pd.read_csv(_mp).set_index('file')
        except Exception:
            _man = None

    def _vintage(files):
        if not files:
            return "via decompose_metric_channels.py (md5 pinned)"
        if _man is None:
            return "UNKNOWN (ge_runs_manifest.csv missing)"
        ids, hors = set(), []
        for f in files:
            b = os.path.basename(f)
            if b in _man.index:
                ids.add(str(_man.loc[b, 'run_id']))
            else:
                hors.append(b)
        ## Lot 86f: the run_id is the name of the machine's run directory
        ## (run_<timestamp>), so writing it made values.tex differ from one
        ## machine to the next. Only the COUNT is written: one run_id means one
        ## vintage, more than one is a mix. The run_ids themselves stay in
        ## output/ge_runs_manifest.csv and values_provenance.tsv.
        out = ("%d run_id for the %d file(s) read (named in ge_runs_manifest.csv)"
               % (len(ids), len(files) - len(hors))) if ids else "no file in the manifest"
        if hors:
            out += " (+%d outside the manifest: %s)" % (len(hors), ", ".join(hors))
        return out

    mixtes = []
    for fam, macros, files, mill, note in RUN_FAMILIES:
        dirs = sorted({os.path.basename(os.path.dirname(f)) for f in files})
        mill = _vintage(files)
        L.append("%% %-16s %s" % (fam, macros))
        L.append("%%   runs_dir : %s"
                 % (", ".join(dirs) if dirs
                    else "via decompose_metric_channels.py (md5 pinned)"))
        L.append("%%   vintage  : %s" % mill)
        if False:
            mixtes.append(fam)
        L.append("%")
    # Le paragraphe « ACKNOWLEDGED AND DECLARED MIX » a ete RETIRE le
    # 11/09/2026. Il decrivait un sc2 du 07/07 remplace par un autre, et trois
    # familles sur quatre lisant l'ancien. Le manifeste dit que tous les runs
    # lus portent un seul run_id : le melange decrit n'existe plus. Le bandeau
    # ci-dessus rend le run_id de chaque famille, mesure a chaque emission,
    # de sorte qu'un melange futur se verra sans qu'on ait a l'annoncer.
    if _man is not None:
        _all = sorted({str(v) for v in _man['run_id']})
        L += ["% " + "-" * 68,
              "%% RUN GRID : %d file(s) in the manifest, %d run_id (named in ge_runs_manifest.csv)"
              % (len(_man), len(_all)),
              "%%   cap_colmass : %s"
              % ", ".join(sorted({str(v) for v in _man['cap_colmass'].fillna('(vide)')})),
              "%%   amenity_spec: %s"
              % ", ".join(sorted({str(v) for v in
                                  _man['amenity_spec'].fillna('(empty = null)')}))]
    L.append("% " + "=" * 68)
    return "\n".join(L) + "\n"


HEADER = _provenance_banner()


##################################################################
# PROVENANCE PAR MACRO
#
# Une macro du CORPS doit sortir d'un run `amenity = none` : le corps du
# papier tourne a la borne basse. Cette regle n'etait ecrite nulle part, et
# trois macros du Nord ont porte les valeurs du bras interieur pendant des
# semaines sans que rien le signale -- au point que S5 imprimait « -0.153
# at the lower bound and -0.153 at the interior estimate ».
#
# Chaque appel a macro() enregistre donc DE QUEL BRAS la valeur sort. Le
# defaut est le bras de tete (celui des fichiers sans suffixe), ce qui est
# le bon defaut : c'est celui qu'on oublie de qualifier.
# tools/check_macros.py croise ce registre avec le .tex et echoue sur toute
# macro du corps issue d'un run interieur.
##################################################################
PROVENANCE = []


# Les CSV producteurs emis par le meme passage que les macros.
TABLES_CSV = {}


def macro(name, val, prov=None):
    """prov : cle de RUN_PROV, ou None pour le bras de tete."""
    # LOT 74 : une macro qui s'imprime « -0.000 » ou « +0.000 » n'est pas emise.
    # Garde d'EMISSION (format), bloquante meme en reprise.
    if _g74.zero_signe(val):
        raise SystemExit('[STOP] %s would print %r: a signed zero. Add '
                         'decimals (n3 does) or do not emit.' % (name, val))
    PROVENANCE.append((name, str(val), prov or 'headline'))
    return f"\\newcommand{{\\{name}}}{{{val}}}\n"


# ======================================================================
# GRID AND RUN CONSTANTS -- derived from the artefacts, never typed
# ======================================================================
# THE DENSITY CUT IS A TRAP AND IT IS DISARMED HERE.
# The manuscript describes a cut at "0.1 inhabitant per km2". That figure
# governs lambda, NOT the cell mask. 02a_calibration_at_1deg.r l.302-305 keeps
# the INTERSECTION of three conditions:
#       productivity > 0  AND  wage > 0  AND  population > 1   (per CELL)
# Measured 11/09/2026: population > 1 alone would keep 14,055 cells; the mask
# keeps 13,990. The binding condition is the positivity of wage and
# productivity, not the population cut. Emitting "0.1" would republish a claim
# the code does not make, so the macros below carry the RULE THAT RUNS, and the
# reconstruction is asserted rather than trusted.
def _grid_constants():
    out = {}
    pop_all = pd.read_csv(os.path.join(INPUTS, 'population0.csv'))
    out['potential'] = len(pop_all)
    out['pop_only'] = int((pop_all['value'] > 1).sum())
    got = None
    paths = {k: os.path.join(INPUTS, f'{k}0.csv') for k in ('wage', 'productivity')}
    if all(os.path.exists(p) for p in paths.values()):
        m = (pop_all.set_index('id')['value'] > 1)
        for k, p in paths.items():
            d = pd.read_csv(p).set_index('id')['value']
            m &= (d.reindex(m.index) > 0)
        got = sorted(m.index[m].astype(int))
    out['rebuilt'] = got
    return out


GRID = _grid_constants()
if GRID['rebuilt'] is not None:
    _same = (len(GRID['rebuilt']) == len(kept)
             and all(a == b for a, b in zip(GRID['rebuilt'], sorted(kept))))
    print("[grille] mask re-derived (prod>0 & wage>0 & pop>1): %d cells, "
          "identical to indices_drop.csv: %s"
          % (len(GRID['rebuilt']), "YES" if _same else "*** NO ***"))
    if not _same:
        sys.exit("STOP: the described retention rule (prod>0 & wage>0 & pop>1) "
                 "does not reproduce indices_drop.csv. A macro that DESCRIBES a "
                 "rule must be the one the code applies -- that is the density-cut "
                 "defect, and it is not republished.")
else:
    print("[grille] wage0/productivity0 missing: mask NOT re-derived, "
          "the rule macros stay PROVISOIRE.")
print("[grille] %d cells kept out of %d; the population threshold alone "
      "(>1/cell) would keep %d -- so it is NOT the binding constraint"
      % (N, GRID['potential'], GRID['pop_only']))


def _n_runs():
    p = os.path.join(OUTPUT_DIR, 'ge_runs_manifest.csv')
    if not os.path.exists(p):
        return None
    try:
        return int(len(pd.read_csv(p)))
    except Exception:
        return None


N_RUNS = _n_runs()


def emit_values_tex():
    s = HEADER
    s += "\n% ---- guardrail: any macro not emitted shouts ----\n"
    s += ("\\providecommand{\\PROVISOIRE}{\\textbf{??}"
          "\\PackageWarning{values}{PROVISOIRE: value not emitted}}\n")

    s += "\n% ---- grid and run constants (derived, never typed) ----\n"
    s += macro('GEncells', f"{N:,}".replace(',', '{,}'))
    s += macro('GEncellsPotential',
               f"{GRID['potential']:,}".replace(',', '{,}'))
    s += macro('GEnperiods', str(T + 1))
    s += macro('GEntransitions', str(T))
    # SSP2 IAMC : NOMBRE DE MODELES PAR PAYS (lot 61). 02a moyenne les modeles
    # disponibles pour chaque pays (OECD Env-Growth, IIASA GDP, et PIK GDP-32
    # pour USA seul) ; il ecrit output/ssp2_modeles_par_pays.csv. Comptes sur
    # les pays DES POIDS (dans_les_poids), ceux qui entrent dans la calibration.
    # Garde : la clause dit « un a trois modeles » ; un quatrieme, ou un pays a
    # zero, la rendrait fausse -> arret.
    _ssp = os.path.join(OUTPUT_DIR, 'ssp2_modeles_par_pays.csv')
    if os.path.exists(_ssp):
        _sm = pd.read_csv(_ssp)
        _sm = _sm[_sm['dans_les_poids'].astype(str).str.upper() == 'TRUE']
        if not set(_sm['n_modeles']) <= {1, 2, 3}:
            raise SystemExit("[STOP] ssp2_modeles_par_pays.csv: numbers of models "
                             "%s outside {1, 2, 3}." % sorted(set(_sm['n_modeles'])))
        for _k, _nm in ((1, 'GEnsspOneModel'), (2, 'GEnsspTwoModels'),
                        (3, 'GEnsspThreeModels')):
            s += macro(_nm, str(int((_sm['n_modeles'] == _k).sum())))
    else:
        for _nm in ('GEnsspOneModel', 'GEnsspTwoModels', 'GEnsspThreeModels'):
            s += macro(_nm, '\\PROVISOIRE')
    # The horizon LABEL, so the manuscript stops spelling a year out. It is a
    # label and not a computed year: the terminal row is the transition
    # HORIZON_BLOCK, and the model never evaluates a one-year 2100 window.
    # Source: 00_parameters.R / 00_parameters.py, shared with every producer.
    s += ("% GEhorizonLabel is an ETIQUETTE, not a year the model indexes on:\n"
          "% the terminal row is the transition " + HORIZON_BLOCK + ".\n")
    s += macro('GEhorizonLabel', HORIZON_LABEL)
    s += macro('GEhorizonSpan', HORIZON_SPAN)
    # The four GCMs of the climate projection, from CMIP6_MODELS (00_parameters
    # .py, twin of the R list that 00_temperature_projections_JMP.R checks the
    # directory against). Directory, code and text are thus one list. The file
    # prefixes upper-cased are the CMIP6 source_id spellings.
    _gcm = [m.upper() for m in _PARAMS.CMIP6_MODELS]
    s += macro('GEcmipModels', ', '.join(_gcm[:-1]) + ' and ' + _gcm[-1])
    # THE CUT THAT RUNS, not the one the prose described. See the block above.
    s += ("% GEdensityCut is the cut the CODE applies: population > 1 inhabitant\n"
          "% per CELL (02a_calibration_at_1deg.r l.302-305), intersected with\n"
          "% wage > 0 and productivity > 0. It is NOT the 0.1/km2 of the prose,\n"
          "% which governs lambda and not the cell mask, and it is NOT binding:\n"
          "% the population test alone would keep GEncellsPopOnly cells.\n")
    if GRID['rebuilt'] is not None:
        s += macro('GEdensityCut', '1')
        s += macro('GEdensityCutUnit', 'inhabitant per cell')
        s += macro('GEncellsPopOnly',
                   f"{GRID['pop_only']:,}".replace(',', '{,}'))
    else:
        for nm in ('GEdensityCut', 'GEdensityCutUnit', 'GEncellsPopOnly'):
            s += macro(nm, '\\PROVISOIRE')
    ##########################################################################
    # \GEtheta -- L'ELASTICITE DU COMMERCE, LUE A SA SOURCE UNIQUE.
    # Posee le 22/09/2026, avec la clause de l'equation (5).
    #
    # POURQUOI ELLE EXISTE. theta = 6,5 gouverne DEUX choses : le solveur
    # l'applique dans `trade_weights = (Z/x)^theta`, et l'estimation divise par
    # lui (`SCALE_PROD = 1/THETA_GE`). Depuis la clause du 22/09 il apparait en
    # plus a DEUX endroits du manuscrit -- §2.3 et le supplement A. Trois sites
    # pour une constante, dont deux en prose : c'est la configuration exacte ou
    # deux sites divergent, et c'est ce qui vient de se produire sur
    # l'ATTRIBUTION de cette meme valeur.
    #
    # ⛔ LUE, JAMAIS RECOPIEE. Meme discipline que 00_parameters.R l.111-146,
    # et pour la meme raison : recopier la valeur creerait une seconde source de
    # verite pour la constante dont le mauvais usage etait le defaut.
    # ⛔ AUCUN REPLI. Fichier introuvable ou cle absente => \PROVISOIRE, jamais
    # un nombre invente.
    ##########################################################################
    _theta_csv = os.path.join(HERE, '02c_model', 'specification_study.csv')
    _theta_val = None
    try:
        with open(_theta_csv, encoding='utf-8', errors='replace') as _f:
            _hits = [l for l in _f if l.strip().startswith('param_theta;')]
        if len(_hits) == 1:
            _t = float(_hits[0].split(';')[1])
            if _t > 0:
                # %g : 6.5 reste « 6.5 », un theta entier resterait entier.
                _theta_val = ('%g' % _t)
    except (OSError, ValueError, IndexError):
        _theta_val = None
    s += ("% GEtheta: param_theta of 02c_model/specification_study.csv, the\n"
          "% ONLY source of the constant. The solver applies it in\n"
          "% trade_weights = (Z/x)^theta; the estimation divides by it\n"
          "% (SCALE_PROD = 1/THETA_GE). The .tex must never copy it.\n")
    s += macro('GEtheta', _theta_val if _theta_val else '\\PROVISOIRE')
    ##########################################################################
    # ⚠️ SI alpha DEVIENT UNE MACRO UN JOUR : ELLE DESCEND DE `param_mu`, ET
    #   JAMAIS DE `alpha`. La contrainte est ecrite ICI, a cote de la seule
    #   autre constante de calibration emise, et pas dans une liste.
    #
    # LA COLLISION DE SYMBOLES, MESUREE LE 22/09/2026. Le supplement A ecrit
    # $\alpha=0.6$ pour la PART DU TRAVAIL, l'exposant du salaire dans
    # lambda_0 ∝ Z w^{-alpha*theta} tau^{-theta}. Or specification_study.csv
    # porte DEUX cles :
    #     param_mu;0.6; labor share (Desmet & Rossi-Hansberg, 2015)   <- celle-ci
    #     alpha;0.85;residential land share in consumption expenditure <- PAS celle-ci
    # Confirme dans le solveur : forward.py l.43-48,
    #     unit_price_bundle_hat_t(...) = w_hat * (ratio ** (1 - param_mu)),
    # donc `param_mu` EST la part du travail du coût unitaire.
    #
    # ⇒ Un symbole, deux objets, et le PIEGE est que le mauvais porte le nom
    #   que le papier emploie. C'est la configuration qui a produit le defaut
    #   theta, en plus vicieuse : la, deux sources disaient deux choses sur une
    #   meme constante ; ici une seule source dit deux choses sous un meme nom,
    #   selon qu'on lise le papier ou le code.
    #
    # alpha reste pour l'instant un LITTERAL de calibration dans le supplement
    # A, ce qui est admissible -- c'est un parametre IMPOSE, pas un resultat.
    ##########################################################################
    if N_RUNS is not None:
        s += macro('GEnruns', str(N_RUNS))
    else:
        s += macro('GEnruns', '\\PROVISOIRE')
    # WHAT THE ARM ACTUALLY IS, read from the manifest. See _arm_label_from_manifest.
    s += macro('GEarmLabel', ARM_LABEL)
    # La specification d'amenite du bras headline, MEME clause que les legendes
    # generees de S5. La legende de fig:ge_traj la cite au lieu de la retaper :
    # une figure et les tables qui l'entourent ne doivent pas pouvoir declarer
    # deux bras differents.
    s += macro('GEamenitySpec', AMEN_CLAUSE)
    # LA CLAUSE DU BRAS DU CORPS, distincte de celle du bras de tete. La
    # legende de fig:ge_traj appelait GEamenitySpec et annoncait donc le bras
    # INTERIEUR au-dessus d'un paragraphe qui cite des valeurs a f_a = 0 ;
    # depuis le 17/09/2026 la figure sort du bras du corps et sa legende le dit.
    s += macro('GEamenitySpecBody', BODY_CLAUSE)

    # ---- initial calibrated state (companion CSV, never the 867 MB matrix) ---
    s += "\n% ---- initial calibrated state (lambda_0, mu_0) ----\n"
    s += ("% Source: output/ge_calib_diagnostic.csv, written once by\n"
          "% 02_GE/ge_calibration_diagnostic.py, which IMPORTS the solver's own\n"
          "% loaders (load_h5, load_migration). The emitter never reads lambda.h5:\n"
          "% 867 MB and 37 s measured, for statistics that describe the INPUT\n"
          "% matrices and do not depend on the run.\n")
    _cal_p = os.path.join(OUTPUT_DIR, 'ge_calib_diagnostic.csv')
    if not os.path.exists(_cal_p):
        s += ("% *** NOT EMITTED: ge_calib_diagnostic.csv is missing. ***\n"
              "% Regenerate: python 02_GE/ge_calibration_diagnostic.py\n")
        print("[warn] ge_calib_diagnostic.csv missing -> 8 macros left "
              "PROVISOIRE (the calls in the .tex will stay undefined, which is "
              "the intended polarity).")
        for nm in ('GElamDiagMedian', 'GElamDiagMean', 'GElamDiagShareHigh',
                   'GEimportShareMean',
                   'GEclearMedian', 'GEclearGap', 'GEclearLo', 'GEclearHi',
                   'GEmigSenders', 'GEmigReceivers',
                   'GEdomCellsOut', 'GEdomPopOutM', 'GEdomPopOutShare', 'GEdomCtyOutN',
                   'GEdomCtyOutPopM', 'GEdomCtyPartialN', 'GEdomMigLostM'):
            s += macro(nm, '\\PROVISOIRE')
    else:
        _cal = pd.read_csv(_cal_p).set_index('quantity')['value']

        def _cal_get(k, fmt=None):
            """Rend PROVISOIRE quand le diagnostic dit NON CALCULE.

            Le CSV ecrit l'absence en toutes lettres plutot que de laisser une
            case vide : une case vide se lit comme un zero au premier read_csv
            distrait."""
            if k not in _cal.index:
                return None
            v = str(_cal.loc[k])
            if v.strip().upper().startswith('NON CALCULE'):
                return None
            return fmt % float(v) if fmt else v

        def _clear_gap():
            """|1 - mediane| en pour cent. None quand la mediane n'est pas
            calculee : la macro sort alors PROVISOIRE, comme ses voisines."""
            v = _cal_get('clear_median')
            return None if v is None else '%.0f' % abs(100.0 * (1.0 - float(v)))

        _vals = [
            ('GElamDiagMedian',    _cal_get('lam_diag_median', '%.2f')),
            ('GElamDiagMean',      _cal_get('lam_diag_mean', '%.2f')),
            ('GElamDiagShareHigh', _cal_get('lam_diag_share_gt_099', '%.1f')),
            # The mean import share, 100 x (1 - mean lambda_rr), in per cent and
            # unweighted. Computed unrounded by the producer: the complement of
            # the rounded mean would carry the rounding of lam_diag_mean.
            # One decimal: the grid is close to autarky and no more is claimed.
            ('GEimportShareMean',  _cal_get('import_share_mean_pct', '%.1f')),
            ('GEclearMedian',      _cal_get('clear_median', '%.2f')),
            # L'ECART AU CLEARING, en pour cent, derive de la MEME valeur
            # brute que la mediane affichee a cote de lui : les deux ne
            # peuvent donc pas se contredire dans la meme phrase.
            ('GEclearGap',         _clear_gap()),
            # Les deux bornes de l'inversion Eaton-Kortum : le paquet ne la
            # contient pas, cf. l'entete de ge_calibration_diagnostic.py.
            ('GEclearLo',          _cal_get('clear_q05')),
            ('GEclearHi',          _cal_get('clear_q95')),
            # L'OBJET DE CES DEUX LISTES (lot 74). mu_0 TEL QUE LE SOLVEUR LE LIT :
            # mig_flows.h5 via load_migration, restreint aux cellules ACTIVES du
            # run (indices_drop.csv), probabilites x population0 de la cellule
            # d'origine, cellules reparties entre pays par part_de_la_cellule.
            # Ce n'est PAS l'objet du temoin de 02b (STEP 15b, mig_witness_by_
            # country.csv) : flux bruts sur les 17 048 ids, avant la restriction
            # aux cellules actives. Les deux classements peuvent donc differer
            # (01/10/2026 : ARE ici, RUS au temoin, cinquieme recepteur).
            ('GEmigSenders',       _cal_get('mig_senders')),
            ('GEmigReceivers',     _cal_get('mig_receivers')),
            # LE DOMAINE (lot 68). Ce que la grille ne couvre pas : cellules
            # habitees absentes de la grille, leur population et sa part, les pays
            # sans aucune population active, les pays partiels, et les flux
            # internationaux que 02b STEP 13 retire avec ces cellules.
            ('GEdomCellsOut',      _cal_get('dom_cells_out')),
            ('GEdomPopOutM',       _cal_get('dom_pop_out_M', '%.1f')),
            ('GEdomPopOutShare',   _cal_get('dom_pop_out_share_pct', '%.1f')),
            ('GEdomCtyOutN',       _cal_get('dom_cty_out_n')),
            ('GEdomCtyOutPopM',    _cal_get('dom_cty_out_pop_M', '%.1f')),
            ('GEdomCtyPartialN',   _cal_get('dom_cty_partial_n')),
            ('GEdomMigLostM',      _cal_get('dom_mig_lost_M', '%.1f')),
        ]
        _prov = [nm for nm, v in _vals if v is None]
        # Provenance 'invariant' (25/09/2026): these statistics describe the
        # input matrices lambda0 and mu0 and descend from no run, so no amenity
        # arm applies. They carried the default 'headline', which maps to the
        # interior arm; check_macros control 5 refused GEimportShareMean once
        # the body cited it.
        for nm, v in _vals:
            s += macro(nm, v if v is not None else '\\PROVISOIRE', 'invariant')
        if _prov:
            s += ("%% %d macro(s) NOT EMITTED for lack of a producer: %s\n"
                  % (len(_prov), ', '.join(_prov)))
            print("[calib] %d macro(s) left PROVISOIRE: %s"
                  % (len(_prov), ', '.join(_prov)))
        print("[calib] lambda_rr median=%s mean=%s >0.99=%s%% | sales/income "
              "median=%s" % (_cal.get('lam_diag_median'), _cal.get('lam_diag_mean'),
                             _cal.get('lam_diag_share_gt_099'),
                             _cal.get('clear_median')))
    # ---- E4 (lot 74) : la masse transfrontaliere du bloc interne de mu_0 -------
    # Producteur : tools/mig_internal_cross.R (T1), qui recalcule STEP 10-11 de 02b
    # et s'assure contre 02b (A1 : masse STEP 11 ; A2 : soldes par pays de
    # mig_witness_by_country.csv a 1 personne pres, connu-faux Syrie). Objet :
    # masse hors diagonale de W'FW, niveau ZONE, avant domaine et plafond.
    # Provenance : le CSV porte le md5 de la migration_matrix.rds qu'il a lue ; une
    # matrice reecrite depuis rend le CSV perime, et l'emission s'arrete.
    _mx_p = os.path.join(INPUTS, 'mig_internal_cross.csv')
    if not os.path.exists(_mx_p):
        print("[warn] mig_internal_cross.csv missing -> GEmigInternalCross* PROVISOIRE "
              "(Rscript tools/mig_internal_cross.R --masse-step11 <M>).")
        for nm in ('GEmigInternalCrossM', 'GEmigInternalCrossShare'):
            s += macro(nm, '\\PROVISOIRE', 'invariant')
    else:
        _mx = pd.read_csv(_mx_p).set_index('quantity')
        _mx_md5 = str(_mx['source_md5'].iloc[0])
        _mm_p = os.path.join(INPUTS, 'migration_matrix.rds')
        _h = hashlib.md5()
        with open(_mm_p, 'rb') as _fh:
            for _blk in iter(lambda: _fh.read(1 << 24), b''):
                _h.update(_blk)
        if _h.hexdigest() != _mx_md5:
            raise SystemExit("[STOP] mig_internal_cross.csv was computed on a "
                             "migration_matrix.rds of md5 %s; the one on disk is %s. "
                             "Re-run tools/mig_internal_cross.R." % (_mx_md5, _h.hexdigest()))
        s += macro('GEmigInternalCrossM',
                   '%.2f' % (float(_mx.loc['internal_cross_gross', 'value']) / 1e6), 'invariant')
        s += macro('GEmigInternalCrossShare',
                   '%.1f' % float(_mx.loc['internal_cross_share_pct', 'value']), 'invariant')
        print("[mu0] cross-border mass of the internal block: %.2f M, %.1f %% "
              "(mig_internal_cross.csv, matrix %s)"
              % (float(_mx.loc['internal_cross_gross', 'value']) / 1e6,
                 float(_mx.loc['internal_cross_share_pct', 'value']), _mx_md5[:8]))
    s += macro('GEcapColmass', ARM_CAP if ARM_CAP else '\\PROVISOIRE')
    # NO PRODUCER IN THIS PACKAGE: the share of border cells needs the grid
    # construction, and no archived artefact carries it. Left PROVISOIRE so the
    # absence shouts at compile time instead of being filled with a guess.
    s += ("% GEborderShare: NOT EMITTED -- no archived artefact carries the\n"
          "% border-cell classification, and inventing a definition here would\n"
          "% be a new result, not an emission.\n")
    s += macro('GEborderShare', '\\PROVISOIRE')

    s += "\n% ---- North aggregate >= 45N (headline) ----\n"
    # ⚠️ CES TROIS MACROS ONT PORTE LE BRAS INTERIEUR JUSQU'AU 16/09/2026.
    # Le corps les presente comme la borne basse, S5 aussi, et elles valaient
    # +0,063 / -0,153 / -0,091 -- les colonnes q50 de tab_ge_north_sens. La
    # borne basse donne +0,065 / -0,109 / -0,044. Rien ne l'a signale parce
    # qu'aucune macro ne portait sa provenance ; c'est ce que PROVENANCE et
    # tools/check_macros.py ferment.
    # The three Northern margins print at four decimals together (see n4()):
    # the manuscript states daily extra = daily total - annual total, and at
    # three decimals -0.006 - (+0.010) = -0.016 while the extra printed -0.017.
    s += macro('GEnorthAnnual', n4(BODY['north']['annual']), 'null')
    s += macro('GEnorthDailyExtra', n4(BODY['north']['extra']), 'null')
    s += macro('GEnorthNet', n4(BODY['north']['daily']), 'null')
    s += macro('GEnorthPopLoss', n1(BODY['north_pop_loss']), 'null')
    s += macro('GEnorthPopM', f"{BODY['north_pop_M']:.0f}", 'null')
    # Les memes quantites sur le bras interieur, sous des noms qui le disent.
    # S5 et « Limits » les consomment ; le corps ne doit pas.
    s += macro('GEnorthAnnualQfifty', n3(north['annual']), 'headline')
    s += macro('GEnorthNetQfifty', n3(north['daily']), 'headline')
    # Lot 74 : §3.4 (« the annual arm itself projects a loss for the aggregate
    # above 45N ») et le supplement (amenite interieure) citent ces deux valeurs
    # comme NEGATIVES. Garde de phrase.
    _ok, _det = _g74.garde_amenite_ge(float(north['annual']), float(north['extra']))
    _garde_phrase(_ok, '§3.4 and supplement: at f_a^q50 the annual arm projects a loss for the North',
                  'annual < 0 and extra < 0', _det,
                  "[STOP] interior amenity: %s. The sentences of §3.4 and of the supplement "
                  "(\"projects a loss\") no longer hold." % _det)
    # ⚠️ PAS de GEnorthExtraQfifty ICI : la boucle de bande l'emet deja sous
    # `mac + 'Qfifty'` pour l'unite `north`. Une seconde definition faisait
    # echouer la compilation sur « Command already defined » -- ce qui est la
    # bonne polarite, mais la bonne reponse est de n'emettre qu'une fois.
    # Les deux valeurs sont le meme objet : le q50 de la bande est controle a
    # l'identite flottante contre les agregats du module.
    # EMITTED from lat_rows, that is from tab:ge_latitude.
    # The threshold grid is DISCRETE: the daily total is negative at one
    # threshold and positive at the next, so the crossing is BRACKETED, not
    # located. The bound macro carries the first positive threshold, which is a
    # bound and not a measurement. Both macros are emitted so that the sentence
    # can say so.
    # GEnorthLatBound RETIREE (11/09/2026). Elle portait la premiere borne
    # POSITIVE de la grille de seuils ; il n'y en avait pas, et elle etait emise
    # a `nan` -- un « nan » rendu dans le manuscrit. GEnorthLatAllNeg dit si la
    # grille est negative a tous les seuils.
    # ⚠️ LOT 74 (01/10/2026) : ELLE NE L'EST PLUS. Depuis la course du 01/10 (mu_0
    # corrige), le run central vaut +0,000076 a >= 66,5N (-0,000076 avant) :
    # GEnorthLatAllNeg = false, GEnorthLatLastNeg = 60. Le seuil 66,5 est « within
    # rounding of zero » (§3.2), ce que la garde ci-dessous verifie, et sa
    # population est GEnorthArcticPopM.
    _blr = BODY['lat_rows']
    _bneg = [r for r in _blr if r[4] < 0]
    _bpos = [r for r in _blr if r[4] >= 0]
    _blast = max(r[0] for r in _bneg) if _bneg else float('nan')
    s += macro('GEnorthLatLastNeg', f"{_blast:g}", 'null')
    s += macro('GEnorthLatAllNeg', 'true' if not _bpos else 'false', 'null')
    # Population au dernier seuil negatif : le manuscrit ecrivait « some two
    # million inhabitants » en dur a cote d'une macro qui valait nan.
    _lr_last = [r for r in _blr if r[0] == _blast]
    s += macro('GEnorthLatLastNegPopM',
               f"{_lr_last[0][1]:.1f}" if _lr_last else '\\PROVISOIRE', 'null')
    # The daily total AT the last negative threshold, four decimals (25/09/2026).
    # At 66.5 N it is about -2.5e-4, which tab:ge_latitude prints as 0.000 at
    # three: the sentence that calls the threshold negative cites this value so
    # that the table's zero does not read as a contradiction. Same lat_rows
    # row as GEnorthLatLastNeg, so the two cannot drift apart.
    if _lr_last:
        if not _lr_last[0][4] < 0:
            raise SystemExit('[STOP] GEnorthLatLastNegDaily: the last negative '
                             'threshold does not carry a negative total.')
        s += macro('GEnorthLatLastNegDaily', n4(_lr_last[0][4]), 'null')
    else:
        s += macro('GEnorthLatLastNegDaily', '\\PROVISOIRE', 'none')
    s += ("%% GEnorthLatAllNeg=%s : %s\n"
          % ('true' if not _bpos else 'false',
             'the daily total is negative at EVERY threshold of the latitude grid.'
             if not _bpos else
             'the daily total of the central run is non-negative at %s N.'
             % ', '.join('%g' % r[0] for r in _bpos)))
    # \GEnorthArcticPopM (lot 74) : population >= 66,5N, depuis la ligne 66,5 de
    # lat_rows -- meme passage que GEnorthLatLastNegPopM. Et la garde de §3.2 :
    # « negative at every latitude threshold of the grid, and within rounding of
    # zero above 66.5N » (lot 85q : clause « every threshold » = GEnorthLatAllNeg).
    _lr_arc = [r for r in _blr if abs(r[0] - 66.5) < 1e-9]
    if len(_lr_arc) != 1:
        raise SystemExit('[STOP] lat_rows: %d row(s) at threshold 66.5 (1 expected).' % len(_lr_arc))
    s += macro('GEnorthArcticPopM', f"{_lr_arc[0][1]:.1f}", 'null')
    _ok, _det = _g74.garde_seuil_arctique(float(_lr_arc[0][4]), not _bpos)
    _garde_phrase(_ok, '§3.2: negative at every threshold, and "within rounding of zero" >= 66.5N',
                  'all negative; < 0.0005', _det,
                  "[STOP] arctic threshold: %s. The sentence of §3.2 no longer holds." % _det)
    # ---- counts behind the tab_ge_north caption, emitted so the body can cite
    # them without retyping a country name that a re-run may invalidate.
    _nc = BODY['counts']
    s += macro('GEnorthNcty', str(_nc['n']), 'null')
    s += macro('GEnorthNGain', str(len(_nc['gain'])), 'null')
    s += macro('GEnorthNExtraNeg', str(_nc['n'] - len(_nc['extra_nonneg'])), 'null')
    s += macro('GEnorthExceptions',
               _enum(_nc['loser']) if _nc['loser'] else 'none', 'null')

    # ---- aggregation gap decomposition (evaluation / shape channels) ----
    # decompose_metric_channels refuses to run unless its three welfare matrices
    # come from one run grid, checked against output/ge_runs_manifest.csv.
    # Assembling them from runs that do not belong together reverses the
    # conclusion, so these values are never recomputed by hand here.
    # north['extra'] is passed as the closure check: the module re-derives the
    # same aggregate through its own mask, weights and cumulation, and must land
    # on this value.
    # ⛔ LA DECOMPOSITION N'EXISTE PAS A f_a = 0.
    # Elle a besoin de sc4 (la branche annuelle lisant la courbe journaliere) ;
    # le manifeste ne porte qu'un seul sc4, `welfare_sc4_q50_dailycurve.csv`,
    # a `amenity_spec = baseline`. Emettre sa valeur sous un nom nu remettrait
    # une valeur interieure dans le corps -- exactement le defaut qu'on ferme.
    # Elles sortent donc en PROVISOIRE, et la valeur interieure est emise a
    # part, sous un nom qui dit son bras.
    s += "\n% ---- aggregation gap decomposition (sc4) ----\n"
    from decompose_metric_channels import compute as _metric_channels
    # Le bras INTERIEUR, sous des noms qui le disent.
    _ch = _metric_channels(cross_check=north['extra'])
    s += macro('GEnorthChanEvalQfifty',  n3(_ch['evaluation']), 'headline')
    s += macro('GEnorthChanShapeQfifty', n3(_ch['shape']), 'headline')
    # Le bras DU CORPS, f_a = 0. Le controle croise est le daily extra du
    # corps : decompose_metric_channels re-derive le total par son propre
    # masque et ses propres poids, sans partager de code avec cet emetteur,
    # donc leur accord n'est pas une tautologie. Il refuse au-dela de sa
    # tolerance, et ce refus est le test que le brief demande.
    # ⚠️ L'ORDRE DES SUFFIXES : le bras vient AVANT « _dailycurve », parce que
    # main.py compose _sc<label><q> + <amenite> + <cap> + <drf>. Ce test
    # cherchait « _dailycurve_NULLamen » et ne trouvait donc jamais le run,
    # meme une fois celui-ci sur le disque. On demande le nom au producteur du
    # triplet plutot que de le reconstruire ici une seconde fois.
    from decompose_metric_channels import triplet_for as _triplet_for
    _sc4null = os.path.join(RUNS_DIR, _triplet_for(NULL_TAG)[2])
    if HAVE_NULL and os.path.exists(_sc4null):
        _chn = _metric_channels(cross_check=BODY['north']['extra'],
                                tag=NULL_TAG)
        s += macro('GEnorthChanEval',  n4(_chn['evaluation']), 'null')   # n4 : see n4()
        s += macro('GEnorthChanShape', n4(_chn['shape']), 'null')
        ##############################################################
        # LE TERME DE JENSEN : UN RAPPORT ET UN SIGNE, PAS UN ADJECTIF
        #
        # Le manuscrit ecrivait « of opposite sign and an order of magnitude
        # smaller ». Les deux moities dependent du BRAS et changent avec lui :
        # au bras interieur l'evaluation vaut -0,003 contre -0,150, donc MEME
        # signe que le total et facteur 50 ; a la borne elle vaut +0,017 contre
        # -0,125, donc signe OPPOSE et facteur 7,4. Un adjectif ne peut pas
        # suivre ca. On emet le rapport et le signe, et la phrase les lit.
        _r_ev, _r_sh = _chn['evaluation'], _chn['shape']
        if _r_ev == 0:
            raise SystemExit("[STOP] zero Jensen term: the ratio diverges.")
        s += macro('GEnorthChanRatio', f"{abs(_r_sh) / abs(_r_ev):.1f}", 'null')
        s += macro('GEnorthChanSameSign',
                   'of the same sign as' if _r_ev * _r_sh > 0
                   else 'of opposite sign to', 'null')
        s += macro('GEworldChanEval',  n3(_chn['world_evaluation']), 'null')
        s += macro('GEworldChanShape', n3(_chn['world_shape']), 'null')
        # La somme des deux canaux EST le daily extra, par telescopage exact.
        # On l'asserte ici aussi : le controle croise de decompose porte sur le
        # total re-derive, celui-ci porte sur la macro que le .tex imprime a
        # cote. Deux phrases differentes, deux controles.
        _sum = _chn['evaluation'] + _chn['shape']
        if abs(_sum - BODY['north']['extra']) > 5e-4:
            raise SystemExit(
                f"[STOP] decomposition at f_a=0: evaluation + shape = "
                f"{_sum:+.4f}, daily extra of the body = "
                f"{BODY['north']['extra']:+.4f}.")
    else:
        s += ("% PROVISOIRE: no sc4 run at amenity=none in the manifest.\n")
        for _nm in ('GEnorthChanEval', 'GEnorthChanShape',
                    'GEworldChanEval', 'GEworldChanShape',
                    'GEnorthChanRatio', 'GEnorthChanSameSign'):
            s += macro(_nm, '\\PROVISOIRE', 'none')

    s += "\n% ---- world aggregate, three horizons ----\n"
    s += ("% Thirty=2030, Fifty=2050, Hundred=" + HORIZON_LABEL + " (digits are\n"
          "% barred from a macro name). 'Hundred' is a NAME, kept stable across\n"
          "% every relabelling so no call site in the manuscript has to move.\n")
    for label in HORIZON_ROW:
        for m, tag in (('daily', 'Daily'), ('annual', 'Annual'), ('extra', 'Extra')):
            s += macro(f'GEagg{tag}{HORIZON_WORD[label]}',
                       n3(BODY['horizon'][label][m]), 'null')
    s += macro('GEampRatio', n1(BODY['amp_ratio']), 'null')
    # annual effect at 2030 in per cent, to match the literature comparator
    s += macro('GEaggAnnualThirtyPct',
               n1(100 * (np.exp(BODY['horizon']['2030']['annual']) - 1)), 'null')
    for m, tag in (('daily', 'Daily'), ('annual', 'Annual'), ('extra', 'Extra')):
        s += macro(f'GEpopLoss{tag}', n1(BODY['pop_loss'][m]), 'null')
    # complements: the share of world population that GAINS under each margin.
    # That is the object of the reclassification claim, not the losing share.
    for m, tag in (('daily', 'Daily'), ('annual', 'Annual')):
        s += macro(f'GEpopGain{tag}', n1(100 - BODY['pop_loss'][m]), 'null')

    s += "\n% ---- countries: reversals / amplifications / France ----\n"
    for iso, _ in REVERSAL + AMPLIFICATION + FRANCE:
        r = BODY['cty'].loc[iso]
        nm = MACRO_NAME[iso]
        s += macro(f'GEcty{nm}Net', n3(r['daily']), 'null')
        s += macro(f'GEcty{nm}Annual', n3(r['annual']), 'null')
        s += macro(f'GEcty{nm}DailyExtra', n3(r['extra']), 'null')

    s += "\n% ---- 1990 GDP per capita quartiles ----\n"
    for q in ('Q1', 'Q2', 'Q3', 'Q4'):
        qn = {'Q1': 'Qone', 'Q2': 'Qtwo', 'Q3': 'Qthree', 'Q4': 'Qfour'}[q]
        for m, tag in (('daily', 'Daily'), ('annual', 'Annual'), ('extra', 'Extra')):
            s += macro(f'GEtier{qn}{tag}', n3(BODY['tiers'][q][m]), 'null')

    # Q4 over Q1 ratio on the daily-extra margin. The multiplier is a result,
    # not a turn of phrase, so it is a macro: otherwise the text lies when the
    # run changes.
    s += macro('GEtierExtraRatio',
               f"{abs(BODY['tiers']['Q4']['extra']) / abs(BODY['tiers']['Q1']['extra']):.1f}",
               'null')

    s += "\n% ---- endogenous income channel by latitude band ----\n"
    s += "% The only publishable channel object: robust to the hat-path orientation.\n"
    s += "% Neither the amplification ratio (not robust) nor the quartile (wrong axis).\n"
    # Les hat-paths existent pour les trois scenarios du bras null, donc le
    # canal SE CALCULE a f_a = 0 : pas de PROVISOIRE ici.
    for _tg in ('Trop', 'Mid', 'Sub', 'Arc'):
        s += macro(f'GEchanIncome{_tg}', n3(BODY['chan_lat'][_tg]), 'null')

    # ---- source decomposition by exposure masks (25/09/2026) -----------------
    # Written by 02_GE/02i_exposure_mask_decomposition.py from four masked runs
    # (02g, 02h). A mask freezes exposure at 2015 outside its kept cells, so
    # warming is switched off there; general equilibrium is not. The macro names
    # say what is switched on: NorthWarming = warming at y >= 45 N only,
    # NonNorthWarming = warming below 45 N only. f_a = 0, central draw, one draw.
    # Emitted: what the text cites, plus the interaction residual (the largest
    # absolute share over the cited rows).
    s += "\n% ---- source decomposition by exposure masks (f_a = 0, central draw) ----\n"
    _mdp = os.path.join(OUTPUT_DIR, 'exposure_mask_decomposition.csv')
    _mcp = os.path.join(OUTPUT_DIR, 'exposure_mask_cells.csv')
    _MASK_NAMES = ('GEmaskNorthWarmingExtraShare', 'GEmaskNonNorthWarmingExtraShare',
                   'GEmaskNorthWarmingIncomeSubShare', 'GEmaskNorthWarmingIncomeArcShare',
                   'GEmaskNorthWarmingAnnual', 'GEmaskNorthWarmingDaily',
                   'GEmaskNonNorthWarmingDailyShare', 'GEmaskInteractionMaxShare',
                   'GEmaskNorthCells', 'GEmaskNonNorthCells')
    if not (os.path.exists(_mdp) and os.path.exists(_mcp)):
        s += ("% *** NOT EMITTED: exposure_mask_decomposition.csv or\n"
              "% exposure_mask_cells.csv is missing. Run the steps masks,\n"
              "% maskruns and maskdec of run_all.sh. ***\n")
        print("[warn] mask decomposition missing -> %d macros PROVISOIRE"
              % len(_MASK_NAMES))
        for _nm in _MASK_NAMES:
            s += macro(_nm, '\\PROVISOIRE', 'none')
    else:
        _md = pd.read_csv(_mdp)

        def _mg(reg, q, col):
            _row = _md[(_md['region'] == reg) & (_md['quantity'] == q)]
            if len(_row) != 1:
                raise SystemExit(f"[STOP] exposure_mask_decomposition.csv: "
                                 f"{len(_row)} row(s) for {reg}/{q}")
            return float(_row[col].iloc[0])

        # The full column must BE the body, not resemble it: a decomposition
        # computed on another vintage of the headline runs is refused.
        for _q, _k in (('daily_extra', 'extra'), ('daily_total', 'daily'),
                       ('annual_total', 'annual')):
            _gap = _mg('north_45', _q, 'full') - BODY['north'][_k]
            if abs(_gap) > 1e-9:
                raise SystemExit(f"[STOP] mask decomposition: full {_q} "
                                 f"differs from the body by {_gap:+.2e}. CSV of another "
                                 f"vintage: re-run 02i.")
        for _tg, _reg in (('Sub', 'band_45_60'), ('Arc', 'band_60_90')):
            _gap = _mg(_reg, 'income_part', 'full') - BODY['chan_lat'][_tg]
            if abs(_gap) > 1e-9:
                raise SystemExit(f"[STOP] mask decomposition: income "
                                 f"component {_reg} differs from GEchanIncome{_tg} "
                                 f"by {_gap:+.2e}.")
        _cited = (('north_45', 'daily_extra'), ('north_45', 'daily_total'),
                  ('north_45', 'annual_total'), ('band_45_60', 'income_part'),
                  ('band_60_90', 'income_part'))
        _imax = max(abs(_mg(r, q, 'share_interaction')) for r, q in _cited)
        _cells = pd.read_csv(_mcp).set_index('mask')['n_cells_kept']
        if int(_cells['north']) + int(_cells['south']) != N:
            raise SystemExit(f"[STOP] masks: {int(_cells['north'])} + "
                             f"{int(_cells['south'])} cells != {N}.")
        _thou = lambda v: f"{int(v):,}".replace(',', '{,}')
        _mv = [
            ('GEmaskNorthWarmingExtraShare',
             n1(_mg('north_45', 'daily_extra', 'share_north_warming_only'))),
            ('GEmaskNonNorthWarmingExtraShare',
             n1(_mg('north_45', 'daily_extra', 'share_nonnorth_warming_only'))),
            ('GEmaskNorthWarmingIncomeSubShare',
             n1(_mg('band_45_60', 'income_part', 'share_north_warming_only'))),
            ('GEmaskNorthWarmingIncomeArcShare',
             n1(_mg('band_60_90', 'income_part', 'share_north_warming_only'))),
            ('GEmaskNorthWarmingAnnual',
             n3(_mg('north_45', 'annual_total', 'north_warming_only'))),
            ('GEmaskNorthWarmingDaily',
             n3(_mg('north_45', 'daily_total', 'north_warming_only'))),
            ('GEmaskNonNorthWarmingDailyShare',
             n1(_mg('north_45', 'daily_total', 'share_nonnorth_warming_only'))),
            ('GEmaskInteractionMaxShare', f"{_imax:.2f}"),
            ('GEmaskNorthCells', _thou(_cells['north'])),
            ('GEmaskNonNorthCells', _thou(_cells['south'])),
        ]
        for _nm, _v in _mv:
            s += macro(_nm, _v, 'null')
        print("[masques] %s" % ' | '.join(f"{a}={b}" for a, b in _mv))
        ## G9 (lot 79, R7) : le paragraphe de §3.4 sur la composante revenu.
        ## Connu-faux : tools/test_gardes_lot78.py.
        import gardes_lot79 as _g79   # noqa: E402
        _ok, _det = _g79.garde_g9(
            float(BODY['chan_lat']['Sub']), float(BODY['chan_lat']['Arc']),
            float(BODY['chan_lat']['Mid']),
            _mg('band_45_60', 'income_part', 'share_north_warming_only'),
            _mg('band_60_90', 'income_part', 'share_north_warming_only'))
        print('[G9] %s' % _det)
        ## LOT 85o : re-ancree sur le nouveau §3.4 (bande temperee = \GEchanIncomeMid).
        _garde_phrase(_ok, '§3.4 (R7): income loss growing with latitude, North alone > 90 %',
                      'true', _det,
                      "[STOP] §3.4 (R7): %s. The sentence \"an income loss that grows with latitude "
                      "... the daily extra is <Mid> in the temperate band, <Sub> between 45 and "
                      "60N and <Arc> above 60N ... Nearly all of it is produced by warming in the "
                      "North itself\" is false." % _det)

    # The five CALIB_FROM_NOTE literals are no longer emitted (25/09/2026); see
    # the note where the dictionary used to be.
    # CALIBcapMagLo / CALIBcapMagHi RETIREES le 11/09/2026, avec la branche qui
    # les portait. Elles comparaient un bras C=1,5 a un bras C=1,2 ; le
    # manifeste dit cap_colmass=none sur les onze runs, donc aucun des deux
    # n'existe. `cap_moves` etait vide en dur et la branche ne pouvait pas
    # s'executer : une macro qu'aucun chemin n'emet et que personne n'appelle.

    # ---- diagnostic de colmass lambda0 / mu0 (BLOC 3, 04/08) ---------------
    # Produced from lambda.h5 and mig_flows.h5 -- NOT from a welfare run: these
    # statistics describe the input matrices, so they are invariant to the run
    # grid. The producer RECORDS the md5 of the two matrices in the CSV rather
    # than pinning them: the package regenerates both, so a pin certified
    # against another artefact would refuse every legitimate run.
    s += "\n% ---- colmass lambda0 / mu0 (app:ge_robustness) ----\n"
    if CAP_COLMASS_DIAG is None:
        # A mechanism that can do nothing must say that it did nothing. NO
        # fallback value is emitted: the absence must show up as an
        # "Undefined control sequence" at compile time rather than as a plausible
        # number.
        s += ("% *** NOT EMITTED: output/cap_colmass_diagnostic.csv is missing. ***\n"
              "% Regenerate: python 02_GE/cap_colmass_diagnostic.py\n")
        print("[warn] cap_colmass_diagnostic.csv missing -> 4 CALIB*Colmass macros "
              "NOT emitted (the calls in the .tex will stay undefined, which is "
              "the intended polarity).")
    else:
        d = CAP_COLMASS_DIAG
        lam = d[d['matrix'] == 'lambda0'].iloc[0]
        mu = d[d['matrix'] == 'mu0'].iloc[0]
        # Thousands separator `{,}`, the dominant convention of the manuscript
        # (13{,}990, 3{,}806{,}447). L3284 portait `2\,379`, seule occurrence en
        # thin space, harmonised here
        # Les deux MAXIMA decrivent les matrices et restent emis.
        # CALIBlamNcolsExceed / CALIBmuNcolsExceed sont RETIREES : elles
        # comptaient les colonnes au-dessus de C, et sans C la quantite n'a pas
        # de definition. Zero appel dans le manuscrit apres la reecriture de S2.
        s += macro('CALIBlamMaxColmass', f"{float(lam['max_colmass']):,.0f}"
                   .replace(',', '{,}'))
        s += macro('CALIBmuMaxColmass', f"{float(mu['max_colmass']):.1f}")

    s += "\n% ---- amenity range (lower bound g=0 against the interior g=q50) ----\n"
    if HAVE_NULL:
        s += "% NET North aggregate (daily margin), bounded between the two amenity specs.\n"
        s += "% Lo/Hi = numerical min/max (both bounds are net losses).\n"
        s += macro('GEsensNorthLo', n3(sens_north_lo))
        s += macro('GEsensNorthHi', n3(sens_north_hi))
        s += "% Share of the AGGREGATE North magnitude preserved under g=0 (min/max over the three margins).\n"
        s += macro('GEsensPreservedLo', f"{sens_pres_lo:.0f}\\%")
        s += macro('GEsensPreservedHi', f"{sens_pres_hi:.0f}\\%")
        s += "% Same AT COUNTRY LEVEL (extra margin, 6 cold economies): lower than the aggregate.\n"
        s += macro('GEsensPreservedCtyLo', f"{sens_pres_cty_lo:.0f}\\%")
        s += macro('GEsensPreservedCtyHi', f"{sens_pres_cty_hi:.0f}\\%")
        s += ('% How much deeper the INTERIOR arm is than the lower\n'
              '% bound, North >= 45N. Ratio of magnitudes of the same sign.\n')
        s += macro('GEsensNetRatio', '%.1f' % sens_net_ratio)
        s += macro('GEsensExtraRatio', '%.1f' % sens_extra_ratio)
        # ---- LOWER BOUND f_a = 0 against the interior f_a^q50 --------------
        # Les six macros que S5 appelle deja. Elles sortent des MEMES agregats
        # que tab_ge_sensitivity et tab_ge_north_sens -- NULLagg et les
        # agregats du module -- donc elles ne peuvent pas diverger des tables
        # qui les entourent. Elles remplacent les litteraux perimes de V10.
        s += "\n% ---- lower bound f_a=0 vs interior f_a^q50 (S5) ----\n"
        # ⚠️ LES NOMS SANS SUFFIXE SONT CEUX QUE LE MANUSCRIT APPELLE.
        # L'emetteur ecrivait GEworldDailyLB / GEworldAnnualLB ; le .tex appelle
        # GEworldDaily / GEworldAnnual (note NEW-MACRO l.666 : « f_a=0 »). Les
        # deux ont coexiste sans se rencontrer : les macros suffixees sortaient
        # et n'etaient jamais lues, les appels restaient indefinis. On emet les
        # noms APPELES. Ce sont les memes valeurs, et ce sont celles de la
        # Table 2 du corps, qui tourne a f_a = 0.
        s += macro('GEworldDaily',        n3(NULLagg['world']['daily']), 'null')
        s += macro('GEworldAnnual',       n3(NULLagg['world']['annual']), 'null')
        # La troisieme marge manquait : elle est la difference des deux, et elle
        # est calculee ICI plutot que soustraite dans le .tex.
        s += macro('GEworldExtra',        n3(NULLagg['world']['extra']), 'null')
        s += macro('GEnorthExtraLB',      n3(NULLagg['north']['extra']), 'null')
        # Same number as GEnorthNet (BODY = NULLagg) at three decimals, for the
        # sentence that sets it beside GEnorthNetQfifty (see n4()).
        s += macro('GEnorthNetLB',        n3(NULLagg['north']['daily']), 'null')
        s += macro('GEworldDailyQfifty',  n3(world['daily']), 'headline')
        s += macro('GEworldAnnualQfifty', n3(world['annual']), 'headline')
        # GEnorthNetQfifty / GEnorthAnnualQfifty / GEnorthExtraQfifty sont
        # emises avec le bloc du Nord, a cote des macros du corps qu'elles
        # doublent : c'est la seule facon de voir d'un coup d'oeil que les
        # deux jeux ne portent pas le meme bras.
    else:
        s += "% No null runs in this arm: macros left PROVISIONAL (a LaTeX warning).\n"
        for nm in ('GEsensNorthLo', 'GEsensNorthHi',
                   'GEsensPreservedLo', 'GEsensPreservedHi',
                   'GEsensPreservedCtyLo', 'GEsensPreservedCtyHi',
                   'GEworldDaily', 'GEworldAnnual', 'GEworldExtra',
                   'GEnorthExtraLB', 'GEnorthNetLB',
                   'GEworldDailyQfifty', 'GEworldAnnualQfifty',
                   'GEsensNetRatio', 'GEsensExtraRatio'):
            s += macro(nm, '\\PROVISOIRE', 'none')

    # ---- quantiles of the MC band on the dose-responses -------------------
    s += "\n% ---- MC band on the dose-responses: q05 / q50 / q95 ----\n"
    s += ("% WHAT Qlo AND Qhi NAME. Qlo is the aggregate evaluated on the q05\n"
          "% COLUMN of the dose-response tables, Qhi on the q95 column. They name\n"
          "% the quantile of the CURVE, not a sorted bound on this aggregate: the\n"
          "% columns are POINTWISE quantiles across MC draws (lib_drf.R,\n"
          "% mc_summarise: quantile(y,.05) by x), so no single draw traces them and\n"
          "% a functional of such an envelope is not a quantile of that\n"
          "% functional's distribution. The median is NOT guaranteed to lie\n"
          "% between Qlo and Qhi. Do not caption these as a confidence band.\n"
          "% MARGIN: daily extra (sc2/sc3) -- the only margin whose numerator and\n"
          "% denominator are evaluated at the same severity level, so the severity\n"
          "% assumption differences out of the comparison between metrics.\n"
          "% Same path as the median for every column: aggregates_for(), same cell\n"
          "% mask, same weights, same cumulation, same terminal period; the q50\n"
          "% column is checked against the module aggregates to fp identity.\n")
    if BAND_OK:
        s += ("%% AMENITY OF THE BAND RUNS, read from the manifest: %s.\n"
              % (BAND_AMENITY or 'UNKNOWN'))
        for unit, mac in BAND_UNITS:
            s += macro(mac + 'Qfifty', n3(_band_val(BAND['q50'], unit, BAND_MARGIN)))
            s += macro(mac + 'Qlo', n3(_band_val(BAND['q05'], unit, BAND_MARGIN)))
            s += macro(mac + 'Qhi', n3(_band_val(BAND['q95'], unit, BAND_MARGIN)))
        s += "% Share of world population GAINING on the daily-extra margin.\n"
        s += macro('GEpopGainExtra', n1(BAND['q50']['gain'][BAND_MARGIN]))
        s += macro('GEpopGainExtraQlo', n1(BAND['q05']['gain'][BAND_MARGIN]))
        s += macro('GEpopGainExtraQhi', n1(BAND['q95']['gain'][BAND_MARGIN]))
        # ---- the question the six runs were run to answer ------------------
        s += ("% RANGE of the Northern DAILY TOTAL across the three solves. Not a\n"
              "% band: the extremes of three runs, reported as an extent.\n")
        s += macro('GEnorthDailyRunMin', n3(NORTH_DAILY_LO))
        s += macro('GEnorthDailyRunMax', n3(NORTH_DAILY_HI))
        s += macro('GEnorthDailyRunAllNeg',
                   'negative in all three' if NORTH_DAILY_ALLNEG
                   else 'NOT negative in all three')
        # ---- units whose median falls outside [q05, q95] --------------------
        # Stated, never hidden: it is the observable signature of the fact that
        # these columns are envelopes rather than draws.
        if BAND_NONMONO:
            s += ("%% %d unit(s) with the median OUTSIDE [q05, q95]: %s\n"
                  % (len(BAND_NONMONO), ", ".join(u for u, *_ in BAND_NONMONO)))
            s += macro('GEbandNonMono',
                       ", ".join(u for u, *_ in BAND_NONMONO))
        else:
            s += macro('GEbandNonMono', 'none')
    else:
        s += ("%% NOT EMITTED: %s\n" % BAND_WHY)
        for unit, mac in BAND_UNITS:
            for suf in ('Qfifty', 'Qlo', 'Qhi'):
                s += macro(mac + suf, '\\PROVISOIRE')
        for nm in ('GEpopGainExtra', 'GEpopGainExtraQlo', 'GEpopGainExtraQhi',
                   'GEnorthDailyRunMin', 'GEnorthDailyRunMax',
                   'GEnorthDailyRunAllNeg', 'GEbandNonMono'):
            s += macro(nm, '\\PROVISOIRE')

    # ---- bande sur tirages classes ---------------------------------------
    if RANKED_OK:
        _rr, _rg = _ranked_rows()
        _by = {r['label']: r for r in _rr}
        _nr = _by[r'North $\geq 45^\circ$N']
        _wr = _by['World']
        s += "\n% ---- band on RANKED draws (rq05/rq50/rq95, f_a = 0) ----\n"
        s += ("% One actual draw per column, ranked by the Northern change in\n"
              "% partial equilibrium. Each arm in ITS own specification, so the\n"
              "% extra margin is BRACKETED (independent arms), not drawn.\n")
        _W = {'rq05': 'RqFive', 'rq50': 'RqFifty', 'rq95': 'RqNinetyFive'}
        for _p, _w in _W.items():
            s += macro(f'GEnorthDaily{_w}', n3(_nr['daily_' + _p]), 'null')
            s += macro(f'GEnorthAnnual{_w}', n3(_nr['annual_' + _p]), 'null')
            s += macro(f'GEworldDaily{_w}', n3(_wr['daily_' + _p]), 'null')
        s += macro('GEnorthExtraRqMin', n3(_nr['extra_min']), 'null')
        s += macro('GEnorthExtraRqFifty', n3(_nr['extra_q50']), 'null')
        s += macro('GEnorthExtraRqMax', n3(_nr['extra_max']), 'null')
        s += macro('GEpopGainDailyRqFive', n1(_rg['rq05']), 'null')
        s += macro('GEpopGainDailyRqNinetyFive', n1(_rg['rq95']), 'null')

        # LE SIGNE, DIT PAR LA MACRO ET NON PAR LA PROSE.
        _neg = [p for p in RANKED_Q if _nr['daily_' + p] < 0]
        if len(_neg) == 3:
            _txt = 'negative at all three percentiles'
        else:
            _pos = [p for p in RANKED_Q if _nr['daily_' + p] >= 0]
            _txt = ('not negative at ' + _enum([p.replace('rq', 'q')
                                                for p in _pos]))
        s += macro('GEnorthDailyRankedAllNeg', _txt, 'null')

        ##############################################################
        # TESTS D'ACCEPTATION DE LA BANDE
        ##############################################################
        # 1. L'identite des marges. Elle est vraie PAR CONSTRUCTION sur la
        #    colonne q50 (extra_q50 est pose comme la difference), donc ce
        #    controle ne teste pas l'arithmetique : il teste qu'on n'a pas
        #    interverti deux colonnes en chemin.
        for _r in _rr:
            _g = abs(_r['daily_rq50'] - (_r['annual_rq50'] + _r['extra_q50']))
            if _g > 5e-4:
                raise SystemExit(
                    f"[STOP] {_r['label']}: daily q50 != annual q50 + extra "
                    f"q50 (gap {_g:.2e}).")
            # Les etiquettes « min » et « max » doivent etre vraies SUR CHAQUE
            # LIGNE, y compris celles que le classement n'ordonne pas.
            if not (_r['extra_min'] <= _r['extra_q50'] <= _r['extra_max']):
                raise SystemExit(
                    f"[STOP] {_r['label']}: the bracket of the extra margin "
                    f"does not bracket its median "
                    f"({_r['extra_min']:+.4f}, {_r['extra_q50']:+.4f}, "
                    f"{_r['extra_max']:+.4f}).")
        # 2. LE CLASSEMENT DOIT SE VOIR EN EQUILIBRE GENERAL.
        #    Les tirages sont ordonnes par le changement nord en equilibre
        #    PARTIEL. Si l'ordre ne survit pas au GE, la colonne « q05 » ne
        #    nomme plus rien et la table serait un classement faux.
        _d = [_nr['daily_' + p] for p in RANKED_Q]
        ## Lot 57 : garde de PHRASE / d'etiquette de table (« q05 » nomme un rang),
        ## pas de calcul -- un GE peut legitimement ne pas conserver l'ordre partiel.
        ## LOT 85o (decision d'auteur, 07/10/2026) : la condition complete
        ## q05 <= q50 <= q95 devient un AVERTISSEMENT. Les macros \GEnorthDailyRq*
        ## qu'elle protegeait ne sont appelees ni par le .tex ni par une table, et le
        ## supplement (supp:ge_band_ranked) dit deja que l'ordre partiel ne survit
        ## qu'imparfaitement. Le tirage central peut tomber hors de [q05, q95] (course
        ## Cholesky : q05=-0.0121 q50=-0.0034 q95=-0.0044). Seul q05 <= q95 reste
        ## BLOQUANT : deux bornes inversees diraient que les colonnes sont interverties.
        ## Ni SystemExit ni ecriture dans output/ pour l'avertissement.
        if not (_d[0] <= _d[2]):
            raise SystemExit(
                "[HOLD] the Northern daily total has its ranked bounds reversed: "
                f"q05={_d[0]:+.4f} > q95={_d[2]:+.4f}. The ranking in partial "
                "equilibrium does not survive the move to general equilibrium.")
        if not (_d[0] <= _d[1] <= _d[2]):
            print("[WARNING] ranked table: the central draw of the Northern daily "
                  f"total leaves [q05, q95] in GE: q05={_d[0]:+.4f} q50={_d[1]:+.4f} "
                  f"q95={_d[2]:+.4f}. \\GEnorthDailyRq* not called; non-blocking.")
        # 3. LE CENTRE DE LA BANDE CONTRE LE CORPS -- RAPPORTE, PAS BLOQUANT.
        #    Le q50 CLASSE est un tirage ; les macros du corps sortent du run
        #    _NULLamen, qui lit la colonne de quantiles PONCTUELS. Les deux ne
        #    sont pas le meme objet et rien ne garantit qu'ils coincident : un
        #    ecart est une mesure de la distance entre une enveloppe ponctuelle
        #    et un tirage median, pas une anomalie. Il est donc EMIS.
        _cmp = (('GEnorthNet', _nr['daily_rq50'], BODY['north']['daily']),
                ('GEnorthAnnual', _nr['annual_rq50'], BODY['north']['annual']),
                ('GEnorthDailyExtra', _nr['extra_q50'], BODY['north']['extra']))
        _gmax = max(abs(a - b) for _, a, b in _cmp)
        print("[bande] ranked centre against the body (pointwise quantiles):")
        for _nm, _a, _b in _cmp:
            print(f"        {_nm:20s} ranked {_a:+.4f}  body {_b:+.4f}  "
                  f"gap {_a - _b:+.4f}")
        s += ("% Largest gap between the centre of the ranked band and the\n"
              "% macros of the body, which come from the column of POINTWISE\n"
              "% quantiles. The two objects differ by construction.\n")
        s += macro('GEnorthRankedCentreGap', n3(_gmax), 'null')
    else:
        # Les six runs ne sont pas la. Les macros sortent en PROVISOIRE plutot
        # que de manquer : une macro absente donne « Undefined control
        # sequence » et AUCUN PDF, donc plus rien n'est relisible ; un
        # PROVISOIRE donne « ?? » a l'endroit exact du trou. La table, elle,
        # reste absente et son \input arrete la compilation -- c'est voulu, et
        # c'est la polarite deja retenue pour tab_ge_band.
        s += "\n% ---- ranked-draw band: NOT EMITTED (runs rq* missing) ----\n"
        for _w in ('RqFive', 'RqFifty', 'RqNinetyFive'):
            for _b in ('GEnorthDaily', 'GEnorthAnnual', 'GEworldDaily'):
                s += macro(_b + _w, '\\PROVISOIRE', 'none')
        for _nm in ('GEnorthExtraRqMin', 'GEnorthExtraRqFifty',
                    'GEnorthExtraRqMax', 'GEpopGainDailyRqFive',
                    'GEpopGainDailyRqNinetyFive', 'GEnorthDailyRankedAllNeg',
                    'GEnorthRankedCentreGap'):
            s += macro(_nm, '\\PROVISOIRE', 'none')

    # ---- arm B : la distribution des sorties GE --------------------------
    if ARMB_OK:
        # LE CSV PRODUCTEUR, REMIS AU FLUX. Une ligne par tirage et par unite,
        # ecrite ICI, dans le meme passage que les macros, pour que les deux ne
        # puissent pas decrire deux etats de la grille. L'emission existait
        # jusqu'au lot 6 du 16/09 ; la refonte du bloc arm B l'a emportee, et
        # le fichier est reste sur disque a trente tirages -- livre tel quel.
        TABLES_CSV['ge_band_ranked_daily.csv'] = \
            assert_armb_csv(_armb_frame()).to_csv(index=False, lineterminator='\n')
        _st = _armb_stats()
        _inv, _rho, _nd = _armb_rank_corr()
        _n = _st['North>=45N']
        _w = _st['World']
        s += "\n% ---- arm B: distribution of the GE outputs (daily arm) ----\n"
        s += ("% These macros describe the DISTRIBUTION OF THE OUTPUTS, not three\n"
              "% chosen draws. The Rq* of the selected draws still exist\n"
              "% but no longer belong in the body: the ranking in partial\n"
              "% equilibrium carries over to them only weakly, and reading\n"
              "% them as percentiles of the outputs would be wrong.\n"
              "% check_macros.py refuses their call in the body.\n")
        s += macro('GEarmBnDraws', str(_nd), 'null')
        s += macro('GEarmBrankInversions', str(_inv), 'null')
        s += macro('GEarmBrankCorr', f"{_rho:.2f}", 'null')

        s += macro('GEnorthDailyRqMed', n3(_n['med']), 'null')
        s += macro('GEnorthDailyRqLo',  n3(_n['lo']),  'null')
        s += macro('GEnorthDailyRqHi',  n3(_n['hi']),  'null')
        s += macro('GEnorthDailyRqMin', n3(_n['min']), 'null')
        s += macro('GEnorthDailyRqMax', n3(_n['max']), 'null')
        s += macro('GEnorthDailyShareNeg', of_neg(_n), 'null')

        # Par unite : la part negative ET la bande, sous le meme suffixe.
        for _lab, _mac, _ in ARMB_UNITS:
            _r = _st[_lab]
            s += macro('GEshareNeg' + _mac, of_neg(_r), 'null')
            s += macro('GEdailyMed' + _mac, n3(_r['med']), 'null')
            s += macro('GEdailyLo' + _mac,  n3(_r['lo']),  'null')
            s += macro('GEdailyHi' + _mac,  n3(_r['hi']),  'null')

        # Parts de population gagnante : monde, Nord, Arctique.
        s += macro('GEpopGainDailyMed', n1(_w['pmed']), 'null')
        s += macro('GEpopGainDailyLo',  n1(_w['plo']),  'null')
        s += macro('GEpopGainDailyHi',  n1(_w['phi']),  'null')
        s += macro('GEpopGainNorthMed', n1(_n['pmed']), 'null')
        s += macro('GEpopGainNorthLo',  n1(_n['plo']),  'null')
        s += macro('GEpopGainNorthHi',  n1(_n['phi']),  'null')
        _a = _st['North>=66.5N']
        s += macro('GEpopGainArcticMed', n1(_a['pmed']), 'null')
        s += macro('GEpopGainArcticMax', n1(_a['pmax']), 'null')

        # ---- quartiles de revenu : la bande de la figure du corps ----
        # Le meme objet que ci-dessus, sur un autre masque. Les bornes emises
        # sont celles que la figure dessine, et elles sortent du meme appel :
        # une figure et une phrase qui liraient deux passages distincts
        # pourraient decrire deux grilles.
        _ts = _armb_tier_stats()
        s += "\n% ---- arm B: distribution by 1990 GDPpc quartile ----\n"
        s += ("% The \"Daily total\" bar of fig_ge_tiers shows BODY['tiers'],\n"
              "% the solve on the median curve; the segment shows the bounds\n"
              "% below, the distribution of the OUTPUTS. Two objects, and the\n"
              "% figure nests one in the other -- asserted below.\n")
        for _t in ('Q1', 'Q4'):
            _m = TIER_MACRO[_t]
            s += macro('GEtier%sDailyLo' % _m, n3(_ts[_t]['lo']), 'null')
            s += macro('GEtier%sDailyHi' % _m, n3(_ts[_t]['hi']), 'null')
        for _t in TIER_KEYS:
            s += macro('GEtierShareNeg%s' % TIER_MACRO[_t],
                       of_neg(_ts[_t]), 'null')

        ##############################################################
        # ASSERTIONS
        ##############################################################
        # 0. LA BARRE EST DANS SA BANDE, QUARTILE PAR QUARTILE.
        #    C'est le controle que la figure du corps ne peut pas faire
        #    d'elle-meme : elle dessine une barre et un segment, et un
        #    segment place a cote de sa barre se lit comme un segment.
        for _t in TIER_KEYS:
            _r, _bar = _ts[_t], BODY['tiers'][_t]['daily']
            if not (_r['min'] <= _bar <= _r['max']):
                raise SystemExit(
                    "[STOP] quartile %s: the daily total of the body "
                    "(%+.4f) falls outside the range of the %d draws "
                    "[%+.4f, %+.4f]." % (_t, _bar, _r['n'], _r['min'], _r['max']))
            if not (_r['min'] <= _r['lo'] <= _r['med'] <= _r['hi'] <= _r['max']):
                raise SystemExit(
                    "[STOP] quartile %s: min/q05/med/q95/max not ordered."
                    % _t)
            if _r['n'] != _nd:
                raise SystemExit("[STOP] quartile %s: %d draws against %d "
                                 "elsewhere." % (_t, _r['n'], _nd))
        # LA BANDE EST-ELLE NEGATIVE A CHAQUE QUARTILE ? C'est une PROPRIETE
        # du resultat, pas une condition de validite : elle est RAPPORTEE,
        # jamais imposee. Le corps porte une phrase qui la suppose ; si elle
        # tombe, c'est la phrase qui doit bouger, pas le controle qui doit
        # ceder. (Meme polarite que la garde d'ordre de la bande quantile.)
        _pos = [t for t in TIER_KEYS if _ts[t]['hi'] >= 0]
        print("[armB] quartiles : " + " | ".join(
            "%s med %+.4f [%+.4f, %+.4f] neg %.0f%%"
            % (t, _ts[t]['med'], _ts[t]['lo'], _ts[t]['hi'], _ts[t]['shareneg'])
            for t in TIER_KEYS))
        if _pos:
            print("[armB] (!) BAND NOT NEGATIVE in %s: the sentence of the body "
                  "\"The band on the daily total is negative at every quartile\""
                  " does NOT hold." % ', '.join(_pos))
        else:
            print("[armB] band negative at all four quartiles (q95 < 0 everywhere)")

        # LE CSV DE LA FIGURE. La figure NE RECALCULE RIEN : elle trace ce
        # fichier. C'est ce qui rend impossible qu'une barre et sa macro
        # divergent -- il n'y a qu'un seul chemin de calcul, celui-ci.
        _fig = pd.DataFrame([{
            'tier': _t,
            'daily':  BODY['tiers'][_t]['daily'],
            'annual': BODY['tiers'][_t]['annual'],
            'extra':  BODY['tiers'][_t]['extra'],
            'band_lo': _ts[_t]['lo'], 'band_hi': _ts[_t]['hi'],
            'band_med': _ts[_t]['med'], 'share_neg': _ts[_t]['shareneg'],
            'n_draws': _ts[_t]['n'],
        } for _t in TIER_KEYS])
        TABLES_CSV['fig_ge_tiers.csv'] = _fig.to_csv(index=False, lineterminator='\n')
        TABLES_CSV['ge_band_tiers_daily.csv'] = \
            _armb_tier_frame().to_csv(index=False, lineterminator='\n')
        # 1. Chaque unite a sa bande : une cellule vide dans la Table 2
        #    passerait inapercue a la relecture.
        for _lab, _mac, _ in ARMB_UNITS:
            _r = _st[_lab]
            if not np.isfinite([_r['med'], _r['lo'], _r['hi']]).all():
                raise SystemExit(f"[STOP] {_lab}: band not finite.")
            if not (_r['min'] <= _r['lo'] <= _r['med'] <= _r['hi'] <= _r['max']):
                raise SystemExit(
                    f"[STOP] {_lab}: min/q05/med/q95/max not ordered "
                    f"({_r['min']:+.4f} {_r['lo']:+.4f} {_r['med']:+.4f} "
                    f"{_r['hi']:+.4f} {_r['max']:+.4f}).")
            if _r['n'] != _nd:
                raise SystemExit(
                    f"[STOP] {_lab}: {_r['n']} draws against {_nd} elsewhere.")
        # 2. LA BANDE N'EST PAS LA COURBE MEDIANE. Le corps solve la courbe
        #    mediane du Monte Carlo ; la bande est la distribution des sorties.
        #    Les deux different, et la Table 2 les met cote a cote : on verifie
        #    que la valeur du corps tombe DANS la bande, sans quoi les deux
        #    colonnes raconteraient deux histoires.
        if not (_n['min'] <= BODY['north']['daily'] <= _n['max']):
            raise SystemExit(
                f"[STOP] the daily total of the body "
                f"({BODY['north']['daily']:+.4f}) falls outside the range of the "
                f"{_nd} draws [{_n['min']:+.4f}, {_n['max']:+.4f}].")
        print(f"[armB] North >=45N: median {_n['med']:+.4f} | "
              f"[{_n['lo']:+.4f}, {_n['hi']:+.4f}] | negative in "
              f"{_n['shareneg']:.1f} % | body (median curve) "
              f"{BODY['north']['daily']:+.4f}")
        print(f"[armB] ranking PE -> GE: {_inv} inversions out of {_nd - 1} "
              f"pairs, correlation {_rho:.2f}")
    else:
        s += "\n% ---- arm B: NOT EMITTED (grid missing or not covering) ----\n"
        for _nm in ('GEarmBnDraws', 'GEarmBrankInversions', 'GEarmBrankCorr',
                    'GEnorthDailyRqMed', 'GEnorthDailyRqLo', 'GEnorthDailyRqHi',
                    'GEnorthDailyRqMin', 'GEnorthDailyRqMax',
                    'GEnorthDailyShareNeg',
                    'GEpopGainDailyMed', 'GEpopGainDailyLo', 'GEpopGainDailyHi',
                    'GEpopGainNorthMed', 'GEpopGainNorthLo', 'GEpopGainNorthHi',
                    'GEpopGainArcticMed', 'GEpopGainArcticMax'):
            s += macro(_nm, '\\PROVISOIRE', 'none')
        for _lab, _mac, _ in ARMB_UNITS:
            for _p in ('GEshareNeg', 'GEdailyMed', 'GEdailyLo', 'GEdailyHi'):
                s += macro(_p + _mac, '\\PROVISOIRE', 'none')
        # Les quartiles suivent la meme polarite : les bornes sortent en
        # PROVISOIRE, et fig_ge_tiers.csv n'est PAS ecrit -- la figure ne se
        # produit donc pas, et son includegraphics arrete la compilation a
        # l'endroit exact du trou. C'est la polarite deja retenue pour
        # tab_ge_band : une figure muette se relit comme un resultat.
        for _t in ('Q1', 'Q4'):
            for _b in ('Lo', 'Hi'):
                s += macro('GEtier%sDaily%s' % (TIER_MACRO[_t], _b),
                           '\\PROVISOIRE', 'none')
        for _t in TIER_KEYS:
            s += macro('GEtierShareNeg%s' % TIER_MACRO[_t], '\\PROVISOIRE', 'none')
    return s


def tabular(caption, label, header, rows, colspec, note=None, arm=None):
    # The PATH, not the bare name. The package holds two files with similar
    # names, one of them differing only by case. A "Source:" mention names a
    # path; a bare name designates nothing, and the intended reader is the data
    # editor.
    s = "% GENERATED by 02_GE/emit_values.py -- DO NOT EDIT BY HAND\n"
    # L'ETIQUETTE DE BRAS, EN TETE DU FICHIER EMIS.
    # Elle est ce qui rend le controle 8 de tools/check_macros.py possible :
    # sans elle, « de quel bras sort cette table ? » ne se repond qu en
    # relisant l emetteur, et une table peut deriver sous un texte qui ne
    # bouge pas. Le commentaire LaTeX ne se rend pas ; il se relit.
    if arm is not None:
        s += "%% amenity = %s\n" % arm
    s += "\\begin{table}[!ht]\n\\centering\n\\small\n"
    s += f"\\caption{{{caption}}}\n\\label{{{label}}}\n"
    s += f"\\begin{{tabular}}{{{colspec}}}\n\\toprule\n"
    s += header + "\n\\midrule\n"
    # addlinespace is a spacing command, not a row: it takes no row separator.
    s += "".join(r + "\n" if r.strip() == "\\addlinespace" else r + " \\\\\n"
                 for r in rows)
    s += "\\bottomrule\n"
    if note:
        ncol = len(colspec.replace('|', ''))
        s += ("\\multicolumn{%d}{p{\\dimexpr\\linewidth-2\\tabcolsep\\relax}}"
              "{\\footnotesize\n" % ncol)
        s += "\n".join(note) + "}\n"
    s += "\\end{tabular}\n\\end{table}\n"
    return s


def emit_tab_aggregate():
    H = BODY['horizon']
    rows = [f"{lab} & {f3(H[lab]['daily'])} & {f3(H[lab]['annual'])}"
            f" & {f3(H[lab]['extra'])}" for lab in HORIZON_ROW]
    return tabular(
        "Aggregate welfare effects, population-weighted cumulative "
        "$\\Delta\\log U$, " + BODY_CLAUSE + ". Daily extra is the difference "
        "of the two totals, computed before rounding, so the last digit of the "
        "three columns need not add up.",
        "tab:ge_aggregate",
        "Horizon & Daily total & Annual total & Daily extra \\\\\n"
        "        & (daily arm vs.~baseline) & (annual arm vs.~baseline) & (daily vs.~annual) \\\\",
        rows, "lrrr", arm=BODY_TAG)


def emit_tab_tiers():
    names = {'daily': 'Daily total', 'annual': 'Annual total', 'extra': 'Daily extra'}
    rows = []
    for m in ('daily', 'annual', 'extra'):
        cells = " & ".join(f3(BODY['tiers'][q][m]) for q in ('Q1', 'Q2', 'Q3', 'Q4'))
        rows.append(f"{names[m]:20} & {cells} & {f3(BODY['world'][m])}")
    return tabular(
        "Pop-weighted cumulative $\\Delta\\log U$ by " + HORIZON_LABEL
        + ", across 1990 GDPpc quartiles, " + BODY_CLAUSE + ".",
        "tab:ge_tiers",
        "Margin & Q1 (low) & Q2 & Q3 & Q4 (high) & ALL \\\\",
        rows, "lrrrrr", arm=BODY_TAG)


# ----------------------------------------------------------------------
# LA CLAUSE DE LEGENDE DE tab_ge_north, DERIVEE DES LIGNES DU TABLEAU
# ----------------------------------------------------------------------
# Aucun nom de pays en dur. Les comptes et les exceptions sortent du meme
# data frame `cty` que les lignes, donc ils ne peuvent pas diverger d'elles.
# (COLD_ORDER_ISO est defini plus haut, avec les categories de pays.)


def _north_counts():
    """Comptes et exceptions sur les economies froides du tableau.

    Lus sur BODY, c'est-a-dire sur LE MEME BRAS que les lignes qu ils
    resument. Avant le 17/09/2026 ils sortaient de `cty` (bras interieur)
    pendant que la legende accompagnait un texte a f_a = 0."""
    rows = [(nm, BODY['cty'].loc[iso]) for iso, nm in COLD_ORDER_ISO]
    n = len(rows)
    gain   = [nm for nm, r in rows if r['annual'] > 0]
    loser  = [nm for nm, r in rows if r['annual'] <= 0]
    exneg  = [nm for nm, r in rows if r['extra'] >= 0]
    return {'n': n, 'gain': gain, 'loser': loser, 'extra_nonneg': exneg}


def _enum(names):
    """« A », « A and B », « A, B and C » -- jamais « A, B, C »."""
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " and " + names[-1]


def _north_caption_clause():
    """UNE PHRASE PAR FAIT, et chaque fait derive des lignes du tableau.

    Ce qui a ete retire le 17/09/2026, et pourquoi :

    - « The reversal deepens continuously with latitude ». Le mot `reversal`
      est banni du manuscrit : le gain projete du Nord est une propriete du
      METRIQUE d'exposition, pas un fait identifie qu'on renverserait. La
      phrase disait en outre une CONTINUITE que ce tableau-ci ne montre pas :
      il liste des pays, pas des bandes de latitude.
    - « no category threshold is implied ». C'est une denegation : elle
      repond a une taxonomie que le tableau n'affiche plus depuis qu'il est
      une liste continue. Nier une categorie absente la reintroduit.
    - les noms de pays en dur. Ils survivent a la valeur qui les justifiait.
    """
    c = _north_counts()
    n = c['n']
    out = []
    # Fait 1 : la marge daily-extra.
    if not c['extra_nonneg']:
        out.append("The daily-extra margin is negative for all %d economies "
                   "shown here." % n)
    else:
        out.append("The daily-extra margin is negative for %d of the %d "
                   "economies shown here; it is not for %s."
                   % (n - len(c['extra_nonneg']), n, _enum(c['extra_nonneg'])))
    # Fait 2 : le gain a la moyenne annuelle, COMPTE puis exception nommee.
    ng = len(c['gain'])
    if not c['loser']:
        out.append("All %d gain under the annual mean." % n)
    else:
        out.append("%d of the %d gain under the annual mean; %s %s not."
                   % (ng, n, _enum(c['loser']),
                      'does' if len(c['loser']) == 1 else 'do'))
    # Fait 3 : ce que porte la derniere ligne.
    out.append("The bottom row aggregates all model cells above "
               "$" + BS + "geq 45^" + BS + "circ$N.")
    # Fait 4 : ou lire le gradient de latitude, sans affirmer sa forme.
    out.append("Table~" + BS + "ref{tab:ge_latitude} reports the same "
               "margins by latitude band.")
    return " ".join(out)


def emit_tab_north():
    """Liste continue des economies froides, dans l'ORDRE DU PANEL C DE LA TABLE 1
    (tab:ge) : REVERSAL+AMPLIFICATION+FRANCE. Plus de categorie de signe
    (Reversal/Amplification) : sous le sc3 corrige, CAN/DEU gagnent sous l'annuel
    et la taxonomie devient caduque (§6). LOT 74 (E5) : cet ordre N'EST PAS
    decroissant en gain annuel (01/10 : FIN +0,029 devant RUS +0,026 en gain, mais
    apres elle dans la liste) ; la docstring et la legende le disaient. On ne trie
    pas -- l'ordre commun aux tables vaut mieux -- on dit lequel c'est. Le
    gradient continu robuste est porte par tab_ge_latitude."""
    COLD_ORDER = REVERSAL + AMPLIFICATION + FRANCE
    rows = []
    for iso, name in COLD_ORDER:
        r = BODY['cty'].loc[iso]
        rows.append(f"{name:10} & {f3(r['daily'])} & {f3(r['annual'])} & {f3(r['extra'])}")
    rows.append("\\addlinespace")
    # LA PHRASE EST DERIVEE DES LIGNES, PAS ECRITE A COTE D'ELLES.
    # L'ancienne legende nommait France en dur sur deux enonces, tous deux faux
    # sur les valeurs courantes : la marge daily-extra est negative pour les SEPT
    # economies (France comprise, -0,063), et les perdants a la moyenne annuelle
    # sont DEUX (Canada -0,003 et France). Un nom de pays en dur dans une legende
    # survit a la valeur qui l'avait justifie.
    # bold in math mode: \mathbf, because \textbf does not cross $...$
    def b3(v):
        return f"$\\mathbf{{{s3(v)}}}$"
    _nb = BODY['north']
    rows.append(f"\\textbf{{North $\\geq 45^\\circ$N}} & {b3(_nb['daily'])}"
                f" & {b3(_nb['annual'])} & {b3(_nb['extra'])}")
    return tabular(
        "Cumulative $\\Delta\\log U$ by " + HORIZON_LABEL
        + " for cold-climate economies, in the order of Table~\\ref{tab:ge}, "
        + BODY_CLAUSE + ". " + _north_caption_clause(),
        "tab:ge_north",
        "Country & Daily total & Annual total & Daily extra \\\\",
        rows, "lrrr")


def emit_tab_latitude():
    rows = [f"$\\geq {t:g}^\\circ$N & {p:,.0f} & {f3(a)} & {f3(e)} & {f3(d)}"
            for t, p, a, e, d in BODY['lat_rows']]
    return tabular(
        "Pop-weighted cumulative $\\Delta\\log U$ by " + HORIZON_LABEL
        + " across latitude bands, "
        + BODY_CLAUSE + ".",
        "tab:ge_latitude",
        "Threshold & Pop. (M) & Annual total & Daily extra & Daily total \\\\",
        rows, "lrrrr", arm=BODY_TAG)


# ---- amenity tables (lower bound g=0 against the interior g=q50) ----
# These exist only when the null runs are present in the headline arm.
def emit_tab_sensitivity():
    """Agregat MONDE, pop-pondere, sous les deux specs d'amenite.
    (!) La marge Daily extra CHANGE DE SIGNE (base <0, null >0) : l'objet
    robuste est le Nord, pas le monde. Reporte tel quel, honnetement."""
    rows = [
        f"Interior estimate $\\hat f_a^{{\\,q50}}$ & {f3(world['daily'])}"
        f" & {f3(world['annual'])} & {f3(world['extra'])}",
        f"Lower bound $f_a = 0$ & {f3(NULLagg['world']['daily'])}"
        f" & {f3(NULLagg['world']['annual'])} & {f3(NULLagg['world']['extra'])}",
    ]
    return tabular(
        "Sensitivity of global pop-weighted cumulative $\\Delta\\log U$ by "
        + HORIZON_LABEL + " to the amenity specification",
        "tab:ge_sensitivity",
        "Specification & Daily total & Annual total & Daily extra \\\\",
        rows, "lrrr", arm=HEAD_TAG + '+' + BODY_TAG)


def emit_tab_north_sens():
    """Pays froids, groupes par categorie, sous les deux specs. Colonnes :
    Annual total {q50 | LB} et Daily extra {q50 | LB}."""
    def b3(v):
        return f"$\\mathbf{{{s3(v)}}}$"
# Continuous list ordered by decreasing annual gain, with no sign category.
    COLD_ORDER = REVERSAL + AMPLIFICATION + FRANCE
    rows = []
    for iso, name in COLD_ORDER:
        rb = cty[cty['GID_0'] == iso].iloc[0]
        rz = NULLagg['cty'].loc[iso]
        rows.append(f"{name:10} & {f3(rb['annual'])} & {f3(rz['annual'])}"
                    f" & {f3(rb['extra'])} & {f3(rz['extra'])}")
    rows.append("\\addlinespace")
    nz = NULLagg['north']
    rows.append(f"\\textbf{{North $\\geq 45^\\circ$N}} & {b3(north['annual'])}"
                f" & {b3(nz['annual'])} & {b3(north['extra'])} & {b3(nz['extra'])}")
    return tabular(
        "Cumulative $\\Delta\\log U$ by " + HORIZON_LABEL
        + " for cold-climate economies under the "
        "two amenity specifications (interior $\\hat f_a^{\\,q50}$ vs.\\ lower "
        # 22/09/2026 : la legende annoncait « ordered by annual-mean gain ».
        # Elle ne l'est pas : l'ordre affiche est celui du Tableau 1,
        # panel C, ou la Finlande devance la Russie et la Suede la Norvege.
        # On NE TRIE PAS -- l'ordre commun aux deux tableaux vaut mieux
        # qu'un tri qui les desaligne -- on dit lequel c'est.
        # ⚠️ ANTISLASH DOUBLE, ET C'EST OBLIGATOIRE. Ecrit `"...Table~\ref{...}"`
        # dans une chaine NORMALE, Python lit `\r` comme un RETOUR CHARIOT : la
        # legende rendue portait « Table~ef{tab:ge} » et LaTeX imprimait
        # « ef{tab:ge} ». Mesure du 22/09, attrapee par la relecture du rendu.
        "bound $f_a = 0$), in the order of Table~\\ref{tab:ge}.",
        "tab:ge_north_sens",
        "        & \\multicolumn{2}{c}{Annual total} & \\multicolumn{2}{c}{Daily extra} \\\\\n"
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n"
        "Country & $\\hat f_a^{\\,q50}$ & $f_a = 0$ "
        "& $\\hat f_a^{\\,q50}$ & $f_a = 0$ \\\\",
        rows, "lrrrr", arm=HEAD_TAG + '+' + BODY_TAG)


def emit_tab_ge():
    """Table 2 du corps : monde, latitudes, pays, sous f_a = 0.

    Lit UNIQUEMENT NULLagg, c'est-a-dire aggregates_for(NULL_TAG). Aucune
    agregation n'est redefinie ici."""
    A = NULLagg
    ##################################################################
    # LA BANDE EST A COTE, PAS DANS, LA COLONNE « Daily total ».
    #
    # « Daily total » reste le modele solve sur la COURBE MEDIANE du Monte
    # Carlo -- le run publie. Les deux colonnes ajoutees decrivent la
    # DISTRIBUTION DES SORTIES sur les tirages de l'arm B. Les deux objets ne
    # coincident pas (mesure : -0,0431 contre une mediane de tirages a
    # -0,0451) et les confondre dans une seule cellule rendrait le tableau
    # illisible : une mediane de tirages n'est pas le tirage median, et aucun
    # des deux n'est la courbe mediane.
    # Les cles d'unite du cadre arm B ne sont pas celles de NULLagg : la
    # correspondance est explicite, et une unite sans bande est une ERREUR,
    # pas une cellule vide.
    ##################################################################
    _BAND_KEY = {
        'World': 'World', r'$\geq 45^\circ$N': 'North>=45N',
        r'$\geq 55^\circ$N': 'North>=55N',
        r'$\geq 66.5^\circ$N': 'North>=66.5N',
        'Russia': 'Russia', 'Finland': 'Finland', 'Norway': 'Norway',
        'Sweden': 'Sweden', 'Canada': 'Canada', 'Germany': 'Germany',
        'France': 'France',
    }
    _BST = _armb_stats() if ARMB_OK else None
    if ARMB_OK:
        assert_armb_csv(_armb_frame())

    def _band_row(label):
        """La SECONDE ligne de l'unite : la bande, en retrait, sous son total.

        Elle ne partage pas les colonnes des marges : la bande decrit le total
        journalier seul, et l'aligner sous « Annual total » ou « Daily extra »
        suggererait une bande sur ces marges-la, qui n'existe pas."""
        if _BST is None:
            return None
        k = _BAND_KEY.get(label)
        if k is None or k not in _BST:
            raise SystemExit(
                "[STOP] tab_ge: no band for row %r. An empty cell "
                "would pass a proofread; its absence must stop the run." % label)
        r = _BST[k]
        return (r"\quad{\footnotesize 5th--95th pct.} & "
                r"{\footnotesize [%s,\,%s]} & "
                r"\multicolumn{2}{l}{\footnotesize negative in %s draws}"
                % (fb(r['lo']), fb(r['hi']), of_neg(r)))

    def _line(label, d):
        ecart = abs(d["daily"] - (d["annual"] + d["extra"]))
        if ecart > 1e-9:
            raise SystemExit(
                "[STOP] tab_ge, row %s: Daily total (%.6f) != Annual (%.6f) "
                "+ Extra (%.6f), gap %.2e. This is an IDENTITY of logarithms; "
                "it fails only if the three columns do not come from the same "
                "triplet of runs."
                % (label, d["daily"], d["annual"], d["extra"], ecart))
        return "%s & %s & %s & %s" % (label, f3(d["daily"]), f3(d["annual"]),
                                      f3(d["extra"]))

    def _lines(label, d):
        """La ligne de marges, puis sa ligne de bande quand elle existe."""
        out = [_line(label, d)]
        b = _band_row(label)
        if b is not None:
            out.append(b)
        return out

    def _panel(title):
        return r"\multicolumn{4}{l}{\textit{" + title + "}}"

    rows = [_panel("Panel A. World")] + _lines("World", A["world"]) + [
            r"\addlinespace", _panel("Panel B. By latitude")] \
        + _lines(r"$\geq 45^\circ$N", A["north"])
    for tag, thr in sorted(BAND_LAT.items(), key=lambda kv: kv[1]):
        rows += _lines(r"$\geq %g^\circ$N" % thr, A["lat"][tag])
    rows += [r"\addlinespace", _panel("Panel C. Cold-climate economies")]
    # FRANCE COMPRISE. Le corps ecrit « GEnorthNGain of GEnorthNcty gain
    # under the annual mean », soit 6 de 7, et l'exhibit en regard n'en montrait
    # que SIX : le lecteur comptait six lignes, toutes gagnantes, et lisait
    # « 6 de 7 » a cote. Le septieme pays existait, dans la Table du supplement.
    # Une ligne de plus ici, et le compte se verifie sur l'exhibit lui-meme.
    for iso, name in REVERSAL + AMPLIFICATION + FRANCE:
        if iso in A["cty"].index:
            r = A["cty"].loc[iso]
            rows += _lines(name, {m: float(r[m]) for m in MARGINS})
    # Part de la population MONDIALE dans des cellules qui GAGNENT : A["gain"]
    # sort de la meme passe que les marges, sur la meme grille et les memes poids.
    rows += [r"\addlinespace",
             "Population in gaining cells (\\%%) & %s & %s & ---"
             % (n1(A["gain"]["daily"]), n1(A["gain"]["annual"]))]

    ## LEGENDE REECRITE 17/09/2026 (lot 8). L'etiquette d'horizon passe par la
    ## MACRO et non plus par la constante : la table devient une fonction de
    ## values.tex seule, et check_horizon n'a rien de nouveau a surveiller ici.
    return tabular(
        "General-equilibrium welfare effects by \\GEhorizonLabel{}, cumulative "
        "population-weighted $\\Delta\\log U$, amenities at the lower bound "
        "$f_a=0$. Daily total and Annual total are each arm against the "
        "no-climate baseline; Daily extra is their difference, computed before "
        "rounding. Second line under each unit: 5th--95th percentile range and "
        "share of negative outcomes of the daily total across \\GEarmBnDraws{} "
        "Monte Carlo draws of the daily dose-response, each solved in general "
        "equilibrium (Table~\\ref{tab:ge_band_ranked}).",
        "tab:ge",
        "& Daily total & Annual total & Daily extra \\\\",
        rows, "lrrr", arm=BODY_TAG)


BAND_LABEL = {
    'world':  'World',
    'north':  '$\\geq 45^\\circ$N',
    'latFF':  '$\\geq 55^\\circ$N',
    'latA':   '$\\geq 66.5^\\circ$N',
    'ctyRus': 'Russia',
    'ctyFin': 'Finland',
    'ctyCan': 'Canada',
    'ctyDeu': 'Germany',
}


def _band_cross_check():
    """La ligne du Nord contre le producteur versionne 02e.

    02e lit les welfare, refait son propre masque et sa propre ponderation.
    C'est un chemin entierement distinct d'aggregates_for(). L'egalite est donc
    EXIGIBLE : deux chemins qui calculent le meme objet doivent rendre le meme
    nombre, sinon l'un des deux decrit autre chose."""
    p = os.path.join(OUTPUT_DIR, 'ge_north_contrast_band.csv')
    if not os.path.exists(p):
        print("[bande] ge_north_contrast_band.csv missing: cross-check NOT "
              "RUN. This is a gap, not a success.")
        return
    d = pd.read_csv(p).set_index('quantity')['value_raw']
    worst, lignes = 0.0, []
    for qq in BAND_Q:
        k = 'corrected_' + qq
        if k not in d.index:
            continue
        ref = float(d.loc[k])
        got = _band_val(BAND[qq], 'north', 'extra')
        worst = max(worst, abs(ref - got))
        lignes.append((qq, ref, got, abs(ref - got)))
    print("[bande] Northern cross-check against 02e: max gap %.2e over %d quantiles"
          % (worst, len(lignes)))
    if worst > 1e-6:
        for qq, ref, got, dd in lignes:
            print("   %s  02e=%+.6f  emitter=%+.6f  gap=%.2e" % (qq, ref, got, dd))
        sys.exit("STOP: the two paths do not give the same Northern contrast. "
                 "One of them does not describe the object the table announces.")


def emit_tab_ge_band():
    """Les huit unites, marge daily extra, aux trois quantiles."""
    _band_cross_check()
    rows = []
    for unit, _mac in BAND_UNITS:
        v = [_band_val(BAND[qq], unit, BAND_MARGIN) for qq in BAND_Q]
        rows.append("%s & %s & %s & %s"
                    % (BAND_LABEL[unit], fb(v[0]), f3(v[1]), fb(v[2])))
    amen = (BAND_AMENITY or 'UNKNOWN')
    amen_txt = ("the interior estimate $\\hat f_a^{\\,q50}$" if amen == 'baseline'
                else "the lower bound $f_a = 0$" if amen in ('', 'nan')
                else "amenity specification \\texttt{%s}" % amen)
    note = [
        "\\textit{Note.} Columns evaluate the same general equilibrium on the q05, q50",
        "and q95 columns of the Monte Carlo dose-response tables, the same",
        "percentile for the productivity and the amenity arm. Those columns are",
        "POINTWISE quantiles across draws, taken independently at each",
        "temperature: no single draw traces them, and the four Monte Carlo runs",
        "behind them use four distinct seeds. What the three solves share is a",
        "level of assumed severity, not a draw index, so no draw-level",
        "covariance between the curves is carried. A functional of a pointwise",
        "envelope is not a quantile of that functional's distribution: this is",
        "a range of assumed severity, NOT a confidence interval, and the median",
        "is not guaranteed to lie between the bounds -- for the world aggregate",
        "it does not. The three runs use " + amen_txt + ".",
        "The daily extra is the daily arm against the annual arm, so both",
        "numerator and",
        "denominator are read at the same percentile and the severity assumption",
        "differences out of the comparison between metrics.",
    ]
    return tabular(
        "Monte Carlo band on the daily-extra contrast by " + HORIZON_LABEL
        + ", cumulative "
        "pop-weighted $\\Delta\\log U$, by unit.",
        "tab:ge_band",
        "Unit & $q_{05}$ & $q_{50}$ & $q_{95}$ \\\\",
        rows, "lrrr", note, arm=HEAD_TAG)


##################################################################################
# LA BANDE SUR TIRAGES CLASSES (rq05 / rq50 / rq95)
#
# CE QU'ELLE EST, ET CE QU'ELLE N'EST PAS.
# Les six runs `_q05/_q50/_q95` deja au manifeste evaluent le GE sur les
# COLONNES de quantiles ponctuels des tables de dose-reponse : ces colonnes ne
# sont tracees par AUCUN tirage, et une fonctionnelle d'une enveloppe ponctuelle
# n'est pas un quantile de la distribution de cette fonctionnelle. Les runs
# `_rq*` evaluent le GE sur de VRAIS TIRAGES, classes par le changement nord
# projete en equilibre partiel (04_contrast/ct5_rank_mc_draws.R).
#
# CHAQUE BRAS EST TIRE DANS SA PROPRE SPECIFICATION -- journalier spline4,
# annuel spline3, ce que les poids d'Akaike de 01a designent. Il n'existe donc
# PAS de tirage commun aux deux bras, et la marge « extra » ne se calcule pas
# comme sc2/sc3 a un tirage donne : ce serait apparier deux tirages
# independants et publier leur difference comme si elle etait tiree. Elle est
# ENCADREE, en traitant les bras comme independants :
#      min = daily q05 - annual q95      max = daily q95 - annual q05
# C'est le bon sens de l'erreur : la correlation mesuree entre bras est
# POSITIVE, donc la vraie dispersion de la difference est PLUS PETITE.
##################################################################################
RANKED_Q = ('rq05', 'rq50', 'rq95')
RANKED_FILES = [f'welfare_sc{sc}_{q}{NULL_TAG}.csv'
                for sc in (2, 3) for q in RANKED_Q]
RANKED_OK = all(os.path.exists(os.path.join(RUNS_DIR, f)) for f in RANKED_FILES)

##########################################################################
# LE MILLESIME DE CE BLOC, ET SON PERIMETRE. Ajoute le 21/09/2026.
#
# Ce bloc ne testait que l'EXISTENCE des six fichiers, la ou le bras de tete
# (l.479), le bras null (l.861) et la bande ponctuelle (l.1064) passent tous
# par assert_one_arm. Mesure du 21/09 : les six rq* dataient du 16/09, donc
# d'AVANT la division par theta, pendant que le corps venait d'etre refait.
# L'emetteur rendait alors
#     GEnorthNet        classe -0.0402   corps -0.0067
#     GEnorthDailyExtra classe -0.1030   corps -0.0164
# et \GEnorthRankedCentreGap passait de 0,006 a 0,087 -- un ecart qui n'est
# pas une mesure de la distance entre une enveloppe ponctuelle et un tirage
# median, mais le facteur 6,5 de la division. Publiable sans un avertissement.
#
# LA LISTE PORTE NEUF FICHIERS, PAS SIX, et c'est le point. Les six rq* entre
# eux PARTAGENT une empreinte perimee : un assert_one_arm pose sur eux seuls
# serait passe au VERT le 21/09 au matin. Ce que ce bloc met cote a cote, ce
# sont les tirages classes ET le corps -- _ranked_frame() lit sc1 null comme
# baseline, et \GEnorthRankedCentreGap compare le centre classe a BODY, qui
# sort de sc2/sc3 null. Le perimetre d'une assertion est l'ensemble des runs
# qu'une meme phrase rapproche, pas l'ensemble des runs qu'une meme boucle
# ouvre.
#
# Polarite : ARRET, comme les trois autres. Deux conventions mises en regard
# ne font pas une comparaison imprecise, elles font une comparaison d'autre
# chose.
##########################################################################
RANKED_ASSERT_FILES = RANKED_FILES + [
    f'welfare_sc{sc}_{QUANTILE}{NULL_TAG}.csv' for sc in (1, 2, 3)]
if RANKED_OK:
    RANKED_ARM = assert_one_arm(RANKED_ASSERT_FILES, 'ranked band + body')
    print("[bande] six ranked-draw runs present | vintage shared "
          "with the body: %s"
          % (RANKED_ARM['inputs_md5'][:8] if RANKED_ARM else 'NOT CHECKED'))
else:
    _miss = [f for f in RANKED_FILES
             if not os.path.exists(os.path.join(RUNS_DIR, f))]
    print(f"[bande] ranked draws MISSING ({len(_miss)}/6) -> "
          f"table and macros not emitted")

# Les unites de la table, dans l'ordre. Le selecteur est une FONCTION du cadre
# pixel : aucun nom de pays en dur dans la boucle d'emission.
RANKED_UNITS = [
    ('World',                      lambda d: d),
    (r'North $\geq 45^\circ$N',    lambda d: d[d['y'] >= 45.0]),
    (r'North $\geq 55^\circ$N',    lambda d: d[d['y'] >= 55.0]),
    (r'North $\geq 66.5^\circ$N',  lambda d: d[d['y'] >= 66.5]),
    ('Russia',                     lambda d: _pays(d, 'RUS')),
    ('Finland',                    lambda d: _pays(d, 'FIN')),
    ('Canada',                     lambda d: _pays(d, 'CAN')),
    ('Germany',                    lambda d: _pays(d, 'DEU')),
]


def _load_w(sc, q, tag):
    f = os.path.join(RUNS_DIR, f'welfare_sc{sc}_{q}{tag}.csv')
    d = pd.read_csv(f)
    d = d.loc[:, ~d.columns.str.contains('^Unnamed')]
    if d.values.shape != (T, N):
        raise SystemExit(f"[STOP] {os.path.basename(f)}: shape "
                         f"{d.values.shape}, ({T}, {N}) expected.")
    return d.values


def _ranked_frame():
    """Un cadre pixel portant, par percentile, le total journalier (sc2 contre
    la baseline) et le total annuel (sc3 contre la MEME baseline).

    La baseline sc1 n'est PAS rejouee : sans rechauffement les deux fonctions
    de dose-reponse sont evaluees a exposition inchangee, donc elle ne depend
    ni du tirage ni de la specification d'amenite -- welfare_sc1_q50.csv et
    welfare_sc1_q50_NULLamen.csv ont le meme md5. Elle est solvee une fois et
    partagee ; le manifeste porte le meme baseline_run_id partout."""
    W1 = _load_w(1, QUANTILE, NULL_TAG)
    q = pix[['id', 'y', 'x', 'pop']].copy()
    for p_ in RANKED_Q:
        q['daily_' + p_] = np.cumsum(
            np.log(_load_w(2, p_, NULL_TAG) / W1), axis=0)[TERMINAL_ROW, :]
        q['annual_' + p_] = np.cumsum(
            np.log(_load_w(3, p_, NULL_TAG) / W1), axis=0)[TERMINAL_ROW, :]
    return q


def _ranked_rows():
    q = _ranked_frame()
    rows = []
    for label, sel in RANKED_UNITS:
        d = sel(q)
        r = {'label': label}
        for p_ in RANKED_Q:
            r['daily_' + p_] = wmean(d, 'daily_' + p_)
            r['annual_' + p_] = wmean(d, 'annual_' + p_)
        # ⚠️ L'ENCADREMENT SE PREND SUR LES EXTREMES OBSERVES, PAS SUR q05/q95.
        # Une premiere version posait min = daily_q05 - annual_q95 et
        # max = daily_q95 - annual_q05, ce qui suppose que chaque marge CROIT
        # avec le percentile. C'est vrai du Nord -- le classement porte sur lui
        # -- et FAUX ailleurs : sur la ligne World, daily valait
        # (-0,096 ; -0,086 ; -0,103), donc « min » sortait a -0,033 et « max »
        # a -0,073, l'un plus grand que l'autre. Les deux etiquettes mentaient.
        _d = [r['daily_' + p_] for p_ in RANKED_Q]
        _a = [r['annual_' + p_] for p_ in RANKED_Q]
        r['extra_min'] = min(_d) - max(_a)
        r['extra_q50'] = r['daily_rq50'] - r['annual_rq50']
        r['extra_max'] = max(_d) - min(_a)
        rows.append(r)
    gain = {p_: 100.0 - 100.0 * q.loc[q['daily_' + p_] < 0, 'pop'].sum()
            / q['pop'].sum() for p_ in RANKED_Q}
    return rows, gain


def assert_armb_csv(frame):
    """Le CSV de bande decrit-il la grille que le manuscrit annonce ?

    lignes = tirages x unites. Un CSV plus court est une grille PARTIELLE, et
    une grille partielle de tirages ordonnes par percentile n est pas un
    echantillon : c est le bas de la distribution. Un CSV plus long est un
    melange de deux grilles. Les deux se publient sans bruit -- celui qui est
    reste sur disque du 16 au 17/09 en porte trente quand le texte en dit
    soixante.

    Polarite : ARRET. Ce n est pas une propriete du resultat, c est la
    coherence entre un fichier livre et le nombre imprime a cote de lui."""
    att = len(ARMB_COLS) * len(ARMB_UNITS)
    if len(frame) != att:
        raise SystemExit(
            "[STOP] ge_band_ranked_daily.csv: %d row(s), %d expected\n"
            "  (%d arm B draws x %d units).\n"
            "  The table and \\GEarmBnDraws announce %d draws; a CSV that\n"
            "  does not carry as many describes another grid than the text."
            % (len(frame), att, len(ARMB_COLS), len(ARMB_UNITS),
               len(ARMB_COLS)))
    if frame['draw_id'].nunique() != len(ARMB_COLS):
        raise SystemExit(
            "[STOP] ge_band_ranked_daily.csv: %d distinct draw(s) for "
            "%d column(s) read.\n  Two columns carry the same draw_id, or "
            "a column is missing." % (frame['draw_id'].nunique(), len(ARMB_COLS)))
    print("[armB] band CSV: %d rows = %d draws x %d units (asserted)"
          % (len(frame), len(ARMB_COLS), len(ARMB_UNITS)))
    return frame


TIER_KEYS  = ('Q1', 'Q2', 'Q3', 'Q4')
TIER_MACRO = {'Q1': 'Qone', 'Q2': 'Qtwo', 'Q3': 'Qthree', 'Q4': 'Qfour'}
_ARMB_TIER_CACHE = {}


def _armb_tier_frame():
    """Une ligne par tirage et par QUARTILE de PIB par tete 1990 : le total
    journalier, agrege sur les memes cellules et avec les memes poids que la
    barre « Daily total » de la figure des quartiles.

    L'AGREGATION SE FAIT SUR LES SORTIES PAR CELLULE, PAS SUR UN RUN NEUF.
    Chaque colonne d'arm B est deja un solve complet du modele ; il n'en manque
    aucun et on n'en lance aucun. Ce qui change par rapport a _armb_frame(),
    c'est le seul masque : le quartile a la place de l'unite geographique.

    L'APPARTENANCE AU QUARTILE NE DEPEND PAS DU TIRAGE. Elle depend du PIB par
    tete de 1990, qui est un intrant. On reprend donc `tier_map`, l'affectation
    du module, pour que la bande et la barre decrivent les MEMES groupes de
    pays -- une re-partition par tirage donnerait une bande sur une frontiere
    mouvante, ce qui n'est pas la question posee.

    MEMOISE, pour la meme raison que _armb_frame : chaque construction relit
    soixante matrices de welfare."""
    if 'f' in _ARMB_TIER_CACHE:
        return _ARMB_TIER_CACHE['f']
    side = pd.read_csv(os.path.join(RUNS_DIR, 'ge_runs_band_armB.csv'))
    side = side.set_index('column')
    tier_map = dict(zip(cty['GID_0'], cty['tier']))
    W1 = _load_w(1, QUANTILE, NULL_TAG)
    out = []
    for c in ARMB_COLS:
        cum = np.cumsum(np.log(_load_w(2, c, NULL_TAG) / W1),
                        axis=0)[TERMINAL_ROW, :]
        q = pix[['id', 'y', 'x', 'pop']].copy()
        q['daily'] = cum
        q = par_pays(q)
        q['tier'] = q['GID_0'].map(tier_map)
        for t in TIER_KEYS:
            d = q[q['tier'] == t]
            out.append({
                'draw_id': int(side.loc[c, 'draw_id']),
                'draw_rank': int(side.loc[c, 'draw_rank']),
                'column': c, 'percentile': float(side.loc[c, 'percentile']),
                'tier': t,
                'daily_total': wmean(d, 'daily'),
                'n_countries': int(d['GID_0'].nunique()),
            })
    _ARMB_TIER_CACHE['f'] = pd.DataFrame(out)
    return _ARMB_TIER_CACHE['f']


def _armb_tier_stats():
    """Par quartile : la distribution du total journalier sur les sorties GE."""
    ab = _armb_tier_frame()
    out = {}
    for t in TIER_KEYS:
        v = ab.loc[ab['tier'] == t, 'daily_total'].values
        out[t] = {'n': len(v),
                  'med': float(np.median(v)),
                  'lo': float(np.quantile(v, 0.05)),
                  'hi': float(np.quantile(v, 0.95)),
                  'min': float(v.min()), 'max': float(v.max()),
                  'shareneg': 100.0 * float((v < 0).mean()),
                  'nneg': int((v < 0).sum())}
    return out


def frac_neg(r):
    """PARTS DE TIRAGES EN FRACTIONS (lot 34, 25/09/2026) : « 59/60 » et non
    « 98.3 % ». Exact, auto-documente, et « 60/60 » dit « tous » sans suggerer
    une precision au dixieme sur soixante tirages. Le compte vient du meme
    vecteur que la part ; on verifie qu'ils s'accordent."""
    k, n = r['nneg'], r['n']
    if abs(100.0 * k / n - r['shareneg']) > 1e-9:
        raise SystemExit('[STOP] part negative %.6f != %d/%d' % (r['shareneg'], k, n))
    return '%d/%d' % (k, n)


def of_neg(r):
    """LA MEME FRACTION, EN PROSE (lot 35, 26/09/2026) : « 59 of 60 ». Une macro
    lue dans une phrase ne peut pas porter « 59/60 » : « negative in 59/60 of
    draws » n'est pas de l'anglais. Les tables gardent « 59/60 » (frac_neg) sous
    un en-tete qui le dit ; les phrases appellent la macro HORS mode math, sinon
    « of » s'imprimerait en italique et sans espaces."""
    frac_neg(r)                     # meme garde d'accord part / compte
    return '%d of %d' % (r['nneg'], r['n'])


def _armb_stats():
    """Par unite, la distribution du total journalier et de la part gagnante
    sur les sorties GE. Une seule lecture, partagee par la table, par les
    colonnes de la Table 2 et par les macros : les trois ne peuvent pas
    decrire trois etats differents de la grille."""
    ab = _armb_frame()
    out = {}
    for lab, mac, _ in ARMB_UNITS:
        g = ab[ab['unit'] == lab]
        v = g['daily_total'].values
        p = g['pop_gain_share'].values
        out[lab] = {
            'macro': mac, 'n': len(v),
            'med': float(np.median(v)),
            'lo': float(np.quantile(v, 0.05)),
            'hi': float(np.quantile(v, 0.95)),
            'min': float(v.min()), 'max': float(v.max()),
            'shareneg': 100.0 * float((v < 0).mean()),
            'nneg': int((v < 0).sum()),
            'pmed': float(np.median(p)),
            'plo': float(np.quantile(p, 0.05)),
            'phi': float(np.quantile(p, 0.95)),
            'pmax': float(p.max()),
        }
    return out


def _avg_ranks(v):
    """Rangs moyens, ex aequo compris. numpy n a pas de rankdata et scipy n est
    pas une dependance du paquet ; six lignes valent mieux qu un import de
    plus dans une chaine qui doit tourner sur la machine du data editor."""
    order = np.argsort(v, kind='mergesort')          # stable
    r = np.empty(len(v), dtype=float)
    r[order] = np.arange(1, len(v) + 1, dtype=float)
    # Moyenne des rangs sur chaque groupe d ex aequo.
    uniq, inv_idx, counts = np.unique(v, return_inverse=True, return_counts=True)
    if (counts > 1).any():
        sums = np.zeros(len(uniq))
        np.add.at(sums, inv_idx, r)
        r = (sums / counts)[inv_idx]
    return r


def _armb_rank_corr():
    """CORRELATION DE RANGS (Spearman) entre le classement du tirage en
    equilibre PARTIEL et son classement en equilibre GENERAL, Nord >= 45N.

    CE QUI A CHANGE LE 17/09/2026, ET POURQUOI.
    E.7 ecrit « the rank correlation between the two is ... ». La fonction
    rendait un PEARSON entre le rang PE (1..n) et le NIVEAU GE -- une
    correlation entre un rang et une valeur, qui n est pas une correlation de
    rangs et qui n est pas invariante a une transformation monotone du
    niveau. Le mot de la phrase et la statistique ne designaient pas le meme
    objet. C est desormais Spearman : rang contre rang, ce que la phrase dit.

    `inv` ne change pas : il compte les inversions de l ordre PE sur les
    paires ADJACENTES, et il se lit sur les niveaux comme sur les rangs."""
    ab = _armb_frame()
    n = ab[ab['unit'] == 'North>=45N'].sort_values('draw_rank')
    v = n['daily_total'].values
    inv = int(sum(1 for i in range(len(v) - 1) if v[i] > v[i + 1]))
    rho = float(np.corrcoef(np.arange(1, len(v) + 1, dtype=float),
                            _avg_ranks(v))[0, 1])
    return inv, rho, len(v)


def emit_tab_ge_band_ranked():
    """LA DISTRIBUTION DES SORTIES, ET ELLE SEULE.

    (Assertion de grille en tete : la table et le CSV livre decrivent la meme
    chose, ou rien n est emis.)

    Les trois tirages selectionnes (rq05/rq50/rq95) ne figurent plus ici. Ils
    ont servi a construire la bande, et le classement en equilibre partiel qui
    les designe ne se transmet que faiblement au GE : les afficher en colonnes
    « q05/q50/q95 » a cote d'une vraie distribution invitait a les lire comme
    des percentiles de sortie, ce qu'ils ne sont pas. Ils restent decrits en
    S5.6b comme des tirages solves.
    """
    assert_armb_csv(_armb_frame())
    st = _armb_stats()
    L = ["%% GENERATED by 02_GE/emit_values.py -- DO NOT EDIT BY HAND.",
         # Les 60 tirages de l'arm B sont solves a f_a = 0 : cette table sort
         # du bras du corps, et le controle 8 doit pouvoir le lire.
         "%% amenity = " + BODY_TAG,
         r"\begin{table}[!ht]\centering\footnotesize",
         r"\setlength{\tabcolsep}{3pt}",
         r"\caption{Distribution of the general-equilibrium daily total by "
         + HORIZON_LABEL + r" across Monte Carlo draws of the daily "
         r"dose-response, cumulative population-weighted $\Delta\log U$.}",
         r"\label{tab:ge_band_ranked}",
         r"\begin{tabular}{lcccccccc}", r"\toprule",
         r"& \multicolumn{6}{c}{Daily total} "
         r"& \multicolumn{2}{c}{Population in gaining cells (\%)} \\",
         r"\cmidrule(lr){2-7}\cmidrule(lr){8-9}",
         # « % neg. » jusqu'au lot 33 ; la cellule porte une FRACTION de
         # tirages depuis le lot 34 (frac_neg), l'en-tete le dit (lot 35).
         r"Unit & median & q05 & q95 & min & max & neg.\ draws "
         r"& median & q95 \\", r"\midrule"]
    for lab, _mac, _ in ARMB_UNITS:
        r = st[lab]
        L.append(r"%s & %s & %s & %s & %s & %s & %s & %.1f & %.1f \\" % (
            lab.replace('>=45N', r' $\geq 45^\circ$N')
               .replace('>=55N', r' $\geq 55^\circ$N')
               .replace('>=66.5N', r' $\geq 66.5^\circ$N'),
            f3(r['med']), fb(r['lo']), fb(r['hi']),
            fb(r['min']), fb(r['max']),
            frac_neg(r), r['pmed'], r['phi']))
    L += [r"\bottomrule",
          r"\multicolumn{9}{p{\dimexpr\linewidth-2\tabcolsep\relax}}{\footnotesize",
          r"\textit{Note.} Distribution of the daily total across "
          r"$\GEarmBnDraws$ Monte Carlo draws of the daily dose-response",
          r"(specification \EMPfigDrfSpecDaily), each solved in general",
          r"equilibrium with amenities at the lower bound and the no-climate",
          r"baseline common to all runs. Draws are a stratified sample of the",
          r"Monte Carlo distribution, taken at evenly spaced percentiles of the",
          r"projected Northern change in partial equilibrium. Percentiles and",
          r"shares are those of the general-equilibrium outputs.}",
          r"\end{tabular}", r"\end{table}"]
    return "\n".join(L) + "\n"

##################################################################################
# ARM B — LA DISTRIBUTION DES SORTIES, PAS UN ENCADREMENT
#
# Trente tirages du SEUL bras journalier. Le total journalier et la part
# gagnante ne dependent que de lui (daily arm contre la baseline sans climat),
# donc trente runs donnent une vraie distribution de ces deux quantites.
#
# ⚠️ DEUX OBJETS DANS LA MEME TABLE, ET ILS NE SE CONFONDENT PAS.
#   - rq05 / rq50 / rq95 : TROIS TIRAGES SELECTIONNES, solves un par un.
#   - Lo / Med / Hi      : les quantiles de la DISTRIBUTION DES TRENTE SORTIES.
# Le premier dit « voici le resultat au tirage qui est 5e sur cent » ; le second
# dit « voici ou tombe le resultat quand on tire trente fois ». La note doit le
# dire, sans quoi un lecteur lira deux fois le meme intervalle.
##################################################################################
##################################################################################
# LA GRILLE EST CELLE QUI A TOURNE, PAS CELLE QUI ETAIT PREVUE.
#
# Deux lots de l'arm B ont ete tues par le systeme faute de memoire. Exiger les
# trente aurait rendu l'emetteur muet alors que douze tirages exploitables
# etaient sur le disque -- et un emetteur muet ne dit pas combien il lui
# manque, il dit seulement qu'il ne sort rien.
# On prend donc les colonnes PRESENTES ET AU MANIFESTE, et le nombre effectif
# est EMIS (\GEarmBnDraws) puis imprime dans la note de la table : une
# distribution sur douze tirages est une distribution sur douze tirages, et le
# lecteur doit le lire, pas le deviner.
# ⚠️ Present au manifeste, pas seulement sur le disque : run_one_ge.sh
# n'inscrit qu'apres rc=0, donc un run tue en cours d'ecriture laisse un CSV
# que sa seule existence ne disqualifie pas.
##################################################################################
ARMB_MIN = 10
_armb_man = pd.read_csv(os.path.join(RUNS_DIR, 'ge_runs_manifest.csv'),
                        keep_default_na=False)
_armb_mani = set(_armb_man['file'])

##################################################################################
# ASSERTION B1 -- LA TAILLE ATTENDUE NE VIENT PAS DE L'ARTEFACT QUE ct5 ECRIT.
#
# Elle en venait. `ARMB_ALL` etait la liste des colonnes de ge_runs_band_armB.csv,
# c'est-a-dire du sidecar que ct5 avait ecrit juste avant. Un ct5 qui rendait la
# moitie de la grille handait donc a l'emetteur la grille COURTE **et**
# l'attente qui lui correspondait : tous les gardes en aval passaient sur trente
# colonnes -- couverture 5..95, ARMB_MIN, appartenance au manifeste -- et une
# bande de trente tirages dont vingt-neuf d'avant la division etait publiable
# sans un avertissement. C'est ce qui a failli arriver le 20/09/2026.
#
# L'attente vient desormais de ARMB_N, declare dans 00_parameters.R et repete
# dans son jumeau Python, c'est-a-dire d'un endroit ou ct5 n'ecrit pas.
##################################################################################
ARMB_N = int(_PARAMS.ARMB_N)
ARMB_ALL = ['rb%02d' % i for i in range(1, ARMB_N + 1)]
_armb_side = pd.read_csv(os.path.join(RUNS_DIR, 'ge_runs_band_armB.csv'))
if len(_armb_side) != ARMB_N or set(_armb_side['column']) != set(ARMB_ALL):
    raise SystemExit(
        "[STOP] ge_runs_band_armB.csv carries %d row(s) (%s...) while "
        "00_parameters.R declares ARMB_N = %d.\n"
        "  The band grid was written against another N than the declared "
        "one. Re-run 04_contrast/ct5_rank_mc_draws.R; N is changed in "
        "00_parameters.R and nowhere else."
        % (len(_armb_side), ', '.join(sorted(_armb_side['column'])[:4]), ARMB_N))

ARMB_COLS = [c for c in ARMB_ALL
             if f'welfare_sc2_{c}{NULL_TAG}.csv' in _armb_mani
             and os.path.exists(os.path.join(
                 RUNS_DIR, f'welfare_sc2_{c}{NULL_TAG}.csv'))]

##################################################################################
# ASSERTION B2 -- CHAQUE TIRAGE DE LA BANDE PORTE LA MEME CONVENTION QUE SA
# BASELINE, ET DONC LA MEME QUE LE RESTE DU JEU PUBLIE.
#
# CE QUE « POSTERIEUR A LA DIVISION » VEUT DIRE OPERATIONNELLEMENT. Depuis le
# 20/09/2026, `inputs_md5` (run_one_ge.sh) hache aussi la convention d'echelle
# du canal productivite, SCALE_PROD = 1/THETA_GE. Un run d'avant la division et
# un run d'apres ne peuvent donc plus partager une empreinte -- avant ce jour ils
# la partageaient (930ff8e0 pour les deux millesimes) et seul l'horodatage les
# separait. Comparer les empreintes est donc STRICTEMENT PLUS FORT qu'un filtre
# de date : la date dit quand le run a tourne, l'empreinte dit sous quelle
# convention. C'est aussi la cle que 02e_north_contrast_band.py et
# decompose_metric_channels.py utilisent deja pour dire « un seul millesime ».
#
# L'ANCRE EST LA BASELINE. Le total journalier d'un tirage est
# cumsum(log(welfare_sc2_rbXX / welfare_sc1_q50_NULLamen)) : un tirage qui ne
# partage pas la convention de sa baseline produit un rapport entre deux mondes.
##################################################################################
_armb_anchor = f'welfare_sc1_{QUANTILE}{NULL_TAG}.csv'
ARMB_VINTAGE_OK = True
ARMB_VINTAGE_WHY = ''
if ARMB_COLS:
    _need = [_armb_anchor] + [f'welfare_sc2_{c}{NULL_TAG}.csv' for c in ARMB_COLS]
    _miss = [n for n in _need if n not in _armb_mani]
    if _miss:
        ARMB_VINTAGE_OK = False
        ARMB_VINTAGE_WHY = ("missing from the manifest: %s" % ', '.join(_miss[:4]))
    else:
        _sub = _armb_man.set_index('file').loc[_need]
        _fp = sorted({str(v) for v in _sub['inputs_md5']})
        _ts = sorted(str(v) for v in _sub['timestamp'])
        if len(_fp) != 1:
            ARMB_VINTAGE_OK = False
            ARMB_VINTAGE_WHY = (
                "the band draws and their baseline carry %d input "
                "fingerprints: %s. A band made of two conventions is not "
                "a distribution, it is a mix. Re-run the draws that "
                "carry the old one." % (len(_fp), ', '.join(f[:8] for f in _fp)))
        else:
            print("[armB] single vintage: inputs_md5=%s | %d run(s) | "
                  "from %s to %s" % (_fp[0][:8], len(_need), _ts[0], _ts[-1]))
if not ARMB_VINTAGE_OK:
    raise SystemExit("[STOP] arm B band: " + ARMB_VINTAGE_WHY)
##################################################################################
# ⛔ UN PREFIXE N'EST PAS UN ECHANTILLON.
#
# Les colonnes rb01..rbNN sont ORDONNEES par percentile du classement. Les
# douze premieres couvrent les percentiles 5 a 39 : leur mediane est celle du
# bas de la grille, et « 100 % des tirages sont negatifs » y est vrai par
# construction, pas par mesure. Accepter une grille partielle sur le seul
# COMPTE publierait un biais systematique sous le nom de distribution.
# On exige donc que la grille COUVRE SA PLAGE : au moins un tirage dans le
# premier decile et un dans le dernier. Un prefixe echoue immediatement.
##################################################################################
## Le percentile d'une colonne est LU dans le sidecar que ct5 ecrit, jamais
## recalcule d'une formule sur l'indice. Une formule « 5 + 90*(i-1)/(n-1) »
## suppose une grille reguliere de bout en bout ; elle devient fausse des qu'on
## INTERCALE des tirages (rb31.. aux percentiles intermediaires), et elle le
## devient en silence -- le garde de couverture continuerait a passer en
## decrivant les mauvais percentiles.
## (!) Le sidecar reste l'autorite du PERCENTILE de chaque colonne ; il n'est
## plus l'autorite de leur NOMBRE. Les deux roles etaient confondus.
_armb_pct = dict(zip(_armb_side['column'], _armb_side['percentile']))
_missing_pct = [c for c in ARMB_ALL if c not in _armb_pct]
if _missing_pct:
    raise SystemExit(
        "[STOP] ge_runs_band_armB.csv does not carry the percentile of: %s.\n"
        "  Re-run 04_contrast/ct5_rank_mc_draws.R." % ', '.join(_missing_pct))
_pcts = [_armb_pct[c] for c in ARMB_COLS]
ARMB_SPANS = bool(_pcts) and min(_pcts) <= 10.0 and max(_pcts) >= 90.0
ARMB_OK = len(ARMB_COLS) >= ARMB_MIN and ARMB_SPANS
if ARMB_OK:
    print(f"[armB] {len(ARMB_COLS)}/{len(ARMB_ALL)} usable draws, "
          f"range {min(_pcts):.0f}-{max(_pcts):.0f} -> distribution emitted")
    if len(ARMB_COLS) < len(ARMB_ALL):
        print(f"[armB] (!) INCOMPLETE grid: {len(ARMB_ALL) - len(ARMB_COLS)} "
              f"draw(s) missing. The count emitted is the ACTUAL count.")
elif not ARMB_SPANS and ARMB_COLS:
    print(f"[armB] REFUSED: {len(ARMB_COLS)} draw(s) covering only the "
          f"percentiles {min(_pcts):.0f}-{max(_pcts):.0f}. A prefix of the "
          f"ranking is not a sample: its median and shares would be "
          f"biased by construction. Distribution not emitted.")
else:
    print(f"[armB] only {len(ARMB_COLS)} draw(s) (< {ARMB_MIN}) "
          f"-> distribution not emitted")

# Onze unites pour le CSV, huit pour la table. Le suffixe est le nom de macro.
ARMB_UNITS = [
    ('World',            'World',           lambda d: d),
    ('North>=45N',       'North',           lambda d: d[d['y'] >= 45.0]),
    ('North>=55N',       'NorthFiftyFive',  lambda d: d[d['y'] >= 55.0]),
    ('North>=66.5N',     'NorthArctic',     lambda d: d[d['y'] >= 66.5]),
    ('Russia',           'Rus',             lambda d: _pays(d, 'RUS')),
    ('Finland',          'Fin',             lambda d: _pays(d, 'FIN')),
    ('Canada',           'Can',             lambda d: _pays(d, 'CAN')),
    ('Germany',          'Deu',             lambda d: _pays(d, 'DEU')),
    ('France',           'Fra',             lambda d: _pays(d, 'FRA')),
    ('Norway',           'Nor',             lambda d: _pays(d, 'NOR')),
    ('Sweden',           'Swe',             lambda d: _pays(d, 'SWE')),
]


_ARMB_FRAME_CACHE = {}


def _armb_frame():
    """Une ligne par tirage et par unite : total journalier et part gagnante.

    MEMOISE. Trois consommateurs le demandent -- _armb_stats, _armb_rank_corr
    et l emission du CSV -- et chaque construction relit soixante matrices de
    welfare. Le cache ne change RIEN au resultat : la fonction est pure, ses
    entrees sont des fichiers que l emetteur ne reecrit pas.

    Le rang du tirage vient de ge_runs_band_armB.csv, ecrit par ct5 AVANT les
    runs : c'est le classement en equilibre partiel, et c'est lui qu'on veut
    comparer a l'ordre observe en equilibre general."""
    if 'f' in _ARMB_FRAME_CACHE:
        return _ARMB_FRAME_CACHE['f']
    side = pd.read_csv(os.path.join(RUNS_DIR, 'ge_runs_band_armB.csv'))
    side = side.set_index('column')
    W1 = _load_w(1, QUANTILE, NULL_TAG)
    out = []
    for c in ARMB_COLS:
        cum = np.cumsum(np.log(_load_w(2, c, NULL_TAG) / W1),
                        axis=0)[TERMINAL_ROW, :]
        q = pix[['id', 'y', 'x', 'pop']].copy()
        q['daily'] = cum
        for lab, _mac, sel in ARMB_UNITS:
            d = sel(q)
            out.append({
                'draw_id': int(side.loc[c, 'draw_id']),
                'draw_rank': int(side.loc[c, 'draw_rank']),
                'column': c, 'percentile': float(side.loc[c, 'percentile']),
                'unit': lab,
                'daily_total': wmean(d, 'daily'),
                'pop_gain_share': 100.0 - 100.0 * d.loc[d['daily'] < 0, 'pop'].sum()
                                  / d['pop'].sum(),
            })
    _ARMB_FRAME_CACHE['f'] = pd.DataFrame(out)
    return _ARMB_FRAME_CACHE['f']


# ======================================================================
# 5ter. L'ASSERTION TEXTE <-> TABLE POUR S5.1-S5.3
# ======================================================================
# Le defaut que ce controle ferme : le texte de S5.1-S5.3 cite des macros
# emises depuis BODY (f_a = 0) pendant que les tables en regard sortaient du
# bras interieur. Le lecteur avait sous les yeux, sur la meme page, une phrase
# et une cellule qui ne decrivaient pas le meme modele -- et rien, ni la
# compilation, ni les sept controles, ne pouvait le voir : les deux nombres
# etaient corrects, chacun dans son bras.
#
# Ce que le controle compare est ce qui est REELLEMENT ECRIT : la chaine du
# \newcommand contre la chaine de la cellule. Comparer les objets Python en
# amont ne verrait pas une table qui se remettrait a lire l autre bras.
#
# Polarite : ARRET. Une divergence ici n est pas une propriete du resultat,
# c est une page qui se contredit.

BS_MID = '\\midrule'
BS_BOT = '\\bottomrule'
BS_NL  = '\\\\'
BS_ADD = '\\addlinespace'
BS     = '\\'


def _table_cells(tex):
    """rowlabel -> [cellules], depuis le corps d un tabular emis."""
    out = {}
    body = tex.split(BS_MID)[-1].split(BS_BOT)[0] if BS_MID in tex else ''
    for line in body.split(BS_NL):
        line = line.strip()
        if line.startswith(BS_ADD):
            # La commande d espacement addlinespace est emise SANS
            # separateur de ligne : apres le decoupage elle reste collee en
            # tete de la ligne suivante. La retirer, et non sauter la ligne,
            # sinon la ligne d agregat du Nord sortait du controle.
            line = line[len(BS_ADD):].strip()
        if not line or line.startswith('%'):
            continue
        parts = [c.strip() for c in line.split('&')]
        if len(parts) < 2:
            continue
        lab = parts[0]
        for junk in (BS + 'textbf{', BS + 'mathbf{', '$', '}'):
            lab = lab.replace(junk, '')
        out[lab.strip()] = [_clean_cell(c) for c in parts[1:]]
    return out


def _clean_cell(c):
    c = c.strip()
    for junk in (BS + 'textbf{', BS + 'mathbf{', '$', '{', '}'):
        c = c.replace(junk, '')
    return c.strip()


def _macros_of(payload):
    """nom -> valeur, depuis le values.tex qu on s apprete a ecrire."""
    out = {}
    for line in payload.split(chr(10)):
        line = line.strip()
        if not line.startswith(BS + 'newcommand{' + BS):
            continue
        try:
            nm = line.split(BS + 'newcommand{' + BS, 1)[1].split('}', 1)[0]
            val = line.split('}{', 1)[1].rsplit('}', 1)[0]
        except (IndexError, ValueError):
            continue
        out[nm] = val.strip()
    return out


# (macro, fichier de table, etiquette de ligne, index de colonne 0-based)
# Les etiquettes de ligne des horizons viennent de HORIZON_ROW, donc l etiquette
# terminale suit HORIZON_LABEL sans etre epelee ici.
def _body_text_table_pairs():
    _h = list(HORIZON_ROW)                       # ['2030', '2050', <label>]
    pairs = []
    for lab in _h:
        w = HORIZON_WORD[lab]
        for col, tag in enumerate(('Daily', 'Annual', 'Extra')):
            pairs.append((f'GEagg{tag}{w}', 'tab_ge_aggregate.tex', lab, col))
    # Les trois agregats mondiaux du corps SONT la ligne terminale de la table.
    for col, tag in enumerate(('Daily', 'Annual', 'Extra')):
        pairs.append((f'GEworld{tag}', 'tab_ge_aggregate.tex', _h[-1], col))
    ## LES QUARTILES NE SONT PLUS ICI, et la garantie n'a pas disparu avec
    ## eux : elle a change de forme. Tant que tab_ge_tiers existait, la macro
    ## et la cellule etaient deux ECRITURES a comparer. Depuis le lot 10 le
    ## chiffre se montre dans une FIGURE, et une figure n'a pas de cellule a
    ## relire : la comparaison de chaines n'a plus d'objet. Ce qui la remplace
    ## est plus fort -- il n'y a plus qu'un seul chemin de calcul.
    ## emit_figures.py ne recalcule AUCUNE barre : il trace fig_ge_tiers.csv,
    ## ecrit ici, dans le meme passage que les macros.
    for iso, name in REVERSAL + AMPLIFICATION + FRANCE:
        nm = MACRO_NAME[iso]
        for col, tag in enumerate(('Net', 'Annual', 'DailyExtra')):
            pairs.append((f'GEcty{nm}{tag}', 'tab_ge_north.tex', name, col))
    return pairs


def assert_body_text_matches_tables(payload, tables):
    mac = _macros_of(payload)
    cells = {k: _table_cells(v) for k, v in tables.items()}
    bad, checked, skipped = [], 0, 0
    for name, fil, row, col in _body_text_table_pairs():
        if name not in mac or fil not in cells:
            skipped += 1
            continue
        if mac[name].startswith(BS + 'PROVISOIRE'):
            skipped += 1
            continue
        r = cells[fil].get(row)
        if r is None or col >= len(r):
            bad.append((name, fil, row, col, mac[name], '<row/column missing>'))
            continue
        checked += 1
        if mac[name] != r[col]:
            bad.append((name, fil, row, col, mac[name], r[col]))
    if bad:
        msg = ["STOP [text<->table S5.1-S5.3]: %d divergence(s)." % len(bad),
               "  A macro of the body and the cell it accompanies do not carry",
               "  the same value: tables on the interior arm under a text",
               "  at f_a = 0.",
               "  %-26s %-24s %-14s %-9s %-9s" % ('macro', 'table', 'row',
                                                  'macro', 'cell')]
        for name, fil, row, col, mv, cv in bad:
            msg.append("  %-26s %-24s %-14s %-9s %-9s"
                       % (BS + name, fil, row + '/c%d' % col, mv, cv))
        sys.exit(chr(10).join(msg))
    print("[texte<->table] %d cell(s) of S5.1-S5.3 equal to their macro "
          "(%d not applicable)" % (checked, skipped))


# ======================================================================
# 6. Writing
# ======================================================================
##################################################################################
# tab_ge_tiers.tex N'EST PLUS EMISE (18/09/2026, lot 10).
#
# L'auteur a retire la sous-section « Income quartiles » et son input : les
# quartiles passent au CORPS, sous la forme de fig_ge_tiers.pdf. La table
# n'avait donc plus de lecteur, et une table sans lecteur qui continue de
# s'ecrire finit par etre citee.
#
# emit_tab_tiers() EST CONSERVEE, non appelee : elle porte la mise en forme
# exacte de l'exhibit retire, et la reactiver est une ligne. Le fichier
# output/tab_ge_tiers.tex du 17/09 est conserve dans le snapshot du lot.
# ⚠️ Il n'est PAS supprime d'output/ par cet emetteur : rien ici n'efface, et
#    un fichier orphelin qui traine est un probleme de menage, pas de sortie.
##################################################################################
TABLES = {
    'tab_ge_aggregate.tex': emit_tab_aggregate(),
    'tab_ge_north.tex':     emit_tab_north(),
    'tab_ge_latitude.tex':  emit_tab_latitude(),
}
if HAVE_NULL:
    TABLES['tab_ge_sensitivity.tex'] = emit_tab_sensitivity()
    TABLES['tab_ge_north_sens.tex']  = emit_tab_north_sens()
    # Table 2 du corps : elle n'existe que si le bras f_a = 0 est present. Son
    # absence est ANNONCEE plutot que comblee au q50.
    TABLES['tab_ge.tex'] = emit_tab_ge()
if BAND_OK:
    # La bande n'existe que si les six runs de quantile sont la. Son absence est
    # ANNONCEE : sans elle le manuscrit s'arrete sur un \input introuvable, ce
    # qui est la bonne polarite -- mieux qu'une table muette.
    TABLES['tab_ge_band.tex'] = emit_tab_ge_band()
if RANKED_OK:
    TABLES['tab_ge_band_ranked.tex'] = emit_tab_ge_band_ranked()

# L'ASSERTION TOURNE AVANT L'ECRITURE, ET AUSSI EN DRY-RUN.
# En dry-run surtout : c'est le mode ou l'on verifie un emetteur qu'on vient
# de modifier, donc le mode ou une divergence doit se voir. Un controle qui ne
# s'execute que sur le chemin d'ecriture ne protege pas la relecture.
_PAYLOAD = emit_values_tex()
assert_body_text_matches_tables(_PAYLOAD, TABLES)

if DRY_RUN:
    print("\n--- DRY RUN: nothing is written ---\n")
    print(_PAYLOAD)
    for k, v in TABLES.items():
        print(f"\n===== {k} =====\n{v}")
else:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ------------------------------------------------------------------
    # CALCULER D'ABORD, ECRIRE ENSUITE, REMPLACER ATOMIQUEMENT.
    # ------------------------------------------------------------------
    # Le code precedent etait :
    #     with open(VALUES_TEX, 'w') as fh:
    #         fh.write(emit_values_tex())
    # Python evalue `open(..., 'w')` AVANT l'appel : le fichier est tronque a
    # zero avant que le contenu n'existe. Toute garde qui s'arrete dans
    # emit_values_tex() -- le refus md5 de decompose_metric_channels, par
    # exemple -- detruisait donc values.tex au lieu de le laisser tel quel.
    # Mesure du 11/09/2026 (test known-bad) : exit 2, values.tex a 0 octet,
    # 104 macros perdues. Un arret doit laisser l'artefact PRECEDENT intact.
    # Deja calcule plus haut, pour que l'assertion texte<->table porte sur
    # EXACTEMENT la chaine qu'on ecrit, et non sur une seconde evaluation.
    _payload = _PAYLOAD                   # <- calcul AVANT toute ouverture

    ## Lot 86d : LF explicite sur tout systeme (le mode texte ecrivait CRLF sous
    ## Windows, et \r\r\n pour les CSV, dont to_csv() rend deja des \r\n).
    def _write_atomic(path, text):
        tmp = path + '.tmp'
        with open(tmp, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(text)
        os.replace(tmp, path)             # atomique sur le meme volume

    ##################################################################
    # values_provenance.tsv -- DE QUEL RUN SORT CHAQUE MACRO
    #
    # Sans ce fichier, « cette macro est-elle a f_a = 0 ? » ne se repond qu'en
    # relisant l'emetteur. C'est ce qui a permis a trois macros du Nord de
    # porter le bras interieur sous un nom que le corps presentait comme la
    # borne basse. tools/check_macros.py le croise avec le .tex.
    #
    # Les etiquettes de run viennent du MANIFESTE, jamais d'une constante :
    # ge_runs_manifest.csv est la seule source des `amenity` et `cap_colmass`.
    ##################################################################
    # ⚠️ keep_default_na=False : l'etiquette d'amenite du bras borne basse est
    # la CHAINE « null », et pandas la lit comme NaN par defaut. La provenance
    # sortait alors « amenity=nan » sur toutes les macros du corps, et le
    # controle les refusait toutes -- un faux rouge qui ressemble a un vrai.
    _mani = pd.read_csv(os.path.join(RUNS_DIR, 'ge_runs_manifest.csv'),
                        keep_default_na=False)

    ##################################################################
    # LES CONSTANTES N'ONT PAS DE BRAS.
    # La grille, le nombre de periodes, l'etiquette d'horizon : elles ne
    # descendent d'aucun run et les qualifier d'un bras serait faux dans les
    # deux sens. Elles sont marquees `invariant`, et le controle les ignore.
    ##################################################################
    INVARIANT = {
        'GEncells', 'GEncellsPotential', 'GEncellsPopOnly', 'GEnperiods',
        'GEntransitions', 'GEhorizonLabel', 'GEhorizonSpan', 'GEdensityCut',
        'GEdensityCutUnit', 'GEnruns', 'GEarmLabel', 'GEcapColmass',
        'GEnorthLat', 'GEmigSenders', 'GEmigReceivers', 'GEclearGap',
        'GEamenitySpecBody',
        # theta est une CALIBRATION, pas un resultat : elle est lue dans
        # specification_study.csv et ne descend d'aucun run. La qualifier d'un
        # bras serait faux dans les deux sens, exactement comme la grille.
        'GEtheta',
        # la liste des GCM : une donnee d'entree, pas un resultat de run.
        'GEcmipModels',
        # le domaine (lot 68) : la grille, la population et Abel-Cohen, pas un run.
        'GEdomCellsOut', 'GEdomPopOutM', 'GEdomPopOutShare', 'GEdomCtyOutN',
        'GEdomCtyOutPopM', 'GEdomCtyPartialN', 'GEdomMigLostM',
        # SSP2 : comptes de modeles par pays, ecrits par 02a -- une donnee d'entree.
        'GEnsspOneModel', 'GEnsspTwoModels', 'GEnsspThreeModels',
    }
    PROVENANCE[:] = [(nm, val, 'invariant' if nm in INVARIANT else prov)
                     for nm, val, prov in PROVENANCE]

    def _run_row(tag):
        """La ligne de manifeste du scenario 2 de ce bras : elle porte le
        run_id, l'amenite et le cap qui valent pour tout le bras."""
        cible = f'welfare_sc2_{QUANTILE}{tag}.csv'
        hit = _mani[_mani['file'] == cible]
        if hit.empty:
            return dict(run_id='?', amenity='?', cap='?')
        r = hit.iloc[0]
        return dict(run_id=str(r['run_id']), amenity=str(r['amenity_spec']),
                    cap=str(r['cap_colmass']))

    RUN_PROV = {'headline': _run_row(CAP_TAG),
                'null': _run_row(NULL_TAG) if HAVE_NULL else
                        dict(run_id='-', amenity='-', cap='-'),
                'invariant': dict(run_id='-', amenity='invariant', cap='-'),
                'none': dict(run_id='-', amenity='-', cap='-')}

    _lignes = ['macro\tvalue\trun_id\tamenity\tcap_colmass\tdraw']
    for _nm, _val, _prov in PROVENANCE:
        _r = RUN_PROV.get(_prov, RUN_PROV['none'])
        _lignes.append('\t'.join((_nm, _val, _r['run_id'], _r['amenity'],
                                  _r['cap'], '-')))
    _prov_txt = '\n'.join(_lignes) + '\n'

    _write_atomic(VALUES_TEX, _payload)
    print(f"\n[write] {VALUES_TEX}")
    _write_atomic(os.path.join(OUTPUT_DIR, 'values_provenance.tsv'), _prov_txt)
    _n_null = sum(1 for _, _, p in PROVENANCE if p == 'null')
    print(f"[write] output/values_provenance.tsv: {len(PROVENANCE)} macros "
          f"({_n_null} at amenity=none)")
    for k, v in TABLES.items():
        _write_atomic(os.path.join(OUTPUT_DIR, k), v)
        print(f"[write] output/{k}")
    for k, v in TABLES_CSV.items():
        _write_atomic(os.path.join(OUTPUT_DIR, k), v)
        print(f"[write] output/{k}")

if HAVE_NULL:
    print("\n[note] tab_ge_sensitivity / tab_ge_north_sens: EMITTED (null g=0 runs present).")
else:
    print("\n[note] tab_ge_sensitivity / tab_ge_north_sens: NOT emitted "
          "(null AMENITY_SPEC=null runs missing from this arm).")
