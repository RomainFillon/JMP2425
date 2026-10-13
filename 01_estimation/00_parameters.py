# -*- coding: utf-8 -*-
r"""
00_parameters.py -- the Python twin of 00_parameters.R, for the horizon labels.

It carries ONLY what a Python producer needs and an R producer already has, so
that a single decision is not taken twice in two languages. It deliberately
does NOT mirror the estimation parameters: those are read by R alone, and a
second copy of them would be a second thing to keep in step.

TWO GROUPS OF CONSTANTS LIVE HERE, AND FOR THE SAME REASON: the horizon labels,
and -- since 2026-09-18 -- the CONTRAST GRID. A Python producer counts the grid
(emit_values_empirics.py fills \EMPgridProdNspec and the three \...Wide macros;
ct3 names the reference specification), and it cannot source 00_parameters.R.
The duplication is paid for in tools/check_macros.py, control 9, which reads
both files and fails if the two lists differ by one entry.

THE LABEL IS NOT THE COMPUTED YEAR, AND THE TWO ARE DELIBERATELY DIFFERENT.

The quinquennial grid labels each block by its FIRST year, so the last FULL
block is {2095-2099} and carries the label 2095 inside the computation. A block
labelled 2100 would hold a single year and be degenerate, so the model never
evaluates one: TERMINAL_ROW, HORIZON_ROW and every FUT_YEAR stay where they
are, and none of them is touched by what follows.

What the READER is shown is the END of that terminal block, 2100. These four
constants are the ONLY place that decision is taken on the Python side, and
they are checked against the R side by tools/check_horizon.py.

Import idiom, from any producer directory:

    import os, sys
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        '01_estimation'))
    from importlib import import_module
    _P = import_module('00_parameters')      # the name starts with a digit
    HORIZON_LABEL = _P.HORIZON_LABEL
"""

HORIZON_LABEL = "2100"          # a level AT the horizon: "by 2100"
HORIZON_SPAN = "2015--2100"     # a trajectory: en dash, LaTeX source form
HORIZON_SPAN_PLAIN = "2015-2100"  # same, for plot labels and file names
HORIZON_BLOCK = "2095-2099"     # the terminal block itself, when named

# ---------------------------------------------------------------------------
# THE CONTRAST GRID -- the Python twin of GRID_SPECS in 00_parameters.R
# ---------------------------------------------------------------------------
# Author's decision of 2026-09-18: fifteen specifications, not twenty-four.
# Polynomials of degree two to four and natural splines with three and four
# degrees of freedom, each at three symmetric winsorisation thresholds. The
# natural splines with five to seven degrees of freedom are out: they oscillate
# in the DENSE part of the support, 20-28 C.
#
# This is NOT the grid of the body. The body's grid -- 01a, 01b, 01c, the Akaike
# weights, the GE, the CRH appendix -- has four forms and twelve specifications
# and is SPEC_CONFIGS in 00_parameters.R. It is unchanged. The one below is the
# grid of the annual/daily contrast alone: ct1, ct2, ct4, ct6.
#
# The label format is the one every artefact on disk already carries,
# "spline4_wl0.5wu0.5": %g on the thresholds, so 1.0 prints as "1".
GRID_FORMS = (("polynomial", 2), ("polynomial", 3), ("polynomial", 4),
              ("spline", 3), ("spline", 4))
GRID_WINSOR = (0.5, 1.0, 2.5)


def grid_spec_label(spec, degree, w_lo, w_hi):
    """The label of one specification. Same construction as ct1's sprintf."""
    return "%s%d_wl%gwu%g" % (spec, int(degree), w_lo, w_hi)


GRID_SPECS = tuple(
    {"spec": spec, "degree": deg, "winsor_lo": w, "winsor_hi": w,
     "spec_label": grid_spec_label(spec, deg, w, w)}
    for spec, deg in GRID_FORMS for w in GRID_WINSOR)

GRID_NSPEC = len(GRID_SPECS)
GRID_LABELS = tuple(z["spec_label"] for z in GRID_SPECS)

# The reference specification of the contrast, designated before any result:
# ns df = 4 at 0.5% winsorisation, pre-engagement (P4) of ct2. It is the
# specification of the published paper, and 0.5% is the least clipping of the
# three thresholds. NOT the "REF_SPEC" of 00_parameters.R, which is the
# polynomial of the identification diagnostics -- a different object, hence a
# different name.
GRID_REF_SPEC = grid_spec_label("spline", 4, 0.5, 0.5)

# The count is the product, and the labels are distinct. Checked at import, so
# a consumer never counts a grid that lost an entry to a typo.
assert GRID_NSPEC == len(GRID_FORMS) * len(GRID_WINSOR)
assert len(set(GRID_LABELS)) == GRID_NSPEC
assert GRID_REF_SPEC in GRID_LABELS


# The two clamp points of the CRH sensitivity, the twin of the R declaration.
# There is NO clamp in the published projection: project_channel() runs with
# reliable_hi = Inf, and these two are the sensitivities reported against it.
#   CRH_CLAMP_SUPPORT is the edge of the well-populated support, the threshold
#     the manuscript names \EMPcrhSupportHi;
#   CRH_CLAMP_SEVERE is the stricter one asked for in lot 10.
# They are floats here because they index a CSV column written by the R side:
# the lookup is what checks that the two twins still agree.
CRH_CLAMP_SUPPORT = 40.0
CRH_CLAMP_SEVERE = 30.0

# ---- Bande GE "arm B": the number of draws, the twin of the R declaration ----
# ARMB_N is declared in 00_parameters.R and NOWHERE ELSE. It is repeated here
# for the Python consumers exactly as GRID_SPECS is, and check_macros control 9
# refuses a divergence between the two files.
#
# WHY A PYTHON TWIN AT ALL. emit_values.py used to read the SIZE of the band
# grid out of ge_runs_band_armB.csv -- the very artefact ct5 writes. A ct5 that
# wrote half the grid therefore handed the emitter a smaller grid AND the
# expectation that matched it: every downstream guard passed on 30 columns.
# The expectation has to come from somewhere ct5 cannot write to.
ARMB_BASE_N = 30
ARMB_N = 60
assert ARMB_N in (ARMB_BASE_N, 2 * ARMB_BASE_N)

# ---- Omega and the CMIP6 ensemble: twins of the R declarations (lot 38) ------
# OMEGA is the migration elasticity of the estimation, 00_parameters.R `Omega`
# (SCALE_AMEN = Omega). emit_values_empirics.py used to type 0.5 itself and cite
# data_inputs/tab_calib_params.tex, a file withdrawn on 27/09/2026.
# CMIP6_MODELS: the four ISIMIP3b GCMs as the file names spell them; the GE
# emitter prints them as \GEcmipModels. check_macros control 9 refuses a
# divergence from the R side for both.
OMEGA = 0.50
CMIP6_MODELS = ("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0")


__all__ = ["HORIZON_LABEL", "HORIZON_SPAN", "HORIZON_SPAN_PLAIN",
           "HORIZON_BLOCK",
           "GRID_FORMS", "GRID_WINSOR", "grid_spec_label", "GRID_SPECS",
           "GRID_NSPEC", "GRID_LABELS", "GRID_REF_SPEC",
           "CRH_CLAMP_SUPPORT", "CRH_CLAMP_SEVERE",
           "ARMB_BASE_N", "ARMB_N", "OMEGA", "CMIP6_MODELS"]
