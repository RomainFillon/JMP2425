# -*- coding: utf-8 -*-
"""
main.py — Job market paper, Romain Fillon
"""
import os
import time
import sys
import importlib.util
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix

## The root is derived from the location of this file, never from a development
## path: main.py sits in 02_GE/02c_model/, so the root is two levels up.
## JMP_ROOT wins when set, and the orchestrator exports it.
_HERE = os.path.dirname(os.path.abspath(__file__))          # .../02_GE/model/
JMP_ROOT = os.environ.get('JMP_ROOT') or os.path.dirname(os.path.dirname(_HERE))

data_folder = os.environ.get('JMP_DATA',
                             os.path.join(JMP_ROOT, 'intermediate_outputs'))
out_folder  = os.environ.get('JMP_OUT',
                             os.path.join(JMP_ROOT, 'output'))
run_folder  = _HERE          # ou vivent data_cleaning.py et forward.py

## No os.chdir. It would make every relative write depend on the current
## directory, so a write call would no longer say where it writes. The paths
## above are explicit instead.
sys.path.append(run_folder)

for _n, _p in (('data_folder', data_folder), ('out_folder', out_folder)):
    if not os.path.isdir(_p):
        raise SystemExit(f"[STOP] {_n} absent : {_p}")

SCENARIO     = int(os.environ.get('SCENARIO', '2'))   # env-overridable : 1=no CC | 2=daily DRF | 3=annual DRF

##############################################################################
# LE CONSOMMATEUR VALIDE SON ENTREE. Ajoute le 21/09/2026.
#
# « sc4 » est une ETIQUETTE DE SORTIE, pas une valeur d'entree : elle est
# calculee plus bas par
#     _sclabel = 4 if (SCENARIO == 3 and ANNUAL_DRF_SOURCE == 'dailycurve')
# et un run sc4 se lance donc avec SCENARIO=3 + ANNUAL_DRF_SOURCE=dailycurve.
#
# Le 20/09/2026 run_lot11_complement.sh exportait SCENARIO=4. Aucune des trois
# branches de dommage ne mordait, damage_amenity_diff n'etait jamais liee, et
# le solveur mourait en UnboundLocalError l.412 -- QUARANTE SECONDES plus tard,
# apres avoir charge et diagnostique toutes les donnees. Le pilote s'arretant
# au premier echec, ce seul mot a annule douze runs GE et quatre producteurs.
#
# Le defaut n'est pas que le script se soit trompe : c'est que le consommateur
# ne regardait pas son entree. Trois lignes ici valent quarante secondes et un
# lot.
#
# Polarite : ARRET, avant tout chargement.
##############################################################################
if SCENARIO not in (1, 2, 3):
    raise SystemExit(
        f"[STOP] SCENARIO={SCENARIO} : valeurs acceptees 1, 2, 3.\n"
        "  1 = pas de rechauffement | 2 = dose-reponse journaliere | "
        "3 = dose-reponse annuelle.\n"
        "  Le « sc4 » des noms de fichiers est une etiquette CALCULEE : pour "
        "le produire,\n"
        "  poser SCENARIO=3 et ANNUAL_DRF_SOURCE=dailycurve.")

QUANTILE = os.environ.get('QUANTILE', 'q50')

# -- per-channel quantiles --------------------------------------------------
# Setting both channels to q50 reproduces the baseline run exactly. They can be
# decoupled to run one channel at a pessimistic quantile while the other stays
# central.

QUANTILE_PROD = os.environ.get('QUANTILE_PROD', QUANTILE)
QUANTILE_AMENITY = os.environ.get('QUANTILE_AMENITY', QUANTILE)
_qtag = f"_{QUANTILE_PROD}" if QUANTILE_PROD == QUANTILE_AMENITY \
        else f"_pp{QUANTILE_PROD}_pa{QUANTILE_AMENITY}"


print(f"[quant] QUANTILE_PROD={QUANTILE_PROD}  QUANTILE_AMENITY={QUANTILE_AMENITY}  _qtag='{_qtag}'")

# -- amenity DRF specification ------------------------------------------------
# AMENITY_SPEC selects which amenity dose-response lookup table is read and the
# suffix appended to the output filenames. The harness is fully general: any spec
# is just a CSV f(T) over the -100..100 C grid (same schema as the baseline file),
# read and integrated over the daily/annual distribution by forward.py. To add a
# new amenity DRF variant, register it here -- no other code change is needed.
#   'baseline' -> MC_damage_amenity_{daily,annual}_quantiles.csv       (estimated DRF)
#   'null'     -> MC_damage_amenity_{daily,annual}_quantiles_NULL.csv  (g=0 lower bound)
#   'bell'     -> MC_damage_amenity_{daily,annual}_quantiles_BELL.csv    (C1: bell, Albouy unclamped)
#   'bell_cap' -> MC_damage_amenity_{daily,annual}_quantiles_BELLCAP.csv (C1: bell clamped |g|<=|f_z|)
AMENITY_SPECS = {
    'baseline': {'file_tag': '',        'out_tag': ''},
    'null':     {'file_tag': '_NULL',   'out_tag': '_NULLamen'},
    'bell':     {'file_tag': '_BELL',   'out_tag': '_BELLamen'},
    'bell_cap': {'file_tag': '_BELLCAP','out_tag': '_BELLCAPamen'},
}

AMENITY_SPEC = os.environ.get('AMENITY_SPEC', 'baseline')   # env-overridable

# -- source of the annual dose-response --------------------------------------
# Applies to SCENARIO 3 only. It selects which f(T) table is evaluated at the
# ANNUAL MEAN temperature; the aggregation itself does not change:
#   'annual'     -> MC_damage_{prod,amenity}_annual_quantiles*.csv
#   'dailycurve' -> MC_damage_{prod,amenity}_daily_quantiles*.csv
# Both tables share the same schema (temperature_C, q05..q95, y_mean; 201 points
# over [-100, 100]), so the annual damage function consumes either unchanged.
ANNUAL_DRF_SOURCES = {
    'annual':     {'grid': 'annual', 'out_tag': ''},
    'dailycurve': {'grid': 'daily',  'out_tag': '_dailycurve'},
}
ANNUAL_DRF_SOURCE = os.environ.get('ANNUAL_DRF_SOURCE', 'annual')
_ADRF   = ANNUAL_DRF_SOURCES[ANNUAL_DRF_SOURCE]
_drftag = _ADRF['out_tag']
print(f"[drf] ANNUAL_DRF_SOURCE={ANNUAL_DRF_SOURCE}  grid='{_ADRF['grid']}'  _drftag='{_drftag}'")

# -- welfare dump at each outer iteration (convergence diagnostic) -----------
# DUMP_OUTER=1 writes welfare_outer<k><suffix>.csv at each iteration k. The
# default of '0' writes nothing and changes no behaviour.
DUMP_OUTER = os.environ.get('DUMP_OUTER', '0').strip() == '1'
print(f"[dump] DUMP_OUTER={DUMP_OUTER}")

#AMENITY_SPEC = 'baseline'   # 'baseline' | 'null' | 'bell' | 'bell_cap'
_AMEN        = AMENITY_SPECS[AMENITY_SPEC]
NULL_AMENITY = (AMENITY_SPEC == 'null')   # backward-compat alias


# -- trade concentration cap (invariance test) -------------------------------
# CAP_COLMASS bounds the column mass of lambda0, that is, the market share a
# single seller can hold. It is driven by an environment variable so the runs of
# the sweep need no code edit. The default is the headline value; the control
# arm sets CAP_COLMASS=none.
_cap_env    = os.environ.get('CAP_COLMASS', '1.5').strip().lower()
CAP_COLMASS = None if _cap_env in ('none', 'off', '') else float(_cap_env)
_captag = '' if CAP_COLMASS is None else '_cap' + str(CAP_COLMASS).replace('.', 'p')
print(f"[cap] CAP_COLMASS={CAP_COLMASS}  _captag='{_captag}'")

# load modules
def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
        print(f"{name} loaded successfully")
    except Exception as e:
        print(f"Failed to load {name}: {e}")
        raise
    return mod

data_cleaning = _load_module("data_cleaning", os.path.join(run_folder, "data_cleaning.py"))
forward       = _load_module("forward",       os.path.join(run_folder, "forward.py"))

from data_cleaning import (
    load_study_parameters, load_population, load_growth_paths, load_migration,
    load_h5, data_storage)
from forward import (
    compute_migrations_populations, compute_amenity_hat, compute_productivity_hat, damage_function_exogenous_daily, 
	damage_function_exogenous_annual, compute_annual_temperature,
    solve_temporary_equilibrium) 
from forward import (compute_productivity_hat, compute_amenity_hat_t,
                     step_migration, solve_temporary_equilibrium_t)
 
def print_quantiles_sparse(name, mat):

    import numpy as np

    print(f"\n--- {name} quantiles ---")

    data = mat.data  # non-zero entries only

    qs = [0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1]

    for q in qs:
        print(f"q{q:.3f}: {np.quantile(data, q):.3e}")

    print(f"mean: {data.mean():.3e}")
    print(f"max: {data.max():.3e}")

def check_row_sums(name, mat):

    row_sums = np.array(mat.sum(axis=1)).flatten()

    print(f"\n--- {name} row sums ---")
    print(f"min: {row_sums.min():.3e}")
    print(f"max: {row_sums.max():.3e}")
    print(f"mean: {row_sums.mean():.3e}")

    print(f"q01: {np.quantile(row_sums, 0.01):.3e}")
    print(f"q50: {np.quantile(row_sums, 0.50):.3e}")
    print(f"q99: {np.quantile(row_sums, 0.99):.3e}")

    print(f"zeros: {(row_sums == 0).sum()}")
    print(f">10: {(row_sums > 10).sum()}")

# =========================
# FULL DATA DIAGNOSTICS
# =========================
def run_diagnostics(population, migration0, lambda0,
                    growth_pop, growth_prod,
                    utility_hat_guess):

    import numpy as np

    print("\n" + "="*60)
    print("?? FULL DATA DIAGNOSTICS")
    print("="*60)

    def check_array(name, arr):
        print(f"\n--- {name} ---")
        print(f"shape: {arr.shape}")
        print(f"dtype: {arr.dtype}")

        # NaN / inf
        print(f"NaN count: {np.isnan(arr).sum()}")
        print(f"Inf count: {np.isinf(arr).sum()}")

        # stats
        print(f"min: {np.nanmin(arr):.3e}")
        print(f"max: {np.nanmax(arr):.3e}")
        print(f"mean: {np.nanmean(arr):.3e}")

        # zeros / negatives
        print(f"zeros: {(arr == 0).sum()}")
        print(f"negatives: {(arr < 0).sum()}")

    # =========================
    # BASIC ARRAYS
    # =========================
    check_array("population", population)
    check_array("growth_pop", growth_pop)
    check_array("growth_prod", growth_prod)
    check_array("utility_hat_guess", utility_hat_guess)

    # =========================
    # MIGRATION MATRIX
    # =========================
    print("\n--- migration0 ---")
    print(f"type: {type(migration0)}")
    print(f"shape: {migration0.shape}")

    if hasattr(migration0, "nnz"):
        print(f"non-zero entries: {migration0.nnz}")
        print(f"sparsity: {migration0.nnz / (migration0.shape[0]**2):.3e}")

    # row sums
    row_sums = np.array(migration0.sum(axis=1)).flatten()
    print(f"row sum min: {row_sums.min():.3e}")
    print(f"row sum max: {row_sums.max():.3e}")

    # suspicious rows
    print(f"rows with sum=0: {(row_sums == 0).sum()}")
    print(f"rows with sum>10: {(row_sums > 10).sum()}")

    # =========================
    # LAMBDA (trade shares)
    # =========================
    print("\n--- lambda0 ---")
    print(f"type: {type(lambda0)}")
    print(f"shape: {lambda0.shape}")

    row_sums = np.array(lambda0.sum(axis=1)).flatten()
    print(f"row sum min: {row_sums.min():.6f}")
    print(f"row sum max: {row_sums.max():.6f}")
    print(f"mean row sum: {row_sums.mean():.6f}")

    print(f"rows far from 1: {(np.abs(row_sums - 1) > 1e-2).sum()}")

    # =========================
    # CRITICAL ECONOMIC CHECKS
    # =========================
    print("\n--- ECONOMIC CHECKS ---")

    # population zeros
    print(f"locations with zero population: {(population[:,0] == 0).sum()}")

    # extreme growth
    print(f"growth_pop > 2: {(growth_pop > 2).sum()}")
    print(f"growth_pop < 0: {(growth_pop < 0).sum()}")

    print(f"growth_prod > 2: {(growth_prod > 2).sum()}")
    print(f"growth_prod < 0: {(growth_prod < 0).sum()}")

    # =========================
    # DIMENSION CONSISTENCY
    # =========================
    print("\n--- DIMENSION CHECK ---")

    N = population.shape[0]

    print(f"N (population): {N}")
    print(f"migration0: {migration0.shape}")
    print(f"lambda0: {lambda0.shape}")

    if migration0.shape[0] != N:
        print("? migration0 mismatch")

    if lambda0.shape[0] != N:
        print("? lambda0 mismatch")

    print("\n? DIAGNOSTICS DONE")
    print("="*60)


def main():
    indices_drop_id = pd.read_csv(
        os.path.join(data_folder, 'indices_drop.csv')
    ).values.flatten()
    indices_drop = indices_drop_id - 1   # 0-based

    # -- parameters & storage --------------------------------------------------
    params = load_study_parameters(run_folder)

    T            = params['T']
    nb_locations = len(indices_drop_id)


    utility_hat_guess, w_hat, w, income_hat = data_storage(
        T, nb_locations, data_folder, indices_drop_id)

    # -- data ------------------------------------------------------------------
    population = load_population(data_folder, T, indices_drop_id)
    growth_population_exogenous, growth_productivity_exogenous = load_growth_paths(data_folder, T)





    migration0 = load_migration(
        os.path.join(data_folder, 'mig_flows.h5'),
        indices_drop)


    lambda0 = load_h5(
        os.path.join(data_folder, 'lambda.h5'),
        indices_drop, "lambda", cap_colmass=CAP_COLMASS)

    income0   = w[0,:] * population[:,0]
    w_implied = lambda0.T.dot(income0) / np.maximum(population[:,0], 1e-8)
    w_implied = w_implied / np.mean(w_implied)
    ratio     = w_implied / np.maximum(w[0,:], 1e-8)
    print(f"ratio w_implied/w[0]: median={np.median(ratio):.4f}, "
          f"min={ratio.min():.4f}, max={ratio.max():.4f}")


    mask_ok = (ratio > 0.01) & (ratio < 100)
    print(f"ratio sans outliers: median={np.median(ratio[mask_ok]):.4f}, "
          f"n_outliers={(~mask_ok).sum()}, n_ok={mask_ok.sum()}")

    # Identify the outliers
    outlier_idx = np.where(~mask_ok)[0]
    print(f"Top outliers:")
    for i in outlier_idx[:10]:
        print(f"  loc {i}: w={w[0,i]:.4f}, w_implied={w_implied[i]:.4f}, "
              f"ratio={ratio[i]:.2f}, pop={population[i,0]:.0f}")


    # Solve w* = lambda0.T * (w* * L) / L, which is an eigenvalue problem on
    # lambda0.T weighted by population.

    lambda0_csr   = lambda0.tocsr()
    diag          = np.array(lambda0_csr.diagonal())
    zero_idx      = np.where(diag == 0)[0]

    zero_orig_ids = indices_drop[zero_idx[:10]] + 1  # retour en 1-based
    print(f"IDs 1-based : {zero_orig_ids}")

    indices_drop_csv = pd.read_csv(os.path.join(data_folder, 'indices_drop.csv')).values.flatten()
    print(f"Dans indices_drop.csv : {np.isin(zero_orig_ids, indices_drop_csv).sum()} / 10")

    import h5py
    with h5py.File(os.path.join(data_folder, 'lambda.h5'), 'r') as f:
        for orig_id in zero_orig_ids[:5]:
            row = f['lambda'][orig_id - 1, :]
            print(f"ID {orig_id}: diag={row[orig_id-1]:.8f}, sum={row.sum():.6f}")

    # -- climate damages -------------------------------------------------------
    if SCENARIO == 1:
        N = population.shape[0]
        damage_amenity_diff      = np.ones((N, T - 1))
        damage_productivity_diff = np.ones((N, T-1))

    elif SCENARIO == 2:
        tmean = pd.read_csv(os.path.join(data_folder, 'tas_ssp2_daily.csv'), header=0)
        tmean = tmean[tmean['id'].isin(indices_drop_id.tolist())]
        print(f"tmean colonnes: {tmean.columns.tolist()}")
        print(f"tmean shape: {tmean.shape}")
        print(tmean.head(3))
    

        amenity_file = f"MC_damage_amenity_daily_quantiles{_AMEN['file_tag']}.csv"
        damage_amenity = pd.read_csv(os.path.join(out_folder, amenity_file), header=0)
        damage_prod = pd.read_csv(os.path.join(out_folder, 'MC_damage_prod_daily_quantiles.csv'), header=0)
        print(f"sc2 amenity DRF file: {amenity_file}")

        damage_amenity_diff = damage_function_exogenous_daily(tmean, damage_amenity, T, quantile=QUANTILE_AMENITY)
        damage_productivity_diff = damage_function_exogenous_daily(tmean, damage_prod, T, quantile=QUANTILE_PROD)

        print(f"damage_prod_diff: min={damage_productivity_diff.min():.4f}, "f"max={damage_productivity_diff.max():.4f}, "f"mean={damage_productivity_diff.mean():.4f}")
# Expected range: roughly 0.90 to 1.05

    elif SCENARIO == 3:
        # tas_ssp2_annual.csv is in LONG format: 3 cols (id, year, temp).
        # Pivot to wide (id index, year columns) before passing to damage_function_exogenous_annual,
        # which expects a wide-format DataFrame indexed by id with year columns.
        tmean = pd.read_csv(os.path.join(data_folder, 'tas_ssp2_annual.csv'), header=0)
        tmean = tmean[tmean['id'].isin(indices_drop_id.tolist())]

        # Auto-detect year column (numeric values) and value column (the remaining one)
        non_id_cols = [c for c in tmean.columns if c != 'id']
        # year column = the one whose values look like years (4-digit ints)
        year_col, value_col = None, None
        for c in non_id_cols:
            sample = tmean[c].dropna().iloc[0] if len(tmean[c].dropna()) > 0 else None
            try:
                if sample is not None and 1900 <= int(sample) <= 2200:
                    year_col = c
                    break
            except (ValueError, TypeError):
                pass
        if year_col is None:
            # fallback: assume conventional naming
            year_col  = 'year' if 'year' in non_id_cols else non_id_cols[0]
        value_col = [c for c in non_id_cols if c != year_col][0]
        print(f"sc3: long-format tas_ssp2_annual.csv, year_col='{year_col}', value_col='{value_col}'")

        annual_temp = tmean.pivot(index='id', columns=year_col, values=value_col)
        annual_temp = annual_temp.sort_index(axis=1)
        print(f"sc3: pivoted annual_temp shape={annual_temp.shape}, year cols={list(annual_temp.columns)[:5]}...")

        _grid        = _ADRF['grid']
        amenity_file = f"MC_damage_amenity_{_grid}_quantiles{_AMEN['file_tag']}.csv"
        prod_file    = f"MC_damage_prod_{_grid}_quantiles.csv"
        damage_amenity = pd.read_csv(os.path.join(out_folder, amenity_file), header=0)
        damage_prod = pd.read_csv(os.path.join(out_folder, prod_file), header=0)
        print(f"sc3[{ANNUAL_DRF_SOURCE}] amenity DRF file: {amenity_file}")
        print(f"sc3[{ANNUAL_DRF_SOURCE}] prod    DRF file: {prod_file}")

        damage_amenity_diff = damage_function_exogenous_annual(annual_temp, damage_amenity, T, quantile=QUANTILE_AMENITY)
        damage_productivity_diff = damage_function_exogenous_annual(annual_temp, damage_prod, T,  quantile=QUANTILE_PROD)
        print(f"damage_prod_diff: min={damage_productivity_diff.min():.4f}, "f"max={damage_productivity_diff.max():.4f}, "f"mean={damage_productivity_diff.mean():.4f}")

    damage_amenity_diff      = np.nan_to_num(damage_amenity_diff,      nan=1.0)
    damage_productivity_diff = np.nan_to_num(damage_productivity_diff, nan=1.0)

    # Restrict exogenous growth paths to kept locations
    growth_pop_kept  = growth_population_exogenous[indices_drop, :]   # (14463, 36)
    growth_prod_kept = growth_productivity_exogenous[indices_drop, :]  # (14463, 36)

    run_diagnostics(population,migration0,lambda0,growth_pop_kept,growth_prod_kept,utility_hat_guess)


    print_quantiles_sparse("lambda0", lambda0)
    print_quantiles_sparse("migration0", migration0)

    # Output suffix, hoisted BEFORE the loop because DUMP_OUTER needs it.
    # The scenario label differs from SCENARIO only when the annual branch reads
    # the daily curve; otherwise the suffix is unchanged.
    _sclabel = 4 if (SCENARIO == 3 and ANNUAL_DRF_SOURCE == 'dailycurve') else SCENARIO
    suffix = f"_sc{_sclabel}{_qtag}" + _AMEN['out_tag'] + _captag + _drftag

    # -- outer loop ------------------------------------------------------------
    iter_outer = 1
    i_outer    = params['tol_outer'] + 1.0


    income = w[0,:] * population[:,0]
    inflow0 = lambda0.dot(income)
    ratio0  = inflow0 / income
    print(f"sum(income)       = {income.sum():.6e}")
    print(f"sum(lambda.dot)   = {lambda0.dot(income).sum():.6e}") 
    print(f"sum(lambda.T.dot) = {lambda0.T.dot(income).sum():.6e}")
    print(f"ratio inflow/income: median={np.median(ratio0):.3e}, min={ratio0.min():.3e}, max={ratio0.max():.3e}")
    print(f"chi*ratio median={params['param_chi'] * np.median(ratio0):.3e}")


    income0  = w[0,:] * population[:,0]
    inflow_row = lambda0.dot(income0)      # the live variant
    inflow_col = lambda0.T.dot(income0)    # transposed variant

    ratio_row = inflow_row / income0
    ratio_col = inflow_col / income0

    print(f"Row version - ratio: median={np.median(ratio_row):.3e}, max={ratio_row.max():.3e}")
    print(f"Col version - ratio: median={np.median(ratio_col):.3e}, max={ratio_col.max():.3e}")

    # If w[0] is at the lambda0 equilibrium, then lambda0.T.dot(w*L) / L is w.
    income0   = w[0,:] * population[:,0]
    w_implied = lambda0.T.dot(income0) / np.maximum(population[:,0], 1)
    w_implied = w_implied / np.mean(w_implied)

    ratio = w_implied / np.maximum(w[0,:], 1e-8)
    print(f"ratio w_implied/w[0]: median={np.median(ratio):.4f}, "
        f"min={ratio.min():.4f}, max={ratio.max():.4f}")

    worst = np.argsort(np.abs(ratio - 1))[-5:]
    for i in worst:
        print(f"  loc {i}: w={w[0,i]:.4f}, w_implied={w_implied[i]:.4f}, "
            f"ratio={ratio[i]:.2f}, pop={population[i,0]:.0f}")


    # --- local fixed point settings ------------------------------------------
    # tol_pop is a criterion on pop_{t+1}, as a root mean square relative gap.
    TOL_POP   = params.get('tol_pop', 1e-8)
    MAX_POP   = params.get('draw_pop', 60)
    DAMP_POP  = params.get('damp_pop', 0.5)
    
    amenity_hat_path = np.ones((nb_locations, T - 1))

    productivity_hat = compute_productivity_hat(damage_productivity_diff, growth_prod_kept)
    
    mu_t = csr_matrix(migration0)
    lambda_t_path = lambda0.copy()
    w[1:, :] = 1.0
    w_hat[:, :] = 1.0
    
    n_slow = 0
    t_start_all = time.time()
    
    for t in range(T - 1):
        t0 = time.time()
    
        # Initial guess for pop_{t+1}: exogenous growth alone, without
        # migration. It is the closest starting point available without
        # computing anything; starting at pop_t would be further away.
        pop_t = population[:, t].copy()
        pop_next = pop_t * growth_pop_kept[:, t]
        pop_next[pop_next <= 0] = 0.0
    
        i_pop = TOL_POP + 1.0
        k = 0
        inner_iters = 0
    
        # ---- local fixed point on pop_{t+1} ------------------------------
        while i_pop > TOL_POP and k < MAX_POP:
            amen_t = compute_amenity_hat_t(pop_t, pop_next,
                                          params['param_lambda'],
                                          damage_amenity_diff[:, t])
    
            inc_t, lam_next, it_in, err_in, conv_in = solve_temporary_equilibrium_t(
                t, w, w_hat, pop_t, pop_next, productivity_hat[:, t],
                lambda_t_path, params)
            inner_iters += it_in
    
            # A temporary equilibrium that has not converged injects noise
            # into the local fixed point. Report it rather than absorb it.
            if not conv_in:
                print(f"  [t={t} k={k}] equilibre temporaire NON converge "
                      f"(err={err_in:.2e}, {it_in} iter)")
    
            u_hat_t = amen_t * inc_t
            mu_try, pop_try = step_migration(mu_t, u_hat_t, pop_t,
                                            growth_pop_kept[:, t],
                                            params['Omega'])
    
            i_pop = np.mean(((pop_try - pop_next)
                            / np.maximum(np.abs(pop_next), 1.0)) ** 2)
            pop_next = DAMP_POP * pop_try + (1.0 - DAMP_POP) * pop_next
            k += 1
    
        if i_pop > TOL_POP:
            n_slow += 1
            print(f"  [t={t}] point fixe local NON converge apres {k} iter "
                  f"(i_pop={i_pop:.2e})")
    
        # ---- the period is closed: freeze it and move on ------------------
        # One last evaluation at the converged pop_next, so that the stored
        # quantities are those of the fixed point and not of the iteration
        # before it.
        amen_t = compute_amenity_hat_t(pop_t, pop_next, params['param_lambda'],
                                      damage_amenity_diff[:, t])
        inc_t, lam_next, it_in, err_in, conv_in = solve_temporary_equilibrium_t(
            t, w, w_hat, pop_t, pop_next, productivity_hat[:, t],
            lambda_t_path, params)
    
        income_hat[t, :] = inc_t
        amenity_hat_path[:, t] = amen_t
        lambda_t_path = lam_next
        mu_t, population[:, t + 1] = step_migration(
            mu_t, amen_t * inc_t, pop_t, growth_pop_kept[:, t], params['Omega'])
    
        # Clearing check: the value of production must equal income.
        ratio = lambda_t_path.dot(w[t + 1, :] * population[:, t + 1])
        denom = w[t + 1, :] * population[:, t + 1]
        mask = denom > 0
        ratio_check = np.median(ratio[mask] / denom[mask])
        flag = "" if abs(ratio_check - 1) <= 0.1 else "  !! APUREMENT NON SATISFAIT"
        print(f"  [t={t}] {k} iter locales | {inner_iters} iter internes | "
              f"i_pop={i_pop:.2e} | apurement={ratio_check:.4f} | "
              f"{time.time() - t0:.1f}s{flag}")
    
        # Warm start: the wage at t is the best starting point for t+1.
        if t < T - 2:
            w_hat[t + 1, :] = (params['damp_tstep'] * w_hat[t, :]
                              + (1.0 - params['damp_tstep']) * w_hat[t + 1, :])
    
    u_hat = amenity_hat_path.T * income_hat
    
    print(f"\n=== trajectoire resolue en {time.time() - t_start_all:.0f}s ===")
    # A mechanism that can do nothing must say whether it acted.
    if n_slow:
        print(f"[!] {n_slow}/{T-1} periodes non convergees a {TOL_POP:.0e}. "
              f"Les resultats en heritent : ne pas publier sans trancher.")
    else:
        print(f"[ok] les {T-1} periodes convergees a {TOL_POP:.0e}.")
    
    ## NO INDEX COLUMN. The emitters call pd.read_csv(f).values, so an index
    ## column would be read as data and shift the whole matrix. The expected
    ## shape of welfare is (16, N): 16 TRANSITIONS between the 17 periods, not
    ## 17 levels.
    os.makedirs(out_folder, exist_ok=True)

    pd.DataFrame(u_hat).to_csv(
        os.path.join(out_folder, f'welfare{suffix}.csv'), index=False)
    pd.DataFrame(amenity_hat_path.T).to_csv(
        os.path.join(out_folder, f'amenity_hat_path{suffix}.csv'), index=False)
    pd.DataFrame(productivity_hat.T).to_csv(
        os.path.join(out_folder, f'productivity_hat_path{suffix}.csv'), index=False)
    print(f"[write] welfare{suffix}.csv  {u_hat.shape}")
    print(f"[write] amenity_hat_path{suffix}.csv  {amenity_hat_path.T.shape}")
    print(f"[write] productivity_hat_path{suffix}.csv  {productivity_hat.T.shape}")

if __name__ == "__main__":
    main()

