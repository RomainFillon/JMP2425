# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from scipy.interpolate import interp1d

# =============================================================================
# forward.py -- the per-period equilibrium blocks of the solver.
#
# Holds two families of functions:
#   A. the per-period evaluations the local fixed point needs;
#   B. update_lambda, price_index_hat and their shared weight, where most of
#      the solver time is spent.
# =============================================================================

import numpy as np
from scipy.sparse import csr_matrix


# =============================================================================
# A. PER-PERIOD VERSIONS
#
# These evaluate ONE period, taking pop_{t+1} as an argument rather than reading
# it from a global trajectory array. pop_{t+1} is the unknown of the local fixed
# point, so it has to be a parameter.
# =============================================================================

def compute_amenity_hat_t(pop_t, pop_tplus, param_lambda, damage_amenity_t):
    """Amenite de la periode t.

    ATTENTION -- C'EST ICI QUE VIT LA CIRCULARITE. L'amenite depend du RATIO
    pop_{t+1}/pop_t : la congestion fait baisser la valeur d'un lieu quand sa
    population croit. Or pop_{t+1} decoule de u_hat_t, qui decoule de l'amenite.
    Le modele n'est donc PAS sequentiel : il y a un point fixe a resoudre a
    l'interieur de chaque periode. Ce qui est evitable, c'est de le resoudre sur
    toute la trajectoire a la fois.
    """
    ratio = np.ones_like(pop_t, dtype=float)
    mask = (pop_t > 0) & (pop_tplus > 0)
    ratio[mask] = pop_tplus[mask] / pop_t[mask]
    return (ratio ** (-param_lambda)) * damage_amenity_t


def unit_price_bundle_hat_t(w_hat_t, pop_t, pop_tplus, param_mu):
    """Meme remarque : depend du ratio de populations, donc de l'inconnue."""
    ratio = np.ones_like(pop_t, dtype=float)
    mask = (pop_t > 0) & (pop_tplus > 0)
    ratio[mask] = pop_tplus[mask] / pop_t[mask]
    return w_hat_t * (ratio ** (1.0 - param_mu))


def step_migration(mu_t, u_hat_t, pop_t, growth_pop_t, Omega):
    """Un pas de migration myope : mu_t ne depend que de u_hat_t.

    Rend (mu_{t+1}, pop_{t+1}). C'est la moitie sequentielle du modele -- aucune
    anticipation, aucun terme en t+1 du cote migration.
    """
    u_pow = u_hat_t ** (1.0 / Omega)
    mu_new = mu_t.multiply(u_pow)

    row_sums = np.asarray(mu_new.sum(axis=1)).ravel()
    inv = np.ones_like(row_sums)
    nz = row_sums != 0
    inv[nz] = 1.0 / row_sums[nz]
    mu_new = mu_new.multiply(inv[:, None]).tocsr()

    if growth_pop_t.mean() != 1:
        pop_tplus = pop_t * growth_pop_t
        flows = mu_new.multiply(pop_t[:, None])
        net = (np.asarray(flows.sum(axis=0)).ravel()
               - np.asarray(flows.sum(axis=1)).ravel())
        pop_tplus = pop_tplus + net
    else:
        pop_tplus = pop_t.copy()

    pop_tplus[pop_tplus <= 0] = 0.0
    return mu_new, pop_tplus


# =============================================================================
# B. THE HOT PATH
#
# update_lambda is the most expensive function of the solver: it multiplies a
# 17048 x 17048 sparse matrix on every inner iteration, and each multiply
# ALLOCATES a new matrix, since scipy does not scale in place. The shared weight
# below removes one of those multiplies, which is the same vector computed twice
# per iteration.
# =============================================================================

def trade_weights(productivity_hat_t, x_hat, param_theta):
    """(Z/x)^theta, calcule UNE FOIS et partage.

    price_index_hat calculait (x/Z)^(-theta) et update_lambda (Z/x)^theta :
    c'est le meme vecteur, calcule deux fois par iteration.
    """
    x_safe = np.where(x_hat > 1e-10, x_hat, 1e-10)
    z_safe = np.where(productivity_hat_t > 1e-10, productivity_hat_t, 1e-10)
    return np.clip((z_safe / x_safe) ** param_theta, 1e-10, 1e10)


def price_index_hat_fast(weights, lambda_t, param_theta):
    """P_hat = (lambda . weights)^(-1/theta), sur des poids deja calcules."""
    P_hat = lambda_t.dot(weights)
    P_hat = np.nan_to_num(P_hat, nan=1.0, posinf=1.0, neginf=1.0)
    P_hat = np.where(P_hat > 0, P_hat, 1.0)
    return np.power(P_hat, -1.0 / param_theta)


def update_lambda_fast(lambda_t, weights):
    """lambda' = normalise_lignes(lambda . diag(weights)).

    DEUX ECONOMIES, aucune approximation.

    1. LE DENOMINATEUR P^theta EST SUPPRIME. L'original faisait
         lambda * num[colonnes] / denom[lignes],  puis normalisait les lignes.
       Or denom est CONSTANT PAR LIGNE : la normalisation qui suit l'elimine
       exactement. Le calculer, l'appliquer, puis le diviser a nouveau etait
       une passe complete sur 17048^2 pour un resultat inchange au bit pres.

    2. LA NORMALISATION SE FAIT EN PLACE sur L.data. L'original appelait un
       troisieme multiply, donc une troisieme allocation. Ici on divise le
       tableau de donnees par les sommes de lignes repetees selon indptr --
       aucune matrice intermediaire.

    Passe de trois allocations a une.
    """
    L = lambda_t.multiply(weights[np.newaxis, :]).tocsr()
    row_sums = np.asarray(L.sum(axis=1)).ravel()

    zero_rows = row_sums == 0
    if zero_rows.any():
        # An empty row means a destination has NO source left. That is not an
        # edge case: it signals that the weights have blown up.
        print(f"  warning: {zero_rows.sum()} lignes nulles dans lambda "
              f"-> reset diagonale")

    safe = np.where(row_sums > 0, row_sums, 1.0)
    L.data /= np.repeat(safe, np.diff(L.indptr))

    if zero_rows.any():
        from scipy.sparse import diags
        L = (L + diags(zero_rows.astype(float))).tocsr()
    return L


# =============================================================================
# C. TEMPORARY EQUILIBRIUM, with pop_{t+1} as an argument
#
# pop_{t+1} is the unknown of the local fixed point, so it is passed in rather
# than read from a global trajectory array.
# =============================================================================

def solve_temporary_equilibrium_t(t, w, w_hat, pop_t, pop_tplus,
                                  productivity_hat_t, lambda_t, params):
    tol_inner  = params['tol_inner']
    draw_inner = params['draw_inner']
    param_mu   = params['param_mu']
    param_theta = params['param_theta']
    param_chi  = params['param_chi']

    autarkic = np.asarray(lambda_t.sum(axis=1)).ravel() <= 0

    i_inner = tol_inner + 1.0
    iter_inner = 1
    lambda_tplus = lambda_t
    x_hat = None

    while i_inner > tol_inner and iter_inner <= draw_inner:
        x_hat = unit_price_bundle_hat_t(w_hat[t, :], pop_t, pop_tplus, param_mu)
        x_hat = np.clip(x_hat, 1e-4, 1e4)

        # One weight computation shared by both uses.
        weights = trade_weights(productivity_hat_t, x_hat, param_theta)
        P_hat = price_index_hat_fast(weights, lambda_t, param_theta)
        lambda_tplus = update_lambda_fast(lambda_t, weights)

        w_hat_new = equilibrium(w[t, :], w_hat[t, :], pop_t,
                                lambda_tplus, param_chi)
        w_hat_new[autarkic] = 1.0

        i_inner = np.mean(((w_hat_new - w_hat[t, :])
                           / np.maximum(w_hat[t, :], 1e-8)) ** 2)
        if i_inner < tol_inner:
            break

        # There is deliberately no early exit on a loose error threshold. A mean
        # relative squared error of the order of 0.05 is large: the loop would
        # leave NOT CONVERGED and the outer fixed point would inherit that
        # noise. An exit threshold that hides non-convergence is worse than a
        # slow loop.

        if i_inner > 1e4:
            damp = 0.05
        elif i_inner > 1e2:
            damp = 0.1
        elif i_inner > 1.0:
            damp = 0.3
        else:
            damp = 0.5

        w_hat[t, :] = damp * w_hat_new + (1.0 - damp) * w_hat[t, :]
        iter_inner += 1

    converged = i_inner <= tol_inner

    w_hat[t, :] = np.clip(w_hat[t, :], 0.1, 10)
    w_hat[t, :] = w_hat[t, :] / np.mean(w_hat[t, :])
    w[t + 1, :] = np.maximum(w_hat[t, :] * w[t, :], 1e-6)
    w[t + 1, :] = w[t + 1, :] / np.mean(w[t + 1, :])

    weights = trade_weights(productivity_hat_t, x_hat, param_theta)
    P_hat = price_index_hat_fast(weights, lambda_t, param_theta)
    income_hat = np.ones_like(pop_t, dtype=float)
    mask = P_hat > 0
    income_hat[mask] = w_hat[t, mask] / P_hat[mask]

    return income_hat, lambda_tplus, iter_inner, i_inner, converged

# =========================
# MIGRATION
# =========================
def compute_migrations_populations(utility_hat_guess, m0, population,
                                   growth_population_exogenous, T, Omega, nb_locations):

    mu_t = csr_matrix(m0)

    for t in range(T - 1):

        # avoids the diagonal
        u_t_hat = utility_hat_guess[t] ** (1.0 / Omega)
        mu_t = mu_t.multiply(u_t_hat)

        # row normalisation
        row_sums = np.array(mu_t.sum(axis=1)).ravel()
        inv = np.ones_like(row_sums)
        mask = row_sums != 0
        inv[mask] = 1.0 / row_sums[mask]
        mu_t = mu_t.multiply(inv[:, None])

        if growth_population_exogenous[:, t].mean() != 1:
            population[:, t + 1] = population[:, t] * growth_population_exogenous[:, t]

            flows = mu_t.multiply(population[:, t][:, None])
            net = np.array(flows.sum(axis=0)).ravel() - np.array(flows.sum(axis=1)).ravel()
            population[:, t + 1] += net
        else:
            population[:, t + 1] = population[:, t]

        population[population[:, t + 1] <= 0, t + 1] = 0

    return population


# =========================
# AMENITIES / PRODUCTIVITY
# =========================
def compute_amenity_hat(population, param_lambda, damage_amenity_diff):

    ratio = np.ones_like(damage_amenity_diff)

    mask = (population[:, :-1] > 0) & (population[:, 1:] > 0)
    ratio[mask] = population[:, 1:][mask] / population[:, :-1][mask]

    return (ratio ** (-param_lambda)) * damage_amenity_diff


def compute_productivity_hat(damage_productivity_diff, growth_productivity_exogenous):
    return growth_productivity_exogenous * damage_productivity_diff


def compute_annual_temperature(tmean):
    df = tmean.copy()

    # year columns
    year_cols = [c for c in df.columns if c.isdigit()]

    # numerator: value * frequency
    for col in year_cols:
        df[col] = df['value'] * df[col]

    # sum by id
    sum_weighted = df.groupby('id')[year_cols].sum()

    # denominator: sum of the frequencies
    freq = tmean.groupby('id')[year_cols].sum()

    # mean temperature
    annual_temp = sum_weighted / freq

    return annual_temp

def damage_function_exogenous_daily(daily_path, damage_consolidated, T,
                                     quantile='q50'):
    """
    Vraie integration daily : pour chaque (id, year),
        f_year = sum_b omega_b * f_daily(T_b)
    daily_path : DataFrame avec colonnes
        - 'id'    : location id
        - 'value' : T du bin (degres C, entier)
        - colonnes annees '2015', '2020', ... contenant la frequence omega_b
    """
    assert quantile in damage_consolidated.columns

    # 1. interpolate f on the temperature bins
    lookup = damage_consolidated[['temperature_C', quantile]].dropna()
    interp = interp1d(lookup['temperature_C'], lookup[quantile],
                      bounds_error=False, fill_value='extrapolate')

    df    = daily_path.copy()
    f_bin = interp(df['value'].to_numpy())  # f(T_b) for each bin

    # 2. year columns, sorted numerically
    year_cols = sorted([c for c in df.columns if str(c).isdigit()],
                       key=lambda c: int(c))

    # NaN -> 0 in the frequencies
    df[year_cols] = df[year_cols].fillna(0.0)

    # 3. sum_b omega_b * f(T_b) by (id, year)
    weighted = df[year_cols].mul(f_bin, axis=0)
    sum_f    = weighted.groupby(df['id']).sum()

    # 4. normalise by sum_b omega_b, giving E_b[f(T)]
    n_b    = df.groupby('id')[year_cols].sum().replace(0.0, np.nan)
    f_year = sum_f / n_b
    f_year = f_year.ffill(axis=1).bfill(axis=1)

    # 5. diff inter-periodes + exp -> facteur multiplicatif
    f_year = f_year.sort_index(axis=1)
    out    = f_year.diff(axis=1).apply(np.exp).iloc[:, 1:]

    arr   = out.to_numpy()
    n_out = T - 1
    if arr.shape[1] > n_out:
        print(f"damage_function_exogenous_daily: truncating from {arr.shape[1]} to {n_out} cols")
        arr = arr[:, :n_out]
    pad   = np.ones((arr.shape[0], max(0, n_out - arr.shape[1])))
    if pad.shape[1] > 0:
        print(f"damage_function_exogenous_daily: padding {pad.shape[1]} cols with 1s "
              f"(after timestep {arr.shape[1]-pad.shape[1]} climate damage is frozen)")
    return np.hstack((arr, pad))

def damage_function_exogenous_daily_past(temperature_path, damage_consolidated, T,
                             quantile='q50'):

    assert quantile in damage_consolidated.columns, \
        f"Colonne '{quantile}' absente â€” colonnes disponibles : {list(damage_consolidated.columns)}"

    # 1. Interpolation
    damage_lookup = damage_consolidated[['temperature_C', quantile]].dropna()
    interp = interp1d(damage_lookup['temperature_C'],
                      damage_lookup[quantile],
                      bounds_error=False,
                      fill_value='extrapolate')

    temperature_path = temperature_path.copy()
    temperature_path['_damage'] = interp(temperature_path['temp_mean'])

    # 2. apply the damage

    temperature_path['damaged'] = temperature_path['_damage']  # f(T) only, not T*f(T)


    # 3. Pivot
    merged = temperature_path.pivot(index='id', columns='five_year', values='damaged')
    merged = merged.sort_index(axis=1)

    # 4. diff + exp
    merged = merged.diff(axis=1).map(np.exp)

    merged = merged.iloc[:, 1:]
    arr = merged.to_numpy()

    n_out = T - 1
    pad = np.ones((arr.shape[0], max(0, n_out - arr.shape[1])))

    return np.hstack((arr, pad))


def damage_function_exogenous_annual(annual_temp, damage_consolidated,
                                      T, quantile='q50'):
    assert quantile in damage_consolidated.columns

    lookup = damage_consolidated[['temperature_C', quantile]].dropna()
    interp = interp1d(lookup['temperature_C'], lookup[quantile],
                      bounds_error=False, fill_value='extrapolate')

    damaged = annual_temp.copy()
    for col in damaged.columns:
        damaged[col] = interp(damaged[col].to_numpy())  # f(T), not T*f(T)

    damaged = damaged.sort_index(axis=1)
    damaged = damaged.diff(axis=1).apply(np.exp)
    damaged = damaged.iloc[:, 1:]

    arr   = damaged.to_numpy()
    n_out = T - 1
    if arr.shape[1] > n_out:
        print(f"damage_function_exogenous_annual: truncating from {arr.shape[1]} to {n_out} cols")
        arr = arr[:, :n_out]
    pad   = np.ones((arr.shape[0], max(0, n_out - arr.shape[1])))
    if pad.shape[1] > 0:
        print(f"damage_function_exogenous_annual: padding {pad.shape[1]} cols with 1s "
              f"(after timestep {arr.shape[1]-pad.shape[1]} climate damage is frozen)")
    return np.hstack((arr, pad))


def damage_function_exogenous_annual_past(annual_temp, damage_consolidated,
                                     T, quantile='q50'):

    assert quantile in damage_consolidated.columns, \
        f"Colonne '{quantile}' absente â€” colonnes disponibles : {list(damage_consolidated.columns)}"

    # 1. Interpolation
    damage_lookup = damage_consolidated[['temperature_C', quantile]].dropna()
    interp = interp1d(damage_lookup['temperature_C'],
                      damage_lookup[quantile],
                      bounds_error=False,
                      fill_value='extrapolate')

    # 2. apply directly on the annual temperature
    damaged = annual_temp.copy()
    for col in damaged.columns:
        damaged[col] = interp(damaged[col].to_numpy())  # f(T) only, not T*f(T)

    # 3. diff + exp, as above
    damaged = damaged.sort_index(axis=1)
    damaged = damaged.diff(axis=1).apply(np.exp)

    damaged = damaged.iloc[:, 1:]
    arr = damaged.to_numpy()

    # padding
    n_out = T - 1
    pad = np.ones((arr.shape[0], max(0, n_out - arr.shape[1])))

    return np.hstack((arr, pad))


# =========================
# CORE ECON
# =========================
def unit_price_bundle_hat(t, w_hat_t, population, param_mu):

    ratio = np.ones(population.shape[0])
    mask = (population[:, t] > 0) & (population[:, t + 1] > 0)

    ratio[mask] = population[mask, t + 1] / population[mask, t]

    return w_hat_t * (ratio ** (1.0 - param_mu))

def price_index_hat(productivity_hat_t, x_hat, lambda_t, param_theta):
    # guard x_hat
    x_hat_safe = np.where(x_hat > 1e-10, x_hat, 1e-10)
    prod_safe   = np.where(productivity_hat_t > 1e-10, productivity_hat_t, 1e-10)
    
    effective_cost = x_hat_safe / prod_safe
    weights = effective_cost ** (-param_theta)
    
    # clip before the dot product, to avoid an overflow
    weights = np.clip(weights, 0, 1e10)
    
    P_hat = lambda_t.dot(weights)
    
    # guard against NaN AND zero
    P_hat = np.nan_to_num(P_hat, nan=1.0, posinf=1.0, neginf=1.0)
    P_hat = np.where(P_hat > 0, P_hat, 1.0)
    
    return np.power(P_hat, -1.0 / param_theta)


def update_lambda(lambda_t, productivity_hat_t, P_hat, x_hat, param_theta):
    
    # guard the inputs
    x_hat_safe  = np.where(x_hat > 1e-10, x_hat, 1e-10)
    prod_safe   = np.where(productivity_hat_t > 1e-10, productivity_hat_t, 1e-10)
    P_hat_safe  = np.where(P_hat > 1e-10, P_hat, 1e-10)
    
    # numerator: (Z(s) / x(s))^theta, scaling the columns (sources s)
    num = (prod_safe / x_hat_safe) ** param_theta
    num = np.clip(num, 1e-10, 1e10)
    
    # denominator: P(r)^theta, scaling the rows (destinations r)
    denom = P_hat_safe ** param_theta
    denom = np.clip(denom, 1e-10, 1e10)
    
    # update
    lambda_t = lambda_t.multiply(num[np.newaxis, :])   # colonnes = sources
    lambda_t = lambda_t.multiply(1.0 / denom[:, np.newaxis])  # lignes = destinations
    
    # renormalisation
    row_sums = np.array(lambda_t.sum(axis=1)).ravel()
    
    # detect the empty rows
    zero_rows = row_sums == 0
    if zero_rows.sum() > 0:
        print(f"  warning: {zero_rows.sum()} lignes nulles dans lambda ? reset diagonale")
    
    inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
    lambda_t = lambda_t.multiply(inv[:, np.newaxis])
    
    # reset the diagonal on the empty rows
    if zero_rows.sum() > 0:
        from scipy.sparse import diags
        lambda_t = lambda_t + diags(zero_rows.astype(float))
    
    return lambda_t


def equilibrium(w_t, w_hat_t, population_t, lambda_tplus, param_chi):
    w_tplus  = w_hat_t * w_t
    income   = w_tplus * population_t
    

#    weighted = param_chi * lambda_tplus.T.dot(income)

    weighted = param_chi * lambda_tplus.T.dot(income)
    
    w_new = np.ones_like(w_tplus)
    mask  = population_t > 0
    w_new[mask] = weighted[mask] / population_t[mask]
    
    w_hat_new = np.ones_like(w_hat_t)
    mask_w = w_t > 1e-8
    w_hat_new[mask_w] = w_new[mask_w] / w_t[mask_w]
    w_hat_new[~mask_w] = 1.0  # force to 1 when w is too small
    w_hat_new = np.clip(w_hat_new, 0.1, 10)
    w_hat_new = w_hat_new / np.mean(w_hat_new)

    return w_hat_new

# =========================
# SOLVER
# =========================

def solve_temporary_equilibrium(t, w, w_hat, population, productivity_hat,
                                 lambda_t, params):

    param_mu      = params['param_mu']
    param_theta   = params['param_theta']
    param_chi     = params['param_chi']
    draw_inner    = params['draw_inner']
    tol_inner     = params['tol_inner']
    param_epsilon = params['epsilon']

    # Warm start entre periodes
    if t == 0:
        w_hat[t, :] = np.ones(population.shape[0])
    else:
        w_hat[t, :] = w_hat[t - 1, :]

    # Precomputations, hoisted out of the while loop
    lambda_t_csr = lambda_t.tocsr()
    diag         = np.array(lambda_t_csr.diagonal())
    autarkic     = diag > 0.99

    i_inner    = tol_inner + 1.0
    iter_inner = 1

    while i_inner > tol_inner and iter_inner <= draw_inner:

        x_hat = unit_price_bundle_hat(t, w_hat[t, :], population, param_mu)
        x_hat = np.clip(x_hat, 1e-4, 1e4)

        P_hat = price_index_hat(productivity_hat[:, t], x_hat, lambda_t, param_theta)

        lambda_tplus = update_lambda(lambda_t, productivity_hat[:, t], P_hat, x_hat, param_theta)

        w_hat_new = equilibrium(w[t, :], w_hat[t, :], population[:, t], lambda_tplus, param_chi)
        w_hat_new[autarkic] = 1.0

        i_inner = np.mean(((w_hat_new - w_hat[t, :]) / np.maximum(w_hat[t, :], 1e-8)) ** 2)

        # Critere arret
        if i_inner < 1e-6:
            break
        if iter_inner > 50 and i_inner < 0.05:
            break

        # Damping adaptatif
        if i_inner > 1e4:
            damp = 0.05
        elif i_inner > 1e2:
            damp = 0.1
        elif i_inner > 1.0:
            damp = 0.3
        else:
            damp = 0.5

        w_hat[t, :] = damp * w_hat_new + (1.0 - damp) * w_hat[t, :]
        iter_inner += 1

    # Print minimal
    if t % 10 == 0 or t == 0:
        print(f"[t={t}] iter={iter_inner} | error={i_inner:.3e}")

    # Update economy
    w_hat[t, :] = np.clip(w_hat[t, :], 0.1, 10)
    w_hat[t, :] = w_hat[t, :] / np.mean(w_hat[t, :])
    w[t + 1, :] = w_hat[t, :] * w[t, :]
    w[t + 1, :] = np.maximum(w[t + 1, :], 1e-6)
    w[t + 1, :] = w[t + 1, :] / np.mean(w[t + 1, :])

    P_hat = price_index_hat(productivity_hat[:, t], x_hat, lambda_t, param_theta)

    income_hat = np.ones(population.shape[0])
    mask = P_hat > 0
    income_hat[mask] = w_hat[t, mask] / P_hat[mask]

    return income_hat, lambda_tplus


def solve_temporary_equilibrium_past(t, w, w_hat, population, productivity_hat,
                                 lambda_t, params):

    param_mu      = params['param_mu']
    param_theta   = params['param_theta']
    param_chi     = params['param_chi']
    draw_inner    = params['draw_inner']
    tol_inner     = params['tol_inner']
    damp_initial  = params['damp_inner_initial']
    damp_final    = params['damp_inner_final']
    param_epsilon = params['epsilon']

    # =========================
    # WARM START OR NOT

    if t > 0:
        w_hat[t, :] = w_hat[t - 1, :]

    i_inner = tol_inner + 1.0
    iter_inner = 1

    # opening print
    print(f"[t={t}] start inner loop")

  # before the while loop of solve_temporary_equilibrium
    x_hat_init = unit_price_bundle_hat(t, w_hat[t,:], population, param_mu)
    print(f"x_hat: min={x_hat_init.min():.3e}, max={x_hat_init.max():.3e}, NaN={np.isnan(x_hat_init).sum()}")
    weights_init = (x_hat_init / productivity_hat[:,t]) ** (-param_theta)


    income0  = w[0,:] * population[:,0]
    # test both conventions
    recv_row = lambda_t.dot(income0)    # if lambda[r,s] is the purchase share of r
    recv_col = lambda_t.T.dot(income0)  # if lambda[r,s] is the sales share of r

    if t == 0:
        print(f"w[0,:] median = {np.median(w[0,:]):.3e}")
        print(f"w[0,:] min = {w[0,:].min():.3e}, max = {w[0,:].max():.3e}")
        print(f"population[:,0] median = {np.median(population[:,0]):.3e}")
        print(f"population[:,0] min = {np.min(population[:,0]):.3e}")
        recv_row = lambda_t.dot(income0)
        recv_col = lambda_t.T.dot(income0)
        print(f"Row: recv/income mediane = {np.median(recv_row/income0):.4f}")
        print(f"Col: recv/income mediane = {np.median(recv_col/income0):.4f}")


    lambda_t_csr = lambda_t.tocsr() 
    diag = np.array(lambda_t_csr.diagonal())
    autarkic = diag > 0.99

    x_hat = unit_price_bundle_hat(t, w_hat[t, :], population, param_mu)

    while i_inner > tol_inner and iter_inner <= draw_inner:



        x_hat = unit_price_bundle_hat(t, w_hat[t, :], population, param_mu)
        x_hat = np.clip(x_hat, 1e-4, 1e4)   # ? ajouter ce clip

        P_hat = price_index_hat(productivity_hat[:, t], x_hat, lambda_t, param_theta)

        if iter_inner <= 3:
            print(f"  P_hat: min={P_hat.min():.3e}, max={P_hat.max():.3e}, NaN={np.isnan(P_hat).sum()}, Inf={np.isinf(P_hat).sum()}")

        lambda_tplus = update_lambda(lambda_t, productivity_hat[:, t], P_hat, x_hat, param_theta)

        if iter_inner <= 3:
            rs = np.array(lambda_tplus.sum(axis=1)).ravel()
            print(f"  lambda_tplus row sums: min={rs.min():.3e}, max={rs.max():.3e}, NaN={np.isnan(rs).sum()}")


        # POPULATION AT t, NOT t+1. The wage block clears on the population in
        # place at the start of the period; taking t+1 would let this period
        # migration feed back into this period wage.
        w_hat_new = equilibrium(w[t, :], w_hat[t, :], population[:, t],lambda_tplus, param_chi)

        # An autarkic cell trades with nobody, so w_hat is not determined by
        # trade and is pinned at 1.

        w_hat_new[autarkic] = 1.0

        if iter_inner <= 3:
            print(f"  w_hat_new: min={w_hat_new.min():.3e}, max={w_hat_new.max():.3e}, NaN={np.isnan(w_hat_new).sum()}, Inf={np.isinf(w_hat_new).sum()}")


        # RELATIVE criterion, not an absolute sum of squares: with wages
        # spanning several orders of magnitude an absolute residual would be
        # dominated by the largest cells alone.
        i_inner = np.mean(((w_hat_new - w_hat[t, :]) / np.maximum(w_hat[t, :], 1e-8)) ** 2)


        # ---- diagnostics ----

        # Progress trace.
        if iter_inner % 100 == 0:
            print(f"[t={t}] iter={iter_inner} | error={i_inner:.3e}")

        # Convergence proper.
        if i_inner < 1e-6:
            print(f"[t={t}] ? convergence fine: {i_inner:.3e} (iter {iter_inner})")


        # SECOND, LOOSER ACCEPTANCE. Past fifty iterations an error of 1e-4 is
        # accepted and the loop breaks. This is an escape hatch, not the
        # tolerance: a run that leaves through here has NOT reached 1e-6, and
        # the residual it stops at is reported as a calibration diagnostic.
        if iter_inner > 50 and i_inner < 1e-4:
            print(f"[t={t}] convergence acceptable: {i_inner:.3e}")
            break

        if i_inner > 1e4:
            damp = 0.05
        elif i_inner > 1e2:
            damp = 0.1
        elif i_inner > 1.0:
            damp = 0.3
        else:
            damp = 0.5

        w_hat[t, :] = damp * w_hat_new + (1.0 - damp) * w_hat[t, :]

        if w_hat[t,:].max() > 100:
            print(f"ALERTE: w_hat max={w_hat[t,:].max():.1f}: localisation du pb:")
            worst = np.argmax(w_hat[t,:])
            print(f"ALERTE: w_hat max={w_hat[t,:].max():.1f}: loc={worst}, "f"pop={population[worst,t]:.0f}, w={w[t,worst]:.4f}")

        iter_inner += 1

    print(f"[t={t}] end inner loop | final error={i_inner:.3e} | iter={iter_inner}")

    # =========================
    # UPDATE ECONOMY
    # =========================

    w_hat[t, :] = np.clip(w_hat[t, :], 0.1, 10)
    w_hat[t, :] = w_hat[t, :] / np.mean(w_hat[t, :])
    w[t + 1, :] = w_hat[t, :] * w[t, :]
    w[t + 1, :] = np.maximum(w[t + 1, :], 1e-6)
    w[t + 1, :] = w[t + 1, :] / np.mean(w[t + 1, :])

    y_t = population[:, t] * w[t, :]
    y_tp1 = population[:, t + 1] * w[t + 1, :]

    P_hat = price_index_hat(productivity_hat[:, t], x_hat, lambda_t, param_theta)

    income_hat = np.ones(population.shape[0])

    mask = P_hat > 0
    income_hat[mask] = w_hat[t, mask] / P_hat[mask]


    print(f"w_hat[t,:]: min={w_hat[t,:].min():.3e}, max={w_hat[t,:].max():.3e}")
    print(f"P_hat_final: min={P_hat.min():.3e}, max={P_hat.max():.3e}")
    print(f"income_hat: min={income_hat.min():.3e}, max={income_hat.max():.3e}")

    return income_hat, lambda_tplus

