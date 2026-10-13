# -*- coding: utf-8 -*-
import os
import csv
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, diags
import pyreadr


def load_study_parameters(run_folder):
    params = {}
    with open(os.path.join(run_folder, 'specification_parameters.csv'), 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            params[row[0]] = int(row[1])
    with open(os.path.join(run_folder, 'specification_study.csv'), 'r') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            params[row[0]] = float(row[1])
    return params


def _load_csv_ordered(filepath, id_col='id', value_cols=None):
    """
    Load a CSV, sort it on id_col, and return (ids, values).
    value_cols: columns to extract; if None, every column except id_col.
    """
    df = pd.read_csv(filepath)
    assert id_col in df.columns, f"Colonne '{id_col}' absente dans {filepath}"
    df = df.sort_values(id_col).reset_index(drop=True)
    ids = df[id_col].values
    if value_cols is None:
        value_cols = [c for c in df.columns if c != id_col]
    return ids, df[value_cols].values


# Checks that the id VECTORS match, not that the ids denote the same cells.
# Two vintages of identifier_xy_final.csv both number 1..17048, so this passes
# even when id 1 denotes a different cell in each. Cabling it as a guard against
# misalignment would give a green light it cannot earn.
def check_id_alignment(ids_reference, **others):
    """
    Check that every id vector equals ids_reference, element by element.
    Call: check_id_alignment(ids_ref, wage0=ids_w, growth_pop=ids_gp)
    """
    for name, ids in others.items():
        if not np.array_equal(ids_reference, ids):
            raise ValueError(
                f"IDs désalignés : population0 vs {name}\n"
                f"  ex. pop[0]={ids_reference[0]}, {name}[0]={ids[0]}"
            )
    print("✓ Alignement IDs vérifié pour tous les fichiers")


def load_population(data_folder, T, indices_drop_id):
    ids, vals = _load_csv_ordered(
        os.path.join(data_folder, 'population0.csv'),
        value_cols=['value'])
    
    # Filter by ID explicitly, in the right order
    id_to_val = dict(zip(ids, vals[:, 0]))
    
    # Check all IDs are present
    missing = [i for i in indices_drop_id if i not in id_to_val]
    print(f"Missing IDs in population0.csv: {len(missing)}")
    
    pop_filtered = np.array([id_to_val.get(i, 0.0) for i in indices_drop_id])
    
    population = np.empty((len(indices_drop_id), T))
    population[:, 0] = pop_filtered
    
    print(f"population t=0: min={population[:,0].min():.3e}, max={population[:,0].max():.3e}, zeros={(population[:,0]==0).sum()}")
    
    return population


def load_growth_paths(data_folder, T):
    growth_pop  = pd.read_csv(os.path.join(data_folder, 'data_final_population.csv'))
    growth_prod = pd.read_csv(os.path.join(data_folder, 'data_final_productivity.csv'))

    def _pad(arr, T, name=''):
        # Drop a leading 'id' column if present (defensive)
        if 'id' in arr.columns:
            arr = arr.drop(columns=['id'])
        arr = arr.to_numpy()
        n_out = T - 1
        if arr.shape[1] > n_out:
            print(f"load_growth_paths[{name}]: truncating from {arr.shape[1]} to {n_out} cols")
            arr = arr[:, :n_out]
        pad = np.ones((arr.shape[0], max(0, n_out - arr.shape[1])))
        if pad.shape[1] > 0:
            print(f"load_growth_paths[{name}]: padding {pad.shape[1]} cols with 1s")
        return np.hstack((arr, pad))

    return _pad(growth_pop, T, 'pop'), _pad(growth_prod, T, 'prod')

def data_storage(T, nb_locations, data_folder, indices_drop_id):
    ids_w, vals = _load_csv_ordered(
            os.path.join(data_folder, 'wage0.csv'),
            value_cols=['value'])
    
    # Filter to kept locations in the right order
    ids_w = ids_w.flatten()
    id_to_val = dict(zip(ids_w, vals[:, 0]))
    wage_filtered = np.array([id_to_val[i] for i in indices_drop_id])
    
    w      = np.ones((T, nb_locations))
    w[0,:] = wage_filtered / 1000.0
    w[0, :] = w[0, :] / np.mean(w[0, :])#normalize wage


    print(f"w[0,:] after filter: min={w[0,:].min():.3e}, max={w[0,:].max():.3e}, zeros={(w[0,:]==0).sum()}")
    
    utility_hat_guess = np.ones((T - 1, nb_locations))
    w_hat             = np.ones((T - 1, nb_locations))
    income_hat        = np.ones((T - 1, nb_locations))
    return utility_hat_guess, w_hat, w, income_hat

def load_h5(file_path, indices_drop, name, cap_colmass=None):
    with h5py.File(file_path, 'r') as f:
        data = f[name][:].T   # R (hdf5r) writes row-stochastic in column-major
                              # order, so h5py reads it transposed.
    data = np.nan_to_num(data, nan=0.0)
    data = csr_matrix(data)
    
    # 1. filter FIRST
    data = data[indices_drop, :][:, indices_drop]
    
    # 2. drop the very small values AFTER filtering
    data.data[data.data < 1e-15] = 0.0
    data.eliminate_zeros()
    
    # 3. renormalise the rows
    row_sums = np.array(data.sum(axis=1)).ravel()
    inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
    data = data.multiply(inv[:, np.newaxis]).tocsr()

    # 3b. Orientation probe (always) and column-mass cap (optional).
    #     colmass(s) = sum_r lambda[r, s] is the market share of seller s, that
    #     is, the weight of s in the propagation operator lambda.T used by
    #     equilibrium(). After the transpose, axis=0 indexes sellers; if the axis
    #     were the wrong way round the maximum would sit near 1, which is what
    #     the probe detects. The cap bounds that mass at C, renormalises the
    #     rows, and iterates. It touches no population, wage or productivity.
    colmass0 = np.array(data.sum(axis=0)).ravel()   # axis=0 = columns = sellers
    print(f"{name} colmass AVANT cap: max={colmass0.max():.3f}, "
          f"n_cols>1.5={(colmass0 > 1.5).sum()}  (croiser diagnostic: ~14 / ~511)")
    if cap_colmass is not None:
        C = float(cap_colmass)
        for _ in range(100):
            colmass = np.array(data.sum(axis=0)).ravel()
            if colmass.max() <= C * (1.0 + 1e-6):
                break
            scale = np.minimum(1.0, C / np.maximum(colmass, 1e-300))
            data = data.multiply(scale[np.newaxis, :]).tocsr()   # scale the columns
            row_sums = np.array(data.sum(axis=1)).ravel()
            inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
            data = data.multiply(inv[:, np.newaxis]).tocsr()      # renormalise the rows
        colmass = np.array(data.sum(axis=0)).ravel()
        print(f"{name} colmass APRES cap C={C}: max={colmass.max():.3f}, "
              f"n_cols>C={(colmass > C*(1+1e-6)).sum()}")

    # 4. diagonal check
    diag = np.array(data.diagonal())
    print(f"diagonale avant fix: min={diag.min():.4f}, median={np.median(diag):.4f}")
    
    # 5. recalibrate the diagonal for market clearing
    income  = None  # passed from main when needed
    rs = np.array(data.sum(axis=1)).ravel()
    print(f"{name} row sums: min={rs.min():.6f}, max={rs.max():.6f}")
    
    return data

def load_migration(file_path, indices_drop):
    with h5py.File(file_path, 'r') as f:
        data = f["mig_flows"][:].T  # transposed: R saves the matrix row by row
    
    data = np.nan_to_num(data, nan=0.0)
    data = csr_matrix(data)
    data.eliminate_zeros()
    
    print(f"migration raw: row sums min={data.sum(axis=1).min():.4f}, max={data.sum(axis=1).max():.4f}")
    
    # Filter to the active locations
    data = data[indices_drop, :][:, indices_drop]
    
    # Renormalise the rows after filtering
    row_sums = np.array(data.sum(axis=1)).ravel()
    inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
    data = data.multiply(inv[:, np.newaxis])
    
    # Repair the empty rows: force the diagonal to 1, meaning stay in place
    row_sums_after = np.array(data.sum(axis=1)).ravel()
    zero_rows = row_sums_after == 0
    if zero_rows.sum() > 0:
        print(f"migration: {zero_rows.sum()} lignes nulles ? diagonale 1")
        data = data + diags(zero_rows.astype(float))
    
    # Final check
    row_sums_final = np.array(data.sum(axis=1)).ravel()
    print(f"migration row sums after renorm: min={row_sums_final.min():.6f}, max={row_sums_final.max():.6f}")
    print(f"migration diag: min={np.array(data.diagonal()).min():.4f}, median={np.median(np.array(data.diagonal())):.4f}")
    

    # Diagnostic on t = 0 only



    return data

