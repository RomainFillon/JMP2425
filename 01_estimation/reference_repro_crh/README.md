# Frozen witnesses — the CRH inference chain

These files are frozen witnesses: no step of the package writes them. They hold the published
results of the Cruz and Rossi-Hansberg triangulation, against which the inference chain
(`03_inference/step1_gate.R` … `step7_period_sensitivity.R`) checks that it rebuilds the frame
the published curves came from. Pointing the checks at `intermediate_outputs/crh/` instead
would compare `01d` with itself and verify nothing.

| File | Read by | What is compared |
|---|---|---|
| `daily_amen_bins.csv`, `daily_prod_bins.csv` | `step1_gate.R` | the refitted daily-bin coefficients and bin centres, per channel |
| `crh_daily_projection_cells.csv` | `step2_maps.R` | the per-cell daily and annual-mean projections, reproduced by the linear maps |
| `tab_crh_northern_reversal_daily.csv` | `step4_mc.R` | the country levels (daily and mean) and their cell counts |

Each check prints the largest absolute difference (`max|diff|`) and **does not stop the
chain**: there is no tolerance. A non-zero gap is expected, because `01d` regenerates the inputs
of the chain. On the current tree it is about 8.8e-02 on the productivity bins and about 2.2e-04
on the amenity bins. The published intervals come from the chain itself (`step4_mc.R`), not from
these files.
