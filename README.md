# Revisiting the Global Distribution of Climate Damages — replication package

Romain Fillon, job market paper (2026).

## Overview

This package reproduces the paper `JMP_V11.tex`. The chain estimates the damage functions
(trade and migration, then their robustness and the Cruz and Rossi-Hansberg triangulation),
calibrates and runs the general-equilibrium model, and emits every number, table and figure the
paper prints: the macros of `output/values.tex` and `output/values_empirics.tex`, 19 tables and
9 figures. `run_all.sh` runs it end to end, in about 25 and a half hours on a 4-core laptop.

Two steps that process the raw climate files run on a server and are not replayed by the chain;
their three outputs are deposited on Zenodo (<https://doi.org/10.5281/zenodo.23318496>).

A fresh clone of this repository (full run of 9–12 October 2026, Windows 11, R 4.3.1)
reproduced every macro of `values.tex` and `values_empirics.tex` and all tables and figures of
the paper.

Generative AI (Claude, Anthropic) assisted with additional robustness checks, code cleaning, data handling and copy-editing; the design, estimation and conclusions are the author's.

## Data availability and provenance

`DATA.md` lists every input: product, version, download location, licence, expected path,
size and md5. `DATA_AVAILABILITY.md` holds the data availability statement and the author's
statement of rights.

- Public sources: CEPII BACI and Gravity, the WITS concordance, Abel and Cohen bilateral
  migration, EM-DAT, JRC GHSL, GADM, CIESIN GPWv4, Kummu et al., ILOSTAT, World Bank WDI, the
  Wittgenstein Centre, the SSP database, the Cruz and Rossi-Hansberg replication package
  (Zenodo), ERA5 (Copernicus Climate Data Store) and ISIMIP3b.
- GLORIA input–output tables (release 057) are obtained on request from the Global MRIO Lab and
  are not redistributed.
- The three derived files produced on the server are on Zenodo, doi:10.5281/zenodo.23318496,
  under CC BY 4.0; a fourth, `intermediate_outputs/tas_ssp2_annual.csv`, is in this repository.

## Computational requirements

| | |
|---|---|
| Tested on | Windows 11 Pro, Intel Core i7-1165G7 (4 cores), 32 GB of RAM, Git Bash |
| R | 4.3.1, with Rtools 4.3 on Windows; packages pinned in `renv.lock` |
| Python | 3.14.3; packages pinned in `requirements.txt` |
| LaTeX | `pdflatex` and `bibtex` (MiKTeX), to compile the paper |
| GDAL / `ogr2ogr` | not required |
| Memory | `01c` peaks at 20.4 GB; `02b` needs more than 19 GB; each general-equilibrium run needs 14 GB free |
| Disk | about 40 GB: raw inputs about 28 GB, `intermediate_outputs/` about 9.5 GB, `output/` about 1 GB, R packages about 0.7 GB |
| Time | about 25 and a half hours (below) |

Time by block, measured on the machine above (step by step in `ENVIRONMENT.md`):

| Block | Steps | Time |
|---|---|---|
| Population weights and local climate | `matrice`, `poids`, `climat_local` | 48 min |
| Damage functions and robustness | `01a` … `01c_reemit` | 2 h 31 |
| CRH triangulation and inference | `01d`, `crh1`–`crh7`, `01d_bis`, `crhclamp`, `crhfig` | 6 h 19 |
| Diagnostics and annual/daily contrast | `01e`, `ct1`–`ct9` | 34 min |
| General-equilibrium calibration | `02a`, `02b`, `mig_cross`, `colmass`, `calibdiag`, `miglat` | 1 h 51 |
| General-equilibrium runs | `ge` (78 solver runs) | 12 h 26 |
| Exposure masks, band, emission, closing record | `masks` … `releve_pa1` | 1 h 04 |
| **Total** | | **about 25 h 30** |

Random seeds, all fixed in the code:

| File | Seed | Used for |
|---|---|---|
| `01_estimation/00_parameters.R` (`MC_SEED`) | 44, 45, 42, 43 | Monte Carlo draws of the damage functions: productivity daily, productivity annual, amenity daily, amenity annual |
| `run_all.sh` (`SEED_OH_2`, passed as `JMP_MC_SEED_PROD_DAILY`) | 20260930 | second seed of step `01c_graine2` |
| `03_inference/step4_mc.R`, `step5_boot.R` | 20260728 | CRH interval draws and bootstrap |
| `04_contrast/ct5_rank_mc_draws.R` (`SEED_RANK`) | 20260916 (+1 for the annual arm) | ranking of the Monte Carlo draws |
| `02_GE/02a_calibration_at_1deg.r` | 42 | five locations printed as a check; no output depends on it |

The general-equilibrium solver (Python) draws no random numbers.

## Description of programs

```
run_all.sh          the whole chain, in order
01_estimation/      population weights, damage functions (01a trade, 01b migration),
                    robustness (01c), CRH triangulation (01d), diagnostics (01e, 01f)
02_GE/              general-equilibrium calibration (02a, 02b), solver (02c_model/),
                    exposure masks (02g–02i), emitters (02d and emit_*.py)
03_inference/       inference on the Cruz and Rossi-Hansberg triangulation (step1–step7)
04_contrast/        the annual-mean versus daily contrast (ct1–ct9)
00_fetch/           md5 manifests of every input
tools/              the checks the chain runs, compilation, comparison (tools/README.md)
data_inputs/        raw inputs, not shipped (DATA.md), except two author-built tables: gloria_to_isic3.csv and theta_k.csv
intermediate_outputs/  derived files written by the chain (not shipped)
output/             what the paper reads: values*.tex, tab_*.tex, figures
```

`00_fetch/` also holds two one-time extraction scripts, not steps of `run_all.sh`:
`extract_netmgr_2010_2014.R` reads the Niva et al. (2023) net migration raster
(`raster_netMgr_2000_2019_5yrSum.tif`) and writes its 2010–2014 layer to
`data_inputs/Kummu/raster_netMgr_2010_2014_5yrSum.tif`, and `oh_extract_hdi.R` reads the gridded
HDI of Kummu et al. (2018) and writes the shipped `output/oh_hdi_country_year.csv`.

The steps of `run_all.sh`, in order (`bash run_all.sh --list`):

| Step | Script | What it produces |
|---|---|---|
| `matrice` | `01_estimation/00c_matrice_population.py` | cell × country population matrix from GHSL and GADM 4.1 |
| `poids` | `01_estimation/00d_poids_population.py` | population weights per cell and per country |
| `climat_local` | `01_estimation/00e_era5_pays_local.py` | country temperature series from the shipped ERA5 grid |
| `01a` | `01_estimation/01a_productivity_drf.R` | productivity damage functions (trade), their tables and figures |
| `01a_reemit` | `01_estimation/01a_reemit_tab_absorption.R` | Tables 2–3, re-emitted from their CSV files |
| `01f` | `01_estimation/01f_oh_surviving_share.R` | share of the climate variation that survives the fixed effects |
| `01b` | `01_estimation/01b_amenity_drf.R` | amenity damage functions (migration), their tables and figures |
| `01c` | `01_estimation/01c_robustness.R` | robustness: Oster bounds, omitted-heterogeneity schemes |
| `01c_graine2` | `tools/run_01c_graine2.sh` | `01c` at a second seed, compared with frozen witnesses at 1e-8 |
| `01c_reemit` | `01_estimation/01c_reemit_tab_oh_bounds.R` | Table 5, re-emitted from its CSV files |
| `01d` | `01_estimation/01d_crh_triangulation.R` | CRH triangulation, first pass: inputs of the inference chain |
| `crh1`–`crh7` | `03_inference/step1_gate.R` … `step7_period_sensitivity.R` | reproduction checks, maps, covariances, analytical intervals (the published ones), bootstrap witness, deliverable, period sensitivity |
| `01d_bis` | `01_estimation/01d_crh_triangulation.R` | second pass: Table 19, with the intervals of the chain |
| `crhclamp` | `01_estimation/01d_crh_clamp_sensitivity.R` | sensitivity of the CRH results to the clamp |
| `crhfig` | `01_estimation/01d_fig_crh_responses.R` | Figure 9 |
| `01e` | `01_estimation/01e_diagnostics.R` | diagnostics: cross-weighting table, maps |
| `ct1`, `ct2` | `04_contrast/ct1_delta_north_grid.R`, `ct2_between_arm_test.R` | Northern contrast over the specification grid, test between arms |
| `ct8`, `ct9` | `04_contrast/ct8_binning_variants.R`, `ct9_binning_mechanism.R` | binning variants and their mechanism |
| `ct3` | `04_contrast/ct3_weighted_curve_gap.py` | decomposition of the weighted curve gap |
| `ct4` | `04_contrast/ct4_emit_tab_contribA_grid.R` | Tables 7–8, Figure 2 |
| `ct6` | `04_contrast/ct6_fig_drf_paper.R` | Figure 1 |
| `ct5` | `04_contrast/ct5_rank_mc_draws.R` | ranked Monte Carlo draws for the general-equilibrium band |
| `ct7` | `04_contrast/ct7_cold_tail_variants.R` | cold-tail variants |
| `02a` | `02_GE/02a_calibration_at_1deg.r` | one-degree calibration; Figures 6–7 |
| `02b` | `02_GE/02b_mig_flows.R` | bilateral migration matrix (17,048 × 17,048) |
| `mig_cross` | `tools/run_mig_cross.sh` | cross-border share of internal migration |
| `colmass`, `calibdiag` | `02_GE/cap_colmass_diagnostic.py`, `ge_calibration_diagnostic.py` | calibration diagnostics |
| `miglat` | `02_GE/02f_fig_mig_latitude.py` | Figure 5 |
| `ge` | `02_GE/02c_model/main.py` | the 78 general-equilibrium runs |
| `masks`, `maskruns`, `maskdec` | `02_GE/02g_exposure_masks.py`, `02h_exposure_mask_runs.py`, `02i_exposure_mask_decomposition.py` | exposure masks, four masked runs, decomposition |
| `band` | `02_GE/02e_north_contrast_band.py` | Northern contrast at the quantiles of the damage functions |
| `02d` | `02_GE/02d_emit_latex.R` | runs the three emitters (`emit_values_empirics.py`, `emit_values.py`, `emit_figures.py`): `values*.tex`, the general-equilibrium tables, Figures 3–4 |
| `releve_pa1` | `tools/etape_releve_pa1.py` | non-blocking record: re-emission and the nine checks of `tools/README.md` |
| `targets` | (in `run_all.sh`) | closing inventory of what the paper loads and the chain did not produce |

## Instructions to replicators

1. Restore the environments, from this directory (Windows: install Rtools 4.3 first,
   <https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html>):
   ```
   Rscript -e "renv::restore()"
   pip install -r requirements.txt
   ```
   Do not call `renv::snapshot()`. `ENVIRONMENT.md` says what `renv::restore()` installs.
2. Download the three derived files from Zenodo (<https://doi.org/10.5281/zenodo.23318496>) and
   place them in `intermediate_outputs/`.
3. Download the raw inputs listed in `DATA.md` and place each at its expected path.
4. Check the data: `python tools/check_data.py` (add `--presence` for a quick check).
5. Run the chain, from this directory, in Git Bash on Windows or bash elsewhere:
   ```
   bash run_all.sh                  # everything, in order; checks the data first
   bash run_all.sh --list           # step names
   bash run_all.sh --from 02a       # resume at a step, after a stop or a fix
   bash run_all.sh --until 01d      # stop cleanly after a step
   ```
   Launch from the package root: R reads `.Rprofile`, which activates renv and sets
   `JMP_ROOT`. Keep the machine awake for the whole run (Windows:
   `powercfg /change standby-timeout-ac 0`, and restore it afterwards).
6. Compile the paper: `python tools/compiler.py` (pdflatex, bibtex, pdflatex ×3; writes
   `JMP_V11.pdf` and `tools/COMPILATION.tsv`; rc=0 when clean).
7. Check:
   ```
   python tools/reference_phase.py      # every published value against the reference: 0 differences
   md5sum -c MANIFEST_md5.txt           # every shipped file
   ```
   The reference holds byte-level md5 of the shipped outputs as written on Windows. Files
   written by the R steps carry CRLF line endings there and LF on Linux or macOS, so on those
   systems compare them after normalising line endings; files written by the Python emitters
   are LF on every platform.

The run's framing (header, data check, stops, final checks) is printed in English; the logs of the estimation and calibration steps themselves are in French.

## List of tables and figures

Numbers as in the compiled paper. Tables 17 and 18 are typed in `JMP_V11.tex`; their entries
are macros of `values.tex` and `values_empirics.tex`.

| Paper object | Script (step) | Output file |
|---|---|---|
| Table 1 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge.tex` |
| Table 2 | `01_estimation/01a_reemit_tab_absorption.R` (`01a_reemit`) | `output/tab_fe_absorption.tex` |
| Table 3 | `01_estimation/01a_reemit_tab_absorption.R` (`01a_reemit`) | `output/tab_ctrl_absorption.tex` |
| Table 4 | `01_estimation/01c_robustness.R` (`01c`) | `output/tab_oster_delta.tex` |
| Table 5 | `01_estimation/01c_reemit_tab_oh_bounds.R` (`01c_reemit`) | `output/tab_oh_bounds.tex` |
| Table 6 | `02_GE/emit_values_empirics.py` (`02d`) | `output/tab_binning.tex` |
| Table 7 | `04_contrast/ct4_emit_tab_contribA_grid.R` (`ct4`) | `output/tab_robustness_prod_joint.tex` |
| Table 8 | `04_contrast/ct4_emit_tab_contribA_grid.R` (`ct4`) | `output/tab_contrast_grid.tex` |
| Table 9 | `01_estimation/01b_amenity_drf.R` (`01b`) | `output/tab_exogeneity_mrt_migration.tex` |
| Table 10 | `01_estimation/01b_amenity_drf.R` (`01b`) | `output/tab_robustness_amenity_formatted.tex` |
| Table 11 | `01_estimation/01b_amenity_drf.R` (`01b`) | `output/tab_robustness_amenity_annual_formatted.tex` |
| Table 12 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge_aggregate.tex` |
| Table 13 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge_latitude.tex` |
| Table 14 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge_sensitivity.tex` |
| Table 15 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge_north_sens.tex` |
| Table 16 | `02_GE/emit_values.py` (`02d`) | `output/tab_ge_band_ranked.tex` |
| Table 17 | `02_GE/emit_values.py`, `02_GE/emit_values_empirics.py` (`02d`) | macros in `output/values.tex`, `output/values_empirics.tex` |
| Table 18 | `02_GE/emit_values.py` (`02d`) | macros in `output/values.tex` |
| Table 19 | `01_estimation/01d_crh_triangulation.R` (`01d_bis`) | `output/tab_crh_northern.tex` |
| Figure 1 | `04_contrast/ct6_fig_drf_paper.R` (`ct6`) | `output/MC_prod_daily_fan_paper.pdf`, `output/MC_prod_annual_fan_paper.pdf` |
| Figure 2 | `04_contrast/ct4_emit_tab_contribA_grid.R` (`ct4`) | `output/fig_contrast_grid.pdf` |
| Figure 3 | `02_GE/emit_figures.py` (`02d`) | `output/map_welfare_2100_country.pdf` |
| Figure 4 | `02_GE/emit_figures.py` (`02d`) | `output/fig_ge_tiers.pdf` |
| Figure 5 | `02_GE/02f_fig_mig_latitude.py` (`miglat`) | `output/fig_M1_net_by_latitude.png` |
| Figure 6 | `02_GE/02a_calibration_at_1deg.r` (`02a`) | `output/myplot_migration.pdf` |
| Figure 7 | `02_GE/02a_calibration_at_1deg.r` (`02a`) | `output/no_activity_map_smooth_light.pdf` |
| Figure 8 | `01_estimation/01b_amenity_drf.R` (`01b`) | `output/MC_amenity_daily_fan_paper.pdf`, `output/MC_amenity_annual_fan_paper.pdf` |
| Figure 9 | `01_estimation/01d_fig_crh_responses.R` (`crhfig`) | `output/fig_crh_responses.pdf` |
| Numbers in the text | `02_GE/emit_values.py`, `02_GE/emit_values_empirics.py` (`02d`) | `output/values.tex`, `output/values_empirics.tex` |

## Known non-determinism

- Sixteen PDF figures in `output/` carry their creation date in their metadata: fourteen
  written by R (`pdf()`, `ggsave`) and two by cairo. A re-run of `01a`, `01b`, `01c`, `01d`,
  `01e`, `ct4` or `02a` changes their md5 but not their content: compare their content.
  (The figures of `02d` are written without a date and are byte-identical across runs.)
  To compare, copy the shipped `output/` aside before the run, then
  `python tools/cmp_rerun.py --ref output_shipped --new output`: a PDF whose only differences
  are `/CreationDate`, `/ModDate` and `/ID` is classed `horodatage` (timestamp), not divergent.
- `python tools/check_fraicheur.py` compares file timestamps, so its return code depends on the
  order in which a machine wrote the files; it is informative and does not enter the reference
  (`tools/reference_phase.py` runs and prints it, but does not pin it). In the author's run it
  reports two fan charts, `MC_prod_{daily,annual}_fan_paper.pdf`, as older than their inputs:
  `01a` writes them about two seconds before the CSV files they plot. This is the order of
  writing, not a stale figure.

## Citation

Paper: Fillon, Romain (2026). *Revisiting the Global Distribution of Climate Damages*. Job
market paper.

Package: Replication package, GitHub `RomainFillon/JMP2425`; derived inputs: Zenodo,
doi:10.5281/zenodo.23318496.

## Licence

Code: MIT (see `LICENSE`). Derived data on Zenodo: CC BY 4.0. Third-party inputs keep their
own terms (`DATA.md`).
