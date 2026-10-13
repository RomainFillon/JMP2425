# Environment

## Machine that produced the results

| | |
|---|---|
| OS | Windows 11 Pro, build 10.0.26200 |
| CPU | Intel Core i7-1165G7, 4 cores, 8 threads |
| RAM | 32 GB (31.7 GB usable) |
| R | 4.3.1 (2023-06-16, ucrt); packages pinned in `renv.lock` |
| Python | 3.14.3; packages pinned in `requirements.txt` |

Line endings: files written by the Python emitters (`values.tex`, `values_empirics.tex`, `tab_binning.tex`, `tab_ge*.tex`) use LF on every platform; files written by R steps follow the platform convention, and on Windows every shipped text file they write (the other `tab_*.tex`) is CRLF.

## System tools

| Tool | Used for | Version on the machine above |
|---|---|---|
| Git Bash (bash, `md5sum`, `cygpath`) | `run_all.sh`, the checks | Git for Windows |
| A LaTeX distribution with `pdflatex`, `bibtex` | compiling `JMP_V11.tex` (`tools/compiler.py`) | MiKTeX 26.5, pdfTeX 4.27 |
| `texcount` | word count of the manuscript | MiKTeX |
| Rtools 4.3 (Windows only) | building from source the R packages that `renv::restore()` installs | Rtools43, GCC 12.3.0 |

GDAL is not required: `run_all.sh` exports `OGR2OGR`, but no step calls it.

## Restoring the libraries

```
Rscript -e "renv::restore()"        # from the package root
pip install -r requirements.txt
```

On Windows, install Rtools 4.3 first:
<https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html>. Most versions locked in
`renv.lock` are no longer available as CRAN binaries for R 4.3, so `renv::restore()` has to
compile them from source. Without Rtools, 68 of the 122 packages fail to install (21 at
compilation, 47 through a failed dependency) and `renv::restore()` stops. With Rtools 4.3, on
the machine above and starting from an empty renv cache, `renv::restore()` took 29 minutes
(1,727 s): 122 packages, 25 installed as binaries, 97 built from source, none failed. No other
system library is needed: Rtools 4.3 provides the libraries the source builds link against,
GDAL 3.8.2, PROJ 9.3.1 and GEOS 3.11.2 (`sf`, `terra`; GEOS for `exactextractr`), udunits2
(`units`), HDF5 1.12.1 (`hdf5r`) and netCDF 4.9.2 (`ncdf4`); `rhdf5` is a Bioconductor
binary that brings its own HDF5 1.10.7 (`Rhdf5lib`).

`renv.lock` records the R packages the code loads and their dependencies (121 packages, one
from Bioconductor: `rhdf5`), at the versions installed on the machine above. Launch R and
`run_all.sh` from the package root: `.Rprofile` activates renv there. `_dependencies.R`
declares what renv cannot detect: the packages that `02a` loads through a vector, and
`ncdf4`, which `terra` needs to write netCDF but only suggests. `requirements.txt` lists the Python packages the scripts
import, at the installed versions; `geopandas` also installs its reading engine `pyogrio`.

Packages added to `_dependencies.R` because they are loaded at run time:

- `rnaturalearthdata` 1.0.0: loaded at run time by `rnaturalearth::ne_countries(scale = "medium")` (`01c`, `01d`, `01e`, `02a`, `03_inference/step3`, `step5`, `step7`), invisible to `renv::dependencies()`.

## What `renv::restore()` does on a fresh clone

`renv.lock` describes exactly the library that produced the results: the 121 packages at the
versions installed on the machine above (do not call `renv::snapshot()`, which would replace
this record with whatever is installed on your machine). renv 1.2.2 resolves version
constraints from the current CRAN metadata rather than from the lock, and therefore installs
neighbouring versions of some packages. What it installs depends on whether Rtools is present
and on the state of CRAN on the day `renv::restore()` is run.

### With the renv cache already populated

Measured on 03/10/2026 on the machine above, without Rtools, with a renv cache that already
held the packages. `renv::restore()` prints eight constraint warnings that do not apply to the
locked versions:

```
'tidyterra' requires 'dplyr >= 1.2.0'       'units'    requires 'Rcpp >= 1.1.0'
'vctrs'     requires 'rlang >= 1.1.7'       'dplyr'    requires 'rlang >= 1.1.7'
'dreamerr'  requires 'stringmagic >= 1.2.0' 'fixest'   requires 'stringmagic >= 1.2.0'
'lubridate' requires 'timechange >= 0.4.0'  'dplyr'    requires 'vctrs >= 0.7.1'
```

It then reports that "the dependency tree was repaired": it installs a newer version of eight
locked packages, and seven packages that are not in the lock at all.

| Package | `renv.lock` | Installed by `renv::restore()` |
|---|---|---|
| cli | 3.6.1 | 3.6.4 |
| cpp11 | 0.4.7 | 0.5.5 |
| curl | 5.0.2 | 6.2.2 |
| dplyr | 1.1.3 | 1.1.4 |
| lifecycle | 1.0.4 | 1.0.5 |
| Rcpp | 1.0.11 | 1.0.14 |
| rlang | 1.1.1 | 1.1.5 |
| vctrs | 0.6.3 | 0.6.5 |
| clipr | not in the lock | 0.8.1 |
| forcats | not in the lock | 1.0.1 |
| fs | not in the lock | 2.1.0 |
| haven | not in the lock | 2.5.5 |
| readr | not in the lock | 2.2.0 |
| tzdb | not in the lock | 0.5.0 |
| vroom | not in the lock | 1.7.1 |

### From an empty cache with Rtools 4.3

Measured on 04/10/2026 on the machine above, with Rtools 4.3 and an empty renv cache.
`renv::restore()` prints no constraint warning and takes 29 minutes (1,727 s). It reports that
"the dependency tree was repaired": it installs a newer version of ten locked packages, and
seven packages that are not in the lock at all.

| Package | `renv.lock` | Installed by `renv::restore()` |
|---|---|---|
| cli | 3.6.1 | 3.6.6 |
| cpp11 | 0.4.7 | 0.5.5 |
| curl | 5.0.2 | 8.0.0 |
| dplyr | 1.1.3 | 1.2.1 |
| lifecycle | 1.0.4 | 1.0.5 |
| Rcpp | 1.0.11 | 1.1.2 |
| rlang | 1.1.1 | 1.3.0 |
| stringmagic | 1.1.2 | 1.3.0 |
| timechange | 0.3.0 | 0.4.0 |
| vctrs | 0.6.3 | 0.7.3 |
| clipr | not in the lock | 0.8.1 |
| forcats | not in the lock | 1.0.1 |
| fs | not in the lock | 2.1.0 |
| haven | not in the lock | 2.5.5 |
| readr | not in the lock | 2.2.0 |
| tzdb | not in the lock | 0.5.0 |
| vroom | not in the lock | 1.7.1 |

A fresh clone restored with `renv::restore()` (9–12 October 2026) reproduced every macro of
`values.tex` and `values_empirics.tex` byte for byte; the point estimates, AIC weights and
counts are identical under the locked library and under the versions `renv::restore()`
installs.

### The second seed of step `01c_graine2` depends on the library

- Until 5 October 2026, the second-seed draw of `01c_graine2` depended on the R library: for the importer scheme, q95 is 1.2227 under the renv project library and 1.2441 under the user library (04/10/2026); the first seed and the other schemes do not move.
- The mechanism was identified on 5 October 2026. The coefficient covariance `V_hat` that enters the Monte Carlo differs between two runs at the order of 1e-13: the library, and the number of threads `fixest` sums over, change the order of floating-point summation. The draw was taken by `MASS::mvrnorm`, whose root is spectral (`eigen`); the sign and order of the eigenvectors are discontinuous in `V_hat`, so a 1e-13 difference can flip a sign and move a draw by O(0.1). This gave two exact states rather than continuous noise. Observed directly on the exporter scheme's `V_hat`: 200 symmetric perturbations at 1e-13 move one `mvrnorm` draw by 0.111, and a Cholesky draw by at most 1.3e-12.
- Remedy: every Monte Carlo draw now uses a Cholesky root, with the same normals in the same order. These are `mc_run` (`01_estimation/lib_drf.R::mc_draw_beta`), `04_contrast/ct5_rank_mc_draws.R` and `03_inference/step4_mc.R`. Negative eigenvalues are truncated to zero there, `V+ = U diag(max(lambda, 0)) U'`, then a pivoted Cholesky is taken. A `V_hat` that is not positive definite stops `mc_run` and `ct5`, with no jitter. `01c` estimates on a single thread.
- The package's values are those of the first run under the Cholesky draw (5 October 2026): `01_estimation/reference_repro_01c_graine2/`. The macros `\EMPohSeed*` are re-emitted from that state by step `02d`.

## Cost of each step of `run_all.sh`

Wall-clock time of each step on a fresh clone of the repository (9–10 October 2026, machine
above), read in its run logs: `run_all.sh --from 01a` up to `02a`,
then `--from 02b` to the end. `matrice`, `poids` and `climat_local` were not run on that clone;
their times are those of an earlier run. Memory is given where a log measured it; elsewhere it
was not measured.

| Step | Time | Peak memory |
|---|---|---|
| `matrice` | 48 min (earlier run) | not measured |
| `poids` | 3 s (earlier run) | not measured |
| `climat_local` | 22 s (earlier run) | not measured |
| `01a` | 24 min | not measured |
| `01a_reemit` | 17 s | not measured |
| `01f` | 1 min 27 s | not measured |
| `01b` | 3 min | not measured |
| `01c` | 1 h 03 (single thread) | 20.4 GB (memory probe of 17/09/2026) |
| `01c_graine2` | 59 min | not measured (same script as `01c`) |
| `01c_reemit` | 9 s | not measured |
| `01d` | 2 h 59 | not measured |
| `crh1`–`crh7` | 22 min in all, of which `crh5` 14 min | not measured |
| `01d_bis` | 2 h 57 | not measured |
| `crhclamp`, `crhfig` | 22 s, 8 s | not measured |
| `01e` | 34 s | not measured |
| `ct1`, `ct2` | 11 min, 10 min | not measured |
| `ct8`, `ct9` | 8 min, 2 min | not measured |
| `ct3`–`ct7` | under 20 s each | not measured |
| `02a` | 1 h 17 | not measured |
| `02b` | 29 min | above 19 GB (a run was stopped at 19 GB on 30/09/2026) |
| `mig_cross` | 2 min 41 s (added 11/10/2026; timed in the working copy, not on the clone) | not measured |
| `colmass`, `calibdiag`, `miglat` | 1 min 14 s, 57 s, 21 s | not measured |
| `ge` | 12 h 26 for 78 solver runs | 14 GB of free memory required per run |
| `masks` | 1 min 35 s | not measured |
| `maskruns` | 40 min | not measured |
| `maskdec`, `band` | 12 s, 5 s | not measured |
| `02d` | 9 min | not measured |
| `releve_pa1` | 13 min after `02d` (11/10/2026, working copy); 40 s on the clone, where it ran after `ct7` and stopped early on a missing input | not measured |

`crh5` measures the cost of its bootstrap before running it and skips a channel whose
extrapolated cost exceeds 75 min. On the clone the amenity channel was estimated at 83 min
and skipped, so `crh5` took 14 min; when that channel runs, `crh5` takes about 1 h 20. The
skipped channel is a control and enters no published number.

In all, about 13 hours of steps plus 12 h 26 of general equilibrium: about 25 and a half hours
(the sum of the table above).
A machine with 32 GB of RAM is needed; close other programs during `01c`, `02b` and `ge`.

## Disk space

About 40 GB, measured on the machine above: raw inputs in `data_inputs/` about 28 GB,
`intermediate_outputs/` about 9.5 GB after a full run, `output/` about 1 GB, and the R
packages of `renv.lock` about 0.7 GB.
