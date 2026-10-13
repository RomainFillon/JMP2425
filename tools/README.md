# `tools/` — checks, compilation, comparison

Every check returns `rc=0` when clean, `rc=1` when it finds something, and `rc=2` when it
cannot conclude (a missing input is not a pass). Run them from the package root.

## Checks the chain runs

| Tool | What it checks | Called by |
|---|---|---|
| `check_data.py` | presence and md5 of every input of `DATA.md`, against `00_fetch/MD5_*.tsv`; `--presence` for presence only | first thing in `run_all.sh`; stops the run on failure |
| `run_01c_graine2.sh` | step `01c_graine2`: `01c` at a second seed, then a cell-by-cell comparison with `01_estimation/reference_repro_01c_graine2/` at 1e-8 | `run_all.sh` |
| `run_mig_cross.sh`, `mig_internal_cross.R` | step `mig_cross`: cross-border share of internal migration, from the matrix written by `02b` | `run_all.sh` |
| `valide_matrice_population.R`, `temoin_population_pays.R` | the population implied by the weights, against world and World Bank totals | `00d_poids_population.py` |
| `temoin_rattachement.py`, `anomalies_rattachement_connues.csv` | the attachment of populated pixels outside any polygon, against the list of known anomalies | `00c_matrice_population.py` |
| `temoins_poids.py` | each population weight sums to one over its own denominator (cell or country); `02_GE/lib_poids.R` is its R twin, used by `02a` and `02b` | `00d_poids_population.py`, `ge_calibration_diagnostic.py`, the emitters |
| `gardes_lot74.py` … `gardes_lot82.py` | pure predicates: each one checks that a sentence of the paper still holds on the data | the emitters (`02d`) |
| `known_bad_emission.py` | an emission refused half-way leaves the previous artefacts intact (run in a temporary copy of `output/`) | `02d` (blocking) |
| `check_macros.py` | every macro the paper calls is defined, none is `\PROVISOIRE`, and the consistency of macros, tables and producers (specification grid, arms, bound pairs) | `02d` (non-blocking) |
| `reprise.py` | the declared-resumption switch: when `tools/REPRISE_DELIBEREE.tsv` exists, non-regression guards report instead of stopping | the emitters |

## The nine checks of `releve_pa1`

After `02d`, step `releve_pa1` (`etape_releve_pa1.py`, then `releve_pa1.py`) runs the nine
checks listed in `reference_phase.py` (`CONTROLES`) and writes a record under `logs/`. It never
stops the chain.

| Tool | What it checks |
|---|---|
| `check_macros.py` | as above |
| `check_horizon.py` | the horizon labels of `00_parameters.R` and `00_parameters.py` agree |
| `scan_litteraux.py` | lists the numbers written by hand in the body of the paper (information, `rc=0`) |
| `known_bad_emission.py` | as above |
| `check_refs.py` | every `\ref` has a `\label`, every citation a `.bib` entry |
| `scan_octets_controle.py` | no control character in a generated `.tex` |
| `check_fraicheur.py` | no published artefact is older than the data it reads (informative, not pinned in the reference; rc=1 on the author's machine (two fan charts written seconds before their CSV), rc=0 on a fresh clone: see `README.md`, "Known non-determinism") |
| `check_edits_differes.py` | the workshop's register of deferred edits; not shipped, so the check prints "register not shipped; check skipped" and returns `rc=0` |
| `check_labels.py` | no internal scenario number in a generated artefact |

## Reference, compilation, comparison

| Command | What it does |
|---|---|
| `python tools/reference_phase.py` | every published value (macros, artefact md5s, check codes, compilation) against `REFERENCE_PHASE.tsv`; expected: 0 differences, `rc=0` |
| `python tools/compiler.py` | compiles `JMP_V11.tex` (pdflatex, bibtex, pdflatex ×3) and writes `COMPILATION.tsv`, which `reference_phase.py` reads; `rc=0` when clean |
| `python tools/cmp_rerun.py --ref A --new B` | compares two output directories file by file and classes each difference: identical, timestamp (PDF metadata only), rounding, divergent |

`REFERENCE_PHASE.tsv` and `COMPILATION.tsv` are written by these tools; do not edit them.
Their headers and the console output of the tools are in French.
