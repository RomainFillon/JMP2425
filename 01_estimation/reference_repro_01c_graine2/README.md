# Frozen witnesses of step `01c_graine2` — no step of the package writes them.

Origin: step `01c_graine2` of a run of `run_all.sh` on 5 October 2026 (Cholesky Monte Carlo draw, single-threaded estimation), output `output_01c_graine2/`, copied byte for byte on 6 October 2026.
md5: oh_hotslope_schemes.csv b4df649f4e973c61caa585b001606265; oh_scheme_effects.csv fe9ff0050d1a842162f1da8f923f4ae0; oster_delta_full.csv da21b3c220bc9186725fa9ed8825a564.
Library: produced under the renv project library (`ENVIRONMENT.md`, "The second seed of step `01c_graine2` depends on the library").
Seed: `JMP_MC_SEED_PROD_DAILY=20260930` (productivity, daily arm); everything else as in step `01c`.
Check: `tools/run_01c_graine2.sh` compares `output_01c_graine2/` with these three CSV files through `tools/gardes_lot74.py comparer`, cell by cell, at 1e-8 (absolute).
What it compares: same columns, same shape, identical text and booleans, numbers within 1e-8; a missing file is reported ABSENT and stops the step.
