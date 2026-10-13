# Frozen witnesses — the annual/daily contrast chain

These files are frozen witnesses: no step of the package writes them, which is what gives them
their value. The same arrangement serves the CRH chain in `01_estimation/reference_repro_crh/`.

## `weighted_curve_gap_20260806.csv`

The eight published values of `ct3_weighted_curve_gap.py`, measured on 6 August 2026 on the
dose-response curves of that date.

They are deliberately **non-blocking**. The curves they summarise have since been re-estimated,
and a legitimate re-run of the package moves them again; a refusal on these values would be a
check that ends up switched off. The risk a refusal would guard against is a cross-matching of
inputs, and that one is refused: `ct3` recomputes the provenance fingerprint of its six inputs
and stops if it differs from the one `ct1` recorded.

`ct3` therefore **prints the gap** to these eight values and writes it to
`output/pe_weighted_gap_north_scalars.csv`; it does not stop. A non-zero gap is a measurement,
not a failure: it says how far the re-estimation has moved the decomposition since 6 August,
and the ten `\EMPgap*` macros of the paper derive from it.

Do not refresh this file to make a gap disappear. A refresh should say, in this README, what the
new values were checked against.
