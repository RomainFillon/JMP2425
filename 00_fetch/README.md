# `00_fetch/` — md5 manifests of the inputs

- `MD5_data_inputs.tsv`: the full md5 of each of the 103 raw input files that produced the
  results, with its path relative to `data_inputs/`.
- `MD5_intrants_derives.tsv`: the md5 of the derived inputs (the four server outputs, three of
  them on Zenodo, and `output/oh_hdi_country_year.csv`), with paths relative to the package root.
- `verifier_md5.py`: checks a filled `data_inputs/` against `MD5_data_inputs.tsv`, file by file.
  `python 00_fetch/verifier_md5.py` checks everything; `python 00_fetch/verifier_md5.py GLORIA`
  checks one source directory. Returns 0 if everything matches, 1 if a file is missing or differs.
  It downloads nothing and modifies nothing.

`tools/check_data.py` reads both manifests and is what `run_all.sh` runs first. Where to obtain
each input is in `DATA.md`.
