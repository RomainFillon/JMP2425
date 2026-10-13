# Data and Code Availability Statement

## Statement of rights

The author certifies that he has legitimate access to and permission to use all data used in
this paper, and permission to redistribute the derived files deposited on Zenodo.

## The statement

> **Data and code availability.** All code needed to reproduce the results in
> this paper is available in the replication package (GitHub
> `RomainFillon/JMP2425`). The derived climate files that the package starts
> from are deposited on Zenodo, doi:10.5281/zenodo.23318496, under CC BY 4.0.
> The package contains the estimation, general equilibrium and inference code,
> and every input that may be redistributed.
>
> Bilateral trade data are from the CEPII BACI database, version V202401,
> distributed under the Etalab Open Licence 2.0 (Gaulier and Zignago, 2010).
> Climate reanalysis data are from ERA5, produced by the Copernicus Climate
> Change Service (C3S) and ECMWF, and future projections are from ISIMIP3b.
> The country-level temperature series derived from the ERA5 grid are
> aggregated by population weights and therefore constitute modified
> Copernicus information: *Contains modified Copernicus Climate Change Service
> information 2026*. Neither the European Commission nor ECMWF is responsible
> for any use that may be made of the Copernicus information or data it
> contains. The bias-adjusted climate input data are the ISIMIP3b dataset of
> Lange and Büchner (2021), doi:10.48364/ISIMIP.842396.1, released under CC0 1.0
> with the request that data providers be credited by citing the DOI. The four
> GCMs used are a choice of this paper; the dataset also carries UKESM1-0-LL.
> The calibrated geography of the model builds on the replication
> package of Cruz and Rossi-Hansberg (2024), distributed under CC BY 4.0.
>
> The GLORIA multi-regional input-output tables (release 057, covering 1995 to
> 2019) are distributed by the Industrial Ecology Lab at the University of
> Sydney. Access is requested on the IELab GLORIA page
> (https://ielab.info/labs/ielab-gloria) and the tables are downloaded manually.
> The database is licensed CC BY-NC-SA 4.0 and is therefore not redistributed
> here. Release 057 and the years 1995–2019 are read from the file names in
> `data_inputs/GLORIA/`. The package documents the file names, their checksums
> and the release used.
> The raw ERA5 and ISIMIP3b files are likewise not redistributed, for size
> rather than licensing reasons: 28.2 GB of ERA5 daily fields (Climate Data
> Store, free account, registration required) and 69.0 GB of ISIMIP3b daily
> fields, 97.2 GB in all. The package documents the files used and their
> checksums, and the derived files built from them are on Zenodo.

## Why each sentence is there

| Clause | Reason |
|---|---|
| Gaulier and Zignago (2010) | Etalab Open Licence 2.0 requires the citation. The key `gaulier2010baci` is in `References_JMP.bib` but is not cited in `JMP_V11.tex`, so this statement carries it. |
| *"Contains modified Copernicus…"* | The Copernicus licence prescribes this wording for modified products. The derived series are modified twice: spatial mean to one degree (`01_estimation/00_era5_distribution.R`), then population-share aggregation to countries. |
| The liability disclaimer | Required alongside the attribution by the same licence. |
| Cruz and Rossi-Hansberg, CC BY 4.0 | Attribution is the condition of CC BY. |
| GLORIA not redistributed | CC BY-NC-SA 4.0 permits non-commercial redistribution under the same licence; the package does not redistribute the tables, which are obtained on request from IELab. |

The Copernicus sentence is contractual; the scientific credit for ERA5 is Hersbach et al.
(2020), in the bibliography. The package owes both.

## Sources

The licences above were read at their sources on 22 September 2026, except GLORIA (12 October
2026) and ISIMIP3b (13 October 2026). `DATA.md` carries the
licence of each input; the code licence and the licence of the Zenodo deposit are in the
"Licence" section of `README.md`.
