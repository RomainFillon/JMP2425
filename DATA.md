# Data

The package does not redistribute any third-party data. Each input below is identified by
its exact product and version, a download location, and the md5 of the file we used.
`tools/check_data.py` verifies every file against `00_fetch/MD5_data_inputs.tsv` (public
inputs, one md5 per file) and `00_fetch/MD5_intrants_derives.tsv` (derived inputs).
Attributions that licences impose (CRH package, BACI, the Copernicus sentence for ERA5-derived
files) are in `DATA_AVAILABILITY.md`.

## Raw inputs read by the chain

For a directory, the md5 is that of the set: the md5 of the sorted lines `md5  path`, one
per file, taken from the manifest. Per-file md5s are in `00_fetch/MD5_data_inputs.tsv`.

| Data | Product and version | Download | Licence | Expected path | Size | md5 |
|---|---|---|---|---|---|---|
| Bilateral trade | CEPII BACI, HS92, release `V202401`, years 1995–2022 (28 annual files), with `country_codes_V202401.csv` | `https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37` (direct download, no registration; HS92 release used: see Version column) | Etalab Open Licence 2.0; cite Gaulier and Zignago (2010), CEPII WP 2010-23 | `data_inputs/BACI_CEPII_trade/` | 11.8 GB | `96be4e440172809f10047a8e9b5e1408` |
| ↳ HS–ISIC concordance | WITS concordance H0 → ISIC Rev.3 (`WITS_Concordance_H0_to_I3.CSV`, same directory) | `https://wits.worldbank.org/product_concordance.html`, file `Concordance_H0_to_I3.zip` (HS 1988/92 → ISIC Rev. 3) | World Bank WITS usage conditions (free for research use) | (in the BACI directory) | 0.5 MB | in the set above |
| ↳ Sectoral trade elasticities | θ^k of Caliendo and Parro (2015), `caliendo_parro_2015_thetak/theta_k.csv` | Caliendo and Parro (2015), Table 1, transcribed into `data_inputs/BACI_CEPII_trade/caliendo_parro_2015_thetak/theta_k.csv` (shipped) | Published estimates; no data licence | (in the BACI directory; shipped) | < 0.1 MB | in the set above |
| ↳ GLORIA–ISIC concordance | `gloria_to_isic3.csv` | Built for this paper: 120 GLORIA 057 sectors mapped to 48 ISIC Rev. 3 divisions (the GLORIA ReadMe provides ISIC Rev. 4 concordances only); shipped as `data_inputs/BACI_CEPII_trade/gloria_to_isic3.csv` | Author's construction, MIT (as the code) | (in the BACI directory; shipped) | < 0.1 MB | in the set above |
| Domestic absorption | GLORIA multi-region input–output, release `057`, build `20230314`, 120-sector mother table, years 1995–2019 (25 files) and `GLORIA_ReadMe_057.xlsx` | IELab GLORIA page, <https://ielab.info/labs/ielab-gloria> (access on request, manual download) | CC BY-NC-SA 4.0; not redistributed | `data_inputs/GLORIA/` | 3.2 GB | `2de417de936b398e3db51af14d763d1f` |
| Gravity controls | CEPII Gravity, release `V202211` (12 files; the chain reads `Gravity_V202211.csv`) | `https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8`, version 202211 (`Gravity_V202211.csv`), direct download | Etalab Open Licence 2.0; cite Conte, Cotterlaz and Mayer (2022), CEPII WP 2022-05 | `data_inputs/gravity_CEPII/` | 1.3 GB | `602acfaa19068e0584d84e92be349a99` |
| Bilateral migration | Abel and Cohen, five-year bilateral flows (`bilat_mig.csv`) | figshare collection `https://doi.org/10.6084/m9.figshare.c.4470464` (Abel and Cohen 2019, *Scientific Data* 6:82), file of origin–destination flow estimates | CC BY 4.0 (figshare) | `data_inputs/abel_cohen_migration/bilat_mig.csv` | 9.7 MB | `e5f4a0079a8102aaacb120eb00b5256b` |
| Disasters (control) | EM-DAT public custom request of 2025-10-28 | `https://public.emdat.be` (free registration for academic users); public custom request, EM-DAT version 2025-10-24, file created 28 October 2025: all countries, start years 1990–2020, technological disasters and the geophysical, biological and extra-terrestrial natural subgroups — no meteorological, hydrological or climatological events; 9,040 records. Climate-related events are excluded by design, so that the disaster control does not itself depend on temperature | EM-DAT terms of use: non-commercial research use, no redistribution of the data; cite per `https://doc.emdat.be/docs/legal/citation-policy/` | `data_inputs/EMDAT/public_emdat_custom_request_2025-10-28_6d0f12d0-20cf-4795-b590-f503f33b2473.xlsx` | 2.0 MB | `80add6f00b26ba1bc0c438aad24601ca` |
| Population grids | JRC GHSL `GHS_POP`, release `R2023A`, epochs 1985–2015 (7 GeoTIFFs), WGS84, 30 arc-seconds | GHS-POP R2023A, `https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE` (JRC); 7 epoch tif files, see Path | European Commission open data: free reuse with attribution; cite Schiavina et al. (2023) | `data_inputs/GHSL_POP/` | 2.6 GB | `6e05f8fe9a8739a623e10ead3d4bde44` |
| Administrative boundaries | GADM 4.1 (`gadm_410.gpkg`, and the archive it was unzipped from) and a simplified GADM 3.6 layer (`gadm_simplified.*`, `gadm_countries.gpkg`) | `https://gadm.org/download_world.html`, GADM 4.1 (`gadm_410.gpkg`). `gadm_simplified.*` is a simplified copy of the GADM 3.6 shapefile (same 339,127 features and 59 fields), made on 17 April 2024; the simplification tool and tolerance were not recorded | GADM licence: free for academic and other non-commercial use; redistribution not allowed (`https://gadm.org/license.html`) | `data_inputs/gadm36_shp/` | 6.2 GB | `4e4bcadfec9d2526aa872db325125ea3` |
| Land area | CIESIN GPWv4 revision 11, land area, 1 degree | `https://doi.org/10.7927/H4Z60M4Z` (GPWv4.11 Land and Water Area, NASA SEDAC) | CC BY 4.0 | `data_inputs/GPW/gpw_v4_land_water_area_rev11_landareakm_1_deg.tif` | 0.1 MB | `b55ec373dba77e31cdac27fe719ddf97` |
| Gridded GDP and net migration | Kummu et al.: total GDP 1990–2022 at 30 arc-minutes; net migration 2010–2014 (5-year sum) | GDP: `https://doi.org/10.5281/zenodo.16741980` (Kummu et al. 2025), file `rast_gdpTot_1990_2022_30arcmin.tif` (33 annual bands, GDP PPP in 2017 international dollars). Net migration: `https://doi.org/10.5281/zenodo.7997134` (Niva et al. 2023), file `raster_netMgr_2000_2019_5yrSum.tif`; the shipped `raster_netMgr_2010_2014_5yrSum.tif` is its third layer (`2010_2014`) at 5 arc-minutes, extracted by `00_fetch/extract_netmgr_2010_2014.R` (identical values on all valid cells; only the nodata code differs) | CC BY 4.0 (both) | `data_inputs/Kummu/` | 20 MB | `1c95c8f8422711ca49d81c378e9fff03` |
| Labour share | ILOSTAT `LAP_2GDP_NOC_RT_A` (labour income share of GDP), extracted 2024-09-22, with ILO country groupings | `https://ilostat.ilo.org/data/`, indicator `LAP_2GDP_NOC_RT` (Labour income share as a percent of GDP, ILO modelled estimates) | CC BY 4.0 (ILO open data policy, datasets released from May 2023; earlier vintages: ILO copyright, research use permitted) | `data_inputs/ILO_labor_share_income/` | 0.9 MB | `2e0fd798feada96f52c28d42ea483a4e` |
| GDP by country | World Bank, GDP per capita, constant 2015 US$ (`NY.GDP.PCAP.KD`), with UN codes | `https://data.worldbank.org/indicator/NY.GDP.PCAP.KD` (series code read in the file's Indicator Code column) | CC BY 4.0 | `data_inputs/WB_gdppercapita_country/` | 0.3 MB | `adcb24807be62531039185233f89469e` |
| Population by country | World Bank population (`population_data.csv`) and a 1990–2020 extract | `https://data.worldbank.org/indicator/SP.POP.TOTL` | CC BY 4.0 | `data_inputs/WB_pop_country/` | 0.4 MB | `71ae1d0abd18576f78ac1274998c91af` |
| Population projections | Wittgenstein Centre, SSP2, `V13` (with and without migration), with UNSD methodology table | `https://dataexplorer.wittgensteincentre.org/wcde-v3/`, SSP2 population by country; cite K.C. et al. (2024), IIASA WP-24-003 | IIASA / Wittgenstein Centre terms of use (free for non-commercial research with attribution) | `data_inputs/Wittgenstein_population_SSP2/` | 193 MB | `cfef30dc8a1ff740aebdddc9cba42bb7` |
| Productivity growth | SSP2 GDP projections, IAMC format (`iamc_db.xlsx`) | IIASA SSP Public Database, Version 2.0 (`https://tntcat.iiasa.ac.at/SspDb`), export generated 12 April 2024, file `iamc_db.xlsx`: SSP2, variable `GDP\|PPP` (billion US$2005/yr), models OECD Env-Growth (Dellink et al. 2017), IIASA GDP, and PIK GDP-32 (Leimbach et al. 2017; United States only) | IIASA SSP Database terms of use (free for non-commercial research with attribution; cite as on `https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=citation`) | `data_inputs/SSP2_productivity_growth/iamc_db.xlsx` | 0.1 MB | `65eb9cb3900d73ec137e15e7aa510b6a` |
| CRH replication package | Cruz and Rossi-Hansberg, *The Economic Geography of Global Warming*, replication package | Zenodo, doi:10.5281/zenodo.7681424 | CC BY 4.0 | `data_inputs/EGGW_Replication_Package.zip` | 2.3 GB | `3d1c6bce20bded3c27944f2c0073a21e` |
| Daily temperature (raw) | ERA5 single levels (not ERA5-Land), 2 m temperature `t2m`, daily mean, 0.25°, 1985–2019, one netCDF per year | Copernicus Climate Data Store; one NetCDF file of daily mean 2 m temperature (`t2m`, K, 0.25°, 1440 × 721, longitudes 0–360) per year, 1985–2019, 35 files; the CDS dataset name and request were not recorded | Copernicus C3S licence (free, attribution: « Contains modified Copernicus Climate Change Service information 2026 ») | processed on a server, not in the package (see derived inputs) | about 97 GB | per-file md5 recorded for 19 of the 35 files (1985–1999, 2005–2008) in `intermediate_outputs/era5_inventory.csv`; 2000–2004 and 2009–2019 not recorded |
| Projected temperature (raw) | ISIMIP3b, bias-adjusted against W5E5, SSP2-4.5, daily `tas`, models GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 (nine decadal files each) | `https://doi.org/10.48364/ISIMIP.842396.1` (Lange and Büchner 2021), daily `tas`, GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, ssp245, 2015–2100; the dataset also carries UKESM1-0-LL, not used | CC0 1.0; cite the DOI | processed on a server, not in the package (see derived inputs) | 69.0 GB (36 files: 4 GCMs × 9 decadal files, daily `tas`, ssp245, 2015–2100) | not recorded for the raw files |
| Human Development Index | `iso3, year, HDI`, 5,746 rows; used by one robustness task | Kummu et al. (2018), Dryad `https://doi.org/10.5061/dryad.dk1j0`, file `HDI_1990_2015_v2.nc` (5 arc-minutes, 1990–2015), averaged to countries by area-weighted zonal means over national polygons (`00_fetch/oh_extract_hdi.R`, one-time extraction; its output `output/oh_hdi_country_year.csv` is shipped); inputs of the extraction: `HDI_1990_2015_v2.nc` from the Dryad deposit and `polyg_adm0_dataNetMgr.gpkg` from the Niva et al. Zenodo deposit (doi:10.5281/zenodo.7997134) | CC0 1.0 (Dryad) | `output/oh_hdi_country_year.csv` | 0.2 MB | `3bd664e3ce03e941ec71e9fa06f5820e` |

## Derived inputs produced on a server

The raw climate files are too large to process locally. Two server steps,
`01_estimation/00_era5_distribution.R` and `02_GE/00_temperature_projections_JMP.R`, turn them
into the files below; the local chain starts from these. They cannot be downloaded from
the original providers in this form. `tas_ssp2_annual.csv` (8.3 MB) is shipped in the
repository; the three others exceed what a repository can hold and are deposited on Zenodo,
<https://doi.org/10.5281/zenodo.23318496> (CC BY 4.0). Download the three files and place
them in `intermediate_outputs/`; their md5s are in the table below and in
`00_fetch/MD5_intrants_derives.tsv`, which `tools/check_data.py` reads.

| File | Produced by | Size | md5 |
|---|---|---|---|
| `intermediate_outputs/ERA_tas_distribution_grid_annual.csv` | `00_era5_distribution.R`, from ERA5 | 622 MB | `33e28df127942c0ea9243e1fb5e4d0ff` |
| `intermediate_outputs/CMIP6_tas.csv` | `00_temperature_projections_JMP.R`, from ISIMIP3b | 462 MB | `b698f887c6eea4257958ab2b7b524b71` |
| `intermediate_outputs/tas_ssp2_daily.csv` | `00_temperature_projections_JMP.R`, from ISIMIP3b | 224 MB | `aaf4697057bb82e284a5104fb65d409d` |
| `intermediate_outputs/tas_ssp2_annual.csv` | `00_temperature_projections_JMP.R`, from ISIMIP3b | 8.3 MB | `6b8d68d55e11d56f5688114da7e882f7` |

## Obtaining the data

1. Download each raw input from the location in the table and place it at its expected path.
   Directory names and file names must match exactly.
2. `data_inputs/gadm36_shp/`: unzip `gadm_410-gpkg.zip` into `gadm_410.gpkg` in the same directory.
3. Download the three large server-derived files from Zenodo
   (<https://doi.org/10.5281/zenodo.23318496>) and place them in `intermediate_outputs/`;
   `tas_ssp2_annual.csv` is already there.
4. `output/oh_hdi_country_year.csv` is shipped with the repository.
5. ERA5 and ISIMIP3b raw files are needed only to re-run the two server steps, not the chain.
6. Check everything with `python tools/check_data.py` (about 30 GB are read; add
   `--presence` for a quick check of presence only).
7. The script exits 1 and lists each missing or different file.
8. `bash run_all.sh` runs this check first and stops before any step if it fails.
9. Totals: about 28 GB of raw inputs and 1.3 GB of derived inputs.
10. Per-file md5s: `00_fetch/MD5_data_inputs.tsv` and `00_fetch/MD5_intrants_derives.tsv`.
