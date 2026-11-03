# <Manuscript title>

**Andrew J. Temple, Kaitlyn A. O'Toole, Michael L. Berumen**

Code and data accompanying the manuscript submitted to *<Journal>*.
<!-- On publication, replace with: Code and data for <full citation>. https://doi.org/<DOI> -->

## Overview

How is research effort distributed across the vertebrate tree of life, and what causes some species to be studied far more than others? We assembled per-species research effort (Scopus publication counts) for 64,707 vertebrate species across nine classes with a global IUCN Red List assessment, together with species traits, geography, national wealth and conservation status. We mapped where research effort falls, resolved its taxonomic structure, and estimated the causal effect of each driver identified in a directed acyclic graph (DAG) using Bayesian negative binomial models with taxonomic random effects, both across all vertebrates and within each of the six major classes. We then used these causal estimates to quantify, by counterfactual prediction, the research each species forgoes because of each bias.

The analyses:

1. **DAG and conditional independence tests:** the hypothesised causal structure and its consistency with the data (`01_dag_development.R`)
2. **Taxonomic structure:** expected publications per species by class, order, family and genus (`02_taxonomic_analysis.R`)
3. **Global causal effects:** one model per driver, adjusted for its backdoor set, across all vertebrates, with E-values (`03_causal_models_global.R`)
4. **Class-level causal effects:** the same models fitted within each of the six major classes (`04_causal_models_by_class.R`)
5. **Research deficits:** counterfactual per-species and total deficits for each driver and for class membership (`05_counterfactual_deficits.R`)
6. **Geography and concentration of research:** global maps of richness and research effort, and the ten most-studied species per class (`06_maps_top10.R`)

## Repository structure

```
vertebrate-research-bias/
├── README.md
├── LICENSE
├── .gitignore
├── session_info.txt                 # Computing environment (from 99_session_info.R)
├── data/                            # Processed, analysis-ready data
│   ├── species_data_final.csv       # Species-level research effort, traits, geography, IUCN status
│   ├── species_top10.csv            # Ten most-studied species per major class
│   ├── global_richness.tif          # Species richness grid
│   ├── global_papersum.tif          # Summed publications grid
│   └── global_mean_papers.tif       # Mean publications per species grid
├── analysis/                        # Run in numbered order from the repository root
│   ├── 01_dag_development.R
│   ├── 02_taxonomic_analysis.R
│   ├── 03_causal_models_global.R
│   ├── 04_causal_models_by_class.R
│   ├── 05_counterfactual_deficits.R
│   ├── 06_maps_top10.R
│   └── 99_session_info.R
└── output/                          # Generated figures, tables and objects (not tracked)
    └── .gitkeep
```

## 1. System requirements

### Software

- R (≥ 4.3.0); RStudio optional
- CmdStan (installed through `cmdstanr`; used as the `brms` backend)
- R packages (versions used to produce the reported results): `ape` (5.8-1), `brms` (2.22.0), `cmdstanr` (0.9.0), `cowplot` (1.1.3), `dagitty` (0.3-4), `data.tree` (1.1.0), `EValue` (4.1.3), `ggforce` (0.5.0), `ggridges` (0.5.6), `ggtree` (3.10.1), `grid` (4.3.0), `lme4` (1.1-37), `patchwork` (1.3.1), `performance` (0.15.0), `phangorn` (2.12.1), `posterior` (1.6.1), `rphylopic` (1.5.0), `rstan` (2.32.7), `scales` (1.4.0), `terra` (1.8-54), `tidybayes` (3.0.7), `tidyterra` (0.7.2), `tidyverse` (2.0.0; including `ggplot2` 3.5.2, `dplyr` 1.1.4, `tidyr` 1.3.1, `readr` 2.1.5); CmdStan 2.39.0
- Full environment: [`session_info.txt`](session_info.txt)
- Internet access when running scripts 02, 04, 05 and 06, which download PhyloPic silhouettes with `rphylopic`

### Tested on

- macOS 15.8.1 (x86_64), R 4.3.0

### Hardware

No non-standard hardware is required, but the Bayesian models are computationally heavy. Each `brm()` fit runs 4 chains with 3 threads per chain (12 cores). `03_causal_models_global.R` `04_causal_models_by_class.R` fit a large number of models and required 2-3 days of run time on a 16-core, 64GB RAM workstation, plan accordingly.

## 2. Installation guide

1. Install [R](https://cran.r-project.org/) and, optionally, [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install the required packages:

```r
install.packages(c("tidyverse", "dagitty", "lme4", "performance",
                   "tidybayes", "brms", "rstan", "posterior", "ggforce",
                   "cowplot", "patchwork", "data.tree", "phangorn", "ape",
                   "rphylopic", "EValue", "ggridges", "terra", "tidyterra",
                   "scales"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree"))

# cmdstanr (Stan R-universe) and CmdStan
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

Typical install time: <N–M> minutes (longer if `rstan`/`brms` compile from source, plus CmdStan).

## 3. Demo and reproduction

All analyses run on the data provided in `data/`; no separate demo dataset is needed.

From the repository root, run the scripts in order, each in a fresh R session:

```bash
for f in analysis/[0-9][0-9]_*.R; do Rscript "$f"; done
```

Alternatively, open the scripts in RStudio with the working directory set to the repository root and source them in order. Scripts 03 and 04 must run before 05 (which loads their saved models), and 03 before 04 (which reads the global effect sizes); 02 must run before 05.

### Expected output

| Script | Output (in `output/`) | Manuscript item | Description |
|---|---|---|---|
| `01_dag_development.R` | — (console output) | Table 3 | Implied conditional independencies, likelihood-ratio tests and partial R² |
| `02_taxonomic_analysis.R` | `Taxa_Bias.pdf` | Fig. 2 | (a) taxonomic tree with research-effort rings, (b) species vs publication share by class, (c) expected publications per species by order |
| `02_taxonomic_analysis.R` | `class_effects.csv`, `order_effects.csv`, `family_effects.csv`, `genus_effects.csv` | tables S1-S4 | Expected publications per species by taxonomic level |
| `02_taxonomic_analysis.R` | `m_taxa.rds` | — | Taxonomic model (used by 05) |
| `03_causal_models_global.R` | `global_model_effectsize.csv`, `pairwise_biome.csv`, `pairwise_status.csv` | Tables S5-S7 | Global causal effects (IRRs) and pairwise contrasts |
| `03_causal_models_global.R` | `global_model_evalues.csv` | Table S11 | E-values, global models |
| `03_causal_models_global.R` | `m_<driver>.rds` (8 files) | — | Global causal models |
| `04_causal_models_by_class.R` | `Causal_Effects.pdf` | Fig. 3 | (a) global IRRs, (b) class-level IRRs, (c) heat map of IRRs |
| `04_causal_models_by_class.R` | `class_model_effectsize.csv`, `class_model_pairwise_biome.csv`, `class_model_pairwise_status.csv` | Table S5-S7 | Class-level causal effects and pairwise contrasts |
| `04_causal_models_by_class.R` | `class_model_evalues.csv` | Table S11 | E-values, class-level models |
| `04_causal_models_by_class.R` | `m_<driver>_<class>.rds` | — | Class-level causal models (used by 05) |
| `05_counterfactual_deficits.R` | `Deficit.pdf` | Fig. 4 | (a) per-species deficits by driver and class, (b) per-species taxonomic deficit relative to Mammalia |
| `05_counterfactual_deficits.R` | `summary_all_deficits.csv`, `summary_class_deficits.csv`, `class_deficit_summary.csv` | Table S8 | Total deficits by driver, by driver and class, and by class (taxonomic) |
| `05_counterfactual_deficits.R` | `species_deficit_summary.csv`, `species_taxa_deficit_summary.csv` | Table S9 | Per-species deficit distributions |
| `05_counterfactual_deficits.R` | `reference_levels_summary.csv` | Counterfactual anchor for each driver and class |
| `06_maps_top10.R` | `Map.pdf` | Fig. 1 | (a) species richness, (b) total publications, (c) mean publications per species, (d) ten most-studied species per class |
| `99_session_info.R` | `session_info.txt` (repository root) | — | Computing environment |

Under `Rscript`, diagnostic plots (histograms, trace plots, posterior predictive checks and unsaved panels) are written to `Rplots.pdf`.

### Expected run times

Measured on <CPU cores, RAM, OS>.

| Script | Approximate time |
|---|---|
| `01_dag_development.R` | Minutes |
| `02_taxonomic_analysis.R` | 1 day |
| `03_causal_models_global.R` | 1+ days |
| `04_causal_models_by_class.R` | 1 day |
| `05_counterfactual_deficits.R` | 1 hour |
| `06_maps_top10.R` | Minutes |

Random seeds are fixed (`seed <- 2009` in each script). Results are reproducible on the same platform. Small Monte Carlo differences can occur across platforms and compilers.

## 4. Instructions for use

To apply the analysis to other taxa or an updated dataset, replace `data/species_data_final.csv` with a file that has the same columns (see Data dictionary). Taxonomic columns must be nested (each genus in one family, each family in one order, each order in one class). The main settings are:

- the DAG (`dagitty()` string), which sets each driver's adjustment set; it is defined in scripts 01–04 and must be changed in all of them;
- `classes` and `min_n` (minimum species per biome or status level within a class), in scripts 04 and 05;
- the priors (`priors_re`, `priors`, `priors_cls`) in scripts 02–04;
- `n_draws_use` (posterior draws used for deficits) and `anchor_percentile` (continuous counterfactual anchor) in script 05.

## Data dictionary

All files are comma-separated, UTF-8, with a single header row. Missing values are empty cells (read as `NA` by `read.csv()`).

### `data/species_data_final.csv`

64,707 rows, one per species. Read by scripts 01–05. No variable has missing values.

| Variable | Description | Type / unit |
|---|---|---|
| `phylum` | Phylum (IUCN taxonomy) | Categorical |
| `class` | Class (IUCN taxonomy) | Categorical |
| `order` | Order (IUCN taxonomy) | Categorical |
| `family` | Family (IUCN taxonomy) | Categorical |
| `genus` | Genus (IUCN taxonomy) | Categorical |
| `sci_name` | Accepted scientific name (IUCN) | Categorical |
| `status` | IUCN Red List category | Categorical |
| `year` | Year of description, from the taxonomic authority | Numeric |
| `papers` | Scopus-indexed publications referring to the species (accepted name and synonyms, title–abstract–keywords; articles, reviews, book chapters, letters, notes) | Numeric |
| `median_gdp_bil` | Summed GDP of the countries whose land or EEZ overlaps the species' range, using each country's median GDP over the last 10 years | Numeric |
| `median_lat` | Median latitude of the species' range (signed) | Numeric |
| `range_km2` | Area of the species' range | Numeric |
| `body_mass_g` | Body mass, from species trait databases | Numeric |
| `biome` | Biome(s) occupied; species in more than one biome form combined levels | Categorical |
| `dist_min` | Minimum elevation (terrestrial) or depth (aquatic) | Numeric |
| `dist_max` | Maximum elevation or depth | Numeric |
| `dist_extent` | Extent between minimum and maximum elevation or depth | Numeric |

**Missing data:** none.

### `data/species_top10.csv`

60 rows, one per species: the ten most-studied species in each of the six major classes. Read by `06_maps_top10.R`.

| Variable | Description | Type / unit |
|---|---|---|
| `class` | Class | Categorical |
| `sci_name` | Scientific name | Categorical |
| `common_name` | Common name used as the plot label. Leading and trailing spaces offset the label from the point. An asterisk marks species with current or historical taxonomic uncertainty whose publication counts are combined under the historical concept (10 species). | Categorical |
| `papers` | Scopus publications (as in `species_data_final.csv`) | Numeric |
| `rank` | Rank within class by `papers` | Numeric |

**Missing data:** none.

### `data/global_*.tif`

Global grids (GeoTIFF) of species ranges and research effort, read by `06_maps_top10.R`. Cells with zero or infinite values are treated as `NA`.

| File | Description | Size |
|---|---|---|
| `global_richness.tif` | Number of species whose range intersects each cell | 31 MB |
| `global_papersum.tif` | Sum of `papers` across the species in each cell | 47 MB |
| `global_mean_papers.tif` | Mean `papers` per species in each cell | 63 MB |

Marine and terrestrial/freshwater species were gridded on separate layers aligned to the coastline and then combined.

## Citation

<At publication: the full citation and DOI of the paper and of the archived code release.>

## License

Code in `analysis/` is licensed under the [MIT License](LICENSE). Data in `data/` are licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). See [LICENSE](LICENSE).
