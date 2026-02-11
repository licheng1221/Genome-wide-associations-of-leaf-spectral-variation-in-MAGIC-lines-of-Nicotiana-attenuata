# Genome-wide associations of leaf spectral variation in MAGIC lines of Nicotiana attenuata

This repository contains the code and data for the publication "Genome-wide associations of leaf spectral variation in MAGIC lines of Nicotiana attenuata". 

## Data

The spectral measurement data underlying the results presented in this paper are available as a published dataset at SPECCHIO http://sc22.geo.uzh.ch:8080/SPECCHIO_Web_Interface/search, with Keyword: UZH_SG_Nicotiana_attenuata_ASD
OR
at Zenodo https://zenodo.org/records/10700859 

The processed data can be found under the folder "Data", which are the input files for the downstream analyses.

Place all downloaded datasets and GAPIT-generated Manhattan/QQ plots into the data/ folder before running the scripts.

## Code

There are six R scripts in this repository:

- `Spectra_corr.R`
- `hsc.R`
- `gapit_indices.R`
- `gapit_sw.R`
- `gapit_hsc.R`
- `GWAS_plots.R`

## Dependencies

The results and figures were generated with R (version 4.3.0) and the following R packages:

- tidyverse
- spectrolab
- cowplot
- magick
- data.table
- RColorBrewer
