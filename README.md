# Association study of genetic variation and spectroscopic imaging variants

This repository contains the code and data for the publication "Association study of genetic variation and spectroscopic imaging variants". 

## Data

The spectral measurement data underlying the results presented in this paper are available as a published dataset at SPECCHIO http://sc22.geo.uzh.ch:8080/SPECCHIO_Web_Interface/search, with Keyword: UZH_SG_Nicotiana_attenuata_ASD
OR
at Zenodo https://zenodo.org/records/10700859 

The processed data can be found under the folder "Data", which are the input files for the downstream analyses.

## Code

There are six R scripts in this repository:

- `hsc.R`
- `adjust_spec.R`
- `gapit_indices.R`
- `gapit_sw.R`
- `gapit_hsc.R`
- `plots.R`

## Dependencies

The results and figures were generated with R (version 4.3.0 (2023-04-21)) and the following R packages:

- spectrolab (version 0.0.10)
- tidyverse (version 1.3.2)
