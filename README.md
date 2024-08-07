# Human adolescent brain network development is different for paralimbic versus neocortical zones

[![DOI](https://zenodo.org/badge/800231232.svg)](https://zenodo.org/doi/10.5281/zenodo.12735272)

This repository contains the code for the main analyses of the manuscript "Human adolescent brain network development is different for paralimbic versus neocortical zones" by Lena Dorfschmidt, František Váša, Simon R. White, Rafael Romero-García, Manfred G. Kitzbichler, Aaron Alexander-Bloch, Matthew Cieslak, Kahini Mehta, Theodore D. Satterhwaite , the NSPN consortium, Richard A. Bethlehem, Jakob Seidlitz, Petra E. Vértes, Edward T. Bullmore. 

For details behind these analyses refer to the manuscript: https://doi.org/10.1101/2023.09.17.558126

# Data
All data required to run these analyses can be found at: . Download data from Zenodo and place it into a folder `data/`.

# Requirements
To run all analyses in this publication, you will additionally need to download the code published by [Váša et al. (2018)](https://github.com/frantisekvasa/rotate_parcellation) to estimate the spherical permutation p-values. Download the code here and place them in a folder `scripts/external/`.

# How to Run
To run this code, first download the required data and external scripts (see above). Most of the scripts read in ouputs from other scripts, so the order in which you run them is essential.

1. Generate main results using `scripts/01.morphometric.development.R`. This is the key script, generating all main results.
2. Estimate within-sample replication using `scripts/02.within.sample.replication.R`. This takes *a while*.
3. Estimate structure-function coupling using `scripts/03.structure.function.coupling.R`
4. Estimate functional network metrics using `scripts/04.functional.network.metrics.R`. This takes *a long time*.
5. External replication in HCP-D was performed using, however we cannot include the processed data to run this script. It is included for completeness/so you can run it with your own processing. The outputs of this script are included in the data release, so you will be able to run the replication statistics script in the next step `scripts/05.external.replication.R`
6. Estimate correspondance between results generated in the NSPN and HCP-D sample using `scripts/06.replication.statistics.R`
   






