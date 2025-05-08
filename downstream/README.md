# Project Overview

This repository contains a collection of scripts, datasets, and analysis outputs focused on evolutionary conservation and functional analysis of genomic elements across mammals and birds. The repository is structured to facilitate exploration of conserved coding and noncoding elements, gene enrichment analysis, and the identification of accelerated regions. All the results in this repository can be reproduced by following the workflow outlined in the files `AnalysisElements_part1.Rmd`, `AnalysisElements_part2.Rmd`, `AnalysisElements_part3.Rmd`.

## Repository Structure

### `code/` Directory

Contains R scripts and auxiliary files used for plotting and data analysis in the `.Rmd` files.

### `data/` Directory

Contains subdirectories with input data used for the analyses. The input data from sections 1-6 was obtained from the generation of accelerated elements, outlined in the repository: <https://github.com/paulati/acc_regions_mammals_aves>.

1.  `1. mammals_phastCons/`:
    -   Contains three subfolders with PhastCons data for mammalian conserved elements and their intersections:
        -   `intersection_modneu_phastCons100way_sarcopterygii_mammals/`
        -   `modneu_phastCons100way_mammals/`
        -   `modneu_phastCons100way_sarcopterygii/`
2.  `2. aves_phastCons/`:
    -   Contains conserved element data for birds, including:
        -   `intersection/`: PhastCons data of aves and their intersection with sarcopterygii.
3.  `5. mammals accelerated candidate functional regions/`:
    -   Contains candidate functional element data for accelerated regions in mammals.
    -   The file `mammals_acc_candidate_func_elements_all.txt` is a concatenate of `mammal_join_filtered_elements_norm_cod_nonCod_oneacczerogap_candidatefunctionalregions.csv` and `mammals_join_filtered_elements_norm_cod_oneacczerogap_candidatefunctionalregions.csv`.
4.  `6. aves accelerated functional regions/`:
    -   Contains functional element data for accelerated regions in aves.
    -   The file `aves_acc_candidate_func_elements_all.txt.txt` is a concatenate of `aves_join_filtered_elements_norm_cod_nonCod_oneacczerogap_candidatefunctionalregions.csv` and `aves_join_filtered_elements_norm_cod_oneacczerogap_candidatefunctionalregions.csv`.
5.  `7. human chicken transcripts/`:
    -   Contains transcript data for human and chicken species.
6.  `8. TADs/`:
    -   Contains Topologically Associated Domain (TAD) data for different species.
7.  `9. liftover/`:
    -   Contains chain files for liftover between species, including *galGal6ToHg38*.
8.  `10. lola/`:
    -   Data used for LOLA enrichment analysis, including regions and databases such as cistrome and epigenome.

### `outputs/` Directory

Contains outputs in the form of plots, tables, and enriched data from the analysis. All output files are generated based on the workflow in `AnalysisElements_part1.Rmd`, `AnalysisElements_part2.Rmd`, `AnalysisElements_part3.Rmd`.

## System Requirements

-   R version: 4.4.1
-   Dependencies are specified at `ncomms-24-47757/setup.R` and `ncomms-24-47757/downstream/setup.R`
