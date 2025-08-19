# ncomms-24-47757

This repository contains the analysis code for the study:

"Comparative genomics sheds light on mammalian and avian gene regulation and phenotypic evolution"

(Nature Communications, 2024; ncomms-24-47757)

------------------------------------------------------------------------

## Setup

Dependencies are listed here

<https://github.com/paulati/ncomms-24-47757/blob/main/setup.R>

And can be installed by running this script.

------------------------------------------------------------------------

## Repository Content

-   `accelerated_elements` contains code and example alignments to run the pipeline that calculates accelerated elements for mammals and aves.

-   `downstream` contains code that process the lists of accelerated elements (computed at the previous step) to generate the datasets shared at `downstream/outputs`:

    -   Coding / non-coding analysis

    -   Clustering

    -   LOLA (Genomic Locus Overlap Enrichment Analysis)

    -   `phastbias` analysis on accelerated elements regions

-   `go_analysis` contains code to evaluate GO enrichment on the accelerated elements regions.

-   `lola` contains code to filter `downstream` LOLA results based on statistical significance.

-   `phylogenetic_tree_plot` and `figure5` contain code used to generate parts of the figures presented in this publication.

------------------------------------------------------------------------

## Data Availability

`downstream/outputs`

- `go_analysis`



`LOLA` 
- https://github.com/paulati/ncomms-24-47757/blob/main/lola/data/lola_filtered_mammals.csv
- https://github.com/paulati/ncomms-24-47757/blob/main/lola/data/lola_filtered_aves.csv
- https://github.com/paulati/ncomms-24-47757/blob/main/lola/data/lola_filtered_aves_FAANG.csv




