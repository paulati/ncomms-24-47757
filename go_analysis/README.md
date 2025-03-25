## Go term analysis

In this document you will follow the pipeline used to compute go term enrichment for non-coding accelerated regions.

### Software

rGREAT <https://github.com/jokergoo/rGREAT>

### Input data

#### Mammals (ncMARs)

Coordinates: hg38

Conserved elements: `./data/input/10_mammals_conserved_noncoding_elements.bed`

Accelerated elements: `./data/input/13_non_coding_acc_regions_mammals_ranked.bed`

#### Aves (ncAvARs)

Coordinates: galGal6

Conserved elements: `./data/input/12_aves_conserved_noncoding_elements.bed`

Accelerated elements: `./data/input/14_non_coding_acc_regions_aves_ranked.bed`

#### Common accelerated regions (CARs)

hg38 coordinates

Conserved elements: `./data/input/10_mammals_conserved_noncoding_elements.bed`

Accelerated elements: `./data/input/supp14_ncCAR_hg38.bed`

### Pipeline

#### Includes
``` r
source('/u01/home/pbeati/2024/lucia/paper_acelerados/go_analysis/rGreat_base.R')
source('/u01/home/pbeati/2024/lucia/paper_acelerados/go_analysis/rGreat_simulation.R')
```

#### Step 1

``` r
    go_ontology <- 'bp'

    nc_cons_mammals_file_name <- '10_mammals_conserved_noncoding_elements.bed'
    nc_cons_mammals_file_path <- file.path(data_base_path, 'input', nc_cons_mammals_file_name)
    file.exists(nc_cons_mammals_file_path)
    nc_cons_mammals <- read.delim(nc_cons_mammals_file_path, sep = ' ', header = TRUE)

    # 1. simulation data:
    nc_acc_mammals_file_name <- '13_non_coding_acc_regions_mammals_ranked.bed'
    nc_acc_mammals_file_path <- file.path(data_base_path, 'input', nc_acc_mammals_file_name)
    file.exists(nc_acc_mammals_file_path)
    nc_acc_mammals <- read.delim(nc_acc_mammals_file_path, sep = ' ', header = TRUE)

    hits_stats_sim <- simulation_stats_clade_go_terms(nc_cons_mammals,
                                        nc_acc_mammals, go_ontology, 'mammals')
```



#### Step 2
