#library(rGREAT)


enrichment_clade_go_terms_base <- function(regions_clade, go_ontology, clade,
                                      file_name_prefix, use_biomart) {

    regions_clade_gr <- GenomicRanges::makeGRangesFromDataFrame(regions_clade)

    if(clade == 'mammals') {
        go_terms_regions_clade <- great_hg38_whole_genome(regions_clade_gr,
                                                    go_ontology, use_biomart)
    } else if(clade == 'aves') {
        go_terms_regions_clade <- great_galgal6_whole_genome(regions_clade_gr,
                                                    go_ontology, use_biomart)
    } else if(clade == 'cars') {
        go_terms_regions_clade <- great_hg38_whole_genome(regions_clade_gr,
                                                    go_ontology, use_biomart)

    }

    if(use_biomart) {
        folder_name <- 'biomart'
    } else {
        folder_name <- 'default'
    }

    # go_terms_regions_clade_file_name <- paste0('nc_acc_', clade, '_go',
    #                                           go_ontology, '_all.csv')
    go_terms_regions_clade_file_name <- paste0(file_name_prefix, clade, '_go',
                                               go_ontology, '_all.csv')
    go_terms_regions_clade_base_path <- file.path(data_base_path, 'output',
                        'go_terms', clade, 'rgreat_whole_genome', folder_name)
    if(! dir.exists(go_terms_regions_clade_base_path)) {
        dir.create(go_terms_regions_clade_base_path, recursive = TRUE,
                   showWarnings = FALSE)
    }
    go_terms_regions_clade_file_path <- file.path(
        go_terms_regions_clade_base_path, go_terms_regions_clade_file_name)

    write.table(go_terms_regions_clade,
                file = go_terms_regions_clade_file_path,
                sep = ' \t', quote = FALSE,
                col.names = TRUE, row.names = FALSE)

    return(go_terms_regions_clade)

}

enrichment_clade_go_terms <- function(regions_clade, go_ontology, clade,
                                      file_name_prefix) {

    result_default <- enrichment_clade_go_terms_base(regions_clade, go_ontology, clade,
                                  file_name_prefix, use_biomart = FALSE)

    result_biomart <- enrichment_clade_go_terms_base(regions_clade, go_ontology, clade,
                                   file_name_prefix, use_biomart = TRUE)

    result <- list('default' = result_default,
                   'biomart' = result_biomart)

    return(result)

}

# go_terms_regions_clade_ls <- go_terms_nc_acc_cars_all
# hits_stats_sim_ls <- hits_stats_sim
# clade <- 'cars'
# file_name_prefix <- 'nc_acc_'
# use_biomart <- FALSE
significant_clade_go_terms_base <- function(go_terms_regions_clade_ls,
                            hits_stats_sim_ls, go_ontology, clade,
                            file_name_prefix, use_biomart) {


    for(item_name  in names(go_terms_regions_clade_ls)) {

        print(item_name)

        go_terms_regions_clade <- go_terms_regions_clade_ls[[item_name]]
        hits_stats_sim <- hits_stats_sim_ls[[item_name]]


        # regions
        go_terms_stat_values <- go_terms_regions_clade$observed_region_hits
        names(go_terms_stat_values) <- go_terms_regions_clade$id
        go_terms_distibution <- hits_stats_sim$observed_regions_hits_summary
        go_terms_regions_clade$regions_empirical_p_value <-
            rgreat_empirical_pval(go_terms_stat_values,
                                  go_terms_distibution)
        go_terms_regions_clade$regions_empirical_p_adjust <-
            p.adjust(go_terms_regions_clade$regions_empirical_p_value,
                     method="BH")

        ?p.adjust

        # min(go_terms_regions_clade$regions_empirical_p_adjust)

        # genes
        go_terms_stat_values <- go_terms_regions_clade$observed_gene_hits
        names(go_terms_stat_values) <- go_terms_regions_clade$id
        go_terms_distibution <-hits_stats_sim$observed_genes_hits_summary
        go_terms_regions_clade$genes_empirical_p_value <- rgreat_empirical_pval(
            go_terms_stat_values,
            go_terms_distibution)
        go_terms_regions_clade$genes_empirical_p_adjust <-
            p.adjust(go_terms_regions_clade$genes_empirical_p_value,
                     method="BH")

        out_file_name_base <- paste0(file_name_prefix, clade, '_go', go_ontology)
        if(use_biomart) {
            folder_name <- 'biomart'
        } else {
            folder_name <- 'default'
        }
        out_base_path <- file.path(data_base_path, 'output', 'go_terms', clade,
                                   'rgreat_sim_significant', folder_name)
        if(!dir.exists(out_base_path)) {
            dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
        }

        go_terms_regions_clade_file_path <- file.path(out_base_path,
                                                      paste0(out_file_name_base,
                                                             '_all_empiricalp.csv'))
        write.table(go_terms_regions_clade,
                    file = go_terms_regions_clade_file_path,
                    sep = ' \t', quote = FALSE,
                    col.names = TRUE, row.names = FALSE)

        #hist(go_terms_regions_clade$regions_empirical_p_adjust)
        #hist(go_terms_regions_clade$genes_empirical_p_adjust)

        out_file_name <- paste0(out_file_name_base,
                                '_regions_empirical_p_adjust.csv')
        out_file_path <- file.path(out_base_path, out_file_name)
        mask_sign <- go_terms_regions_clade$regions_empirical_p_adjust < 0.05
        # sum(mask_sign)
        result_regions <- go_terms_regions_clade[mask_sign, ]
        write.table(result_regions, file = out_file_path,
                    sep = '\t', quote = FALSE,
                    col.names = TRUE, row.names = FALSE)


        out_file_name <- paste0(out_file_name_base,
                                '_genes_empirical_p_adjust.csv')
        out_file_path <- file.path(out_base_path, out_file_name)
        mask_sign <- go_terms_regions_clade$genes_empirical_p_adjust < 0.05
        # sum(mask_sign, na.rm = TRUE)
        result_regions <- go_terms_regions_clade[mask_sign, ]
        write.table(result_regions, file = out_file_path,
                    sep = '\t', quote = FALSE,
                    col.names = TRUE, row.names = FALSE)

    }




}

significant_clade_go_terms <- function(go_terms_regions_clade_ls,
                                       hits_stats_sim_ls, go_ontology,
                                       clade, file_name_prefix) {

    significant_clade_go_terms_base(go_terms_regions_clade_ls,
                                    hits_stats_sim_ls, go_ontology, clade,
                                    file_name_prefix, use_biomart = FALSE)

    significant_clade_go_terms_base(go_terms_regions_clade_ls,
                                    hits_stats_sim_ls, go_ontology, clade,
                                    file_name_prefix, use_biomart = TRUE)


}



# nc_cons_clade <- nc_cons_mammals
# nc_acc_clade <- nc_acc_mammals
# go_ontology <- 'bp'
# clade <- 'mammals'

simulation_stats_clade_go_terms <- function(nc_cons_clade, nc_acc_clade, go_ontology, clade) {

    universe_gr <- GenomicRanges::makeGRangesFromDataFrame(nc_cons_clade,
                                            keep.extra.columns = TRUE)
    target_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(nc_acc_clade,
                                                  keep.extra.columns = TRUE)

    # tables with observed regions hits / observed gene hits:
    statitics_sim <- background_elements_distribution(target_regions_gr,
                                                universe_gr,
                                                go_ontology, clade)


    return(statitics_sim)

}




main_mammals <- function() {

    #go_ontology <- 'mf'
    # go_ontology <- 'bp'
    #go_ontology <- 'cc'

    #for(go_ontology in c('bp', 'mf', 'cc')) {

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

        # 2. conserved elements go terms
        #go_terms_nc_cons_mammals_all <- conserved_nc_mammals_go_terms(nc_cons_mammals, go_ontology, filter = 'all')
        #go_terms_nc_cons_mammals_sign <- conserved_nc_mammals_go_terms(nc_cons_mammals, go_ontology, filter = 'bin005')


        # 3. calculate rgreat for acc elements
        go_terms_nc_acc_mammals_all <- enrichment_clade_go_terms(nc_acc_mammals,
                                                go_ontology, 'mammals', 'nc_acc_')

        # go_terms_nc_acc_mammals_all <- accelerated_nc_clade_go_terms(nc_acc_mammals,
        #                                                 go_ontology, 'mammals')
        #hist(go_terms_nc_acc_mammals_all$p_adjust)
        #hist(go_terms_nc_acc_mammals_all$p_adjust_hyper)

        # 4. empirical p values for acc elements
        significant_clade_go_terms(go_terms_nc_acc_mammals_all, hits_stats_sim,
                                   go_ontology, 'mammals', 'nc_acc_')
        # sum(mask_sign, na.rm = TRUE)
    # }

}

main_aves <- function() {

    # go_ontology <- 'bp'

    #for(go_ontology in c('bp', 'mf', 'cc')) {

        go_ontology <- 'mf'

        nc_cons_aves_file_name <- '12_aves_conserved_noncoding_elements.bed'
        nc_cons_aves_file_path <- file.path(data_base_path, 'input', nc_cons_aves_file_name)
        file.exists(nc_cons_aves_file_path)
        nc_cons_aves <- read.delim(nc_cons_aves_file_path, sep = ' ', header = TRUE)

        # 2. simulation data:
        nc_acc_aves_file_name <- '14_non_coding_acc_regions_aves_ranked.bed'
        nc_acc_aves_file_path <- file.path(data_base_path, 'input', nc_acc_aves_file_name)
        file.exists(nc_acc_aves_file_path)
        nc_acc_aves <- read.delim(nc_acc_aves_file_path, sep = ' ', header = TRUE)

        hits_stats_sim <- simulation_stats_clade_go_terms(nc_cons_aves, nc_acc_aves,
                                                go_ontology, 'aves')

        # 3. calculate rgreat for acc elements
        go_terms_nc_acc_aves_all <- enrichment_clade_go_terms(nc_acc_aves,
                                                        go_ontology, 'aves',
                                                        'nc_acc_')


        # 4. empirical p values for acc elements
        significant_clade_go_terms(go_terms_nc_acc_aves_all, hits_stats_sim,
                                   go_ontology, 'aves', 'nc_acc_')

#    }


}


main_cars <- function() {

    #go_ontology <- 'mf'
    # go_ontology <- 'bp'
    #go_ontology <- 'cc'

    #for(go_ontology in c('bp', 'mf', 'cc')) {

    go_ontology <- 'cc'

    nc_cons_mammals_file_name <- '10_mammals_conserved_noncoding_elements.bed'
    nc_cons_mammals_file_path <- file.path(data_base_path, 'input', nc_cons_mammals_file_name)
    file.exists(nc_cons_mammals_file_path)
    nc_cons_mammals <- read.delim(nc_cons_mammals_file_path, sep = ' ', header = TRUE)

    # 1. simulation data:
    nc_acc_cars_file_name <- 'supp14_ncCAR_hg38.bed'
    nc_acc_cars_file_path <- file.path(data_base_path, 'input', nc_acc_cars_file_name)
    file.exists(nc_acc_cars_file_path)
    nc_acc_cars <- read.delim(nc_acc_cars_file_path, sep = '\t', header = TRUE)

    hits_stats_sim <- simulation_stats_clade_go_terms(nc_cons_mammals,
                                                      nc_acc_cars, go_ontology, 'cars')

    # 3. calculate rgreat for acc elements
    go_terms_nc_acc_cars_all <- enrichment_clade_go_terms(nc_acc_cars,
                                            go_ontology, 'cars', 'nc_acc_')


    # 4. empirical p values for acc elements
    significant_clade_go_terms(go_terms_nc_acc_cars_all, hits_stats_sim,
                               go_ontology, 'cars', 'nc_acc_')
    # }


}

# global:
# data_base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/go_analysis/data'

# main_mammals()

# main_aves()

