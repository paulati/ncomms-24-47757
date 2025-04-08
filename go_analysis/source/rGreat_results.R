# combine results from biomart and default annotations
# add genes info (annotations)
# format column names

combine_results_base <- function(clade, ontology, file_name_biomart,
                                 file_name_default, out_file_name) {

    data_base_path_default <- file.path(go_terms_data_base_path,
                                    clade, 'rgreat_sim_significant', 'default')
    data_base_path_biomart <- file.path(go_terms_data_base_path,
                                    clade, 'rgreat_sim_significant', 'biomart')

    file_path_default <- file.path(
        data_base_path_default, file_name_default)

    file.exists(file_path_default)

    file_path_biomart <- file.path(data_base_path_biomart,
                                               file_name_biomart)

    file.exists(file_path_biomart)

    data_default <- read.delim(file_path_default, sep = '\t', header = TRUE)

    data_biomart <- read.delim(file_path_biomart, sep = '\t', header = TRUE)


    data_all <- merge(data_default, data_biomart,
                      all.x = TRUE, all.y = TRUE,
                      by = 'id',
                      suffixes = c('_default_annot', '_biomart_annot'))

    if(nrow(data_all) > 0) {

        intersection_mask <- ! is.na(data_all$description_default_annot) &
            ! is.na(data_all$description_biomart_annot)
        data_all$source[intersection_mask] <- 'default_and_biomart_annot'

        tss_mask <- is.na(data_all$description_biomart_annot) &
            ! is.na(data_all$description_default_annot)
        data_all$source[tss_mask] <- 'default_annot'

        biomart_mask <- is.na(data_all$description_default_annot) &
            ! is.na(data_all$description_biomart_annot)
        data_all$source[biomart_mask] <- 'biomart_annot'


        data_all$description <- apply(data_all, MARGIN = 1, function(row) {

            if(is.na(row["description_default_annot"]))  {
                result <- row["description_biomart_annot"]
            } else {
                result <- row["description_default_annot"]
            }
            return(result)
        })

        data_all$fold_enrichment <- apply(data_all, MARGIN = 1, function(row) {
            result <- max(c(row["fold_enrichment_default_annot"],
                            row["fold_enrichment_biomart_annot"]), na.rm = TRUE)
            return(result)
        })

        data_all$p_adjust <- apply(data_all, MARGIN = 1, function(row) {
            result <- min(c(row["p_adjust_default_annot"],
                            row["p_adjust_biomart_annot"]), na.rm = TRUE)
            return(result)
        })

        data_all$regions_empirical_p_adjust <- apply(data_all, MARGIN = 1,
                                                     function(row) {
            result <- min(c(row["regions_empirical_p_adjust_default_annot"],
                            row["regions_empirical_p_adjust_biomart_annot"]),
                          na.rm = TRUE)
            return(result)
        })


        result_col_names <- unique(c("id", "description", "fold_enrichment",
                                     "source", "p_adjust",
                                     "regions_empirical_p_adjust",
                                     colnames(data_all)))

        result <- data_all[, result_col_names]
    } else {

        result <- data_all

    }

    out_base_path <- file.path(go_terms_data_base_path, clade, 'rgreat_sim_significant',
                               'union')
    if(!dir.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    out_file_path <- file.path(out_base_path, out_file_name)

    write.table(result, out_file_path, sep = '\t',
                col.names = TRUE, row.names = FALSE,
                quote = FALSE)


    return(out_file_path)
}


combine_results_regions <- function(clade, ontology) {

    file_name_default <- paste0('nc_acc_', clade, '_go',
                                    ontology, '_regions_empirical_p_adjust.csv')
    file_name_biomart <- paste0('nc_acc_', clade, '_go',
                                    ontology, '_regions_empirical_p_adjust.csv')

    out_file_name <- paste0('nc_acc_', clade, '_go', ontology,
                            '_regions_empirical_p_adjust_union.csv')

    out_file_path <- combine_results_base(clade, ontology, file_name_biomart,
                                          file_name_default, out_file_name)

    return(out_file_path)
}

combine_results_annotations <- function(clade, ontology) {

    file_name_default <- paste0('nc_acc_', clade, '_go',
                            ontology, '_regions_empirical_p_adjust_annotated.csv')
    file_name_biomart <- paste0('nc_acc_', clade, '_go',
                            ontology, '_regions_empirical_p_adjust_annotated.csv')

    out_file_name <- paste0('nc_acc_', clade, '_go', ontology,
                            '_regions_empirical_p_adjust_annotated_union.csv')

    out_file_path <- combine_results_base(clade, ontology, file_name_biomart,
                                          file_name_default, out_file_name)

    return(out_file_path)
}


main_base <- function(clade) {

    go_terms_data_base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/go_analysis/data/output/go_terms'

    for (ontology in c('bp', 'mf', 'cc')) {
        print(ontology)
        data_file_path_1 <- combine_results_regions(clade, ontology)
        data_1 <- read.delim(data_file_path_1, sep = '\t', header = TRUE)

        data_file_path_2 <- combine_results_annotations(clade, ontology)
        data_2 <- read.delim(data_file_path_2, sep = '\t', header = TRUE)

        print(nrow(data_1))
        print(nrow(data_2))
        print(nrow(data_1) == nrow(data_2))
    }
}


main_mammals <- function() {

    clade <- 'mammals'
    main_base(clade)

}

main_aves <- function() {


    clade <- 'aves'
    main_base(clade)

}

