

main_base <- function(clade) {

    go_data_base_path <- file.path(data_base_path, 'output', 'go_terms',
                                  clade, 'rgreat_sim_significant', 'union')

    for(go_ontology in c('bp', 'mf', 'cc')) {

        file_name <- paste0('nc_acc_', clade, '_go', go_ontology,
                            '_regions_empirical_p_adjust_annotated_union.csv')
        file_path <- file.path(go_data_base_path, file_name)

        data <- read.delim(file_path, sep = '\t', header = TRUE)

        exclude_cols <- c("genes_in_term_default_annot",
                          "genes_in_term_and_region_default_annot",
                          "genes_in_term_biomart_annot",
                          "genes_in_term_and_region_biomart_annot")
        keep_cols <- setdiff(colnames(data), exclude_cols)

        result <- data[order(data$fold_enrichment, decreasing = TRUE),
                       keep_cols]

        out_file_name <- stringr::str_replace(file_name, pattern = "\\.csv",
                                              replacement = "_paper\\.csv")
        out_file_path <- file.path(go_data_base_path, out_file_name)
        write.table(result, file = out_file_path,
                    sep = '\t', quote = FALSE,
                    col.names = TRUE, row.names = FALSE)

    }

}

main <- function() {

    data_base_path <<- "/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/go_analysis/data"

    clade <- 'mammals'

    main_base(clade)

    clade <- 'aves'

    main_base(clade)


}
