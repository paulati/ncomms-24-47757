
save_plot <- function(p, go_ontology, clade) {

    plot_base_path <- file.path(go_terms_data_base_path, clade, 'dot_plot')
    if(! dir.exists(plot_base_path)) {
        dir.create(plot_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    plot_file_name <- paste0('acc_nc_', clade, '_go', go_ontology, '_top',
                             max_terms_count, 'fe')
    svg_plot_file_name <- paste0(plot_file_name, '.svg')
    pdf_plot_file_name <- paste0(plot_file_name, '.pdf')
    ggplot2::ggsave(file = file.path(plot_base_path, svg_plot_file_name),
                    plot = p, width=29.7, height=21, units = "cm")
    ggplot2::ggsave(file = file.path(plot_base_path, pdf_plot_file_name),
                    plot = p, width=29.7, height=21, units = "cm")

}

dot_plot_clade_ontology <- function(go_terms_data_base_path, go_ontology, clade) {

    data_base_path <- file.path(go_terms_data_base_path, clade, 'rgreat_sim_significant',
                                'union')
    file_name <- paste0('nc_acc_', clade, '_go', go_ontology,
                        '_regions_empirical_p_adjust_annotated_union.csv')
    file_path <- file.path(data_base_path, file_name)
    file.exists(file_path)
    data <- read.delim(file_path, sep = '\t', header = TRUE)

    plot_data <- data.frame('name' = character(nrow(data)),
                            stringsAsFactors = FALSE)

    plot_data$name <- paste0(data$id, ' - ', data$description)
    plot_data$fold_enrichment <- data$fold_enrichment
    plot_data$regions_empirical_p_adjust <- apply(data,
                                                  MARGIN = 1,
                                                  function(row){
        result <-     max(
            as.numeric(row['regions_empirical_p_adjust_biomart_annot']),
            as.numeric(row['regions_empirical_p_adjust_default_annot'])
            , na.rm = TRUE)
        return(result)
    })

    plot_data$genes_prop <- apply(data, MARGIN = 1, function(row){
        result <-     max(
            as.numeric(row['genes_prop_biomart_annot']),
            as.numeric(row['genes_prop_default_annot']), na.rm = TRUE)
        return(result)
    })

    plot_data_sort <- plot_data[order(plot_data$fold_enrichment,
                                      decreasing = TRUE), ]

    #plot_data_sort$order <- seq(1, nrow(plot_data_sort))

    max_terms_count <<- 30

    if(nrow(plot_data_sort) < max_terms_count) {
        plot_data_top <- plot_data_sort

    } else {
        plot_data_top <- plot_data_sort[seq(1, max_terms_count), ]
    }


    font_size <- 9

    plot <- ggplot2::ggplot(data = plot_data_top,
                                    ggplot2::aes(x = fold_enrichment, y = name,
                                            color = regions_empirical_p_adjust,
                                            size = genes_prop)) +
        ggplot2::geom_point() +
        ggplot2::labs(x = "regions fold enrichment",
                      y = "",
                      color = "regions empirical p-adjust",
                      size = "genes in term proportion") +
        ggplot2::scale_color_gradient(low = "blue", high = "red") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = font_size),
            legend.title = ggplot2::element_text(size = font_size),
            axis.text.x = ggplot2::element_text(size = font_size),
            axis.text.y = ggplot2::element_text(size = font_size),
            axis.title.x = ggplot2::element_text(size = font_size)

        )


    return(plot)

}

main_base <- function(clade) {

    for(go_ontology in c('bp', 'cc', 'mf')) {

        plot_to_save <- dot_plot_clade_ontology(go_terms_data_base_path, go_ontology, clade)
        save_plot(plot_to_save, go_ontology, clade)
    }

}


main_mamals <- function() {

    go_terms_data_base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/go_analysis/data/output/go_terms'

    clade <- 'mammals'

    main_base(clade)

}

main_aves <- function() {

    go_terms_data_base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/go_analysis/data/output/go_terms'

    clade <- 'aves'

    main_base(clade)

}




# https://stackoverflow.com/questions/76395098/dotplot-of-enrichgo-results-with-all-of-the-ontology-terms-on-same-plot-for-comp
#        aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
#     geom_point() +
#     facet_grid(Category ~ .) +
#     scale_color_gradient(low = "blue", high = "red") +
#     ylab("")


