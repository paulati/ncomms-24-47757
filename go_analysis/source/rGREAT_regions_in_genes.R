#script_base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/go_analysis'
#source(file.path(script_base_path, 'rGreat_base.R'))


regions_in_genes <- function(clade, ontology, go_annotation_type) {

    # https://jokergoo.github.io/rGREAT/reference/getRegionGeneAssociations-GreatObject-method.html


    if(clade == 'mammals') {

        nc_acc_mammals_file_name <- '13_non_coding_acc_regions_mammals_ranked.bed'
        nc_acc_mammals_file_path <- file.path(data_base_path, 'input',
                                              nc_acc_mammals_file_name)
        file.exists(nc_acc_mammals_file_path)
        nc_acc_mammals <- read.delim(nc_acc_mammals_file_path, sep = ' ',
                                     header = TRUE)
        input_regions <- GenomicRanges::makeGRangesFromDataFrame(nc_acc_mammals)

        if(go_annotation_type == 'default') {
            use_biomart <- FALSE
            x <- org.Hs.eg.db::org.Hs.eg.db
        } else if(go_annotation_type == 'biomart') {
            use_biomart <- TRUE
        }
        rgreat_output <- great_hg38_whole_genome(input_regions, ontology,
                                            use_biomart)
    }  else if(clade == 'aves') {

        nc_acc_aves_file_name <- '14_non_coding_acc_regions_aves_ranked.bed'
        nc_acc_aves_file_path <- file.path(data_base_path, 'input',
                                           nc_acc_aves_file_name)
        file.exists(nc_acc_aves_file_path)
        nc_acc_aves <- read.delim(nc_acc_aves_file_path, sep = ' ',
                                  header = TRUE)

        input_regions <- GenomicRanges::makeGRangesFromDataFrame(nc_acc_aves)

        if(go_annotation_type == 'default') {

            use_biomart <- FALSE
            x <- org.Gg.eg.db::org.Gg.eg.db

        } else if(go_annotation_type == 'biomart') {
            use_biomart <- TRUE
        }

        rgreat_output <- great_galgal6_whole_genome(input_regions, ontology,
                                                 use_biomart)

    } else if(clade == 'cars') {


        # TODO
    }

    rgreat_obj <- rgreat_output$obj

    region_gene_association <- rGREAT:::getRegionGeneAssociations(rgreat_obj)

    result <- data.frame(region_gene_association)

    if(go_annotation_type == 'default') {
        # convert gene symbol to ensembl id

        gene_symbol <- unique(do.call(c, result$annotated_genes))

        ensembl_ids <-  AnnotationDbi::mapIds(x, keys = gene_symbol,
                                              keytype="SYMBOL",
                                              column = "ENSEMBL")

        symbol_ensembl_map <- data.frame('gene_symbol' = gene_symbol,
                                         'ensembl_ids' = ensembl_ids,
                                         stringsAsFactors = FALSE)


        #result <- data.frame('id' = names(xx), stringsAsFactors = FALSE)

        genes_in_region <- lapply(seq(1, nrow(result)), function(i) {

            print(i)

            gene_symbol <- result$annotated_genes[[i]]
            #entrez_ids <- as.character(xx[[i]])
            mask <- symbol_ensembl_map$gene_symbol %in% gene_symbol

            result <- symbol_ensembl_map$ensembl_ids[mask]
            return(result)

        })

        result$annotated_genes <- genes_in_region

    } else {

        # do nothing, genes id in ENSEMBL
        #print("completar!!!")
    }

    return(result)


}

go_terms_gene_sets <- function(clade, ontology, go_annotation_type) {

    if(go_annotation_type == 'default') {

        if(clade == 'aves') {
            # https://support.bioconductor.org/p/9145022/

            xx <- as.list(org.Gg.eg.db::org.Gg.egGO2ALLEGS)
            x <- org.Gg.eg.db::org.Gg.eg.db

        } else if(clade == 'mammals') {

            xx <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
            x <- org.Hs.eg.db::org.Hs.eg.db

        } else if (clade == 'cars') {

                print("TODO")
        }

        entrez_ids_lst <- lapply(seq(1, length(xx)), function(i) {
            as.character(xx[[i]])
        })
        entrez_ids <- unique(do.call(c, entrez_ids_lst))
        ensembl_ids <-  AnnotationDbi::mapIds(x, keys = entrez_ids,
                                              keytype="ENTREZID",
                                              column = "ENSEMBL")

        entrez_ensembl_map <- data.frame('entrez_ids' = entrez_ids,
                                         'ensembl_ids' = ensembl_ids,
                                         stringsAsFactors = FALSE)

        result <- data.frame('id' = names(xx), stringsAsFactors = FALSE)

        genes_in_term <- lapply(seq(1, nrow(result)), function(i) {

            print(i)

            entrez_ids <- as.character(xx[[i]])
            mask <- entrez_ensembl_map$entrez_ids %in% entrez_ids

            result <- entrez_ensembl_map$ensembl_ids[mask]
            return(result)
        })

        result$genes_in_term <- genes_in_term


    } else if(go_annotation_type == 'biomart') {

        if(clade == 'aves') {

            biomart_dataset <- 'ggallus_gene_ensembl'

        } else if(clade == 'mammals') {

            biomart_dataset <- 'hsapiens_gene_ensembl'

        } else if (clade == 'cars') {
            print("TODO")
        }

        gene_sets <- rGREAT::getGeneSetsFromBioMart(biomart_dataset,
                                                      ontology)


        result <- data.frame('id' = character(length(gene_sets)),
                             stringsAsFactors = FALSE)

        result$id <- unlist(lapply(seq(1, length(gene_sets)), function(i) {
            names(gene_sets[i])
        }))

        result$genes_in_term <- lapply(seq(1, length(gene_sets)), function(i) {
           gene_sets[[i]]
        })

    }

    return(result)


}

format_data_result <- function(genes_in_term_sign) {

    data_types <- sapply(genes_in_term_sign, class)
    list_mask <- data_types =='list'
    to_flat <- data_types[list_mask]

    for(col_name in names(to_flat)) {
        print(col_name)
        flat_value <- unlist(lapply(genes_in_term_sign[, col_name],
                             function(x) {
                                paste(x, collapse = " ")}))
        # length(flat_value)
        # genes_in_term_sign[ , paste0(col_name, '_tmp')] <- flat_value
        genes_in_term_sign[ , col_name] <- flat_value
    }

    return(genes_in_term_sign)

}

bind_genes_data <- function(clade, ontology, go_annotation_type) {

    file_name_base <- paste0('nc_acc_', clade, '_go', ontology,
                             '_regions_empirical_p_adjust')
    # if(mode == 'union') {
    #     file_name_base <- paste0(file_name_base, '_union')
    # }
    file_name <- paste0(file_name_base, '.csv')

    file_path <- file.path(go_terms_data_base_path, clade,
                           'rgreat_sim_significant', go_annotation_type,
                           file_name)
    data <- read.delim(file_path, sep = '\t', header = TRUE)

    genes_in_term <- go_terms_gene_sets(clade, ontology, go_annotation_type)

    genes_in_term_sign <- merge(data, genes_in_term, by = 'id',
                                all.x = FALSE, all.y = FALSE)

    genes_and_regions <- regions_in_genes(clade, ontology, go_annotation_type)


    # keep only genes from genes_and_regions associated to significant terms
    genes_in_regions_unique <- unique(unlist(genes_and_regions$annotated_genes))
    #length(genes_in_regions_unique)

    # exists a region associated to the gene, and this gene is associated with
    # a singnificant term
    genes_in_term_sign$genes_in_term_and_region <- lapply(
        genes_in_term_sign$genes_in_term, function(x) {
            result <- intersect(genes_in_regions_unique, x)
            return(result)

        })


    genes_in_term_len <- unlist(lapply(genes_in_term_sign$genes_in_term,
                                       length))
    genes_in_term_and_region_len <- unlist(lapply(
        genes_in_term_sign$genes_in_term_and_region, length))
    genes_in_term_sign$genes_prop <- genes_in_term_and_region_len /
                                        genes_in_term_len
    #hist(genes_in_term_sign$genes_prop)


    result <- format_data_result(genes_in_term_sign)


    out_file_name <- paste0(file_name_base, '_annotated.csv')
    out_file_path <- file.path(go_terms_data_base_path, clade,
                           'rgreat_sim_significant', go_annotation_type,
                           out_file_name)
    write.table(result, out_file_path, sep = '\t',
                col.names = TRUE, row.names = FALSE)

    file.exists(out_file_path)

    return(out_file_path)
}

main_base <- function(clade) {

    # biomart_dataset: NULL, 'hsapiens_gene_ensembl'.  'ggallus_gene_ensembl'
    # The rGREAT ‘biomart_dataset’ parameter allows us to specify different GO gene set annotations.
    # If this parameter is set to null, the default annotations databases used are org.Gg.egGO2ALLEGS
    # (DOI:10.18129/B9.bioc.org.Gg.eg.db) and org.Hs.egGO2ALLEGS (DOI:10.18129/B9.bioc.org.Hs.eg.db)
    # If ‘biomart_dataset’ is set to 'hsapiens_gene_ensembl' or  'ggallus_gene_ensembl' annotations
    # are retrived from the Bioconductor package BioMartGOGeneSets

    data_base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/go_analysis/data'
    go_terms_data_base_path <<- file.path(data_base_path, 'output/go_terms')

    for(ontology in c('bp', 'mf', 'cc')) {

        go_annotation_type <- 'biomart'
        result_file_path <- bind_genes_data(clade, ontology, go_annotation_type)
        result_data_biomart <- read.delim(result_file_path, sep = '\t', header = TRUE)

        go_annotation_type <- 'default'
        result_file_path <- bind_genes_data(clade, ontology, go_annotation_type)
        result_data_default <- read.delim(result_file_path, sep = '\t', header = TRUE)

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


