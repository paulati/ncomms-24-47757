great_whole_genome <- function(input_regions, gene_sets,
                               tss_source, biomart_dataset = NA,
                               ncores = 2) {

    if(is.na(biomart_dataset)) {
        results <- rGREAT::great(
            gr = input_regions,
            background = NULL,
            gene_sets = gene_sets,
            tss_source = tss_source,
            cores = ncores)
    } else {
        results <- rGREAT::great(
            gr = input_regions,
            background = NULL,
            gene_sets = gene_sets,
            tss_source = tss_source,
            biomart_dataset = biomart_dataset,
            cores = ncores)

    }


    # Extract significant terms
    results_table <- rGREAT::getEnrichmentTables(results)

    result <- list('obj' = results,
                   'enrichment' = results_table)

    return(result)


}

# go_ontology bp cc mf
great_hg38_whole_genome <- function(input_regions, go_ontology,
                                    use_biomart = FALSE) {

    gene_sets <- switch(go_ontology,
               'bp' = "GO:BP",
               'mf' = "GO:MF",
               'cc' = "GO:CC")

    tss_source <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    biomart_dataset <- 'hsapiens_gene_ensembl'

    if(use_biomart) {
        result <- great_whole_genome(input_regions, gene_sets, tss_source,
                                     biomart_dataset)
    } else {
        result <- great_whole_genome(input_regions, gene_sets, tss_source,
                                     biomart_dataset = NA)
    }

    return(result)
}

great_galgal6_whole_genome <- function(input_regions, go_ontology, use_biomart) {

    gene_sets <- switch(go_ontology,
                        'bp' = "GO:BP",
                        'mf' = "GO:MF",
                        'cc' = "GO:CC")

    tss_source <- "TxDb.Ggallus.UCSC.galGal6.refGene"
    biomart_dataset <- 'ggallus_gene_ensembl'

    if(use_biomart) {
        result <- great_whole_genome(input_regions, gene_sets, tss_source,
                                     biomart_dataset)
    } else {
        result <- great_whole_genome(input_regions, gene_sets, tss_source,
                                     biomart_dataset = NA)
    }

    return(result)
}




# filter can be 'bin005' 'hyp005' 'and'
# significant_terms <- function(go_terms, filter = 'and') {
#
#     mask1 <- go_terms$p_adjust < 0.05
#     mask2 <- go_terms$p_adjust_hyper < 0.05
#
#     if(filter == 'and') {
#         mask3 <- mask1 & mask2
#         go_terms_sign <- go_terms[mask3, ]
#     } else if (filter == "bin005") {
#         go_terms_sign <- go_terms[mask1, ]
#     } else if (filter == "hyp005") {
#         go_terms_sign <- go_terms[mask2, ]
#     }
#
#     return(go_terms_sign)
#
# }
