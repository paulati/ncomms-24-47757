
main <- function() {

    base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/lola'
    data_base_path <<- file.path(base_path, 'data', 'liftover')

    file_name <- '14_non_coding_acc_regions_aves_ranked.bed'
    file_path <- file.path(data_base_path, file_name)

    data <- read.delim(file_path, sep = ' ', header = TRUE)

    data_gr_galgal6 <- GenomicRanges::makeGRangesFromDataFrame(
                                                data, keep.extra.columns = TRUE)

    chain_file_name <- 'galGal6ToHg38.over.chain'
    chain_file_path <- file.path(data_base_path, chain_file_name)
    chain <- rtracklayer::import.chain(chain_file_path)

    ?rtracklayer::liftOver

    data_hg38 <- rtracklayer::liftOver(data_gr_galgal6, chain)

    data_hg38_df <- data.frame(data_hg38)

    out_file_name <- '14_non_coding_acc_regions_aves_ranked_liftover_hg38.csv'
    out_file_path <- file.path(data_base_path, out_file_name)
    write.table(data_hg38_df, out_file_path,
                sep = '\t', quote = FALSE,
                col.names = TRUE, row.names = FALSE)



    galgal6_elements_count <- length(unique(
        data_gr_galgal6@elementMetadata$element_accelerated_ranking_name))

    missing_elements_count <- galgal6_elements_count -
        length(unique(data_hg38_df$element_accelerated_ranking_name))

    missing_elements_prop <- missing_elements_count / galgal6_elements_count

}
