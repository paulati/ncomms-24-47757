library(dplyr)

sample_regions <- function(target_regions, universe, nrep, ncores = 4) {

    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl, varlist = c('target_regions', 'universe'),
                            envir = environment())

    result <- parallel::parLapply(cl, seq(1, nrep), function(i) {

        simulated_regions <-
            regioneR::resampleRegions(target_regions,
                                      universe = universe,
                                      per.chromosome = TRUE)
        simulated_regions_sort <- sort(simulated_regions)

        return(simulated_regions_sort)
    })

    parallel::stopCluster(cl)
    return(result)
}


load_data <- function() {

    aves_acc_file_path <- file.path(data_base_path, 'input',
                                    '14_non_coding_acc_regions_aves_ranked.bed')
    aves_acc_data <- read.delim(aves_acc_file_path, sep = ' ', header = TRUE)

    mammals_acc_file_path <- file.path(data_base_path, 'input',
                                '13_non_coding_acc_regions_mammals_ranked.bed')
    mammals_acc_data <- read.delim(mammals_acc_file_path, sep = ' ',
                                   header = TRUE)

    cars_acc_file_path <- file.path(data_base_path, 'input',
                                    'supp14_ncCAR_hg38.bed')
    cars_acc_data <- read.delim(cars_acc_file_path, sep = '\t', header = TRUE)


    aves_cons_file_path <- file.path(data_base_path, 'input',
                                    '12_aves_conserved_noncoding_elements.bed')
    aves_cons_data <- read.delim(aves_cons_file_path, sep = ' ', header = TRUE)


    mammals_cons_file_path <- file.path(data_base_path, 'input',
                                '10_mammals_conserved_noncoding_elements.bed')
    mammals_cons_data <- read.delim(mammals_cons_file_path, sep = ' ',
                                    header = TRUE)


    result <- list(
        'aves_acc' = GenomicRanges::makeGRangesFromDataFrame(aves_acc_data,
                                                keep.extra.columns = TRUE),
        'mammals_acc' = GenomicRanges::makeGRangesFromDataFrame(
            mammals_acc_data, keep.extra.columns = TRUE),
        'cars_acc' = GenomicRanges::makeGRangesFromDataFrame(cars_acc_data,
                                                keep.extra.columns = TRUE),
        'aves_cons' = GenomicRanges::makeGRangesFromDataFrame(aves_cons_data,
                                                keep.extra.columns = TRUE),
        'mammals_cons' = GenomicRanges::makeGRangesFromDataFrame(
            mammals_cons_data, keep.extra.columns = TRUE)
    )

    return(result)
}


liftover_regions <- function(data_gr) {

    path <- file.path(data_base_path, 'input', 'galGal6ToHg38.over.chain')
    chain <-  rtracklayer::import.chain(path)

    result <- rtracklayer::liftOver(data_gr, chain)

    return(result)
}

main <- function() {


    base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/car'
    data_base_path <- file.path(base_path, 'data')

    data <- load_data()

    # length(data$cars_acc)

    # hist(GenomicRanges::width(data$mammals_cons))
    # hist(GenomicRanges::width(data$mammals_acc))

    #------------------------------------------------------
    # long_acc_aves <- data$aves_acc[GenomicRanges::width(data$aves_acc) >= 3000]
    # long_cons_aves <- data$aves_cons[GenomicRanges::width(data$aves_cons) >= 3000]

    # GenomicRanges::seqnames(long_acc_aves)
    # acc_mammals_mask <- GenomicRanges::seqnames(data$mammals_acc) == 'chr7'
    # sum(acc_mammals_mask)
    # acc_mammals_to_intersect <- data$mammals_acc[acc_mammals_mask]
    #
    # long_regions_overlap <- GenomicRanges::findOverlaps(
    #     acc_mammals_to_intersect, long_acc_aves)
    #
    # common_regions_in_mammals_long <- acc_mammals_to_intersect[
    #     S4Vectors::queryHits(long_regions_overlap),] %>%
    #     unique()
    #
    # tmp <- data.frame(acc_mammals_to_intersect)
    #
    # result <- length(common_regions_in_mammals_long)
    #------------------------------------------------------


    # cons_filtered <- data$aves_cons[GenomicRanges::width(data$aves_cons) < 3000]
    # length(cons_filtered)
    # length(data$aves_cons)
    # hist(GenomicRanges::width(cons_filtered))
    #
    # hist(GenomicRanges::width(data$aves_acc))
    # hist(GenomicRanges::width(data$aves_cons))
    # hist(GenomicRanges::width(aves_sample[[1000]]))

    sim_count <- 50000

    sample_data_file_path <- file.path(data_base_path, "samples.RData")

    if(file.exists(sample_data_file_path)) {

        load(sample_data_file_path)

    } else {

        mammals_sample <- sample_regions(target_regions = data$mammals_acc,
                                         universe = data$mammals_cons,
                                         nrep = sim_count, ncores = 10)

        # liftover galgal6 to hg38
        aves_acc_hg38_lst <- liftover_regions(data$aves_acc)

        # length(aves_acc_hg38)
        # length(data$aves_acc)
        #tmp <- lapply(aves_acc_hg38, function(x) length(x))
        #hist(unlist(tmp))
        #table(tmp)

        aves_acc_hg38 <- aves_acc_hg38_lst |>
            unlist() |>
            GenomicRanges::reduce(min.gapwidth = 10L)

        aves_cons_hg38_lst <- liftover_regions(data$aves_cons)
        aves_cons_hg38 <- aves_cons_hg38_lst |>
            unlist() |>
            GenomicRanges::reduce(min.gapwidth = 10L)

        aves_sample <- sample_regions(target_regions = aves_acc_hg38,
                                      universe = aves_cons_hg38,
                                      nrep = sim_count, ncores = 10)

        save(mammals_sample, aves_sample,
             file = sample_data_file_path)

    }

    intersection_count <- lapply(seq(1, sim_count), function(i) {
        gr_set1 <- mammals_sample[[i]]
        gr_set2 <- aves_sample[[i]]

        regions_overlap <- GenomicRanges::findOverlaps(gr_set1, gr_set2)

        common_regions_in_mammals <- gr_set1[
            S4Vectors::queryHits(regions_overlap),] %>%
            unique()

        result <- length(common_regions_in_mammals)

        return(result)

    })


    #hist(unlist(intersection_count))
    plot_data <- data.frame('x' = unlist(intersection_count))
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x=x)) +
        ggplot2::geom_histogram(binwidth=5, color="black", fill="white") +
        ggplot2::labs(x = "Common conserved regions count in random sets") +
        ggplot2::theme_minimal()
    # p

    mask <- unlist(intersection_count) > length(data$cars_acc)

    empirical_p <- (sum(mask) + 1) / (length(mask) + 1)

    plot_file_name <- 'sample_common_conserved_count_distribution'
    svg_plot_file_name <- paste0(plot_file_name, '.svg')
    pdf_plot_file_name <- paste0(plot_file_name, '.pdf')
    ggplot2::ggsave(file = file.path(base_path, svg_plot_file_name),
                    plot = p, width=29.7, height=21, units = "cm")
    ggplot2::ggsave(file = file.path(base_path, pdf_plot_file_name),
                    plot = p, width=29.7, height=21, units = "cm")


}
