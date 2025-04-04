
load_raw_data <- function() {

    mammals_file_name <- 'lola_raw_mammals.csv'
    aves_file_name <- 'lola_raw_aves.csv'
    aves_faang_file_name <- 'lola_raw_aves_FAANG.csv'

    mammals_file_path <- file.path(data_base_path, mammals_file_name)
    aves_file_path <- file.path(data_base_path, aves_file_name)
    aves_faang_file_path <- file.path(data_base_path, aves_faang_file_name)

    file.exists(aves_file_path)
    file.exists(mammals_file_path)
    file.exists(aves_faang_file_path)

    data_mammals <- read.delim(mammals_file_path, sep = '\t', header = TRUE)
    data_aves <- read.delim(aves_file_path, sep = '\t', header = TRUE)
    data_aves_faang <- read.delim(aves_faang_file_path, sep = '\t',
                                  header = TRUE)

    result <- list('mammals' = data_mammals,
                   'aves' = data_aves,
                   'aves_faang' = data_aves_faang)

    return(result)

}

filter_raw_data <- function(data, support_threshold, p_value_threshold,
                            odds_ratio_threshold) {

    # data_support_counts <- table(data$support)
    # print(data_support_counts)
    # since I want a general idea of the most frequent enrichment, based on the
    # support counts I will exclude those records with support < 6

    exclusion_mask1 <- data$support < support_threshold

    log_p_value_threshold <- -log(p_value_threshold, base = 10)
    exclusion_mask2 <- data$pValueLog <= log_p_value_threshold

    exclusion_mask3 <- data$oddsRatio < odds_ratio_threshold

    exclusion_mask <- exclusion_mask1 | exclusion_mask2 | exclusion_mask3

    data_sign <- data[! exclusion_mask, ]

    result <- data_sign[order(data_sign$meanRnk, decreasing = FALSE), ]

    return(result)

}


main <- function() {

    base_path <- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/lola'
    data_base_path <<- file.path(base_path, 'data')

    data <- load_raw_data()

    #----

    support_threshold_mammals <- 6
    p_value_threshold <- 0.05
    odds_ratio_threshold <- 1

    data_mammals_sign <- filter_raw_data(data$mammals, support_threshold_mammals,
                                    p_value_threshold, odds_ratio_threshold)

    out_file_name <- 'lola_filtered_mammals.csv'
    out_file_path <- file.path(data_base_path, out_file_name)
    write.table(data_mammals_sign, out_file_path,
                sep = '\t', quote = FALSE,
                col.names = TRUE, row.names = FALSE)

    #----

    support_threshold_aves <- 4
    data_aves_sign <- filter_raw_data(data$aves, support_threshold_aves,
                                         p_value_threshold, odds_ratio_threshold)

    out_file_name <- 'lola_filtered_aves.csv'
    out_file_path <- file.path(data_base_path, out_file_name)
    write.table(data_aves_sign, out_file_path,
                sep = '\t', quote = FALSE,
                col.names = TRUE, row.names = FALSE)


    #----

    support_threshold_aves_faang <- 4
    data_aves_faang_sign <- filter_raw_data(data$aves_faang, support_threshold_aves_faang,
                                      p_value_threshold, odds_ratio_threshold)

    # remove first column (row index)
    data_aves_faang_sign <- data_aves_faang_sign[, -1]

    out_file_name <- 'lola_filtered_aves_FAANG.csv'
    out_file_path <- file.path(data_base_path, out_file_name)
    write.table(data_aves_faang_sign, out_file_path,
                sep = '\t', quote = FALSE,
                col.names = TRUE, row.names = FALSE)


}

