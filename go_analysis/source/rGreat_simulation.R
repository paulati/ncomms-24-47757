
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

# go_ontology: bp cc mf
compute_great <- function(input_regions_sim, go_ontology, clade,
                          use_biomart = FALSE, ncores = 4) {

    if(use_biomart) {
        sim_base_path <- file.path(data_base_path, 'output',
                                   'simulation', clade, 'whole_genome', 'biomart',
                                   paste0('go', go_ontology))

    } else {
        sim_base_path <- file.path(data_base_path, 'output',
                'simulation', clade, 'whole_genome', 'default',
                paste0('go', go_ontology))
    }


    if(!dir.exists(sim_base_path)) {
        dir.create(sim_base_path, recursive = TRUE, showWarnings = FALSE)

    }

    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl, varlist = c('input_regions_sim',
                                            'go_ontology',
                                            'data_base_path',
                                            'great_hg38_whole_genome',
                                            'great_galgal6_whole_genome',
                                            'great_whole_genome'),
                            envir = environment())

    result <- parallel::parLapply(cl, seq(1, length(input_regions_sim)),
                                  function(i) {


        input_region <- input_regions_sim[[i]]

        go_terms_input_region_file_name <- paste0('sim_input_region_', i,
                                                  '_go', go_ontology, '.csv')
        go_terms_input_region_file_path <- file.path(sim_base_path,
                                            go_terms_input_region_file_name)

        if(! file.exists(go_terms_input_region_file_path)) {
            if(clade == 'mammals') {
                go_terms_input_region <- great_hg38_whole_genome(input_region,
                                                    go_ontology, use_biomart)
            } else if(clade == 'aves') {
                go_terms_input_region <- great_galgal6_whole_genome(
                                        input_region, go_ontology, use_biomart)
            } else if(clade == 'cars') {
                go_terms_input_region <- great_hg38_whole_genome(input_region,
                                                    go_ontology, use_biomart)
            }

            # keep all results do not filter by significance!
            # go_terms_input_region_sign <- significant_terms(
            #     go_terms_input_region, go_terms_input_region_file_path)

            write.table(go_terms_input_region$enrichment,
                        go_terms_input_region_file_path,
                        sep = '\t', quote = FALSE,
                        col.names = TRUE, row.names = FALSE)

            rm()
            gc()

        } else {

            go_terms_input_region <- read.delim(
                go_terms_input_region_file_path, sep = '\t', header = TRUE)
        }

        return(go_terms_input_region)
    })

    parallel::stopCluster(cl)
    return(result)


}


base_hits_summary <- function(files_paths, stat_col_name, out_file_name, clade,
                              use_biomart = FALSE) {

    if(use_biomart) {

        out_base_path <- file.path(data_base_path, 'output', 'simulation',
                                   clade, 'whole_genome_summary', 'biomart')

    } else {

        out_base_path <- file.path(data_base_path, 'output', 'simulation',
                                   clade, 'whole_genome_summary', 'default')

    }

    if(! dir.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    out_file_path <- file.path(out_base_path, out_file_name)
    out_file_path_gz <- file.path(out_base_path, paste0(out_file_name, '.gz'))

    if(! file.exists(out_file_path)) {

        # TODO; paralelizar esto
        rows_from_files <- lapply(seq(1, length(files_paths)), function(i) {

            #print(i)
            file_path <- files_paths[i]
            #data <- read.delim(file_path)
            # read compressed files:
            data <- read.table(file_path, sep = '\t', header = TRUE, quote = "\"")

            sim_id <- as.integer(unlist(stringr::str_split(file_path, pattern = '_'))[7])

            row <- c(sim_id, data[, stat_col_name])
            names(row) <- c('sim_id', data$id)
            rm()
            gc()
            return(row)
        })

        go_terms_all <- unique(unlist(lapply(rows_from_files, function(row) {
            go_terms_row <- names(row)[-1] # exclude sim_id
        })))

        result_mtx <- matrix(0L,
                             nrow = length(rows_from_files),
                             ncol = length(go_terms_all) + 1)
        result <- data.frame(result_mtx)
        colnames(result) <- c('sim_id', go_terms_all)
        for (i in seq(1, length(rows_from_files))) {
            #print(i)
            row_data <- rows_from_files[[i]]
            result[i, names(row_data)] <- row_data
        }

        # length(rows_from_files[[i]])
        # ncol(result) - length(rows_from_files[[i]])
        # sum(is.na(result[i, ]))

        write.table(result, out_file_path, sep = '\t',
                    col.names = TRUE, row.names = FALSE)

        R.utils::gzip(filename = out_file_path, destname = out_file_path_gz,
                      ext="gz", remove = TRUE, overwrite = TRUE)

    } else {

        result <- read.delim(out_file_path_gz, sep = '\t', header = TRUE)
    }

    return(result)

}


observed_genes_hits_summary <- function(files_paths, go_ontology, clade,
                                        use_biomart = FALSE) {

    if(use_biomart) {
        folder_name <- 'biomart'
    } else {
        folder_name <- 'default'
    }

    stat_col_name <- 'observed_gene_hits'
    out_file_name <- paste0('go_terms_observed_genes_hits_', go_ontology, '.csv')
    out_file_path <- file.path(data_base_path, 'output', 'simulation',
                               clade, 'whole_genome_summary', folder_name,
                               paste0(out_file_name, '.gz'))

    print(out_file_path)

    if(! file.exists(out_file_path)) {
        if(anyNA(files_paths)) {
            result <- data.frame()
        } else {
            result <- base_hits_summary(files_paths, stat_col_name, out_file_name,
                                        clade, use_biomart)
        }
    } else {
        result <- read.delim(out_file_path, sep = '\t', header = TRUE)
    }
    return(result)
}

observed_regions_hits_summary <- function(files_paths, go_ontology, clade,
                                          use_biomart = FALSE) {

    if(use_biomart) {
        folder_name <- 'biomart'
    } else {
        folder_name <- 'default'
    }

    stat_col_name <- 'observed_region_hits'
    out_file_name <- paste0('go_terms_observed_region_hits_', go_ontology, '.csv')
    out_file_path <- file.path(data_base_path, 'output', 'simulation',
                               clade, 'whole_genome_summary', folder_name,
                               paste0(out_file_name, '.gz'))

    print(out_file_path)

    if(! file.exists(out_file_path)) {
        if(anyNA(files_paths)) {
            result <- data.frame()
        } else {
            # TODO check if sim files exists, if not run simulation
            result <- base_hits_summary(files_paths, stat_col_name, out_file_name,
                                    clade, use_biomart)
        }
    } else {

        result <- read.delim(out_file_path, sep = '\t', header = TRUE)
    }

    return(result)
}


rgreat_empirical_pval <- function(go_terms_stat_values, go_terms_distibution) {

    #fix:
    colnames(go_terms_distibution) <- stringr::str_replace(
        colnames(go_terms_distibution),
        pattern = "\\.",
        replacement = "\\:")

    empirical_pvalues <- unlist(lapply(seq(1, length(go_terms_stat_values)), function(i) {

        go_id <- names(go_terms_stat_values)[i]
        go_value <- go_terms_stat_values[i]
        if(go_id %in% colnames(go_terms_distibution)) {
            distribution <- go_terms_distibution[, go_id]
            result <- (sum(go_value <= distribution) + 1)/
                (length(distribution) + 1)
        } else {
            print(go_id)
            result <- NA
        }
        return(result)
    }))

    return(empirical_pvalues)

}

background_elements_distribution_base <- function(target_regions_gr, universe_gr,
                                                  go_ontology, clade,
                                                  use_biomart) {

    observed_genes_hits <- observed_genes_hits_summary(files_paths = NA,
                                                       go_ontology, clade,
                                                       use_biomart)
    observed_regions_hits <- observed_regions_hits_summary(files_paths = NA,
                                                           go_ontology, clade,
                                                           use_biomart)

    nrep <- 5000

    run_simulation <-  ! (nrow(observed_genes_hits) == nrep &
                              nrow(observed_regions_hits) == nrep)

    if(use_biomart) {
        folder_name <- 'biomart'
    } else {
        folder_name <- 'default'
    }


    if(run_simulation) {

        go_terms_input_region_base_path <- file.path(data_base_path, 'output',
                                                     'simulation', clade,
                                                     'whole_genome', folder_name,
                                                     paste0('go', go_ontology))
        files <- list.files(go_terms_input_region_base_path)

        if(length(files) == 0) {

            input_regions_sim <- sample_regions(target_regions_gr, universe_gr, nrep)
            go_terms_sim <- compute_great(input_regions_sim, go_ontology, clade,
                                          use_biomart, ncores = 10)

        } else {

            print('sim folder should be empty')
            print(go_terms_input_region_base_path)
        }

        #pongo algo con un wait??
        files <- list.files(go_terms_input_region_base_path, full.names = TRUE)
        if(length(files) < nrep){
            print('wait rgreat evaluation ends')
        } else {

            observed_genes_hits <- observed_genes_hits_summary(files,
                                                go_ontology, clade, use_biomart)
            observed_regions_hits <- observed_regions_hits_summary(files,
                                                go_ontology, clade, use_biomart)
        }

    }

    result <- list(
        'observed_genes_hits_summary' = observed_genes_hits,
        'observed_regions_hits_summary' = observed_regions_hits
    )

    return(result)

}


# build a table where each row is the results of a simulation and
# each column is a go term
# the value of the elements is the calculated statistics ()
#cons_elements_distribution <- function(target_regions_gr, universe_gr,

background_elements_distribution <- function(target_regions_gr, universe_gr,
                                       go_ontology, clade) {


    result_default <- background_elements_distribution_base(target_regions_gr, universe_gr,
                                         go_ontology, clade,
                                         use_biomart = FALSE)

    result_biomart <- background_elements_distribution_base(target_regions_gr, universe_gr,
                                                            go_ontology, clade,
                                                            use_biomart = TRUE)

    result <- list('default' = result_default,
                   'biomart' = result_biomart)

    return(result)

}


