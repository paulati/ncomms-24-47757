clean_env <- function() {

    if(exists("user_base_path")) {
        rm(user_base_path, envir = globalenv())
    }

}

copy_example_alignment <- function(alignment_id, base_path) {

    # copy example data
    data_base_path <- file.path(base_path, 'data')

    if(exists("user_base_path")) {
        out_base_path <- file.path(user_base_path, alignment_id, 'raw')
    } else {
        tmp_base_path <-  cladeAcc:::pkg_data_tmp_base_path()
        out_base_path <- file.path(tmp_base_path, alignment_id, 'raw')
    }

    if(! file.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    if(alignment_id == "100_way") {
        mammals_data_base_path <- file.path(data_base_path, '100way')
        file_name <- 'chr14sample.maf.gz'
        file_path <- file.path(mammals_data_base_path, file_name)
        out_file_name <- 'chr114.maf.gz'
    } else if(alignment_id == "77_way") {
        aves_data_base_path <- file.path(data_base_path, '77way')
        file_name <- 'chr11sample.maf.gz'
        file_path <- file.path(aves_data_base_path, file_name)
        out_file_name <- 'chr111.maf.gz'
    }

    out_file_path <- file.path(out_base_path, out_file_name)
    if(! file.exists(out_file_path)) {
        file.copy(file_path, out_file_path)
    }


}
