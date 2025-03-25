cran_dependencies <- c(
    'BiocManager',
    'devtools',
    'factoextra',
    'ggplot2',
    'aws.s3',
    'R.utils',
    'ape',
    'rstudioapi'
    )


bioc_dependencies <- c(
    'GenomicRanges',
    'AnnotationDbi',
    'org.Gg.eg.db',
    'BioMartGOGeneSets',
    'TxDb.Ggallus.UCSC.galGal6.refGene',
    'TxDb.Hsapiens.UCSC.hg38.knownGene'
    )

for(lib_name in cran_dependencies) {

    if (!require(lib_name, quietly = TRUE))
        install.packages(lib_name)
}

for(lib_name in bioc_dependencies) {

    if (!require(lib_name, quietly = TRUE))
        BiocManager::install(lib_name)
}

if (!require(rphast, quietly = TRUE))
    devtools::install_github("CshlSiepelLab/RPHAST")

if (!require(rGREAT, quietly = TRUE))
    devtools::install_github("jokergoo/rGREAT")

