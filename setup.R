cran_dependencies <- c(
    'BiocManager',
    'devtools',
    'factoextra',
    'ggplot2',
    'aws.s3',
    'R.utils',
    'ape',
    'rstudioapi',
    'rphylopic',
    'ggimage'

    )


bioc_dependencies <- c(
    'GenomicRanges',
    'AnnotationDbi',
    'org.Gg.eg.db',
    'BioMartGOGeneSets',
    'TxDb.Ggallus.UCSC.galGal6.refGene',
    'TxDb.Hsapiens.UCSC.hg38.knownGene',
    'treeio',
    'rtracklayer'
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


if (!require(simpleCache, quietly = TRUE))
    install.packages("simpleCache")

if (!require(LOLA, quietly = TRUE))
    BiocManager::install("LOLA")

if (!require(rphylopic, quietly = TRUE))
    install.packages("rphylopic")


if (!require(ggthemr, quietly = TRUE))
    devtools::install_github(repo = 'https://github.com/Mikata-Project/ggthemr')

if (!require("cladeAcc", quietly = TRUE))
    devtools::install_github(repo = 'https://github.com/paulati/cladeAcc')


# markdown pdf:
# install.packages(c('tinytex', 'rmarkdown'))
# tinytex::install_tinytex()
