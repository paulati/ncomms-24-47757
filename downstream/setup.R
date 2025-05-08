

if (!require(simpleCache, quietly = TRUE))
    install.packages("simpleCache")

if (!require(LOLA, quietly = TRUE))
    BiocManager::install("LOLA")

if (!require(ggthemr, quietly = TRUE))
    devtools::install_github(repo = 'https://github.com/Mikata-Project/ggthemr')

