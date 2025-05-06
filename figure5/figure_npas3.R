# https://bernatgel.github.io/karyoploter_tutorial/
# https://bernatgel.github.io/karyoploter_tutorial/Tutorial/PlotRegions/PlotRegions.html
# https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotDensity/PlotDensity.html

# plot params
# https://bernatgel.github.io/karyoploter_tutorial/Tutorial/PlotParams/PlotParams.html

# colors
# https://convertingcolors.com/hex-color-013220.html?search=Hex(013220)

# BiocManager::install("Rhtslib")
# BiocManager::install("karyoploteR")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(dplyr)
library(stringr)
library(karyoploteR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)


citation("karyoploteR")

base_directory <- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/figure5'
file.exists(base_directory)

# tyes: M, Av
get_nc_elements <- function(type = "M") {

    ar_base_directory <- file.path(base_directory, 'data')

    header <- TRUE
    chr_column_name <- 'seqnames'
    star_column_name <- 'start'
    end_column_name <- 'end'


    if(type == "M") {
        file_name <- 'ncMAR.csv'
    } else if(type == "Av") {
        file_name <- 'ncAvAR.csv'
    }

    file_path <- file.path(ar_base_directory, file_name)

    data <- read.table(file_path, sep = '\t', header = header)

    return(data)
}

get_c_elements <- function(type = "M") {

    ar_base_directory <- file.path(base_directory, 'data')

    header <- TRUE
    chr_column_name <- 'seqnames'
    star_column_name <- 'start'
    end_column_name <- 'end'


    if(type == "M") {
        file_name <- 'mammals_acc_coding.csv'
    } else if(type == "Av") {
        file_name <- 'aves_acc_coding.csv'
    }

    file_path <- file.path(ar_base_directory, file_name)

    data <- read.table(file_path, sep = '\t', header = header)

    return(data)
}

get_all_genes <- function(chrs = NA) {

    data_base_directory <- file.path(base_directory, 'public_data')
    file_name <- 'hg38.refGene.gtf.gz'

    header <- FALSE
    chr_column_name <- 'V1'
    star_column_name <- 'V4'
    end_column_name <- 'V5'
    strand_column_name <- 'V7'
    gene_id_columns_name <- 'gene_id'

    file_path <- file.path(data_base_directory, file_name)

    raw_data <- read.table(file_path, sep = '\t', header = header)

    data_info_part_lst <- str_split(raw_data$V9, pattern = ";")

    gene_id <- lapply(data_info_part_lst, function(x) {
        str_replace(string = x[[1]], pattern = "gene_id ", replacement = "")
    })

    raw_data$gene_id <- unlist(gene_id)

    data <- raw_data %>% filter(V3 == 'transcript')

    result_columns <- c(chr_column_name, star_column_name, end_column_name,
                        strand_column_name, gene_id_columns_name)

    if(is.na(chrs)) {

        data_chr <- data  %>%
            select(all_of(result_columns))

    } else {

        if(length(chrs) == 0) {

            data_chr <- data %>%
                select(all_of(result_columns))

        } else if(length(chrs) == 1) {

            data_chr <- data %>% filter(V1 == chrs[1]) %>%
                select(all_of(result_columns))

        } else {

            data_chr <- data %>% filter(V1 %in% chrs) %>%
                select(all_of(result_columns))

        }

    }

    colnames(data_chr) <- c("seqnames", "start", "end", "strand", "gene_id")
    data_chr_gr <- GRanges(data_chr)

    #length(data_chr_gr)
    result <- reduce(data_chr_gr)
    # length(result)
    # length(unique(data_chr$gene_id))

    return(result)

}

get_figure_genes <- function() {

    result <- GRanges(seqnames = "chr14",
                      ranges = IRanges(start = 32934396, end = 33820863 ),
                      strand = "*",
                      hgnc_symbol = "NPAS3" #tiene que ser factor??
                      )
    return(result)
}

get_figure_data <- function(chrs = NA) {

    genome <- BSgenome.Hsapiens.UCSC.hg38

    nc_elems <- get_nc_elements()
    nc_elem_regions <- GRanges(seqnames = nc_elems$seqnames,
                               IRanges(start = nc_elems$start, end = nc_elems$end))

    c_elems <- get_c_elements()
    c_elem_regions <- GRanges(seqnames = c_elems$seqnames,
                               IRanges(start = c_elems$start, end = c_elems$end))

    # chrs <- paste0('chr', seq(1, 22))

    markers <- get_figure_genes()

    #all_genes <- get_all_genes()

    all_genes <- get_all_genes(chrs)

    result <- list(
        'genome' = genome,
        'nc_elems' = nc_elem_regions,
        'c_elems' = c_elem_regions,
        'markers' = markers,
        'genes' = all_genes
    )

    return(result)

}

#kp <- plotKaryotype(chromosomes = chrs)

#kpPlotDensity(kp, data=nc_elem_regions)

#kpPlotRegions(kp, data=nc_elem_regions)


#----

#kp <- plotKaryotype(plot.type=2, chromosomes = chrs)
# ?plotKaryotype
# ?kpPlotRegions
# ? kpPlotMarkers

plot_figure <- function(data, chrs, plot_file_name, zoom = NULL) {

    height_unit <- 0.2
    custom_border <- 0.1

    nc_density_track_base <- 0
    nc_density_track_height <- nc_density_track_base + 1
    nc_elem_track_base <- 0
    nc_elem_track_height <- nc_elem_track_base +  height_unit
    c_elem_track_base <- nc_elem_track_height
    c_elem_track_height <- c_elem_track_base +  height_unit
    genes_track_base <- c_elem_track_height
    genes_track_height <- genes_track_base +  height_unit
    c_density_track_base <- genes_track_height
    c_density_track_height <- c_density_track_base +  2

    plot_params <- getDefaultPlotParams(plot.type=2)
    plot_params$leftmargin <- 0.2
    plot_params$rightmargin <- 0.2
    plot_params$topmargin <- 0.2
    plot_params$bottommargin <- 0.2
    plot_params$data1height <- 200
    plot_params$data2height <- 400
    plot_params$ideogramheight <- 20
    plot_params$data2max <- ceiling(c_density_track_height)

    plot_file_path <- file.path(base_directory, plot_file_name)
    png(plot_file_path, width = 1500, height = 700, units = "px")
    # pdf(plot_file_path, width = 1500, height = 700)


    #kp <- plotKaryotype(plot.type=2, cytobands = GRanges(), chromosomes = "chr14")
    kp <- plotKaryotype(plot.type=2, chromosomes = chrs,
                        plot.params = plot_params,
                        zoom = zoom)

    # kpDataBackground(kp,
    #                  r0=nc_density_track_base, r1=nc_density_track_height,
    #                  col="#FFDDDD")



    #kp <- plotKaryotype(plot.type=2, chromosomes = chrs)
    nc_density_plot <- kpPlotDensity(kp,
                                     data=data$nc_elems,
                                     window.size = 1e6,
                                     col = "#99CAB1",
                                     border = "#013220",
                                     r0=nc_density_track_base,
                                     r1=nc_density_track_height)

    kpAxis(kp,
           ymax = nc_density_plot$latest.plot$computed.values$max.density,
           r0=nc_density_track_base, r1=nc_density_track_height,
           cex=0.8, side = "right")
    kpAbline(kp,
             h = mean(nc_density_plot$latest.plot$computed.values$density),
             lty = 2,
             ymax = nc_density_plot$latest.plot$computed.values$max.density,
             r0=nc_density_track_base, r1=nc_density_track_height)


    kpAddLabels(kp, labels = "density noncod",
                label.margin = 0, side = "left",
                r0=nc_density_track_base, r1=nc_density_track_height, pos = NULL)

    #-----------------------------------

    kpPlotRegions(kp,
                  data = data$nc_elems,
                  data.panel = 2,
                  col = "#99CAB1", border = "#013220",
                  r0 = nc_elem_track_base, r1 = nc_elem_track_height,
    )

    kpAddLabels(kp, labels = "noncod elem regions", data.panel=2,
                label.margin = 0, side = "left",
                r0 = nc_elem_track_base, r1 = nc_elem_track_height, pos = NULL)

    #-----------------------------------

    kpPlotRegions(kp, data = data$c_elems, data.panel = 2,
                  col="#FE6F5E", border = "#FE6F5E",
                  r0=c_elem_track_base, r1=c_elem_track_height - custom_border)

    kpAddLabels(kp, labels = "cod elem regions",
                data.panel=2,
                r0=c_elem_track_base, r1=c_elem_track_height, pos = NULL,
                side = "left")

    #-----------------------------------

    kpPlotRegions(kp, data = data$genes, data.panel=2,
                  col="#72A0C1", border="#205472",
                  r0=genes_track_base, r1=genes_track_height - custom_border)

    kpAddLabels(kp, labels = "genes", data.panel=2,
                r0=genes_track_base, r1=genes_track_height, pos = NULL,
                side = "left")

    #-----------------------------------

    c_density_plot <- kpPlotDensity(kp, data = data$c_elems, window.size = 1e6,
                                    data.panel = 2,
                                    col = "#FF8B77", border = "#FE6F5E",
                                    r0=c_density_track_base, r1=c_density_track_height)

    kpAxis(kp,
           ymax = c_density_plot$latest.plot$computed.values$max.density,
           r0=c_density_track_base, r1=c_density_track_height,
           cex=0.8,
           data.panel=2,
           side = "right")
    kpAbline(kp,
             h = mean(c_density_plot$latest.plot$computed.values$density),
             lty = 2,
             ymax = c_density_plot$latest.plot$computed.values$max.density,
             r0=c_density_track_base, r1=c_density_track_height,
             data.panel=2)


    kpAddLabels(kp, labels = "density cod", data.panel=2,
                r0=c_density_track_base, r1=c_density_track_height, pos = NULL,
                side = "left")

    #-----------------------------------

    kpPlotMarkers(kp, data = data$markers, labels = data$markers$hgnc_symbol,
                  label.margin = 50, text.orientation = "horizontal")

    kpPlotMarkers(kp, data = data$markers, labels = data$markers$hgnc_symbol,
                  data.panel = 2,
                  label.margin = 20, text.orientation = "horizontal",
                  r0=c_density_track_base, r1=c_density_track_height + 5
    )

    #Close the device
    dev.off()

}



plot_chr14 <- function() {

    #?getDefaultPlotParams

    plot_file_name <- "chr14.png"
    # plot_file_name <- "chr14.pdf"

    chrs <- c("chr14")

    data <- get_figure_data(chrs)

    plot_figure(data, chrs, plot_file_name)

}

plot_chrs <- function() {

    #?getDefaultPlotParams

    all_chrs <- paste0("chr", seq(1, 22))

    for(chr in all_chrs) {

        plot_file_name <- paste0(chr, ".png")

        chrs <- c(chr)

        data <- get_figure_data(chrs)

        plot_figure(data, chrs, plot_file_name)

    }

}


plot_NPAS3 <- function() {

    npas3_region <- get_figure_genes()

    chrs <- c("chr14")

    data <- get_figure_data(chrs)

    plot_file_name <- "npas3.png"

    plot_figure(data, chrs, plot_file_name, zoom = npas3_region)

}

# plot_NPAS3()
plot_chr14()
