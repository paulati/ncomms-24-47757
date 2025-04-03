# references:
# https://yulab-smu.top/treedata-book/chapter2.html
# https://bioconnector.github.io/workshops/r-ggtree.html
# https://cran.r-project.org/web/packages/rphylopic/vignettes/a-getting-started.html

save_images <- function() {

    data_raw <- plot_info()
    data <- data_raw$mammals

    for(i in seq(1, nrow(data))) {

        print(i)

        img <- rphylopic::get_phylopic(uuid = data$phylopic_id[i])
        out_img_path <- file.path(data_base_path, 'images', '100way',
                                  paste0(data$code[i], '.png'))
        # Save image
        rphylopic::save_phylopic(img = img, path = out_img_path)
    }

    data <- data_raw$aves

    for(i in seq(1, nrow(data))) {

        print(i)

        img <- rphylopic::get_phylopic(uuid = data$phylopic_id[i])
        out_img_path <- file.path(data_base_path, 'images', '77way',
                                  paste0(data$code[i], '.png'))
        # Save image
        rphylopic::save_phylopic(img = img, path = out_img_path)
    }



}


plot_info <- function() {

    data_path <- file.path(data_base_path, '100way_labels.csv')
    data_mammals <- read.delim(data_path, sep = '\t', header = TRUE)

    data_path <- file.path(data_base_path, '77way_labels.csv')
    data_aves <- read.delim(data_path, sep = '\t', header = TRUE)

    result <- list('mammals' = data_mammals,
                   'aves' = data_aves)

    return(result)

}


plot_tree_base <- function(clade, nhx_reduced, data_labels) {

    if(clade == 'mammals') {
        xlim_max <- 35
        image_folder_name <- '100way'
        label_ingroup <- 'Mammals'
        label_outgroup <- 'Amniotes'
        clade_offset_ingroup <- 15
        clade_offset_outgroup <- 12
        extend_ingroup <- 8.5
        extend_outgroup <- 12.5
        plot_file_name <- 'tree_mammals'
        font_size <- 5
        image_size <- .06
        ingroup_tree_node <- 7
        outgroup_tree_node <- 12


    } else if(clade == 'aves') {

        xlim_max <- 70
        image_folder_name <- '77way'
        label_ingroup <- 'Aves'
        label_outgroup <- 'Amniotes'
        clade_offset_ingroup <- 30
        clade_offset_outgroup <- 23
        extend_ingroup <- 26
        extend_outgroup <- 30
        plot_file_name <- 'tree_aves'
        font_size <- 3
        image_size <- .02
        ingroup_tree_node <- 23
        outgroup_tree_node <- 19
    }

    p <- ggtree::ggtree(nhx_reduced, ladderize = TRUE, branch.length="none") +
        xlim(NA, xlim_max) +
        ggtree::geom_tiplab(geom="label", offset = 2,
                            fontface = 4, label.size = 0,
                            family = "sans", size = font_size,
                            mapping = aes(
                                color = data_labels[label, 'text_color'],
                                label = data_labels[label, 'name'])
        ) +
        ggtree::scale_color_manual(values=c(gray_code, green_code)) +
        ggtree::geom_tiplab(aes(
            image = file.path(data_base_path, 'images', image_folder_name,
                              paste0(label, '.png'))),
            geom="image", offset=1, align=FALSE, size=image_size,
            color = gray_code) +
        theme(legend.position="none")

    p_ingroup <- p + geom_cladelabel(geom = "label", label.size = 0,
                                  node=ingroup_tree_node,
                                  label = label_ingroup,
                                  offset = clade_offset_ingroup, align=T,
                                  color=black_code, angle = 90,
                                  extend = extend_ingroup,
                                  offset.text = 0.5, fontface = 4)


    p_ingroup_outgroup <- p_ingroup + geom_cladelabel(node = outgroup_tree_node,
                                label = label_outgroup,
                                offset = clade_offset_outgroup, align=T,
                                color=black_code, angle = 90,
                                extend = extend_outgroup,
                                offset.text = 0.5, fontface = 4)

    plot_base_path <- file.path(base_path, 'plots')

    svg_plot_file_name <- paste0(plot_file_name, '.svg')
    pdf_plot_file_name <- paste0(plot_file_name, '.pdf')
    ggplot2::ggsave(file = file.path(plot_base_path, svg_plot_file_name),
                    plot = p_ingroup_outgroup, width=29.7, height=21,
                    units = "cm")
    ggplot2::ggsave(file = file.path(plot_base_path, pdf_plot_file_name),
                    plot = p_ingroup_outgroup, width=29.7, height=21,
                    units = "cm")

}



plot_tree_aves <- function() {

    nwk_file_name <- 'galGal6.phastCons77way.labeled.mod'
    nwk_file_path <- file.path(data_base_path, nwk_file_name)
    file.exists(nwk_file_path)

    tree <- ape::read.tree(nwk_file_path)
    all_species <- tree$tip.label

    nodes_77way_aves <- c('galGal6', 'cotJap2', 'melGal5', 'tytAlb1', 'bucRhi1',
                          'anaPla1', 'apaVit1', 'calAnn1', 'cucCan1', 'chaVoc2',
                          'fulGla1', 'tauEry1', 'opiHoa1', 'phoRub1', 'colLiv1',
                          'lepDis1', 'merNub1', 'pelCri1', 'phaCar1', 'phaLep1',
                          'pteGut1', 'nipNip1', 'egrGar1', 'pygAde1', 'aptFor1',
                          'carCri1', 'mesUni1', 'eurHel1', 'balPav1', 'chlUnd1',
                          'falChe1', 'falPer1', 'aquChr2', 'halAlb1', 'halLeu1',
                          'corBra1', 'corCor1', 'acaChl1', 'ficAlb2', 'serCan1',
                          'zonAlb1', 'geoFor1', 'taeGut2', 'pseHum1', 'gavSte1',
                          'capCar1', 'melUnd1', 'amaVit1', 'araMac1', 'colStr1',
                          'picPub1', 'strCam1', 'tinGut2')
    nodes_77way_sarcopterygii <- c('allMis1', 'cheMyd1', 'chrPic2', 'pelSin1',
                                   'apaSpi1', 'anoCar2', 'xenTro9')
    nodes_77way <- c(nodes_77way_aves, nodes_77way_sarcopterygii)

    length(nodes_77way)

    to_drop <- setdiff(all_species, nodes_77way)
    length(to_drop)

    nhx_reduced <- treeio::drop.tip(tree, to_drop)

    data_labels <- plot_info()
    data_labels_77_way <- data_labels$aves


    text_color <- unlist(lapply(data_labels_77_way$code, function(x) {
        if(x %in% nodes_77way_aves) {
            result <- green_code
        } else {
            result <- black_code
        }

    }))
    data_labels_77_way$text_color <- as.factor(text_color)
    data_labels_77_way$image_path <- file.path(data_base_path, 'images',
                            '77way', paste0(data_labels_77_way$code, '.png'))

    rownames(data_labels_77_way) <- data_labels_77_way$code


    plot_tree_base('aves', nhx_reduced, data_labels_77_way)


}


plot_tree_mammals <- function() {

    nwk_file_name <- 'hg38.phastCons100way.labeled.mod'
    nwk_file_path <- file.path(data_base_path, nwk_file_name)
    file.exists(nwk_file_path)

    tree <- ape::read.tree(nwk_file_path)
    all_species <- tree$tip.label


    nodes_100way_mammals <- c('hg38', 'otoGar3', 'rheMac3', 'mm10', 'oryCun2',
                              'ochPri3', 'susScr3', 'turTru2', 'bosTau8',
                              'felCat8', 'myoLuc2', 'loxAfr3', 'echTel2',
                              'dasNov3', 'monDom5', 'macEug2', 'ornAna1')
    nodes_100way_sarcopterygii <- c('allMis1', 'cheMyd1', 'chrPic2', 'pelSin1',
                                    'apaSpi1', 'anoCar2', 'xenTro7', 'latCha1')
    nodes_100way <- c(nodes_100way_mammals, nodes_100way_sarcopterygii)

    length(nodes_100way)

    to_drop <- setdiff(all_species, nodes_100way)
    length(to_drop)

    nhx_reduced <- treeio::drop.tip(tree, to_drop)
    length(nhx_reduced$node.label)
    length(nhx_reduced$tip.label)


    data_labels <- plot_info()
    data_labels_100_way <- data_labels$mammals


    text_color <- unlist(lapply(data_labels_100_way$code, function(x) {
        if(x %in% nodes_100way_mammals) {
            result <- green_code
        } else {
            result <- black_code
        }

    }))
    data_labels_100_way$text_color <- as.factor(text_color)
    data_labels_100_way$image_path <- file.path(data_base_path, 'images',
                        '100way', paste0(data_labels_100_way$code, '.png'))

    rownames(data_labels_100_way) <- data_labels_100_way$code


    plot_tree_base('mammals', nhx_reduced, data_labels_100_way)
}


main <- function() {

    base_path <<- '/u01/home/pbeati/2024/lucia/paper_acelerados/ncomms-24-47757/phylogenetic_tree'
    data_base_path <<- file.path(base_path, 'data')

    # save_images()

    green_code <<- '#558658'
    gray_code <<- '#808080'
    black_code <<- '#262626'

    plot_tree_aves()

    plot_tree_mammals()

}

