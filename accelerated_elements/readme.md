
The strategy is to evaluate acceleration in internal branches (mammals and aves independently) on genomic regions conserved in amniotes and mammals or conserved in amniotes and aves to identify shifts that can be functional region candidates.

# Data

## 1. Multiple Species Alignments

At the time of this study (2020) there were no multiple alignments that included a large number of mammals, birds, and other amniotes, so we decided to conduct our search on two large independent individual multiple alignments: one ("100 pathways", <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/>) containing primarily mammals along with other amniotes and another containing primarily birds and other amniotes ("77 pathways", <http://hgdownload.cse.ucsc.edu/goldenPath/galGal6/multiz77way/>).

## 2. Species

### a. Amniotes

The species extracted from the 100-way multiple alignment of vertebrate genomes to compute conservation of amniotes in the mammals' analysis are: hg38 (human, Homo sapiens, as reference species), allMis1 (alligator, Alligator mississippiensis), cheMyd1 (sea turtle, Chelonia mydas), chrPic2 (painted turtle, Chrysemys picta belli), pelSin1 (Chinese turtle, Pelodiscus sinensis), apaSpi1 (spiny turtle, Apalone spinifera), anoCar2 (lizard, Anolis carolinesis), xenTro7 (frog, Xenopus tropicalis) and latCha1 (coelacanth, Latimeria chalumnae)

These species are listed in the configuration file at line 17 of <https://github.com/paulati/cladeAcc/blob/master/inst/config.yml>

The species extracted from the 77-way multiple alignment of vertebrate genomes to compute conservation of amniotes in the aves analysis are: galGal6 (chicken, Gallus gallus as reference species), allMis1 (alligator, Alligator mississippiensis), cheMyd1 (sea turtle, Chelonia mydas), chrPic2 (painted turtle, Chrysemys picta belli), pelSin1 (Chinese turtle, Pelodiscus sinensis), apaSpi (spiny turtle, Apalone spinifera), anoCar2 (lizard, Anolis carolinensis) and xenTro9 (frog, Xenopus tropicalis).

These species are listed in the configuration file at line 19 of <https://github.com/paulati/cladeAcc/blob/master/inst/config.yml>

We would like to remark that there is a core of seven shared amniotes between the 77-way and the 100-way alignment involved in the computation of conserved elements.

### b. Mammals

We extracted data for mammalian species in the publicly available 100-way multiple alignment and sought a balance in the representation of different branches of mammalian diversity. The species selected for the analysis were: hg38 (human, Homo sapiens), otoGar3 (Garnett's galago, Otolemur garnetti), rheMac3 (rhesus monkey, Macaca mulatta), mm10 (mouse, Mus musculus), oryCun2 (rabbit, Oryctolagus cuniculus), ochPri3 (American pika, Ochotona princeps), susScr3 (pig, Sus scrofa), turTru2 (dolphin, Tursiops truncatus), bosTau8 (cow, Bos taurus), felCat8 (cat, Felis catus), myoLuc2 (bat, Myotis lucifugus), loxAfr3 (elephant, Loxodonta africana), echTel2 (hedgehog, Echinops telfairi), dasNov3 (armadillo, Dasypus novemcinctus), monDom5 (gray siskin, Monodelphis domestica), macEug2 (wallaby, Notamacropus eugenii), ornAna1 (platypus, Ornithorhynchus anatinus). The selected species were chosen to represent major mammal groups and adaptations to diverse habitats.

These species are listed in the configuration file at line 16 of <https://github.com/paulati/cladeAcc/blob/master/inst/config.yml>

### c. Aves

We selected the following birds species from the 77-way multiple alignment of vertebrate genomes that includes mainly birds and other amniotes: galGal6 (chicken, Gallus gallus as reference species), cotJap2 (Japanese quail, Coturnix japonica), melGal5 (turkey, Meleagris gallopavo), tytAlb1 (owl, Tyto alba), bucRhi1 (rhinoceros hornbill, Rhinoceros hornbill), anaPla1 (duck, Anas platyrhynchos), apaVit1 (mountain trogon, Apaloderma vittatum), calAnn1 (hummingbird, Calypte anna), cucCan1 (cuckoo, Cuculus canorus), chaVoc2 (plover, Charadrius vociferus), fulGla1 (northern fulmar, Fulmarus glacialis), tauEry1 (red-crested turaco, Tauraco erythrolophus), opiHoa1 (hoatzin, Opisthocomus hoazin), phoRub1 (flamingo, Phoenicopterus roseus), colLiv1 (pigeon, Columba livia), lepDis1 (carraca curol, Leptosomus discolor), merNub1 (crimson honey bee-eater, Merops nubicus), pelCri1 (pelican, Pelecanus crispus), phaCar1 (cormorant, Phalacrocorax carbo), phaLep1 (common yellowtail, Phaethon lepturus), pteGut1 (red-crowned sandgrouse, Pterocles gutturalis), nipNip1 (Japanese ibis, Nipponia nippon), egrGar1 (little egret, Egretta garzetta), pygAde1 (Adelie penguin, Pygoscelis adeliae), aptFor1 (emperor penguin, Aptenodytes forsteri), carCri1 (crested caryama, Cariama cristata), mesUni1 (mesito unicolor, Mesitornis unicolor), eurHel1 (tigana, Eurypyga helias), balPav1 (black-necked crowned crane, Balearica regulorum), chlUnd1 (MacQueen's Houbara, Chlamydotis undulata), falChe1 (saker falcon, Falco cherrug), falPer1 (peregrine falcon, Falco peregrinus), aquChr2 (golden eagle, Aquila chrysaetos), halAlb1 (white-tailed eagle, Haliaeetus albicilla), halLeu1 (bald eagle, Haliaeetus leucocephalus), corBra1 (American crow, Corvus brachyrhynchos), corCor1 (carrion crow, Corvus cornix), acaChl1 (greenish acanthisita, Acanthisitta chloris), ficAlb2 (collared flycatcher, Ficedula albicollis), serCan1 (canary, Serinus canaria), zonAlb1 (white-throated chingolo, Zonotrichia albicollis), geoFor1 (Darwin's finch, Geospiza fortis), taeGut2 (zebra finch, Taeniopygia guttata), pseHum1 (great tit, Pseudopodoces humilis), gavSte1 (little grebe, Gavia stellata), capCar1 (nightjar, Antrostomus carolinensis), melUnd1 (parakeet, Melopsittacus undulatus), amaVit1 (Puerto Rican parakeet, Amazona vittata), araMac1 (scarlet macaw, Ara macao), colStr1 (common mousebird, Colius striatus), picPub1 (lesser downy woodpecker, Picoides pubescens), strCam1 (ostrich, Struthio camelus), tinGut2 (white throated tinamou, Tinamus guttatus).

These species are listed in the configuration file at line 18 of <https://github.com/paulati/cladeAcc/blob/master/inst/config.yml>

## 3. Neutral model

The neutral model for 100-way alignment was downloaded from <https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.mod> This all-species tree model for this track was generated using the phyloFit program from the PHAST package (REV model, EM algorithm, medium precision) using multiple alignments of 4-fold degenerate sites extracted from the 100-way alignment (msa_view). The 4d sites were derived from the RefSeq (Reviewed+Coding) gene set, and filtered to select single-coverage long transcripts.

The neutral model for 77-way alignment was downloaded from <https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/phastCons77way/galGal6.phastCons77way.mod> This all-species tree model was generated using the phyloFit program from the PHAST package (REV model, EM algorithm, medium precision) using multiple alignments of 4-fold degenerate sites extracted from the 77-way alignment (msa_view). The 4d sites were derived from the NCBI RefSeq gene set, and filtered to select single-coverage long transcripts.

We added "MAMMALS" label to Newick tree in hg38.phastCons100way.mod and "AVES" label to the tree in galGal6.phastCons77way.mod to be able to specify the internal branch where acceleration would be evaluated.

Labeling is performed by functions in the file <https://github.com/paulati/cladeAcc/blob/master/R/acceleration_neutral_model_delegates.R> The method that calculates acceleration takes the labeling function as a parameter, allowing these methods to be overridden with custom code in analyses involving a different set of species.

# Pipeline

## 1. Conserved elements computation

### a. Amniotes

We used the phastCons function from the RPHAST package (<https://github.com/CshlSiepelLab/RPHAST>) to identify conserved elements in a set of amniotes. The input alignments used for this calculation are the "100 pathways" and "77 pathways", filtered by the species listed in section 2.a.

First, we identified conserved regions in amniotes using `phastCons` function with the following parameters: expected.length = 45, target.coverage = 0.3, rho = 0. 3, viterbi=TRUE, for each of the two sub-alignments with the species listed in Data 2.a. `phastCons` function returns a list containing parameter estimates, including an object (most.conserved) of type feat which describes conserved elements detected by the Viterbi algorithm.

This computation is performed at line 10 of <https://github.com/paulati/cladeAcc/blob/master/R/conservation.R>

These conserved elements were then filtered according to our criteria for informative regions A region was considered informative when alligator, lizard, and at least a turtle were present in the alignment and had no gaps in these species.

This condition is implemented in the functions `required_species_features_sarcopterygii_100way` and `required_species_features_sarcopterygii_77way` in <https://github.com/paulati/cladeAcc/blob/master/R/conservation_delegates.R>

The method that computes conservation takes required-species-features function as a parameter, allowing it to be overridden with custom code in analyses involving different conditions.

The reported amniotes conserved elements resulted from the intersection of the conserved elements detected by the Viterbi algorithm and the informative regions in the alignment.

This computation is performed at line 15 of <https://github.com/paulati/cladeAcc/blob/master/R/conservation.R>

### b. Mammals

We identified conserved regions in mammals using `phastCons` function with the same parameters values used for amniotes conserved elements computation (expected.length = 45, target.coverage = 0.3, rho = 0. 3, viterbi=TRUE), for a sub-alignment of "100 pathways" containing the species listed in Data 2.b. phastCons `most.conserved` results (conserved elements detected by the Viterbi algorithm) were then filtered according to our criteria for informative regions.

A region was considered informative when the associated alignment fulfilled the following three conditions:

-   hg38, otoGar3, mm10, oryCun2, ochPri3, susScr3, turTru2, bosTau8, felCat8, myoLuc2, ornAna1 were present in the alignment and these species had no gaps,

-   At least one of loxAfr3, echTel2 was present in the alignment and this specia had no gaps,

-   At least one of dasNov3, monDom5, macEug2 was present in the alignment and this specie had no gaps.

These conditions is implemented in the function `required_species_features_mammals_100way` in <https://github.com/paulati/cladeAcc/blob/master/R/conservation_delegates.R>

The method that computes conservation takes required-species-features function as a parameter, allowing it to be overridden with custom code in analyses involving different conditions.

The reported mammals' conserved elements resulted from the intersection of the conserved elements detected by the Viterbi algorithm and the informative regions in the alignment.

This computation is performed at line 15 of <https://github.com/paulati/cladeAcc/blob/master/R/conservation.R>

### c. Aves

We identified conserved regions in aves using `phastCons` function with the same parameters values used for amniotes conserved elements computation (expected.length = 45, target.coverage = 0.3, rho = 0. 3, viterbi=TRUE), for a sub-alignment of "77 pathways" containing the species listed in Data 2.c. phastCons `most.conserved` results (conserved elements detected by the Viterbi algorithm) were then filtered according to our criteria for informative regions. 

A region was considered informative when the associated alignment fulfilled the following four conditions:

-   At least one of falChe1, falPer1, halLeu1 was present in the alignment and this specie had no gaps,
-   At least one of taeGut2, geoFor1 was present in the alignment and this specie had no gaps,
-   At least one of strCam1, tinGut2 was present in the alignment and this specie had no gaps,
-   At least eleven species were present in the alignment and these species had no gaps 

These conditions is implemented in the function `required_species_features_aves_77way` in <https://github.com/paulati/cladeAcc/blob/master/R/conservation_delegates.R>

The method that computes conservation takes required-species-features function as a parameter, allowing it to be overridden with custom code in analyses involving different conditions.

The reported aves' conserved elements resulted from the intersection of the conserved elements detected by the Viterbi algorithm and the informative regions in the alignment.

This computation is performed at line 15 of <https://github.com/paulati/cladeAcc/blob/master/R/conservation.R>


## 2. Accelerated elements computation

The set of candidate regions to evaluate acceleration in mammals or aves was the result of the intersection between the amniotes conserved regions (output from Pipeline 1.a) and the mammals conserved regions (output from Pipeline 1.b.) or aves conserved regions (output from Pipeline 1.c.) 

This intersection is performed by teh function `conserved_elements_in_common` at <https://github.com/paulati/cladeAcc/blob/master/R/conservation.R>

The two intersections of these sets resulted in conserved elements in both amniotes and mammals or amniotes and aves. This intersection is useful for ruling out those elements only conserved in amniotes (which could have been lost in mammals) or only conserved in mammals or aves (probably arisen de novo). 

These candidate regions were then split into 25 bp and 50 bp chunks and acceleration was evaluated in these chunks. 

This step is implemented in the function `split_candidate_regions` in <https://github.com/paulati/cladeAcc/blob/master/R/acceleration.R>

We used the `phyloP` function program from the RPHAST package (<https://github.com/CshlSiepelLab/RPHAST>) with the internal branch test to identify mammalian (parameter branches = "MAMMALS") or avian (parameter branches = "AVES") accelerated elements in a set of chunks mentioned above. 

This step is implemented in the function `compute_observed_phyloP` of <https://github.com/paulati/cladeAcc/blob/master/R/acceleration.R>

The neutral model in `phyloP` function was the same one used for the conserved elements computation. 

To assess the significance of accelerated elements obtained as the result of `phyloP` function, we computed empirical p-values using non-parametric simulations instead of relying on the assumption of a chi-square null distribution. (<http://compgen.cshl.edu/rphast/vignette2.pdf>) 

This step in implemented at <https://github.com/paulati/cladeAcc/blob/master/R/acceleration_non_parametric_stats.R>

We began by extracting the regions from the original alignment that correspond to the conserved elements set (candidate regions to evaluate acceleration, as mentioned above), using the `extract.feature.msa` function. This original alignment for mammals is the "100 pathway" filtered by species listed in sections 2.a and 2.b, while for aves it is the "77 pathway" filtered by the species in sections 2.1 and 2.c.

We then generated 100,000 synthetic alignments by sampling with replacement from this set of regions and ran phyloP on these alignments to obtain our null distribution of log-likelihood ratios, applying phyloP in the same way as for the real data. 

This step in implemented in the function `calculate_non_parametric_distribution_phyloP` of <https://github.com/paulati/cladeAcc/blob/master/R/acceleration_non_parametric_stats.R>

We calculated empirical p-values for the result of phyloP function on real data and adjusted them for multiple comparisons using p.adjust with Benjamini-Hochberg method to compute approximate false discovery rates (FDRs) for each element. 

This step in implemented in the function `calculate_non_parametric_significance` of <https://github.com/paulati/cladeAcc/blob/master/R/acceleration_non_parametric_stats.R>


## 3.  Filtering


To filter the results from the previous step, we implemented a filtering pipeline composed of the following steps:

- Step 1: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step1.R>
 
At line 87, we extracted elements with FDR \< 0.05 and treated them as accelerated elements in mammals or aves and putative functional regions.    

We then used these significative regions to extract the associated MAFs:

- MAFs for ingroup species (mammals / aves)

- MAFs for outgroup species (arcopterygii)

These alignments are required in the next steps of the filtering pipeline.

- Step 2: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step2.R>

We calculated the consensus sequences for the ingroup and outgroup for each MAF representing a significative accelerated region.

- Step 3: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step3.R>

We define a "shift sequence" as the sequence in a MAF associated to the specie where we want to see the shift (ornAna1 / tinGut2)
We define a "shift position" as the location in a MAF where ingroup (same for outgroup) consensus sequence is different of the shift sequence

Based on the consensus sequences obtained in the previous step, we computed these metrics for each region:
- length of the element
- gaps count
- shift count, the number of positions in the MAF that are not shift positions in the ingroup but are shift positions in the outgroup (i.e. for position p, shift sequence[p] == ingroup_consensus[p] and shift sequence[p] != outgroup_consensus[p])
- hamming distance between ingroup consensus sequence and the shift sequence (hamming_in)
- hamming distance between outgroup consensus sequence and the shift sequence (hamming_out)

- Step 4: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step4.R>

Filter elements based on the metrics calculated in the previous step:
The condition is (hamming_in < hamming_out) and (shift count > 0) 
i.e. there are more differences between the outgroup and the shift sequence than between the ingroup and shift sequence.

- Step 5: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step5.R>

We associate each significant accelerated element obtained from the previous step to a conserved element obtained from step 1 of the pipeline (Conserved elements computation).

- Step 6: <https://github.com/paulati/cladeAcc/blob/master/R/custom_filtering_step6.R>

We filter the set of conserved element set obtained in the previous step, keeping only those that overlap with at least one significantly accelerated region.

The result of the filtering process is a set of conserved regions (as defined in Step 1 of the pipeline) that overlaps significant accelerated regions.


## 4.  Reporting

We defined a "candidate functional element" as a conserved region (i.e., the set of candidate regions where acceleration was evaluated, Pipeline) unified by a distance of less than 20 bp from each other and with a minimum size of 100 bp.


# Examples
