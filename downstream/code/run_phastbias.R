# Author: Anabella Trigila
# ------


library(rphast)


run_bgc_analysis <- function(mod_file, maf_file, branch_name, out_dir = "outputs_bgc") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Read MAF file and neutral model
  align <- read.msa(filename = maf_file, format = "MAF")
  neutralMod <- read.tm(mod_file)

  # Classify mutations and run nucleotide tests
  mutation_counts <- classify.muts.bgc(align, neutralMod, branch = branch_name)
  nuc_results <- bgc.nucleotide.tests(align, neutralMod, branch_name)
  print(nuc_results)

  # Write nucleotide test results
  out_file_nuc <- file.path(out_dir, paste0("bgc.nucleotide.tests_results_", branch_name, ".txt"))
  write.table(nuc_results, file = out_file_nuc, sep = "\t", row.names = TRUE, quote = FALSE)

  # Run phastBias (Capra et al. 2013) for GC bias detection
  hmm_results <- phastBias(align, neutralMod, foreground = branch_name)
  print(hmm_results$tracts)

  # Write phastBias tracts
  out_file_hmm <- file.path(out_dir, paste0("hmm_results_track_", branch_name, ".txt"))
  write.table(hmm_results$tracts, file = out_file_hmm, sep = "\t", row.names = TRUE, quote = FALSE)

  return(list(nucleotide_tests = nuc_results, phast_bias_tracks = hmm_results$tracts))
}

# # Example:
# # Mammals analysis:
# mammals_results <- run_bgc_analysis(
#   mod_file = here::here("data", "13. MAFs", "hg38.phastCons100way.mod_branch_mammals.mod"),
#   maf_file = here::here("data", "outputs", "ncMAR-2845.maf"),
#   branch_name = "mammals"
# )
#
# # AVES analysis:
# aves_results <- run_bgc_analysis(
#   mod_file = here::here("data", "13. MAFs", "galGal6.phastCons77way_branch_aves.mod"),
#   maf_file = here::here("data", "13. MAFs", "ncAvARs_maf.maf"),
#   branch_name = "AVES"
# )
