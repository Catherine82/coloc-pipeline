# Run one locus end to end: every dataset at both windows, then the summaries.
# All the logic lives in coloc_lib.R (run_locus); this file only collects the settings.
#
# From RStudio: edit `params`, then click Source.
# From the command line:
#   Rscript code/run_locus.R <label> <chr> <lead_pos_hg19> [pan|eur] [nfilter]
#   Rscript code/run_locus.R L1_chr1_50.0Mb 1 50000000
#   Rscript code/run_locus.R L1_chr1_50.0Mb 1 50000000 eur
#
# Datasets: GTEx v8 cortex, aorta, coronary artery, tibial artery and tibial nerve; PsychENCODE; the
# eight Bryois cell types; and, for comparison only, GTEx v7 cortex and any eQTL Catalogue slice on
# disk that holds the locus. Windows: 500 kb (primary) and 250 kb (comparison).
# Nothing is downloaded: the artery and nerve regions for a locus must already have been fetched
# with fetch_gtex_v8_regions.R, or those datasets are reported as failed with the reason.
#
# Output in results/coloc/: one file per dataset and window, the run logs, and
# <label>__summary_500kb.csv (the primary result) and <label>__summary_250kb.csv. Runs with the
# European file or the sample-size filter end in __eur, __pan_nfilter or __eur_nfilter.

params <- list(
  label = "L1_chr1_50.0Mb",
  chr = "1",
  lead_pos = 50000000,         # hg19 position of the lead SNP (GenomicRiskLoci.txt)
  gwas_version = "pan",        # "pan" or "eur"
  sample_size_filter = FALSE   # TRUE keeps only SNPs whose sample size is within 10% of the lead SNP's
)

script_dir <- function() {
  from_rscript <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(from_rscript)) return(dirname(normalizePath(from_rscript[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    sourced <- sys.frame(i)$ofile
    if (!is.null(sourced)) return(dirname(normalizePath(sourced)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    open_file <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(open_file)) return(dirname(normalizePath(open_file)))
  }
  getwd()
}
root <- script_dir()
while (!file.exists(file.path(root, "code", "coloc_lib.R"))) {
  if (dirname(root) == root) stop("could not find the project folder (the one holding code/coloc_lib.R)", call. = FALSE)
  root <- dirname(root)
}
source(file.path(root, "code", "coloc_lib.R"))

cli <- commandArgs(trailingOnly = TRUE)
if (!interactive()) {
  if (length(cli) < 3 || length(cli) > 5) {
    stop("usage: Rscript code/run_locus.R <label> <chr> <lead_pos_hg19> [pan|eur] [nfilter]", call. = FALSE)
  }
  params <- list(label = cli[1], chr = cli[2], lead_pos = as.numeric(cli[3]),
                 gwas_version = if (length(cli) >= 4) cli[4] else "pan",
                 sample_size_filter = length(cli) >= 5 && cli[5] == "nfilter")
}

done <- run_locus(params$label, params$chr, params$lead_pos, project_paths(root),
                  gwas_version = params$gwas_version, sample_size_filter = params$sample_size_filter)

primary <- fread(done$summaries[[as.character(PRIMARY_WINDOW_KB)]])
message(sprintf("\n%s, %s: %d genes at %s kb; replicated %d, single resource %d -> %s", params$label,
                gwas_label(params$gwas_version, params$sample_size_filter), nrow(primary), PRIMARY_WINDOW_KB,
                sum(primary$support == "replicated"), sum(primary$support == "single resource"),
                done$summaries[[as.character(PRIMARY_WINDOW_KB)]]))
shown <- intersect(c("symbol", "gene_id", "support", "support_strict_prior", "n_datasets_shared",
                     "fraction_in_segmental_duplication", "datasets_two_signals"), names(primary))
print(primary[support != "none", ..shown])
if (nrow(done$failures)) {
  message(sprintf("\n%d dataset runs FAILED:", nrow(done$failures)))
  print(done$failures)
  quit(status = 1)
}
