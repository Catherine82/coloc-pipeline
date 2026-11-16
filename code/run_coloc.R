# Run colocalization (coloc.abf) for one GWAS locus against one eQTL source.
# All the logic lives in coloc_lib.R; this file only collects the settings and calls run_coloc().
#
# From RStudio: edit the `params` block below, then click Source.
#
# From the command line:
#   Rscript code/run_coloc.R <label> <chr> <lead_pos_hg19> <window_kb> <source> [out.csv] [dataset_label] [pan|eur] [nfilter]
#   Rscript code/run_coloc.R L1_chr1_50.0Mb 1 50000000 250 gtex_v7
#   Rscript code/run_coloc.R L1_chr1_50.0Mb 1 50000000 250 bryois:Excitatory.neurons
#   Rscript code/run_coloc.R L1_chr1_50.0Mb 1 50000000 500 \
#     catalogue:data/eqtl-catalogue/QTD000171_brain_cortex_chr16_27.3-29.9Mb_hg38.tsv "" GTEx_brain_cortex
#
# <source> is one of:
#   gtex_v7              GTEx v7 Brain Cortex, the dataset ezQTL uses
#   gtex_v8              GTEx v8 Brain Cortex, GTEx's own results (GRCh38 file, read with tabix)
#   gtex_v8:<Tissue>     another GTEx v8 tissue, for example gtex_v8:Artery_Tibial; its region for
#                        this locus must first be downloaded with fetch_gtex_v8_regions.R
#   psychencode          PsychENCODE full cis-eQTL set
#   bryois:<CellType>    Bryois et al. 2022 single-nucleus eQTL (file <CellType>.<chr>.gz)
#   catalogue:<file>     a region slice of an eQTL Catalogue file; needs a dataset label
#
# Output: results/coloc/<label>__<dataset>__<window>kb.csv unless out.csv is given, plus
# run_logs/<same name>.skipped.csv (every gene not tested, and why) and
# run_logs/<same name>.settings.tsv (every setting, count and approximation behind the numbers).
#
# The GWAS summary statistics are read from disk and never leave this machine.

# ---- Settings used when this file is sourced in RStudio ----
params <- list(
  label = "L1_chr1_50.0Mb",      # free text, used in file names
  chr = "1",
  lead_pos = 50000000,        # hg19 position of the lead SNP (GenomicRiskLoci.txt)
  window_kb = 500,            # half-width in kb; 500 is the primary window, 250 the ezQTL comparison
  source = "gtex_v7",
  out_file = NULL,            # NULL = default name under results/coloc/
  dataset_label = NULL,       # NULL = default for the source; required for catalogue:<file>
  gwas_version = "pan",       # "pan" (pan-ancestry file) or "eur" (European-only file)
  sample_size_filter = FALSE  # TRUE keeps only SNPs whose sample size is within 10% of the lead SNP's
)

# ---- Find the project folder, whichever way this file was started ----
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

find_project_root <- function(start) {
  dir <- normalizePath(start)
  repeat {
    if (file.exists(file.path(dir, "code", "coloc_lib.R"))) return(dir)
    parent <- dirname(dir)
    if (parent == dir) stop("could not find the project folder (the one holding code/coloc_lib.R)", call. = FALSE)
    dir <- parent
  }
}

root <- find_project_root(script_dir())
source(file.path(root, "code", "coloc_lib.R"))

# ---- Command-line arguments replace the block above ----
cli <- commandArgs(trailingOnly = TRUE)
if (!interactive()) {
  if (length(cli) < 5 || length(cli) > 9) {
    stop("usage: Rscript code/run_coloc.R <label> <chr> <lead_pos_hg19> <window_kb> <source> [out.csv] [dataset_label] [pan|eur] [nfilter]",
         call. = FALSE)
  }
  params <- list(label = cli[1], chr = cli[2], lead_pos = as.numeric(cli[3]), window_kb = as.numeric(cli[4]),
                 source = cli[5], out_file = if (length(cli) >= 6) cli[6] else NULL,
                 dataset_label = if (length(cli) >= 7) cli[7] else NULL,
                 gwas_version = if (length(cli) >= 8) cli[8] else "pan",
                 sample_size_filter = length(cli) >= 9 && cli[9] == "nfilter")
}

run <- run_coloc(label = params$label, chr = params$chr, lead_pos = params$lead_pos,
                 window_kb = params$window_kb, source = params$source, out_file = params$out_file,
                 dataset_label = params$dataset_label, paths = project_paths(root),
                 gwas_version = params$gwas_version, sample_size_filter = params$sample_size_filter)

if (nrow(run$results)) {
  print(run$results[seq_len(min(12, .N)),
                    .(symbol, gene_id, n_snps, min_p_eqtl = signif(min_p_eqtl, 3),
                      H3 = round(PP.H3_two_signals, 3), H4 = round(PP.H4_shared, 3),
                      H4_strict_prior = round(PP.H4_shared_strict_prior, 3), outcome, top_snp)])
}
if (nrow(run$skipped)) {
  message(sprintf("%d genes skipped; reasons:", nrow(run$skipped)))
  print(run$skipped[, .N, by = reason])
}
