# Run every locus of the FUMA job end to end, then write the all-loci tables.
# All the logic lives in coloc_lib.R (run_locus, overview_all_loci, supported_genes_all_loci).
#
# From the command line:
#   Rscript code/run_all_loci.R [pan|eur] [nfilter] [only-missing] [tables-only]
#   Rscript code/run_all_loci.R                       # pan-ancestry file, every locus
#   Rscript code/run_all_loci.R eur only-missing      # European-only file, skipping loci already run
#   Rscript code/run_all_loci.R pan nfilter           # pan file with the sample-size filter
#   Rscript code/run_all_loci.R eur tables-only       # just rebuild the two tables from what is on disk
#
# Loci come from FUMA's GenomicRiskLoci.txt. Each locus is run against the fixed dataset panel at
# 500 kb and 250 kb (see run_locus.R). Nothing is downloaded: fetch the artery and nerve regions
# for a locus first with fetch_gtex_v8_regions.R, or those datasets are reported as failed.
#
# Output in results/coloc/ (file names end __eur, __pan_nfilter or __eur_nfilter for those settings):
#   ALL_LOCI__overview_500kb.csv          one row per locus: signal, testable tests, genes by label
#   ALL_LOCI__supported_genes_500kb.csv   one row per gene with a shared call at any locus
#   ALL_LOCI__failures.csv                every dataset run that failed, with the reason

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(script)) dirname(dirname(normalizePath(script[1]))) else getwd()
source(file.path(root, "code", "coloc_lib.R"))

cli <- commandArgs(trailingOnly = TRUE)
known <- c("pan", "eur", "nfilter", "only-missing", "tables-only")
if (length(setdiff(cli, known))) stop("usage: Rscript code/run_all_loci.R [pan|eur] [nfilter] [only-missing] [tables-only]", call. = FALSE)
gwas_version <- if ("eur" %in% cli) "eur" else "pan"
sample_size_filter <- "nfilter" %in% cli
setting <- gwas_label(gwas_version, sample_size_filter)
paths <- project_paths(root)
loci <- locus_table(paths)
all_file <- function(what) file.path(paths$results, sprintf("ALL_LOCI__%s%s.csv", what, gwas_suffix(setting)))

failures <- list()
if (!"tables-only" %in% cli) {
  for (i in seq_len(nrow(loci))) {
    locus <- loci[i]
    if ("only-missing" %in% cli && file.exists(summary_file(paths$results, locus$label, PRIMARY_WINDOW_KB, setting))) {
      message(sprintf("[%d/%d] %s (%s): already run, skipped", i, nrow(loci), locus$label, setting))
      next
    }
    started <- Sys.time()
    done <- tryCatch(suppressMessages(run_locus(locus$label, locus$chr, locus$lead_pos, paths, gwas_version = gwas_version,
                                                sample_size_filter = sample_size_filter)),
                     error = function(e) list(failures = data.table(window_kb = NA_real_, dataset = "WHOLE LOCUS", error = conditionMessage(e))))
    if (nrow(done$failures)) failures[[locus$label]] <- cbind(label = locus$label, done$failures)
    message(sprintf("[%d/%d] %s (%s): done in %.0f s, %d dataset runs failed", i, nrow(loci), locus$label, setting,
                    as.numeric(difftime(Sys.time(), started, units = "secs")), nrow(done$failures)))
  }
  failures <- rbindlist(failures, fill = TRUE)
  if (!nrow(failures)) failures <- data.table(label = character(), window_kb = numeric(), dataset = character(), error = character())
  fwrite(failures, all_file("failures"))
}

overview <- overview_all_loci(paths, loci, setting)
fwrite(overview, all_file(sprintf("overview_%skb", PRIMARY_WINDOW_KB)))
supported <- supported_genes_all_loci(paths, loci, setting)
fwrite(supported, all_file(sprintf("supported_genes_%skb", PRIMARY_WINDOW_KB)))
message(sprintf("\n%s: %d of %d loci run; %d genes with a shared call; tables written to %s", setting, sum(overview$status == "run"),
                nrow(overview), nrow(supported), paths$results))
