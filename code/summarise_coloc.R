# Build one table per locus and window from every result file run_coloc.R has written for it.
#
# From RStudio: edit `params`, then click Source.
# From the command line:  Rscript code/summarise_coloc.R <label> <window_kb> [gwas setting]
#   The GWAS setting is pan (default), eur, pan_nfilter or eur_nfilter: whichever the runs used.
#   Rscript code/summarise_coloc.R L1_chr1_50.0Mb 250
#
# Output: results/coloc/<label>__summary_<window>kb.csv, one row per gene, replicated genes first.
# The reporting rule (which genes count as replicated) is defined once, in coloc_lib.R.

params <- list(label = "L1_chr1_50.0Mb", window_kb = 500, gwas_setting = "pan")

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
  if (!length(cli) %in% 2:3) stop("usage: Rscript code/summarise_coloc.R <label> <window_kb> [gwas setting]", call. = FALSE)
  params <- list(label = cli[1], window_kb = as.numeric(cli[2]), gwas_setting = if (length(cli) == 3) cli[3] else "pan")
}

paths <- project_paths(root)
results_dir <- paths$results
# The primary-window summary also shows each gene's label at the comparison window, when that was run.
comparison <- NULL
if (params$window_kb == PRIMARY_WINDOW_KB) {
  comparison_exists <- tryCatch({
    summarise_locus(results_dir, params$label, COMPARISON_WINDOW_KB, gwas_setting = params$gwas_setting)
    TRUE
  }, error = function(e) FALSE)
  if (comparison_exists) comparison <- COMPARISON_WINDOW_KB
}
summary_table <- summarise_locus(results_dir, params$label, params$window_kb,
                                 annotation = load_duplication_annotation(paths), comparison_window_kb = comparison,
                                 gwas_setting = params$gwas_setting)
out_file <- summary_file(results_dir, params$label, params$window_kb, params$gwas_setting)
fwrite(summary_table, out_file)
message(sprintf("%d genes across %d datasets -> %s", nrow(summary_table),
                sum(startsWith(names(summary_table), "H4_")), out_file))
shown <- intersect(c("symbol", "gene_id", "support", "support_strict_prior", "n_datasets_shared",
                     "fraction_in_segmental_duplication", "max_H4_primary", "best_primary_dataset"), names(summary_table))
print(summary_table[seq_len(min(15, .N)), ..shown])
