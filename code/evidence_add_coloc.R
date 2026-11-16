# Add colocalization columns to the gene evidence table for one locus that has been run with run_locus.R.
#
# From the command line:
#   Rscript code/evidence_add_coloc.R <label> <FUMA locus number>
#   Rscript code/evidence_add_coloc.R L1_chr1_50.0Mb 1
#
# Reads the gene evidence table named in the study settings and the locus's primary summary
# (results/coloc/<label>__summary_500kb.csv, pan-ancestry file), plus the European-only and
# sample-size-filtered summaries when they exist. Fills these columns for the genes of that locus
# and leaves every other column, and every other locus's rows, exactly as they were:
#   coloc_status                       "tested", "not tested in any dataset", or "locus not run"
#   coloc_support, coloc_support_strict_prior    the reporting rule under each prior (500 kb, pan file)
#   coloc_n_datasets_shared, coloc_datasets_shared, coloc_datasets_two_signals
#   coloc_n_datasets_testable          primary datasets in which the gene could be tested
#   coloc_max_H4, coloc_best_dataset   highest PP.H4 in a primary dataset, and where
#   coloc_support_at_250kb, coloc_support_european_only, coloc_support_sample_size_filter
#   fraction_in_segmental_duplication
# The existing evidence counts (evidence_lines, functional_lines) are not changed.

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(script)) dirname(dirname(normalizePath(script[1]))) else getwd()
source(file.path(root, "code", "coloc_lib.R"))

cli <- commandArgs(trailingOnly = TRUE)
if (length(cli) != 2) stop("usage: Rscript code/evidence_add_coloc.R <label> <FUMA locus number>", call. = FALSE)
label <- cli[1]
locus_number <- as.character(cli[2])

paths <- project_paths(root)
evidence_file <- file.path(root, "results", STUDY_FILES$gene_evidence)
evidence <- fread(evidence_file, colClasses = "character", na.strings = NULL)
original_columns <- names(evidence)[!startsWith(names(evidence), "coloc_") & names(evidence) != "fraction_in_segmental_duplication"]

primary <- fread(summary_file(paths$results, label, PRIMARY_WINDOW_KB, "pan"))
other_setting <- function(setting) {
  file <- summary_file(paths$results, label, PRIMARY_WINDOW_KB, setting)
  if (!file.exists(file)) return(NULL)
  fread(file)[, .(gene_id, support)]
}
european <- other_setting("eur")
filtered <- other_setting("pan_nfilter")
label_from <- function(table, ids, missing_file = "not run", missing_gene = "not tested") {
  if (is.null(table)) return(rep(missing_file, length(ids)))
  value <- table$support[match(ids, table$gene_id)]
  ifelse(is.na(value), missing_gene, value)
}

new_columns <- c("coloc_status", "coloc_support", "coloc_support_strict_prior", "coloc_n_datasets_shared",
                 "coloc_datasets_shared", "coloc_datasets_two_signals", "coloc_n_datasets_testable", "coloc_max_H4",
                 "coloc_best_dataset", "coloc_support_at_250kb", "coloc_support_european_only",
                 "coloc_support_sample_size_filter", "fraction_in_segmental_duplication")
for (column in new_columns) if (!column %in% names(evidence)) evidence[, (column) := ""]
evidence[coloc_status == "", coloc_status := "locus not run"]

rows <- which(evidence$locus == locus_number)
if (!length(rows)) stop(sprintf("no genes with locus %s in the evidence table", locus_number), call. = FALSE)
hit <- match(evidence$ensg[rows], primary$gene_id)
tested <- !is.na(hit)
text <- function(x) ifelse(is.na(x), "", as.character(x))
at_250 <- paste0("support_at_", COMPARISON_WINDOW_KB, "kb")
set(evidence, rows, "coloc_status", ifelse(tested, "tested", "not tested in any dataset"))
set(evidence, rows, "coloc_support", text(primary$support[hit]))
set(evidence, rows, "coloc_support_strict_prior", text(primary$support_strict_prior[hit]))
set(evidence, rows, "coloc_n_datasets_shared", text(primary$n_datasets_shared[hit]))
set(evidence, rows, "coloc_datasets_shared", text(primary$datasets_at_threshold[hit]))
set(evidence, rows, "coloc_datasets_two_signals", text(primary$datasets_two_signals[hit]))
set(evidence, rows, "coloc_n_datasets_testable", text(primary$n_primary_datasets_testable[hit]))
set(evidence, rows, "coloc_max_H4", text(round(primary$max_H4_primary[hit], 3)))
set(evidence, rows, "coloc_best_dataset", text(primary$best_primary_dataset[hit]))
set(evidence, rows, "coloc_support_at_250kb", if (at_250 %in% names(primary)) text(primary[[at_250]][hit]) else "not run")
set(evidence, rows, "coloc_support_european_only", ifelse(tested, label_from(european, evidence$ensg[rows]), ""))
set(evidence, rows, "coloc_support_sample_size_filter", ifelse(tested, label_from(filtered, evidence$ensg[rows]), ""))
set(evidence, rows, "fraction_in_segmental_duplication",
    if ("fraction_in_segmental_duplication" %in% names(primary)) text(primary$fraction_in_segmental_duplication[hit]) else "")

setcolorder(evidence, c(original_columns, new_columns))
fwrite(evidence, evidence_file)
message(sprintf("%s (locus %s): %d genes in the evidence table, %d tested by colocalization, %d not tested in any dataset -> %s",
                label, locus_number, length(rows), sum(tested), sum(!tested), evidence_file))
print(evidence[rows][coloc_support %in% c("replicated", "single resource"),
                     .(symbol, ensg, functional_lines, coloc_support, coloc_support_strict_prior, coloc_n_datasets_shared,
                       coloc_support_european_only, fraction_in_segmental_duplication)])
