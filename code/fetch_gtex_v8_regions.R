# Download the region around one locus from official GTEx v8 tissue files on the EBI server, so that
# run_coloc.R can use sources such as gtex_v8:Artery_Tibial without the whole 3.6 to 4 GB files.
#
# From RStudio: edit `params`, then click Source.
# From the command line:
#   Rscript code/fetch_gtex_v8_regions.R <chr> <lead_pos_hg19> <window_kb> <Tissue> [<Tissue> ...]
#   Rscript code/fetch_gtex_v8_regions.R 1 50000000 500 Artery_Aorta Artery_Coronary Artery_Tibial Nerve_Tibial
#
# What it does: finds where the hg19 window sits in GRCh38 (by matching the window's SNPs by rsID in
# the local GTEx v8 Brain Cortex file), then downloads that span plus 1 Mb on each side from each
# tissue through its tabix index and adds it to data/eqtl-full/gtex_v8_regions/<Tissue>.hg38.tsv.gz.
# Give the widest window you plan to run; narrower windows are then covered too.
# Only genome coordinates are sent to the server. No study data leaves this machine.
# Tissue names are the file names at ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/.

params <- list(chr = "1", lead_pos = 50000000, window_kb = 500,
               tissues = c("Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Nerve_Tibial"))

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
  if (length(cli) < 4) {
    stop("usage: Rscript code/fetch_gtex_v8_regions.R <chr> <lead_pos_hg19> <window_kb> <Tissue> [<Tissue> ...]", call. = FALSE)
  }
  params <- list(chr = cli[1], lead_pos = as.numeric(cli[2]), window_kb = as.numeric(cli[3]), tissues = cli[-(1:3)])
}

paths <- project_paths(root)
window <- locus_window(params$lead_pos, params$window_kb)
span <- grch38_span_of_window(paths, as.character(params$chr), window$lo, window$hi)
from <- max(1, span[1] - REGION_FETCH_MARGIN_BP)
to <- span[2] + REGION_FETCH_MARGIN_BP
message(sprintf("hg19 %s:%s-%s is GRCh38 %s-%s; fetching %s-%s", params$chr, format_position(window$lo),
                format_position(window$hi), format_position(span[1]), format_position(span[2]),
                format_position(from), format_position(to)))
for (tissue in params$tissues) {
  got <- fetch_gtex_v8_region(tissue, as.character(params$chr), from, to, paths$gtex_v8_regions,
                              hg19_window = c(window$lo, window$hi))
  message(sprintf("%s: %s%s", tissue, got$status, if (is.na(got$rows)) "" else sprintf(", %d rows", got$rows)))
}
