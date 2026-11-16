# Functions for colocalization of one GWAS locus against one eQTL dataset.
#
# Method: coloc.abf (Giambartolomei et al. 2014; R package coloc), one gene at a time, which assumes
# at most one causal variant per trait in the window.
#
# This file only defines functions and constants. Run an analysis with run_coloc.R; run the tests
# with `Rscript code/tests/test_coloc_lib.R`.
#
# Every eQTL reader returns the same thing, so the matching and coloc code is shared:
#   $data      data.table with columns gene_id, gene_name, pos, rsid, ak, beta, se, p, maf
#   $spec      how to use it: match_by, n (sample size) or sd_y, approximations (plain text)
#   $counts    named integer vector: rows read and rows dropped by each filter
#
# coloc.abf uses beta^2 / variance only, so which allele an effect refers to does not matter; SNPs
# are matched on an unordered allele pair ("ak") wherever the eQTL source reports alleles.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# One data.table thread, on purpose. With worker threads, R can deadlock at exit when its output is
# piped into something that stops reading early (for example `Rscript code/run_coloc.R ... | head -1`):
# the broken-pipe signal lands on a worker thread and R never finishes shutting down. Seen twice on
# November 13, 2026. Run time is the same either way (it is dominated by unzipping the GWAS file).
setDTthreads(1)

# ---- Study settings -----------------------------------------------------------------------------

# Everything that belongs to one study (sample sizes, GWAS file names, the FUMA job folder, fixed
# locus labels) lives in code/study_config.R, which git does not track. To set one up, copy
# code/study_config.example.R to code/study_config.R and fill it in. It is read here when present;
# the environment variable COLOC_STUDY_CONFIG can point at a different file.
STUDY <- NULL                          # list(n_cases, n_controls) of the primary GWAS
GWAS_VERSIONS <- NULL                  # per GWAS version: file name, n_cases, n_controls
LOCUS_LABEL_OVERRIDES <- character()   # optional fixed file-name labels, named by FUMA locus number
STUDY_FILES <- NULL                    # study-specific file and folder names
.study_config <- Sys.getenv("COLOC_STUDY_CONFIG", "")
if (!nzchar(.study_config)) {
  for (.frame in rev(seq_len(sys.nframe()))) {       # the folder this file is being sourced from
    .this_file <- sys.frame(.frame)$ofile
    if (!is.null(.this_file)) {
      .study_config <- file.path(dirname(normalizePath(.this_file)), "study_config.R")
      break
    }
  }
}
if (nzchar(.study_config) && file.exists(.study_config)) source(.study_config, local = TRUE)
suppressWarnings(rm(.study_config, .frame, .this_file))

#' Stop with instructions when the study settings have not been loaded.
need_study_config <- function() {
  if (is.null(STUDY) || is.null(GWAS_VERSIONS) || is.null(STUDY_FILES)) {
    stop("study settings are not set: copy code/study_config.example.R to code/study_config.R and fill it in", call. = FALSE)
  }
}

# Optional sample-size filter: keep only SNPs whose sample size is within this share of the lead
# SNP's. In the pan file SNPs measured in only some ancestry groups have smaller sample sizes.
SAMPLE_SIZE_TOLERANCE <- 0.10

MIN_SNPS_DEFAULT <- 50

GWAS_COLUMNS <- c("CHR", "BP", "A1", "A2", "rsID", "P", "Beta", "SE", "N")

GTEX_V7_COLUMNS <- c("gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf",
                     "pval_nominal", "slope", "slope_se")

# QTLtools nominal layout, as documented in data/eqtl-full/README.md. The file has no header.
PSYCHENCODE_COLUMNS <- c("gene", "gene_chr", "gene_start", "gene_end", "strand", "n_variants_tested",
                         "distance", "variant", "variant_chr", "variant_start", "variant_end",
                         "p", "slope", "top_variant")

BRYOIS_COLUMNS <- c("gene", "rsid", "tss_distance", "p", "beta")
BRYOIS_SNP_COLUMNS <- c("SNP", "SNP_id_hg38", "SNP_id_hg19", "effect_allele", "other_allele", "MAF")

# eQTL Catalogue "all associations" layout. Region slices made with tabix have no header.
CATALOGUE_COLUMNS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
                       "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
                       "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

# Official GTEx v8 all-pairs file as imported (not re-processed) by the eQTL Catalogue. Has a header.
GTEX_V8_COLUMNS <- c("variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id", "maf",
                     "gene_id", "median_tpm", "beta", "se", "an", "ac", "chromosome", "position", "ref",
                     "alt", "type", "rsid")

# GRCh38 files are read over the hg19 window plus this much on each side, because the same SNP can
# sit up to a few Mb away in the other build. SNPs are then matched by rsID, so the extra rows are
# harmless; run_coloc() stops if the matched SNPs come within MIN_EDGE_MARGIN_BP of the region's edge.
BUILD_PAD_MB_DEFAULT <- 5
MIN_EDGE_MARGIN_BP <- 100000

# Where the official GTEx v8 files for other tissues live, and how far past the matched SNPs a
# downloaded region reaches on each side (see fetch_gtex_v8_region).
GTEX_V8_TISSUE_URL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/%s.tsv.gz"
REGION_FETCH_MARGIN_BP <- 1e6

# coloc priors. p1 and p2 stay at coloc's defaults (1e-4 each). The prior that one SNP drives both
# traits (p12) is run two ways in every analysis: coloc's default, and the stricter value from
# Wallace 2020 (PLoS Genet, doi 10.1371/journal.pgen.1008720): "the coloc default of p12 = 10^-5 may
# be overly liberal ... p12 = 5 x 10^-6 may be a more generally robust choice".
COLOC_P12_DEFAULT <- 1e-5
COLOC_P12_STRICT <- 5e-6

# Reporting rule, fixed before any locus other than the first was run.
# A gene is "replicated" when PP.H4 >= SUPPORT_THRESHOLD in two or more independent resources,
# "single resource" when in one, otherwise "none". Datasets from the same people count once:
# every GTEx tissue and release is one resource, and the Bryois cell types are one resource.
# The datasets below are sensitivity checks and never count towards the rule.
# GTEx v7 cortex joined that list on November 14, 2026: it is the earlier release of the same GTEx
# collection as v8 and is kept only for comparing with ezQTL, which uses it.
SUPPORT_THRESHOLD <- 0.8
SENSITIVITY_ONLY_DATASETS <- c("GTEx_v8_brain_cortex_reprocessed", "CommonMind_DLPFC", "ROSMAP_DLPFC",
                               "BrainSeq_DLPFC", "GTEx_v7_brain_cortex")

# What makes a gene test "testable", fixed before the first locus was rerun under them:
#   1. the GWAS file has a migraine signal at the locus: the smallest p within SIGNAL_SEARCH_KB of
#      the lead SNP is below GWAS_SIGNAL_P;
#   2. the gene has an expression signal in that dataset: the smallest eQTL p among the tested SNPs
#      is below EXPRESSION_SIGNAL_P;
#   3. the tested SNPs keep at least SIGNAL_KEPT_MIN of the migraine signal, on the log10 p scale.
# A test that fails one of these is "not testable" with the reason, never a negative result.
GWAS_SIGNAL_P <- 1e-5
SIGNAL_SEARCH_KB <- 250
EXPRESSION_SIGNAL_P <- 1e-5
SIGNAL_KEPT_MIN <- 0.9
# A testable test is "shared" (PP.H4 >= SUPPORT_THRESHOLD), "two signals" (PP.H3 >= TWO_SIGNALS_MIN),
# "leaning" (PP.H4 >= LEANING_MIN) or "undecided".
TWO_SIGNALS_MIN <- 0.8
LEANING_MIN <- 0.5

# File-name labels are L<locus>_chr<chr>_<lead position in Mb>Mb (see locus_table). A study can fix
# the label of a locus in its settings (LOCUS_LABEL_OVERRIDES), so files named earlier stay valid.

# The datasets a locus is run against (see locus_sources).
# Three arteries, tibial nerve, two adipose depots and pituitary are run next to brain. Their
# regions must be downloaded first (fetch_gtex_v8_regions.R).
GTEX_V8_REGION_TISSUES <- c("Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Nerve_Tibial",
                            "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Pituitary")
BRYOIS_CELL_TYPES <- c("Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons",
                       "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes")
CATALOGUE_SLICE_LABELS <- c(QTD000171 = "GTEx_v8_brain_cortex_reprocessed", QTD000075 = "CommonMind_DLPFC",
                            QTD000434 = "ROSMAP_DLPFC", QTD000051 = "BrainSeq_DLPFC")

# Window (lead SNP plus and minus this many kb). The primary window is the one results are reported
# from; the comparison window is ezQTL's. A FUMA locus that reaches further from its lead SNP than
# the window is cut by it; every run log records whether the window covers the FUMA locus.
PRIMARY_WINDOW_KB <- 500
COMPARISON_WINDOW_KB <- 250

# UCSC hg19 tables used to flag genes inside duplicated sequence (layouts from the .sql files at
# hgdownload.soe.ucsc.edu/goldenPath/hg19/database/). Neither file has a header.
UCSC_SEGDUP_COLUMNS <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "otherChrom",
                         "otherStart", "otherEnd", "otherSize", "uid", "posBasesHit", "testResult", "verdict",
                         "chits", "ccov", "alignfile", "alignL", "indelN", "indelS", "alignB", "matchB",
                         "mismatchB", "transitionsB", "transversionsB", "fracMatch", "fracMatchIndel", "jcK", "k2K")
UCSC_ENSGENE_COLUMNS <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
                          "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

# Sample sizes that the files themselves do not carry.
# Wang et al. 2018 (Science, doi 10.1126/science.aat8464, PMC6413328), Figure 1 legend: 1362 PsychENCODE + 25 GTEx individuals
# = 1387 with genotypes matched to RNA-seq for QTL analysis; same number in the Gerstein lab FAQ.
PSYCHENCODE_N <- 1387
BRYOIS_N <- 192         # Bryois et al. 2022.

MATCH_MODES <- c("pos_alleles", "pos", "rsid_alleles")

# ---- Small helpers ------------------------------------------------------------------------------

#' Unordered allele pair, so A/G and G/A give the same key.
allele_key <- function(a1, a2) {
  ifelse(a1 < a2, paste(a1, a2), paste(a2, a1))
}

#' Standard error implied by a two-sided p-value and an effect size.
se_from_p <- function(beta, p) {
  abs(beta) / abs(qnorm(p / 2))
}

#' Stop with a readable message unless a table has exactly the expected column names.
check_columns <- function(actual, expected, what) {
  if (!identical(as.character(actual), expected)) {
    stop(sprintf("%s: unexpected columns.\n  expected: %s\n  found:    %s", what,
                 paste(expected, collapse = ", "), paste(actual, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless a file exists.
check_file <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path), call. = FALSE)
  invisible(TRUE)
}

#' Whole number as plain digits (never 1e+05), for building shell commands.
format_position <- function(x) {
  format(round(x), scientific = FALSE, trim = TRUE)
}

#' Shell command that streams a file to stdout, decompressing if it ends in .gz.
stream_command <- function(path) {
  if (grepl("\\.gz$", path)) paste("gzip -dc", shQuote(path)) else paste("cat", shQuote(path))
}

#' fread() on a command, returning an empty table (not a warning) when the command prints nothing.
fread_allow_empty <- function(cmd, ...) {
  withCallingHandlers(fread(cmd = cmd, ...), warning = function(w) {
    if (grepl("has size 0", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

#' Path to tabix (htslib). RStudio started from the Dock does not always have Homebrew on its PATH.
find_tabix <- function() {
  candidates <- c(Sys.which("tabix"), "/opt/homebrew/bin/tabix", "/usr/local/bin/tabix")
  found <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(found)) stop("tabix not found; install htslib (brew install htslib)", call. = FALSE)
  unname(found[1])
}

#' Window around a lead SNP, as whole base pairs.
locus_window <- function(lead_pos, window_kb) {
  if (!is.numeric(lead_pos) || length(lead_pos) != 1 || is.na(lead_pos) || lead_pos <= 0) {
    stop("lead_pos must be one positive number", call. = FALSE)
  }
  if (!is.numeric(window_kb) || length(window_kb) != 1 || is.na(window_kb) || window_kb <= 0) {
    stop("window_kb must be one positive number", call. = FALSE)
  }
  half <- round(window_kb * 1000)
  list(lo = max(0, round(lead_pos) - half), hi = round(lead_pos) + half)
}

#' Empty eQTL table in the standard layout.
empty_eqtl <- function() {
  data.table(gene_id = character(), gene_name = character(), pos = integer(), rsid = character(),
             ak = character(), beta = numeric(), se = numeric(), p = numeric(), maf = numeric())
}

# Small in-memory store for the last large file read (see read_eqtl_bryois).
READ_CACHE <- new.env()

# ---- GWAS ---------------------------------------------------------------------------------------

#' Read the GWAS rows inside a window (hg19). The file is about 1 GB compressed, so it is filtered
#' in a stream and never loaded whole. Nothing here sends data anywhere.
#'
#' @return list(data = data.table(pos, rsid, ak, gwas_beta, gwas_se, gwas_p, gwas_n), counts = ...)
read_gwas_window <- function(gwas_file, chr, lo, hi) {
  check_file(gwas_file, "GWAS file")
  cmd <- sprintf("%s | awk -F'\\t' -v c=%s -v lo=%s -v hi=%s 'NR==1 || ($1==c && $2>=lo && $2<=hi)'",
                 stream_command(gwas_file), shQuote(as.character(chr)),
                 format_position(lo), format_position(hi))
  raw <- fread(cmd = cmd, colClasses = list(character = c("CHR", "A1", "A2", "rsID")))
  check_columns(names(raw), GWAS_COLUMNS, "GWAS file")
  usable <- raw[!is.na(Beta) & !is.na(SE) & SE > 0]
  data <- usable[, .(pos = as.integer(BP), rsid = rsID, ak = allele_key(A1, A2),
                     gwas_beta = Beta, gwas_se = SE, gwas_p = P, gwas_n = N)]
  list(data = data,
       counts = c(gwas_rows_in_window = nrow(raw), gwas_rows_with_beta_and_se = nrow(usable)))
}

#' The file for a GWAS version. `paths$gwas` is a named vector (pan = ..., eur = ...); a single
#' unnamed path is taken as the pan file.
gwas_file_for <- function(paths, version) {
  if (is.null(GWAS_VERSIONS)) need_study_config()
  if (!version %in% names(GWAS_VERSIONS)) {
    stop(sprintf("unknown GWAS version '%s' (use %s)", version, paste(names(GWAS_VERSIONS), collapse = " or ")), call. = FALSE)
  }
  if (!is.null(names(paths$gwas))) {
    if (!version %in% names(paths$gwas)) stop(sprintf("no file given for GWAS version '%s'", version), call. = FALSE)
    return(paths$gwas[[version]])
  }
  if (version == "pan") return(paths$gwas)
  stop(sprintf("no file given for GWAS version '%s'", version), call. = FALSE)
}

#' Name of a GWAS setting ("pan", "eur", "pan_nfilter", "eur_nfilter") and the file-name ending it
#' gives results. The pan file without the filter has no ending, so those file names stay as before.
gwas_label <- function(version, sample_size_filter = FALSE) {
  paste0(version, if (isTRUE(sample_size_filter)) "_nfilter" else "")
}
gwas_suffix <- function(label) {
  if (label == "pan") "" else paste0("__", label)
}

#' Keep SNPs whose sample size is within `tolerance` of the lead SNP's (of the smallest-p SNP near
#' the lead position when the lead itself is not in the file).
#' @return list(data, reference_n, dropped)
filter_sample_size <- function(gwas, lead_pos, tolerance = SAMPLE_SIZE_TOLERANCE) {
  reference <- gwas[pos == round(lead_pos), gwas_n][1]
  if (is.na(reference)) {
    near <- gwas[abs(pos - round(lead_pos)) <= SIGNAL_SEARCH_KB * 1000]
    reference <- near[which.min(gwas_p), gwas_n][1]
  }
  if (is.na(reference)) stop("sample-size filter: no SNP near the lead position to take a reference from", call. = FALSE)
  keep <- !is.na(gwas$gwas_n) & abs(gwas$gwas_n - reference) <= tolerance * reference
  list(data = gwas[keep], reference_n = reference, dropped = sum(!keep))
}

#' The migraine signal at the locus in this GWAS file: the smallest p within SIGNAL_SEARCH_KB of
#' the lead position (NA when the file has no SNP there).
migraine_reference_p <- function(gwas, lead_pos) {
  near <- gwas[abs(pos - round(lead_pos)) <= SIGNAL_SEARCH_KB * 1000, gwas_p]
  if (!length(near)) NA_real_ else min(near)
}

# ---- eQTL readers -------------------------------------------------------------------------------

#' Rows of the pre-sliced GTEx v7 Brain Cortex file for one chromosome and window, unfiltered.
read_gtex_v7_rows <- function(file, chr, lo, hi) {
  check_file(file, "GTEx v7 slice")
  cmd <- sprintf("awk -F'\\t' -v pre=%s 'NR==1 || index($2, pre)==1' %s",
                 shQuote(paste0(chr, "_")), shQuote(file))
  raw <- fread(cmd = cmd)
  check_columns(names(raw), GTEX_V7_COLUMNS, "GTEx v7 slice")
  if (!nrow(raw)) return(raw[, `:=`(pos = integer(), ak = character())])
  parts <- tstrsplit(raw$variant_id, "_", fixed = TRUE)
  if (length(parts) != 5) stop("GTEx v7 slice: variant_id is not chr_pos_ref_alt_build", call. = FALSE)
  raw[, `:=`(pos = as.integer(parts[[2]]), ak = allele_key(parts[[3]], parts[[4]]))]
  raw[pos >= lo & pos <= hi]
}

#' GTEx v7 Brain Cortex (the dataset ezQTL uses). Has slope, slope_se and maf.
#' Expression is treated as unit variance (sdY = 1); GTEx normalises expression before mapping.
read_eqtl_gtex_v7 <- function(file, chr, lo, hi) {
  raw <- read_gtex_v7_rows(file, chr, lo, hi)
  kept <- raw[!is.na(slope_se) & slope_se > 0]
  data <- kept[, .(gene_id = sub("\\..*", "", gene_id), gene_name = NA_character_, pos,
                   rsid = NA_character_, ak, beta = slope, se = slope_se, p = pval_nominal, maf)]
  list(data = data,
       spec = list(match_by = "pos_alleles", n = NULL, sd_y = 1,
                   approximations = "expression treated as unit variance (sdY = 1)"),
       counts = c(eqtl_rows_in_window = nrow(raw), eqtl_rows_kept = nrow(kept)))
}

#' PsychENCODE full cis-eQTL set (QTLtools nominal output: p-value and slope only).
#' No alleles, no SE, no MAF in the file, so: SE is recovered from p and slope, SNPs are matched on
#' position only, and MAF is borrowed from the GTEx v7 file at the same position.
read_eqtl_psychencode <- function(file, gtex_v7_file, chr, lo, hi, n = PSYCHENCODE_N) {
  check_file(file, "PsychENCODE slice")
  cmd <- sprintf("awk -v c=%s -v lo=%s -v hi=%s '$9==c && $10>=lo && $10<=hi' %s",
                 shQuote(paste0("chr", chr)), format_position(lo), format_position(hi), shQuote(file))
  raw <- fread_allow_empty(cmd, header = FALSE)
  if (!nrow(raw)) {
    return(list(data = empty_eqtl(),
                spec = list(match_by = "pos", n = n, sd_y = NULL, approximations = ""),
                counts = c(eqtl_rows_in_window = 0L, eqtl_rows_kept = 0L)))
  }
  if (ncol(raw) != length(PSYCHENCODE_COLUMNS)) {
    stop(sprintf("PsychENCODE slice: expected %d columns, found %d",
                 length(PSYCHENCODE_COLUMNS), ncol(raw)), call. = FALSE)
  }
  setnames(raw, PSYCHENCODE_COLUMNS)
  layout_ok <- all(grepl("^ENSG", raw$gene)) && all(grepl("^chr", raw$variant_chr)) &&
    is.numeric(raw$variant_start) && is.numeric(raw$p) && is.numeric(raw$slope) &&
    all(raw$p >= 0 & raw$p <= 1) && all(raw$variant_start == raw$variant_end)
  if (!layout_ok) stop("PsychENCODE slice: columns do not look like the QTLtools nominal layout", call. = FALSE)

  usable <- raw[p > 0 & p < 1 & slope != 0]
  maf_by_pos <- read_gtex_v7_rows(gtex_v7_file, chr, lo, hi)[, .(maf = maf[1]), by = pos]
  usable <- usable[, .(gene_id = sub("\\..*", "", gene), gene_name = NA_character_,
                       pos = as.integer(variant_start), rsid = NA_character_, ak = NA_character_,
                       beta = slope, se = se_from_p(slope, p), p)]
  data <- merge(usable, maf_by_pos, by = "pos")
  setcolorder(data, names(empty_eqtl()))
  list(data = data,
       spec = list(match_by = "pos", n = n, sd_y = NULL,
                   approximations = paste("SE recovered from p-value and slope;",
                                          "matched on position only (file has no alleles);",
                                          "MAF borrowed from GTEx v7 at the same position")),
       counts = c(eqtl_rows_in_window = nrow(raw), eqtl_rows_p_and_slope_usable = nrow(usable),
                  eqtl_rows_kept = nrow(data)))
}

#' Bryois et al. 2022 single-nucleus eQTL, one cell type, one chromosome
#' (Zenodo 10.5281/zenodo.6104982). The file has p-value and beta; SE is recovered from them.
#' hg19 position, alleles and MAF come from the study's snp_pos table.
read_eqtl_bryois <- function(dir, cell_type, chr, lo, hi, n = BRYOIS_N) {
  file <- file.path(dir, sprintf("%s.%s.gz", cell_type, chr))
  snp_file <- file.path(dir, "snp_pos.txt.gz")
  check_file(file, sprintf("Bryois file for %s chromosome %s", cell_type, chr))
  check_file(snp_file, "Bryois snp_pos table")

  # A whole-chromosome file takes seconds to read and a locus needs it once per window, so the
  # last file read (and the last chromosome of the SNP table) is kept in memory.
  if (!identical(READ_CACHE$bryois_file, file)) {
    raw <- fread(cmd = stream_command(file), header = FALSE)
    if (ncol(raw) != length(BRYOIS_COLUMNS)) {
      stop(sprintf("Bryois file: expected %d columns, found %d", length(BRYOIS_COLUMNS), ncol(raw)), call. = FALSE)
    }
    setnames(raw, BRYOIS_COLUMNS)
    if (!all(grepl("_ENSG[0-9]+$", raw$gene)) || !is.numeric(raw$p) || !is.numeric(raw$beta)) {
      stop("Bryois file: columns do not look like gene, rsID, distance, p-value, beta", call. = FALSE)
    }
    READ_CACHE$bryois_file <- file
    READ_CACHE$bryois_raw <- raw
  }
  raw <- READ_CACHE$bryois_raw

  snp_key <- paste(snp_file, chr)
  if (!identical(READ_CACHE$bryois_snp_key, snp_key)) {
    cmd <- sprintf("%s | awk -F'\\t' -v pre=%s 'NR==1 || index($3, pre)==1'",
                   stream_command(snp_file), shQuote(paste0("chr", chr, ":")))
    chromosome_snps <- fread(cmd = cmd)
    check_columns(names(chromosome_snps), BRYOIS_SNP_COLUMNS, "Bryois snp_pos table")
    chromosome_snps[, pos := as.integer(sub(".*:", "", SNP_id_hg19))]
    READ_CACHE$bryois_snp_key <- snp_key
    READ_CACHE$bryois_snps <- chromosome_snps
  }
  snps <- READ_CACHE$bryois_snps
  snps <- snps[pos >= lo & pos <= hi,
               .(rsid = SNP, pos, ak = allele_key(effect_allele, other_allele), maf = MAF)]

  in_window <- merge(raw, snps, by = "rsid")
  kept <- in_window[p > 0 & p < 1 & beta != 0 & maf > 0 & maf < 1]
  data <- kept[, .(gene_id = sub(".*_", "", gene), gene_name = sub("_ENSG[0-9]+$", "", gene),
                   pos, rsid, ak, beta, se = se_from_p(beta, p), p, maf)]
  list(data = data,
       spec = list(match_by = "rsid_alleles", n = n, sd_y = NULL,
                   approximations = "SE recovered from p-value and beta"),
       counts = c(eqtl_rows_in_window = nrow(in_window), eqtl_rows_kept = nrow(data)))
}

#' Shared by the two readers below: filter and standardise rows in eQTL Catalogue column names.
standardise_catalogue_rows <- function(raw) {
  kept <- raw[!is.na(beta) & !is.na(se) & se > 0 & maf > 0 & maf < 1 & !is.na(rsid) & rsid != "NA"]
  kept[, .(gene_id, gene_name = NA_character_, pos = as.integer(position), rsid,
           ak = allele_key(ref, alt), beta, se, p = pvalue, maf)]
}

#' A region slice of an eQTL Catalogue "all associations" file (GRCh38, no header, 19 columns).
#' The build differs from the GWAS (hg19), so SNPs are matched on rsID plus allele pair and the
#' slice is not filtered by position here. Sample size is read from the allele number column.
read_eqtl_catalogue <- function(file) {
  check_file(file, "eQTL Catalogue slice")
  raw <- fread(file, header = FALSE)
  if (ncol(raw) != length(CATALOGUE_COLUMNS)) {
    stop(sprintf("eQTL Catalogue slice: expected %d columns, found %d",
                 length(CATALOGUE_COLUMNS), ncol(raw)), call. = FALSE)
  }
  setnames(raw, CATALOGUE_COLUMNS)
  if (!all(grepl("^ENSG", raw$gene_id)) || !all(grepl("^chr", raw$variant)) || !is.numeric(raw$an)) {
    stop("eQTL Catalogue slice: columns do not look like the 'all associations' layout", call. = FALSE)
  }
  data <- standardise_catalogue_rows(raw)
  list(data = data,
       spec = list(match_by = "rsid_alleles", n = max(raw$an, na.rm = TRUE) / 2, sd_y = NULL,
                   positions_build = "GRCh38", region_read = range(raw$position), coverage = NULL,
                   approximations = "positions are GRCh38; matched to the hg19 GWAS by rsID and alleles"),
       counts = c(eqtl_rows_in_slice = nrow(raw), eqtl_rows_kept = nrow(data)))
}

#' Regions present in a region file made by fetch_gtex_v8_region(); NULL for a whole-genome file.
#' The list sits next to the data file: <Tissue>.hg38.tsv.gz -> <Tissue>.hg38.regions.tsv.
read_region_coverage <- function(file) {
  manifest <- sub("\\.tsv\\.gz$", ".regions.tsv", file)
  if (!file.exists(manifest)) return(NULL)
  coverage <- fread(manifest, colClasses = list(character = "chromosome"))
  check_columns(names(coverage), c("chromosome", "start", "end"), "region list")
  coverage
}

#' The hg19 windows a region file was downloaded for (written by fetch_gtex_v8_region), or NULL.
#' It sits next to the data file: <Tissue>.hg38.tsv.gz -> <Tissue>.hg38.fetched_windows.tsv.
read_fetched_windows <- function(file) {
  listing <- sub("\\.tsv\\.gz$", ".fetched_windows.tsv", file)
  if (!file.exists(listing)) return(NULL)
  windows <- fread(listing, colClasses = list(character = "chromosome"))
  if (!all(c("chromosome", "hg19_lo", "hg19_hi") %in% names(windows))) {
    stop("fetched-windows list: columns chromosome, hg19_lo and hg19_hi are needed", call. = FALSE)
  }
  windows
}

#' Official GTEx v8, one tissue, every SNP-gene pair (GTEx's own analysis, as imported by the eQTL
#' Catalogue; GRCh38; bgzip-compressed with a tabix index next to it). Read by region with tabix.
#' Brain Cortex is a whole-genome file; other tissues are region files holding only the downloaded
#' loci, with a region list next to them that run_coloc() checks the locus against.
#' The window is given in hg19, so the region read is the window plus `pad_mb` on each side and SNPs
#' are matched to the GWAS by rsID plus allele pair. Expression is treated as unit variance
#' (sdY = 1), the same convention as the GTEx v7 reader; the sample size is recorded but not used.
read_eqtl_gtex_v8 <- function(file, chr, lo, hi, pad_mb = BUILD_PAD_MB_DEFAULT) {
  check_file(file, "GTEx v8 file")
  check_file(paste0(file, ".tbi"), "GTEx v8 tabix index")
  con <- gzfile(file)
  header <- tryCatch(strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]], finally = close(con))
  check_columns(header, GTEX_V8_COLUMNS, "GTEx v8 file")

  pad <- round(pad_mb * 1e6)
  region <- c(max(1, lo - pad), hi + pad)
  cmd <- sprintf("%s %s %s", shQuote(find_tabix()), shQuote(file),
                 shQuote(sprintf("%s:%s-%s", chr, format_position(region[1]), format_position(region[2]))))
  text_columns <- c("variant", "molecular_trait_object_id", "molecular_trait_id", "gene_id", "chromosome",
                    "ref", "alt", "type", "rsid")
  raw <- fread_allow_empty(cmd, header = FALSE, col.names = GTEX_V8_COLUMNS,
                           colClasses = list(character = match(text_columns, GTEX_V8_COLUMNS)))
  coverage <- read_region_coverage(file)
  fetched <- read_fetched_windows(file)
  spec <- list(match_by = "rsid_alleles", n = NULL, sd_y = 1, positions_build = "GRCh38", region_read = region,
               coverage = if (is.null(coverage)) NULL else coverage[chromosome == as.character(chr)],
               # TRUE when this exact hg19 window is known to have been downloaded for this tissue
               window_downloaded = !is.null(fetched) && nrow(fetched[chromosome == as.character(chr) & hg19_lo <= lo & hg19_hi >= hi]) > 0,
               approximations = paste("expression treated as unit variance (sdY = 1);",
                                      "positions are GRCh38; matched to the hg19 GWAS by rsID and alleles;",
                                      sprintf("region read is the hg19 window plus %s Mb each side", pad_mb)))
  if (!nrow(raw)) {
    return(list(data = empty_eqtl(), spec = spec, counts = c(eqtl_rows_in_region_read = 0L, eqtl_rows_kept = 0L)))
  }
  if (!all(grepl("^ENSG", raw$gene_id)) || !is.numeric(raw$an) || !is.numeric(raw$beta) || !is.numeric(raw$se)) {
    stop("GTEx v8 file: columns do not look like the eQTL Catalogue import layout", call. = FALSE)
  }
  data <- standardise_catalogue_rows(raw)
  spec$n <- max(raw$an, na.rm = TRUE) / 2
  list(data = data, spec = spec,
       counts = c(eqtl_rows_in_region_read = nrow(raw), eqtl_rows_kept = nrow(data)))
}

# ---- Matching -----------------------------------------------------------------------------------

#' Add the column `key` used to join GWAS and eQTL rows, and drop rows that cannot be keyed.
add_match_key <- function(dt, match_by) {
  match_by <- match.arg(match_by, MATCH_MODES)
  dt <- copy(dt)
  if (match_by == "pos_alleles") {
    dt[, key := paste(pos, ak)]
  } else if (match_by == "pos") {
    dt[, key := as.character(pos)]
  } else {
    dt <- dt[!is.na(rsid) & rsid != "." & rsid != ""]
    dt[, key := paste(rsid, ak)]
  }
  dt[]
}

#' Drop every row whose key occurs more than once (both copies), since they cannot be told apart.
drop_ambiguous_keys <- function(dt) {
  dt[!key %in% key[duplicated(key)]]
}

#' For eQTL data in another build: keep only rows whose SNP is in the GWAS window, and make sure the
#' data that were available reach well past those SNPs on both sides (otherwise the window was cut).
#' `coverage` is the region list of a region file (NULL when the whole genome is in the file).
#' @return list(eqtl = restricted table, info = named list for the settings file)
restrict_to_gwas_window <- function(eqtl, gwas, region_read, coverage = NULL, min_margin = MIN_EDGE_MARGIN_BP,
                                    window_downloaded = FALSE) {
  in_window <- eqtl$key %in% gwas$key
  matched_pos <- eqtl$pos[in_window]
  how_to_fetch <- if (is.null(coverage)) "" else
    "; this eQTL file holds only downloaded regions, so fetch this locus with code/fetch_gtex_v8_regions.R"
  if (!length(matched_pos) && isTRUE(window_downloaded)) {
    # The window was downloaded for this dataset and still nothing matches: the dataset tests no gene
    # within reach of the window (a gene-poor stretch for this tissue). That is an empty result, not an error.
    return(list(eqtl = eqtl[0],
                info = list(eqtl_region_read = sprintf("%s-%s", format_position(region_read[1]), format_position(region_read[2])),
                            eqtl_region_available = "window downloaded for this dataset",
                            eqtl_span_of_matched_snps = "none: this dataset tests no gene with a SNP in the GWAS window",
                            eqtl_rows_for_snps_outside_gwas_window = nrow(eqtl),
                            genes_with_no_snp_in_gwas_window = uniqueN(eqtl$gene_id))))
  }
  if (!length(matched_pos)) {
    stop("no SNP in the GWAS window was found by rsID in the eQTL region read; ",
         "the region may not reach this locus in the eQTL file's build", how_to_fetch, call. = FALSE)
  }
  available <- region_read
  if (!is.null(coverage)) {
    holding <- coverage[start <= min(matched_pos) & end >= max(matched_pos)]
    if (!nrow(holding)) {
      stop(sprintf("the matched SNPs (%s-%s) are not inside one downloaded region of this eQTL file%s",
                   format_position(min(matched_pos)), format_position(max(matched_pos)), how_to_fetch), call. = FALSE)
    }
    available <- c(max(region_read[1], holding$start[1]), min(region_read[2], holding$end[1]))
  }
  low_margin <- if (available[1] <= 1) Inf else min(matched_pos) - available[1]
  high_margin <- available[2] - max(matched_pos)
  if (min(low_margin, high_margin) < min_margin) {
    stop(sprintf(paste("matched SNPs (%s-%s) come within %s bp of the edge of the eQTL data available (%s-%s);",
                       "the window is probably cut off, so read a wider region%s"),
                 format_position(min(matched_pos)), format_position(max(matched_pos)),
                 format_position(min_margin), format_position(available[1]), format_position(available[2]),
                 how_to_fetch),
         call. = FALSE)
  }
  list(eqtl = eqtl[in_window],
       info = list(eqtl_region_read = sprintf("%s-%s", format_position(region_read[1]), format_position(region_read[2])),
                   eqtl_region_available = sprintf("%s-%s", format_position(available[1]), format_position(available[2])),
                   eqtl_span_of_matched_snps = sprintf("%s-%s", format_position(min(matched_pos)),
                                                       format_position(max(matched_pos))),
                   eqtl_rows_for_snps_outside_gwas_window = sum(!in_window),
                   genes_with_no_snp_in_gwas_window = uniqueN(eqtl$gene_id) - uniqueN(eqtl$gene_id[in_window])))
}

# ---- Region files for GTEx v8 tissues ------------------------------------------------------------

#' Merge overlapping or touching intervals, per chromosome.
merge_intervals <- function(intervals) {
  if (!nrow(intervals)) return(intervals[, .(chromosome, start, end)])
  x <- copy(intervals)[order(chromosome, start, end)]
  x[, group := cumsum(c(TRUE, start[-1] > cummax(end)[-.N] + 1)), by = chromosome]
  x[, .(start = min(start), end = max(end)), by = .(chromosome, group)][, .(chromosome, start, end)]
}

#' Add rows (a headerless file in GTEX_V8_COLUMNS layout) to a region file, keeping it sorted,
#' free of duplicate rows, bgzip-compressed and indexed, and record the region in its region list.
add_rows_to_region_file <- function(file, new_rows_file, chromosome, start, end) {
  tabix <- find_tabix()
  bgzip <- file.path(dirname(tabix), "bgzip")
  if (!file.exists(bgzip)) stop("bgzip not found next to tabix; install htslib", call. = FALSE)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  body <- tempfile("region_body_")
  building <- paste0(file, ".building")
  on.exit(unlink(c(body, building)), add = TRUE)

  existing <- if (file.exists(file)) sprintf("gzip -dc %s | tail -n +2; ", shQuote(file)) else ""
  sort_cmd <- sprintf("{ %scat %s; } | LC_ALL=C sort -u | LC_ALL=C sort -t '\t' -k13,13 -k14,14n > %s",
                      existing, shQuote(new_rows_file), shQuote(body))
  if (system(sort_cmd) != 0) stop("could not sort the rows for ", file, call. = FALSE)
  pack_cmd <- sprintf("{ printf '%%s\\n' %s; cat %s; } | %s -c > %s",
                      shQuote(paste(GTEX_V8_COLUMNS, collapse = "\t")), shQuote(body), shQuote(bgzip), shQuote(building))
  if (system(pack_cmd) != 0) stop("could not compress ", file, call. = FALSE)
  if (!file.rename(building, file)) stop("could not write ", file, call. = FALSE)
  if (system2(tabix, c("-f", "-s", "13", "-b", "14", "-e", "14", "-S", "1", shQuote(file))) != 0) {
    stop("could not index ", file, call. = FALSE)
  }
  old <- read_region_coverage(file)
  if (is.null(old)) old <- data.table(chromosome = character(), start = numeric(), end = numeric())
  coverage <- merge_intervals(rbind(old, data.table(chromosome = as.character(chromosome), start = start, end = end)))
  fwrite(coverage, sub("\\.tsv\\.gz$", ".regions.tsv", file), sep = "\t")
  invisible(coverage)
}

#' Download one GRCh38 region of an official GTEx v8 tissue file from the EBI server (through its
#' tabix index, so only that region travels) and add it to the tissue's region file.
#' Only genome coordinates are sent; no study data leaves this machine.
#' @return list(file, rows, status) where status is "downloaded" or "already present"
fetch_gtex_v8_region <- function(tissue, chr, from, to, dir, hg19_window = NULL) {
  file <- file.path(dir, paste0(tissue, ".hg38.tsv.gz"))
  note_window <- function() {
    if (is.null(hg19_window)) return(invisible())
    listing <- sub("\\.tsv\\.gz$", ".fetched_windows.tsv", file)
    row <- data.table(chromosome = as.character(chr), hg19_lo = hg19_window[1], hg19_hi = hg19_window[2],
                      grch38_from = from, grch38_to = to, noted = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    fwrite(row, listing, sep = "\t", append = file.exists(listing))
  }
  coverage <- read_region_coverage(file)
  if (!is.null(coverage) && nrow(coverage[chromosome == as.character(chr) & start <= from & end >= to])) {
    note_window()
    return(list(file = file, rows = NA_integer_, status = "already present"))
  }
  url <- sprintf(GTEX_V8_TISSUE_URL, tissue)
  index_dir <- file.path(dir, "idx")   # tabix drops the remote index into the working directory
  dir.create(index_dir, showWarnings = FALSE, recursive = TRUE)
  new_rows <- tempfile("region_new_")
  on.exit(unlink(new_rows), add = TRUE)
  region <- sprintf("%s:%s-%s", chr, format_position(from), format_position(to))
  cmd <- sprintf("cd %s && %s %s %s > %s", shQuote(index_dir), shQuote(find_tabix()), shQuote(url),
                 shQuote(region), shQuote(new_rows))
  if (system(cmd) != 0) stop(sprintf("download of %s %s failed", tissue, region), call. = FALSE)

  head_rows <- fread_allow_empty(sprintf("head -200 %s", shQuote(new_rows)), header = FALSE)
  if (!nrow(head_rows)) stop(sprintf("the server returned no rows for %s %s", tissue, region), call. = FALSE)
  layout_ok <- ncol(head_rows) == length(GTEX_V8_COLUMNS) && all(grepl("^ENSG", head_rows[[7]])) &&
    all(as.character(head_rows[[13]]) == as.character(chr)) && all(head_rows[[14]] >= from & head_rows[[14]] <= to)
  if (!layout_ok) stop(sprintf("rows returned for %s %s are not in the expected layout", tissue, region), call. = FALSE)
  rows <- as.integer(system(sprintf("wc -l < %s", shQuote(new_rows)), intern = TRUE))

  add_rows_to_region_file(file, new_rows, chr, from, to)
  note_window()
  log_file <- file.path(dir, "FETCH_LOG.tsv")
  fwrite(data.table(fetched = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), tissue = tissue, region_grch38 = region,
                    rows = rows, url = url), log_file, sep = "\t", append = file.exists(log_file))
  list(file = file, rows = rows, status = "downloaded")
}

#' Where an hg19 window sits in GRCh38, found by matching the window's SNPs by rsID in the local
#' whole-genome GTEx v8 Brain Cortex file. @return c(lowest, highest) matched GRCh38 position.
grch38_span_of_window <- function(paths, chr, lo, hi) {
  gwas <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file_for(paths, "pan"), chr, lo, hi)$data,
                                            "rsid_alleles"))
  cortex <- read_eqtl_gtex_v8(paths$gtex_v8, chr, lo, hi)
  matched <- restrict_to_gwas_window(add_match_key(cortex$data, "rsid_alleles"), gwas,
                                     cortex$spec$region_read, cortex$spec$coverage)
  range(matched$eqtl$pos)
}

# ---- Duplicated sequence --------------------------------------------------------------------------

#' Read a headerless UCSC table and name its columns, checking the column count.
read_ucsc_table <- function(file, columns, what) {
  check_file(file, what)
  raw <- fread(cmd = stream_command(file), header = FALSE, sep = "\t")
  if (ncol(raw) != length(columns)) {
    stop(sprintf("%s: expected %d columns, found %d", what, length(columns), ncol(raw)), call. = FALSE)
  }
  setnames(raw, columns)
  raw[grepl("^chr([0-9]+|X|Y|M)$", chrom)]   # main chromosomes only, not alternate haplotypes
}

#' One row per Ensembl gene (hg19): from the first to the last base of any of its transcripts.
#' Coordinates are UCSC's (start counted from 0), the same as the duplication table.
read_gene_spans <- function(file) {
  raw <- read_ucsc_table(file, UCSC_ENSGENE_COLUMNS, "UCSC ensGene table")
  if (!all(grepl("^ENSG", raw$name2))) stop("UCSC ensGene table: gene IDs not found in column 13", call. = FALSE)
  raw[, .(start = min(txStart), end = max(txEnd)), by = .(gene_id = name2, chromosome = sub("^chr", "", chrom))]
}

#' Segmental duplications (hg19): one row per duplication record, with the identity between copies.
read_segmental_duplications <- function(file) {
  raw <- read_ucsc_table(file, UCSC_SEGDUP_COLUMNS, "UCSC genomicSuperDups table")
  raw[, .(chromosome = sub("^chr", "", chrom), start = chromStart, end = chromEnd, identity = fracMatch)]
}

#' Both tables, or NULL (with a message) when either file is missing.
load_duplication_annotation <- function(paths) {
  if (is.null(paths$gene_spans) || is.null(paths$segdups) || !file.exists(paths$gene_spans) || !file.exists(paths$segdups)) {
    message("duplication annotation not found (data/annotation/); the duplication columns will be left out")
    return(NULL)
  }
  list(gene_spans = read_gene_spans(paths$gene_spans), duplications = read_segmental_duplications(paths$segdups))
}

#' Total length covered by a set of intervals (start counted from 0, end exclusive), overlaps counted once.
union_length <- function(start, end) {
  if (!length(start)) return(0)
  o <- order(start, end)
  start <- start[o]
  end <- end[o]
  reach <- cummax(end)
  block <- cumsum(c(TRUE, start[-1] > reach[-length(reach)]))
  sum(tapply(end, block, max) - tapply(start, block, min))
}

#' For each gene: its hg19 span, the fraction of that span inside a segmental duplication, and the
#' highest identity of any overlapping duplication. Genes missing from the gene table get NA
#' (unknown); genes with no overlap get 0.
duplication_by_gene <- function(gene_ids, annotation) {
  gene_ids <- unique(gene_ids)
  spans <- annotation$gene_spans[gene_id %in% gene_ids]
  spans <- spans[, .SD[1], by = gene_id]   # a gene listed on two chromosomes (X and Y) is taken once
  dups <- annotation$duplications[chromosome %in% spans$chromosome]
  setkey(dups, chromosome, start, end)
  hits <- if (nrow(spans) && nrow(dups)) {
    foverlaps(spans, dups, by.x = c("chromosome", "start", "end"), nomatch = NULL)
  } else {
    data.table(gene_id = character(), start = numeric(), end = numeric(), i.start = numeric(), i.end = numeric(),
               identity = numeric())
  }
  per_gene <- if (nrow(hits)) {
    hits[, .(duplicated_bp = union_length(pmax(start, i.start), pmin(end, i.end)),
             max_duplication_identity = max(identity)), by = gene_id]
  } else {
    data.table(gene_id = character(), duplicated_bp = numeric(), max_duplication_identity = numeric())
  }
  out <- merge(spans, per_gene, by = "gene_id", all.x = TRUE)
  out[is.na(duplicated_bp), duplicated_bp := 0]
  out <- out[, .(gene_id, gene_span_hg19 = sprintf("%s:%s-%s", chromosome, format_position(start + 1), format_position(end)),
                 fraction_in_segmental_duplication = round(duplicated_bp / (end - start), 3), max_duplication_identity)]
  merge(data.table(gene_id = gene_ids), out, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# ---- Gene names ---------------------------------------------------------------------------------

#' Ensembl ID to symbol, from FUMA's MAGMA gene file, with a count of how many IDs share each
#' symbol. Two different genes can carry one symbol.
load_symbol_map <- function(magma_file) {
  if (is.null(magma_file) || !file.exists(magma_file)) {
    return(data.table(gene_id = character(), symbol = character(), ids_sharing_symbol = integer()))
  }
  magma <- fread(magma_file)
  if (!all(c("GENE", "SYMBOL") %in% names(magma))) {
    stop("MAGMA gene file: GENE and SYMBOL columns not found", call. = FALSE)
  }
  map <- unique(magma[, .(gene_id = GENE, symbol = SYMBOL)])
  map[, ids_sharing_symbol := .N, by = symbol]
  map[]
}

# ---- coloc --------------------------------------------------------------------------------------

#' Run coloc.abf for one gene on already matched rows, under the default and the strict p12 prior.
#' @return list(summary, summary_strict = named numeric or NULL, top = one-row data.frame or NULL,
#'              error = message or NA, warnings = character vector from the default-prior run)
coloc_one_gene <- function(matched, study, spec, p12 = COLOC_P12_DEFAULT, p12_strict = COLOC_P12_STRICT) {
  gwas <- list(beta = matched$gwas_beta, varbeta = matched$gwas_se^2, snp = matched$key, type = "cc",
               s = study$n_cases / (study$n_cases + study$n_controls),
               N = study$n_cases + study$n_controls)
  eqtl <- list(beta = matched$beta, varbeta = matched$se^2, snp = matched$key, type = "quant")
  if (!is.null(spec$sd_y)) {
    eqtl$sdY <- spec$sd_y
  } else {
    eqtl$N <- spec$n
    eqtl$MAF <- pmin(matched$maf, 1 - matched$maf)
  }
  warnings_seen <- character()
  error <- NA_character_
  run_abf <- function(prior, keep_warnings) {
    fit <- NULL
    invisible(capture.output(
      fit <- tryCatch(
        withCallingHandlers(
          suppressMessages(coloc.abf(gwas, eqtl, p12 = prior)),
          warning = function(w) {
            if (keep_warnings) warnings_seen <<- c(warnings_seen, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
        error = function(e) {
          error <<- conditionMessage(e)
          NULL
        })
    ))
    fit
  }
  fit <- run_abf(p12, keep_warnings = TRUE)
  strict <- if (is.null(fit)) NULL else run_abf(p12_strict, keep_warnings = FALSE)
  if (is.null(fit) || is.null(strict)) {
    return(list(summary = NULL, summary_strict = NULL, top = NULL, error = error, warnings = warnings_seen))
  }
  list(summary = fit$summary, summary_strict = strict$summary,
       top = fit$results[which.max(fit$results$SNP.PP.H4), ], error = NA_character_, warnings = warnings_seen)
}

#' The outcome of one gene test in words (see the constants above GWAS_SIGNAL_P).
#' `gwas_signal` and `signal_kept` may be NA when no reference p-value was given; they then do not gate.
outcome_label <- function(h3, h4, expression_signal, signal_kept = NA_real_, gwas_signal = NA) {
  if (isFALSE(gwas_signal)) return("not testable: no migraine signal in this GWAS file")
  if (!isTRUE(expression_signal)) return("not testable: no expression signal")
  if (!is.na(signal_kept) && signal_kept < SIGNAL_KEPT_MIN) return("not testable: migraine signal not in the tested SNPs")
  if (h4 >= SUPPORT_THRESHOLD) "shared" else if (h3 >= TWO_SIGNALS_MIN) "two signals" else if (h4 >= LEANING_MIN) "leaning" else "undecided"
}

#' Run coloc for every gene in an eQTL table against one GWAS window.
#'
#' @param gwas,eqtl tables that already carry `key` (see add_match_key) with ambiguous GWAS keys dropped.
#' @param gwas_reference_p the migraine signal at the locus (see migraine_reference_p); used to decide
#'   whether each test is testable. NA skips the two checks that need it.
#' @return list(results, skipped, warnings). Every gene in `eqtl` appears in exactly one of
#'   `results` or `skipped`; nothing is dropped silently.
coloc_locus <- function(gwas, eqtl, spec, study = STUDY, min_snps = MIN_SNPS_DEFAULT,
                        label = "", dataset_label = "", lead_pos = NA, symbol_map = load_symbol_map(NULL),
                        gwas_reference_p = NA_real_, gwas_version = "pan") {
  if (is.null(study)) need_study_config()
  lead_key <- gwas[pos == lead_pos, key][1]
  gwas_signal <- if (is.na(gwas_reference_p)) NA else gwas_reference_p < GWAS_SIGNAL_P
  results <- list()
  skipped <- list()
  all_warnings <- character()

  for (gid in unique(eqtl$gene_id)) {
    gene_rows <- drop_ambiguous_keys(eqtl[gene_id == gid])
    source_name <- gene_rows$gene_name[1]
    matched <- merge(gwas, gene_rows[, .(key, beta, se, p, maf)], by = "key")

    if (nrow(matched) < min_snps) {
      skipped[[gid]] <- data.table(gene_id = gid, n_snps = nrow(matched),
                                   reason = sprintf("fewer than %d matched SNPs", min_snps))
      next
    }
    fit <- coloc_one_gene(matched, study, spec)
    all_warnings <- c(all_warnings, fit$warnings)
    if (is.null(fit$summary)) {
      skipped[[gid]] <- data.table(gene_id = gid, n_snps = nrow(matched),
                                   reason = paste("coloc.abf error:", fit$error))
      next
    }

    expression_signal <- min(matched$p) < EXPRESSION_SIGNAL_P
    best_tested_p <- min(matched$gwas_p)
    signal_kept <- if (is.na(gwas_reference_p) || gwas_reference_p >= 1) {
      NA_real_
    } else if (best_tested_p <= gwas_reference_p) {
      1
    } else {
      log10(best_tested_p) / log10(gwas_reference_p)
    }

    magma_row <- symbol_map[gene_id == gid]
    magma_symbol <- if (nrow(magma_row)) magma_row$symbol[1] else NA_character_
    has_source_name <- !is.na(source_name) && nzchar(source_name)
    results[[gid]] <- data.table(
      locus = label, eqtl_dataset = dataset_label, gwas_version = gwas_version, gene_id = gid,
      symbol = if (has_source_name) source_name else if (!is.na(magma_symbol)) magma_symbol else "",
      symbol_from = if (has_source_name) "eqtl_source" else if (!is.na(magma_symbol)) "fuma_magma" else "",
      fuma_symbol = if (!is.na(magma_symbol)) magma_symbol else "",
      fuma_symbol_shared_by_n_gene_ids = if (nrow(magma_row)) magma_row$ids_sharing_symbol[1] else NA_integer_,
      n_snps = nrow(matched), lead_snp_included = !is.na(lead_key) && lead_key %in% matched$key,
      min_p_gwas = min(matched$gwas_p), min_p_eqtl = min(matched$p),
      PP.H0 = fit$summary[["PP.H0.abf"]], PP.H1_gwas_only = fit$summary[["PP.H1.abf"]],
      PP.H2_eqtl_only = fit$summary[["PP.H2.abf"]], PP.H3_two_signals = fit$summary[["PP.H3.abf"]],
      PP.H4_shared = fit$summary[["PP.H4.abf"]],
      PP.H3_two_signals_strict_prior = fit$summary_strict[["PP.H3.abf"]],
      PP.H4_shared_strict_prior = fit$summary_strict[["PP.H4.abf"]],
      expression_signal = expression_signal, migraine_signal_kept = round(signal_kept, 3),
      outcome = outcome_label(fit$summary[["PP.H3.abf"]], fit$summary[["PP.H4.abf"]],
                              expression_signal, signal_kept, gwas_signal),
      outcome_strict_prior = outcome_label(fit$summary_strict[["PP.H3.abf"]], fit$summary_strict[["PP.H4.abf"]],
                                           expression_signal, signal_kept, gwas_signal),
      top_snp = matched$rsid[match(fit$top$snp, matched$key)], top_snp_PP.H4 = fit$top$SNP.PP.H4,
      n_coloc_warnings = length(fit$warnings))
  }

  results <- rbindlist(results)
  if (nrow(results)) setorder(results, -PP.H4_shared)
  skipped <- rbindlist(skipped)
  if (!nrow(skipped)) skipped <- data.table(gene_id = character(), n_snps = integer(), reason = character())
  list(results = results, skipped = skipped, warnings = all_warnings)
}

#' Reduce coloc's warnings to their first line with numbers masked, so they can be counted by kind.
normalise_warnings <- function(messages) {
  first_line <- sub("\n.*", "", messages)
  gsub("[0-9][0-9.]*(e-?[0-9]+)?", "<n>", first_line)
}

# ---- One whole run ------------------------------------------------------------------------------

#' Default input locations under the project folder.
project_paths <- function(root) {
  need_study_config()
  fuma_job <- file.path(root, "results", "fuma-downloads", STUDY_FILES$fuma_job_dir)
  list(
    gwas = c(pan = file.path(root, "data", "fuma-input", GWAS_VERSIONS$pan$file),
             eur = file.path(root, "data", "fuma-input", GWAS_VERSIONS$eur$file)),
    catalogue_dir = file.path(root, "data", "eqtl-catalogue"),
    gtex_v7 = file.path(root, "data", "eqtl-full", STUDY_FILES$gtex_v7_slice),
    gtex_v8 = file.path(root, "data", "eqtl-full", "GTEx_v8_Brain_Cortex.allpairs.hg38.tsv.gz"),
    gtex_v8_regions = file.path(root, "data", "eqtl-full", "gtex_v8_regions"),
    fuma_loci = file.path(fuma_job, "GenomicRiskLoci.txt"),
    ancestry = file.path(root, "results", STUDY_FILES$ancestry_table),
    gene_spans = file.path(root, "data", "annotation", "ensGene.hg19.txt.gz"),
    segdups = file.path(root, "data", "annotation", "genomicSuperDups.hg19.txt.gz"),
    psychencode = file.path(root, "data", "eqtl-full", STUDY_FILES$psychencode_slice),
    bryois_dir = file.path(root, "data", "eqtl-full", "bryois2022"),
    magma = file.path(fuma_job, "magma.genes.out"),
    results = file.path(root, "results", "coloc"))
}

#' Split a source string into its kind and its argument.
#' Accepted: "gtex_v7", "gtex_v8" (Brain Cortex), "gtex_v8:<Tissue>" (a region file, for example
#' gtex_v8:Artery_Tibial), "psychencode", "bryois:<CellType>", "catalogue:<path to slice>".
parse_source <- function(source) {
  if (!is.character(source) || length(source) != 1 || is.na(source)) stop("source must be one string", call. = FALSE)
  kind <- sub(":.*", "", source)
  arg <- if (grepl(":", source, fixed = TRUE)) sub("^[^:]*:", "", source) else ""
  if (!kind %in% c("gtex_v7", "gtex_v8", "psychencode", "bryois", "catalogue")) {
    stop(sprintf("unknown source '%s' (use gtex_v7, gtex_v8, gtex_v8:<Tissue>, psychencode, bryois:<CellType> or catalogue:<file>)",
                 source),
         call. = FALSE)
  }
  if (kind %in% c("bryois", "catalogue") && !nzchar(arg)) {
    stop(sprintf("source '%s' needs an argument after the colon", kind), call. = FALSE)
  }
  list(kind = kind, arg = arg)
}

#' Name used in output file names and the eqtl_dataset column when none is given.
default_dataset_label <- function(source) {
  parsed <- parse_source(source)
  switch(parsed$kind,
         gtex_v7 = "GTEx_v7_brain_cortex",
         gtex_v8 = if (nzchar(parsed$arg)) paste0("GTEx_v8_", parsed$arg) else "GTEx_v8_brain_cortex",
         psychencode = "PsychENCODE",
         bryois = paste0("Bryois_", gsub("[^A-Za-z0-9.]+|\\.{2,}", "_", parsed$arg)),
         catalogue = stop("a catalogue source needs dataset_label", call. = FALSE))
}

#' Read the eQTL data for a source.
read_eqtl <- function(source, paths, chr, lo, hi) {
  parsed <- parse_source(source)
  switch(parsed$kind,
         gtex_v7 = read_eqtl_gtex_v7(paths$gtex_v7, chr, lo, hi),
         gtex_v8 = read_eqtl_gtex_v8(
           if (nzchar(parsed$arg)) file.path(paths$gtex_v8_regions, paste0(parsed$arg, ".hg38.tsv.gz")) else paths$gtex_v8,
           chr, lo, hi),
         psychencode = read_eqtl_psychencode(paths$psychencode, paths$gtex_v7, chr, lo, hi),
         bryois = read_eqtl_bryois(paths$bryois_dir, parsed$arg, chr, lo, hi),
         catalogue = read_eqtl_catalogue(parsed$arg))
}

#' What is already known about a locus that decides how far its coloc result can be trusted:
#' how many independent signals FUMA found (coloc.abf assumes one) and how strong the lead SNP is
#' in the European-only GWAS (the eQTL data are mostly European). Looked up by chromosome and lead
#' position; entries say "not found" when a file or the locus is missing.
locus_context <- function(paths, chr, lead_pos, window_kb = NULL) {
  out <- list(fuma_lead_snps_in_locus = "not found", fuma_independent_significant_snps = "not found",
              one_signal_assumption = "not checked", fuma_locus_reach_from_lead_kb = "not found",
              window_covers_fuma_locus = "not checked", gwas_lead_p_european_only = "not found",
              european_signal = "not checked")
  locus_chr <- as.character(chr)
  locus_pos <- round(lead_pos)
  if (!is.null(paths$fuma_loci) && file.exists(paths$fuma_loci)) {
    loci <- fread(paths$fuma_loci)
    hit <- loci[as.character(chr) == locus_chr & pos == locus_pos]
    if (nrow(hit)) {
      out$fuma_lead_snps_in_locus <- hit$nLeadSNPs[1]
      out$fuma_independent_significant_snps <- hit$nIndSigSNPs[1]
      if (all(c("start", "end") %in% names(hit))) {
        reach_kb <- max(locus_pos - hit$start[1], hit$end[1] - locus_pos) / 1000
        out$fuma_locus_reach_from_lead_kb <- round(reach_kb)
        if (!is.null(window_kb)) out$window_covers_fuma_locus <- reach_kb <= window_kb
      }
      out$one_signal_assumption <- if (hit$nLeadSNPs[1] > 1) {
        sprintf("doubtful: %d lead SNPs", hit$nLeadSNPs[1])
      } else if (hit$nIndSigSNPs[1] > 1) {
        sprintf("check: one lead SNP but %d independent significant SNPs", hit$nIndSigSNPs[1])
      } else {
        "ok: one independent significant SNP"
      }
    }
  }
  if (!is.null(paths$ancestry) && file.exists(paths$ancestry)) {
    ancestry <- fread(paths$ancestry)
    hit <- ancestry[as.character(chr) == locus_chr & lead_pos_hg19 == locus_pos]
    if (nrow(hit)) {
      p_eur <- hit$EUR_lead_p[1]
      out$gwas_lead_p_european_only <- if (is.na(p_eur)) "lead SNP absent from the European-only file" else signif(p_eur, 3)
      out$european_signal <- if (is.na(p_eur)) {
        "not testable with European eQTL data"
      } else if (p_eur < 5e-8) {
        "genome-wide significant in Europeans"
      } else if (p_eur < 1e-5) {
        "suggestive in Europeans (5e-8 to 1e-5)"
      } else {
        "weak in Europeans (above 1e-5)"
      }
    }
  }
  out
}

#' Default result file for a run.
result_file <- function(results_dir, label, dataset_label, window_kb, gwas_setting = "pan") {
  file.path(results_dir, sprintf("%s__%s__%skb%s.csv", label, dataset_label,
                                 format(window_kb, scientific = FALSE, trim = TRUE), gwas_suffix(gwas_setting)))
}

#' Run one locus against one eQTL source and write the result, the skipped genes and the settings.
#'
#' Writes <out_file>, and next to it run_logs/<name>.skipped.csv and run_logs/<name>.settings.tsv.
#' @param gwas_version "pan" or "eur" (see GWAS_VERSIONS); the case and control counts follow it
#'   unless `study` is given.
#' @param sample_size_filter TRUE keeps only SNPs whose sample size is within SAMPLE_SIZE_TOLERANCE
#'   of the lead SNP's.
#' @param gwas_read the result of read_gwas_window() for this window, to avoid reading the GWAS
#'   file again when many datasets are run for one locus.
#' @return (invisibly) list(results, skipped, settings, out_file)
run_coloc <- function(label, chr, lead_pos, window_kb, source, out_file = NULL, dataset_label = NULL,
                      paths, study = NULL, min_snps = MIN_SNPS_DEFAULT, gwas_version = "pan",
                      sample_size_filter = FALSE, gwas_read = NULL) {
  chr <- as.character(chr)
  window <- locus_window(lead_pos, window_kb)
  gwas_file <- gwas_file_for(paths, gwas_version)
  gwas_setting <- gwas_label(gwas_version, sample_size_filter)
  if (is.null(study)) study <- GWAS_VERSIONS[[gwas_version]][c("n_cases", "n_controls")]
  if (is.null(dataset_label) || !nzchar(dataset_label)) dataset_label <- default_dataset_label(source)
  if (is.null(out_file) || !nzchar(out_file)) {
    out_file <- result_file(paths$results, label, dataset_label, window_kb, gwas_setting)
  }

  if (is.null(gwas_read)) gwas_read <- read_gwas_window(gwas_file, chr, window$lo, window$hi)
  gwas_data <- gwas_read$data
  size_filter_info <- list(sample_size_filter = "off")
  if (isTRUE(sample_size_filter)) {
    filtered <- filter_sample_size(gwas_data, lead_pos)
    gwas_data <- filtered$data
    size_filter_info <- list(sample_size_filter = sprintf("on: within %s%% of the lead SNP's", 100 * SAMPLE_SIZE_TOLERANCE),
                             sample_size_reference = filtered$reference_n,
                             gwas_rows_dropped_by_sample_size = filtered$dropped)
  }
  reference_p <- migraine_reference_p(gwas_data, lead_pos)

  eqtl_read <- read_eqtl(source, paths, chr, window$lo, window$hi)
  spec <- eqtl_read$spec
  gwas_keyed <- add_match_key(gwas_data, spec$match_by)
  gwas <- drop_ambiguous_keys(gwas_keyed)
  eqtl <- add_match_key(eqtl_read$data, spec$match_by)
  other_build_info <- list()
  if (identical(spec$positions_build, "GRCh38")) {
    restricted <- restrict_to_gwas_window(eqtl, gwas, spec$region_read, spec$coverage,
                                          window_downloaded = isTRUE(spec$window_downloaded))
    eqtl <- restricted$eqtl
    other_build_info <- restricted$info
  }
  lead_rsid <- gwas_data[pos == round(lead_pos), rsid][1]
  lead_key <- gwas[pos == round(lead_pos), key][1]

  fit <- coloc_locus(gwas, eqtl, spec, study = study, min_snps = min_snps, label = label,
                     dataset_label = dataset_label, lead_pos = round(lead_pos),
                     symbol_map = load_symbol_map(paths$magma), gwas_reference_p = reference_p,
                     gwas_version = gwas_setting)

  count_outcome <- function(pattern) if (nrow(fit$results)) sum(grepl(pattern, fit$results$outcome)) else 0L
  warning_counts <- table(normalise_warnings(fit$warnings))
  settings <- c(
    list(run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), locus = label, chr = chr,
         lead_pos_hg19 = format_position(lead_pos), lead_rsid = lead_rsid, window_kb = window_kb,
         window_role = if (window_kb == PRIMARY_WINDOW_KB) "primary" else if (window_kb == COMPARISON_WINDOW_KB) {
           "comparison (ezQTL's window)"
         } else "other",
         window_lo = format_position(window$lo), window_hi = format_position(window$hi),
         source = source, eqtl_dataset = dataset_label, gwas_version = gwas_setting, gwas_file = gwas_file,
         method = "coloc.abf, one gene at a time; p1 = p2 = 1e-4 (coloc defaults)",
         prior_p12_default = COLOC_P12_DEFAULT, prior_p12_strict = COLOC_P12_STRICT,
         n_cases = study$n_cases, n_controls = study$n_controls,
         case_fraction = signif(study$n_cases / (study$n_cases + study$n_controls), 6),
         eqtl_n = if (is.null(spec$n)) "not in file" else spec$n,
         eqtl_sd_y = if (is.null(spec$sd_y)) "estimated by coloc from N and MAF" else spec$sd_y,
         matched_by = spec$match_by, approximations = spec$approximations, min_snps = min_snps),
    locus_context(paths, chr, lead_pos, window_kb),
    as.list(gwas_read$counts), size_filter_info,
    list(gwas_rows_keyed = nrow(gwas_keyed), gwas_rows_after_dropping_ambiguous_keys = nrow(gwas),
         lead_snp_usable = !is.na(lead_rsid),
         migraine_signal_p = if (is.na(reference_p)) "no SNP near the lead position" else signif(reference_p, 4),
         migraine_signal_in_this_gwas_file = if (is.na(reference_p)) FALSE else reference_p < GWAS_SIGNAL_P),
    as.list(eqtl_read$counts), other_build_info,
    list(lead_snp_in_eqtl_data = !is.na(lead_key) && lead_key %in% eqtl$key,
         gwas_snps_with_an_eqtl_match = sum(gwas$key %in% eqtl$key),
         genes_in_eqtl_data = uniqueN(eqtl$gene_id), genes_tested = nrow(fit$results),
         genes_skipped = nrow(fit$skipped),
         tests_testable = count_outcome("^(shared|leaning|two signals|undecided)$"),
         tests_shared = count_outcome("^shared$"), tests_leaning = count_outcome("^leaning$"),
         tests_two_signals = count_outcome("^two signals$"), tests_undecided = count_outcome("^undecided$"),
         tests_not_testable = count_outcome("^not testable"),
         coloc_warnings = if (length(warning_counts)) {
           paste(sprintf("%s (x%d)", names(warning_counts), as.integer(warning_counts)), collapse = " | ")
         } else "none",
         r_version = R.version.string, coloc_version = as.character(packageVersion("coloc")),
         data.table_version = as.character(packageVersion("data.table"))))

  write_run(out_file, fit$results, fit$skipped, settings)
  message(sprintf("%s vs %s (%s): %d GWAS SNPs usable, lead %s, %d genes tested (%d testable, %d shared), %d skipped -> %s",
                  label, dataset_label, gwas_setting, nrow(gwas), lead_rsid, nrow(fit$results),
                  settings$tests_testable, settings$tests_shared, nrow(fit$skipped), out_file))
  invisible(list(results = fit$results, skipped = fit$skipped, settings = settings, out_file = out_file))
}

#' Write a result table plus its two side files.
write_run <- function(out_file, results, skipped, settings) {
  log_dir <- file.path(dirname(out_file), "run_logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  stem <- sub("\\.csv$", "", basename(out_file))
  if (ncol(results)) fwrite(results, out_file) else file.create(out_file)   # no gene tested: an empty file
  fwrite(skipped, file.path(log_dir, paste0(stem, ".skipped.csv")))
  fwrite(data.table(setting = names(settings), value = vapply(settings, function(x) as.character(x)[1], "")),
         file.path(log_dir, paste0(stem, ".settings.tsv")), sep = "\t")
  invisible(out_file)
}

# ---- Summary across datasets --------------------------------------------------------------------

#' Which independent resource a dataset belongs to: datasets drawn from the same people count once.
resource_of <- function(dataset) {
  fifelse(startsWith(dataset, "GTEx"), "GTEx", fifelse(startsWith(dataset, "Bryois"), "Bryois2022", dataset))
}

#' The reporting rule in words, from the number of independent resources at the threshold.
support_label <- function(n_resources) {
  if (n_resources >= 2) "replicated" else if (n_resources == 1) "single resource" else "none"
}

#' Read one result file written by run_coloc; NULL when no gene was tested (the file is empty).
read_result_file <- function(file, ...) {
  if (!file.exists(file) || file.size(file) == 0) return(NULL)
  one <- fread(file, ...)
  if (!nrow(one)) NULL else one
}

#' Default summary file for a locus, window and GWAS setting.
summary_file <- function(results_dir, label, window_kb, gwas_setting = "pan") {
  file.path(results_dir, sprintf("%s__summary_%skb%s.csv", label,
                                 format(window_kb, scientific = FALSE, trim = TRUE), gwas_suffix(gwas_setting)))
}

#' One row per gene across every result file for a locus, window and GWAS setting.
#'
#' Reads <results_dir>/<label>__<dataset>__<window>kb[__<gwas setting>].csv (files written by
#' run_coloc) and applies the reporting rule (see SUPPORT_THRESHOLD) to the primary datasets, counting
#' only tests whose outcome is "shared", so a test that was not testable can never count:
#' `support` under coloc's default prior and `support_strict_prior` under the stricter one.
#' Also: how many primary datasets were testable and shared, where the gene reads as two separate
#' signals, and whether it is shared in one dataset and two signals in another. Then, per dataset,
#' the outcome, PP.H4 under both priors, PP.H3, SNP count and top SNP. Replicated genes come first.
#'
#' @param annotation from load_duplication_annotation(); adds the duplicated-sequence columns.
#' @param comparison_window_kb if given, the support labels at that window are added as columns
#'   (`support_at_<n>kb`, `support_strict_prior_at_<n>kb`; "not tested" when a gene is absent there).
#' @param gwas_setting which GWAS setting's files to read ("pan", "eur", "pan_nfilter", ...).
summarise_locus <- function(results_dir, label, window_kb, threshold = SUPPORT_THRESHOLD,
                            sensitivity_only = SENSITIVITY_ONLY_DATASETS, annotation = NULL,
                            comparison_window_kb = NULL, gwas_setting = "pan") {
  suffix <- sprintf("__%skb%s.csv", format(window_kb, scientific = FALSE, trim = TRUE), gwas_suffix(gwas_setting))
  prefix <- paste0(label, "__")
  files <- list.files(results_dir, pattern = paste0(gsub(".", "[.]", suffix, fixed = TRUE), "$"))
  files <- files[startsWith(files, prefix) & !grepl("__summary_|coloc_on_ezQTL_inputs", files)]
  if (!length(files)) {
    stop(sprintf("no result files for %s at %s kb (%s) in %s", label, window_kb, gwas_setting, results_dir), call. = FALSE)
  }

  long <- rbindlist(lapply(files, function(f) {
    one <- read_result_file(file.path(results_dir, f))
    if (is.null(one)) return(NULL)
    if (!all(c("PP.H4_shared_strict_prior", "outcome", "outcome_strict_prior") %in% names(one))) {
      stop(sprintf("%s was written before outcomes were recorded; rerun it with run_coloc.R or run_locus.R", f), call. = FALSE)
    }
    one[, .(dataset = substr(f, nchar(prefix) + 1, nchar(f) - nchar(suffix)), gene_id,
            symbol = as.character(symbol), symbol_from,
            ids_sharing = as.integer(fuma_symbol_shared_by_n_gene_ids),
            H4 = PP.H4_shared, H4strict = PP.H4_shared_strict_prior, H3 = PP.H3_two_signals,
            outcome = as.character(outcome), outcome_strict = as.character(outcome_strict_prior), n_snps, top_snp)]
  }))
  if (!nrow(long)) {
    # Every dataset was run and none could test a gene (for example, too few SNPs left after a filter).
    return(data.table(symbol = character(), gene_id = character(), support = character(), support_strict_prior = character(),
                      n_datasets_shared = integer(), n_datasets_shared_strict_prior = integer(), datasets_at_threshold = character(),
                      datasets_two_signals = character(), shared_and_two_signals = logical(), n_primary_datasets_tested = integer(),
                      n_primary_datasets_testable = integer(), max_H4_primary = numeric(), best_primary_dataset = character()))
  }
  long[, `:=`(primary = !dataset %in% sensitivity_only, resource = resource_of(dataset),
              testable = !startsWith(outcome, "not testable"))]

  # One name per gene: the eQTL source's own name if any dataset gives one, else FUMA's, else blank.
  best_symbol <- function(symbol, symbol_from) {
    from_source <- symbol[symbol_from == "eqtl_source" & !is.na(symbol) & symbol != ""]
    any_name <- symbol[!is.na(symbol) & symbol != ""]
    if (length(from_source)) from_source[1] else if (length(any_name)) any_name[1] else ""
  }
  genes <- long[, {
    hit <- primary & outcome == "shared"
    hit_strict <- primary & outcome_strict == "shared"
    two <- primary & outcome == "two signals"
    .(symbol = best_symbol(symbol, symbol_from),
      fuma_symbol_shared_by_n_gene_ids = ids_sharing[!is.na(ids_sharing)][1],
      support = support_label(uniqueN(resource[hit])),
      support_strict_prior = support_label(uniqueN(resource[hit_strict])),
      n_datasets_shared = sum(hit), n_datasets_shared_strict_prior = sum(hit_strict),
      resources_at_threshold = paste(sort(unique(resource[hit])), collapse = "; "),
      datasets_at_threshold = paste(dataset[hit][order(-H4[hit])], collapse = "; "),
      datasets_two_signals = paste(dataset[two][order(-H3[two])], collapse = "; "),
      shared_and_two_signals = any(hit) && any(two),
      n_primary_datasets_tested = sum(primary), n_primary_datasets_testable = sum(primary & testable),
      max_H4_primary = if (any(primary)) max(H4[primary]) else NA_real_,
      best_primary_dataset = if (any(primary)) dataset[primary][which.max(H4[primary])] else "",
      n_datasets_tested = .N, n_H4_ge_0.5 = sum(H4 >= 0.5), n_H4_ge_0.8 = sum(H4 >= 0.8),
      n_H3_ge_0.5 = sum(H3 >= 0.5), max_H4 = max(H4), best_dataset = dataset[which.max(H4)])
  }, by = gene_id]

  wide <- dcast(long, gene_id ~ dataset, value.var = c("outcome", "H4", "H4strict", "H3", "n_snps", "top_snp"))
  per_dataset <- unlist(lapply(sort(unique(long$dataset)), function(d) {
    paste(c("outcome", "H4", "H4strict", "H3", "n_snps", "top_snp"), d, sep = "_")
  }))
  extra <- character()
  if (!is.null(comparison_window_kb)) {
    other <- summarise_locus(results_dir, label, comparison_window_kb, threshold, sensitivity_only, gwas_setting = gwas_setting)
    names_at <- paste0(c("support_at_", "support_strict_prior_at_"),
                       format(comparison_window_kb, scientific = FALSE, trim = TRUE), "kb")
    other <- other[, .(gene_id, support, support_strict_prior)]
    setnames(other, c("support", "support_strict_prior"), names_at)
    genes <- merge(genes, other, by = "gene_id", all.x = TRUE)
    for (column in names_at) set(genes, which(is.na(genes[[column]])), column, "not tested")
    extra <- c(extra, names_at)
  }
  if (!is.null(annotation)) {
    genes <- merge(genes, duplication_by_gene(genes$gene_id, annotation), by = "gene_id", all.x = TRUE)
    extra <- c(extra, "gene_span_hg19", "fraction_in_segmental_duplication", "max_duplication_identity")
  }
  out <- merge(genes, wide, by = "gene_id")
  lead_columns <- c("symbol", "gene_id", "support", "support_strict_prior", extra)
  setcolorder(out, c(lead_columns, setdiff(names(genes), lead_columns), per_dataset))
  tier <- match(out$support, c("replicated", "single resource", "none"))
  out[order(tier, -n_datasets_shared, -max_H4_primary, na.last = TRUE)]
}

# ---- One whole locus ----------------------------------------------------------------------------

#' The datasets a locus is run against, as a table of source strings and dataset labels.
#' Primary: GTEx v8 cortex, three arteries, tibial nerve, two adipose depots and pituitary,
#' PsychENCODE, the eight Bryois cell types. Comparison only: GTEx v7 cortex, and any eQTL
#' Catalogue slice on disk whose file name (..._chr<chr>_<from>-<to>Mb_hg38.tsv) holds the locus
#' with 1 Mb to spare on each side.
locus_sources <- function(paths, chr, lead_pos) {
  sources <- data.table(
    source = c("gtex_v8", paste0("gtex_v8:", GTEX_V8_REGION_TISSUES), "psychencode",
               paste0("bryois:", BRYOIS_CELL_TYPES), "gtex_v7"),
    dataset_label = NA_character_)
  if (!is.null(paths$catalogue_dir) && dir.exists(paths$catalogue_dir)) {
    pattern <- "^(QTD[0-9]+)_.*_chr([0-9XY]+)_([0-9.]+)-([0-9.]+)Mb_hg38[.]tsv$"
    for (f in list.files(paths$catalogue_dir, pattern = pattern)) {
      part <- regmatches(f, regexec(pattern, f))[[1]]
      holds_locus <- part[3] == as.character(chr) && lead_pos / 1e6 >= as.numeric(part[4]) + 1 &&
        lead_pos / 1e6 <= as.numeric(part[5]) - 1
      if (holds_locus && part[2] %in% names(CATALOGUE_SLICE_LABELS)) {
        sources <- rbind(sources, data.table(source = paste0("catalogue:", file.path(paths$catalogue_dir, f)),
                                             dataset_label = CATALOGUE_SLICE_LABELS[[part[2]]]))
      }
    }
  }
  sources
}

#' Run one locus end to end: every dataset at every window, then the summaries.
#'
#' A dataset that fails does not stop the others: its error is recorded, and any older result file
#' for it is removed so that a stale result cannot be summarised with the new ones.
#' Nothing is downloaded here; a GTEx v8 tissue whose region is missing fails with the reason.
#' @return list(failures = table of window, source and error; summaries = named summary files)
run_locus <- function(label, chr, lead_pos, paths, gwas_version = "pan", sample_size_filter = FALSE,
                      windows = c(PRIMARY_WINDOW_KB, COMPARISON_WINDOW_KB),
                      sources = locus_sources(paths, chr, lead_pos)) {
  chr <- as.character(chr)
  gwas_setting <- gwas_label(gwas_version, sample_size_filter)
  failures <- list()
  gwas_reads <- lapply(windows, function(window_kb) {
    window <- locus_window(lead_pos, window_kb)
    read_gwas_window(gwas_file_for(paths, gwas_version), chr, window$lo, window$hi)
  })
  for (i in seq_len(nrow(sources))) {     # dataset by dataset, so both windows reuse one read of a large eQTL file
    source <- sources$source[i]
    dataset_label <- if (is.na(sources$dataset_label[i])) default_dataset_label(source) else sources$dataset_label[i]
    for (w in seq_along(windows)) {
      window_kb <- windows[w]
      gwas_read <- gwas_reads[[w]]
      error <- tryCatch({
        suppressMessages(run_coloc(label, chr, lead_pos, window_kb, source, dataset_label = dataset_label, paths = paths,
                                   gwas_version = gwas_version, sample_size_filter = sample_size_filter,
                                   gwas_read = gwas_read))
        NA_character_
      }, error = function(e) conditionMessage(e))
      if (!is.na(error)) {
        unlink(result_file(paths$results, label, dataset_label, window_kb, gwas_setting))
        failures[[length(failures) + 1]] <- data.table(window_kb = window_kb, dataset = dataset_label, error = error)
      }
      message(sprintf("%4s kb  %-36s %s", window_kb, dataset_label, if (is.na(error)) "done" else paste("FAILED:", error)))
    }
  }
  annotation <- load_duplication_annotation(paths)
  summaries <- character()
  for (window_kb in windows) {
    comparison <- if (window_kb == PRIMARY_WINDOW_KB && COMPARISON_WINDOW_KB %in% windows) COMPARISON_WINDOW_KB else NULL
    table <- summarise_locus(paths$results, label, window_kb, annotation = annotation,
                             comparison_window_kb = comparison, gwas_setting = gwas_setting)
    file <- summary_file(paths$results, label, window_kb, gwas_setting)
    fwrite(table, file)
    summaries[as.character(window_kb)] <- file
  }
  failures <- rbindlist(failures)
  if (!nrow(failures)) failures <- data.table(window_kb = numeric(), dataset = character(), error = character())
  list(failures = failures, summaries = summaries)
}

# ---- All loci -----------------------------------------------------------------------------------

#' The loci of the FUMA job: locus number, file-name label, chromosome, lead position, lead SNP, p.
locus_table <- function(paths) {
  check_file(paths$fuma_loci, "FUMA GenomicRiskLoci file")
  loci <- fread(paths$fuma_loci)
  needed <- c("GenomicLocus", "rsID", "chr", "pos", "p")
  if (!all(needed %in% names(loci))) {
    stop("FUMA GenomicRiskLoci file: columns GenomicLocus, rsID, chr, pos and p are needed", call. = FALSE)
  }
  out <- loci[, .(locus = as.integer(GenomicLocus), chr = as.character(chr), lead_pos = pos, lead_snp = rsID, gwas_p = p)]
  out[, label := sprintf("L%d_chr%s_%.1fMb", locus, chr, lead_pos / 1e6)]
  named <- match(as.character(out$locus), names(LOCUS_LABEL_OVERRIDES))
  out[!is.na(named), label := LOCUS_LABEL_OVERRIDES[named[!is.na(named)]]]
  out[]
}

#' One row per locus: what is known about it and what its colocalization run found (primary window).
#' A locus whose summary file does not exist is listed with status "not run".
overview_all_loci <- function(paths, loci = locus_table(paths), gwas_setting = "pan", window_kb = PRIMARY_WINDOW_KB) {
  primary_expected <- locus_sources(list(), "1", 1)[source != "gtex_v7", .N]
  name <- function(symbol, gene_id) ifelse(is.na(symbol) | symbol == "", gene_id, symbol)
  rows <- lapply(seq_len(nrow(loci)), function(i) {
    locus <- loci[i]
    context <- locus_context(paths, locus$chr, locus$lead_pos, window_kb)
    base <- data.table(locus = locus$locus, label = locus$label, chr = locus$chr, lead_snp = locus$lead_snp,
                       lead_pos_hg19 = locus$lead_pos, gwas_p_pan = locus$gwas_p,
                       gwas_lead_p_european_only = as.character(context$gwas_lead_p_european_only),
                       european_signal = as.character(context$european_signal),
                       one_signal_assumption = as.character(context$one_signal_assumption),
                       window_covers_fuma_locus = as.character(context$window_covers_fuma_locus))
    file <- summary_file(paths$results, locus$label, window_kb, gwas_setting)
    if (!file.exists(file)) return(cbind(base, status = "not run"))
    genes <- fread(file)
    suffix <- sprintf("__%skb%s.csv", format(window_kb, scientific = FALSE, trim = TRUE), gwas_suffix(gwas_setting))
    prefix <- paste0(locus$label, "__")
    files <- list.files(paths$results, pattern = paste0(gsub(".", "[.]", suffix, fixed = TRUE), "$"))
    files <- files[startsWith(files, prefix) & !grepl("__summary_|coloc_on_ezQTL_inputs", files)]
    datasets <- substr(files, nchar(prefix) + 1, nchar(files) - nchar(suffix))
    is_primary <- !datasets %in% SENSITIVITY_ONLY_DATASETS
    tests <- rbindlist(lapply(files[is_primary], function(f) read_result_file(file.path(paths$results, f), select = "outcome")), fill = TRUE)
    outcome <- if (nrow(tests)) as.character(tests$outcome) else character()
    listed <- function(keep) paste(name(genes$symbol[keep], genes$gene_id[keep]), collapse = "; ")
    best <- if (nrow(genes) && any(!is.na(genes$max_H4_primary))) which.max(genes$max_H4_primary) else NA
    cbind(base, data.table(
      status = "run", primary_datasets_run = sum(is_primary), primary_datasets_expected = primary_expected,
      genes_tested = nrow(genes), gene_tests = length(outcome), testable_tests = sum(!startsWith(outcome, "not testable")),
      shared_tests = sum(outcome == "shared"), two_signal_tests = sum(outcome == "two signals"),
      leaning_tests = sum(outcome == "leaning"),
      n_replicated = sum(genes$support == "replicated"), genes_replicated = listed(genes$support == "replicated"),
      n_replicated_strict_prior = sum(genes$support_strict_prior == "replicated"),
      genes_replicated_strict_prior = listed(genes$support_strict_prior == "replicated"),
      n_single_resource = sum(genes$support == "single resource"), genes_single_resource = listed(genes$support == "single resource"),
      genes_shared_and_two_signals = listed(genes$shared_and_two_signals %in% TRUE),
      top_gene = if (is.na(best)) "" else name(genes$symbol[best], genes$gene_id[best]),
      top_gene_H4 = if (is.na(best)) NA_real_ else round(genes$max_H4_primary[best], 3),
      top_gene_dataset = if (is.na(best)) "" else genes$best_primary_dataset[best]))
  })
  rbindlist(rows, fill = TRUE)
}

#' Every gene with a shared call at any locus (support is not "none"), one row per gene, loci in order.
supported_genes_all_loci <- function(paths, loci = locus_table(paths), gwas_setting = "pan", window_kb = PRIMARY_WINDOW_KB) {
  wanted <- c("symbol", "gene_id", "support", "support_strict_prior", paste0("support_at_", COMPARISON_WINDOW_KB, "kb"),
              "n_datasets_shared", "datasets_at_threshold", "datasets_two_signals", "shared_and_two_signals",
              "n_primary_datasets_testable", "max_H4_primary", "best_primary_dataset", "gene_span_hg19",
              "fraction_in_segmental_duplication", "max_duplication_identity")
  rows <- lapply(seq_len(nrow(loci)), function(i) {
    file <- summary_file(paths$results, loci$label[i], window_kb, gwas_setting)
    if (!file.exists(file)) return(NULL)
    genes <- fread(file)[support != "none"]
    if (!nrow(genes)) return(NULL)
    cbind(data.table(locus = loci$locus[i], label = loci$label[i], lead_snp = loci$lead_snp[i]), genes[, intersect(wanted, names(genes)), with = FALSE])
  })
  out <- rbindlist(rows, fill = TRUE)
  if (!nrow(out)) out <- data.table(locus = integer(), label = character(), lead_snp = character(), symbol = character(), gene_id = character(), support = character())
  out
}
