# Tests for coloc_lib.R. Plain R, no test package needed.
#   Rscript code/tests/test_coloc_lib.R
# Uses small made-up files in a temporary folder; never touches the real data.

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- if (length(script)) dirname(dirname(normalizePath(script[1]))) else file.path(getwd(), "code")
source(file.path(code_dir, "coloc_lib.R"))

# Made-up study settings, so the tests never depend on code/study_config.R.
STUDY <- list(n_cases = 21000, n_controls = 84000)
GWAS_VERSIONS <- list(pan = list(file = "gwas_pan.tsv.gz", n_cases = STUDY$n_cases, n_controls = STUDY$n_controls),
                      eur = list(file = "gwas_eur.tsv.gz", n_cases = 15000, n_controls = 60000))
LOCUS_LABEL_OVERRIDES <- c(`7` = "L7_named")
STUDY_FILES <- list(fuma_job_dir = "FUMA_job", gene_evidence = "gene-evidence.csv", ancestry_table = "ancestry.csv",
                    gtex_v7_slice = "gtex_v7_slice.tsv", psychencode_slice = "psychencode_slice.txt")

passed <- 0L
skipped_tests <- 0L
failed <- character()
test <- function(name, expr) {
  ok <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) {
    failed <<- c(failed, sprintf("%s: %s", name, conditionMessage(e)))
    FALSE
  })
  if (ok) passed <<- passed + 1L
  cat(sprintf("%s  %s\n", if (ok) "ok  " else "FAIL", name))
}
skip_test <- function(name, why) {
  skipped_tests <<- skipped_tests + 1L
  cat(sprintf("skip  %s (%s)\n", name, why))
}
expect_error <- function(expr, pattern) {
  msg <- tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(msg)) stop("expected an error, got none")
  if (!grepl(pattern, msg)) stop(sprintf("error did not match '%s': %s", pattern, msg))
}

# data.table on this machine cannot write .gz itself, so write plain text and compress it.
write_gz <- function(table, path, ...) {
  plain <- sub("\\.gz$", "", path)
  fwrite(table, plain, ...)
  if (system2("gzip", c("-f", shQuote(plain))) != 0) stop("gzip failed for ", plain)
}

tmp <- file.path(tempdir(), "coloc_lib_tests")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)

# ---- Made-up locus: 200 SNPs on chromosome 16, one causal SNP for the GWAS -----------------------
set.seed(1)
n_snp <- 200
positions <- 28000000L + seq_len(n_snp) * 1000L
alleles_1 <- sample(c("A", "C"), n_snp, replace = TRUE)
alleles_2 <- sample(c("G", "T"), n_snp, replace = TRUE)
rsids <- paste0("rs", 1000 + seq_len(n_snp))
maf <- runif(n_snp, 0.1, 0.5)
ld_decay <- function(causal) exp(-abs(seq_len(n_snp) - causal) / 4)

gwas_se <- rep(0.01, n_snp)
gwas_z <- 8 * ld_decay(100) + rnorm(n_snp, sd = 0.3)
gwas_table <- data.table(CHR = "16", BP = positions, A1 = alleles_1, A2 = alleles_2, rsID = rsids,
                         P = 2 * pnorm(-abs(gwas_z)), Beta = gwas_z * gwas_se, SE = gwas_se, N = 648172)
gwas_file <- file.path(tmp, "gwas.tsv.gz")
write_gz(gwas_table, gwas_file, sep = "\t")

eqtl_z <- function(causal) 7 * ld_decay(causal) + rnorm(n_snp, sd = 0.3)
z_shared <- eqtl_z(100)    # same causal SNP as the GWAS
z_distinct <- eqtl_z(30)   # a different causal SNP
eqtl_se <- 0.1

gtex_file <- file.path(tmp, "gtex.tsv")
gtex_rows <- function(gene, z) {
  data.table(gene_id = gene, variant_id = paste("16", positions, alleles_1, alleles_2, "b37", sep = "_"),
             tss_distance = 0L, ma_samples = 50L, ma_count = 60L, maf = maf,
             pval_nominal = 2 * pnorm(-abs(z)), slope = z * eqtl_se, slope_se = eqtl_se)
}
fwrite(rbind(gtex_rows("ENSG00000000001.5", z_shared), gtex_rows("ENSG00000000002.7", z_distinct),
             gtex_rows("ENSG00000000003.1", z_shared)[1:10]), gtex_file, sep = "\t")

psychencode_file <- file.path(tmp, "psychencode.txt")
fwrite(data.table("ENSG00000000001.5", "chr16", 28100000L, 28100000L, "+", 200L, 0L,
                  paste0("16:", positions), "chr16", positions, positions,
                  2 * pnorm(-abs(z_shared)), z_shared * eqtl_se, 0L),
       psychencode_file, sep = " ", col.names = FALSE)

bryois_dir <- file.path(tmp, "bryois")
dir.create(bryois_dir)
write_gz(data.table(paste0("GENEA_ENSG00000000001"), rsids, 0L, 2 * pnorm(-abs(z_shared)), z_shared * eqtl_se),
       file.path(bryois_dir, "Testcells.16.gz"), sep = " ", col.names = FALSE)
write_gz(data.table(SNP = rsids, SNP_id_hg38 = paste0("chr16:", positions + 5L),
                  SNP_id_hg19 = paste0("chr16:", positions), effect_allele = alleles_2,
                  other_allele = alleles_1, MAF = maf),
       file.path(bryois_dir, "snp_pos.txt.gz"), sep = "\t")

catalogue_file <- file.path(tmp, "catalogue.tsv")
flank_positions <- c(27000000L + 1:5 * 1000L, 29500000L + 1:5 * 1000L)   # far outside the GWAS window
catalogue_rows <- function(gene, pos, a1, a2, z, ids, maf_values) {
  data.table(gene, "16", pos, a1, a2, paste("chr16", pos, a1, a2, sep = "_"), 50L, maf_values,
             2 * pnorm(-abs(z)), z * eqtl_se, eqtl_se, "SNP", 60L, 410L, NA, gene, gene, 10.5, ids)
}
fwrite(rbind(catalogue_rows("ENSG00000000001", positions + 5L, alleles_1, alleles_2, z_shared, rsids, maf),
             catalogue_rows("ENSG00000000008", flank_positions, "A", "G", rep(1, 10), paste0("rs9", 1:10), rep(0.3, 10))),
       catalogue_file, sep = "\t", col.names = FALSE)

# GTEx v8 look-alike: GRCh38 positions 1.2 Mb away from hg19, header row, bgzip + tabix index.
build_shift <- 1200000L
tabix_bin <- tryCatch(find_tabix(), error = function(e) "")
bgzip_bin <- file.path(dirname(tabix_bin), "bgzip")
have_htslib <- nzchar(tabix_bin) && file.exists(bgzip_bin)
gtex_v8_file <- file.path(tmp, "gtex_v8.tsv.gz")
v8_rows <- function(gene, pos, a1, a2, z, ids, maf_values) {
  data.table(variant = paste("chr16", pos, a1, a2, sep = "_"), r2 = NA, pvalue = 2 * pnorm(-abs(z)),
             molecular_trait_object_id = gene, molecular_trait_id = gene, maf = maf_values, gene_id = gene,
             median_tpm = 5.5, beta = z * eqtl_se, se = eqtl_se, an = 410L, ac = 60L, chromosome = "16",
             position = pos, ref = a1, alt = a2, type = "SNP", rsid = ids)
}
if (have_htslib) {
  v8 <- rbind(v8_rows("ENSG00000000001", positions + build_shift, alleles_1, alleles_2, z_shared, rsids, maf),
              v8_rows("ENSG00000000002", positions + build_shift, alleles_1, alleles_2, z_distinct, rsids, maf),
              v8_rows("ENSG00000000007", 32000000L + 1:60 * 1000L, "A", "G", rep(1, 60), paste0("rs7", 1:60), rep(0.3, 60)))
  setorder(v8, position)
  plain <- sub("[.]gz$", "", gtex_v8_file)
  fwrite(v8, plain, sep = "\t")
  stopifnot(system2(bgzip_bin, c("-f", shQuote(plain))) == 0,
            system2(tabix_bin, c("-s", "13", "-b", "14", "-e", "14", "-S", "1", shQuote(gtex_v8_file))) == 0)
}

magma_file <- file.path(tmp, "magma.genes.out")
fwrite(data.table(GENE = c("ENSG00000000001", "ENSG00000000002", "ENSG00000000009"),
                  SYMBOL = c("TWIN", "SOLO", "TWIN")), magma_file, sep = "\t")

fuma_loci_file <- file.path(tmp, "GenomicRiskLoci.txt")
fwrite(data.table(GenomicLocus = 1:3, rsID = c("rs1100", "rsB", "rsC"), p = c(1e-15, 2e-9, 3e-8), chr = c(16L, 16L, 2L), pos = c(positions[100], 5000000L, 7000000L),
                  start = c(positions[100] - 100000L, 4990000L, 6990000L), end = c(positions[100] + 300000L, 5010000L, 7010000L),
                  nIndSigSNPs = c(1L, 4L, 2L), nLeadSNPs = c(1L, 1L, 2L)), fuma_loci_file, sep = "\t")
ancestry_file <- file.path(tmp, "ancestry.csv")
fwrite(data.table(locus = 1:3, chr = c(16L, 16L, 2L), lead_pos_hg19 = c(positions[100], 5000000L, 7000000L),
                  EUR_lead_p = c(2e-7, 3e-9, NA)), ancestry_file)
regions_dir <- file.path(tmp, "gtex_v8_regions")

paths <- list(gwas = gwas_file, gtex_v7 = gtex_file, gtex_v8 = gtex_v8_file, gtex_v8_regions = regions_dir,
              fuma_loci = fuma_loci_file, ancestry = ancestry_file, psychencode = psychencode_file,
              bryois_dir = bryois_dir, magma = magma_file, results = file.path(tmp, "results"))
lo <- 28000000
hi <- 28300000

# ---- Helpers ------------------------------------------------------------------------------------
test("allele_key ignores allele order", {
  stopifnot(identical(allele_key("A", "G"), allele_key("G", "A")),
            allele_key("A", "G") != allele_key("A", "T"),
            identical(allele_key(c("T", "C"), c("A", "CA")), c("A T", "C CA")))
})
test("se_from_p inverts a two-sided p-value", {
  beta <- c(0.2, -0.5)
  se <- c(0.05, 0.1)
  stopifnot(isTRUE(all.equal(se_from_p(beta, 2 * pnorm(-abs(beta / se))), se)))
})
test("check_columns stops on a wrong header", {
  check_columns(GWAS_COLUMNS, GWAS_COLUMNS, "x")
  expect_error(check_columns(rev(GWAS_COLUMNS), GWAS_COLUMNS, "GWAS file"), "GWAS file: unexpected columns")
})
test("locus_window validates and handles fractional kb", {
  stopifnot(identical(locus_window(1000000, 0.5), list(lo = 999500, hi = 1000500)),
            identical(locus_window(1000, 250)$lo, 0))
  expect_error(locus_window(NA_real_, 250), "lead_pos")
  expect_error(locus_window(1000, -1), "window_kb")
})
test("format_position never uses scientific notation", {
  stopifnot(identical(format_position(100000), "100000"), identical(format_position(28500001), "28500001"))
})
test("parse_source and default labels", {
  stopifnot(identical(parse_source("bryois:OPCs...COPs"), list(kind = "bryois", arg = "OPCs...COPs")),
            identical(parse_source("gtex_v8"), list(kind = "gtex_v8", arg = "")),
            identical(parse_source("gtex_v8:Artery_Tibial"), list(kind = "gtex_v8", arg = "Artery_Tibial")),
            identical(default_dataset_label("gtex_v8:Artery_Tibial"), "GTEx_v8_Artery_Tibial"),
            identical(default_dataset_label("gtex_v8"), "GTEx_v8_brain_cortex"),
            identical(default_dataset_label("bryois:OPCs...COPs"), "Bryois_OPCs_COPs"),
            identical(default_dataset_label("bryois:Endothelial.cells"), "Bryois_Endothelial.cells"),
            identical(default_dataset_label("gtex_v7"), "GTEx_v7_brain_cortex"))
  expect_error(parse_source("gtex_v9"), "unknown source")
  expect_error(parse_source("bryois"), "needs an argument")
  expect_error(default_dataset_label("catalogue:x.tsv"), "dataset_label")
})

test("normalise_warnings groups warnings that differ only by a number", {
  w <- c("minimum p value is: 0.00017192\nIf this is what you expected", "minimum p value is: 1.2e-05\nmore", "other")
  stopifnot(identical(normalise_warnings(w), c("minimum p value is: <n>", "minimum p value is: <n>", "other")),
            identical(normalise_warnings(character()), character()))
})

# ---- Readers ------------------------------------------------------------------------------------
test("read_gwas_window keeps only the window and checks the header", {
  g <- read_gwas_window(gwas_file, "16", 28050000, 28060000)
  stopifnot(nrow(g$data) == 11, all(g$data$pos >= 28050000 & g$data$pos <= 28060000),
            identical(names(g$data), c("pos", "rsid", "ak", "gwas_beta", "gwas_se", "gwas_p", "gwas_n")),
            nrow(read_gwas_window(gwas_file, "17", 28050000, 28060000)$data) == 0)
  bad <- file.path(tmp, "bad_gwas.tsv")
  fwrite(setnames(copy(gwas_table), "Beta", "BETA"), bad, sep = "\t")
  expect_error(read_gwas_window(bad, "16", lo, hi), "unexpected columns")
  expect_error(read_gwas_window(file.path(tmp, "nope.tsv.gz"), "16", lo, hi), "not found")
})
test("read_gwas_window drops rows without a usable SE", {
  broken <- copy(gwas_table)
  broken[1:3, SE := c(0, NA, -1)]
  f <- file.path(tmp, "gwas_broken.tsv")
  fwrite(broken, f, sep = "\t")
  g <- read_gwas_window(f, "16", lo, hi)
  stopifnot(g$counts[["gwas_rows_in_window"]] == 200, g$counts[["gwas_rows_with_beta_and_se"]] == 197)
})
test("read_eqtl_gtex_v7 standardises columns and strips gene versions", {
  e <- read_eqtl_gtex_v7(gtex_file, "16", lo, hi)
  stopifnot(identical(names(e$data), names(empty_eqtl())), e$spec$match_by == "pos_alleles",
            e$spec$sd_y == 1, setequal(unique(e$data$gene_id), sprintf("ENSG0000000000%d", 1:3)),
            nrow(read_eqtl_gtex_v7(gtex_file, "16", 28050000, 28060000)$data) == 11 * 2)
})
test("read_eqtl_psychencode recovers SE and borrows MAF", {
  e <- read_eqtl_psychencode(psychencode_file, gtex_file, "16", lo, hi)
  stopifnot(nrow(e$data) == 200, isTRUE(all.equal(e$data$se, rep(eqtl_se, 200))),
            isTRUE(all.equal(e$data[order(pos)]$maf, maf)), e$spec$match_by == "pos", e$spec$n == PSYCHENCODE_N,
            nrow(read_eqtl_psychencode(psychencode_file, gtex_file, "5", lo, hi)$data) == 0)
  bad <- file.path(tmp, "bad_psychencode.txt")
  fwrite(fread(psychencode_file, header = FALSE)[, 1:13], bad, sep = " ", col.names = FALSE)
  expect_error(read_eqtl_psychencode(bad, gtex_file, "16", lo, hi), "expected 14 columns")
})
test("read_eqtl_bryois takes names, positions and alleles from the right places", {
  e <- read_eqtl_bryois(bryois_dir, "Testcells", "16", lo, hi)
  stopifnot(nrow(e$data) == 200, all(e$data$gene_id == "ENSG00000000001"), all(e$data$gene_name == "GENEA"),
            setequal(e$data$pos, positions), e$spec$n == BRYOIS_N, e$spec$match_by == "rsid_alleles")
  expect_error(read_eqtl_bryois(bryois_dir, "Nocells", "16", lo, hi), "not found")
})
test("read_eqtl_catalogue reads N from the allele number", {
  e <- read_eqtl_catalogue(catalogue_file)
  stopifnot(nrow(e$data) == 210, e$spec$n == 205, e$spec$match_by == "rsid_alleles",
            e$spec$positions_build == "GRCh38", identical(e$spec$region_read, c(27001000L, 29505000L)))
  bad <- file.path(tmp, "bad_catalogue.tsv")
  fwrite(fread(catalogue_file, header = FALSE)[, 1:18], bad, sep = "\t", col.names = FALSE)
  expect_error(read_eqtl_catalogue(bad), "expected 19 columns")
})

if (have_htslib) {
  test("read_eqtl_gtex_v8 reads a padded region through the index and checks the header", {
    e <- read_eqtl_gtex_v8(gtex_v8_file, "16", lo, hi)
    stopifnot(nrow(e$data) == 460, identical(names(e$data), names(empty_eqtl())), e$spec$sd_y == 1, e$spec$n == 205,
              e$spec$match_by == "rsid_alleles", identical(e$spec$region_read, c(lo - 5e6, hi + 5e6)),
              setequal(e$data[gene_id == "ENSG00000000001", pos], positions + build_shift),
              nrow(read_eqtl_gtex_v8(gtex_v8_file, "16", lo, hi, pad_mb = 0.5)$data) == 0,
              nrow(read_eqtl_gtex_v8(gtex_v8_file, "5", lo, hi)$data) == 0)
    bad <- file.path(tmp, "bad_v8.tsv.gz")
    writeLines(c(paste(rev(GTEX_V8_COLUMNS), collapse = "\t"), "x"), sub("[.]gz$", "", bad))
    system2("gzip", c("-f", shQuote(sub("[.]gz$", "", bad))))
    file.copy(paste0(gtex_v8_file, ".tbi"), paste0(bad, ".tbi"))
    expect_error(read_eqtl_gtex_v8(bad, "16", lo, hi), "GTEx v8 file: unexpected columns")
    expect_error(read_eqtl_gtex_v8(file.path(tmp, "nope.tsv.gz"), "16", lo, hi), "not found")
  })
  test("restrict_to_gwas_window keeps window SNPs and refuses a region that cuts the window", {
    g <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file, "16", lo, hi)$data, "rsid_alleles"))
    e <- add_match_key(read_eqtl_gtex_v8(gtex_v8_file, "16", lo, hi)$data, "rsid_alleles")
    r <- restrict_to_gwas_window(e, g, c(lo - 5e6, hi + 5e6))
    span <- range(positions + build_shift)
    stopifnot(nrow(r$eqtl) == 400, r$info$genes_with_no_snp_in_gwas_window == 1,
              r$info$eqtl_rows_for_snps_outside_gwas_window == 60,
              r$info$eqtl_span_of_matched_snps == sprintf("%d-%d", span[1], span[2]),
              nrow(restrict_to_gwas_window(e, g, c(1, hi + 5e6))$eqtl) == 400)   # a region starting at 1 is not an edge
    expect_error(restrict_to_gwas_window(e, g, c(span[1] - 10, span[2] + 5e6)), "edge of the eQTL data available")
    expect_error(restrict_to_gwas_window(e, g, c(span[1] - 5e6, span[2] + 10)), "edge of the eQTL data available")
    expect_error(restrict_to_gwas_window(e[gene_id == "ENSG00000000007"], g, c(lo - 5e6, hi + 5e6)), "no SNP in the GWAS window")
  })
  test("run_coloc on the GTEx v8 format finds the shared signal and logs the cross-build facts", {
    r <- run_coloc("V", "16", positions[100], 150, "gtex_v8", paths = paths)
    stopifnot(r$results[gene_id == "ENSG00000000001", PP.H4_shared] > 0.9,
              r$results[gene_id == "ENSG00000000002", PP.H4_shared] < 0.1, nrow(r$skipped) == 0,
              basename(r$out_file) == "V__GTEx_v8_brain_cortex__150kb.csv",
              r$settings$genes_with_no_snp_in_gwas_window == 1, r$settings$gwas_snps_with_an_eqtl_match == 200,
              r$settings$eqtl_n == 205, r$settings$eqtl_sd_y == 1, r$settings$genes_tested == 2)
  })
} else {
  skip_test("GTEx v8 reader tests", "bgzip/tabix not found")
}

# ---- Matching -----------------------------------------------------------------------------------
test("add_match_key and drop_ambiguous_keys", {
  d <- data.table(pos = c(10L, 10L, 20L, 30L), rsid = c("rs1", "rs1", ".", "rs3"), ak = c("A G", "A T", "A G", "C T"))
  stopifnot(identical(add_match_key(d, "pos_alleles")$key, c("10 A G", "10 A T", "20 A G", "30 C T")),
            identical(drop_ambiguous_keys(add_match_key(d, "pos"))$key, c("20", "30")),
            identical(add_match_key(d, "rsid_alleles")$key, c("rs1 A G", "rs1 A T", "rs3 C T")),
            !"key" %in% names(d))
  expect_error(add_match_key(d, "nonsense"), "arg")
})
test("load_symbol_map counts gene IDs that share a symbol", {
  m <- load_symbol_map(magma_file)
  stopifnot(m[gene_id == "ENSG00000000001", ids_sharing_symbol] == 2, m[symbol == "SOLO", ids_sharing_symbol] == 1,
            nrow(load_symbol_map(NULL)) == 0)
})

# ---- coloc --------------------------------------------------------------------------------------
test("coloc_locus finds the shared signal, rejects the distinct one, logs the skipped gene", {
  g <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file, "16", lo, hi)$data, "pos_alleles"))
  e_read <- read_eqtl_gtex_v7(gtex_file, "16", lo, hi)
  fit <- coloc_locus(g, add_match_key(e_read$data, "pos_alleles"), e_read$spec, lead_pos = positions[100],
                     symbol_map = load_symbol_map(magma_file))
  shared <- fit$results[gene_id == "ENSG00000000001"]
  distinct <- fit$results[gene_id == "ENSG00000000002"]
  stopifnot(shared$PP.H4_shared > 0.9, distinct$PP.H4_shared < 0.1, distinct$PP.H3_two_signals > 0.9,
            shared$n_snps == 200, shared$lead_snp_included,
            shared$symbol == "TWIN", shared$fuma_symbol_shared_by_n_gene_ids == 2,
            nrow(fit$skipped) == 1, fit$skipped$gene_id == "ENSG00000000003", fit$skipped$n_snps == 10,
            grepl("fewer than 50", fit$skipped$reason),
            isTRUE(all.equal(rowSums(fit$results[, .(PP.H0, PP.H1_gwas_only, PP.H2_eqtl_only,
                                                     PP.H3_two_signals, PP.H4_shared)]), c(1, 1))),
            # the stricter prior can only lower the shared-signal probability, and must match a direct call
            all(fit$results$PP.H4_shared_strict_prior < fit$results$PP.H4_shared),
            all(c("PP.H3_two_signals_strict_prior", "PP.H4_shared_strict_prior") %in% names(fit$results)))
  m <- merge(g, add_match_key(e_read$data[gene_id == "ENSG00000000001"], "pos_alleles")[, .(key, beta, se)], by = "key")
  invisible(capture.output(direct <- suppressWarnings(coloc.abf(
    list(beta = m$gwas_beta, varbeta = m$gwas_se^2, snp = m$key, type = "cc", s = STUDY$n_cases / (STUDY$n_cases + STUDY$n_controls), N = STUDY$n_cases + STUDY$n_controls),
    list(beta = m$beta, varbeta = m$se^2, snp = m$key, type = "quant", sdY = 1), p12 = 5e-6))))
  stopifnot(isTRUE(all.equal(shared$PP.H4_shared_strict_prior, direct$summary[["PP.H4.abf"]])),
            COLOC_P12_DEFAULT == 1e-5, COLOC_P12_STRICT == 5e-6)
})
test("coloc_locus records a coloc error instead of hiding it", {
  g <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file, "16", lo, hi)$data, "pos_alleles"))
  e_read <- read_eqtl_gtex_v7(gtex_file, "16", lo, hi)
  bad_spec <- modifyList(e_read$spec, list(sd_y = NULL, n = NULL))   # neither sdY nor N: coloc must refuse
  fit <- coloc_locus(g, add_match_key(e_read$data[gene_id == "ENSG00000000001"], "pos_alleles"), bad_spec)
  stopifnot(nrow(fit$results) == 0, nrow(fit$skipped) == 1, grepl("^coloc.abf error:", fit$skipped$reason))
})
test("run_coloc agrees across all four source formats and writes its side files", {
  runs <- list(
    run_coloc("T", "16", positions[100], 150, "gtex_v7", paths = paths),
    run_coloc("T", "16", positions[100], 150, "psychencode", paths = paths),
    run_coloc("T", "16", positions[100], 150, "bryois:Testcells", paths = paths),
    run_coloc("T", "16", positions[100], 150, paste0("catalogue:", catalogue_file), dataset_label = "Cat", paths = paths))
  h4 <- vapply(runs, function(r) r$results[gene_id == "ENSG00000000001", PP.H4_shared], 0)
  stopifnot(all(h4 > 0.9), runs[[3]]$results$symbol == "GENEA", runs[[3]]$results$symbol_from == "eqtl_source",
            runs[[1]]$results[gene_id == "ENSG00000000001", symbol_from] == "fuma_magma")
  out <- runs[[1]]$out_file
  stopifnot(basename(out) == "T__GTEx_v7_brain_cortex__150kb.csv", file.exists(out),
            file.exists(file.path(dirname(out), "run_logs", "T__GTEx_v7_brain_cortex__150kb.skipped.csv")))
  settings <- fread(file.path(dirname(out), "run_logs", "T__GTEx_v7_brain_cortex__150kb.settings.tsv"))
  value <- function(name) settings[setting == name, value]
  stopifnot(value("n_cases") == as.character(STUDY$n_cases), value("window_lo") == "27950000", value("genes_tested") == "2",
            value("genes_skipped") == "1", value("lead_rsid") == "rs1100", value("matched_by") == "pos_alleles",
            value("prior_p12_default") == "1e-05", value("prior_p12_strict") == "5e-06", value("window_role") == "other",
            value("fuma_locus_reach_from_lead_kb") == "300", value("window_covers_fuma_locus") == "FALSE",
            value("one_signal_assumption") == "ok: one independent significant SNP",
            value("european_signal") == "suggestive in Europeans (5e-8 to 1e-5)")
  written <- fread(out)
  stopifnot(nrow(written) == 2, written$PP.H4_shared[1] >= written$PP.H4_shared[2])
})

test("summarise_locus gives one row per gene across datasets", {
  for (src in c("gtex_v7", "psychencode", "bryois:Testcells")) run_coloc("S", "16", positions[100], 150, src, paths = paths)
  s <- summarise_locus(paths$results, "S", 150)
  top <- s[gene_id == "ENSG00000000001"]
  stopifnot(nrow(s) == 2, s$gene_id[1] == "ENSG00000000001", top$n_datasets_tested == 3, top$n_H4_ge_0.8 == 3,
            top$symbol == "GENEA", top$fuma_symbol_shared_by_n_gene_ids == 2,
            s[gene_id == "ENSG00000000002", n_datasets_tested] == 1,
            all(c("H4_GTEx_v7_brain_cortex", "H3_PsychENCODE", "top_snp_Bryois_Testcells",
                  "H4strict_PsychENCODE") %in% names(s)),
            top$support == "replicated", top$resources_at_threshold == "Bryois2022; PsychENCODE",   # GTEx v7 does not vote
            top$n_datasets_shared == 2, top$n_primary_datasets_testable == 2, !top$shared_and_two_signals,
            "outcome_PsychENCODE" %in% names(s), top$outcome_GTEx_v7_brain_cortex == "shared",
            s[gene_id == "ENSG00000000002", support] == "none")
  expect_error(summarise_locus(paths$results, "S", 999), "no result files")
})

test("summarise_locus names a gene from whichever dataset names it", {
  dir <- file.path(tmp, "naming"); dir.create(dir)
  row <- function(symbol, symbol_from) {
    data.table(locus = "N", eqtl_dataset = "x", gene_id = "ENSG00000000042", symbol = symbol, symbol_from = symbol_from,
               fuma_symbol = "", fuma_symbol_shared_by_n_gene_ids = NA_integer_, n_snps = 100L,
               PP.H3_two_signals = 0.1, PP.H4_shared = 0.6, PP.H4_shared_strict_prior = 0.5, outcome = "leaning",
               outcome_strict_prior = "leaning", top_snp = "rs1")
  }
  fwrite(row("", ""), file.path(dir, "N__Unnamed__100kb.csv"))
  fwrite(row("NAMED", "eqtl_source"), file.path(dir, "N__Named__100kb.csv"))
  stopifnot(summarise_locus(dir, "N", 100)$symbol == "NAMED")
  fwrite(row("FUMANAME", "fuma_magma"), file.path(dir, "N__Unnamed__100kb.csv"))
  stopifnot(summarise_locus(dir, "N", 100)$symbol == "NAMED")
  file.remove(file.path(dir, "N__Named__100kb.csv"))
  stopifnot(summarise_locus(dir, "N", 100)$symbol == "FUMANAME")
})

# Hand-made result files for the summary tests. The outcome follows the probabilities unless given.
result_rows <- function(locus, dataset, genes, h4, h4_strict = h4 - 0.05, h3 = rep(0.05, length(genes)), outcome = NULL,
                        gwas_version = "pan") {
  label <- function(h4, h3) ifelse(h4 >= 0.8, "shared", ifelse(h3 >= 0.8, "two signals", ifelse(h4 >= 0.5, "leaning", "undecided")))
  data.table(locus = locus, eqtl_dataset = dataset, gwas_version = gwas_version, gene_id = genes, symbol = genes,
             symbol_from = "fuma_magma", fuma_symbol = genes, fuma_symbol_shared_by_n_gene_ids = 1L, n_snps = 100L,
             PP.H3_two_signals = h3, PP.H4_shared = h4, PP.H4_shared_strict_prior = h4_strict,
             outcome = if (is.null(outcome)) label(h4, h3) else outcome, outcome_strict_prior = label(h4_strict, h3), top_snp = "rs1")
}
test("the reporting rule counts independent resources, only testable tests, and flags conflicts", {
  dir <- file.path(tmp, "rule"); dir.create(dir)
  put <- function(dataset, ...) fwrite(result_rows("R", dataset, c("A", "B", "C", "D"), ...), file.path(dir, sprintf("R__%s__250kb.csv", dataset)))
  #                             A: two GTEx tissues   B: two resources   C: sensitivity sets only   D: high but not testable
  put("GTEx_v8_brain_cortex",  h4 = c(0.90, 0.90, 0.10, 0.95), h4_strict = c(0.85, 0.85, 0.05, 0.90),
      outcome = c("shared", "shared", "undecided", "not testable: no expression signal"))
  put("GTEx_v8_Artery_Tibial", h4 = c(0.95, 0.02, 0.10, 0.10), h4_strict = c(0.90, 0.01, 0.05, 0.05), h3 = c(0.01, 0.95, 0.05, 0.05))
  put("PsychENCODE",           h4 = c(0.30, 0.85, 0.10, 0.92), h4_strict = c(0.20, 0.70, 0.05, 0.85),
      outcome = c("undecided", "shared", "undecided", "not testable: migraine signal not in the tested SNPs"))
  put("CommonMind_DLPFC",      h4 = c(0.99, 0.99, 0.99, 0.10))
  put("GTEx_v7_brain_cortex",  h4 = c(0.99, 0.99, 0.99, 0.10))
  s <- summarise_locus(dir, "R", 250)
  row <- function(g) s[gene_id == g]
  stopifnot(identical(s$gene_id[1:2], c("B", "A")),
            row("A")$support == "single resource", row("A")$resources_at_threshold == "GTEx", row("A")$n_datasets_shared == 2,
            row("A")$datasets_at_threshold == "GTEx_v8_Artery_Tibial; GTEx_v8_brain_cortex",
            row("B")$support == "replicated", row("B")$support_strict_prior == "single resource",
            row("B")$datasets_two_signals == "GTEx_v8_Artery_Tibial", row("B")$shared_and_two_signals, !row("A")$shared_and_two_signals,
            row("C")$support == "none", row("C")$n_H4_ge_0.8 == 2, row("C")$max_H4_primary == 0.10,       # sensitivity sets and GTEx v7 never count
            row("D")$support == "none", row("D")$n_H4_ge_0.8 == 2, row("D")$n_primary_datasets_testable == 1,   # not testable never counts
            row("A")$n_primary_datasets_tested == 3, row("A")$n_primary_datasets_testable == 3,
            "GTEx_v7_brain_cortex" %in% SENSITIVITY_ONLY_DATASETS,
            identical(resource_of(c("GTEx_v8_Nerve_Tibial", "Bryois_Microglia", "PsychENCODE")), c("GTEx", "Bryois2022", "PsychENCODE")))
  # files from another GWAS setting are kept apart
  fwrite(result_rows("R", "PsychENCODE", "A", 0.99, gwas_version = "eur"), file.path(dir, "R__PsychENCODE__250kb__eur.csv"))
  stopifnot(nrow(summarise_locus(dir, "R", 250)) == 4, row("A")$support == "single resource",
            summarise_locus(dir, "R", 250, gwas_setting = "eur")$n_datasets_tested == 1,
            basename(summary_file(dir, "R", 250, "eur")) == "R__summary_250kb__eur.csv",
            basename(summary_file(dir, "R", 250)) == "R__summary_250kb.csv")
  fwrite(fread(file.path(dir, "R__PsychENCODE__250kb.csv"))[, outcome := NULL], file.path(dir, "R__PsychENCODE__250kb.csv"))
  expect_error(summarise_locus(dir, "R", 250), "before outcomes were recorded")
})
test("locus_context reports signals and the European tier", {
  a <- locus_context(paths, "16", positions[100]); b <- locus_context(paths, 16, 5000000); d <- locus_context(paths, "2", 7000000)
  stopifnot(a$one_signal_assumption == "ok: one independent significant SNP", a$european_signal == "suggestive in Europeans (5e-8 to 1e-5)",
            b$one_signal_assumption == "check: one lead SNP but 4 independent significant SNPs",
            b$european_signal == "genome-wide significant in Europeans",
            d$one_signal_assumption == "doubtful: 2 lead SNPs", d$european_signal == "not testable with European eQTL data",
            a$fuma_locus_reach_from_lead_kb == 300, a$window_covers_fuma_locus == "not checked",
            isTRUE(locus_context(paths, "16", positions[100], 500)$window_covers_fuma_locus),
            isFALSE(locus_context(paths, "16", positions[100], 250)$window_covers_fuma_locus),
            PRIMARY_WINDOW_KB == 500, COMPARISON_WINDOW_KB == 250,
            locus_context(paths, "9", 1)$one_signal_assumption == "not checked",
            locus_context(list(), "16", positions[100])$fuma_lead_snps_in_locus == "not found")
})
# ---- Duplicated sequence: tiny stand-ins for the two UCSC tables ----
gene_spans_file <- file.path(tmp, "ensGene.txt.gz")
transcript <- function(gene, chrom, from, to) {
  data.table(bin = 1L, name = paste0("ENST", from), chrom = chrom, strand = "+", txStart = from, txEnd = to, cdsStart = from,
             cdsEnd = to, exonCount = 1L, exonStarts = paste0(from, ","), exonEnds = paste0(to, ","), score = 0L,
             name2 = gene, cdsStartStat = "cmpl", cdsEndStat = "cmpl", exonFrames = "0,")
}
write_gz(rbind(transcript("ENSG00000000101", "chr16", 1000L, 1500L), transcript("ENSG00000000101", "chr16", 1200L, 2000L),
               transcript("ENSG00000000102", "chr16", 5000L, 6000L), transcript("ENSG00000000103", "chr16", 9000L, 9500L),
               transcript("ENSG00000000102", "chr6_apd_hap1", 1L, 90000L)),
         gene_spans_file, sep = "\t", col.names = FALSE)
segdups_file <- file.path(tmp, "genomicSuperDups.txt.gz")
duplication <- function(on_chrom, from, to, identity) {
  x <- as.data.table(as.list(rep(0, length(UCSC_SEGDUP_COLUMNS))))
  setnames(x, UCSC_SEGDUP_COLUMNS)
  x[, chrom := as.character(chrom)]
  set(x, j = c("chrom", "chromStart", "chromEnd", "fracMatch"), value = list(on_chrom, from, to, identity))
  x
}
write_gz(rbind(duplication("chr16", 900L, 2100L, 0.99), duplication("chr16", 5200L, 5400L, 0.95),
               duplication("chr16", 5300L, 5600L, 0.97), duplication("chr6_apd_hap1", 1L, 90000L, 0.999)),
         segdups_file, sep = "\t", col.names = FALSE)
annotation_paths <- list(gene_spans = gene_spans_file, segdups = segdups_file)

test("duplication_by_gene measures the fraction of each gene inside duplicated sequence", {
  stopifnot(union_length(c(0, 5, 20), c(10, 15, 30)) == 25, union_length(numeric(), numeric()) == 0)
  a <- load_duplication_annotation(annotation_paths)
  stopifnot(nrow(a$gene_spans) == 3, a$gene_spans[gene_id == "ENSG00000000101", end - start] == 1000, nrow(a$duplications) == 3)
  d <- duplication_by_gene(c("ENSG00000000102", "ENSG00000000101", "ENSG00000000103", "ENSG00000009999", "ENSG00000000101"), a)
  stopifnot(identical(d$gene_id, c("ENSG00000000102", "ENSG00000000101", "ENSG00000000103", "ENSG00000009999")),
            identical(d$fraction_in_segmental_duplication, c(0.4, 1, 0, NA)),
            identical(d$max_duplication_identity, c(0.97, 0.99, NA, NA)),
            d$gene_span_hg19[2] == "16:1001-2000", is.na(d$gene_span_hg19[4]),
            is.null(suppressMessages(load_duplication_annotation(list()))))
  bad <- file.path(tmp, "bad_ensGene.txt.gz")
  write_gz(transcript("ENSG00000000101", "chr16", 1L, 2L)[, 1:15], bad, sep = "\t", col.names = FALSE)
  expect_error(read_gene_spans(bad), "expected 16 columns")
})
test("the primary-window summary carries the comparison-window labels and the duplication columns", {
  dir <- file.path(tmp, "windows"); dir.create(dir)
  put <- function(dataset, window, genes, h4) fwrite(result_rows("W", dataset, genes, h4), file.path(dir, sprintf("W__%s__%dkb.csv", dataset, window)))
  put("GTEx_v8_brain_cortex", 500, c("ENSG00000000101", "ENSG00000000102"), c(0.90, 0.90))
  put("PsychENCODE", 500, c("ENSG00000000101", "ENSG00000000102"), c(0.90, 0.10))
  put("GTEx_v8_brain_cortex", 250, "ENSG00000000101", 0.90)
  put("PsychENCODE", 250, "ENSG00000000101", 0.20)
  s <- summarise_locus(dir, "W", 500, annotation = load_duplication_annotation(annotation_paths), comparison_window_kb = 250)
  stopifnot(identical(names(s)[1:9], c("symbol", "gene_id", "support", "support_strict_prior", "support_at_250kb",
                                       "support_strict_prior_at_250kb", "gene_span_hg19",
                                       "fraction_in_segmental_duplication", "max_duplication_identity")),
            identical(s$gene_id, c("ENSG00000000101", "ENSG00000000102")),
            identical(s$support, c("replicated", "single resource")), identical(s$support_at_250kb, c("single resource", "not tested")),
            identical(s$fraction_in_segmental_duplication, c(1, 0.4)),
            !"support_at_250kb" %in% names(summarise_locus(dir, "W", 500)))
})
test("merge_intervals joins overlapping and touching regions per chromosome", {
  x <- data.table(chromosome = c("16", "16", "16", "2", "16"), start = c(100, 150, 301, 100, 1000), end = c(200, 300, 400, 200, 1100))
  stopifnot(identical(merge_intervals(x), data.table(chromosome = c("16", "16", "2"), start = c(100, 1000, 100), end = c(400, 1100, 200))),
            nrow(merge_intervals(x[0])) == 0)
})
if (have_htslib) {
  test("region files: built without duplicates, coverage recorded, and a cut window is refused", {
    region_file <- file.path(regions_dir, "Testtissue.hg38.tsv.gz")
    body <- v8[gene_id %in% c("ENSG00000000001", "ENSG00000000002")]
    span <- range(body$position)
    rows_file <- file.path(tmp, "new_rows.tsv")
    fwrite(body, rows_file, sep = "\t", col.names = FALSE)
    add_rows_to_region_file(region_file, rows_file, "16", span[1] - 1e6, span[2] + 1e6)
    add_rows_to_region_file(region_file, rows_file, "16", span[1] - 5e5, span[2] + 2e6)   # overlapping refetch
    coverage <- read_region_coverage(region_file)
    stopifnot(nrow(coverage) == 1, coverage$chromosome == "16", coverage$start == span[1] - 1e6, coverage$end == span[2] + 2e6,
              nrow(fread(cmd = paste("gzip -dc", shQuote(region_file)))) == 400, is.null(read_region_coverage(gtex_v8_file)))
    e <- read_eqtl_gtex_v8(region_file, "16", lo, hi)
    stopifnot(nrow(e$data) == 400, nrow(e$spec$coverage) == 1)
    r <- run_coloc("C", "16", positions[100], 150, "gtex_v8:Testtissue", paths = paths)
    stopifnot(r$results[gene_id == "ENSG00000000001", PP.H4_shared] > 0.9, basename(r$out_file) == "C__GTEx_v8_Testtissue__150kb.csv",
              r$settings$eqtl_region_available == sprintf("%s-%s", format_position(span[1] - 1e6), format_position(span[2] + 2e6)))
    # a region file that stops in the middle of the window must be refused, not used
    cut_file <- file.path(regions_dir, "Cuttissue.hg38.tsv.gz")
    fwrite(body[position >= span[1] + 50000], rows_file, sep = "\t", col.names = FALSE)
    add_rows_to_region_file(cut_file, rows_file, "16", span[1] + 50000, span[2] + 1e6)
    expect_error(run_coloc("C", "16", positions[100], 150, "gtex_v8:Cuttissue", paths = paths), "edge of the eQTL data available")
    expect_error(run_coloc("C", "16", positions[100], 150, "gtex_v8:Missingtissue", paths = paths), "not found")
    g <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file, "16", lo, hi)$data, "rsid_alleles"))
    expect_error(restrict_to_gwas_window(add_match_key(e$data, "rsid_alleles"), g, e$spec$region_read,
                                         data.table(chromosome = "16", start = 1, end = 1000)), "not inside one downloaded region")
  })
} else {
  skip_test("region file tests", "bgzip/tabix not found")
}

# ---- Testable or not, GWAS versions, the sample-size filter, one whole locus ----
test("outcome_label names each outcome and each reason a test is not testable", {
  stopifnot(outcome_label(0.05, 0.90, TRUE, 1, TRUE) == "shared", outcome_label(0.90, 0.05, TRUE, 1, TRUE) == "two signals",
            outcome_label(0.30, 0.60, TRUE, 1, TRUE) == "leaning", outcome_label(0.30, 0.30, TRUE, 1, TRUE) == "undecided",
            outcome_label(0.05, 0.90, FALSE, 1, TRUE) == "not testable: no expression signal",
            outcome_label(0.05, 0.90, TRUE, 0.89, TRUE) == "not testable: migraine signal not in the tested SNPs",
            outcome_label(0.05, 0.90, TRUE, 1, FALSE) == "not testable: no migraine signal in this GWAS file",
            outcome_label(0.05, 0.90, TRUE, NA, NA) == "shared",   # no reference p-value: the two GWAS checks do not gate
            EXPRESSION_SIGNAL_P == 1e-5, SIGNAL_KEPT_MIN == 0.9, GWAS_SIGNAL_P == 1e-5, SAMPLE_SIZE_TOLERANCE == 0.10)
})
test("coloc_locus labels each test: shared, two signals, and the three kinds of not testable", {
  g <- drop_ambiguous_keys(add_match_key(read_gwas_window(gwas_file, "16", lo, hi)$data, "pos_alleles"))
  e_read <- read_eqtl_gtex_v7(gtex_file, "16", lo, hi)
  e <- add_match_key(e_read$data, "pos_alleles")
  reference <- migraine_reference_p(g, positions[100])
  stopifnot(reference == min(g$gwas_p), is.na(migraine_reference_p(g, 99000000)))
  run <- function(eqtl, ref = reference) coloc_locus(g, eqtl, e_read$spec, lead_pos = positions[100], gwas_reference_p = ref)$results
  full <- run(e)
  stopifnot(full[gene_id == "ENSG00000000001", outcome] == "shared", full[gene_id == "ENSG00000000002", outcome] == "two signals",
            all(full$expression_signal), all(full$migraine_signal_kept == 1), all(full$gwas_version == "pan"),
            full[gene_id == "ENSG00000000001", outcome_strict_prior] == "shared")
  weak <- copy(e)[gene_id == "ENSG00000000001"][, p := 0.5]
  stopifnot(run(weak)$outcome == "not testable: no expression signal", !run(weak)$expression_signal)
  # a gene with its own strong expression signal elsewhere, tested without the SNPs that carry the migraine signal
  far <- e[gene_id == "ENSG00000000002" & abs(pos - positions[100]) > 30000]
  cut <- run(far)
  stopifnot(cut$outcome == "not testable: migraine signal not in the tested SNPs", cut$migraine_signal_kept < 0.9, cut$expression_signal)
  stopifnot(all(run(e, ref = 0.001)$outcome == "not testable: no migraine signal in this GWAS file"))
})
uneven <- copy(gwas_table)[c(1:15, 186:200), N := 400000]   # 30 SNPs measured in fewer people, away from the signal
uneven_file <- file.path(tmp, "gwas_uneven.tsv.gz")
write_gz(uneven, uneven_file, sep = "\t")
eur_table <- copy(gwas_table)[, N := 464672]
eur_file <- file.path(tmp, "gwas_eur.tsv.gz")
write_gz(eur_table, eur_file, sep = "\t")
test("the sample-size filter drops SNPs measured in a different number of people, and says how many", {
  g <- read_gwas_window(uneven_file, "16", lo, hi)$data
  f <- filter_sample_size(g, positions[100])
  stopifnot(f$dropped == 30, f$reference_n == 648172, nrow(f$data) == 170, positions[100] %in% f$data$pos,
            filter_sample_size(g, 28100500)$reference_n == 648172)   # lead position absent: nearest strongest SNP is used
  r <- run_coloc("F", "16", positions[100], 150, "gtex_v7", paths = modifyList(paths, list(gwas = uneven_file)), sample_size_filter = TRUE)
  stopifnot(basename(r$out_file) == "F__GTEx_v7_brain_cortex__150kb__pan_nfilter.csv", r$settings$gwas_rows_dropped_by_sample_size == 30,
            r$settings$gwas_version == "pan_nfilter", all(r$results$gwas_version == "pan_nfilter"), r$results[1, n_snps] == 170,
            run_coloc("F", "16", positions[100], 150, "gtex_v7", paths = paths)$settings$sample_size_filter == "off")
})
test("a run can use the European file with its own case counts, and is labelled", {
  two <- modifyList(paths, list(gwas = c(pan = gwas_file, eur = eur_file)))
  r <- run_coloc("E", "16", positions[100], 150, "gtex_v7", paths = two, gwas_version = "eur")
  stopifnot(basename(r$out_file) == "E__GTEx_v7_brain_cortex__150kb__eur.csv", r$settings$n_cases == 15000, r$settings$n_controls == 60000,
            r$settings$gwas_version == "eur", basename(r$settings$gwas_file) == "gwas_eur.tsv.gz", isTRUE(r$settings$lead_snp_in_eqtl_data),
            isTRUE(r$settings$migraine_signal_in_this_gwas_file), r$settings$tests_shared == 1, r$settings$tests_two_signals == 1,
            run_coloc("E", "16", positions[100], 150, "gtex_v7", paths = two)$settings$n_cases == STUDY$n_cases,
            identical(gwas_label("eur", TRUE), "eur_nfilter"), identical(gwas_suffix("pan"), ""), identical(gwas_suffix("eur"), "__eur"),
            gwas_file_for(paths, "pan") == gwas_file)
  expect_error(run_coloc("E", "16", positions[100], 150, "gtex_v7", paths = paths, gwas_version = "eur"), "no file given for GWAS version 'eur'")
  expect_error(gwas_file_for(two, "asn"), "unknown GWAS version")
})
test("locus_sources lists the fixed panel and finds catalogue slices that hold the locus", {
  catalogue_dir <- file.path(tmp, "catalogue"); dir.create(catalogue_dir)
  for (f in c("QTD000075_CommonMind_DLPFC_chr16_27.3-29.9Mb_hg38.tsv", "QTD000171_brain_cortex_chr16_27.3-29.9Mb_hg38.tsv",
              "QTD000434_ROSMAP_DLPFC_chr2_27.3-29.9Mb_hg38.tsv", "README.md")) writeLines("x", file.path(catalogue_dir, f))
  with_slices <- locus_sources(list(catalogue_dir = catalogue_dir), "16", 28500000)
  stopifnot(nrow(locus_sources(list(), "16", 28500000)) == 18,
            identical(locus_sources(list(), "16", 28500000)$source[c(1, 2, 6, 8, 9, 10, 18)],
                      c("gtex_v8", "gtex_v8:Artery_Aorta", "gtex_v8:Adipose_Subcutaneous", "gtex_v8:Pituitary", "psychencode",
                        "bryois:Astrocytes", "gtex_v7")),
            nrow(with_slices) == 20, setequal(with_slices$dataset_label[19:20], c("CommonMind_DLPFC", "GTEx_v8_brain_cortex_reprocessed")),
            nrow(locus_sources(list(catalogue_dir = catalogue_dir), "16", 29500000)) == 18)   # too close to the slice's edge
})
test("run_locus runs every dataset at both windows, summarises, and reports a failed dataset without keeping its old file", {
  locus_paths <- modifyList(paths, list(results = file.path(tmp, "locus_results")))
  sources <- data.table(source = c("gtex_v7", "psychencode", "bryois:Testcells", "gtex_v8:Missingtissue"), dataset_label = NA_character_)
  dir.create(locus_paths$results)
  stale <- result_file(locus_paths$results, "Z", "GTEx_v8_Missingtissue", 500)
  writeLines("stale", stale)
  done <- suppressMessages(run_locus("Z", "16", positions[100], locus_paths, sources = sources, windows = c(500, 250)))
  primary <- fread(done$summaries[["500"]])
  stopifnot(nrow(done$failures) == 2, all(done$failures$dataset == "GTEx_v8_Missingtissue"), all(grepl("not found", done$failures$error)),
            !file.exists(stale), file.exists(done$summaries[["250"]]),
            basename(done$summaries[["500"]]) == "Z__summary_500kb.csv", "support_at_250kb" %in% names(primary),
            primary[gene_id == "ENSG00000000001", support] == "replicated", primary[gene_id == "ENSG00000000001", n_datasets_shared] == 2,
            length(list.files(locus_paths$results, pattern = "^Z__.*kb[.]csv$")) == 8)   # 3 datasets x 2 windows + 2 summaries
})

test("the Bryois reader gives the same rows whether or not it reuses the last file read", {
  rm(list = ls(READ_CACHE), envir = READ_CACHE)
  wide_first <- read_eqtl_bryois(bryois_dir, "Testcells", "16", lo, hi)$data
  narrow_cached <- read_eqtl_bryois(bryois_dir, "Testcells", "16", 28050000, 28060000)$data     # reuses the file just read
  rm(list = ls(READ_CACHE), envir = READ_CACHE)
  narrow_fresh <- read_eqtl_bryois(bryois_dir, "Testcells", "16", 28050000, 28060000)$data
  wide_cached <- read_eqtl_bryois(bryois_dir, "Testcells", "16", lo, hi)$data
  stopifnot(isTRUE(all.equal(narrow_cached, narrow_fresh)), isTRUE(all.equal(wide_first, wide_cached)), nrow(narrow_fresh) == 11,
            nrow(read_eqtl_bryois(bryois_dir, "Testcells", "16", lo, hi)$data) == 200)
})
test("locus_table labels every locus, and the all-loci tables report run and not-run loci", {
  all_paths <- modifyList(paths, list(results = file.path(tmp, "all_results")))
  dir.create(all_paths$results)
  loci <- locus_table(all_paths)
  stopifnot(identical(loci$label, c("L1_chr16_28.1Mb", "L2_chr16_5.0Mb", "L3_chr2_7.0Mb")), identical(loci$lead_snp, c("rs1100", "rsB", "rsC")),
            LOCUS_LABEL_OVERRIDES[["7"]] == "L7_named")
  named <- file.path(tmp, "loci_named.txt")
  fwrite(data.table(GenomicLocus = 7L, rsID = "rsNamed", chr = 16L, pos = 28500000L, p = 1e-8), named, sep = "\t")
  stopifnot(locus_table(list(fuma_loci = named))$label == "L7_named")
  sources <- data.table(source = c("gtex_v7", "psychencode", "bryois:Testcells"), dataset_label = NA_character_)
  suppressMessages(run_locus(loci$label[1], loci$chr[1], loci$lead_pos[1], all_paths, sources = sources))
  overview <- overview_all_loci(all_paths, loci)
  first <- overview[locus == 1]
  genes <- supported_genes_all_loci(all_paths, loci)
  stopifnot(identical(overview$status, c("run", "not run", "not run")), first$primary_datasets_run == 2, first$primary_datasets_expected == 17,
            first$n_replicated == 1, first$genes_replicated == "GENEA", first$gene_tests == 2, first$testable_tests == 2, first$shared_tests == 2,
            first$european_signal == "suggestive in Europeans (5e-8 to 1e-5)", first$top_gene == "GENEA", first$top_gene_H4 > 0.9,
            nrow(genes) == 1, genes$symbol == "GENEA", genes$locus == 1, genes$support == "replicated",
            nrow(supported_genes_all_loci(all_paths, loci, gwas_setting = "eur")) == 0,
            all(overview_all_loci(all_paths, loci, gwas_setting = "eur")$status == "not run"))
})

if (have_htslib) {
  test("a downloaded window with no gene in reach is an empty result, and an undownloaded one is still an error", {
    desert_file <- file.path(regions_dir, "Deserttissue.hg38.tsv.gz")
    rows_file <- file.path(tmp, "desert_rows.tsv")
    far <- v8[gene_id == "ENSG00000000007"]                       # a gene whose SNPs all lie outside the GWAS window
    fwrite(far, rows_file, sep = "\t", col.names = FALSE)
    add_rows_to_region_file(desert_file, rows_file, "16", min(positions) + build_shift - 1e6, max(far$position) + 1e6)
    expect_error(run_coloc("D", "16", positions[100], 150, "gtex_v8:Deserttissue", paths = paths), "no SNP in the GWAS window")
    window <- locus_window(positions[100], 150)
    fwrite(data.table(chromosome = "16", hg19_lo = window$lo, hg19_hi = window$hi), sub("[.]tsv[.]gz$", ".fetched_windows.tsv", desert_file), sep = "\t")
    r <- run_coloc("D", "16", positions[100], 150, "gtex_v8:Deserttissue", paths = paths)
    stopifnot(nrow(r$results) == 0, nrow(r$skipped) == 0, r$settings$genes_in_eqtl_data == 0, r$settings$genes_with_no_snp_in_gwas_window == 1,
              grepl("^none", r$settings$eqtl_span_of_matched_snps), is.null(read_result_file(r$out_file)),
              isTRUE(read_eqtl_gtex_v8(desert_file, "16", window$lo, window$hi)$spec$window_downloaded),
              !isTRUE(read_eqtl_gtex_v8(desert_file, "16", window$lo - 1, window$hi)$spec$window_downloaded),   # a wider window was not downloaded
              is.null(read_fetched_windows(gtex_v8_file)))
  })
  test("a locus where no dataset can test a gene still gets a summary and an overview row", {
    empty_paths <- modifyList(paths, list(results = file.path(tmp, "empty_results")))
    dir.create(empty_paths$results)
    loci <- locus_table(empty_paths)[1]
    done <- suppressMessages(run_locus(loci$label, loci$chr, loci$lead_pos, empty_paths,
                                       sources = data.table(source = "gtex_v8:Deserttissue", dataset_label = NA_character_), windows = c(150, 100)))
    summary <- fread(done$summaries[["150"]])
    overview <- overview_all_loci(empty_paths, loci, window_kb = 150)
    stopifnot(nrow(done$failures) == 0, nrow(summary) == 0, all(c("symbol", "gene_id", "support") %in% names(summary)),
              overview$status == "run", overview$genes_tested == 0, overview$gene_tests == 0, overview$n_replicated == 0, overview$top_gene == "",
              nrow(supported_genes_all_loci(empty_paths, loci, window_kb = 150)) == 0)
  })
} else {
  skip_test("empty-result tests", "bgzip/tabix not found")
}

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d passed, %d failed, %d skipped\n", passed, length(failed), skipped_tests))
if (length(failed)) {
  cat(paste0("  ", failed, collapse = "\n"), "\n")
  quit(status = 1)
}
