# Study settings for the colocalization pipeline: a template.
# Copy this file to code/study_config.R (which git does not track) and fill in your study's values.
# coloc_lib.R reads code/study_config.R when it is sourced.

# Case and control counts of the primary GWAS.
STUDY <- list(n_cases = 10000, n_controls = 50000)

# GWAS versions a run can use: the file name under data/fuma-input/ and its own counts.
# "pan" is the file the loci were defined in; "eur" is a European-only version of it.
GWAS_VERSIONS <- list(
  pan = list(file = "GWAS_pan.tsv.gz", n_cases = STUDY$n_cases, n_controls = STUDY$n_controls),
  eur = list(file = "GWAS_eur.tsv.gz", n_cases = 7000, n_controls = 35000))

# Optional fixed file-name labels, named by FUMA locus number, for example c(`3` = "L3_my_locus").
LOCUS_LABEL_OVERRIDES <- character()

# Study-specific file and folder names.
STUDY_FILES <- list(
  fuma_job_dir = "FUMA_job_full",                    # under results/fuma-downloads/: GenomicRiskLoci.txt, magma.genes.out
  gene_evidence = "gene-evidence.csv",               # under results/: table that evidence_add_coloc.R fills
  ancestry_table = "per-ancestry-signal-by-locus.csv",   # under results/: lead-SNP p per ancestry group
  gtex_v7_slice = "GTEx_v7_Brain_Cortex.loci_1Mb_hg19.tsv",    # under data/eqtl-full/
  psychencode_slice = "PsychENCODE_Full.loci_1Mb_hg19.txt")    # under data/eqtl-full/
