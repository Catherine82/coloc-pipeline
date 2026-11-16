Colocalization pipeline: GWAS loci against eQTL datasets (coloc.abf)
R 4.6.1 with the packages coloc 5.2.3 and data.table. tabix and bgzip (htslib) for the GTEx v8 files.

What it does
Tests, one gene at a time, whether a GWAS locus and a gene's eQTL share one causal variant
(coloc.abf). Loci come from a FUMA SNP2GENE job. Each locus is run against a fixed panel of eQTL
datasets: GTEx v8 brain cortex, arteries, tibial nerve, adipose and pituitary, PsychENCODE, and
eight brain cell types (Bryois et al. 2022).

Files
code/coloc_lib.R              All functions and all analysis rules. The rules are named constants
                              near the top of the file: the two priors, the reporting rule, what
                              makes a gene test "testable", the outcome labels, the tissue panel
                              and the windows.
code/study_config.example.R   Template for the study settings (sample sizes, file names). Copy it
                              to code/study_config.R, which is not tracked, and fill it in.
code/run_coloc.R              One locus against one eQTL dataset.
code/run_locus.R              One locus against every dataset at both windows, then the summary table.
code/run_all_loci.R           Every locus of the FUMA job, then the all-loci tables.
code/fetch_gtex_v8_regions.R  Downloads the region around a locus from the official GTEx v8 tissue
                              files. Only genome coordinates are sent to the server.
code/summarise_coloc.R        Rebuilds one locus summary from its result files.
code/evidence_add_coloc.R     Adds the colocalization columns to a gene evidence table.
code/tests/test_coloc_lib.R   39 tests on small made-up files; they need no study settings.
                              Run from this folder:  Rscript code/tests/test_coloc_lib.R

Not included
GWAS summary statistics, eQTL files, results, and the study settings. Where the code looks for
the data is set in project_paths() in coloc_lib.R.
