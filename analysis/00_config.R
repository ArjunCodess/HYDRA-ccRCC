set.seed(20260530)

PROJECT_ID <- "TCGA-KIRC"

LOCAL_R_LIB <- file.path("environment", "R-library")
if (!dir.exists(LOCAL_R_LIB)) dir.create(LOCAL_R_LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(LOCAL_R_LIB, winslash = "/", mustWork = TRUE), .libPaths()))

DIRS <- list(
  raw = file.path("data", "raw"),
  processed = file.path("data", "processed"),
  metadata = file.path("data", "metadata"),
  tables = file.path("results", "tables"),
  figures = file.path("results", "figures")
)

for (path in DIRS) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

FILES <- list(
  tcga_se = file.path(DIRS$processed, "tcga_kirc_star_counts_se.rds"),
  tcga_counts = file.path(DIRS$processed, "tcga_kirc_counts.rds"),
  tcga_coldata = file.path(DIRS$metadata, "tcga_kirc_coldata.csv"),
  tcga_clinical = file.path(DIRS$metadata, "tcga_kirc_clinical.csv"),
  tcga_vst = file.path(DIRS$processed, "tcga_kirc_vst.rds"),
  tcga_deg = file.path(DIRS$tables, "tcga_kirc_deseq2_tumor_vs_normal.csv"),
  tcga_survival = file.path(DIRS$tables, "tcga_kirc_cox_models.csv"),
  tcga_enrichment = file.path(DIRS$tables, "tcga_kirc_enrichment_go_bp.csv"),
  tcga_session = file.path("environment", "sessionInfo.txt")
)

THRESHOLDS <- list(
  deg_fdr = 0.05,
  deg_abs_log2fc = 1,
  strict_survival_fdr = 0.05,
  high_confidence_survival_fdr = 0.01,
  ph_min_p = 0.05,
  min_abs_log_hr = log(1.25),
  high_confidence_abs_log_hr = log(1.5),
  min_geo_abs_log2fc = 0.25,
  min_count = 10,
  min_samples = 10
)

RESAMPLING <- list(
  seed = 20260725,
  stability_repeats = 20,
  stability_fraction = 0.80,
  cv_repeats = 20,
  cv_folds = 5,
  pipeline_bootstrap_repeats = 100
)

SOURCE_URLS <- list(
  tcga_gdc = "https://portal.gdc.cancer.gov/projects/TCGA-KIRC",
  geo_gse40435 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40435",
  geo_gse53757 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53757",
  geo_gse29609 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29609",
  emtab1980 = "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1980",
  hpa_single_cell = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip",
  aran2015_purity_study = "https://doi.org/10.1038/ncomms9971",
  aran2015_purity_data = paste0(
    "https://static-content.springer.com/esm/art%3A10.1038%2F",
    "ncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"
  )
)

BIOC_PACKAGES <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "limma",
  "GEOquery",
  "Biobase",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "GO.db",
  "msigdbr",
  "fgsea",
  "ComplexHeatmap",
  "EnhancedVolcano"
)

CRAN_PACKAGES <- c(
  "tidyverse",
  "janitor",
  "here",
  "readxl",
  "survival",
  "survminer",
  "broom",
  "ggrepel",
  "patchwork",
  "pheatmap",
  "UpSetR"
)

TCGA_SAMPLE_TYPES <- c("Primary Tumor", "Solid Tissue Normal")
