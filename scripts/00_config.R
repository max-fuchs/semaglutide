library(tibble)
library(dplyr)

prefix <- "hs_sema"
seed <- 42

proj_dir <- "/project/ag-baer/Semaglutide"
raw_dir <- file.path(proj_dir, "data", "raw")
cr_dir <- file.path(raw_dir, "2026_040595_2", "N8773-N8774", "02_DataProcessing", "02_cellranger_primary_analysis")
rds_dir <- file.path(proj_dir, "data", "rds")
out_dir <- file.path(proj_dir, "output")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

invisible(lapply(c(raw_dir, rds_dir, fig_dir, tab_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

step_dir <- function(step) {
  d <- file.path(out_dir, step)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

save_plot <- function(p, file, width = 8, height = 6, dpi = 150, show = TRUE) {
  ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
  if (show) print(p)
  invisible(p)
}

samples <- tribble(
  ~sample_id, ~run,  ~ocm,  ~patient, ~donor,     ~treatment,
  "S01s1",    "S01", "OB1", "P1",     "H.25.50",  "Ctr",
  "S01s2",    "S01", "OB2", "P2",     "H.25.121", "Ctr",
  "S01s3",    "S01", "OB3", "P4",     "H.25.167", "Ctr",
  "S01s4",    "S01", "OB4", "P5",     "H.25.39",  "Ctr",
  "S02s1",    "S02", "OB1", "P1",     "H.25.50",  "Semaglutide",
  "S02s2",    "S02", "OB2", "P2",     "H.25.121", "Semaglutide",
  "S02s3",    "S02", "OB3", "P4",     "H.25.167", "Semaglutide",
  "S02s4",    "S02", "OB4", "P5",     "H.25.39",  "Semaglutide"
) |>
  mutate(
    sample = paste(patient, treatment, sep = "_"),
    dir = file.path(cr_dir, run, "per_sample_outs", sample_id),
    mtx_dir = file.path(dir, "sample_filtered_feature_bc_matrix"),
    mtx_raw = file.path(dir, "sample_raw_feature_bc_matrix"),
    metrics = file.path(dir, "metrics_summary.csv")
  )

lvl_treatment <- c("Ctr", "Semaglutide")
lvl_patient <- c("P1", "P2", "P4", "P5")
lvl_run <- c("S01", "S02")
lvl_sample <- as.vector(t(outer(lvl_patient, lvl_treatment, paste, sep = "_")))

col_treatment <- c(Ctr = "#7F7F7F", Semaglutide = "#D55E00")
col_patient <- c(P1 = "#0072B2", P2 = "#009E73", P4 = "#CC79A7", P5 = "#E69F00")
col_sample <- setNames(
  unlist(lapply(lvl_patient, function(p) c(colorspace::lighten(col_patient[[p]], 0.45), col_patient[[p]]))),
  lvl_sample
)

markers_overview <- c(
  "TTN", "RYR2", "MYH7", "MYH6", "TNNT2", "NPPA", "NPPB",
  "DCN", "GSN", "PDGFRA", "COL1A1",
  "VWF", "PECAM1", "CDH5", "NPR3", "PROX1",
  "RGS5", "ABCC9", "KCNJ8",
  "MYH11", "ACTA2",
  "CD163", "MRC1", "F13A1",
  "PTPRC", "CD3E", "NKG7",
  "KIT", "CPA3",
  "ADIPOQ", "PLIN1",
  "NRXN1", "NRXN3",
  "MSLN", "MKI67"
)

rds_raw <- file.path(rds_dir, paste0(prefix, "_raw.rds"))
rds_qc <- file.path(rds_dir, paste0(prefix, "_qc.rds"))
rds_clustered <- file.path(rds_dir, paste0(prefix, "_clustered.rds"))
rds_soupx <- file.path(rds_dir, paste0(prefix, "_soupx.rds"))
rds_clean <- file.path(rds_dir, paste0(prefix, "_clean.rds"))
rds_annotated <- file.path(rds_dir, paste0(prefix, "_annotated.rds"))