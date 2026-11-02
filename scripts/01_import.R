source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

out <- step_dir("01_import")

stopifnot(all(dir.exists(samples$mtx_dir)))

counts <- lapply(seq_len(nrow(samples)), function(i) {
  m <- Read10X(samples$mtx_dir[i])
  if (is.list(m)) m <- m[["Gene Expression"]]
  colnames(m) <- paste(samples$sample[i], colnames(m), sep = "_")
  m
})
names(counts) <- samples$sample

stopifnot(all(sapply(counts, function(m) identical(rownames(m), rownames(counts[[1]])))))

cell_sample <- rep(samples$sample, sapply(counts, ncol))

seu <- CreateSeuratObject(counts = do.call(cbind, counts), project = prefix,
                          min.cells = 0, min.features = 0)
rm(counts)
gc()

stopifnot(ncol(seu) == length(cell_sample))

md <- samples[match(cell_sample, samples$sample), ]

seu$sample    <- factor(md$sample, levels = lvl_sample)
seu$sample_id <- md$sample_id
seu$patient   <- factor(md$patient, levels = lvl_patient)
seu$donor     <- md$donor
seu$treatment <- factor(md$treatment, levels = lvl_treatment)
seu$library   <- factor(md$run, levels = lvl_run)
seu$ocm       <- md$ocm

seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

summary_import <- seu@meta.data %>%
  group_by(sample, patient, treatment, library, ocm) %>%
  summarise(
    n_nuclei      = n(),
    median_nCount = median(nCount_RNA),
    median_nFeat  = median(nFeature_RNA),
    median_mt     = median(percent.mt),
    median_ribo   = median(percent.ribo),
    .groups = "drop"
  )

pairs_import <- summary_import %>%
  select(patient, treatment, n_nuclei) %>%
  pivot_wider(names_from = treatment, values_from = n_nuclei) %>%
  mutate(ratio_sema_ctr = round(Semaglutide / Ctr, 2))

metrics_cr <- map2_dfr(samples$metrics, samples$sample, function(f, s) {
  read_csv(f, col_types = cols(.default = "c")) %>% mutate(sample = s, .before = 1)
})

print(summary_import)
print(pairs_import)
print(grep("^MT-", rownames(seu), value = TRUE))
print(table(seu$patient, seu$treatment))

write_csv(summary_import, file.path(out, "summary_import.csv"))
write_csv(pairs_import, file.path(out, "pairs_import.csv"))
write_csv(metrics_cr, file.path(out, "cellranger_metrics_all.csv"))

saveRDS(seu, rds_raw)

sessionInfo()
