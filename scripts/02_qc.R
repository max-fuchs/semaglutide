source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(scDblFinder)
})

out <- step_dir("02_qc")

seu <- readRDS(rds_raw)

seu$log10_nCount   <- log10(seu$nCount_RNA)
seu$log10_nFeature <- log10(seu$nFeature_RNA)
seu$complexity     <- seu$log10_nFeature / seu$log10_nCount

qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "complexity")

violin_vars <- c(log10_nCount = "log10 UMI counts", log10_nFeature = "log10 genes detected",
                 percent.mt = "Mitochondrial reads (%)", percent.ribo = "Ribosomal reads (%)",
                 complexity = "Complexity (log10 genes / log10 UMI)")

qc_quantiles <- function(md) {
  md %>%
    select(sample, all_of(qc_vars)) %>%
    pivot_longer(-sample, names_to = "metric") %>%
    group_by(metric, sample) %>%
    summarise(
      min = min(value),
      q01 = quantile(value, 0.01),
      q05 = quantile(value, 0.05),
      q50 = median(value),
      q95 = quantile(value, 0.95),
      q99 = quantile(value, 0.99),
      max = max(value),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

qc_violin <- function(md, title) {
  md %>%
    select(sample, all_of(names(violin_vars))) %>%
    pivot_longer(-sample, names_to = "metric") %>%
    mutate(metric = factor(violin_vars[metric], levels = violin_vars)) %>%
    ggplot(aes(sample, value, fill = sample)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.2) +
    scale_fill_manual(values = col_sample) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = NULL, title = title)
}

density_plot <- function(md, var, label, log_x = FALSE) {
  p <- ggplot(md, aes(.data[[var]], colour = treatment)) +
    geom_density(linewidth = 0.6) +
    scale_colour_manual(values = col_treatment) +
    facet_wrap(~ patient, ncol = 4) +
    theme_bw() +
    labs(x = label, y = "Density", colour = "Treatment")
  if (log_x) p <- p + scale_x_log10()
  p
}

md_pre <- seu@meta.data %>% as_tibble(rownames = "cell")

quantiles_pre <- qc_quantiles(md_pre)
print(quantiles_pre, n = Inf, width = Inf)
write_csv(quantiles_pre, file.path(out, "qc_quantiles_prefilter.csv"))

save_plot(qc_violin(md_pre, "QC metrics per sample before filtering"),
          file.path(out, "qc_violin_prefilter.png"), width = 8, height = 14)

p_density_pre <- density_plot(md_pre, "nFeature_RNA", "Genes detected", log_x = TRUE) /
  density_plot(md_pre, "nCount_RNA", "UMI counts", log_x = TRUE) /
  density_plot(md_pre, "percent.mt", "Mitochondrial reads (%)") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "QC distributions per patient before filtering, Ctr vs. Semaglutide")

save_plot(p_density_pre, file.path(out, "qc_density_prefilter.png"), width = 14, height = 10)

p_scatter_count <- ggplot(md_pre, aes(nCount_RNA, nFeature_RNA, colour = percent.mt)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_c(limits = c(0, 5), oob = scales::squish) +
  facet_wrap(~ sample, ncol = 4) +
  theme_bw() +
  labs(x = "UMI counts", y = "Genes detected", colour = "MT (%)",
       title = "UMI counts vs. genes detected before filtering")

save_plot(p_scatter_count, file.path(out, "qc_scatter_counts_features_prefilter.png"), width = 16, height = 8)

p_scatter_mt <- ggplot(md_pre, aes(nFeature_RNA, percent.mt)) +
  geom_point(size = 0.2, alpha = 0.4, colour = "grey30") +
  scale_x_log10() +
  facet_wrap(~ sample, ncol = 4) +
  theme_bw() +
  labs(x = "Genes detected", y = "Mitochondrial reads (%)",
       title = "Genes detected vs. mitochondrial reads before filtering")

save_plot(p_scatter_mt, file.path(out, "qc_scatter_features_mt_prefilter.png"), width = 16, height = 8)

qc_thr <- list(min_features = 500, min_counts = 1000, max_mt = 10, min_complexity = 0.8)

seu$qc_pass <- seu$nFeature_RNA >= qc_thr$min_features &
  seu$nCount_RNA >= qc_thr$min_counts &
  seu$percent.mt < qc_thr$max_mt &
  seu$complexity >= qc_thr$min_complexity

md_status <- seu@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(qc_status = factor(if_else(qc_pass, "pass", "fail"), levels = c("pass", "fail")))

p_scatter_status <- ggplot(md_status, aes(nCount_RNA, nFeature_RNA, colour = qc_status)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c(pass = "grey60", fail = "#D55E00")) +
  facet_wrap(~ sample, ncol = 4) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "UMI counts", y = "Genes detected", colour = "QC",
       title = sprintf("Nuclei removed by QC (genes >= %d, UMI >= %d, MT < %d%%, complexity >= %.2f)",
                       qc_thr$min_features, qc_thr$min_counts, qc_thr$max_mt, qc_thr$min_complexity))

save_plot(p_scatter_status, file.path(out, "qc_scatter_status.png"), width = 16, height = 8)

qc_steps <- md_status %>%
  group_by(sample) %>%
  summarise(
    n_raw           = n(),
    fail_features   = sum(nFeature_RNA < qc_thr$min_features),
    fail_counts     = sum(nCount_RNA < qc_thr$min_counts),
    fail_mt         = sum(percent.mt >= qc_thr$max_mt),
    fail_complexity = sum(complexity < qc_thr$min_complexity),
    n_pass          = sum(qc_pass),
    .groups = "drop"
  ) %>%
  mutate(pct_removed_qc = round(100 * (1 - n_pass / n_raw), 1))

seu <- subset(seu, subset = qc_pass)

set.seed(seed)
sce_dbl <- scDblFinder(LayerData(seu, assay = "RNA", layer = "counts"),
                       samples = as.character(seu$sample))

stopifnot(identical(colnames(sce_dbl), colnames(seu)))

seu$scDblFinder_class <- factor(sce_dbl$scDblFinder.class, levels = c("singlet", "doublet"))
seu$scDblFinder_score <- sce_dbl$scDblFinder.score
rm(sce_dbl)
gc()

md_post <- seu@meta.data %>% as_tibble(rownames = "cell")

qc_steps <- qc_steps %>%
  left_join(
    md_post %>%
      group_by(sample) %>%
      summarise(n_doublet = sum(scDblFinder_class == "doublet"), .groups = "drop"),
    by = "sample"
  ) %>%
  mutate(
    pct_doublet       = round(100 * n_doublet / n_pass, 1),
    n_singlet         = n_pass - n_doublet,
    pct_removed_total = round(100 * (1 - n_singlet / n_raw), 1)
  )

print(qc_steps, width = Inf)
write_csv(qc_steps, file.path(out, "qc_steps.csv"))

quantiles_post <- qc_quantiles(md_post)
write_csv(quantiles_post, file.path(out, "qc_quantiles_postfilter.csv"))

save_plot(qc_violin(md_post, "QC metrics per sample after filtering"),
          file.path(out, "qc_violin_postfilter.png"), width = 8, height = 14)

p_density_post <- density_plot(md_post, "nFeature_RNA", "Genes detected", log_x = TRUE) /
  density_plot(md_post, "nCount_RNA", "UMI counts", log_x = TRUE) /
  density_plot(md_post, "percent.mt", "Mitochondrial reads (%)") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "QC distributions per patient after filtering, Ctr vs. Semaglutide")

save_plot(p_density_post, file.path(out, "qc_density_postfilter.png"), width = 14, height = 10)

p_dbl_scatter <- md_post %>%
  arrange(scDblFinder_class) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour = scDblFinder_class)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c(singlet = "grey60", doublet = "#CC0000")) +
  facet_wrap(~ sample, ncol = 4) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "UMI counts", y = "Genes detected", colour = "scDblFinder",
       title = "Predicted doublets after QC filtering (flagged, not removed)")

save_plot(p_dbl_scatter, file.path(out, "qc_doublets_scatter.png"), width = 16, height = 8)

p_dbl_bar <- ggplot(qc_steps, aes(sample, pct_doublet, fill = sample)) +
  geom_col() +
  geom_text(aes(label = n_doublet), vjust = -0.4, size = 3) +
  scale_fill_manual(values = col_sample) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Predicted doublets (%)", title = "Predicted doublet rate per sample")

save_plot(p_dbl_bar, file.path(out, "qc_doublets_rate.png"), width = 7, height = 5)

saveRDS(seu, rds_qc)

sessionInfo()