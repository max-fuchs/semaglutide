source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

out <- step_dir("04_clean_recluster")

remove_clusters <- c(doublet = "11", low_quality_cm_fragments = "4", low_quality = "6")

n_pcs       <- 30
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
res_use     <- 0.4

seu <- readRDS(rds_soupx)

seu$cluster_initial <- seu$seurat_clusters
seu@meta.data <- seu@meta.data %>% select(-starts_with("RNA_snn_res."), -seurat_clusters)

seu$removal_reason <- case_when(
  seu$cluster_initial == remove_clusters[["doublet"]]                  ~ "doublet_cluster",
  seu$scDblFinder_class == "doublet"                                   ~ "scDblFinder_doublet",
  seu$cluster_initial == remove_clusters[["low_quality_cm_fragments"]] ~ "low_quality_cm_fragments",
  seu$cluster_initial == remove_clusters[["low_quality"]]              ~ "low_quality",
  TRUE                                                                 ~ "kept"
)

clean_steps <- seu@meta.data %>%
  dplyr::count(sample, removal_reason) %>%
  pivot_wider(names_from = removal_reason, values_from = n, values_fill = 0) %>%
  mutate(
    n_before    = rowSums(across(where(is.numeric))),
    pct_removed = round(100 * (1 - kept / n_before), 1)
  ) %>%
  relocate(sample, n_before) %>%
  arrange(factor(sample, levels = lvl_sample))

print(clean_steps, width = Inf)
write_csv(clean_steps, file.path(out, "clean_steps.csv"))

seu <- subset(seu, subset = removal_reason == "kept")
seu$removal_reason <- NULL

DefaultAssay(seu) <- "SoupX"
seu <- DietSeurat(seu, assays = c("RNA", "SoupX"), layers = c("counts", "data"))
seu <- NormalizeData(seu, assay = "SoupX")

s_genes   <- intersect(cc.genes.updated.2019$s.genes,   rownames(seu))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(seu))
seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

save_plot(ElbowPlot(seu, ndims = 50) + ggtitle("PCA elbow plot after cleaning (SoupX)"),
          file.path(out, "pca_elbow.png"), width = 7, height = 5)

loadings_tbl <- map_dfr(1:10, function(i) {
  l <- Loadings(seu, "pca")[, i]
  tibble(
    PC  = paste0("PC_", i),
    pos = paste(names(sort(l, decreasing = TRUE))[1:15], collapse = ", "),
    neg = paste(names(sort(l))[1:15], collapse = ", ")
  )
})
print(loadings_tbl, n = Inf, width = Inf)
write_csv(loadings_tbl, file.path(out, "pca_loadings_top15.csv"))

pc_cor <- cor(
  Embeddings(seu, "pca")[, 1:20],
  seu@meta.data %>%
    transmute(
      log10_nCount   = log10(nCount_RNA),
      log10_nFeature = log10(nFeature_RNA),
      complexity, percent.mt, percent.ribo, percent.malat1,
      soupx_removed_frac,
      S.Score, G2M.Score,
      scDblFinder_score
    ) %>%
    as.matrix()
)

p_pc_cor <- as.data.frame(pc_cor) %>%
  rownames_to_column("PC") %>%
  pivot_longer(-PC, names_to = "variable", values_to = "r") %>%
  mutate(PC = factor(PC, levels = paste0("PC_", 1:20))) %>%
  ggplot(aes(PC, variable, fill = r)) +
  geom_tile() +
  geom_text(aes(label = round(r, 2)), size = 2.5) +
  scale_fill_gradient2(limits = c(-1, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = "Pearson r", title = "Correlation of PCs with QC variables after cleaning")

save_plot(p_pc_cor, file.path(out, "pca_cor_qc_vars.png"), width = 12, height = 5)

seu <- FindNeighbors(seu, dims = 1:n_pcs)
seu <- FindClusters(seu, resolution = resolutions)
seu <- RunUMAP(seu, dims = 1:n_pcs)

res_cols <- paste0("SoupX_snn_res.", resolutions)

p_res <- wrap_plots(lapply(res_cols, function(r) {
  DimPlot(seu, group.by = r, label = TRUE, label.size = 3) + NoLegend() +
    ggtitle(paste("Resolution", sub("SoupX_snn_res.", "", r, fixed = TRUE)))
}), ncol = 3)
save_plot(p_res, file.path(out, "umap_resolutions.png"), width = 15, height = 10)

Idents(seu) <- paste0("SoupX_snn_res.", res_use)
seu$seurat_clusters <- Idents(seu)

p_meta <- (DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + NoLegend() +
             ggtitle(paste("Clusters after cleaning, resolution", res_use))) +
  (DimPlot(seu, group.by = "sample", cols = col_sample, shuffle = TRUE) + ggtitle("Sample")) +
  (DimPlot(seu, group.by = "patient", cols = col_patient, shuffle = TRUE) + ggtitle("Patient")) +
  (DimPlot(seu, group.by = "treatment", cols = col_treatment, shuffle = TRUE) + ggtitle("Treatment")) +
  plot_layout(ncol = 2)
save_plot(p_meta, file.path(out, "umap_meta.png"), width = 13, height = 11)

p_split <- DimPlot(seu, group.by = "seurat_clusters", split.by = "patient", label = TRUE, label.size = 2.5, ncol = 4) +
  NoLegend() + ggtitle("Clusters after cleaning split by patient")
save_plot(p_split, file.path(out, "umap_split_patient.png"), width = 18, height = 5)

p_initial <- DimPlot(seu, group.by = "cluster_initial", label = TRUE, label.size = 3) + NoLegend() +
  ggtitle("Initial clusters (03_cluster, resolution 0.4) on new UMAP")
save_plot(p_initial, file.path(out, "umap_initial_clusters.png"), width = 7, height = 6)

p_qc_umap <- FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", "complexity", "percent.mt",
                                           "percent.malat1", "soupx_removed_frac", "scDblFinder_score", "G2M.Score"),
                         ncol = 4, order = TRUE)
save_plot(p_qc_umap, file.path(out, "umap_qc_vars.png"), width = 18, height = 9)

initial_vs_new <- seu@meta.data %>%
  dplyr::count(cluster_initial, seurat_clusters) %>%
  group_by(seurat_clusters) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_map <- ggplot(initial_vs_new, aes(seurat_clusters, cluster_initial, fill = pct)) +
  geom_tile() +
  geom_text(aes(label = ifelse(pct >= 5, round(pct), "")), size = 2.5) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "New cluster", y = "Initial cluster", fill = "% of new cluster",
       title = "Origin of new clusters from initial clusters")
save_plot(p_map, file.path(out, "initial_vs_new_clusters.png"), width = 10, height = 7)

cluster_tbl <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n                 = n(),
    median_nCount     = median(nCount_RNA),
    median_nFeat      = median(nFeature_RNA),
    median_complexity = round(median(complexity), 3),
    median_mt         = round(median(percent.mt), 2),
    median_malat1     = round(median(percent.malat1), 1),
    median_soupx_rm   = round(100 * median(soupx_removed_frac), 1),
    median_dbl_score  = round(median(scDblFinder_score), 3),
    median_G2M        = round(median(G2M.Score), 3),
    .groups = "drop"
  )

cluster_sample <- seu@meta.data %>%
  dplyr::count(seurat_clusters, sample) %>%
  group_by(seurat_clusters) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = sample, values_from = pct, values_fill = 0)

cluster_tbl <- left_join(cluster_tbl, cluster_sample, by = "seurat_clusters")
print(cluster_tbl, n = Inf, width = Inf)
write_csv(cluster_tbl, file.path(out, "cluster_summary.csv"))

prop_sample <- seu@meta.data %>%
  dplyr::count(sample, seurat_clusters) %>%
  group_by(sample) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = sample, values_from = pct, values_fill = 0) %>%
  arrange(seurat_clusters)
print(prop_sample, n = Inf, width = Inf)
write_csv(prop_sample, file.path(out, "cluster_pct_of_sample.csv"))

p_cl_qc <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol = 1, log = TRUE) /
  VlnPlot(seu, features = "percent.mt", pt.size = 0)
save_plot(p_cl_qc, file.path(out, "vln_cluster_qc.png"), width = 12, height = 10)

p_dot <- DotPlot(seu, features = intersect(markers_overview, rownames(seu)), cluster.idents = TRUE) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Overview markers per cluster after cleaning (SoupX)")
save_plot(p_dot, file.path(out, "dotplot_markers_overview.png"), width = 14, height = 8)

saveRDS(seu, rds_clean)

sessionInfo()
