source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

out <- step_dir("03_cluster")

n_pcs       <- 30
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
res_use     <- 0.4

seu <- readRDS(rds_qc)

seu[["percent.malat1"]] <- PercentageFeatureSet(seu, features = "MALAT1")

seu <- NormalizeData(seu)

s_genes   <- intersect(cc.genes.updated.2019$s.genes,   rownames(seu))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(seu))
seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

save_plot(ElbowPlot(seu, ndims = 50) + ggtitle("PCA elbow plot"),
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
  labs(x = NULL, y = NULL, fill = "Pearson r", title = "Correlation of PCs with QC variables")

save_plot(p_pc_cor, file.path(out, "pca_cor_qc_vars.png"), width = 12, height = 4.5)

p_pca <- (DimPlot(seu, reduction = "pca", group.by = "sample", cols = col_sample, shuffle = TRUE, raster = FALSE) +
            ggtitle("PCA by sample")) |
  (DimPlot(seu, reduction = "pca", group.by = "treatment", cols = col_treatment, shuffle = TRUE, raster = FALSE) +
     ggtitle("PCA by treatment"))
save_plot(p_pca, file.path(out, "pca_sample_treatment.png"), width = 12, height = 5)

seu <- FindNeighbors(seu, dims = 1:n_pcs)
seu <- FindClusters(seu, resolution = resolutions)
seu <- RunUMAP(seu, dims = 1:n_pcs)

res_cols <- paste0("RNA_snn_res.", resolutions)

p_res <- wrap_plots(lapply(res_cols, function(r) {
  DimPlot(seu, group.by = r, label = TRUE, label.size = 3) + NoLegend() +
    ggtitle(paste("Resolution", sub("RNA_snn_res.", "", r, fixed = TRUE)))
}), ncol = 3)
save_plot(p_res, file.path(out, "umap_resolutions.png"), width = 15, height = 10)

Idents(seu) <- paste0("RNA_snn_res.", res_use)
seu$seurat_clusters <- Idents(seu)

p_meta <- (DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + NoLegend() +
             ggtitle(paste("Clusters, resolution", res_use))) +
  (DimPlot(seu, group.by = "sample", cols = col_sample, shuffle = TRUE) + ggtitle("Sample")) +
  (DimPlot(seu, group.by = "patient", cols = col_patient, shuffle = TRUE) + ggtitle("Patient")) +
  (DimPlot(seu, group.by = "treatment", cols = col_treatment, shuffle = TRUE) + ggtitle("Treatment")) +
  plot_layout(ncol = 2)
save_plot(p_meta, file.path(out, "umap_meta.png"), width = 13, height = 11)

p_split <- DimPlot(seu, group.by = "seurat_clusters", split.by = "patient", label = TRUE, label.size = 2.5, ncol = 4) +
  NoLegend() + ggtitle("Clusters split by patient")
save_plot(p_split, file.path(out, "umap_split_patient.png"), width = 18, height = 5)

p_qc_umap <- FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", "complexity", "percent.mt",
                                           "percent.malat1", "percent.ribo", "S.Score", "G2M.Score"),
                         ncol = 4, order = TRUE)
save_plot(p_qc_umap, file.path(out, "umap_qc_vars.png"), width = 18, height = 9)

p_dbl_umap <- (DimPlot(seu, group.by = "scDblFinder_class", order = "doublet",
                       cols = c(singlet = "grey85", doublet = "firebrick")) +
                 ggtitle("scDblFinder class")) +
  FeaturePlot(seu, features = "scDblFinder_score", order = TRUE)
save_plot(p_dbl_umap, file.path(out, "umap_doublets.png"), width = 13, height = 5.5)

cluster_tbl <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n                 = n(),
    pct_doublet       = round(100 * mean(scDblFinder_class == "doublet"), 1),
    median_nCount     = median(nCount_RNA),
    median_nFeat      = median(nFeature_RNA),
    median_complexity = round(median(complexity), 3),
    median_mt         = round(median(percent.mt), 2),
    median_malat1     = round(median(percent.malat1), 1),
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

p_cl_qc <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol = 1, log = TRUE) /
  VlnPlot(seu, features = "percent.mt", pt.size = 0)
save_plot(p_cl_qc, file.path(out, "vln_cluster_qc.png"), width = 12, height = 10)

p_dot <- DotPlot(seu, features = intersect(markers_overview, rownames(seu)), cluster.idents = TRUE) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Overview markers per cluster")
save_plot(p_dot, file.path(out, "dotplot_markers_overview.png"), width = 14, height = 8)

sex_genes <- intersect(c("XIST", "RPS4Y1", "UTY", "DDX3Y", "KDM5D"), rownames(seu))

sex_tbl <- FetchData(seu, vars = c("sample", sex_genes), layer = "data") %>%
  group_by(sample) %>%
  summarise(across(all_of(sex_genes), ~ round(100 * mean(.x > 0), 1)), .groups = "drop")
print(sex_tbl)
write_csv(sex_tbl, file.path(out, "sex_check_pct_expressing.csv"))

p_sex <- VlnPlot(seu, features = sex_genes, group.by = "sample", cols = col_sample,
                 pt.size = 0, ncol = length(sex_genes))
save_plot(p_sex, file.path(out, "vln_sex_genes.png"), width = 3 * length(sex_genes), height = 4)

saveRDS(seu, rds_clustered)

sessionInfo()
