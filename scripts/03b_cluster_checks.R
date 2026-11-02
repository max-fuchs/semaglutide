source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

out <- step_dir("03b_cluster_checks")

seu <- readRDS(rds_clustered)
Idents(seu) <- "seurat_clusters"

prop_sample <- seu@meta.data %>%
  dplyr::count(sample, seurat_clusters) %>%
  group_by(sample) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = sample, values_from = pct, values_fill = 0) %>%
  arrange(seurat_clusters)

print(prop_sample, n = Inf, width = Inf)
write_csv(prop_sample, file.path(out, "cluster_pct_of_sample.csv"))

p_prop <- seu@meta.data %>%
  dplyr::count(sample, patient, treatment, seurat_clusters) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(treatment, pct, fill = treatment)) +
  geom_col() +
  scale_fill_manual(values = col_treatment) +
  facet_grid(seurat_clusters ~ patient, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = "Nuclei of sample (%)", fill = "Treatment",
       title = "Cluster proportions per sample, Ctr vs. Semaglutide per patient")

save_plot(p_prop, file.path(out, "cluster_pct_of_sample.png"), width = 9, height = 16)

unclear <- c("4", "5", "6", "15", "18")

mk_unclear <- map_dfr(unclear, function(cl) {
  FindMarkers(seu, ident.1 = cl, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25,
              max.cells.per.ident = 1000, random.seed = seed) %>%
    rownames_to_column("gene") %>%
    slice_head(n = 20) %>%
    mutate(cluster = cl, .before = 1)
}) %>%
  as_tibble()

print(mk_unclear %>% dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2), n = Inf)
write_csv(mk_unclear, file.path(out, "markers_unclear_clusters.csv"))

cm_genes <- intersect(c("MYH7", "TNNT2", "MB", "MYL2", "NPPA", "TTN"), rownames(seu))

ambient_cm <- FetchData(seu, vars = c("sample", "seurat_clusters", cm_genes), layer = "counts") %>%
  filter(seurat_clusters %in% c("0", "1", "2", "9")) %>%
  group_by(seurat_clusters, sample) %>%
  summarise(across(all_of(cm_genes), ~ round(100 * mean(.x > 0), 1)), .groups = "drop")

print(ambient_cm, n = Inf, width = Inf)
write_csv(ambient_cm, file.path(out, "cm_transcripts_in_non_cm_clusters.csv"))

p_ambient <- ambient_cm %>%
  pivot_longer(all_of(cm_genes), names_to = "gene", values_to = "pct") %>%
  ggplot(aes(sample, pct, fill = sample)) +
  geom_col() +
  scale_fill_manual(values = col_sample) +
  facet_grid(gene ~ seurat_clusters, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Nuclei with >= 1 UMI (%)",
       title = "Cardiomyocyte transcripts detected in non-cardiomyocyte clusters (0 fibroblasts, 1 pericytes, 2 endothelial, 9 macrophages)")

save_plot(p_ambient, file.path(out, "cm_transcripts_in_non_cm_clusters.png"), width = 14, height = 12)
