source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(tidyverse)
  library(patchwork)
})

out <- step_dir("03c_soupx")

rho_method <- "auto"

cm_est_genes   <- c("MB", "MYL2", "MYL3", "TNNT2", "MYH7", "CSRP3", "TNNC1", "CKM")
cm_check_genes <- c("MB", "TNNT2", "MYL2", "MYH7", "NPPA", "TTN")
non_cm_clusters <- c("0", "1", "2", "9")

own_markers <- list(
  "0" = c("DCN", "PDGFRA", "COL1A1"),
  "1" = c("RGS5", "ABCC9", "KCNJ8"),
  "2" = c("VWF", "PECAM1", "CDH5"),
  "9" = c("CD163", "MRC1", "F13A1"),
  "3" = c("TTN", "RYR2", "TNNT2"),
  "5" = c("TTN", "RYR2", "MYO18B")
)

seu <- readRDS(rds_clustered)
counts_rna <- LayerData(seu, assay = "RNA", layer = "counts")

soup_res <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples$sample[i]
  message("SoupX: ", s)
  
  cells <- colnames(seu)[seu$sample == s]
  toc <- counts_rna[, cells]
  colnames(toc) <- sub(paste0("^", s, "_"), "", colnames(toc))
  
  tod <- Read10X(samples$mtx_raw[i])
  if (is.list(tod)) tod <- tod[["Gene Expression"]]
  stopifnot(identical(rownames(tod), rownames(toc)), all(colnames(toc) %in% colnames(tod)))
  
  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  rm(tod)
  gc()
  
  sc <- setClusters(sc, setNames(as.character(seu$seurat_clusters[cells]), colnames(toc)))
  
  rho_auto <- tryCatch(
    mean(autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE, verbose = FALSE)$metaData$rho),
    error = function(e) NA_real_
  )
  
  cm_list <- list(CM = intersect(cm_est_genes, rownames(toc)))
  use_cm  <- estimateNonExpressingCells(sc, nonExpressedGeneList = cm_list)
  rho_cm  <- mean(calculateContaminationFraction(sc, cm_list, useToEst = use_cm,
                                                 forceAccept = TRUE, verbose = FALSE)$metaData$rho)
  
  method_used <- if (rho_method == "auto" && !is.na(rho_auto)) "auto" else "cm_markers"
  rho_used    <- if (method_used == "auto") rho_auto else rho_cm
  
  sc  <- setContaminationFraction(sc, rho_used, forceAccept = TRUE)
  adj <- adjustCounts(sc, roundToInt = TRUE, verbose = 0)
  colnames(adj) <- paste(s, colnames(adj), sep = "_")
  
  soup_top <- sc$soupProfile %>%
    rownames_to_column("gene") %>%
    arrange(desc(est)) %>%
    slice_head(n = 15) %>%
    mutate(sample = s)
  
  list(
    adj      = adj,
    rho      = tibble(sample = s, rho_auto = rho_auto, rho_cm = rho_cm,
                      rho_used = rho_used, method_used = method_used),
    soup_top = soup_top
  )
})

rho_tbl <- bind_rows(lapply(soup_res, `[[`, "rho")) %>%
  mutate(across(starts_with("rho"), ~ round(.x, 3)))
print(rho_tbl, width = Inf)
write_csv(rho_tbl, file.path(out, "soupx_rho.csv"))

soup_top <- bind_rows(lapply(soup_res, `[[`, "soup_top"))
write_csv(soup_top, file.path(out, "soupx_soup_profile_top15.csv"))

soup_mat <- do.call(cbind, lapply(soup_res, `[[`, "adj"))[, colnames(seu)]
rm(soup_res)
gc()

seu[["SoupX"]] <- CreateAssay5Object(counts = soup_mat)
seu <- NormalizeData(seu, assay = "SoupX")

seu$soupx_removed_frac <- 1 - Matrix::colSums(soup_mat) / Matrix::colSums(counts_rna[, colnames(seu)])
rm(soup_mat)
gc()

removed_tbl <- seu@meta.data %>%
  group_by(sample) %>%
  summarise(
    median_removed = round(100 * median(soupx_removed_frac), 2),
    q90_removed    = round(100 * quantile(soupx_removed_frac, 0.9), 2),
    .groups = "drop"
  )
print(removed_tbl)
write_csv(removed_tbl, file.path(out, "soupx_removed_per_sample.csv"))

p_rho <- rho_tbl %>%
  pivot_longer(c(rho_auto, rho_cm), names_to = "estimate", values_to = "rho") %>%
  mutate(sample = factor(sample, levels = lvl_sample),
         estimate = recode(estimate, rho_auto = "autoEstCont", rho_cm = "CM marker genes")) %>%
  ggplot(aes(sample, rho, fill = estimate)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(autoEstCont = "grey40", `CM marker genes` = "#0072B2")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Estimated contamination fraction (rho)", fill = "Estimate",
       title = paste("SoupX contamination estimates per sample (used:", rho_method, ")"))

save_plot(p_rho, file.path(out, "soupx_rho.png"), width = 8, height = 5)

soup_genes <- unique(soup_top$gene)

p_soup <- soup_top %>%
  mutate(sample = factor(sample, levels = lvl_sample),
         gene = factor(gene, levels = rev(soup_genes))) %>%
  ggplot(aes(sample, gene, fill = 100 * est)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Share of soup (%)",
       title = "Top 15 genes of the ambient (soup) profile per sample")

save_plot(p_soup, file.path(out, "soupx_soup_profile.png"), width = 8, height = 10)

p_removed <- ggplot(seu@meta.data, aes(seurat_clusters, 100 * soupx_removed_frac, fill = treatment)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = col_treatment) +
  theme_bw() +
  labs(x = "Cluster (resolution 0.4)", y = "Counts removed by SoupX (%)", fill = "Treatment",
       title = "Fraction of counts removed per nucleus by cluster")

save_plot(p_removed, file.path(out, "soupx_removed_per_cluster.png"), width = 12, height = 5)

expr_pct <- function(assay, genes, clusters) {
  m <- LayerData(seu, assay = assay, layer = "counts")
  genes <- intersect(genes, rownames(m))
  seu@meta.data %>%
    as_tibble(rownames = "cell") %>%
    filter(seurat_clusters %in% clusters) %>%
    group_by(seurat_clusters, sample) %>%
    group_modify(function(d, k) {
      tibble(gene = genes, pct = 100 * Matrix::rowMeans(m[genes, d$cell, drop = FALSE] > 0))
    }) %>%
    ungroup() %>%
    mutate(assay = assay)
}

cm_check <- bind_rows(
  expr_pct("RNA", cm_check_genes, non_cm_clusters),
  expr_pct("SoupX", cm_check_genes, non_cm_clusters)
) %>%
  mutate(assay = factor(assay, levels = c("RNA", "SoupX")))

cm_check_wide <- cm_check %>%
  mutate(pct = round(pct, 1)) %>%
  pivot_wider(names_from = assay, values_from = pct)
write_csv(cm_check_wide, file.path(out, "check_cm_genes_in_non_cm_clusters.csv"))

p_cm <- ggplot(cm_check, aes(sample, pct, fill = assay)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(RNA = "grey60", SoupX = "#0072B2")) +
  facet_grid(gene ~ seurat_clusters, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Nuclei with >= 1 UMI (%)", fill = "Counts",
       title = "Cardiomyocyte transcripts in non-cardiomyocyte clusters before and after SoupX")

save_plot(p_cm, file.path(out, "check_cm_genes_in_non_cm_clusters.png"), width = 14, height = 12)

own_check <- imap_dfr(own_markers, function(g, cl) {
  bind_rows(expr_pct("RNA", g, cl), expr_pct("SoupX", g, cl))
}) %>%
  mutate(assay = factor(assay, levels = c("RNA", "SoupX")))

own_check_summary <- own_check %>%
  pivot_wider(names_from = assay, values_from = pct) %>%
  group_by(seurat_clusters, gene) %>%
  summarise(
    pct_RNA   = round(mean(RNA), 1),
    pct_SoupX = round(mean(SoupX), 1),
    retained  = round(100 * mean(SoupX) / mean(RNA), 1),
    .groups = "drop"
  )
print(own_check_summary, n = Inf)
write_csv(own_check_summary, file.path(out, "check_own_markers_retained.csv"))

p_own <- own_check %>%
  mutate(label = paste0("Cluster ", seurat_clusters, ": ", gene)) %>%
  ggplot(aes(sample, pct, fill = assay)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(RNA = "grey60", SoupX = "#0072B2")) +
  facet_wrap(~ label, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Nuclei with >= 1 UMI (%)", fill = "Counts",
       title = "Over-correction check: own markers in their cluster before and after SoupX")

save_plot(p_own, file.path(out, "check_own_markers_retained.png"), width = 14, height = 18)

saveRDS(seu, rds_soupx)

sessionInfo()
