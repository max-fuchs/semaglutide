---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
# snRNA-seq of human heart slices: Semaglutide vs. control

Single-nucleus RNA-seq analysis of human heart slices treated with semaglutide compared with untreated controls.

## Research questions

- Which cell types are present in the slices, and does their composition shift under semaglutide?
- Which genes are altered by semaglutide within individual cell types, per patient and as pseudobulk across all patients?
- Are there treatment-specific subclusters or cell states?

## Design

Four patients, each contributing one untreated (Ctr) and one treated (Semaglutide) slice. Paired design, pseudobulk model `~ patient + treatment`.

| Sample ID | Run | OCM | Patient | Donor    | Treatment   | Sample in object  |
|-----------|-----|-----|---------|----------|-------------|-------------------|
| S01s1     | S01 | OB1 | P1      | H.25.50  | Ctr         | P1_Ctr            |
| S01s2     | S01 | OB2 | P2      | H.25.121 | Ctr         | P2_Ctr            |
| S01s3     | S01 | OB3 | P4      | H.25.167 | Ctr         | P4_Ctr            |
| S01s4     | S01 | OB4 | P5      | H.25.39  | Ctr         | P5_Ctr            |
| S02s1     | S02 | OB1 | P1      | H.25.50  | Semaglutide | P1_Semaglutide    |
| S02s2     | S02 | OB2 | P2      | H.25.121 | Semaglutide | P2_Semaglutide    |
| S02s3     | S02 | OB3 | P4      | H.25.167 | Semaglutide | P4_Semaglutide    |
| S02s4     | S02 | OB4 | P5      | H.25.39  | Semaglutide | P5_Semaglutide    |

### Key limitation

Treatment and 10x run are fully confounded: all controls were run in S01, all semaglutide samples in S02. Treatment effects therefore cannot be separated from run effects (ambient profile, sequencing depth, chip). Consequences:

- no integration across `library`; if integration is used, it is done across `patient` and only for annotation
- differential expression on uncorrected counts

## Data

- Chemistry: 10x GEM-X Universal 3' v4, 4-plex On-Chip Multiplexing (OCM)
- Processing: `cellranger multi` by the core facility (RCUG), order 2026_040595_2, libraries N8773 (S01) and N8774 (S02)
- Reference: 10x GRCh38-2024-A
- Only the demultiplexed per-sample matrices are used (`sample_filtered_feature_bc_matrix`, `sample_raw_feature_bc_matrix`), not the aggregate
- Raw data live in `data/raw/` and are not part of the repository

## Project structure

```
/project/ag-baer/Semaglutide/
├── scripts/      analysis scripts, numbered, one step per script
├── data/
│   ├── raw/      Cell Ranger output (not versioned)
│   └── rds/      Seurat objects per stage (not versioned)
└── output/       plots and tables per script (not versioned)
```

## Pipeline

| Script              | Content                                                         | Output                  |
|---------------------|-----------------------------------------------------------------|-------------------------|
| `00_config.R`       | paths, sample table, factor levels, colours, RDS paths          |                         |
| `01_import.R`       | import of mtx matrices, metadata, percent.mt / percent.ribo     | `hs_sema_raw.rds`       |
| `02_qc.R`           | QC plots, thresholds, filtering, doublet detection (flag only)  | `hs_sema_qc.rds`        |
| `03_cluster.R`      | normalisation, cell cycle, PCA diagnostics, clustering, UMAP    | `hs_sema_clustered.rds` |
| `03b_cluster_checks.R` | cluster proportions per sample, markers of unclear clusters, ambient check | |
| `03c_soupx.R`       | SoupX ambient correction per sample as separate assay, over-correction checks | `hs_sema_soupx.rds` |
| `04_clean_recluster.R` | removal of doublets and low-quality clusters, re-clustering  | `hs_sema_clean.rds`     |
| `05_annotation.R`   | marker panels, FindAllMarkers, cell type annotation             | `hs_sema_annotated.rds` |
| `06_deg.R`          | DE per cell type, per patient and as pseudobulk                 |                         |
| `07_composition.R`  | cell type proportions per sample and treatment                  |                         |

All scripts are run from the project directory and source `scripts/00_config.R` first. Parameters such as QC thresholds and clustering resolution are set only after inspecting the plots and are recorded in the respective script.

## Environment

- R with Seurat v5, tidyverse, colorspace
- HPC with SLURM, RStudio in a Singularity container

## Conventions

- The Seurat object is called `seu` throughout the project
- Ctr is always the first factor level
- Barcodes carry the sample as prefix (`P1_Ctr_<barcode>`)
- HGNC-compliant gene symbols
- Code without comments; decisions are documented in this README
- All plot titles, labels, tables and documentation in English

## Decisions and open issues

### Resolved

- Paired design with patient as block
- Import of demultiplexed per-sample mtx matrices instead of the aggregate
- Reference 10x GRCh38-2024-A; `include-introns` not set, so the Cell Ranger default applies (introns counted)
- QC thresholds, identical for all samples: genes >= 500, UMI >= 1000, MT < 10 %, complexity (log10 genes / log10 UMI) >= 0.80, no upper cut-off
- MT threshold set to 10 % instead of 5 %: one third to one half of nuclei with MT >= 5 % have >= 1500 genes and are likely intact nuclei (risk of selectively losing cardiomyocytes)
- No filtering on percent.ribo
- Doublets detected with scDblFinder per sample and flagged only (`scDblFinder_class`, `scDblFinder_score`); removal decided in `04_clean_recluster.R`
- Ambient RNA correction with SoupX added after clustering because cardiomyocyte transcripts are widespread in non-cardiomyocyte nuclei and higher in Semaglutide samples; run per sample on the QC-filtered nuclei with the per-sample raw matrix and the resolution 0.4 clusters; corrected counts stored as separate assay `SoupX`, assay `RNA` unchanged
- Every analysis using SoupX counts is checked for over-correction (retention of own cell type markers) and compared with the uncorrected counts
- Cleaning in `04_clean_recluster.R`: removal of all scDblFinder doublets and of initial clusters 11 (doublet cluster), 4 (cardiomyocyte fragments) and 6 (low quality); initial cluster 15 (unclear, patient-derived) and 18 (B cells) kept
- Re-clustering, annotation and composition on the `SoupX` assay; `RNA` assay kept for DE comparison
- Normalisation with LogNormalize; raw counts retained for pseudobulk DE
- No integration; to be revisited only if patient effects obscure cell types
- Cell cycle scored for information, not regressed

### Observations

- P4_Ctr has about half the nuclei of P4_Semaglutide with normal per-nucleus quality, and a lower fraction of MT-rich nuclei
- Semaglutide samples show slightly higher percent.mt in all four pairs and lower UMI counts in three of four pairs; confounded with run S02
- P5 has markedly lower percent.mt in both samples (patient effect)
- Low-quality shoulder at 300 to 1500 genes with elevated MT in all samples, independent of treatment
- QC removes 5.7 to 11.0 % of nuclei per sample (27,214 of 29,845 pass), mainly via the UMI and gene thresholds
- scDblFinder flags 4.8 to 7.8 % per sample (1,855 in total), concentrated at high UMI and gene counts as expected
- A small group of nuclei below the diagonal (few genes for their UMI count) remains in P5, mainly P5_Semaglutide; to be checked after clustering
- Sex check (XIST, RPS4Y1, UTY, DDX3Y, KDM5D): P1 male, P2, P4 and P5 female; Ctr and Semaglutide consistent within each patient, supporting the OCM sample assignment
- Without integration, non-cardiomyocyte cell types mix well across patients; cardiomyocytes separate by patient (P5 in its own cluster)
- PC_8 combines cytoplasmic cardiomyocyte transcripts (MB, TNNC1, MYL2, MYL3, DES, CKM) with percent.mt (r = -0.72) and percent.malat1 (r = 0.57), suggesting a cytoplasmic carry-over axis rather than biology
- Proliferating fibroblasts and pericytes (MKI67, PC_9) and an inflammatory/activation programme (IL6, CXCL8, HMGA2, INHBA, PC_10) are present, consistent with slice culture
- T/NK cells and a small PTPRC+ cluster are largely from P4
- Cluster 4: cytoplasmic cardiomyocyte transcripts (CSRP3, CRYAB, DES, MYL2, TPM1) plus ribosomal and MT genes, low gene counts; damaged cardiomyocytes or cytoplasmic fragments
- Cluster 5 (P5): genuine cardiomyocyte nuclei (TTN, MYO18B, NRAP, MYPN, ACTN2, ANKRD1, high MALAT1)
- Cluster 6: no specific markers, ribosomal and housekeeping genes, low gene counts; low quality
- Cluster 15: neuronal/synaptic profile (GRID2, RIMS2, SYT1, UNC5D, RMST) with XACT and IGF2BP3 and very low MALAT1; sex genes follow the patient (Y genes only in P1, XIST only in the female patients), so patient-derived rather than external contamination; nuclear lncRNAs (MALAT1, XIST) generally low; identity to be resolved in annotation
- Cluster 18: B cells (MS4A1, BANK1, PAX5, CD79A), P4 only
- Cardiomyocyte transcripts are detected in 20 to 80 % of fibroblast, pericyte, endothelial and macrophage nuclei (MB, TNNT2, MYL2, TTN), mostly higher in Semaglutide samples; strong ambient / cytoplasmic carry-over that is confounded with treatment
- SoupX: autoEstCont rho 0.064 to 0.106 per sample (used); CM-marker-based rho 0.11 to 0.23 (not used, likely inflated by physically attached cytoplasm). Median 5 to 9 % of counts removed per nucleus. Soup profile dominated by MALAT1, MT genes, FTH1/FTL, TIMP1, B2M and ribosomal genes, not by cardiomyocyte genes. autoEstCont rho is not systematically higher in Semaglutide samples; the higher CM transcript load there reflects a soup richer in CM transcripts
- SoupX over-correction check passed: own cell type markers retained in their clusters (mostly > 95 %, DCN in fibroblasts about 90 %)
- After SoupX, cardiomyocyte transcripts remain detectable in 10 to 45 % of non-cardiomyocyte nuclei, highest in endothelial cells; likely physical cytoplasmic carry-over that SoupX cannot remove. Cardiomyocyte genes must be treated with caution in DE of non-cardiomyocyte cell types
- Cleaning removed 11.9 to 22.7 % of nuclei per sample; 22,325 nuclei kept
- After cleaning and re-clustering on SoupX (resolution 0.4, 17 clusters): non-cardiomyocyte cell types mix well across patients and treatments; cardiomyocytes still separate by sample (P5 in its own cluster, within P5 Ctr vs. Semaglutide; P1_Semaglutide separate within the main cardiomyocyte cluster)
- New cluster 12 (294 nuclei, from initial fibroblast cluster 0): fibroblast and cardiomyocyte markers, median scDblFinder score 0.23 (others <= 0.08), MT 4 %; residual fibroblast/cardiomyocyte doublets or carry-over
- New cluster 9 (cardiomyocytes, from initial cluster 12) is enriched in Semaglutide samples of P1, P2 and P4
- New cluster 14 expresses SOX10, CDH19, NRXN1 (Schwann cells / glia)
- PC_8 still captures the cytoplasmic carry-over axis (FTH1, MT genes, CRYAB, TNNC1, DES, MB; r = 0.7 with percent.mt and percent.ribo)
- First composition hints (resolution 0.4, before cleaning, not tested): cardiomyocyte clusters 3 and 12 higher in Semaglutide in P1, P2 and P4 (P5 cardiomyocytes in cluster 5 unchanged); pericytes (cluster 1) higher in Semaglutide in all four patients; cluster 15 lower in Semaglutide in all four patients. Compositional data, n = 4, one slice per condition and patient, treatment confounded with run

### Open

- Slice preparation: culture duration, semaglutide concentration, heart region
- Nuclei loading numbers for P4_Ctr
- Low-quality cluster check after clustering (nuclei with 5 to 10 % MT)
- Number of PCs and clustering resolution for `03_cluster.R` (currently 30 PCs, resolution 0.4)
- Genotype matching of the pairs using the BAM files (optional)