# snRNA-seq humane Herzschnitte: Semaglutide vs. Kontrolle

Single-nucleus-RNA-seq-Analyse humaner Herzschnitte nach Behandlung mit Semaglutide im Vergleich zur unbehandelten Kontrolle.

## Fragestellung

- Welche Zelltypen sind in den Schnitten vertreten, und verschiebt sich ihre Zusammensetzung unter Semaglutide?
- Welche Gene sind innerhalb einzelner Zelltypen durch Semaglutide verändert, pro Patient und als Pseudobulk über alle Patienten?
- Gibt es behandlungsspezifische Subcluster bzw. Zellzustände?

## Design

Vier Patienten, aus jedem Patienten ein unbehandelter (Ctr) und ein behandelter (Semaglutide) Schnitt. Gepaartes Design, Pseudobulk-Modell `~ patient + treatment`.

| Sample-ID | Run | OCM | Patient | Spender  | Behandlung  | Sample im Objekt  |
|-----------|-----|-----|---------|----------|-------------|-------------------|
| S01s1     | S01 | OB1 | P1      | H.25.50  | Ctr         | P1_Ctr            |
| S01s2     | S01 | OB2 | P2      | H.25.121 | Ctr         | P2_Ctr            |
| S01s3     | S01 | OB3 | P4      | H.25.167 | Ctr         | P4_Ctr            |
| S01s4     | S01 | OB4 | P5      | H.25.39  | Ctr         | P5_Ctr            |
| S02s1     | S02 | OB1 | P1      | H.25.50  | Semaglutide | P1_Semaglutide    |
| S02s2     | S02 | OB2 | P2      | H.25.121 | Semaglutide | P2_Semaglutide    |
| S02s3     | S02 | OB3 | P4      | H.25.167 | Semaglutide | P4_Semaglutide    |
| S02s4     | S02 | OB4 | P5      | H.25.39  | Semaglutide | P5_Semaglutide    |

### Wichtige Einschränkung

Behandlung und 10x-Run sind vollständig konfundiert: Alle Kontrollen liefen in S01, alle Semaglutide-Proben in S02. Behandlungseffekte sind daher nicht von Run-Effekten (Ambient-Profil, Sequenziertiefe, Chip) trennbar. Konsequenzen:

- keine Integration über `library`; falls integriert wird, dann über `patient` und nur für die Annotation
- Ambient-RNA-Korrektur (SoupX) pro Sample
- DE auf unkorrigierten Counts

## Daten

- Chemie: 10x GEM-X Universal 3' v4, 4-plex On-Chip Multiplexing (OCM)
- Prozessierung: `cellranger multi` durch die Core Facility (RCUG), Auftrag 2026_040595_2, Libraries N8773 (S01) und N8774 (S02)
- Verwendet werden ausschließlich die demultiplexten Matrizen pro Sample (`sample_filtered_feature_bc_matrix`, `sample_raw_feature_bc_matrix`), nicht das Aggregat
- Rohdaten liegen unter `data/raw/` und sind nicht im Repository

## Projektstruktur

```
/project/ag-baer/Semaglutide/
├── scripts/      Analyse-Skripte, nummeriert, ein Arbeitsschritt pro Skript
├── data/
│   ├── raw/      Cell-Ranger-Output (nicht versioniert)
│   └── rds/      Seurat-Objekte pro Stufe (nicht versioniert)
└── output/       Plots und Tabellen pro Skript (nicht versioniert)
```

## Pipeline

| Skript              | Inhalt                                                           | Output                 |
|---------------------|------------------------------------------------------------------|------------------------|
| `00_config.R`       | Pfade, Sample-Tabelle, Faktorlevel, Farben, RDS-Pfade            |                        |
| `01_import.R`       | Import der mtx-Matrizen, Metadaten, percent.mt / percent.ribo    | `hs_sema_raw.rds`      |
| `02_qc.R`           | QC-Plots, Schwellen, Filter, Ambient-Korrektur, Doublets         | `hs_sema_qc.rds`       |
| `03_cluster.R`      | Normalisierung, Zellzyklus, PCA-Diagnostik, Clustering, UMAP     | `hs_sema_clustered.rds`|
| `04_annotation.R`   | Marker-Panels, FindAllMarkers, Zelltyp-Annotation                | `hs_sema_annotated.rds`|
| `05_deg.R`          | DE pro Zelltyp, pro Patient und als Pseudobulk                   |                        |
| `06_composition.R`  | Zelltyp-Anteile pro Sample und Behandlung                        |                        |

Alle Skripte werden aus dem Projektverzeichnis gestartet und lesen zuerst `scripts/00_config.R`. Parameter wie QC-Schwellen und Clustering-Resolution werden erst nach Sichtung der Plots festgelegt und im jeweiligen Skript festgehalten.

## Umgebung

- R mit Seurat v5, tidyverse, colorspace
- HPC mit SLURM, RStudio im Singularity-Container

## Konventionen

- Seurat-Objekt heißt projektweit `seu`
- Ctr ist immer das erste Faktorlevel
- Barcodes tragen das Sample als Präfix (`P1_Ctr_<barcode>`)
- HGNC-konforme Gensymbole
- Code ohne Kommentare; Entscheidungen werden in diesem README dokumentiert

## Entscheidungen und offene Punkte

### Abgeschlossen

- Gepaartes Design mit Patient als Block
- Import der demultiplexten mtx-Matrizen pro Sample statt Aggregat

### Offen

- Referenzgenom und Intron-Zählung (`config.csv` prüfen)
- Präparation der Schnitte: Kulturdauer, Semaglutide-Konzentration, Herzregion
- Geschlecht der Spender für den Geschlechts-Check
- Geringere Kernzahl bei P4_Ctr gegenüber P4_Semaglutide prüfen
- Ambient-Korrektur: SoupX vs. decontX
- Integration ja/nein
- Genotyp-Abgleich der Paare über die BAM-Dateien (optional)
