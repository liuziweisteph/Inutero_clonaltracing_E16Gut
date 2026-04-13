# Clonal analysis of the developing mouse gut using *in utero* lineage barcoding and scRNA-seq

Code repository for single-cell transcriptome and clonal relations analysis of the **mouse embryonic day 16.5 (E16.5) gut**, combined with **TREX lineage analysis** to dissect lineage relationships amongst diverse gut cell types. The repository also includes a re-analysis of publicly available **E9.5/E10.5 neural crest** scRNA-seq data from the De Haan & He *et al.* studies to identify the gut-innervating neural crest populations targeted by *in utero* E7.5 amniotic cavity lentivirus injection. Together, these scripts generate all figures in the accompanying manuscript.

## Repository structure

```
.
├── Preprocessing/
│   ├── Assessing_ambinetRNA.R        # Ambient RNA assessment (SoupX)
│   ├── Expt1_Preprocessing.R         # E1: QC → doublet removal → merge → SCT → Harmony
│   ├── Expt2_Preprocessing.R         # E2: SoupX → QC → doublet removal → merge → SCT → Harmony
│   ├── Expt1_Processing.Rmd          # E1: annotation, markers, ENS sub-clustering
│   └── Expt2_Processing.Rmd          # E2: annotation, markers, label transfer from E1
├── Clonal_analysis/
│   ├── Barcode_E1_Clones.Rmd         # Clonal analysis — Experiment 1 (Fig. 3)
│   ├── Barcode_E2_Clones.Rmd         # Clonal analysis — Experiment 2 (Fig. 4)
│   ├── Clonal_Regional.Rmd           # Regional clonal sharing, ENS migration (Fig. 5)
│   ├── E1E2_Comparison.Rmd           # Cross-experiment comparison (Supplementary)
│   ├── Transcriptome_Regional.Rmd    # Merged E1 + E2 transcriptome integration (Supplementary)
│   ├── TREX_output_E1/               # TREX outputs (clones.txt, clone_details.txt, umi_count_matrix.csv, log.txt, …)
│   └── TREX_output_E2/               # TREX outputs for Experiment 2
├── E9E10/
│   └── E9E10_analysis.html           # Reanalysis of E9.5/E10.5 neural crest data (Scanpy/Python)
├── TREX/
│   └── medium_dataset_exclusion_list.csv   # CloneID sequences passed to TREX --filter-cloneids (overrepresentation filter)
├── Helpers/
│   ├── SoupXDoubletfinder.R          # Custom functions for SoupX correction, DoubletFinder, elbow plots
│   └── helpers.R                     # Custom functions for Seurat object creation, QC, batch processing
├── Data/                             # Intermediate outputs and tables (not tracked by git)
├── renv.lock                         # Reproducible R environment snapshot
├── renv/                             # renv infrastructure
├── LICENSE                           # MIT License
└── README.md
```

## Data flow

The diagram below shows how scripts are linked through shared data objects. All paths are **relative to the project root**.

```
E9.5/E10.5 public data (De Haan & He *et al.* DOI: 10.1126/science.adq9248; ArrayExpress E-MTAB-14817)
│
└─ E9E10/E9E10_analysis ─────────────── Neural crest subset (NC_AC) → Fig. 1, Supp. Fig. 2
   Identifies vagal gut-innervating NCC population targeted by TREX barcoding


Raw 10X counts — E16.5 gut (GEO: GSE325733)
│
├─ Assessing_ambinetRNA.R ──────────────── Data/ambient_RNA_QC_summary*.csv
│
├─ Expt1_Preprocessing.R ──────────────── Data/E16_barcodes_E1.rds
│   │                                      Data/mmus_s_genes.txt
│   │                                      Data/mmus_g2m_genes.txt
│   ▼
│   Expt1_Processing.Rmd ─────────────── Data/E16_barcodes_E1.rds  (annotated)
│                                         Data/E16_barcodes_E1_ENSsubset.rds
│
├─ Expt2_Preprocessing.R ──────────────── Data/E16_barcodes_E2.rds
│   │
│   ▼
│   Expt2_Processing.Rmd ─────────────── Data/E16_barcodes_E2.rds     (annotated + label transfer)
│
├─ TREX (run10x) ◄── BAM + cell barcode list + TREX/medium_dataset_exclusion_list.csv (--filter-cloneids)
│       │            Per embryo: if merging regional BAMs, give each region unique cell-barcode IDs before samtools merge; then index merged BAM (see TREX section)
│       │            Drops overrepresented cloneID sequences before graph / clone calling
│       ▼
│   TREX_output_E1/ & TREX_output_E2/  (clones.txt, clone_details.txt, …)
│
├─ Barcode_E1_Clones.Rmd ◄── E16_barcodes_E1.rds + TREX_output_E1/
│       │
│       ▼ overwrites Data/E16_barcodes_E1.rds (adds clone metadata)
│       ▼ Data/E1_df_clonalanalysis_*.csv  → lineage coupling (Bandler *et al.* format)
│
├─ Barcode_E2_Clones.Rmd ◄── E16_barcodes_E2.rds + TREX_output_E2/
│       │
│       ▼ writes Data/E16_barcodes_E2.rds, E16_barcodes_ENS_E2.rds, etc. (adds clone metadata)
│       ▼ Data/E2_df_clonalanalysis_*.csv  → lineage coupling (Bandler *et al.* format)
│
├─ E1E2_Comparison.Rmd ◄──── E16_barcodes_E1.rds + E16_barcodes_E2.rds
│
├─ Clonal_Regional.Rmd ◄──── E16_barcodes_E1.rds + E16_barcodes_E2.rds + E16_barcodes_ENS_E2.rds
│
└─ Transcriptome_Regional.Rmd ◄── E16_barcodes_E1.rds + E16_barcodes_E2.rds
        │
        ▼ Data/E16_barcodes_E1E2_merged.rds (joint Harmony integration)
```

## Raw data access

Raw sequencing and processed archives for the two arms of this project are hosted in public repositories:

| Data | Description | Repository |
|------|-------------|------------|
| **E9.5–E10.5 neural crest** | Whole-embryo scRNA-seq used for the neural-crest reanalysis (`E9E10/`; De Haan & He *et al.*) | **[ArrayExpress E-MTAB-14817](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14817)** |
| **E16.5 barcoded guts** | *In utero* lineage–barcoded gut scRNA-seq (this study: TREX, Seurat, clonal analysis) | **[GEO GSE325733](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325733)** |

Download raw FASTQs / matrices from those accessions before running the E9E10 notebook or the E16 preprocessing and TREX pipelines.

## Execution order

Scripts should be run in the following order to reproduce the full analysis:


| Step | Script                       | Language        | Purpose                                                                                                                                                                      |
| ---- | ---------------------------- | --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0    | `E9E10/E9E10_analysis`       | Python (Scanpy) | Reanalysis of E9.5/E10.5 neural crest — identify gut-innervating vagal NCC                                                                                                   |
| 1    | `Assessing_ambinetRNA.R`     | R               | Quantify ambient RNA contamination                                                                                                                                           |
| 2a   | `Expt1_Preprocessing.R`      | R               | E1 QC, doublet removal, SCTransform, Harmony                                                                                                                                 |
| 2b   | `Expt2_Preprocessing.R`      | R               | E2 SoupX, QC, doublet removal, SCTransform, Harmony                                                                                                                          |
| 3a   | `Expt1_Processing.Rmd`       | R               | E1 cell-type annotation, DEG analysis, ENS sub-clustering                                                                                                                    |
| 3b   | `Expt2_Processing.Rmd`       | R               | E2 cell-type annotation, label transfer from E1                                                                                                                              |
| 3c   | TREX `run10x`                | —               | Per experiment: extract barcodes from BAM with `**--filter-cloneids TREX/medium_dataset_exclusion_list.csv**`; write `Clonal_analysis/TREX_output_E1/` and `TREX_output_E2/` |
| 4a   | `Barcode_E1_Clones.Rmd`      | R               | E1 clonal integration (requires TREX output)                                                                                                                                 |
| 4b   | `Barcode_E2_Clones.Rmd`      | R               | E2 clonal integration (requires TREX output)                                                                                                                                 |
| 6a   | `E1E2_Comparison.Rmd`        | R               | Cross-experiment clone comparison                                                                                                                                            |
| 6b   | `Clonal_Regional.Rmd`        | R               | Regional clonal dispersion, ENS migration                                                                                                                                    |
| 6c   | `Transcriptome_Regional.Rmd` | R               | Merged E1 + E2 transcriptome                                                                                                                                                 |


## Setup and reproducibility

### Prerequisites

- **R 4.5.0** (as recorded in `renv.lock`)
- **Bioconductor 3.21**
- **Python 3** with `scanpy`, `harmonypy`, `numpy`, `pandas`, `matplotlib` (for `E9E10/E9E10_analysis`)
- **Python** with `leidenalg` (accessed via `reticulate`; required for Leiden clustering in R scripts)
- **TREX** ([Ratz *et al.](https://github.com/frisen-lab/TREX)*) — run separately on each experiment’s merged 10x BAM to populate `Clonal_analysis/TREX_output_E1/` and `TREX_output_E2/`. Pass `**TREX/medium_dataset_exclusion_list.csv`** to TREX `**--filter-cloneids**` so highly overrepresented cloneID sequences are ignored before clone calling (see `TREX_output_*/log.txt` for the exact command).
- **Lineage coupling (optional)** — Python environment for [Bandler *et al.* lineage coupling](https://github.com/mayer-lab/Bandler-et-al_lineage) (`Lineage Coupling Analysis/`); see below.

### Installation

```bash
git clone <repository-url>
cd <repository>
```

Open an R session in the project root and restore the exact package environment:

```r
install.packages("renv")
renv::restore()
```

### Key package versions (from `renv.lock`)


| Package        | Version | Source       |
| -------------- | ------- | ------------ |
| Seurat         | 5.4.0   | CRAN         |
| SeuratObject   | 5.3.0   | CRAN         |
| sctransform    | 0.4.3   | CRAN         |
| DoubletFinder  | 2.0.6   | GitHub       |
| SoupX          | 1.6.2   | CRAN         |
| DropletUtils   | 1.28.1  | Bioconductor |
| scran          | 1.36.0  | Bioconductor |
| biomaRt        | 2.64.0  | Bioconductor |
| ComplexHeatmap | 2.24.1  | Bioconductor |
| UpSetR         | 1.4.0   | CRAN         |
| ggplot2        | 4.0.1   | CRAN         |
| ggalluvial     | 0.12.5  | CRAN         |
| dplyr          | 1.2.0   | CRAN         |
| tidyr          | 1.3.2   | CRAN         |
| scales         | 1.4.0   | CRAN         |
| patchwork      | 1.3.2   | CRAN         |
| circlize       | 0.4.17  | CRAN         |
| scCustomize    | 3.2.4   | CRAN         |
| reticulate     | 1.44.1  | CRAN         |
| renv           | 1.1.6   | CRAN         |


All ~200 dependencies are pinned in `renv.lock`. Run `renv::restore()` for a fully reproducible environment.

### Data availability

See **[Raw data access](#raw-data-access)** above: **E-MTAB-14817** (E9.5–E10.5 neural crest) and **GSE325733** (E16.5 barcoded guts).

To skip E16 preprocessing and jump straight to analysis, download the two base Seurat objects from **GSE325733** and place them in `Data/`:


| File                  | Description                                                   |
| --------------------- | ------------------------------------------------------------- |
| `E16_barcodes_E1.rds` | Experiment 1 — merged, QC'd, Harmony-integrated Seurat object |
| `E16_barcodes_E2.rds` | Experiment 2 — merged, QC'd, Harmony-integrated Seurat object |


These are the starting inputs for `Expt1_Processing.Rmd` and `Expt2_Processing.Rmd`, respectively. All downstream scripts consume the annotated objects produced by those two notebooks.

If you wish to **re-run preprocessing from raw counts**, the Cell Ranger output directories are also available from GEO. The preprocessing scripts expect them at a local path set via:

```r
external_data_dir <- "~/scRNA_seq/Trex_ENS/Submission/Data"
```

**Edit this path** to match where you have extracted the raw count matrices.

### E9.5/E10.5 neural crest reanalysis (`E9E10/`)

The `E9E10/E9E10_analysis` notebook is a **Python/Scanpy** reanalysis of publicly available single-cell RNA-seq data from the De Haan and He studies. These datasets profile the **E9.5 and E10.5 mouse neural crest** at whole-embryo level.

**Purpose:** To identify the targeted gut-innervating vagal neural crest cell (NCC) population that is subsequently barcoded by TREX in our E16 experiments. This provides developmental context for the lineage barcoding strategy (Fig. 1).

**Analysis workflow:**

1. Load the full E9.5/E10.5 dataset (`E9E10_all.h5ad`) and the neural/epithelial subset (`E9E10NC.AC_neural_and_epi.h5ad`)
2. Subset the **NC and sensory neuron** population from the full atlas
3. Filter to **samples**  — E9AC1, E9AC21, E9AC22 (E9.5) and E10AC (E10.5)
4. Re-normalise, identify highly variable genes (batch-corrected by `sample_id`), and run PCA
5. Batch correct with **Harmony** across samples
6. Leiden clustering (resolution 0.4 → 13 clusters) and marker gene analysis
7. Characterise clusters using curated marker panels for vagal, cranial, sympathetic, parasympathetic, and ENS identity
8. Hox gene expression analysis to map axial identity and confirm vagal (Hoxb2–Hoxb4-positive) gut-innervating NCCs

**Input data:** The source h5ad files are from the published De Haan & He datasets. The processed NC_AC object is saved as `NC_AC_E9E10.h5ad`.

**Key Python dependencies:** `scanpy`, `harmonypy`, `numpy`, `pandas`, `matplotlib`

### Regional libraries from the same embryo (before TREX)

Gut regions were processed as **separate Cell Ranger runs** per embryo, so each `**outs/possorted_genome_bam.bam`** uses the **same raw 10x barcode space** and barcodes would **collide** across regions if you merged BAMs without relabelling. **Before** `samtools merge` (or equivalent), **rewrite read tags** so every retained cell’s **CB** tag matches the **prefixed** barcode naming used in this project: `**<library_prefix>_` + 16-character 10x barcode**, exactly as in the allowlists. The library prefixes match `**merge(..., add.cell.ids = ...)`** in `Preprocessing/Expt1_Preprocessing.R` and `Expt2_Preprocessing.R`:


| Experiment | `add.cell.ids` (Seurat merge) | Full cell ID pattern in TSV / BAM                                                  |
| ---------- | ----------------------------- | ---------------------------------------------------------------------------------- |
| **E1**     | `e16_s1`, `e16_i1`            | `e16_s1_` + 16 bp barcode and `e16_i1_` + 16 bp barcode                            |
| **E2**     | `e16_s1`, `e16_ji1`, `e16_c1` | `e16_s1_`, `e16_ji1_`, or `e16_c1_`, each followed by the 16-character 10x barcode |



| File                  | Role                                                                                                                      |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `Data/cellids_E1.tsv` | One barcode per line (tab-separated index and full cell ID); only **`e16_s1_`** and **`e16_i1_`** prefixes (Experiment 1) |
| `Data/cellids_E2.tsv` | Same layout for Experiment 2: **`e16_s1_`**, **`e16_ji1_`**, and **`e16_c1_`** |


**Cell identity in the merged BAM must map to these files:** the **CB** strings in the BAM after prefixing must be the **same** strings as column 2 of `cellids_E1.tsv` / `cellids_E2.tsv` (no extra characters). To add a library prefix inside the BAM (here **`e16_s1_`** for one regional library), stream the alignments through **`awk`** and write a new BAM; repeat per region with the matching prefix (`e16_i1_`, `e16_ji1_`, `e16_c1_`) before **`samtools merge`**.

```bash
samtools view -h possorted_genome_bam.bam \
  | awk 'BEGIN{OFS="\t"}
         /^@/ {print; next}
         {
           for (i=12; i<=NF; i++) {
             if ($i ~ /^CB:Z:/) {
               sub(/^CB:Z:/, "", $i);
               $i = "CB:Z:e16_s1_" $i;
             }
           }
           print
         }' \
  | samtools view -b -o possorted_genome_bam_e16s1_prefixed.bam
```

Repeat with the same **`awk`** pattern but **`CB:Z:e16_ji1_`** / **`CB:Z:e16_c1_`** (and output filenames) for the other Experiment 2 libraries, then merge, sort, and index—for example for **Embryo 2** (three regions: **`e16_s1`**, **`e16_ji1`**, **`e16_c1`**):

```bash
samtools merge -@ 32 -o merged_e16_Embryo2.bam \
  possorted_genome_bam_e16s1_prefixed.bam \
  possorted_genome_bam_e16ji1_prefixed.bam \
  possorted_genome_bam_e16c1_prefixed.bam

samtools sort -@ 60 -o merged_e16_Embryo2.sorted.bam merged_e16_Embryo2.bam
samtools index merged_e16_Embryo2.sorted.bam
```

Install the result as **`outs/possorted_genome_bam.bam`** (and **`outs/possorted_genome_bam.bam.bai`**) under the directory you pass to **`trex run10x`**, or symlink/copy `merged_e16_Embryo2.sorted.bam` to that name.

**`-f`** must list the **same** barcodes as in the BAM (e.g. second column of `cellids_E1.tsv` / `cellids_E2.tsv` via **`cut -f2`**, or an equivalent one-column file). Then run **`trex run10x`** on that **`outs/`** directory (plus **`--filter-cloneids`** as below).

### TREX cloneID filtering (`TREX/`)

TREX is run on the merged 10x **`possorted_genome_bam.bam`** with a per-experiment cell-barcode list (see `TREX_output_*/log.txt` for the exact paths used on the cluster). To reduce bias from **overrepresented cloneID sequences** (e.g. dominant or artefactual barcodes), the same **`--filter-cloneids`** file is supplied for both experiments:


| File                                     | Role                                                                                                              |
| ---------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `TREX/medium_dataset_exclusion_list.csv` | One cloneID sequence per line; TREX **ignores** these IDs during analysis so they do not inflate clone statistics |


After filtering, TREX writes **`clones.txt`**, **`clone_details.txt`**, **`cells_filtered.txt`**, **`umi_count_matrix.csv`**, and related files under `Clonal_analysis/TREX_output_E1/` and `TREX_output_E2/`. The clonal R Markdown notebooks merge these calls into the Seurat objects.

### Lineage coupling analysis (Bandler *et al.*)

Clonal coupling heatmaps in this manuscript follow the **lineage coupling** framework published with Bandler *et al.* The reference implementation and Python dependencies are in the Mayer lab repository: **[mayer-lab/Bandler-et-al_lineage](https://github.com/mayer-lab/Bandler-et-al_lineage)** (see the **Lineage Coupling Analysis** folder and `2_lineage_coupling_analysis_wagner_way_5.py`).

`Barcode_E1_Clones.Rmd` and `Barcode_E2_Clones.Rmd` export **CSV tables** suitable as input to the **lineage coupling** Python workflow in [Bandler-et-al_lineage](https://github.com/mayer-lab/Bandler-et-al_lineage): each row is one cell, with columns `**cloneID*`* and `**ident**` (cell-type label), matching the format expected by `**2_lineage_coupling_analysis_wagner_way_5.py**` (`--help` for options). Exports use `write.csv(..., row.names = TRUE)` so **cell barcodes** are preserved in the first column when present. Files are written to `Data/`:


| File                                                  | Contents                                     | Manuscript panels |
| ----------------------------------------------------- | -------------------------------------------- | ----------------- |
| `E1_df_clonalanalysis_subtypes.csv`                   | E1 — subtypes                                | Fig. 3F           |
| `E1_df_clonalanalysis_subtypes_meso.csv`              | E1 — mesodermal-only clones, subtype `ident` | Fig. 3G           |
| `E2_df_clonalanalysis_subtypes_nonepiblast.csv`       | E2 — subtypes, non-epiblast                  | Fig. 4H           |         |
| `E2_df_clonalanalysis_subtypes_meso.csv`              | E2 — mesodermal-label subset, subtypes       | Fig. 4I           |


Clone the Bandler-et-al repository, create the **conda** environment from `**requirements.txt`** in **Lineage Coupling Analysis/** as described in [their README](https://github.com/mayer-lab/Bandler-et-al_lineage), then run the coupling script on the exported CSV. The upstream `**1_generate_input_for_lineage_coupling_analysis_*.R`** in that repository is an alternative way to build compatible inputs from a Seurat object; this project instead emits CSVs directly from annotated metadata.

### Running R Markdown notebooks

All `.Rmd` files include `knitr::opts_knit$set(root.dir = normalizePath(".."))` so that `Data/`, `Figures/`, and `Clonal_analysis/` resolve from the **project root** regardless of whether you knit from the subdirectory or the root.

## Figure mapping


| Figure                             | Primary script(s)                                                                                                                                      |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Fig. 1                             | `E9E10/E9E10_analysis` → `E9E10/Figures/`                                                                                                              |
| Supp. Fig. 1                       | `E9E10/E9E10_analysis` → `E9E10/Figures/`                                                                                                              |
| Fig. 2                             | `Expt1_Processing.Rmd` → `Figures/Fig2/`                                                                                                               |
| Fig. 3                             | `Barcode_E1_Clones.Rmd` → `Figures/Fig3/`                                                                                                              |
| Fig. 4                             | `Expt2_Processing.Rmd`, `Barcode_E2_Clones.Rmd` → `Figures/Fig4/`                                                                                      |
| Fig. 5                             | `Clonal_Regional.Rmd`, `E1E2_Comparison.Rmd`, `Transcriptome_Regional.Rmd` → `Figures/Fig5/`                                                           |
| Fig. 3F–G, 4H–I (lineage coupling) | Exported `Data/E1_df_clonalanalysis_*.csv`, `E2_df_clonalanalysis_*.csv` + [Bandler-et-al_lineage](https://github.com/mayer-lab/Bandler-et-al_lineage) |


## Variable naming convention

Across all scripts, the main Seurat objects follow a consistent naming scheme:


| Variable                  | Content                                        |
| ------------------------- | ---------------------------------------------- |
| `E16_seurat_E1`           | Experiment 1 — annotated E16 gut Seurat object |
| `E16_seurat_E2`           | Experiment 2 — annotated E16 gut Seurat object |
| `E16_seurat_E1_ENSsubset` | E1 ENS-only subset (re-clustered)              |



## License

This project is licensed under the **MIT License** — see `[LICENSE](LICENSE)` for details.

Copyright (c) 2025 Ziwei Liu, Karolinska Institutet, 2026
