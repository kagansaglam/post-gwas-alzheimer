# 🧬 post-gwas-alzheimer

![R](https://img.shields.io/badge/R-4.2+-276DC3?style=flat&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.16+-87B13F?style=flat)
![License](https://img.shields.io/badge/license-MIT-green?style=flat)
![Status](https://img.shields.io/badge/status-active-brightgreen?style=flat)
![Disease](https://img.shields.io/badge/disease-Alzheimer's-blueviolet?style=flat)
![Methods](https://img.shields.io/badge/methods-GWAS%20%7C%20eQTL%20%7C%20Colocalization-orange?style=flat)

A robust, modular Post-GWAS analysis pipeline focused on **Alzheimer's Disease (AD)**, designed to move from GWAS-significant SNPs to biologically interpretable gene sets and pathways.

---

## 📌 Overview

Genome-Wide Association Studies (GWAS) identify statistical associations between genetic variants and traits, but the biological interpretation of these signals requires substantial downstream work. This pipeline bridges that gap for Alzheimer's Disease by:

1. Fetching curated AD risk SNPs from the **NHGRI-EBI GWAS Catalog** (with a robust offline fallback)
2. Mapping SNPs to candidate genes using **eQTL data** filtered for brain tissue
3. Converting gene identifiers and running **multi-tool enrichment analysis**
4. Producing pathway, ontology, and disease association results ready for interpretation

The pipeline is built with a **hybrid server/backup architecture** — if any external API is unavailable, the analysis continues using curated literature-based data, making it reproducible in any environment.

---

## 🗂️ Repository Structure

```
post-gwas-alzheimer/
│
├── README.md                        # This file
├── Post_GWAS_Analysis_new.R              # Script 1: Hybrid (gene-based, with fallback)
├── Post_GWAS_Analysis_bigdata.R          # Script 2: Big Data (full study download)
├── Post_GWAS_Colocalization_HyPrColoc.R  # Script 3: Standalone colocalization
├── Post_GWAS_SummaryStats_Pipeline.R     # Script 4: Real summary stats (full pipeline)
├── Post_GWAS_DrugTarget_OpenTargets.R    # Script 5: Drug target prioritization: Multi-trait colocalization
│
├── methods/                         # Per-tool methodology documentation
│   ├── METHODS_GWAS_Catalog.md
│   ├── METHODS_QTLizer.md
│   ├── METHODS_BigData_Studies.md
│   ├── METHODS_HyPrColoc.md
│   ├── METHODS_ReactomePA.md
│   ├── METHODS_clusterProfiler_GO.md
│   ├── METHODS_gProfiler2.md
│   ├── METHODS_enrichR.md
│   └── METHODS_DOSE.md
│
└── results/                         # Output directory (generated at runtime)
    ├── lead_snps.csv
    ├── qtls_filtered.csv
    ├── gene_mapping.csv
    ├── Reactome_Results.csv
    ├── GO_Results.csv
    ├── gProfiler2_Results.csv
    ├── EnrichR_Results_TopDB.csv
    └── DOSE_Results.csv
```

---

## ⚙️ Pipeline Steps

| Step | Description | Tool |
|------|-------------|------|
| 1 | Fetch GWAS associations for AD risk genes | `gwasrapidd` + GWAS Catalog API |
| 2 | Map SNPs to genes via brain eQTLs | `Qtlizer` |
| 3 | Convert gene symbols to Entrez IDs | `org.Hs.eg.db` + `clusterProfiler::bitr` |
| 4a | Reactome pathway enrichment | `ReactomePA` |
| 4b | GO Biological Process enrichment | `clusterProfiler` |
| 4c | Multi-database enrichment (GO, REAC, KEGG) | `gprofiler2` |
| 4d | Curated library enrichment | `enrichR` |
| 4e | Disease Ontology enrichment | `DOSE` |
| 5 | Multi-trait colocalization (GWAS + brain eQTLs) | `HyPrColoc` |
| 6 | Drug target prioritization & clinical evidence | `OpenTargets API` |

---

## 🚀 Getting Started

### Requirements

- R ≥ 4.2.0
- Bioconductor ≥ 3.16
- Internet connection (optional — pipeline works offline via backup data)

### Installation

```r
# Install CRAN packages
install.packages(c("tidyverse", "gwasrapidd", "httr", "jsonlite",
                   "gprofiler2", "enrichR", "readr", "stringr", "ggplot2"))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Qtlizer", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "ReactomePA", "DOSE", "AnnotationDbi"))
```

### Which Script to Use?

| | `Post_GWAS_Analysis_new.R` | `Post_GWAS_Analysis_bigdata.R` | `Post_GWAS_Colocalization_HyPrColoc.R` |
|---|---|---|---|
| **Purpose** | Enrichment analysis | Enrichment (full studies) | Causal gene evidence |
| **SNP source** | 19 curated AD genes | Bellenguez 2022 + Kunkle 2019 | Same as Script 1 |
| **Offline fallback** | ✅ Yes | ❌ No | ✅ Yes |
| **Speed** | Fast | Slow | Moderate |
| **Use when** | First-pass analysis | Full coverage | After enrichment, to confirm causal genes |

### Run

```r
# Step 1: Enrichment analysis (recommended starting point)
source("Post_GWAS_Analysis_new.R")

# Step 1 (alternative): Full study download
source("Post_GWAS_Analysis_bigdata.R")

# Step 2: Colocalization — run after enrichment to confirm causal genes
source("Post_GWAS_Colocalization_HyPrColoc.R")
```

Results will be saved to `~/Documents/run_results/` by default. Change the `work_dir` parameter at the top of the script to modify this.

### Key Parameters

```r
trait          <- "Alzheimer_Hybrid_Analysis"
pval_threshold <- 5e-8      # GWAS significance threshold
ld_method      <- "r2"      # LD correlation method
ld_corr        <- 0.8       # LD r² cutoff
qtl_tissue     <- "Brain"   # Tissue filter for eQTL mapping
```

---

## 🧪 Target Genes

The pipeline queries associations for 19 well-established AD GWAS risk genes:

`APOE`, `BIN1`, `CLU`, `ABCA7`, `CR1`, `PICALM`, `MS4A6A`, `CD33`, `MS4A4E`, `CD2AP`, `EPHA1`, `INPP5D`, `MEF2C`, `TREM2`, `SORL1`, `PLCG2`, `ADAM10`, `ACE`, `TOMM40`

---

## 🔭 Future Directions

This pipeline is designed as a foundation. The following analyses are planned for future development:

### 1. 🔗 Colocalization Analysis (`coloc` / `HyPrColoc`)
Formally test whether a GWAS signal and a brain eQTL share the same causal variant — providing much stronger evidence of gene causality than LD-proxy matching alone. Multi-trait colocalization with `HyPrColoc` will allow simultaneous testing across multiple AD-related phenotypes.

### 2. ⚖️ Mendelian Randomization (`TwoSampleMR`)
Use brain eQTL variants as instruments to test whether altered expression of candidate genes *causally* affects AD risk. This moves beyond correlation and enables directional inference about gene function in disease.

### 3. 🧠 Cell-Type-Specific Enrichment (`MAGMA Celltyping`)
Integrate GWAS summary statistics with single-cell RNA-seq signatures (e.g. Allen Brain Atlas, Human Cell Atlas) to identify which brain cell types — microglia, astrocytes, excitatory neurons — are most enriched for AD heritability.

### 4. 💊 Drug Target Prioritization (OpenTargets / DGIdb)
Map candidate genes to known drug targets and compounds using OpenTargets and DGIdb. Several AD GWAS genes (`TREM2`, `ADAM10`, `PLCG2`) are already actionable targets — this step will systematically flag druggable candidates and existing clinical compounds.

---

## 📚 Methods

Detailed documentation for each tool used in this pipeline is available in the [`methods/`](methods/) directory.

---

## 📄 License

MIT License — free to use, adapt, and distribute with attribution.

---

## 🙋 Author

Contributions and feedback welcome. Please open an issue or submit a pull request.
