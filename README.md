# 🧠 Post-GWAS Analysis Pipeline — Alzheimer's Disease

![R](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![License](https://img.shields.io/badge/License-MIT-blue)

A robust, hybrid post-GWAS analysis pipeline for **Alzheimer's Disease (AD)**, designed to go beyond raw GWAS hits and extract biological meaning from associated variants. The pipeline maps lead SNPs to candidate genes via brain eQTLs, performs multi-database functional enrichment, and is built with fault-tolerance in mind (server fallback logic included).

---

## 📌 Table of Contents

- [Overview](#overview)
- [Pipeline Structure](#pipeline-structure)
- [Methods Used](#methods-used)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [Future Directions](#future-directions)
- [References](#references)

---

## Overview

Genome-wide association studies (GWAS) identify thousands of disease-associated variants, but translating these SNPs into biological mechanisms remains a major challenge. This pipeline implements a **post-GWAS functional annotation workflow** focused on Alzheimer's Disease, integrating:

- GWAS Catalog API querying for AD-associated SNPs
- SNP-to-gene mapping through brain-specific eQTLs (via Qtlizer)
- Multi-database pathway and disease enrichment analysis

The pipeline is designed to be **robust** — if external APIs are unavailable, it gracefully falls back to curated literature-based SNP and gene lists, ensuring analysis always completes.

**Target trait:** Alzheimer's Disease  
**Genome build:** GRCh38  
**Significance threshold:** p ≤ 5×10⁻⁸  
**LD threshold:** r² ≥ 0.8  
**eQTL tissue filter:** Brain

---

## Pipeline Structure

```
post-gwas-alzheimer/
├── README.md
├── scripts/
│   └── Post_GWAS_Analysis_new.R   # Main pipeline script
├── results/                        # Output files (generated at runtime)
│   ├── lead_snps.csv
│   ├── qtls_filtered.csv
│   ├── gene_mapping.csv
│   ├── Reactome_Results.csv
│   ├── GO_Results.csv
│   ├── gProfiler2_Results.csv
│   ├── EnrichR_Results_TopDB.csv
│   └── DOSE_Results.csv
└── docs/
    └── methods.md                  # Detailed methods descriptions
```

---

## Methods Used

### 1. 🔍 GWAS Data Retrieval — `gwasrapidd`
Queries the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) for associations linked to known AD risk genes. Falls back to a curated list of 46 well-validated AD SNPs from published literature if the API is unavailable.

**Key genes queried:** APOE, BIN1, CLU, ABCA7, CR1, PICALM, TREM2, SORL1, PLCG2, ADAM10, and more.

### 2. 🧬 SNP → Gene Mapping — `Qtlizer`
Maps lead SNPs (and their LD proxies) to target genes using expression quantitative trait loci (eQTL) data, filtered to brain-relevant tissues. Batch processing is implemented to prevent API timeouts.

### 3. 🔄 Gene ID Conversion — `clusterProfiler::bitr` + `org.Hs.eg.db`
Converts gene symbols to Entrez IDs for downstream enrichment tools. Handles antisense and versioned gene names via regex cleaning.

### 4. 📊 Pathway & Disease Enrichment
Multiple complementary enrichment databases are queried in parallel:

| Tool | Database | What it tests |
|------|----------|---------------|
| `ReactomePA` | Reactome | Curated biological pathways |
| `clusterProfiler` | Gene Ontology (BP) | Biological process terms |
| `gprofiler2` | GO + KEGG + Reactome | Cross-database enrichment |
| `enrichR` | DisGeNET + GO 2023 | Disease-gene associations |
| `DOSE` | Disease Ontology | DO-based disease enrichment |

All results use **Benjamini-Hochberg** FDR correction at p < 0.05.

---

## Installation

```r
# Install CRAN packages
install.packages(c("tidyverse", "gwasrapidd", "httr", "jsonlite",
                   "gprofiler2", "enrichR", "ggplot2"))

# Install Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Qtlizer", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "ReactomePA", "DOSE", "AnnotationDbi"))
```

**R version:** ≥ 4.2.0 recommended

---

## Usage

```r
# Clone the repository
# git clone https://github.com/YOUR_USERNAME/post-gwas-alzheimer.git

# Set your working directory and run
source("scripts/Post_GWAS_Analysis_new.R")
```

To change the trait, tissue filter, or significance threshold, edit the parameters block at the top of the script:

```r
trait          <- "Alzheimer_Hybrid_Analysis"
pval_threshold <- 5e-8
qtl_tissue     <- "Brain"   # Change to "all" for no tissue filter
ld_corr        <- 0.8
```

---

## Output Files

| File | Description |
|------|-------------|
| `lead_snps.csv` | Final SNP list used for analysis |
| `qtls_filtered.csv` | SNP-to-gene mapping results from Qtlizer |
| `gene_mapping.csv` | Symbol → Entrez ID conversion table |
| `Reactome_Results.csv` | Reactome pathway enrichment |
| `GO_Results.csv` | GO Biological Process enrichment (simplified) |
| `gProfiler2_Results.csv` | Multi-database enrichment (gProfiler2) |
| `EnrichR_Results_TopDB.csv` | enrichR top database results |
| `DOSE_Results.csv` | Disease Ontology enrichment |

---

## Future Directions

This pipeline is under active development. Planned extensions include:

### 🔬 Stronger Causal Inference
- **Colocalization analysis** (`coloc` / `HyPrColoc`) — formally test whether GWAS and eQTL signals share a causal variant, replacing LD-proxy filtering with probabilistic evidence
- **Mendelian Randomization** (`TwoSampleMR`) — test whether changes in gene expression causally affect AD risk, not just correlate with it

### 🧠 Brain Cell-Type Resolution
- **MAGMA cell-type enrichment** — partition heritability across brain cell types (microglia, astrocytes, oligodendrocytes, neurons) using single-cell RNA-seq reference panels
- **Multi-tissue eQTL comparison** — extend Qtlizer queries across all 13 GTEx brain subregions and compute tissue-specificity scores

### 🎯 Druggability & Translation
- **Drug target prioritization** — integrate candidate genes with OpenTargets and DGIdb to flag druggable targets and existing compounds
- **PheWAS** — query all GWAS Catalog traits for top SNPs to identify pleiotropic effects and shared mechanisms with other diseases

### 📦 Pipeline Infrastructure
- **GWAS summary statistics input** — accept full sumstats (z-scores, betas, SE) to unlock finemapping tools like SuSiE and FINEMAP
- **Automated HTML report** — Quarto/RMarkdown report bundling all plots, tables, and a methods summary
- **Locus zoom plots** — per-locus visualization using `locuszoomr`

---

## References

- Buniello A. et al. (2019) The NHGRI-EBI GWAS Catalog. *Nucleic Acids Research*
- Võsa U. et al. (2021) Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci. *Nature Genetics*
- Krassowski M. et al. — Qtlizer R package
- Yu G. et al. (2012) clusterProfiler. *OMICS*
- Fabregat A. et al. (2018) The Reactome Pathway Knowledgebase. *Nucleic Acids Research*
- Lambert J.C. et al. (2013) Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for AD. *Nature Genetics*
- Jansen I.E. et al. (2019) Genome-wide meta-analysis identifies new loci for AD. *Nature Genetics*

---

## License

MIT License — see [LICENSE](LICENSE) for details.

> **Contact:** Open an issue or pull request for questions and contributions.
