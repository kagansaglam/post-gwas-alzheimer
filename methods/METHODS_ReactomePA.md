# METHODS: Reactome Pathway Enrichment

## Tool
**ReactomePA** — Bioconductor package for Reactome pathway over-representation analysis

## Purpose
Identify biological pathways from the Reactome database that are significantly overrepresented among the candidate genes identified from eQTL mapping.

## How It Works
ReactomePA performs a **hypergeometric test** (over-representation analysis, ORA) comparing the proportion of candidate genes belonging to each Reactome pathway against the expected proportion based on the background genome.

The test asks: *"Given how many genes are in this pathway by chance, is it surprising that this many of our candidate genes fall in it?"*

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `organism` | `"human"` | Species for pathway database |
| `pAdjustMethod` | `"BH"` | Benjamini-Hochberg FDR correction |
| `pvalueCutoff` | 0.05 | Significance threshold (adjusted p-value) |
| `readable` | `TRUE` | Convert Entrez IDs to gene symbols in output |

## Input
Entrez Gene IDs derived from eQTL-mapped candidate genes (via `org.Hs.eg.db` conversion).

## Output
- `Reactome_Results.csv` — significant pathways with gene counts, p-values, and gene lists
- Dotplot visualization (shown during runtime): top 15 pathways ranked by gene ratio

## Notes
Reactome is a curated, peer-reviewed pathway database. It is particularly strong for immune system, metabolism, and cell signalling pathways — all relevant to AD neuroinflammation biology.

## References
- Yu & He (2016). ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. *Molecular BioSystems*. https://doi.org/10.1039/C5MB00663E
- Gillespie et al. (2022). The reactome pathway knowledgebase 2022. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gkab1028
