# METHODS: Disease Ontology Enrichment (DOSE)

## Tool
**DOSE** — Bioconductor package for Disease Ontology Semantic and Enrichment analysis

## Purpose
Identify Disease Ontology (DO) terms significantly associated with the candidate gene set, providing a disease-level biological context beyond pathway and GO analysis.

## How It Works
DOSE performs over-representation analysis against the **Human Disease Ontology (DO)** — a standardized ontology that classifies human diseases and their genetic relationships. This allows the pipeline to ask: *"Beyond Alzheimer's Disease itself, what other diseases share a genetic basis with our candidate gene set?"*

This is particularly valuable for:
- Identifying comorbidity patterns (e.g. shared genetics between AD and cardiovascular disease)
- Validating that the pipeline recovers known disease associations
- Discovering unexpected biological connections

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `pvalueCutoff` | 0.05 | Significance threshold (adjusted p-value) |
| `readable` | `TRUE` | Display gene symbols instead of Entrez IDs |

## Input
Entrez Gene IDs from eQTL-mapped candidate genes.

## Output
- `DOSE_Results.csv` — significant disease ontology terms with gene counts and statistics
- Dotplot visualization: top 15 disease terms

## Interpreting Results
Results should recover neurodegenerative diseases (Alzheimer's, dementia) as top hits — this serves as an internal validation of the pipeline. Additional hits (e.g. lipid metabolism disorders, immune diseases) reveal the pleiotropic nature of AD risk loci.

## References
- Yu et al. (2015). DOSE: an R/Bioconductor package for disease ontology semantic and enrichment analysis. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btu684
- Schriml et al. (2022). The Human Disease Ontology 2022 update. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gkab1063
