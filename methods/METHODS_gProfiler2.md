# METHODS: Multi-Database Enrichment (gProfiler2)

## Tool
**gprofiler2** — R client for the g:Profiler web service

## Purpose
Run simultaneous enrichment analysis across multiple databases (GO, Reactome, KEGG) in a single query, providing a broad functional overview of candidate genes.

## How It Works
g:Profiler performs ordered/unordered over-representation analysis against multiple annotation sources simultaneously. It uses a custom **g:SCS (Set Counts and Sizes)** multiple testing correction that accounts for the non-independence of GO term tests (due to the DAG structure of GO), making it more appropriate than naive FDR correction for ontology data.

## Databases Queried

| Source | Description |
|--------|-------------|
| `GO` | Gene Ontology (BP, MF, CC) |
| `REAC` | Reactome pathways |
| `KEGG` | KEGG metabolic and signalling pathways |

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `organism` | `"hsapiens"` | Homo sapiens |
| `correction_method` | `"fdr"` | FDR multiple testing correction |
| `sources` | `c("GO","REAC","KEGG")` | Annotation databases |

## Input
Gene symbols (HGNC) from eQTL-mapped candidate genes.

## Output
- `gProfiler2_Results.csv` — all significant terms across all databases
- Manhattan-style `gostplot` visualization grouped by source database

## Notes
gProfiler2 is particularly useful for a first-pass multi-database overview. Results can be cross-validated against the individual ReactomePA and clusterProfiler runs for consistency.

## References
- Kolberg et al. (2023). g:Profiler — interoperable web service for functional enrichment analysis and gene identifier mapping. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gkad347
- g:Profiler web server: https://biit.cs.ut.ee/gprofiler/
