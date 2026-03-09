# METHODS: GO Biological Process Enrichment

## Tool
**clusterProfiler** — Bioconductor package for gene ontology and pathway enrichment

## Purpose
Identify overrepresented Gene Ontology (GO) Biological Process terms among candidate genes, with redundancy reduction via semantic similarity clustering.

## How It Works
1. **Over-representation analysis (ORA)** is performed against the GO Biological Process ontology using a hypergeometric test
2. Results are **simplified** using semantic similarity (Wang method) to remove redundant GO terms — closely related terms are merged, keeping only the most significant representative
3. Results visualized as a dotplot showing gene ratio and adjusted p-value

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ont` | `"BP"` | GO sub-ontology: Biological Process |
| `keyType` | `"ENTREZID"` | Input gene ID type |
| `pAdjustMethod` | `"BH"` | Benjamini-Hochberg FDR correction |
| `pvalueCutoff` | 0.05 | Significance threshold |
| `readable` | `TRUE` | Show gene symbols in output |
| `simplify cutoff` | 0.7 | Semantic similarity threshold for term merging |

## Why Simplify?
GO annotations are hierarchical and highly redundant — many terms describe overlapping concepts. Without simplification, results lists are dominated by semantically near-identical terms. The `simplify()` function removes this redundancy based on semantic distance between GO terms.

## Output
- `GO_Results.csv` — simplified GO:BP terms with gene counts and statistics
- Dotplot visualization: top 15 simplified terms

## References
- Wu et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*. https://doi.org/10.1016/j.xinn.2021.100141
- Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics*. https://doi.org/10.1093/genetics/iyad031
