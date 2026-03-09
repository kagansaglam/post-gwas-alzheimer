# METHODS: Curated Library Enrichment (enrichR)

## Tool
**enrichR** — R package interfacing with the Enrichr web service (Ma'ayan Laboratory)

## Purpose
Perform enrichment analysis against a broad collection of curated gene-set libraries, including disease associations (DisGeNET), GO terms, and Reactome pathways from recent annotation releases.

## How It Works
Enrichr uses three scoring approaches:
1. **Fisher exact test** — classic overlap significance
2. **Z-score** — deviation of rank from expected
3. **Combined score** — log(p-value) × z-score, balancing statistical significance with ranking quality

The combined score is generally the most informative metric for ranking results.

## Libraries Used

| Library | Description |
|---------|-------------|
| `GO_Biological_Process_2023` | GO:BP terms (2023 release) |
| `Reactome_Pathways_2024` | Reactome pathways (2024 release) |
| `DisGeNET` | Disease-gene associations from literature and databases |

## Input
HGNC gene symbols from eQTL-mapped candidate genes.

## Output
- `EnrichR_Results_TopDB.csv` — results from the top-ranked database query

## Why Use enrichR Alongside Other Tools?
Each enrichment tool uses slightly different gene-set versions, background sets, and statistical methods. Consistent findings across ReactomePA, clusterProfiler, gProfiler2, and enrichR provide higher confidence in enriched pathways. The **DisGeNET** library also adds disease-level context not available in the other tools.

## References
- Kuleshov et al. (2016). Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gkw377
- Enrichr web server: https://maayanlab.cloud/Enrichr/
