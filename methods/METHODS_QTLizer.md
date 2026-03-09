# METHODS: SNP-to-Gene Mapping via eQTLs

## Tool
**Qtlizer** — R package for querying the Qtlizer REST API (GTEx-based eQTL database)

## Purpose
Map GWAS lead SNPs (and their LD proxies) to candidate genes through expression quantitative trait loci (eQTL) data, filtered for brain tissue.

## How It Works
For each input SNP, Qtlizer identifies:
1. The SNP itself and its LD proxies (r² ≥ 0.8)
2. Any eQTL associations in the GTEx database where those variants regulate gene expression
3. Results are filtered to retain only brain-tissue eQTLs

SNPs are processed in batches of 50 to prevent API timeouts.

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ld_method` | `r2` | Linkage disequilibrium metric |
| `ld_corr` | 0.8 | Minimum LD r² threshold for proxy SNPs |
| `qtl_tissue` | `"Brain"` | Tissue keyword filter applied to results |
| `batch_size` | 50 | SNPs per API request |

## Why Brain Tissue?
Alzheimer's Disease is a neurological condition. Restricting eQTLs to brain-derived tissues (GTEx brain subregions) increases biological relevance and reduces noise from peripheral tissue associations.

## Fallback Strategy
If the Qtlizer API fails, the pipeline falls back to a manually curated gene list of core AD GWAS genes (`APOE`, `BIN1`, `CLU`, `PICALM`, etc.) to allow enrichment steps to proceed.

## Output
- `qtls_filtered.csv` — eQTL mappings with gene names and tissue annotations

## References
- Munz et al. (2020). Qtlizer: comprehensive QTL annotation of GWAS results. *Scientific Reports*. https://doi.org/10.1038/s41598-020-75770-7
- GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory effects. *Science*. https://doi.org/10.1126/science.aaz1776
