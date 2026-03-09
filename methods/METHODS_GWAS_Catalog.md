# METHODS: GWAS Catalog Data Retrieval

## Tool
**gwasrapidd** — R interface to the NHGRI-EBI GWAS Catalog REST API

## Purpose
Fetch genome-wide significant SNP-trait associations for Alzheimer's Disease risk genes from the curated GWAS Catalog database.

## How It Works
The pipeline queries the GWAS Catalog using a list of 19 established AD risk gene names. For each gene, all recorded associations are retrieved, then filtered to retain only genome-wide significant hits (p ≤ 5×10⁻⁸). Risk allele variant IDs (rsIDs) are extracted from the filtered associations.

## Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| `pval_threshold` | 5e-8 | Standard genome-wide significance cutoff |
| `gene_name` | 19 AD risk genes | Query target list |

## Target Gene List
`APOE`, `BIN1`, `CLU`, `ABCA7`, `CR1`, `PICALM`, `MS4A6A`, `CD33`, `MS4A4E`, `CD2AP`, `EPHA1`, `INPP5D`, `MEF2C`, `TREM2`, `SORL1`, `PLCG2`, `ADAM10`, `ACE`, `TOMM40`

## Fallback Strategy
If the GWAS Catalog API is unavailable (timeout, server error), the pipeline activates a **backup SNP list** of 46 well-validated Alzheimer's SNPs curated from published literature. This ensures the pipeline completes in offline or restricted-network environments.

## Output
- `lead_snps.csv` — list of genome-wide significant rsIDs

## References
- Buniello et al. (2019). The NHGRI-EBI GWAS Catalog. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gky1120
- gwasrapidd package: https://rmagno.eu/gwasrapidd/
- GWAS Catalog: https://www.ebi.ac.uk/gwas/
