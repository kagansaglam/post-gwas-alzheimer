# METHODS: Big Data GWAS Study Fetching

## Script
`Post_GWAS_Analysis_bigdata.R`

## Purpose
Fetch ALL genome-wide significant SNPs directly from the two largest Alzheimer's Disease GWAS meta-analyses, rather than querying by individual gene names.

## Target Studies

| Study ID | Reference | Sample Size | Notable Findings |
|----------|-----------|-------------|-----------------|
| `GCST90027158` | Bellenguez et al. (2022) *Nature Genetics* | ~111,000 cases/controls | 75 risk loci, including 42 novel |
| `GCST007320` | Kunkle et al. (2019) *Nature Genetics* | ~94,000 cases/controls | Confirmed EPHA1, PTK2B, and others |

## How It Differs from the Hybrid Pipeline

| Feature | Hybrid (`Post_GWAS_Analysis_new.R`) | Big Data (`Post_GWAS_Analysis_bigdata.R`) |
|---------|--------------------------------------|-------------------------------------------|
| SNP source | 19 curated target genes | Full study download (GCST IDs) |
| SNP count | ~46 (backup) or variable | Potentially hundreds |
| Offline fallback | Yes | No — stops if API fails |
| Speed | Fast | Slow (large API payload) |
| Coverage | Targeted | Comprehensive |

## When to Use This Script
- You want unbiased, comprehensive coverage of all significant loci
- You have a stable internet connection
- You are willing to wait for larger Qtlizer batch processing

## Known Limitation
The GWAS Catalog API may return a **500 timeout error** for very large studies. If this happens, the script stops with a message directing you to use the hybrid version instead. This is intentional — no silent fallback, to avoid analysing incomplete data without knowing it.

## Batch Processing
SNPs are sent to Qtlizer in batches of **100** (larger than the hybrid pipeline's 50) with a 0.2-second pause between batches to reduce server load and prevent crashes.

## References
- Bellenguez et al. (2022). New insights into the genetic etiology of Alzheimer's disease and related dementias. *Nature Genetics*. https://doi.org/10.1038/s41588-022-01024-z
- Kunkle et al. (2019). Genetic meta-analysis of diagnosed Alzheimer's disease. *Nature Genetics*. https://doi.org/10.1038/s41588-019-0358-2
