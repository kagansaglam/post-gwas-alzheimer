# METHODS: Multi-Trait Colocalization (HyPrColoc)

## Script
`Post_GWAS_Colocalization_HyPrColoc.R`

## Purpose
Formally test whether Alzheimer's Disease GWAS signals and GTEx brain eQTL
signals share the same underlying causal variant — providing causal gene
evidence that goes beyond LD-based SNP-to-gene mapping.

---

## What Is Colocalization?

Standard post-GWAS pipelines (like the hybrid pipeline in this repo) link SNPs
to genes using LD proxies and eQTL overlap. This is useful but has a critical
limitation: two signals can appear co-localised simply because they are in
linkage disequilibrium, not because they share a causal variant.

Colocalization analysis asks a more specific question:

> *"Is the variant driving the GWAS association the same variant driving the
> eQTL association, or are they two distinct variants that happen to be
> correlated?"*

---

## Why HyPrColoc Over Standard coloc?

| Feature | `coloc` | `HyPrColoc` |
|---------|---------|-------------|
| Traits per test | 2 (pairwise) | 2 to N simultaneously |
| AD GWAS + multiple eQTL genes | Requires separate runs | Single integrated run |
| Identifies trait clusters | No | Yes |
| Accounts for shared architecture | Partially | Yes |
| Computational cost | Low | Moderate |

`HyPrColoc` (Hypothesis Prioritisation for multi-trait Colocalization) tests
all traits simultaneously, identifying *clusters* of traits that share a causal
variant. For this pipeline, the traits are: (1) AD GWAS, and (2–N) brain eQTL
signals for each candidate gene.

---

## How It Works

### Input
- **Beta matrix** (SNPs × Traits): effect sizes from GWAS and eQTL datasets
- **SE matrix** (SNPs × Traits): standard errors for each effect estimate
- One row per SNP, one column per trait (AD GWAS + one per eQTL gene)

### Algorithm
HyPrColoc uses a Bayesian framework with two key priors:

| Prior | Value Used | Meaning |
|-------|-----------|---------|
| `prior.1` | 1×10⁻⁴ | Probability any SNP is associated with one trait |
| `prior.c` | 0.02 | Conditional probability of sharing given association |

It returns for each cluster:
- **Posterior Probability (PP)**: probability that all traits in the cluster
  share a single causal variant
- **Regional Probability**: probability that the region contains any
  colocalization
- **Candidate SNP**: the most likely shared causal variant

### Output Interpretation

| PP Value | Evidence Strength | Interpretation |
|----------|-----------------|----------------|
| ≥ 0.90 | Strong | Very likely shared causal variant |
| ≥ 0.70 | Moderate | Good evidence — warrants follow-up |
| ≥ 0.50 | Suggestive | Possible sharing — treat cautiously |
| < 0.50 | Weak | Insufficient evidence |

The default threshold in this pipeline is **PP ≥ 0.70**.

---

## Pipeline Steps

1. **Fetch GWAS data** — same hybrid approach as main pipeline (server → backup)
2. **Fetch brain eQTLs** — GTEx brain tissue via Qtlizer (batched, r² ≥ 0.8)
3. **Build input matrices** — align SNPs across GWAS and eQTL datasets
4. **Quality control** — remove traits with >80% missing SNPs, impute missing
   values with null effects
5. **Run HyPrColoc** — multi-trait colocalization across all brain eQTL genes
6. **Filter results** — apply PP and regional probability thresholds
7. **Visualise** — posterior probability bar chart + regional vs posterior scatter
8. **Save outputs** — CSVs and plots to `work_dir`

---

## Important Notes on Summary Statistics

HyPrColoc requires **per-SNP beta and SE values** for every trait. The backup
mode in this pipeline uses placeholder values (beta = 0.1, SE = 0.05) when the
GWAS Catalog API is unavailable. These placeholders allow the pipeline to run
end-to-end for testing, but **will not produce meaningful colocalization
results**.

For publication-quality analysis, replace `gwas_summary_stats.csv` with real
full-locus summary statistics from:
- Bellenguez et al. 2022 (GCST90027158) — available from GWAS Catalog
- Kunkle et al. 2019 (GCST007320) — available from GWAS Catalog

Full summary stats can be downloaded as tab-separated files from:
https://www.ebi.ac.uk/gwas/studies/GCST90027158

---

## Output Files

| File | Description |
|------|-------------|
| `HyPrColoc_Full_Results.csv` | All colocalization clusters with PP values |
| `HyPrColoc_Significant_Hits.csv` | Filtered hits above PP threshold |
| `HyPrColoc_PosteriorProb_Plot.png` | Bar chart of PP per cluster |
| `HyPrColoc_Regional_vs_Posterior.png` | Scatter plot of both probabilities |
| `HyPrColoc_Summary.txt` | Run summary (SNP counts, hits, parameters) |
| `gwas_summary_stats.csv` | GWAS summary stats used as input |
| `qtls_brain_coloc.csv` | Brain eQTL data used as input |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pp4_threshold` | 0.7 | Minimum PP to report as significant hit |
| `regional_prob` | 0.7 | Minimum regional colocalization probability |
| `locus_window_kb` | 500 | Window size around lead SNP (kb) |
| `prior.1` | 1e-4 | Per-SNP association prior |
| `prior.c` | 0.02 | Conditional sharing prior |
| `qtl_tissue` | "Brain" | GTEx tissue filter |
| `ld_corr` | 0.8 | LD r² threshold for proxy SNPs |

---

## References

- Foley et al. (2021). A fast and efficient colocalization algorithm for
  identifying shared genetic risk factors across multiple traits.
  *Nature Communications*. https://doi.org/10.1038/s41467-020-20885-8

- Giambartolomei et al. (2014). Bayesian Test for Colocalisation between
  Pairs of Genetic Association Studies Using Summary Statistics.
  *PLOS Genetics*. https://doi.org/10.1371/journal.pgen.1004383

- GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory
  effects across human tissues. *Science*.
  https://doi.org/10.1126/science.aaz1776

- HyPrColoc GitHub: https://github.com/jrs95/hyprcoloc
