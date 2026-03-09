# Methods Documentation

Detailed descriptions of each analytical step in the Post-GWAS Alzheimer's Disease pipeline.

---

## Step 1: GWAS Data Retrieval

**Tool:** `gwasrapidd` R package  
**Source:** NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/)

The pipeline queries the GWAS Catalog REST API for SNP associations linked to a curated set of 19 established AD risk genes. Only associations meeting genome-wide significance (p ≤ 5×10⁻⁸) are retained. Risk allele variant IDs (rsIDs) are extracted from the `risk_alleles` slot of the returned `associations` S4 object.

**Fallback:** If the API is unreachable, a curated list of 46 validated AD SNPs from published meta-analyses (Lambert 2013, Jansen 2019) is loaded automatically.

---

## Step 2: SNP-to-Gene Mapping via eQTLs

**Tool:** `Qtlizer` R package  
**Method:** LD-based eQTL proxy mapping (r² ≥ 0.8)  
**Tissue filter:** Brain

Lead SNPs are expanded to LD blocks (r² ≥ 0.8) and mapped to gene targets via expression quantitative trait loci (eQTL) data. Qtlizer integrates multiple eQTL databases and returns tissue-specific gene associations.

Requests are batched (50 SNPs/batch) to prevent API timeout errors. Results are filtered to brain-relevant tissues using a case-insensitive regex match on the `tissue` column.

**Why eQTL mapping?**  
Most GWAS SNPs fall in non-coding regions. eQTL mapping links these regulatory variants to the genes whose expression they modulate, providing a biologically interpretable gene set.

---

## Step 3: Gene ID Conversion

**Tools:** `clusterProfiler::bitr`, `org.Hs.eg.db`

Gene symbols are cleaned to remove version suffixes (e.g., `GENE.1`) and antisense annotations (e.g., `GENE-AS1`). Cleaned symbols are mapped to Entrez IDs using the human annotation database `org.Hs.eg.db`. Unmapped genes are dropped (`drop = TRUE`).

Both the symbol and Entrez ID sets are carried forward — symbols for gProfiler2 and enrichR, Entrez IDs for ReactomePA, clusterProfiler, and DOSE.

---

## Step 4a: Reactome Pathway Enrichment

**Tool:** `ReactomePA::enrichPathway`  
**Database:** Reactome (https://reactome.org)  
**Correction:** Benjamini-Hochberg FDR, cutoff p < 0.05

Tests for over-representation of candidate genes in curated Reactome biological pathways. Results are visualised as a dot plot (top 15 pathways by gene ratio).

---

## Step 4b: Gene Ontology Enrichment (Biological Process)

**Tool:** `clusterProfiler::enrichGO`  
**Ontology:** Biological Process (BP)  
**Correction:** BH FDR, cutoff p < 0.05  
**Simplification:** `simplify()` at 0.7 semantic similarity cutoff

GO-BP terms are tested for over-representation. Because GO terms are highly redundant, `simplify()` removes semantically similar terms, retaining only the most significant representative in each cluster.

---

## Step 4c: gProfiler2 Multi-Database Enrichment

**Tool:** `gprofiler2::gost`  
**Databases:** GO (all branches), Reactome, KEGG  
**Correction:** FDR (g:SCS method)

gProfiler2 provides an independent cross-database enrichment with its own multiple testing correction (g:SCS), which is more conservative than standard FDR for gene set testing. The Manhattan-style `gostplot` visualises significance across all databases simultaneously.

---

## Step 4d: enrichR

**Tool:** `enrichR::enrichr`  
**Databases queried:**
- `GO_Biological_Process_2023`
- `Reactome_Pathways_2024`
- `DisGeNET`

enrichR queries the Enrichr web server with the gene symbol list and returns ranked enrichment results. DisGeNET results are particularly informative for cross-disease comparisons.

---

## Step 4e: Disease Ontology Enrichment (DOSE)

**Tool:** `DOSE::enrichDO`  
**Database:** Disease Ontology (DO)  
**Correction:** BH FDR, cutoff p < 0.05

Tests whether candidate genes are enriched in disease ontology terms beyond Alzheimer's disease itself. This is useful for identifying shared genetic architecture with related conditions (e.g., frontotemporal dementia, Parkinson's disease, vascular dementia).

---

## Statistical Notes

- All enrichment tests use a hypergeometric test (over-representation analysis, ORA)
- Background gene set: all annotated human genes in the respective database
- Multiple testing correction: Benjamini-Hochberg FDR throughout
- Significance threshold: adjusted p < 0.05

---

## Limitations

- ORA does not account for effect sizes or LD structure of GWAS hits
- eQTL mapping depends on available tissue coverage in Qtlizer's backend databases
- Functional enrichment is sensitive to gene set size; very small gene lists (<10 genes) may yield no significant results
- The fallback SNP/gene lists are static and may not reflect the latest GWAS findings
