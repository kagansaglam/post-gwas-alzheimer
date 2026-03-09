# METHODS: Drug Target Prioritization (OpenTargets Platform)

## Script
`Post_GWAS_DrugTarget_OpenTargets.R`

## Purpose
Systematically identify which Alzheimer's Disease GWAS candidate genes are
druggable, and whether existing approved or clinical-stage compounds already
target them — directly actionable intelligence for drug discovery.

---

## Why Drug Target Prioritization?

GWAS identifies genetic risk factors, but not all risk genes are equally useful
for drug development. This module answers three critical questions:

1. **Is the gene druggable?** Does it encode a protein that can be targeted by
   a small molecule or biologic?
2. **Does a drug already exist?** Are there approved compounds or clinical
   candidates targeting this gene?
3. **How strongly is it linked to AD?** Does OpenTargets genetic, expression,
   and functional evidence support this gene as a disease driver?

---

## Data Source: OpenTargets Platform

OpenTargets integrates evidence from multiple sources to score gene-disease
associations and map drug-target relationships:

| Evidence Type | Source |
|--------------|--------|
| Genetic associations | GWAS Catalog, ClinVar, UniProt |
| Somatic mutations | COSMIC, intOGen |
| Gene expression | Expression Atlas |
| Pathway evidence | Reactome, SLAPenrich |
| Animal models | PhenoDigm |
| Text mining | EuropePMC |
| Drug-target links | ChEMBL, FDA approvals |

Access method: **GraphQL API** (`api.platform.opentargets.org/api/v4/graphql`)

---

## Algorithm

### Step 1: Gene Symbol → Ensembl ID Resolution
Each candidate gene symbol is resolved to its Ensembl ID via the OpenTargets
search endpoint — required for all subsequent API queries.

### Step 2: AD Association Score Retrieval
For each gene, the overall OpenTargets association score with Alzheimer's
Disease (`MONDO_0004975`) is retrieved. This score (0–1) integrates all
evidence types listed above.

### Step 3: Drug-Gene Interaction Retrieval
Known drugs targeting each gene are fetched, filtered to **Phase II and above**
(clinical-stage compounds only):

| Phase | Included? |
|-------|-----------|
| Pre-clinical | ❌ No |
| Phase I | ❌ No |
| Phase II | ✅ Yes |
| Phase III | ✅ Yes |
| Approved | ✅ Yes |

### Step 4: Priority Scoring
A composite priority score combines both dimensions:

```
Priority Score = (AD Association Score × 0.6) + (Max Drug Phase / 4 × 0.4)
```

This weights biological relevance (AD score) more than druggability, but
ensures genes with advanced clinical compounds rank highly even if their
overall AD score is moderate.

---

## Output Files

| File | Description |
|------|-------------|
| `DrugTarget_Prioritized.csv` | All genes ranked by priority score |
| `DrugTarget_Top20.csv` | Top 20 druggable targets |
| `DrugTarget_AllDrugs.csv` | Full drug-gene interaction table |
| `DrugTarget_ADScores.csv` | OpenTargets AD association scores |
| `DrugTarget_PriorityPlot.png` | Bar chart: genes ranked by priority |
| `DrugTarget_PhasePlot.png` | Drug pipeline phase distribution |
| `DrugTarget_BubblePlot.png` | AD score vs drug phase bubble chart |
| `DrugTarget_Summary.txt` | Plain-text summary report |

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_drug_phase` | 2 | Minimum clinical phase to include |
| `ad_disease_id` | `MONDO_0004975` | Alzheimer's Disease ontology ID |
| `assoc_score_threshold` | 0.1 | Minimum OpenTargets AD association score |

---

## Interpreting Results

**Priority Score > 0.7:** Strong AD evidence + advanced clinical compound —
highest priority for follow-up. Example: APOE, TREM2, ADAM10.

**Priority Score 0.4–0.7:** Good AD evidence but earlier-stage drugs — worth
monitoring for clinical development.

**Approved Drug + High AD Score:** Most actionable — drug repurposing
candidate. Existing safety/PK data dramatically reduces development costs.

---

## Notable AD GWAS Genes in Drug Pipelines

| Gene | Role | Drug Status |
|------|------|-------------|
| `TREM2` | Microglial receptor, neuroinflammation | Phase II antibodies (AL002) |
| `ADAM10` | APP processing, Aβ production | Phase II/III modulators |
| `PLCG2` | Microglial signaling | Emerging target, early stage |
| `APOE` | Lipid transport, Aβ clearance | Multiple Phase II/III trials |
| `ACE` | Renin-angiotensin system | FDA approved (hypertension) |

---

## References

- Ochoa et al. (2021). Open Targets Platform: supporting systematic drug-target
  identification and prioritisation. *Nucleic Acids Research*.
  https://doi.org/10.1093/nar/gkaa1027

- OpenTargets Platform: https://platform.opentargets.org

- Finan et al. (2017). The druggable genome and support for target
  identification and validation in drug development.
  *Science Translational Medicine*.
  https://doi.org/10.1126/scitranslmed.aag1166
