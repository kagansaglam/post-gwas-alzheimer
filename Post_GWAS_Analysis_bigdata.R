################################################################################
# Post-GWAS Analysis Pipeline: SNPs > GENES > PATHWAYS
# METHOD: DIRECT BIG DATA (Targeting Massive Studies)
# TARGETS: Bellenguez et al. (2022) & Kunkle et al. (2019)
################################################################################

# --- 1. SETUP & LIBRARIES ---
# Standard packages
required_pkgs <- c("tidyverse", "gwasrapidd", "httr", "jsonlite", "dplyr", 
                   "gprofiler2", "enrichR", "readr", "stringr", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("Qtlizer", "clusterProfiler", "org.Hs.eg.db", "enrichplot", 
               "ReactomePA", "DOSE", "AnnotationDbi")
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc, update = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(gwasrapidd)
  library(Qtlizer)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(clusterProfiler)
  library(gprofiler2)
  library(enrichR)
  library(DOSE)
  library(enrichplot)
  library(AnnotationDbi)
})

# --- 2. PARAMETERS ---
work_dir       <- "~/Documents/run_results"
trait          <- "Alzheimer_BigData_Attempt"
pval_threshold <- 5e-8
ld_method      <- "r2"
ld_corr        <- 0.8
qtl_tissue     <- "Brain"  # Tissue Filter

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. FETCH GWAS DATA (BIG DATA MODE) ---
message("\n>>> Step 1: Fetching MASSIVE Study Data...")
message("    Targets: GCST90027158 (Bellenguez 2022) & GCST007320 (Kunkle 2019)")

# These are the two largest Alzheimer's GWAS meta-analyses.
# Attempting to fetch ALL associations from them.
target_studies <- c("GCST90027158", "GCST007320")

snp_list <- NULL

tryCatch({
  assoc <- get_associations(study_id = target_studies)
  
  if (nrow(assoc@associations) > 0) {
    message("    Raw Data Rows: ", nrow(assoc@associations))
    
    # Filter by P-Value
    signif_indices <- which(assoc@associations$pvalue <= pval_threshold)
    lead_assoc_ids <- unique(assoc@associations$association_id[signif_indices])
    
    if(length(lead_assoc_ids) > 0) {
      risk_alleles   <- assoc@risk_alleles
      snp_rows       <- which(risk_alleles$association_id %in% lead_assoc_ids)
      raw_snps       <- unique(risk_alleles$variant_id[snp_rows])
      snp_list       <- grep("^rs", raw_snps, value = TRUE) # Filter RSIDs
      snp_list       <- unique(snp_list)
      
      message("    SUCCESS: Fetched ", length(snp_list), " SNPs from Big Data studies.")
    }
  }
}, error = function(e) {
  message("!!! SERVER TIMEOUT: The dataset is too large for the current connection.")
  message("    Error Details: ", e$message)
})

# CRITICAL CHECK: If Big Data failed, stop here (or use backup if you prefer)
if (is.null(snp_list) || length(snp_list) == 0) {
  stop("Big Data fetch failed (500 Error). Please revert to the 'Hybrid' code provided earlier.")
}

# Save Huge List
write.csv(data.frame(SNP = snp_list), "lead_snps_BIGDATA.csv", row.names = FALSE)


# --- 4. SNP TO GENE MAPPING (BATCH MODE) ---
message("\n>>> Step 2: Mapping SNPs to Genes (Qtlizer)")
message("    Target Tissue: ", qtl_tissue)
message("    Note: This will take time due to large SNP count...")

# Important: Batch size reduced to 100 to prevent Qtlizer crash with Big Data
map_snps_batch <- function(snps_vec, batch_size = 100) {
  res_list <- list()
  batches <- split(snps_vec, ceiling(seq_along(snps_vec)/batch_size))
  
  message("    Processing ", length(batches), " batches...")
  
  for(i in seq_along(batches)) {
    if(i %% 5 == 0) message("    Batch ", i, "/", length(batches), "...")
    try({
      tmp <- Qtlizer::get_qtls(batches[[i]], ld_method = ld_method, corr = ld_corr)
      res_list[[i]] <- tmp
    }, silent = TRUE)
    Sys.sleep(0.2) # Small pause
  }
  bind_rows(res_list)
}

qtls_raw <- map_snps_batch(snp_list)

if (is.null(qtls_raw) || nrow(qtls_raw) == 0) stop("Qtlizer returned no results.")

# Apply Tissue Filter
if (qtl_tissue != "all") {
  qtls <- qtls_raw[grepl(qtl_tissue, qtls_raw$tissue, ignore.case = TRUE), ]
} else {
  qtls <- qtls_raw
}

# Fallback if filter is too strict
if (nrow(qtls) == 0) {
  message("! Warning: No Brain QTLs found. Using all tissues.")
  qtls <- qtls_raw
}

write.csv(qtls, "qtls_filtered_BIGDATA.csv", row.names = FALSE)

genes_qtl <- unique(qtls$gene)
message("    Genes Identified: ", length(genes_qtl))


# --- 5. GENE ID CONVERSION ---
message("\n>>> Step 3: Gene ID Conversion")

genes_clean <- unique(trimws(gsub("\\.\\d+$|-AS\\d+$", "", genes_qtl)))

gene_map <- tryCatch({
  bitr(genes_clean, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
}, error = function(e) {
  return(data.frame(SYMBOL=genes_clean, ENTREZID=seq_along(genes_clean)))
})

gene_entrez  <- unique(gene_map$ENTREZID)
gene_symbols <- unique(gene_map$SYMBOL)

message("    Genes ready for analysis: ", length(gene_entrez))
write.csv(gene_map, "gene_mapping.csv", row.names = FALSE)


# --- 6. ENRICHMENT ANALYSIS ---

# A) ReactomePA
message("\n>>> Step 4a: Reactome Analysis")
try({
  ereac <- enrichPathway(gene = gene_entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ereac) && nrow(as.data.frame(ereac)) > 0) {
    p1 <- dotplot(ereac, showCategory = 15, title = "Reactome Pathways")
    print(p1)
    write.csv(as.data.frame(ereac), "Reactome_Results.csv")
  } else { message("    No significant Reactome pathways.") }
})

# B) clusterProfiler (GO)
message("\n>>> Step 4b: GO (Biological Process) Analysis")
try({
  ego <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_simple <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    p2 <- dotplot(ego_simple, showCategory = 15, title = "GO: Biological Process")
    print(p2)
    write.csv(as.data.frame(ego_simple), "GO_Results.csv")
  }
})

# C) gProfiler2
message("\n>>> Step 4c: gProfiler2 Analysis")
try({
  gost_res <- gost(query = gene_symbols, organism = "hsapiens", sources = c("GO", "REAC", "KEGG"), correction_method = "fdr")
  if (!is.null(gost_res$result)) {
    p_gost <- gostplot(gost_res, capped = TRUE, interactive = FALSE)
    print(p_gost)
    clean_gost <- apply(gost_res$result, 2, as.character)
    write.csv(clean_gost, "gProfiler2_Results.csv", row.names = FALSE)
  }
})

# D) enrichR
message("\n>>> Step 4d: enrichR Analysis")
if(length(gene_symbols) > 0) {
  dbs <- c("GO_Biological_Process_2023", "Reactome_Pathways_2024", "DisGeNET")
  try({
    enrichr_res <- enrichr(gene_symbols, dbs)
    if(length(enrichr_res) > 0) write.csv(enrichr_res[[1]], "EnrichR_Results_TopDB.csv")
  }, silent = TRUE)
}

# E) DOSE
message("\n>>> Step 4e: DOSE Analysis")
try({
  edo <- enrichDO(gene = gene_entrez, pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(edo) && nrow(as.data.frame(edo)) > 0) {
    p4 <- dotplot(edo, showCategory = 15, title = "Disease Ontology")
    print(p4)
    write.csv(as.data.frame(edo), "DOSE_Results.csv")
  }
})

message("\n>>> BIG DATA PIPELINE COMPLETE. Check results in: ", work_dir)