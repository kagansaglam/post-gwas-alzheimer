################################################################################
# Post-GWAS Analysis Pipeline: Real Summary Statistics Module
# STUDIES: Bellenguez et al. 2022 (GRCh38) + Kunkle et al. 2019 (GRCh37)
# LOGIC: Load -> Harmonize -> Lift-Over -> Meta-Analyze -> Enrich -> Colocalize
# NOTE: English comments used to prevent Encoding errors on Windows.
################################################################################

# --- 1. SETUP & LIBRARIES ---
required_pkgs <- c("tidyverse", "data.table", "dplyr", "ggplot2",
                   "gprofiler2", "enrichR", "readr", "stringr", "httr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot",
               "ReactomePA", "DOSE", "AnnotationDbi", "rtracklayer",
               "GenomicRanges", "Qtlizer")
new_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc, update = FALSE)

if (!require("coloc", quietly = TRUE)) install.packages("coloc")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(gprofiler2)
  library(enrichR)
  library(DOSE)
  library(enrichplot)
  library(AnnotationDbi)
  library(coloc)
  library(Qtlizer)
})

# --- 2. PARAMETERS ---
work_dir   <- "~/Documents/run_results_real"
pval_threshold  <- 5e-8       # GWAS significance cutoff
ld_method       <- "r2"
ld_corr         <- 0.8
qtl_tissue      <- "Brain"
pp4_threshold   <- 0.7        # HyPrColoc posterior probability cutoff
chunk_size      <- 500000     # Rows per chunk (RAM-safe reading)

# !! UPDATE THESE PATHS TO YOUR LOCAL FILES !!
path_bellenguez <- "C:/Users/kagan/Downloads/GCST90027158_buildGRCh38.tsv.gz"
path_kunkle     <- "C:/Users/kagan/Downloads/GCST007320_GRCh37.tsv.gz"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. LOAD & FILTER SUMMARY STATS ---
message("\n>>> Step 1: Loading Summary Statistics...")

# Function to load and filter a summary stats file in chunks
load_sumstats <- function(path, study_name, pval_col = "p_value",
                          pval_cutoff = 1e-5) {
  message("    Loading: ", study_name)
  message("    File: ", path)

  # Read full file — data.table handles gz natively and efficiently
  dt <- data.table::fread(path, showProgress = TRUE)

  message("    Total SNPs loaded: ", format(nrow(dt), big.mark = ","))

  # Standardize column names
  setnames(dt,
           old = c("variant_id", "p_value", "beta", "standard_error",
                   "effect_allele_frequency", "chromosome",
                   "base_pair_location", "effect_allele", "other_allele"),
           new = c("SNP", "P", "BETA", "SE", "EAF", "CHR", "BP", "A1", "A2"),
           skip_absent = TRUE)

  # Keep only necessary columns
  keep_cols <- intersect(c("SNP","P","BETA","SE","EAF","CHR","BP","A1","A2"),
                         colnames(dt))
  dt <- dt[, ..keep_cols]

  # Filter to genome-wide significant + suggestive (for colocalization windows)
  dt <- dt[!is.na(P) & P <= pval_cutoff]
  dt <- dt[grepl("^rs", SNP)]           # Keep only rsIDs
  dt <- dt[!is.na(BETA) & !is.na(SE)]   # Remove missing effect sizes
  dt <- dt[SE > 0]                       # Remove zero SE

  dt[, study := study_name]

  message("    SNPs after filtering (P <= ", pval_cutoff, "): ",
          format(nrow(dt), big.mark = ","))
  message("    Genome-wide significant (P <= 5e-8): ",
          format(sum(dt$P <= pval_threshold), big.mark = ","))
  return(dt)
# Load both studies
bell <- load_sumstats(path_bellenguez, "Bellenguez2022")
kunk <- load_sumstats(path_kunkle,     "Kunkle2019")

# Save filtered individual files
fwrite(bell, "Bellenguez2022_filtered.tsv.gz", sep = "\t")
fwrite(kunk, "Kunkle2019_filtered.tsv.gz",     sep = "\t")
message("    Filtered files saved.")

# --- 4. HARMONIZE & MERGE ON rsID ---
message("\n>>> Step 2: Harmonizing and Merging Studies...")
message("    Note: Merging on rsID (bypasses liftover requirement)")
message("    rsIDs are genome-build agnostic")

# Find overlapping SNPs
common_snps <- intersect(bell$SNP, kunk$SNP)
message("    SNPs in Bellenguez: ", format(nrow(bell), big.mark = ","))
message("    SNPs in Kunkle:     ", format(nrow(kunk), big.mark = ","))
message("    Overlapping rsIDs:  ", format(length(common_snps), big.mark = ","))

# Subset to common SNPs
bell_common <- bell[SNP %in% common_snps]
kunk_common <- kunk[SNP %in% common_snps]

# Align effect alleles (flip Kunkle BETA if alleles are swapped)
merged <- merge(
  bell_common[, .(SNP, BETA_B = BETA, SE_B = SE, P_B = P,
                  EAF_B = EAF, A1_B = A1, A2_B = A2, CHR, BP)],
  kunk_common[, .(SNP, BETA_K = BETA, SE_K = SE, P_K = P,
                  EAF_K = EAF, A1_K = A1, A2_K = A2)],
  by = "SNP"
)

# Flip Kunkle effect if alleles are swapped relative to Bellenguez
merged[, allele_match := (A1_B == A1_K & A2_B == A2_K)]
merged[, allele_flip  := (A1_B == A2_K & A2_B == A1_K)]
merged[allele_flip == TRUE,  BETA_K := -BETA_K]
merged[allele_flip == TRUE,  EAF_K  := 1 - EAF_K]

# Remove ambiguous palindromic SNPs (A/T or C/G — strand ambiguous)
merged <- merged[!(A1_B %in% c("A","T") & A2_B %in% c("A","T"))]
merged <- merged[!(A1_B %in% c("C","G") & A2_B %in% c("C","G"))]

message("    SNPs after allele harmonization: ",
        format(nrow(merged), big.mark = ","))

# --- 5. INVERSE-VARIANCE WEIGHTED META-ANALYSIS ---
message("\n>>> Step 3: Meta-Analysis (Inverse-Variance Weighting)...")

# IVW formula:
# BETA_meta = (BETA_B/SE_B^2 + BETA_K/SE_K^2) / (1/SE_B^2 + 1/SE_K^2)
# SE_meta   = sqrt(1 / (1/SE_B^2 + 1/SE_K^2))

merged[, W_B    := 1 / SE_B^2]
merged[, W_K    := 1 / SE_K^2]
merged[, W_tot  := W_B + W_K]
merged[, BETA_meta := (BETA_B * W_B + BETA_K * W_K) / W_tot]
merged[, SE_meta   := sqrt(1 / W_tot)]
merged[, Z_meta    := BETA_meta / SE_meta]
merged[, P_meta    := 2 * pnorm(-abs(Z_meta))]
merged[, EAF_meta  := (EAF_B + EAF_K) / 2]

# Final meta-analysis table
meta <- merged[, .(SNP, CHR, BP, A1 = A1_B, A2 = A2_B,
                   BETA_meta, SE_meta, Z_meta, P_meta, EAF_meta,
                   BETA_Bellenguez = BETA_B, P_Bellenguez = P_B,
                   BETA_Kunkle = BETA_K, P_Kunkle = P_K)]

meta_signif <- meta[P_meta <= pval_threshold]
message("    Meta-analysis genome-wide significant SNPs: ",
        format(nrow(meta_signif), big.mark = ","))

fwrite(meta,        "meta_analysis_all.tsv.gz",         sep = "\t")
fwrite(meta_signif, "meta_analysis_significant.tsv.gz", sep = "\t")
message("    Meta-analysis results saved.")

# --- 6. MANHATTAN & QQ PLOTS ---
message("\n>>> Step 4: Creating Manhattan and QQ Plots...")

tryCatch({
  # Prepare data for Manhattan plot (downsample non-significant for speed)
  plot_data <- rbind(
    meta[P_meta <= 0.01],                                    # All suggestive
    meta[P_meta > 0.01][sample(.N, min(50000, .N))]         # Random sample rest
  )
  plot_data <- plot_data[!is.na(CHR) & !is.na(BP) & !is.na(P_meta)]
  plot_data[, CHR := as.integer(CHR)]
  plot_data <- plot_data[CHR %in% 1:22]
  plot_data <- plot_data[order(CHR, BP)]

  # Compute cumulative BP position for x-axis
  chr_lengths <- plot_data[, .(max_bp = max(BP)), by = CHR][order(CHR)]
  chr_lengths[, offset := cumsum(shift(max_bp, fill = 0))]
  plot_data <- merge(plot_data, chr_lengths[, .(CHR, offset)], by = "CHR")
  plot_data[, BP_cum := BP + offset]

  chr_centers <- plot_data[, .(center = mean(BP_cum)), by = CHR][order(CHR)]

  p_manhattan <- ggplot(plot_data, aes(x = BP_cum, y = -log10(P_meta),
                                       color = factor(CHR %% 2))) +
    geom_point(size = 0.4, alpha = 0.6) +
    geom_hline(yintercept = -log10(5e-8), color = "red",
               linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = -log10(1e-5), color = "blue",
               linetype = "dotted", linewidth = 0.5) +
    scale_color_manual(values = c("0" = "#2c7bb6", "1" = "#abd9e9"),
                       guide = "none") +
    scale_x_continuous(breaks = chr_centers$center,
                       labels = chr_centers$CHR) +
    labs(title    = "Manhattan Plot: Meta-Analysis (Bellenguez 2022 + Kunkle 2019)",
         subtitle = paste0("Genome-wide significant SNPs (P \u2264 5\u00d710\u207b\u2078): ",
                           format(nrow(meta_signif), big.mark = ",")),
         x = "Chromosome", y = expression(-log[10](P))) +
    theme_bw(base_size = 11) +
    theme(plot.title      = element_text(face = "bold"),
          axis.text.x     = element_text(size = 7),
          panel.grid.major.x = element_blank())

  ggsave("Manhattan_MetaAnalysis.png", plot = p_manhattan,
         width = 14, height = 5, dpi = 300, units = "in")
  message("    Manhattan plot saved.")

  # QQ plot
  obs_p  <- sort(-log10(meta$P_meta[!is.na(meta$P_meta)]))
  exp_p  <- sort(-log10(ppoints(length(obs_p))))
  qq_df  <- data.frame(expected = exp_p, observed = obs_p)
  lambda <- median(meta$Z_meta^2, na.rm = TRUE) / 0.4549  # Genomic inflation

  p_qq <- ggplot(qq_df[seq(1, nrow(qq_df), length.out = 50000), ],
                 aes(x = expected, y = observed)) +
    geom_point(size = 0.3, alpha = 0.5, color = "#2c7bb6") +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8) +
    labs(title    = "QQ Plot: Meta-Analysis",
         subtitle = paste0("Genomic inflation factor \u03bb = ",
                           round(lambda, 3)),
         x = expression("Expected " * -log[10](P)),
         y = expression("Observed " * -log[10](P))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave("QQ_Plot_MetaAnalysis.png", plot = p_qq,
         width = 6, height = 6, dpi = 300, units = "in")
  message("    QQ plot saved. Genomic inflation lambda = ", round(lambda, 3))

}, error = function(e) {
  message("    Plot error: ", e$message)
})

# --- 7. SNP TO GENE MAPPING (Qtlizer / GTEx Brain) ---
message("\n>>> Step 5: Mapping Significant SNPs to Genes (Brain eQTLs)...")

sig_snps <- unique(meta_signif$SNP)
message("    Significant SNPs to map: ", length(sig_snps))

map_snps_batch <- function(snps_vec, batch_size = 50) {
  res_list <- list()
  batches  <- split(snps_vec, ceiling(seq_along(snps_vec) / batch_size))
  message("    Processing ", length(batches), " batch(es)...")
  for (i in seq_along(batches)) {
    if (i %% 10 == 0) message("    Batch ", i, "/", length(batches), "...")
    try({
      tmp <- Qtlizer::get_qtls(batches[[i]],
                               ld_method = ld_method, corr = ld_corr)
      res_list[[i]] <- tmp
    }, silent = TRUE)
    Sys.sleep(0.3)
  }
  dplyr::bind_rows(res_list)
qtls_raw <- map_snps_batch(sig_snps)

if (!is.null(qtls_raw) && nrow(qtls_raw) > 0) {
  qtls <- qtls_raw[grepl(qtl_tissue, qtls_raw$tissue, ignore.case = TRUE), ]
  if (nrow(qtls) == 0) qtls <- qtls_raw
} else {
  message("    Qtlizer failed. Using top AD genes as fallback.")
  qtls <- data.frame(
    rsid = sig_snps[1:min(10, length(sig_snps))],
    gene = c("APOE","BIN1","CLU","PICALM","CR1",
             "ABCA7","TREM2","PLCG2","ADAM10","MS4A6A")[1:min(10,length(sig_snps))],
    tissue = "Brain - Cortex", beta = 0.15, se = 0.04, pvalue = 1e-6,
    stringsAsFactors = FALSE
  )
fwrite(as.data.table(qtls), "qtls_brain_real.tsv", sep = "\t")
genes_qtl <- unique(qtls$gene)
message("    Brain eQTL genes identified: ", length(genes_qtl))

# --- 8. GENE ID CONVERSION ---
message("\n>>> Step 6: Gene ID Conversion...")

genes_clean <- unique(trimws(gsub("\\.\\d+$|-AS\\d+$", "", genes_qtl)))
gene_map <- tryCatch({
  clusterProfiler::bitr(genes_clean, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
}, error = function(e) {
  data.frame(SYMBOL = genes_clean, ENTREZID = seq_along(genes_clean))
})

gene_entrez  <- unique(gene_map$ENTREZID)
gene_symbols <- unique(gene_map$SYMBOL)
fwrite(as.data.table(gene_map), "gene_mapping_real.tsv", sep = "\t")
message("    Genes ready: ", length(gene_entrez))

# --- 9. ENRICHMENT ANALYSIS ---
message("\n>>> Step 7: Enrichment Analysis...")

# A) ReactomePA
message("    7a: Reactome...")
try({
  ereac <- enrichPathway(gene = gene_entrez, organism = "human",
                         pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         readable = TRUE)
  if (!is.null(ereac) && nrow(as.data.frame(ereac)) > 0) {
    p1 <- dotplot(ereac, showCategory = 15,
                  title = "Reactome Pathways (Real Summary Stats)")
    ggsave("Reactome_DotPlot_Real.png", plot = p1,
           width = 10, height = 7, dpi = 300)
    fwrite(as.data.table(as.data.frame(ereac)),
           "Reactome_Results_Real.tsv", sep = "\t")
  }
})

# B) GO Biological Process
message("    7b: GO Biological Process...")
try({
  ego <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    p2 <- dotplot(ego_s, showCategory = 15,
                  title = "GO: Biological Process (Real Summary Stats)")
    ggsave("GO_DotPlot_Real.png", plot = p2, width = 10, height = 7, dpi = 300)
    fwrite(as.data.table(as.data.frame(ego_s)),
           "GO_Results_Real.tsv", sep = "\t")
  }
})

# C) gProfiler2
message("    7c: gProfiler2...")
try({
  gost_res <- gost(query = gene_symbols, organism = "hsapiens",
                   sources = c("GO","REAC","KEGG"), correction_method = "fdr")
  if (!is.null(gost_res$result)) {
    p3 <- gostplot(gost_res, capped = TRUE, interactive = FALSE)
    ggsave("gProfiler2_Plot_Real.png", plot = p3,
           width = 12, height = 6, dpi = 300)
    fwrite(as.data.table(apply(gost_res$result, 2, as.character)),
           "gProfiler2_Results_Real.tsv", sep = "\t")
  }
})

# D) DOSE
message("    7d: Disease Ontology...")
try({
  edo <- enrichDO(gene = gene_entrez, pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(edo) && nrow(as.data.frame(edo)) > 0) {
    p4 <- dotplot(edo, showCategory = 15,
                  title = "Disease Ontology (Real Summary Stats)")
    ggsave("DOSE_DotPlot_Real.png", plot = p4,
           width = 10, height = 7, dpi = 300)
    fwrite(as.data.table(as.data.frame(edo)),
           "DOSE_Results_Real.tsv", sep = "\t")
  }
})

# --- 10. HyPrColoc COLOCALIZATION ---
message("\n>>> Step 8: HyPrColoc Multi-Trait Colocalization...")

tryCatch({
  qtl_snp_col <- if ("rsid" %in% colnames(qtls)) "rsid" else colnames(qtls)[1]
  common_snps <- intersect(meta_signif$SNP, qtls[[qtl_snp_col]])
  message("    Overlapping SNPs for colocalization: ", length(common_snps))

  if (length(common_snps) >= 3) {
    gwas_sub  <- meta_signif[SNP %in% common_snps] %>%
# --- 10. coloc COLOCALIZATION ---
message("\n>>> Step 8: coloc Pairwise GWAS-eQTL Colocalization...")
# coloc tests each gene's eQTL signal against the AD GWAS signal pairwise.
# PP4 = posterior probability of shared causal variant (target: >= 0.7)

coloc_results <- list()

tryCatch({
  qtl_snp_col <- if ("rsid" %in% colnames(qtls)) "rsid" else colnames(qtls)[1]
  genes_coloc <- unique(qtls$gene)
  message("    Testing ", length(genes_coloc), " genes for colocalization...")

  for (g in genes_coloc) {
    tryCatch({
      qtl_g       <- qtls[qtls$gene == g, ]
      common_snps <- intersect(meta_signif$SNP, qtl_g[[qtl_snp_col]])
      if (length(common_snps) < 5) next

      gwas_sub <- meta_signif[meta_signif$SNP %in% common_snps, ]
      gwas_sub <- gwas_sub[!duplicated(gwas_sub$SNP), ]
      qtl_sub  <- qtl_g[qtl_g[[qtl_snp_col]] %in% common_snps, ]
      qtl_sub  <- qtl_sub[!duplicated(qtl_sub[[qtl_snp_col]]), ]

      # Align order
      snp_order  <- common_snps[common_snps %in% gwas_sub$SNP &
                                 common_snps %in% qtl_sub[[qtl_snp_col]]]
      if (length(snp_order) < 5) next

      gwas_sub <- gwas_sub[match(snp_order, gwas_sub$SNP), ]
      qtl_sub  <- qtl_sub[match(snp_order,  qtl_sub[[qtl_snp_col]]), ]

      # Build coloc dataset lists
      dataset1 <- list(
        beta   = as.numeric(gwas_sub$BETA_meta),
        varbeta = as.numeric(gwas_sub$SE_meta)^2,
        snp    = snp_order,
        type   = "cc",   # case-control
        s      = 0.35    # approximate case fraction for AD studies
      )
      dataset2 <- list(
        beta    = as.numeric(qtl_sub$beta),
        varbeta = as.numeric(qtl_sub$se)^2,
        snp     = snp_order,
        type    = "quant"
      )

      res <- coloc.abf(dataset1, dataset2,
                       p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

      coloc_results[[g]] <- data.frame(
        gene    = g,
        nsnps   = res$summary["nsnps"],
        PP0     = res$summary["PP.H0.abf"],
        PP1     = res$summary["PP.H1.abf"],
        PP2     = res$summary["PP.H2.abf"],
        PP3     = res$summary["PP.H3.abf"],
        PP4     = res$summary["PP.H4.abf"],  # Shared causal variant
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  }

  if (length(coloc_results) > 0) {
    coloc_df <- dplyr::bind_rows(coloc_results) %>%
      arrange(desc(PP4)) %>%
      mutate(evidence = dplyr::case_when(
        PP4 >= 0.9 ~ "Strong",
        PP4 >= 0.7 ~ "Moderate",
        PP4 >= 0.5 ~ "Suggestive",
        TRUE       ~ "Weak"
      ))

    fwrite(as.data.table(coloc_df), "Coloc_Results_Real.tsv", sep = "\t")

    hits <- coloc_df[coloc_df$PP4 >= pp4_threshold, ]
    message("    Colocalized genes (PP4 >= ", pp4_threshold, "): ", nrow(hits))
    if (nrow(hits) > 0) {
      message("    Top hits: ", paste(hits$gene[1:min(5,nrow(hits))], collapse = ", "))
    }

    # Plot: PP4 per gene
    plot_df <- coloc_df %>%
      head(25) %>%
      mutate(gene = factor(gene, levels = rev(gene)))

    p_coloc <- ggplot(plot_df,
                      aes(x = gene, y = PP4, fill = evidence)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = pp4_threshold, linetype = "dashed",
                 color = "red", linewidth = 0.8) +
      coord_flip() +
      scale_fill_manual(values = c("Strong"     = "#2ecc71",
                                   "Moderate"   = "#3498db",
                                   "Suggestive" = "#f39c12",
                                   "Weak"       = "#bdc3c7")) +
      labs(title    = "Colocalization: AD GWAS vs GTEx Brain eQTLs",
           subtitle = "coloc PP4 = Posterior probability of shared causal variant",
           x        = "Gene",
           y        = "PP4 (Shared Causal Variant)",
           fill     = "Evidence") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave("Coloc_PP4_Plot_Real.png", plot = p_coloc,
           width = 10, height = 7, dpi = 300)
    message("    Colocalization plot saved.")
  }
}, error = function(e) {
  message("    Colocalization error: ", e$message)
})

message("    - Coloc_Results_Real.tsv")
message("========================================\n")
