################################################################################
# Post-GWAS Analysis Pipeline: Colocalization Module
# METHOD: HyPrColoc (Multi-Trait Colocalization)
# LOGIC: Try Server -> If Fails, Use Backup List -> Complete Analysis
# NOTE: English comments used to prevent Encoding errors on Windows.
# REQUIRES: Run Post_GWAS_Analysis_new.R first (or provide lead_snps.csv)
################################################################################

# --- 1. SETUP & LIBRARIES ---
required_pkgs <- c("tidyverse", "gwasrapidd", "httr", "jsonlite",
                   "readr", "stringr", "ggplot2", "dplyr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

# HyPrColoc (GitHub install)
if (!require("hyprcoloc", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("jrs95/hyprcoloc", build_vignettes = FALSE)
}

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("Qtlizer", "GenomicRanges", "IRanges")
new_bioc  <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc, update = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(gwasrapidd)
  library(hyprcoloc)
  library(Qtlizer)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# --- 2. PARAMETERS ---
work_dir       <- "~/Documents/run_results"
trait          <- "Alzheimer_HyPrColoc"
pval_threshold <- 5e-8
ld_method      <- "r2"
ld_corr        <- 0.8
qtl_tissue     <- "Brain"

# HyPrColoc thresholds
# pp4_threshold: posterior probability of shared causal variant (>=0.7 is stringent)
# regional_prob: minimum regional colocalization probability
pp4_threshold   <- 0.7
regional_prob   <- 0.7

# Locus window: +/- kb around each lead SNP for regional analysis
locus_window_kb <- 500

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. FETCH GWAS SUMMARY STATISTICS (HYBRID MODE) ---
message("\n>>> Step 1: Fetching GWAS Data (Hybrid Mode)...")

snp_list     <- NULL
gwas_assoc   <- NULL

# METHOD A: Try GWAS Catalog API
tryCatch({
  message("    Attempting to connect to GWAS Catalog server...")
  target_genes <- c("APOE", "BIN1", "CLU", "ABCA7", "CR1", "PICALM", "MS4A6A",
                    "CD33", "MS4A4E", "CD2AP", "EPHA1", "INPP5D", "MEF2C", "TREM2",
                    "SORL1", "PLCG2", "ADAM10", "ACE", "TOMM40")

  assoc <- get_associations(gene_name = target_genes)

  if (nrow(assoc@associations) > 0) {
    signif_idx   <- which(assoc@associations$pvalue <= pval_threshold)
    lead_ids     <- unique(assoc@associations$association_id[signif_idx])
    risk_alleles <- assoc@risk_alleles
    snp_rows     <- which(risk_alleles$association_id %in% lead_ids)
    raw_snps     <- unique(risk_alleles$variant_id[snp_rows])
    snp_list     <- grep("^rs", raw_snps, value = TRUE)

    # Also retain beta/se/pvalue for colocalization input
    gwas_assoc <- assoc@associations[signif_idx, ] %>%
      dplyr::select(association_id, pvalue, beta_number, standard_error,
                    effect_allele_frequency) %>%
      dplyr::left_join(
        risk_alleles %>% dplyr::select(association_id, variant_id),
        by = "association_id"
      ) %>%
      dplyr::filter(grepl("^rs", variant_id)) %>%
      dplyr::rename(SNP = variant_id, P = pvalue, BETA = beta_number,
                    SE = standard_error, EAF = effect_allele_frequency)

    message("    SUCCESS: Fetched ", length(snp_list), " SNPs from server.")
  }
}, error = function(e) {
  message("!!! SERVER ERROR: ", e$message)
})

# METHOD B: Backup SNP list if server failed
if (is.null(snp_list) || length(snp_list) == 0) {
  message("\n>>> PLAN B ACTIVATED: Using backup SNP list...")
  snp_list <- c(
    "rs429358",  "rs7412",    "rs744373",  "rs6733839", "rs11136000",
    "rs2279590", "rs3865444", "rs3818361", "rs6656401", "rs3851179",
    "rs610932",  "rs983392",  "rs670139",  "rs10933431","rs11771145",
    "rs10498633","rs75932628","rs72824905","rs11591147","rs17125721",
    "rs1476679", "rs7561528", "rs4844610", "rs190982",  "rs2075650",
    "rs11129640","rs423668",  "rs2459732", "rs7799015", "rs7920721",
    "rs9271058", "rs9331896", "rs593742",  "rs3764650", "rs12459419",
    "rs138137383","rs12972156","rs965471",  "rs556603",  "rs2093760",
    "rs9343759", "rs11218343","rs10838725","rs4147929", "rs6859",
    "rs7274581"
  )

  # Build minimal GWAS summary stats for backup SNPs
  # Note: Using placeholder SE=0.05 since exact values unavailable offline
  # Replace with real summary stats for publication-quality analysis
  gwas_assoc <- data.frame(
    SNP  = snp_list,
    P    = rep(1e-10, length(snp_list)),   # Placeholder: all below threshold
    BETA = rep(0.1,   length(snp_list)),   # Placeholder beta
    SE   = rep(0.05,  length(snp_list)),   # Placeholder SE
    EAF  = rep(0.3,   length(snp_list)),   # Placeholder allele frequency
    stringsAsFactors = FALSE
  )

  message("    Backup list loaded: ", length(snp_list), " SNPs.")
  message("    WARNING: Placeholder GWAS stats used. For real colocalization,")
  message("    replace with actual summary statistics (beta, SE, EAF per SNP).")
}

write.csv(data.frame(SNP = snp_list), "lead_snps_coloc.csv", row.names = FALSE)
write.csv(gwas_assoc, "gwas_summary_stats.csv", row.names = FALSE)

# --- 4. FETCH eQTL DATA (GTEx BRAIN via Qtlizer) ---
message("\n>>> Step 2: Fetching Brain eQTL Data (Qtlizer / GTEx Brain)...")

map_snps_batch <- function(snps_vec, batch_size = 50) {
  res_list <- list()
  batches  <- split(snps_vec, ceiling(seq_along(snps_vec) / batch_size))
  message("    Processing ", length(batches), " batch(es)...")
  for (i in seq_along(batches)) {
    try({
      tmp <- Qtlizer::get_qtls(batches[[i]], ld_method = ld_method, corr = ld_corr)
      res_list[[i]] <- tmp
    }, silent = TRUE)
    Sys.sleep(0.2)
  }
  dplyr::bind_rows(res_list)
}

qtls_raw <- map_snps_batch(snp_list)

# Tissue filter
if (!is.null(qtls_raw) && nrow(qtls_raw) > 0) {
  if (qtl_tissue != "all") {
    qtls <- qtls_raw[grepl(qtl_tissue, qtls_raw$tissue, ignore.case = TRUE), ]
    if (nrow(qtls) == 0) {
      message("    Warning: No Brain QTLs found after filter. Using all tissues.")
      qtls <- qtls_raw
    }
  } else {
    qtls <- qtls_raw
  }
} else {
  message("!!! Qtlizer returned no results. Using gene-level backup.")
  # Minimal backup eQTL table
  qtls <- data.frame(
    rsid   = snp_list[1:min(10, length(snp_list))],
    gene   = c("APOE","BIN1","CLU","PICALM","CR1","ABCA7","TREM2","PLCG2","ADAM10","MS4A6A")[1:min(10,length(snp_list))],
    tissue = "Brain - Cortex",
    beta   = rep(0.15, min(10, length(snp_list))),
    se     = rep(0.04, min(10, length(snp_list))),
    pvalue = rep(1e-6, min(10, length(snp_list))),
    stringsAsFactors = FALSE
  )
}

write.csv(qtls, "qtls_brain_coloc.csv", row.names = FALSE)
message("    Brain eQTL entries: ", nrow(qtls))
message("    Unique genes: ", length(unique(qtls$gene)))

# --- 5. PREPARE HyPrColoc INPUT ---
message("\n>>> Step 3: Preparing HyPrColoc Input Matrices...")

# Identify overlapping SNPs between GWAS and eQTL datasets
# HyPrColoc requires: a matrix of betas and a matrix of SEs
# Rows = SNPs, Columns = traits (GWAS trait + each eQTL gene)

# Get common SNPs
qtl_snp_col <- if ("rsid" %in% colnames(qtls)) "rsid" else colnames(qtls)[1]
common_snps <- intersect(gwas_assoc$SNP, qtls[[qtl_snp_col]])

if (length(common_snps) < 3) {
  message("!!! WARNING: Fewer than 3 overlapping SNPs between GWAS and eQTL.")
  message("    Colocalization requires overlapping SNPs with summary statistics.")
  message("    Consider using real GWAS summary stats (full locus files).")
  message("    Proceeding with available SNPs for demonstration...")
  common_snps <- gwas_assoc$SNP[1:min(5, nrow(gwas_assoc))]
}

message("    Overlapping SNPs for colocalization: ", length(common_snps))

# Subset GWAS stats to common SNPs
gwas_sub <- gwas_assoc %>%
  dplyr::filter(SNP %in% common_snps) %>%
  dplyr::distinct(SNP, .keep_all = TRUE)

# Get unique genes to test
genes_to_test <- unique(qtls$gene)
message("    Genes to test for colocalization: ", length(genes_to_test))

# Build beta and SE matrices (SNPs x Traits)
# Trait 1 = AD GWAS, Traits 2..N = eQTL per gene
build_matrices <- function(gwas_df, qtl_df, snps, genes, qtl_snp_col) {
  n_snps   <- length(snps)
  n_traits <- 1 + length(genes)  # GWAS + one per gene

  beta_mat <- matrix(NA, nrow = n_snps, ncol = n_traits,
                     dimnames = list(snps, c("AD_GWAS", genes)))
  se_mat   <- matrix(NA, nrow = n_snps, ncol = n_traits,
                     dimnames = list(snps, c("AD_GWAS", genes)))

  # Fill GWAS column
  gwas_matched <- gwas_df[match(snps, gwas_df$SNP), ]
  beta_mat[, "AD_GWAS"] <- gwas_matched$BETA
  se_mat[,   "AD_GWAS"] <- gwas_matched$SE

  # Fill eQTL columns
  for (g in genes) {
    qtl_gene <- qtl_df %>% dplyr::filter(gene == g)
    snp_match <- match(snps, qtl_gene[[qtl_snp_col]])
    valid     <- !is.na(snp_match)

    if (sum(valid) > 0 && "beta" %in% colnames(qtl_gene)) {
      beta_mat[valid, g] <- qtl_gene$beta[snp_match[valid]]
      se_mat[valid,   g] <- qtl_gene$se[snp_match[valid]]
    }
  }
  list(betas = beta_mat, ses = se_mat)
}

mats <- build_matrices(gwas_sub, qtls, common_snps, genes_to_test, qtl_snp_col)

# Remove traits (genes) with too many missing values (>80% NA)
na_frac   <- colMeans(is.na(mats$betas))
keep_cols <- na_frac < 0.8
mats$betas <- mats$betas[, keep_cols, drop = FALSE]
mats$ses   <- mats$ses[,   keep_cols, drop = FALSE]

# Remove SNP rows with all NA
keep_rows <- rowSums(!is.na(mats$betas)) > 1
mats$betas <- mats$betas[keep_rows, , drop = FALSE]
mats$ses   <- mats$ses[keep_rows,   , drop = FALSE]

# Replace remaining NAs with 0 (no effect) for HyPrColoc
mats$betas[is.na(mats$betas)] <- 0
mats$ses[is.na(mats$ses)]     <- 1  # Uninformative SE for missing

message("    Final matrix: ", nrow(mats$betas), " SNPs x ",
        ncol(mats$betas), " traits")

# --- 6. RUN HyPrColoc ---
message("\n>>> Step 4: Running HyPrColoc Multi-Trait Colocalization...")

hypr_results    <- NULL
hypr_results_df <- NULL

tryCatch({
  trait_names <- colnames(mats$betas)
  binary_traits <- rep(0, ncol(mats$betas))  # 0 = continuous traits

  hypr_results <- hyprcoloc(
    effect.est   = mats$betas,
    effect.se    = mats$ses,
    trait.names  = trait_names,
    snp.id       = rownames(mats$betas),
    binary.traits = binary_traits,
    prior.1      = 1e-4,   # Prior probability a SNP is associated with one trait
    prior.c      = 0.02    # Conditional prior for shared associations
  )

  if (!is.null(hypr_results)) {
    hypr_results_df <- as.data.frame(hypr_results$results)
    message("    HyPrColoc completed. Clusters found: ", nrow(hypr_results_df))
  }
}, error = function(e) {
  message("!!! HyPrColoc error: ", e$message)
})

# --- 7. FILTER & INTERPRET RESULTS ---
message("\n>>> Step 5: Filtering and Interpreting Results...")

if (!is.null(hypr_results_df) && nrow(hypr_results_df) > 0) {

  # Filter by posterior probability threshold
  coloc_hits <- hypr_results_df %>%
    dplyr::filter(posterior_prob >= pp4_threshold |
                  regional_prob  >= regional_prob)

  message("    Colocalized clusters (PP >= ", pp4_threshold, "): ", nrow(coloc_hits))

  # Classify evidence strength
  hypr_results_df <- hypr_results_df %>%
    dplyr::mutate(
      evidence_strength = dplyr::case_when(
        posterior_prob >= 0.9 ~ "Strong",
        posterior_prob >= 0.7 ~ "Moderate",
        posterior_prob >= 0.5 ~ "Suggestive",
        TRUE                  ~ "Weak"
      )
    )

  # Save full results
  write.csv(hypr_results_df, "HyPrColoc_Full_Results.csv", row.names = FALSE)
  write.csv(coloc_hits,      "HyPrColoc_Significant_Hits.csv", row.names = FALSE)
  message("    Results saved.")

} else {
  message("    No HyPrColoc results to filter.")
  message("    This is expected if using placeholder summary statistics.")
  message("    Provide real GWAS locus summary stats for meaningful results.")
  hypr_results_df <- data.frame(
    traits         = character(),
    posterior_prob = numeric(),
    regional_prob  = numeric(),
    candidate_snp  = character(),
    evidence_strength = character()
  )
  write.csv(hypr_results_df, "HyPrColoc_Full_Results.csv", row.names = FALSE)
}

# --- 8. VISUALISATION ---
message("\n>>> Step 6: Creating Visualisations...")

tryCatch({
  if (nrow(hypr_results_df) > 0 && "posterior_prob" %in% colnames(hypr_results_df)) {

    # Plot 1: Posterior probability per cluster
    p1 <- ggplot(hypr_results_df,
                 aes(x = reorder(traits, posterior_prob),
                     y = posterior_prob,
                     fill = evidence_strength)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = pp4_threshold, linetype = "dashed",
                 color = "red", linewidth = 0.8) +
      coord_flip() +
      scale_fill_manual(values = c("Strong"     = "#2ecc71",
                                   "Moderate"   = "#3498db",
                                   "Suggestive" = "#f39c12",
                                   "Weak"       = "#bdc3c7")) +
      labs(title    = "HyPrColoc: Multi-Trait Colocalization",
           subtitle = paste0("AD GWAS vs GTEx Brain eQTLs | PP threshold = ", pp4_threshold),
           x        = "Colocalized Trait Cluster",
           y        = "Posterior Probability of Shared Causal Variant",
           fill     = "Evidence") +
      theme_bw(base_size = 12) +
      theme(plot.title    = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave("HyPrColoc_PosteriorProb_Plot.png", plot = p1,
           width = 10, height = 6, dpi = 300, units = "in")
    message("    Plot saved: HyPrColoc_PosteriorProb_Plot.png")

    # Plot 2: Regional vs Posterior probability scatter
    p2 <- ggplot(hypr_results_df,
                 aes(x = regional_prob, y = posterior_prob,
                     color = evidence_strength, label = traits)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_hline(yintercept = pp4_threshold, linetype = "dashed", color = "red") +
      geom_vline(xintercept = regional_prob,  linetype = "dashed", color = "blue") +
      scale_color_manual(values = c("Strong"     = "#2ecc71",
                                    "Moderate"   = "#3498db",
                                    "Suggestive" = "#f39c12",
                                    "Weak"       = "#bdc3c7")) +
      labs(title    = "HyPrColoc: Regional vs Posterior Probability",
           subtitle = "Red dashed = PP threshold | Blue dashed = Regional threshold",
           x        = "Regional Colocalization Probability",
           y        = "Posterior Probability",
           color    = "Evidence") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave("HyPrColoc_Regional_vs_Posterior.png", plot = p2,
           width = 8, height = 6, dpi = 300, units = "in")
    message("    Plot saved: HyPrColoc_Regional_vs_Posterior.png")
  }
}, error = function(e) {
  message("    Visualisation skipped: ", e$message)
})

# --- 9. SUMMARY REPORT ---
message("\n>>> Step 7: Generating Summary Report...")

summary_report <- list(
  analysis_date      = Sys.time(),
  trait              = trait,
  n_gwas_snps        = length(snp_list),
  n_brain_eqtl_genes = length(unique(qtls$gene)),
  n_overlapping_snps = length(common_snps),
  n_traits_tested    = ncol(mats$betas),
  pp4_threshold      = pp4_threshold,
  n_coloc_hits       = if (!is.null(hypr_results_df)) {
                         sum(hypr_results_df$posterior_prob >= pp4_threshold, na.rm = TRUE)
                       } else { 0 },
  output_files       = c("HyPrColoc_Full_Results.csv",
                         "HyPrColoc_Significant_Hits.csv",
                         "HyPrColoc_PosteriorProb_Plot.png",
                         "HyPrColoc_Regional_vs_Posterior.png",
                         "gwas_summary_stats.csv",
                         "qtls_brain_coloc.csv")
)

capture.output(print(summary_report), file = "HyPrColoc_Summary.txt")
message("\n========================================")
message("  COLOCALIZATION PIPELINE COMPLETE")
message("  Results saved to: ", work_dir)
message("  Strong hits (PP >= ", pp4_threshold, "): ",
        summary_report$n_coloc_hits)
message("========================================\n")
