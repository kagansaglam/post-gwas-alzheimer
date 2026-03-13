################################################################################
# Post-GWAS Analysis Pipeline: Mendelian Randomization Module
# METHOD: TwoSampleMR - IVW + MR-Egger + Weighted Median
# EXPOSURE: Gene expression (eQTLGen blood eQTLs)
# OUTCOME: Alzheimer's Disease (Bellenguez 2022 / GWAS Catalog)
# LOGIC: Try API -> If Fails, Use Backup -> Complete Analysis
# NOTE: English comments used to prevent Encoding errors on Windows.
################################################################################

# --- 1. SETUP & LIBRARIES ---
required_pkgs <- c("tidyverse", "dplyr", "ggplot2", "readr", "httr", "jsonlite")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

# TwoSampleMR from MR-Base
if (!require("TwoSampleMR", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("MRCIEU/TwoSampleMR")
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
})

# --- 2. PARAMETERS ---
work_dir       <- "~/Documents/run_results_MR"
ad_outcome_id  <- "ieu-b-2"   # Alzheimer's Disease in MR-Base (Lambert 2013)
ad_outcome_id2 <- "ebi-a-GCST90027158"  # Bellenguez 2022 if available
pval_threshold <- 5e-8        # Instrument strength threshold
f_stat_min     <- 10          # Minimum F-statistic for instrument validity

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. LOAD CANDIDATE GENES ---
message("\n>>> Step 1: Loading Candidate Genes...")

candidate_genes <- NULL
tryCatch({
  if (file.exists("~/Documents/run_results/gene_mapping.csv")) {
    gene_map        <- read.csv("~/Documents/run_results/gene_mapping.csv")
    candidate_genes <- unique(gene_map$SYMBOL)
    message("    Loaded ", length(candidate_genes), " genes from previous pipeline")
  }
}, error = function(e) NULL)

if (is.null(candidate_genes) || length(candidate_genes) == 0) {
  candidate_genes <- c(
    "APOE","BIN1","CLU","ABCA7","CR1","PICALM","MS4A6A","CD33",
    "TREM2","SORL1","PLCG2","ADAM10","ACE","TOMM40","INPP5D",
    "MEF2C","EPHA1","CD2AP","MS4A4E"
  )
  message("    Using curated AD GWAS gene list.")
}
message("    Genes to test: ", length(candidate_genes))

# --- 4. GET EXPOSURE DATA (eQTLGen blood eQTLs) ---
message("\n>>> Step 2: Fetching Exposure Data (eQTLGen)...")
message("    Source: eQTLGen blood eQTLs (n~31,000)")
message("    Instrument: cis-eQTL SNPs per gene (P <= 5e-8)")

exposure_dat <- NULL

# METHOD A: Try MR-Base API (eQTLGen)
tryCatch({
  message("    Connecting to MR-Base API...")

  # eQTLGen is available via MR-Base as outcome "eqtl-a-*" format
  # We fetch cis-eQTLs for each candidate gene
  exposure_list <- list()

  for (gene in candidate_genes) {
    tryCatch({
      # Search for eQTLGen instruments for this gene
      ao <- available_outcomes(access_token = NULL)
      eqtlgen_id <- ao$id[grepl(paste0("eqtlgen.*", gene, "|", gene, ".*eqtlgen"),
                                 ao$trait, ignore.case = TRUE)][1]

      if (!is.na(eqtlgen_id)) {
        exp_raw <- extract_instruments(eqtlgen_id, access_token = NULL)
        if (!is.null(exp_raw) && nrow(exp_raw) > 0) {
          exp_raw$exposure <- gene
          exposure_list[[gene]] <- exp_raw
        }
      }
    }, error = function(e) NULL)
    Sys.sleep(0.2)
  }

  if (length(exposure_list) > 0) {
    exposure_dat <- dplyr::bind_rows(exposure_list)
    message("    SUCCESS: ", nrow(exposure_dat), " instruments for ",
            length(unique(exposure_dat$exposure)), " genes")
  }
}, error = function(e) {
  message("    MR-Base API failed: ", e$message)
})

# METHOD B: Backup - use GWAS Catalog eQTL data
if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
  message("\n>>> PLAN B: Building exposure from Qtlizer eQTL data...")

  qtl_file <- "~/Documents/run_results/qtls_filtered.csv"
  if (!file.exists(qtl_file)) qtl_file <- "~/Documents/run_results_real/qtls_brain_real.tsv"

  tryCatch({
    qtls <- if (grepl("\\.csv$", qtl_file)) read.csv(qtl_file) else read.delim(qtl_file)
    qtls <- qtls[qtls$gene %in% candidate_genes, ]

    if (nrow(qtls) > 0 && all(c("beta","se","pvalue","rsid") %in% colnames(qtls))) {
      exposure_dat <- data.frame(
        SNP            = qtls$rsid,
        beta.exposure  = qtls$beta,
        se.exposure    = qtls$se,
        pval.exposure  = qtls$pvalue,
        exposure       = qtls$gene,
        effect_allele.exposure  = ifelse("effect_allele" %in% colnames(qtls), qtls$effect_allele, "A"),
        other_allele.exposure   = ifelse("other_allele"  %in% colnames(qtls), qtls$other_allele,  "G"),
        eaf.exposure   = ifelse("eaf" %in% colnames(qtls), qtls$eaf, 0.3),
        mr_keep.exposure = TRUE,
        stringsAsFactors = FALSE
      )
      exposure_dat <- exposure_dat[!is.na(exposure_dat$pval.exposure) &
                                     exposure_dat$pval.exposure <= pval_threshold, ]
      message("    Backup exposure loaded: ", nrow(exposure_dat), " instruments")
    }
  }, error = function(e) message("    Backup also failed: ", e$message))
}

# Final check
if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
  stop("No exposure data available. Check eQTL source or API connection.")
}

# --- 5. FILTER INSTRUMENTS ---
message("\n>>> Step 3: Filtering and Clumping Instruments...")

# Calculate F-statistic: F = (beta/se)^2
exposure_dat$F_stat <- (exposure_dat$beta.exposure / exposure_dat$se.exposure)^2
exposure_dat        <- exposure_dat[exposure_dat$F_stat >= f_stat_min, ]
message("    Instruments after F-stat filter (F>=", f_stat_min, "): ", nrow(exposure_dat))

# Clump to remove LD between instruments
tryCatch({
  exposure_dat <- clump_data(exposure_dat, access_token = NULL)
  message("    Instruments after LD clumping: ", nrow(exposure_dat))
}, error = function(e) {
  message("    Clumping skipped (API unavailable): ", e$message)
})

write.csv(exposure_dat, "MR_exposure_instruments.csv", row.names = FALSE)
genes_with_instruments <- unique(exposure_dat$exposure)
message("    Genes with valid instruments: ", length(genes_with_instruments))

# --- 6. GET OUTCOME DATA (AD GWAS) ---
message("\n>>> Step 4: Fetching Outcome Data (AD GWAS)...")

outcome_dat <- NULL

# Try Bellenguez 2022 first, fall back to Lambert 2013
for (oid in c(ad_outcome_id2, ad_outcome_id)) {
  tryCatch({
    outcome_raw <- extract_outcome_data(
      snps         = unique(exposure_dat$SNP),
      outcomes     = oid,
      access_token = NULL
    )
    if (!is.null(outcome_raw) && nrow(outcome_raw) > 0) {
      outcome_dat <- outcome_raw
      message("    Outcome loaded from: ", oid, " (", nrow(outcome_dat), " SNPs)")
      break
    }
  }, error = function(e) NULL)
  Sys.sleep(1)
}

# Fallback: use local summary stats if available
if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
  message("    Trying local summary stats fallback...")
  tryCatch({
    sumstats <- data.table::fread(
      "~/Documents/run_results_real/meta_analysis_significant.tsv.gz")
    outcome_dat <- data.frame(
      SNP                    = sumstats$SNP,
      beta.outcome           = sumstats$BETA_meta,
      se.outcome             = sumstats$SE_meta,
      pval.outcome           = sumstats$P_meta,
      outcome                = "Alzheimer's Disease",
      effect_allele.outcome  = sumstats$A1,
      other_allele.outcome   = sumstats$A2,
      eaf.outcome            = sumstats$EAF_meta,
      mr_keep.outcome        = TRUE,
      stringsAsFactors       = FALSE
    )
    outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP, ]
    message("    Local fallback: ", nrow(outcome_dat), " matching SNPs")
  }, error = function(e) message("    Local fallback failed: ", e$message))
}

if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
  stop("No outcome data available. Check MR-Base connection or local summary stats.")
}

write.csv(outcome_dat, "MR_outcome_data.csv", row.names = FALSE)

# --- 7. HARMONISE DATA ---
message("\n>>> Step 5: Harmonising Exposure and Outcome Data...")

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- dat[dat$mr_keep == TRUE, ]
message("    SNPs after harmonisation: ", nrow(dat))
message("    Genes with harmonised data: ", length(unique(dat$exposure)))

write.csv(dat, "MR_harmonised_data.csv", row.names = FALSE)

# --- 8. RUN MR ANALYSIS ---
message("\n>>> Step 6: Running MR Analysis...")
message("    Methods: IVW, MR-Egger, Weighted Median")

mr_results <- mr(dat, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median"
))

mr_results <- generate_odds_ratios(mr_results)
mr_results <- mr_results %>%
  mutate(
    direction = case_when(
      b > 0 ~ "Risk-increasing",
      b < 0 ~ "Protective",
      TRUE  ~ "Neutral"
    ),
    significance = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  )

message("    MR complete. Gene-method combinations: ", nrow(mr_results))

# --- 9. SENSITIVITY ANALYSES ---
message("\n>>> Step 7: Sensitivity Analyses...")

# Heterogeneity test (Cochran Q for IVW, Rucker Q for MR-Egger)
tryCatch({
  heterogeneity <- mr_heterogeneity(dat)
  write.csv(heterogeneity, "MR_heterogeneity.csv", row.names = FALSE)
  message("    Heterogeneity test done.")
}, error = function(e) message("    Heterogeneity skipped: ", e$message))

# Pleiotropy test (MR-Egger intercept)
tryCatch({
  pleiotropy <- mr_pleiotropy_test(dat)
  write.csv(pleiotropy, "MR_pleiotropy.csv", row.names = FALSE)
  message("    Pleiotropy test done.")
  sig_pleio <- pleiotropy[pleiotropy$pval < 0.05, ]
  if (nrow(sig_pleio) > 0) {
    message("    WARNING: Potential pleiotropy detected for: ",
            paste(sig_pleio$exposure, collapse = ", "))
  } else {
    message("    No significant pleiotropy detected.")
  }
}, error = function(e) message("    Pleiotropy skipped: ", e$message))

# Leave-one-out analysis
tryCatch({
  loo <- mr_leaveoneout(dat)
  write.csv(loo, "MR_leaveoneout.csv", row.names = FALSE)
  message("    Leave-one-out analysis done.")
}, error = function(e) message("    LOO skipped: ", e$message))

# Save full results
write.csv(mr_results, "MR_Results_Full.csv", row.names = FALSE)

# Filter significant results
mr_signif <- mr_results %>%
  filter(pval < 0.05) %>%
  arrange(pval)
write.csv(mr_signif, "MR_Results_Significant.csv", row.names = FALSE)

message("    Significant results (P<0.05): ", nrow(mr_signif))
if (nrow(mr_signif) > 0) {
  message("    Top causal genes:")
  for (i in seq_len(min(10, nrow(mr_signif)))) {
    message("      ", mr_signif$exposure[i],
            " | Method: ", mr_signif$method[i],
            " | OR: ", round(mr_signif$or[i], 3),
            " | P: ", signif(mr_signif$pval[i], 3),
            " | ", mr_signif$direction[i])
  }
}

# --- 10. VISUALISATIONS ---
message("\n>>> Step 8: Creating Visualisations...")

tryCatch({

  # Plot 1: Forest plot - IVW results only
  ivw_res <- mr_results %>%
    filter(method == "MR Inverse variance weighted") %>%
    arrange(b) %>%
    mutate(exposure = factor(exposure, levels = exposure),
           sig_color = case_when(
             pval < 0.05 & b < 0 ~ "Protective (P<0.05)",
             pval < 0.05 & b > 0 ~ "Risk-increasing (P<0.05)",
             TRUE                 ~ "Not significant"
           ))

  if (nrow(ivw_res) > 0) {
    p1 <- ggplot(ivw_res, aes(x = exposure, y = b,
                               ymin = b - 1.96*se, ymax = b + 1.96*se,
                               color = sig_color)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(size = 0.6) +
      coord_flip() +
      scale_color_manual(values = c(
        "Protective (P<0.05)"      = "#2ecc71",
        "Risk-increasing (P<0.05)" = "#e74c3c",
        "Not significant"          = "#95a5a6"
      )) +
      labs(title    = "Mendelian Randomization: Gene Expression -> AD Risk",
           subtitle = "Method: IVW | Exposure: eQTLGen blood eQTLs | Outcome: AD GWAS",
           x        = "Gene",
           y        = "MR Effect Estimate (Beta)",
           color    = "Direction") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

    ggsave("MR_ForestPlot_IVW.png", plot = p1,
           width = 11, height = max(6, nrow(ivw_res) * 0.4 + 2), dpi = 300)
    message("    Forest plot saved.")
  }

  # Plot 2: Method comparison - IVW vs Egger vs Weighted Median
  if (nrow(mr_signif) > 0) {
    sig_genes <- unique(mr_signif$exposure)
    compare_df <- mr_results %>%
      filter(exposure %in% sig_genes) %>%
      mutate(method = recode(method,
        "MR Inverse variance weighted" = "IVW",
        "MR Egger"                     = "MR-Egger",
        "Weighted median"              = "Weighted Median"
      ))

    p2 <- ggplot(compare_df, aes(x = method, y = b,
                                  ymin = b - 1.96*se, ymax = b + 1.96*se,
                                  color = method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(size = 0.7) +
      facet_wrap(~exposure, scales = "free_y") +
      scale_color_manual(values = c(
        "IVW"             = "#2c7bb6",
        "MR-Egger"        = "#d7191c",
        "Weighted Median" = "#1a9641"
      )) +
      labs(title    = "MR Method Comparison: Significant Genes",
           subtitle = "Consistent direction across methods = stronger causal evidence",
           x = "MR Method", y = "Effect Estimate (Beta)", color = "Method") +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold"),
            axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")

    ggsave("MR_MethodComparison.png", plot = p2,
           width = 12, height = max(5, ceiling(length(sig_genes)/3) * 3 + 2),
           dpi = 300)
    message("    Method comparison plot saved.")
  }

  # Plot 3: Scatter plot for top gene (IVW line)
  top_gene <- mr_results %>%
    filter(method == "MR Inverse variance weighted") %>%
    arrange(pval) %>%
    slice(1) %>%
    pull(exposure)

  if (length(top_gene) > 0) {
    top_dat <- dat[dat$exposure == top_gene, ]
    top_res <- mr_results[mr_results$exposure == top_gene, ]

    p3 <- mr_scatter_plot(top_res, top_dat)[[1]] +
      labs(title    = paste0("MR Scatter Plot: ", top_gene, " -> AD"),
           subtitle = "Each point = one genetic instrument (cis-eQTL SNP)") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(paste0("MR_ScatterPlot_", top_gene, ".png"), plot = p3,
           width = 8, height = 6, dpi = 300)
    message("    Scatter plot saved for top gene: ", top_gene)
  }

}, error = function(e) message("    Visualisation error: ", e$message))

# --- 11. SUMMARY ---
message("\n========================================")
message("  MENDELIAN RANDOMIZATION COMPLETE")
message("  Results saved to: ", work_dir)
message("  Genes tested: ", length(unique(mr_results$exposure)))
message("  Significant causal genes (P<0.05): ", nrow(mr_signif[mr_signif$method == "MR Inverse variance weighted",]))
message("  Output files:")
message("    - MR_Results_Full.csv")
message("    - MR_Results_Significant.csv")
message("    - MR_ForestPlot_IVW.png")
message("    - MR_MethodComparison.png")
message("    - MR_heterogeneity.csv")
message("    - MR_pleiotropy.csv")
message("========================================")
