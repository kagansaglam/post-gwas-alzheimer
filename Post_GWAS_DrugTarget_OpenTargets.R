################################################################################
# Post-GWAS Analysis Pipeline: Drug Target Prioritization Module
# METHOD: OpenTargets Platform GraphQL API
# LOGIC: Candidate Genes -> OpenTargets -> Drug/Clinical Evidence -> Report
# NOTE: English comments used to prevent Encoding errors on Windows.
# REQUIRES: Candidate genes from Post_GWAS_Analysis_new.R or SummaryStats pipeline
################################################################################

# --- 1. SETUP & LIBRARIES ---
required_pkgs <- c("tidyverse", "httr", "jsonlite", "dplyr",
                   "ggplot2", "readr", "stringr", "gwasrapidd")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(gwasrapidd)
})

# --- 2. PARAMETERS ---
work_dir       <- "~/Documents/run_results"
ot_api_url     <- "https://api.platform.opentargets.org/api/v4/graphql"

# Drug phase filter: 3 = Phase III, 4 = Approved
# We include Phase II+ to capture clinical-stage compounds
min_drug_phase <- 2

# OpenTargets disease ID for Alzheimer's Disease
ad_disease_id  <- "MONDO_0004975"

# Association score threshold (0-1): how strongly linked gene is to AD
assoc_score_threshold <- 0.1

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. LOAD CANDIDATE GENES ---
message("\n>>> Step 1: Loading Candidate Genes...")

# Try to load genes from previous pipeline output
candidate_genes <- NULL

tryCatch({
  if (file.exists("gene_mapping.csv")) {
    gene_map        <- read.csv("gene_mapping.csv")
    candidate_genes <- unique(gene_map$SYMBOL)
    message("    Loaded ", length(candidate_genes),
            " genes from gene_mapping.csv")
  } else if (file.exists("gene_mapping_real.tsv")) {
    gene_map        <- read.delim("gene_mapping_real.tsv")
    candidate_genes <- unique(gene_map$SYMBOL)
    message("    Loaded ", length(candidate_genes),
            " genes from gene_mapping_real.tsv")
  }
}, error = function(e) {
  message("    Could not load gene mapping file: ", e$message)
})

# Fallback: use established AD GWAS genes
if (is.null(candidate_genes) || length(candidate_genes) == 0) {
  message("    Using curated AD GWAS gene list as fallback.")
  candidate_genes <- c(
    "APOE", "BIN1", "CLU", "ABCA7", "CR1", "PICALM", "MS4A6A", "CD33",
    "MS4A4E", "CD2AP", "EPHA1", "INPP5D", "MEF2C", "TREM2", "SORL1",
    "PLCG2", "ADAM10", "ACE", "TOMM40", "BDNF", "MAPT", "GRN",
    "TARDBP", "FUS", "PINK1", "LRRK2", "APP", "PSEN1", "PSEN2"
  )
}

write.csv(data.frame(gene = candidate_genes),
          "candidate_genes_input.csv", row.names = FALSE)
message("    Total candidate genes: ", length(candidate_genes))

# --- 4. OpenTargets GRAPHQL QUERIES ---
message("\n>>> Step 2: Querying OpenTargets Platform...")

# Helper: run a GraphQL query
run_graphql <- function(query, variables = list()) {
  resp <- httr::POST(
    url    = ot_api_url,
    body   = jsonlite::toJSON(list(query = query, variables = variables),
                              auto_unbox = TRUE),
    encode = "json",
    httr::add_headers("Content-Type" = "application/json"),
    httr::timeout(30)
  )
  if (httr::status_code(resp) != 200) return(NULL)
  jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                     flatten = TRUE)
}

# --- 4A. Get Ensembl IDs for gene symbols ---
message("    4a: Resolving gene symbols to Ensembl IDs...")

gene_id_query <- '
query GeneSearch($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"]) {
    hits {
      id
      name
      entity
    }
  }
}
'

resolve_gene_ids <- function(symbols) {
  results <- list()
  for (sym in symbols) {
    tryCatch({
      res <- run_graphql(gene_id_query, list(symbol = sym))
      if (!is.null(res$data$search$hits) && nrow(res$data$search$hits) > 0) {
        hit <- res$data$search$hits[1, ]
        results[[sym]] <- data.frame(
          symbol     = sym,
          ensembl_id = hit$id,
          name       = hit$name,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
    Sys.sleep(0.1)
  }
  dplyr::bind_rows(results)
}

gene_ids <- resolve_gene_ids(candidate_genes)
message("    Resolved: ", nrow(gene_ids), "/", length(candidate_genes), " genes")

# --- 4B. Get AD association scores ---
message("    4b: Fetching AD association scores...")

assoc_query <- '
query TargetAssociation($ensemblId: String!, $diseaseId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(enableIndirect: true) {
      rows {
        disease { id name }
        score
        datatypeScores { componentId score }
      }
    }
  }
}
'

fetch_association <- function(ensembl_id, symbol, disease_id) {
  tryCatch({
    res <- run_graphql(assoc_query,
                       list(ensemblId = ensembl_id, diseaseId = disease_id))
    rows <- res$data$target$associatedDiseases$rows
    if (is.null(rows) || nrow(rows) == 0) return(NULL)

    # Filter to AD specifically
    ad_row <- rows[grepl("alzheimer|dementia", rows$disease.name,
                         ignore.case = TRUE), ]
    if (nrow(ad_row) == 0) ad_row <- rows[rows$disease.id == disease_id, ]
    if (nrow(ad_row) == 0) return(NULL)

    data.frame(
      symbol        = symbol,
      ensembl_id    = ensembl_id,
      disease       = ad_row$disease.name[1],
      overall_score = ad_row$score[1],
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

assoc_results <- list()
for (i in seq_len(nrow(gene_ids))) {
  res <- fetch_association(gene_ids$ensembl_id[i],
                           gene_ids$symbol[i], ad_disease_id)
  if (!is.null(res)) assoc_results[[i]] <- res
  Sys.sleep(0.15)
}
assoc_df <- dplyr::bind_rows(assoc_results)
message("    Genes with AD association scores: ", nrow(assoc_df))

# --- 4C. Get drug information ---
message("    4c: Fetching drug/clinical evidence...")

drug_query <- '
query KnownDrugs($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    knownDrugs {
      rows {
        drug { id name maximumClinicalTrialPhase}
        mechanismOfAction
        disease { name }
        phase
        status
        urls { url }
      }
    }
  }
}
'

fetch_drugs <- function(ensembl_id, symbol) {
  tryCatch({
    res   <- run_graphql(drug_query, list(ensemblId = ensembl_id))
    rows  <- res$data$target$knownDrugs$rows
    if (is.null(rows) || nrow(rows) == 0) return(NULL)

    df <- data.frame(
      symbol             = symbol,
      ensembl_id         = ensembl_id,
      drug_id            = rows$drug.id,
      drug_name          = rows$drug.name,
      mechanism          = rows$mechanismOfAction,
      disease_indication = rows$disease.name,
      phase              = rows$phase,
      status             = ifelse(is.null(rows$status), NA, rows$status),
      max_phase          = rows$drug.maximumClinicalTrialPhase,
      stringsAsFactors   = FALSE
    )
    df[df$phase >= min_drug_phase, ]
  }, error = function(e) NULL)
}

drug_results <- list()
for (i in seq_len(nrow(gene_ids))) {
  res <- fetch_drugs(gene_ids$ensembl_id[i], gene_ids$symbol[i])
  if (!is.null(res) && nrow(res) > 0) drug_results[[i]] <- res
  Sys.sleep(0.15)
}
drug_df <- dplyr::bind_rows(drug_results)
message("    Drug-gene interactions found: ", nrow(drug_df))

# --- 5. PRIORITIZATION SCORING ---
message("\n>>> Step 3: Prioritizing Drug Targets...")

# Phase labels
phase_labels <- c("1" = "Phase I", "2" = "Phase II",
                  "3" = "Phase III", "4" = "Approved")

# Summarize drugs per gene
drug_summary <- drug_df %>%
  group_by(symbol) %>%
  summarise(
    n_drugs        = n_distinct(drug_name),
    max_phase      = max(phase, na.rm = TRUE),
    approved_drugs = sum(phase >= 4, na.rm = TRUE),
    phase3_drugs   = sum(phase == 3, na.rm = TRUE),
    phase2_drugs   = sum(phase == 2, na.rm = TRUE),
    drug_names     = paste(unique(drug_name), collapse = "; "),
    mechanisms     = paste(unique(mechanism), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    phase_label = dplyr::recode(as.character(max_phase), !!!phase_labels,
                                .default = "Unknown"),
    druggability = dplyr::case_when(
      approved_drugs > 0 ~ "Approved Drug Exists",
      phase3_drugs   > 0 ~ "Phase III",
      phase2_drugs   > 0 ~ "Phase II",
      TRUE               ~ "Early Stage"
    )
  )

# Join with association scores
priority_df <- assoc_df %>%
  left_join(drug_summary, by = "symbol") %>%
  mutate(
    n_drugs    = tidyr::replace_na(n_drugs, 0),
    max_phase  = tidyr::replace_na(max_phase, 0),
    # Priority score: weighted combination of AD association + drug phase
    priority_score = overall_score * 0.6 + (max_phase / 4) * 0.4
  ) %>%
  arrange(desc(priority_score)) %>%
  filter(overall_score >= assoc_score_threshold)

message("    Prioritized targets: ", nrow(priority_df))

# Top targets
top_targets <- priority_df %>%
  filter(n_drugs > 0) %>%
  head(20)

message("    Top druggable targets (top 20):")
if (nrow(top_targets) > 0) {
  for (i in seq_len(min(10, nrow(top_targets)))) {
    message("      ", top_targets$symbol[i],
            " | AD score: ", round(top_targets$overall_score[i], 3),
            " | Max phase: ", top_targets$max_phase[i],
            " | Drugs: ", top_targets$drug_names[i])
  }
}

# Save results
write.csv(priority_df, "DrugTarget_Prioritized.csv",    row.names = FALSE)
write.csv(drug_df,     "DrugTarget_AllDrugs.csv",       row.names = FALSE)
write.csv(assoc_df,    "DrugTarget_ADScores.csv",       row.names = FALSE)
write.csv(top_targets, "DrugTarget_Top20.csv",          row.names = FALSE)
message("    Results saved.")

# --- 6. VISUALISATIONS ---
message("\n>>> Step 4: Creating Visualisations...")

tryCatch({

  # Plot 1: Top targets by priority score, colored by druggability
  if (nrow(priority_df) > 0) {
    plot_df <- priority_df %>%
      head(25) %>%
      mutate(
        druggability = tidyr::replace_na(druggability, "No Drug Data"),
        symbol       = factor(symbol, levels = rev(symbol))
      )

    p1 <- ggplot(plot_df,
                 aes(x = symbol, y = priority_score, fill = druggability)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = ifelse(n_drugs > 0,
                                   paste0(n_drugs, " drug(s)"), "")),
                hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_manual(values = c(
        "Approved Drug Exists" = "#2ecc71",
        "Phase III"            = "#3498db",
        "Phase II"             = "#f39c12",
        "Early Stage"          = "#e74c3c",
        "No Drug Data"         = "#bdc3c7"
      )) +
      scale_y_continuous(limits = c(0, 1.15)) +
      labs(
        title    = "Drug Target Prioritization: Alzheimer's Disease",
        subtitle = "OpenTargets Platform | Priority = AD association score × drug phase",
        x        = "Candidate Gene",
        y        = "Priority Score",
        fill     = "Druggability Status"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title    = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave("DrugTarget_PriorityPlot.png", plot = p1,
           width = 11, height = 8, dpi = 300, units = "in")
    message("    Priority plot saved.")
  }

  # Plot 2: Drug phase distribution across all genes
  if (nrow(drug_df) > 0) {
    phase_count <- drug_df %>%
      mutate(phase_label = dplyr::recode(
        as.character(phase),
        "1" = "Phase I", "2" = "Phase II",
        "3" = "Phase III", "4" = "Approved",
        .default = "Unknown"
      )) %>%
      count(phase_label) %>%
      mutate(phase_label = factor(phase_label,
                                  levels = c("Phase I","Phase II",
                                             "Phase III","Approved","Unknown")))

    p2 <- ggplot(phase_count, aes(x = phase_label, y = n, fill = phase_label)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
      scale_fill_manual(values = c(
        "Phase I"   = "#e74c3c", "Phase II"  = "#f39c12",
        "Phase III" = "#3498db", "Approved"  = "#2ecc71",
        "Unknown"   = "#bdc3c7"
      ), guide = "none") +
      labs(
        title    = "Drug Pipeline Distribution",
        subtitle = "Across all AD GWAS candidate genes (OpenTargets)",
        x        = "Clinical Phase",
        y        = "Number of Drug-Gene Interactions"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave("DrugTarget_PhasePlot.png", plot = p2,
           width = 8, height = 5, dpi = 300, units = "in")
    message("    Phase distribution plot saved.")
  }

  # Plot 3: AD association score vs max drug phase bubble chart
  if (nrow(priority_df) > 0 && sum(priority_df$n_drugs > 0) > 0) {
    bubble_df <- priority_df %>%
      filter(n_drugs > 0) %>%
      mutate(max_phase = as.numeric(max_phase))

    p3 <- ggplot(bubble_df,
                 aes(x = overall_score, y = max_phase,
                     size = n_drugs, color = druggability, label = symbol)) +
      geom_point(alpha = 0.7) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 15) +
      scale_size_continuous(range = c(3, 12), name = "# Drugs") +
      scale_color_manual(values = c(
        "Approved Drug Exists" = "#2ecc71",
        "Phase III"            = "#3498db",
        "Phase II"             = "#f39c12",
        "Early Stage"          = "#e74c3c"
      )) +
      scale_y_continuous(breaks = 1:4,
                         labels = c("Phase I","Phase II","Phase III","Approved")) +
      labs(
        title    = "AD Association Score vs Drug Development Stage",
        subtitle = "Bubble size = number of drugs | Top-right = highest priority",
        x        = "OpenTargets AD Association Score",
        y        = "Maximum Drug Phase",
        color    = "Druggability"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    # ggrepel optional — skip if not installed
    if (!require("ggrepel", quietly = TRUE)) {
      p3 <- p3 + geom_text(size = 3, vjust = -1)
    }

    ggsave("DrugTarget_BubblePlot.png", plot = p3,
           width = 10, height = 7, dpi = 300, units = "in")
    message("    Bubble chart saved.")
  }

}, error = function(e) {
  message("    Visualisation error: ", e$message)
})

# --- 7. SUMMARY REPORT ---
message("\n>>> Step 5: Generating Summary Report...")

approved_genes <- priority_df %>%
  filter(druggability == "Approved Drug Exists") %>%
  pull(symbol)

phase3_genes <- priority_df %>%
  filter(druggability == "Phase III") %>%
  pull(symbol)

summary_lines <- c(
  "========================================",
  "  DRUG TARGET PRIORITIZATION REPORT",
  paste0("  Disease: Alzheimer's Disease (", ad_disease_id, ")"),
  paste0("  Data source: OpenTargets Platform"),
  paste0("  Date: ", Sys.Date()),
  "----------------------------------------",
  paste0("  Candidate genes queried:     ", length(candidate_genes)),
  paste0("  Genes with AD scores:        ", nrow(assoc_df)),
  paste0("  Genes with drug interactions:", nrow(drug_summary)),
  paste0("  Total drug interactions:     ", nrow(drug_df)),
  "----------------------------------------",
  paste0("  Genes with APPROVED drugs:   ", length(approved_genes)),
  if(length(approved_genes) > 0)
    paste0("    -> ", paste(approved_genes, collapse = ", ")),
  paste0("  Genes in Phase III trials:   ", length(phase3_genes)),
  if(length(phase3_genes) > 0)
    paste0("    -> ", paste(phase3_genes, collapse = ", ")),
  "----------------------------------------",
  "  Output files:",
  "    DrugTarget_Prioritized.csv",
  "    DrugTarget_Top20.csv",
  "    DrugTarget_AllDrugs.csv",
  "    DrugTarget_PriorityPlot.png",
  "    DrugTarget_PhasePlot.png",
  "    DrugTarget_BubblePlot.png",
  "========================================"
)

writeLines(summary_lines, "DrugTarget_Summary.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")
