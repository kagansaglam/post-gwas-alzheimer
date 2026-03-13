################################################################################
# Post-GWAS Analysis Pipeline: Drug Target Prioritization Module
# METHOD: OpenTargets Platform API
# NOTE: English comments used to prevent Encoding errors on Windows.
################################################################################

# --- 1. SETUP & LIBRARIES ---
required_pkgs <- c("tidyverse", "httr", "jsonlite", "dplyr", "ggplot2", "readr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# --- 2. PARAMETERS ---
work_dir              <- "~/Documents/run_results"
ot_api_url            <- "https://api.platform.opentargets.org/api/v4/graphql"
min_drug_phase        <- 2
ad_disease_id         <- "MONDO_0004975"
assoc_score_threshold <- 0.1

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
message("Results will be saved to: ", work_dir)

# --- 3. LOAD CANDIDATE GENES ---
message("\n>>> Step 1: Loading Candidate Genes...")

candidate_genes <- NULL
tryCatch({
  if (file.exists("gene_mapping.csv")) {
    gene_map        <- read.csv("gene_mapping.csv")
    candidate_genes <- unique(gene_map$SYMBOL)
    message("    Loaded ", length(candidate_genes), " genes from gene_mapping.csv")
  } else if (file.exists("gene_mapping_real.tsv")) {
    gene_map        <- read.delim("gene_mapping_real.tsv")
    candidate_genes <- unique(gene_map$SYMBOL)
    message("    Loaded ", length(candidate_genes), " genes from gene_mapping_real.tsv")
  }
}, error = function(e) message("    Could not load gene mapping: ", e$message))

if (is.null(candidate_genes) || length(candidate_genes) == 0) {
  candidate_genes <- c(
    "APOE","BIN1","CLU","ABCA7","CR1","PICALM","MS4A6A","CD33",
    "MS4A4E","CD2AP","EPHA1","INPP5D","MEF2C","TREM2","SORL1",
    "PLCG2","ADAM10","ACE","TOMM40","PSEN1","PSEN2","APP"
  )
  message("    Using curated fallback gene list.")
}

write.csv(data.frame(gene = candidate_genes), "candidate_genes_input.csv", row.names = FALSE)
message("    Total candidate genes: ", length(candidate_genes))

# --- 4. QUERY OpenTargets ---
message("\n>>> Step 2: Querying OpenTargets Platform...")

# 4A. Resolve gene symbols to Ensembl IDs
message("    4a: Resolving gene symbols to Ensembl IDs...")

gene_search_url <- "https://api.platform.opentargets.org/api/v4/graphql"

resolve_one <- function(symbol) {
  q <- paste0('{ search(queryString: "', symbol, '", entityNames: ["target"]) { hits { id name } } }')
  tryCatch({
    resp <- POST(gene_search_url,
                 body = toJSON(list(query = q, variables = list()), auto_unbox = TRUE),
                 encode = "json", add_headers("Content-Type" = "application/json"), timeout(20))
    res  <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
    hits <- res$data$search$hits
    if (!is.null(hits) && nrow(hits) > 0) {
      return(data.frame(symbol = symbol, ensembl_id = hits$id[1], stringsAsFactors = FALSE))
    }
    return(NULL)
  }, error = function(e) NULL)
}

gene_ids <- dplyr::bind_rows(lapply(candidate_genes, function(s) {
  Sys.sleep(0.1)
  resolve_one(s)
}))
message("    Resolved: ", nrow(gene_ids), "/", length(candidate_genes), " genes")

# 4B. Get AD association scores (disease-centric direct POST)
message("    4b: Fetching AD association scores...")

assoc_df <- tryCatch({
  q <- '{ disease(efoId: "MONDO_0004975") { associatedTargets(page: {index: 0, size: 500}) { rows { target { approvedSymbol ensemblId } score } } } }'
  resp <- POST(ot_api_url,
               body = toJSON(list(query = q, variables = list()), auto_unbox = TRUE),
               encode = "json", add_headers("Content-Type" = "application/json"), timeout(60))
  result <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
  rows   <- result$data$disease$associatedTargets$rows
  if (is.null(rows) || nrow(rows) == 0) stop("Empty response")
  df <- data.frame(
    symbol        = rows$target.approvedSymbol,
    ensembl_id    = rows$target.ensemblId,
    overall_score = rows$score,
    stringsAsFactors = FALSE
  )
  df <- df[df$symbol %in% candidate_genes, ]
  message("    Genes with AD scores: ", nrow(df))
  df
}, error = function(e) {
  message("    AD score fetch failed: ", e$message)
  data.frame(symbol = character(), ensembl_id = character(),
             overall_score = numeric(), stringsAsFactors = FALSE)
})

# 4C. Get drug interactions
message("    4c: Fetching drug/clinical evidence...")

fetch_drugs_one <- function(ensembl_id, symbol) {
  q <- paste0('{ target(ensemblId: "', ensembl_id, '") { knownDrugs { rows {
    drug { id name maximumClinicalTrialPhase }
    mechanismOfAction disease { name } phase status } } } }')
  tryCatch({
    resp <- POST(ot_api_url,
                 body = toJSON(list(query = q, variables = list()), auto_unbox = TRUE),
                 encode = "json", add_headers("Content-Type" = "application/json"), timeout(20))
    res  <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
    rows <- res$data$target$knownDrugs$rows
    if (is.null(rows) || nrow(rows) == 0) return(NULL)
    df <- data.frame(
      symbol             = symbol,
      ensembl_id         = ensembl_id,
      drug_name          = rows$drug.name,
      mechanism          = rows$mechanismOfAction,
      disease_indication = rows$disease.name,
      phase              = rows$phase,
      max_phase          = rows$drug.maximumClinicalTrialPhase,
      stringsAsFactors   = FALSE
    )
    df[!is.na(df$phase) & df$phase >= min_drug_phase, ]
  }, error = function(e) NULL)
}

drug_results <- list()
for (i in seq_len(nrow(gene_ids))) {
  res <- fetch_drugs_one(gene_ids$ensembl_id[i], gene_ids$symbol[i])
  if (!is.null(res) && nrow(res) > 0) drug_results[[i]] <- res
  Sys.sleep(0.15)
}
drug_df <- dplyr::bind_rows(drug_results)
message("    Drug-gene interactions found: ", nrow(drug_df))

# --- 5. PRIORITIZATION ---
message("\n>>> Step 3: Prioritizing Drug Targets...")

drug_summary <- drug_df %>%
  group_by(symbol) %>%
  summarise(
    n_drugs        = n_distinct(drug_name),
    max_phase      = max(phase, na.rm = TRUE),
    approved_drugs = sum(phase >= 4, na.rm = TRUE),
    phase3_drugs   = sum(phase == 3, na.rm = TRUE),
    phase2_drugs   = sum(phase == 2, na.rm = TRUE),
    drug_names     = paste(unique(drug_name), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(druggability = case_when(
    approved_drugs > 0 ~ "Approved Drug Exists",
    phase3_drugs   > 0 ~ "Phase III",
    phase2_drugs   > 0 ~ "Phase II",
    TRUE               ~ "Early Stage"
  ))

if (nrow(assoc_df) > 0) {
  priority_df <- assoc_df %>%
    left_join(drug_summary, by = "symbol") %>%
    mutate(
      n_drugs        = replace_na(n_drugs, 0),
      max_phase      = replace_na(max_phase, 0),
      priority_score = overall_score * 0.6 + (max_phase / 4) * 0.4,
      druggability   = replace_na(druggability, "No Drug Data")
    ) %>%
    arrange(desc(priority_score)) %>%
    filter(overall_score >= assoc_score_threshold)
} else {
  message("    AD scores unavailable — using drug phase as priority.")
  priority_df <- drug_summary %>%
    mutate(overall_score = NA_real_, priority_score = max_phase / 4) %>%
    arrange(desc(priority_score))
}

message("    Prioritized targets: ", nrow(priority_df))
top_targets <- priority_df %>% filter(n_drugs > 0) %>% head(20)

write.csv(priority_df, "DrugTarget_Prioritized.csv", row.names = FALSE)
write.csv(drug_df,     "DrugTarget_AllDrugs.csv",    row.names = FALSE)
write.csv(assoc_df,    "DrugTarget_ADScores.csv",    row.names = FALSE)
write.csv(top_targets, "DrugTarget_Top20.csv",       row.names = FALSE)
message("    Results saved.")

# --- 6. VISUALISATIONS ---
message("\n>>> Step 4: Creating Visualisations...")

tryCatch({
  # Plot 1: Priority bar chart
  if (nrow(priority_df) > 0) {
    plot_df <- priority_df %>%
      head(25) %>%
      mutate(symbol = factor(symbol, levels = rev(symbol)),
             druggability = replace_na(druggability, "No Drug Data"))

    p1 <- ggplot(plot_df, aes(x = symbol, y = priority_score, fill = druggability)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c(
        "Approved Drug Exists" = "#2ecc71", "Phase III" = "#3498db",
        "Phase II" = "#f39c12", "Early Stage" = "#e74c3c",
        "No Drug Data" = "#bdc3c7")) +
      labs(title    = "Drug Target Prioritization: Alzheimer's Disease",
           subtitle = "OpenTargets Platform | Priority = AD score x drug phase",
           x = "Gene", y = "Priority Score", fill = "Druggability") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

    ggsave("DrugTarget_PriorityPlot.png", plot = p1, width = 11, height = 8, dpi = 300)
    message("    Priority plot saved.")
  }

  # Plot 2: Phase distribution
  if (nrow(drug_df) > 0) {
    phase_count <- drug_df %>%
      mutate(phase_label = recode(as.character(phase),
        "1"="Phase I","2"="Phase II","3"="Phase III","4"="Approved",.default="Unknown")) %>%
      count(phase_label) %>%
      mutate(phase_label = factor(phase_label,
        levels = c("Phase I","Phase II","Phase III","Approved","Unknown")))

    p2 <- ggplot(phase_count, aes(x = phase_label, y = n, fill = phase_label)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
      scale_fill_manual(values = c(
        "Phase I"="#e74c3c","Phase II"="#f39c12",
        "Phase III"="#3498db","Approved"="#2ecc71","Unknown"="#bdc3c7"), guide="none") +
      labs(title = "Drug Pipeline Distribution",
           subtitle = "AD GWAS candidate genes (OpenTargets)",
           x = "Clinical Phase", y = "Number of Drug-Gene Interactions") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave("DrugTarget_PhasePlot.png", plot = p2, width = 8, height = 5, dpi = 300)
    message("    Phase plot saved.")
  }

  # Plot 3: Bubble chart
  if (nrow(drug_summary) > 0) {
    bub <- drug_summary %>% mutate(max_phase = as.numeric(max_phase))

    p3 <- ggplot(bub, aes(x = reorder(symbol, max_phase), y = max_phase,
                          size = n_drugs, fill = druggability)) +
      geom_point(shape = 21, alpha = 0.85, color = "white") +
      geom_text(aes(label = symbol), vjust = -1.3, size = 3.2, fontface = "bold") +
      scale_size_continuous(range = c(6, 16), name = "# Drug interactions") +
      scale_fill_manual(values = c(
        "Approved Drug Exists"="#2ecc71","Phase III"="#3498db",
        "Phase II"="#f39c12","Early Stage"="#e74c3c")) +
      scale_y_continuous(breaks=1:4,
        labels=c("Phase I","Phase II","Phase III","Approved"), limits=c(0.5,4.8)) +
      coord_flip() +
      labs(title    = "Drug Development Stage by Gene",
           subtitle = "Bubble size = number of drug interactions | OpenTargets",
           x = "Gene", y = "Maximum Clinical Phase", fill = "Druggability") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave("DrugTarget_BubblePlot.png", plot = p3, width = 10, height = 7, dpi = 300)
    message("    Bubble chart saved.")
  }
}, error = function(e) message("    Visualisation error: ", e$message))

# --- 7. SUMMARY ---
message("\n>>> Step 5: Summary Report...")

approved_genes <- if("druggability" %in% colnames(priority_df))
  priority_df$symbol[priority_df$druggability == "Approved Drug Exists"] else character()
phase3_genes <- if("druggability" %in% colnames(priority_df))
  priority_df$symbol[priority_df$druggability == "Phase III"] else character()

lines <- c(
  "========================================",
  "  DRUG TARGET PRIORITIZATION REPORT",
  paste0("  Disease: Alzheimer's Disease (", ad_disease_id, ")"),
  paste0("  Date: ", Sys.Date()),
  "----------------------------------------",
  paste0("  Candidate genes queried:      ", length(candidate_genes)),
  paste0("  Genes with AD scores:         ", nrow(assoc_df)),
  paste0("  Genes with drug interactions: ", nrow(drug_summary)),
  paste0("  Total drug interactions:      ", nrow(drug_df)),
  "----------------------------------------",
  paste0("  Genes with APPROVED drugs: ", length(approved_genes)),
  if(length(approved_genes)>0) paste0("    -> ", paste(approved_genes, collapse=", ")),
  paste0("  Genes in Phase III trials: ", length(phase3_genes)),
  if(length(phase3_genes)>0) paste0("    -> ", paste(phase3_genes, collapse=", ")),
  "========================================"
)
writeLines(lines, "DrugTarget_Summary.txt")
cat(paste(lines, collapse="\n"), "\n")
