setwd('/home/sirlibzy/KAREN/final analysis/')


# 0. Install and Load Libraries

# Load libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(STRINGdb)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(vsn)
library(dplyr)
library(tibble)
library(readr)
library(grid)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(apeglm)

# --- Configuration ---
base_dir <- getwd()

file_paths <- c(
  "GSE52194" = "GSE52194.tsv",
  "GSE63582" = "GSE63582.tsv",
  "GSE142731" = "GSE142731.tsv",
  "GSE171957" = "GSE171957.tsv",
  "GSE206998" = "GSE206998.tsv"
)

metadata_files <- c(
  "GSE52194" = "GSE52194_metadata.tsv",
  "GSE63582" = "GSE63582_metadata.tsv",
  "GSE142731" = "GSE142731_metadata.tsv",
  "GSE171957" = "GSE171957_metadata.tsv",
  "GSE206998" = "GSE206998_metadata.tsv"
)

output_dir <- file.path(base_dir, "analysis_results_COMBINED")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
analysis_name <- "CombinedAnalysis_TNBC_vs_Normal"

# --- Set Global Plotting Themes ---

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "white", colour = "black"),
                  plot.background = element_rect(fill = "white", colour = NA),
                  strip.background = element_rect(fill="grey90")))

# --- Helper function to find the condition column ---
find_condition_column <- function(df, conditions_to_find) {
  for (col_name in names(df)) {
    if(!is.null(df[[col_name]]) && !all(is.na(df[[col_name]]))) {
      if (any(conditions_to_find %in% as.character(df[[col_name]]))) {
        return(col_name)
      }
    }
  }
  return(NULL)
}


# --- Data Loading and Initial Preparation Loop ---
all_counts_list <- list()
all_metadata_list <- list()

message("--- Starting Data Loading and Initial Preparation for Each Dataset ---")

for (dataset_name in names(file_paths)) {
  message(paste("\nProcessing Dataset for combination:", dataset_name))
  
  counts_file <- file.path(base_dir, file_paths[dataset_name])
  metadata_file <- file.path(base_dir, metadata_files[dataset_name])
  
  if (!file.exists(counts_file)) {
    message(paste("Counts file not found for", dataset_name, ":", counts_file, ". Skipping."))
    next
  }
  if (!file.exists(metadata_file)) {
    message(paste("Metadata file not found for", dataset_name, ":", metadata_file, ". Skipping."))
    next
  }
  
  raw_counts_df_full <- tryCatch({
    read_tsv(counts_file, col_types = cols(), .name_repair = "minimal")
  }, error = function(e1) {
    tryCatch({
      read.delim(counts_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    }, error = function(e2){
      message(paste("Failed to read counts file for", dataset_name, "with both read_tsv and read.delim. Error 1:", e1$message, "Error 2:", e2$message))
      return(NULL)
    })
  })
  
  if(is.null(raw_counts_df_full) || ncol(raw_counts_df_full) < 2) {
    message(paste("Counts file for", dataset_name, "is empty or has insufficient columns. Skipping."))
    next
  }
  
  gene_ids_col_name <- names(raw_counts_df_full)[1]
  gene_ids <- as.character(raw_counts_df_full[[gene_ids_col_name]])
  
  if(any(duplicated(gene_ids))){
    message(paste("Warning: Duplicate gene IDs found in", dataset_name, ". Making them unique. First duplicate:", gene_ids[duplicated(gene_ids)][1]))
    gene_ids <- make.unique(gene_ids, sep = "_dup")
  }
  
  valid_gene_rows <- !is.na(gene_ids) & gene_ids != ""
  if(sum(valid_gene_rows) == 0) {
    message(paste("No valid gene IDs found in counts file for", dataset_name, ". Skipping."))
    next
  }
  gene_ids <- gene_ids[valid_gene_rows]
  raw_counts_matrix_data <- raw_counts_df_full[valid_gene_rows, -1, drop = FALSE]
  
  if(ncol(raw_counts_matrix_data) == 0) {
    message(paste("No sample columns found in counts file for", dataset_name, "after selecting gene IDs. Skipping."))
    next
  }
  
  original_colnames_counts <- trimws(colnames(raw_counts_matrix_data))
  
  raw_counts_matrix_data_numeric <- suppressWarnings(
    as.data.frame(lapply(raw_counts_matrix_data, function(x) as.numeric(as.character(x))))
  )
  colnames(raw_counts_matrix_data_numeric) <- original_colnames_counts
  
  rownames(raw_counts_matrix_data_numeric) <- gene_ids
  
  raw_counts_matrix_data_numeric[is.na(raw_counts_matrix_data_numeric)] <- 0
  raw_counts_df <- round(raw_counts_matrix_data_numeric)
  
  raw_counts_df <- raw_counts_df[rowSums(raw_counts_df) > 0, , drop=FALSE]
  if (nrow(raw_counts_df) == 0 || ncol(raw_counts_df) == 0) {
    message(paste("No valid gene count data after processing for", dataset_name,". Skipping."))
    next
  }
  message(paste("Initial count matrix dims for", dataset_name, ":", paste(dim(raw_counts_df), collapse="x")))
  message(paste("Sample names from count file header for", dataset_name, "(first 5, trimmed):", paste(head(colnames(raw_counts_df), 5), collapse=", ")))
  
  metadata_df <- tryCatch({
    read_tsv(metadata_file, col_types = cols())
  }, error = function(e) {
    read.delim(metadata_file, sep = "\t", stringsAsFactors = FALSE)
  })
  if(nrow(metadata_df) == 0) {
    message(paste("Metadata file for", dataset_name, "is empty. Skipping."))
    next
  }
  message(paste("Metadata loaded for", dataset_name, ". Dims:", paste(dim(metadata_df), collapse="x")))
  message(paste("Metadata column names for", dataset_name, ":", paste(names(metadata_df), collapse=", ")))
  
  sample_id_col <- NULL
  
  if (dataset_name == "GSE52194") {
    known_id_col_gse52194 <- "geo_accession" 
    if (known_id_col_gse52194 %in% names(metadata_df)) {
      sample_id_col <- known_id_col_gse52194
      message(paste("Manually assigned '", sample_id_col, "' as sample ID column for ", dataset_name, ".", sep=""))
    } else {
      message(paste("Warning: Manually specified ID column '", known_id_col_gse52194, "' not found for ", dataset_name, ". Proceeding with auto-detection.", sep=""))
    }
  }
  
  if (is.null(sample_id_col)) {
    title_like_cols_to_check <- c("title", "Title", "sample_title", "Sample_Title", "sample title", "Sample Title")
    for(col_name_candidate in title_like_cols_to_check) {
      if (col_name_candidate %in% names(metadata_df)) {
        potential_matches <- sum(trimws(as.character(metadata_df[[col_name_candidate]])) %in% colnames(raw_counts_df))
        if (potential_matches > 0.5 * min(nrow(metadata_df), ncol(raw_counts_df)) && potential_matches > 0) {
          sample_id_col <- col_name_candidate
          message(paste("Auto-selected (title-like) '", sample_id_col, "' for ", dataset_name, " (", potential_matches, " matches).", sep=""))
          break
        }
      }
    }
    
    if (is.null(sample_id_col)) {
      title_col_insensitive <- grep("^title$", names(metadata_df), ignore.case = TRUE, value = TRUE)
      if (length(title_col_insensitive) > 0) {
        for (col_name_candidate in title_col_insensitive) {
          potential_matches <- sum(trimws(as.character(metadata_df[[col_name_candidate]])) %in% colnames(raw_counts_df))
          if (potential_matches > 0.5 * min(nrow(metadata_df), ncol(raw_counts_df)) && potential_matches > 0) {
            sample_id_col <- col_name_candidate
            message(paste("Auto-selected (title-like, case-insensitive) '", sample_id_col, "' for ", dataset_name, " (", potential_matches, " matches).", sep=""))
            break
          }
        }
      }
    }
    
    if (is.null(sample_id_col)) {
      message(paste("Could not find 'title'-like column for", dataset_name, ". Falling back to general ID search."))
      possible_id_cols <- c("Sample_ID", "SampleName", "Run", "GSM", "sample_id", "sample", "id", "geo_accession")
      for(col_name_candidate in possible_id_cols) {
        if (col_name_candidate %in% names(metadata_df)) {
          potential_matches <- sum(trimws(as.character(metadata_df[[col_name_candidate]])) %in% colnames(raw_counts_df))
          if (potential_matches > 0.5 * min(nrow(metadata_df), ncol(raw_counts_df)) && potential_matches > 0) {
            sample_id_col <- col_name_candidate
            message(paste("Auto-selected (general ID) '", sample_id_col, "' for ", dataset_name, " (", potential_matches, " matches).", sep=""))
            break
          }
        }
      }
    }
    
    if (is.null(sample_id_col) && length(names(metadata_df)) > 0) {
      first_col_name <- names(metadata_df)[1]
      potential_matches_first_col <- sum(trimws(as.character(metadata_df[[first_col_name]])) %in% colnames(raw_counts_df))
      if (potential_matches_first_col > 0.5 * min(nrow(metadata_df), ncol(raw_counts_df)) && potential_matches_first_col > 0) {
        sample_id_col <- first_col_name
        message(paste("Auto-selected (first column) '", sample_id_col, "' for ", dataset_name, " (", potential_matches_first_col, " matches).", sep=""))
      }
    }
  }
  
  if (is.null(sample_id_col)) {
    message(paste("CRITICAL: Could not determine sample ID column in metadata for", dataset_name, ". Skipping dataset."))
    message(paste("Count sample names (first 5, trimmed):", paste(head(colnames(raw_counts_df),5), collapse=", ")))
    if(length(names(metadata_df)) > 0) {
      message(paste("Metadata values in first column '", names(metadata_df)[1],"' (first 5, trimmed):", paste(head(trimws(as.character(metadata_df[[names(metadata_df)[1]]]))),5, collapse=", ")))
    }
    message("Metadata preview (first 5 rows) for problematic dataset:")
    print(head(metadata_df))
    next
  }
  
  metadata_df[[sample_id_col]] <- trimws(as.character(metadata_df[[sample_id_col]]))
  
  conditions_to_look_for <- c("TNBC tumor", "Normal")
  condition_col_name <- find_condition_column(metadata_df, conditions_to_look_for)
  
  if (is.null(condition_col_name)) {
    message(paste("Could not find a column containing 'TNBC tumor' or 'Normal' conditions in", dataset_name, ". Skipping."))
    message("Available columns and their unique values (first few):")
    for(col_iter in names(metadata_df)){
      if(!is.null(metadata_df[[col_iter]])) {
        message(paste("  Col '", col_iter, "': ", paste(head(unique(as.character(metadata_df[[col_iter]]))), collapse=", ")))
      }
    }
    next
  }
  message(paste("Using column '", condition_col_name, "' for conditions in", dataset_name))
  
  metadata_df_filtered_condition <- metadata_df[which(as.character(metadata_df[[condition_col_name]]) %in% conditions_to_look_for), ]
  
  if (nrow(metadata_df_filtered_condition) == 0) {
    message(paste("No samples with specified conditions ('TNBC tumor', 'Normal') in column '", condition_col_name, "' for ", dataset_name, ". Skipping.", sep=""))
    next
  }
  
  metadata_df_filtered_condition$condition <- factor(metadata_df_filtered_condition[[condition_col_name]], levels = c("Normal", "TNBC tumor"))
  
  if (any(duplicated(metadata_df_filtered_condition[[sample_id_col]]))) {
    message(paste("Warning: Duplicate sample IDs found in chosen metadata column '", sample_id_col, "' for dataset ", dataset_name,
                  " AFTER condition filtering. Making them unique. First duplicate:",
                  metadata_df_filtered_condition[[sample_id_col]][duplicated(metadata_df_filtered_condition[[sample_id_col]])][1]))
    metadata_df_filtered_condition[[sample_id_col]] <- make.unique(metadata_df_filtered_condition[[sample_id_col]], sep="_metaDup")
  }
  
  common_samples_in_dataset <- intersect(colnames(raw_counts_df), metadata_df_filtered_condition[[sample_id_col]])
  
  if (length(common_samples_in_dataset) < 1) {
    message(paste("No common samples found between (trimmed) counts and (trimmed) metadata for", dataset_name, "using metadata column '", sample_id_col, "'. Skipping."))
    message(paste("Count sample names (first 5, trimmed):", paste(head(colnames(raw_counts_df),5), collapse=", ")))
    message(paste("Metadata sample IDs in '", sample_id_col, "' after condition filtering (first 5, trimmed):", paste(head(metadata_df_filtered_condition[[sample_id_col]],5), collapse=", ")))
    next
  }
  message(paste("Found", length(common_samples_in_dataset), "common samples for", dataset_name, ". First 5:", paste(head(common_samples_in_dataset), collapse=", ")))
  
  counts_filtered_ds <- raw_counts_df[, common_samples_in_dataset, drop=FALSE]
  
  metadata_intermediate_ds <- metadata_df_filtered_condition[metadata_df_filtered_condition[[sample_id_col]] %in% common_samples_in_dataset, , drop=FALSE]
  metadata_intermediate_ds <- as.data.frame(metadata_intermediate_ds)
  rownames(metadata_intermediate_ds) <- metadata_intermediate_ds[[sample_id_col]]
  
  metadata_final_ds <- metadata_intermediate_ds[common_samples_in_dataset, , drop=FALSE]
  
  if(!identical(colnames(counts_filtered_ds), rownames(metadata_final_ds))) {
    message(paste("CRITICAL LOCAL MISMATCH for", dataset_name, "BEFORE global ID generation."))
    message("Colnames counts_filtered_ds (first 5): ", paste(head(colnames(counts_filtered_ds),5), collapse=", "))
    message("Rownames metadata_final_ds (first 5): ", paste(head(rownames(metadata_final_ds),5), collapse=", "))
    stop(paste("Stopping due to local sample ID mismatch in dataset:", dataset_name))
  }
  
  if (ncol(counts_filtered_ds) == 0) {
    message(paste("No samples remaining after matching for", dataset_name, ". Skipping."))
    next
  }
  
  globally_unique_sample_ids <- paste(dataset_name, colnames(counts_filtered_ds), sep = "_")
  
  colnames(counts_filtered_ds) <- globally_unique_sample_ids
  rownames(metadata_final_ds) <- globally_unique_sample_ids
  
  metadata_final_ds$UniqueID <- globally_unique_sample_ids
  metadata_final_ds$dataset <- dataset_name
  
  all_counts_list[[dataset_name]] <- counts_filtered_ds
  all_metadata_list[[dataset_name]] <- metadata_final_ds
  message(paste("Successfully processed and stored data for", dataset_name, ". Final sample count for this dataset:", ncol(counts_filtered_ds)))
}

message("--- Finished Data Loading and Initial Preparation ---")

# --- Combine Datasets ---
message("\n--- Combining Datasets ---")

if (length(all_counts_list) == 0) {
  stop("No datasets were successfully loaded and prepared. Exiting.")
}
if (length(all_counts_list) < 2 && length(file_paths) >1 && length(file_paths) > length(all_counts_list) ) {
  message(paste("Warning: Only", length(all_counts_list), "out of", length(file_paths) ,"datasets were successfully processed. Combined analysis will proceed with available data."))
} else if (length(all_counts_list) < 2 && length(file_paths) > 1) {
  message("Warning: Only one dataset was successfully processed. Combined analysis will effectively be a single dataset analysis.")
}

if (length(all_counts_list) > 1) {
  common_genes <- Reduce(intersect, lapply(all_counts_list, rownames))
} else if (length(all_counts_list) == 1) {
  common_genes <- rownames(all_counts_list[[1]])
} else {
  stop("No data in all_counts_list to process after initial loading.")
}

if(length(common_genes) == 0) {
  stop("No common genes found across all successfully processed datasets. Cannot proceed with combined analysis.")
}
message(paste("Found", length(common_genes), "common genes across all datasets for combined analysis."))

combined_counts_list_subset <- lapply(all_counts_list, function(df) df[common_genes, , drop = FALSE])
full_counts_matrix <- do.call(cbind, combined_counts_list_subset)

essential_meta_cols <- c("condition", "dataset", "UniqueID")
combined_metadata_list_for_rbind <- lapply(all_metadata_list, function(df) {
  for(col_iter in essential_meta_cols){
    if(!col_iter %in% names(df)){
      df[[col_iter]] <- NA
    }
  }
  return(df[, essential_meta_cols, drop = FALSE])
})
combined_metadata_df <- do.call(rbind, combined_metadata_list_for_rbind)

if(!all(colnames(full_counts_matrix) %in% rownames(combined_metadata_df))) {
  missing_in_meta <- colnames(full_counts_matrix)[!colnames(full_counts_matrix) %in% rownames(combined_metadata_df)]
  missing_in_counts <- rownames(combined_metadata_df)[!rownames(combined_metadata_df) %in% colnames(full_counts_matrix)]
  message("ERROR: Mismatch detected BEFORE reordering combined metadata.")
  message(paste("Globally unique IDs in counts but NOT in metadata rownames (first 5):", paste(head(missing_in_meta,5), collapse=", ")))
  message(paste("Globally unique IDs in metadata rownames but NOT in counts (first 5):", paste(head(missing_in_counts,5), collapse=", ")))
  stop("Critical Error: Mismatch between sample names in combined counts and combined metadata. Check sample ID uniqueness and prepending logic during per-dataset processing.")
}
combined_metadata_df <- combined_metadata_df[colnames(full_counts_matrix), ]

condition_counts_combined <- table(combined_metadata_df$condition)
if(any(condition_counts_combined < 2)) {
  stop(paste("Combined dataset does not have at least two replicates for each condition. Counts:",
             paste(names(condition_counts_combined), condition_counts_combined, collapse=", ")))
}
dataset_sample_counts <- table(combined_metadata_df$dataset)
if(any(dataset_sample_counts <= 1) && length(unique(combined_metadata_df$dataset)) > 1){
  message("Warning: Some datasets in the combined analysis have only one sample. This might be problematic if 'dataset' is in the design.")
}

# --- Combined DESeq2 Analysis ---
message("\n--- Starting Combined DESeq2 Analysis ---")

message("Pre-filtering low count genes on combined data...")
keep_genes_combined <- rowSums(full_counts_matrix) >= 10
full_counts_matrix_filtered <- full_counts_matrix[keep_genes_combined, , drop = FALSE]
common_genes_filtered <- rownames(full_counts_matrix_filtered)

if (nrow(full_counts_matrix_filtered) == 0) {
  stop("No genes remaining after low count filter on combined data. Exiting.")
}
message(paste("Dimensions of combined count matrix after filtering:", paste(dim(full_counts_matrix_filtered), collapse=" x ")))

combined_metadata_df$dataset <- factor(combined_metadata_df$dataset)
combined_metadata_df$condition <- factor(combined_metadata_df$condition, levels = c("Normal", "TNBC tumor"))

if (length(unique(combined_metadata_df$dataset)) > 1) {
  design_formula <- ~ dataset + condition
  message("Using design: ~ dataset + condition")
} else {
  design_formula <- ~ condition
  message("Only one dataset source detected in combined data. Using design: ~ condition")
}

dds_combined <- DESeqDataSetFromMatrix(countData = full_counts_matrix_filtered,
                                       colData = combined_metadata_df,
                                       design = design_formula)

message("Running DESeq2 on combined data...")
dds_combined <- DESeq(dds_combined)

message(paste("Extracting results for contrast: TNBC tumor vs Normal"))
res_combined <- results(dds_combined, contrast=c("condition","TNBC tumor","Normal"))

message("Applying LFC shrinkage...")
target_coef_sanitized <- make.names(paste0("condition_", "TNBC tumor", "_vs_", "Normal"))
actual_coef_name_for_shrinkage <- NULL

if(target_coef_sanitized %in% resultsNames(dds_combined)){
  actual_coef_name_for_shrinkage <- target_coef_sanitized
} else {
  possible_coef_names <- resultsNames(dds_combined)
  pattern_positive_level <- paste0("condition.*", gsub("[^A-Za-z0-9_]", ".", make.names("TNBC tumor")))
  
  matches <- grep(pattern_positive_level, possible_coef_names, value = TRUE, ignore.case=TRUE)
  
  if (length(matches) == 1) {
    actual_coef_name_for_shrinkage <- matches[1]
  } else if (length(matches) > 1) {
    if(target_coef_sanitized %in% matches) {
      actual_coef_name_for_shrinkage <- target_coef_sanitized
    } else {
      actual_coef_name_for_shrinkage <- matches[which.min(nchar(matches))]
      message(paste("Warning: Multiple potential coefficients found for LFC shrinkage:", paste(matches, collapse=", "), ". Using shortest match:", actual_coef_name_for_shrinkage))
    }
  } else {
    actual_coef_name_for_shrinkage <- resultsNames(dds_combined)[length(resultsNames(dds_combined))]
    message(paste("Warning: Could not automatically determine a specific coefficient for LFC shrinkage. Using the last coefficient:", actual_coef_name_for_shrinkage, "as a fallback."))
  }
}
message(paste("Using coefficient '", actual_coef_name_for_shrinkage, "' for LFC shrinkage with apeglm.", sep=""))

tryCatch({
  res_combined_shrunk <- lfcShrink(dds_combined,
                                   coef = actual_coef_name_for_shrinkage,
                                   type="apeglm",
                                   res = res_combined)
  res_combined <- res_combined_shrunk
}, error = function(e_shrink) {
  message(paste("LFC shrinkage with apeglm and coef '", actual_coef_name_for_shrinkage, "' failed. Error:", e_shrink$message))
  message("Attempting LFC shrinkage with type='normal'.")
  tryCatch({
    res_combined_shrunk_normal <- lfcShrink(dds_combined,
                                            coef = actual_coef_name_for_shrinkage,
                                            type="normal",
                                            res = res_combined)
    res_combined <- res_combined_shrunk_normal
    message("Successfully applied LFC shrinkage with type='normal'.")
  }, error = function(e_normal){
    message(paste("LFC shrinkage with type='normal' also failed. Error:", e_normal$message))
    message("Proceeding with unshrunken LFC values from results().")
  })
})

message("Mapping gene IDs to Symbols for combined results...")
original_ids_combined <- rownames(res_combined)

is_symbol_like_combined <- all(grepl("^[A-Za-z][A-Za-z0-9._-]*$", original_ids_combined[1:min(10, nrow(res_combined))])) &&
  !all(grepl("^[0-9]+(_[a-zA-Z0-9]+)*$", original_ids_combined[1:min(10, nrow(res_combined))]))

if (!is_symbol_like_combined) {
  id_type_guess_combined <- "ENTREZID"
  if (any(startsWith(rownames(res_combined)[1:min(10,nrow(res_combined))], "ENSG"))) id_type_guess_combined <- "ENSEMBL"
  
  suppressMessages({
    symbols_combined <- mapIds(org.Hs.eg.db,
                               keys = original_ids_combined,
                               column = "SYMBOL",
                               keytype = id_type_guess_combined,
                               multiVals = "first")
  })
  res_combined$symbol <- symbols_combined
  if (sum(is.na(res_combined$symbol)) > 0.8 * nrow(res_combined) && id_type_guess_combined == "ENTREZID") {
    message("Mapping to SYMBOL using ENTREZID resulted in many NAs. Original IDs might be symbols or unmappable.")
    res_combined$symbol <- original_ids_combined
  } else {
    res_combined$symbol[is.na(res_combined$symbol)] <- original_ids_combined[is.na(res_combined$symbol)]
  }
} else {
  message("Row names appear to be gene symbols or similar. Using them directly.")
  res_combined$symbol <- original_ids_combined
}

res_ordered_combined <- as.data.frame(res_combined[order(res_combined$padj), ])
# This CSV contains the data for the "Top 10 DEGs table" (user can sort/filter)
write.csv(res_ordered_combined, file.path(output_dir, paste0(analysis_name, "_differential_expression_results.csv")))


# --- Define Top N DEGs for subsequent PPI and GO/KEGG ----
top_n_for_analysis <- 10 # Define the number of top genes you want to focus on (e.g., for the specific heatmap)
top_n_degs_df <- NULL
actual_n_selected_for_analysis <- 0

# Ensure sig_genes_df_comb is defined based on your criteria (e.g., padj < 0.05, |LFC| > 1)
if (!exists("sig_genes_df_comb") || is.null(sig_genes_df_comb)) {
  if (exists("res_ordered_combined") && nrow(res_ordered_combined) > 0) {
    sig_genes_df_comb_temp <- subset(res_ordered_combined, padj < 0.05 & abs(log2FoldChange) > 1)
    if (nrow(sig_genes_df_comb_temp) == 0) {
      message("No significant DEGs (padj < 0.05 & |LFC| > 1) found in res_ordered_combined for defining top N.")
    }
  } else {
    sig_genes_df_comb_temp <- data.frame() 
    message("res_ordered_combined not available or empty for defining significant DEGs for top N.")
  }
} else { # This 'else' might not be hit if sig_genes_df_comb is defined later
  sig_genes_df_comb_temp <- subset(res_ordered_combined, padj < 0.05 & abs(log2FoldChange) > 1) # Re-derive to be sure
}
# If sig_genes_df_comb_temp wasn't created above, try one more time (e.g. if res_ordered_combined just became available)
if (!exists("sig_genes_df_comb_temp")) {
  if (exists("res_ordered_combined") && nrow(res_ordered_combined) > 0) {
    sig_genes_df_comb_temp <- subset(res_ordered_combined, padj < 0.05 & abs(log2FoldChange) > 1)
  } else {
    sig_genes_df_comb_temp <- data.frame()
  }
}


if (exists("sig_genes_df_comb_temp") && nrow(sig_genes_df_comb_temp) > 0) {
  message(paste("Selecting up to Top", top_n_for_analysis, "DEGs for focused analyses and specific heatmap..."))
  sig_genes_df_ordered_for_top_n <- sig_genes_df_comb_temp[order(sig_genes_df_comb_temp$padj, -abs(sig_genes_df_comb_temp$log2FoldChange)), ]
  
  actual_n_selected_for_analysis <- min(nrow(sig_genes_df_ordered_for_top_n), top_n_for_analysis)
  
  if (nrow(sig_genes_df_ordered_for_top_n) == 0) {
    message(paste("No significant DEGs were actually available for selection of top", top_n_for_analysis))
    actual_n_selected_for_analysis <- 0 
  } else if (actual_n_selected_for_analysis < top_n_for_analysis && actual_n_selected_for_analysis > 0) {
    message(paste("Fewer than", top_n_for_analysis, "DEGs found (", actual_n_selected_for_analysis,
                  "). Using all", actual_n_selected_for_analysis ,"available significant DEGs for focused analyses and specific heatmap."))
  } else if (actual_n_selected_for_analysis > 0) {
    message(paste("Using top", actual_n_selected_for_analysis, "DEGs for focused analyses and specific heatmap."))
  }
  
  if (actual_n_selected_for_analysis > 0) {
    top_n_degs_df <- head(sig_genes_df_ordered_for_top_n, actual_n_selected_for_analysis)
    message(paste("Selected", nrow(top_n_degs_df), "genes for top_n_degs_df. First few gene symbols (if available):"))
    print(head(top_n_degs_df$symbol, 5))
  } else {
    message(paste("No DEGs met the criteria to be selected as 'top DEGs'. Specific heatmap and focused analyses might be skipped or reflect no genes."))
    top_n_degs_df <- data.frame() 
    actual_n_selected_for_analysis <- 0
  }
  
} else {
  message(paste("No significant DEGs found (padj < 0.05 & |LFC| > 1). Specific heatmap and focused analyses for top genes will be skipped or reflect no genes."))
  top_n_degs_df <- data.frame() 
  actual_n_selected_for_analysis <- 0
}


# --- Quality Control and Visualization Plots (Combined Data) ---
message("\n--- Generating QC and Visualization plots for Combined Data ---")
vsd_combined <- NULL
tryCatch({
  vsd_combined <- vst(dds_combined, blind = TRUE)
}, error = function(e_vst){
  message(paste("Error during VST transformation:", e_vst$message))
  message("Attempting rlog transformation instead.")
  tryCatch({
    vsd_combined <- rlog(dds_combined, blind = TRUE)
  }, error = function(e_rlog){
    message(paste("Error during rlog transformation:", e_rlog$message))
    message("Proceeding with plots using normalized counts where possible, but some VST/rlog dependent plots might fail or be inaccurate.")
  })
})


if(interactive()) {
  par(mar=c(max(5, ncol(dds_combined) %/% 10), 4, 4, 2) + 0.1)
  boxplot(log2(counts(dds_combined, normalized=FALSE) + 1),
          main=paste("Raw Counts (log2+1) -", analysis_name),
          las=2, col=rainbow(ncol(counts(dds_combined))), ylab="Log2(Raw Count + 1)",
          frame.plot=TRUE, border="black", cex.axis=0.7)
  par(mar=c(5, 4, 4, 2) + 0.1)
}
png(file.path(output_dir, paste0(analysis_name, "_boxplot_raw_counts.png")), width=max(1000, ncol(dds_combined)*20), height=700, bg="white")
par(mar=c(max(8, ncol(dds_combined) %/% 10), 4, 4, 2) + 0.1)
boxplot(log2(counts(dds_combined, normalized=FALSE) + 1),
        main=paste("Raw Counts (log2+1) -", analysis_name),
        las=2, col=rainbow(ncol(counts(dds_combined))), ylab="Log2(Raw Count + 1)",
        frame.plot=TRUE, border="black", cex.axis=0.6)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

if(interactive()){
  par(mar=c(max(5, ncol(dds_combined) %/% 10), 4, 4, 2) + 0.1)
  boxplot(log2(counts(dds_combined, normalized=TRUE) + 1),
          main=paste("Normalized Counts (log2+1) -", analysis_name),
          las=2, col=rainbow(ncol(counts(dds_combined))), ylab="Log2(Normalized Count + 1)",
          frame.plot=TRUE, border="black", cex.axis=0.7)
  par(mar=c(5, 4, 4, 2) + 0.1)
}
png(file.path(output_dir, paste0(analysis_name, "_boxplot_normalized_counts.png")), width=max(1000, ncol(dds_combined)*20), height=700, bg="white")
par(mar=c(max(8, ncol(dds_combined) %/% 10), 4, 4, 2) + 0.1)

boxplot(log2(counts(dds_combined, normalized=TRUE) + 1),
        main=paste("Normalized Counts (log2+1) -", analysis_name),
        las=2, col=rainbow(ncol(counts(dds_combined))), ylab="Log2(Normalized Count + 1)",
        frame.plot=TRUE, border="black", cex.axis=0.6)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

df_norm_counts_comb <- as.data.frame(log2(counts(dds_combined, normalized=TRUE) + 1))
df_norm_counts_comb$gene <- rownames(df_norm_counts_comb)
df_norm_counts_long_comb <- tidyr::pivot_longer(df_norm_counts_comb,
                                                cols = -gene,
                                                names_to = "sample",
                                                values_to = "log2_norm_count")
col_data_for_merge <- as.data.frame(colData(dds_combined))
col_data_for_merge$UniqueID_rownames <- rownames(col_data_for_merge)

df_norm_counts_long_comb <- merge(df_norm_counts_long_comb,
                                  col_data_for_merge[,c("condition", "dataset", "UniqueID_rownames")],
                                  by.x="sample", by.y="UniqueID_rownames")

p_box_norm_gg_comb <- ggplot(df_norm_counts_long_comb, aes(x=sample, y=log2_norm_count, fill=condition)) +
  geom_boxplot() +
  facet_wrap(~dataset, scales="free_x", ncol = min(3, length(unique(df_norm_counts_long_comb$dataset)))) +
  labs(title=paste("Normalized Counts per Sample -", analysis_name),
       x="Sample (Grouped by Dataset)", y="Log2(Normalized Count + 1)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=max(3, 8 - length(unique(df_norm_counts_long_comb$sample)) %/% 50 )))
if(interactive()) print(p_box_norm_gg_comb)
ggsave(file.path(output_dir, paste0(analysis_name, "_boxplot_normalized_ggplot.png")), plot=p_box_norm_gg_comb, width=16, height=9)


if(!is.null(vsd_combined)){
  pca_plot_comb <- plotPCA(vsd_combined, intgroup = c("condition", "dataset")) +
    ggtitle(paste("PCA Plot (VST/rlog transformed) -", analysis_name)) +
    geom_text(aes(label=name), nudge_y = 0.5, size=2, check_overlap = TRUE)
  if(interactive()) print(pca_plot_comb)
  ggsave(file.path(output_dir, paste0(analysis_name, "_PCA_plot_condition_dataset.png")), plot=pca_plot_comb, width=10, height=8)
} else {
  message("PCA plot skipped as VST/rlog transformation failed.")
}

# --- Sample-to-Sample Heatmap (using top 10 most variable genes) ---
if(!is.null(vsd_combined)){
  rv <- rowVars(assay(vsd_combined))
  select_top_variable <- order(rv, decreasing=TRUE)[1:min(10, length(rv))] 
  
  top_variable_assay_data <- assay(vsd_combined)[select_top_variable,]
  
  sample_dists_top_var_comb <- dist(t(top_variable_assay_data))
  sample_dist_matrix_top_var_comb <- as.matrix(sample_dists_top_var_comb)
  
  colors_s2s <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  annot_s2s_comb <- as.data.frame(colData(vsd_combined)[, c("condition", "dataset"), drop = FALSE])
  
  pheatmap_s2s_grob_comb <- pheatmap(sample_dist_matrix_top_var_comb,
                                     clustering_distance_rows = sample_dists_top_var_comb,
                                     clustering_distance_cols = sample_dists_top_var_comb,
                                     col = colors_s2s,
                                     annotation_col = annot_s2s_comb,
                                     main = paste("Sample-to-Sample Distances (Top 10 Var. Genes) -", analysis_name),
                                     border_color = NA,
                                     fontsize_row=max(5, 12 - ncol(dds_combined) %/% 20), 
                                     fontsize_col=max(5, 12 - ncol(dds_combined) %/% 20),
                                     silent = TRUE)
  if(interactive()) { grid.newpage(); grid.draw(pheatmap_s2s_grob_comb$gtable) }
  png(file.path(output_dir, paste0(analysis_name, "_sample_to_sample_heatmap_top10var.png")), width=max(800, ncol(top_variable_assay_data)*25), height=max(700, ncol(top_variable_assay_data)*25), bg="white") 
  grid.draw(pheatmap_s2s_grob_comb$gtable)
  dev.off()
} else {
  message("Sample-to-sample heatmap skipped as VST/rlog transformation failed.")
}


max_abs_lfc <- if(nrow(res_combined[!is.na(res_combined$log2FoldChange),]) > 0) max(abs(res_combined$log2FoldChange), na.rm=T) else 5
if(interactive()) {
  plotMA(res_combined, ylim = c(-max_abs_lfc, max_abs_lfc), main = paste("MA Plot -", analysis_name), colNonSig = "gray60", colSig="red3")
  abline(h=c(-1,1), col="dodgerblue", lwd=2, lty=2)
}
png(file.path(output_dir, paste0(analysis_name, "_MA_plot.png")), width=800, height=600, bg="white")
plotMA(res_combined, ylim = c(-max_abs_lfc, max_abs_lfc), main = paste("MA Plot -", analysis_name), colNonSig = "gray60", colSig="red3")
abline(h=c(-1,1), col="dodgerblue", lwd=2, lty=2)
dev.off()


gene_labels_for_volcano_comb <- if("symbol" %in% names(res_ordered_combined) && !all(is.na(res_ordered_combined$symbol))) {
  make.unique(as.character(res_ordered_combined$symbol))
} else {
  rownames(res_ordered_combined)
}
significant_genes_for_label_df <- res_ordered_combined[which(res_ordered_combined$padj < 0.05 & abs(res_ordered_combined$log2FoldChange) > 1), ]
significant_genes_for_label_df <- significant_genes_for_label_df[order(significant_genes_for_label_df$padj, -abs(significant_genes_for_label_df$log2FoldChange)), ]
num_genes_to_label_comb <- min(20, nrow(significant_genes_for_label_df))
select_labs_volcano <- NULL

if (num_genes_to_label_comb > 0) {
  labels_source_df <- head(significant_genes_for_label_df, num_genes_to_label_comb)
  select_labs_volcano <- if(!is.null(labels_source_df$symbol) && !all(is.na(labels_source_df$symbol))) {
    make.unique(as.character(labels_source_df$symbol))
  } else {
    rownames(labels_source_df)
  }
}

volcano_plot_comb <- EnhancedVolcano(res_ordered_combined,
                                     lab = gene_labels_for_volcano_comb,
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = paste('Volcano Plot -', analysis_name),
                                     pCutoff = 0.05,
                                     FCcutoff = 1.0,
                                     pointSize = 2.0,
                                     labSize = 3.0,
                                     selectLab = select_labs_volcano,
                                     drawConnectors = TRUE,
                                     widthConnectors = 0.5,
                                     colConnectors = 'grey50',
                                     legendPosition = 'right',
                                     subtitle = paste0("Significant (p<0.05, |LFC|>1): ", sum(res_ordered_combined$padj < 0.05 & abs(res_ordered_combined$log2FoldChange) > 1, na.rm=TRUE)),
                                     gridlines.major = FALSE,
                                     gridlines.minor = FALSE)
if(interactive()) print(volcano_plot_comb)
ggsave(file.path(output_dir, paste0(analysis_name, "_volcano_plot.png")), plot=volcano_plot_comb, width=10, height=8)


# --- General DEG Heatmaps (Top 10, 20, 50, etc. by p-adj) ---
if(!is.null(vsd_combined)){
  sig_genes_df_comb <- subset(res_ordered_combined, padj < 0.05 & abs(log2FoldChange) > 1)
  sig_genes_df_comb <- sig_genes_df_comb[order(sig_genes_df_comb$padj), ] 
  sig_genes_names_comb <- rownames(sig_genes_df_comb)
  
  if (length(sig_genes_names_comb) > 1) {
    mat_vsd_comb_all_sig <- assay(vsd_combined)[sig_genes_names_comb, , drop=FALSE]
    annot_df_heatmap_comb <- as.data.frame(colData(vsd_combined)[, c("condition", "dataset"), drop = FALSE])
    
    max_genes_for_all_heatmap <- 200
    deg_counts_to_plot <- c(10, 20, 50, min(max_genes_for_all_heatmap, length(sig_genes_names_comb)))
    if(length(sig_genes_names_comb) > max_genes_for_all_heatmap && length(sig_genes_names_comb) > 50){
    } else if (length(sig_genes_names_comb) <= max_genes_for_all_heatmap && length(sig_genes_names_comb) > 50 && length(sig_genes_names_comb) != 50) {
      deg_counts_to_plot <- c(deg_counts_to_plot, length(sig_genes_names_comb))
    }
    deg_counts_to_plot <- unique(deg_counts_to_plot[deg_counts_to_plot <= length(sig_genes_names_comb)])
    deg_counts_to_plot <- deg_counts_to_plot[deg_counts_to_plot > 1]
    deg_counts_to_plot <- sort(deg_counts_to_plot)
    
    for (n_genes_to_plot in deg_counts_to_plot) {
      current_plot_label <- ifelse(n_genes_to_plot == length(sig_genes_names_comb) && n_genes_to_plot <= max_genes_for_all_heatmap,
                                   paste0("All_", n_genes_to_plot),
                                   as.character(n_genes_to_plot))
      
      top_n_sig_original_ids <- head(sig_genes_names_comb, n_genes_to_plot) 
      mat_subset_comb <- mat_vsd_comb_all_sig[top_n_sig_original_ids, , drop=FALSE]
      
      current_symbols_subset <- if("symbol" %in% names(sig_genes_df_comb) && !all(is.na(sig_genes_df_comb[top_n_sig_original_ids, "symbol"]))) {
        make.unique(as.character(sig_genes_df_comb[top_n_sig_original_ids, "symbol"]))
      } else {
        top_n_sig_original_ids
      }
      rownames(mat_subset_comb) <- current_symbols_subset
      
      if (nrow(mat_subset_comb) < 2 || ncol(mat_subset_comb) < 2) {
        message(paste("Skipping heatmap for top", n_genes_to_plot, "DEGs (by p-adj) due to insufficient dimensions."))
        next
      }
      
      plot_title_heatmap <- ifelse(n_genes_to_plot == length(sig_genes_names_comb) && n_genes_to_plot <= max_genes_for_all_heatmap,
                                   paste("Heatmap of All Sig. DEGs (by p-adj, ", nrow(mat_subset_comb), ") - ", analysis_name, sep=""),
                                   paste("Heatmap of Top", n_genes_to_plot, "DEGs (by p-adj) - ", analysis_name, sep=""))
      
      h_cell <- max(5, min(12, 700 / nrow(mat_subset_comb))) 
      png_height <- min(4000, max(800, nrow(mat_subset_comb)*h_cell + 250)) 
      png_width <- min(4000, max(1000, ncol(mat_subset_comb)*15 + 200))    
      
      pheatmap_deg_grob_comb <- pheatmap(mat_subset_comb,
                                         annotation_col = annot_df_heatmap_comb,
                                         scale = "row",
                                         clustering_distance_rows = "correlation",
                                         main = plot_title_heatmap,
                                         fontsize_row = max(4, h_cell - 1), 
                                         fontsize_col = max(4, 8 - ncol(mat_subset_comb) %/% 20),
                                         cellheight = if(nrow(mat_subset_comb) > 150) NA else h_cell, 
                                         border_color = NA,
                                         silent = TRUE)
      if(interactive()) { grid.newpage(); grid.draw(pheatmap_deg_grob_comb$gtable) }
      
      tryCatch({
        png(file.path(output_dir, paste0(analysis_name, "_heatmap_top_", current_plot_label, "_DEGs_by_padj.png")), width=png_width, height=png_height, bg="white", type="cairo") 
        grid.draw(pheatmap_deg_grob_comb$gtable)
        dev.off()
      }, error = function(e_heat) {
        message(paste("Failed to save heatmap for top", current_plot_label, "DEGs (by p-adj). Error:", e_heat$message))
      })
    }
    if(length(sig_genes_names_comb) > max_genes_for_all_heatmap){
      message(paste("Total significant DEGs found:", length(sig_genes_names_comb), ". Heatmap for 'all' DEGs (by p-adj) was capped at", max_genes_for_all_heatmap, "genes for plotting."))
    }
  } else {
    message(paste("Not enough significant DEGs found for general heatmap generation in", analysis_name))
  }
} else {
  message("General DEG heatmaps skipped as VST/rlog transformation failed.")
}


# --- DETAILED HEATMAP: Top N DEGs (as defined for PPI/GO, e.g., Top 10) ---
# This heatmap specifically uses the `top_n_degs_df` gene list.
message(paste0("\n--- Generating Detailed Heatmap for Top ", actual_n_selected_for_analysis, " DEGs (defined by p-adj & LFC) ---"))
if (!is.null(vsd_combined) && !is.null(top_n_degs_df) && nrow(top_n_degs_df) > 0 && actual_n_selected_for_analysis > 0) {
  
  # Genes are already selected in top_n_degs_df (e.g., top 10 by padj then LFC)
  genes_for_specific_heatmap <- rownames(top_n_degs_df) 
  
  if (length(genes_for_specific_heatmap) >= 2) { # Need at least 2 genes for heatmap
    mat_top_specific_degs_heatmap <- assay(vsd_combined)[genes_for_specific_heatmap, , drop = FALSE]
    
    # Get gene symbols for row labels from top_n_degs_df
    symbols_for_specific_heatmap_rows <- top_n_degs_df[genes_for_specific_heatmap, "symbol"]
    # Fallback for any missing symbols (though top_n_degs_df should have them if mapping worked)
    if (any(is.na(symbols_for_specific_heatmap_rows)) || any(symbols_for_specific_heatmap_rows == "" | duplicated(symbols_for_specific_heatmap_rows))) {
      original_ids_for_fallback <- genes_for_specific_heatmap
      symbols_for_specific_heatmap_rows <- ifelse(is.na(symbols_for_specific_heatmap_rows) | symbols_for_specific_heatmap_rows == "", 
                                                  original_ids_for_fallback, 
                                                  symbols_for_specific_heatmap_rows)
      symbols_for_specific_heatmap_rows <- make.unique(as.character(symbols_for_specific_heatmap_rows))
    }
    rownames(mat_top_specific_degs_heatmap) <- symbols_for_specific_heatmap_rows
    
    annot_df_heatmap_specific_top_n <- as.data.frame(colData(vsd_combined)[, c("condition", "dataset"), drop = FALSE])
    
    heatmap_title_specific_top_n <- paste("Heatmap of Top", length(genes_for_specific_heatmap), "DEGs (p-adj & LFC sorted) -", analysis_name)
    
    h_cell_specific <- max(8, min(25, 800 / length(genes_for_specific_heatmap))) 
    fontsize_row_specific <- max(6, min(12, h_cell_specific -1)) # Adjust font size based on cell height
    
    png_height_specific <- max(600, length(genes_for_specific_heatmap) * h_cell_specific + 250) 
    png_width_specific <- max(800, ncol(mat_top_specific_degs_heatmap) * 20 + 300) # Give more width for sample names if possible
    
    # Cap dimensions to avoid overly large images
    png_height_specific <- min(4000, png_height_specific)
    png_width_specific <- min(4000, png_width_specific)
    
    pheatmap_specific_top_n_grob <- pheatmap(mat_top_specific_degs_heatmap,
                                             annotation_col = annot_df_heatmap_specific_top_n,
                                             scale = "row", 
                                             clustering_distance_rows = "correlation",
                                             clustering_distance_cols = "euclidean", 
                                             main = heatmap_title_specific_top_n,
                                             fontsize_row = fontsize_row_specific,
                                             fontsize_col = max(6, 10 - ncol(mat_top_specific_degs_heatmap) %/% 15),
                                             cellheight = if(length(genes_for_specific_heatmap) > 40) NA else h_cell_specific, # auto if many rows
                                             cellwidth = if(ncol(mat_top_specific_degs_heatmap) > 50) NA else max(10, 700 / ncol(mat_top_specific_degs_heatmap)), # auto if many cols
                                             border_color = "grey60",
                                             show_colnames = (ncol(mat_top_specific_degs_heatmap) <= 60), 
                                             silent = TRUE)
    
    if(interactive()) { 
      grid.newpage()
      grid.draw(pheatmap_specific_top_n_grob$gtable) 
    }
    
    png_file_specific_top_n_heatmap <- file.path(output_dir, paste0(analysis_name, "_heatmap_DETAILED_Top", length(genes_for_specific_heatmap), "_DEGs.png"))
    tryCatch({
      png(png_file_specific_top_n_heatmap, width = png_width_specific, height = png_height_specific, bg = "white", type = "cairo")
      grid.draw(pheatmap_specific_top_n_grob$gtable)
      dev.off()
      message(paste("Saved detailed heatmap for Top", length(genes_for_specific_heatmap), "DEGs to:", png_file_specific_top_n_heatmap))
    }, error = function(e_heat_specific) {
      message(paste("Failed to save detailed heatmap for Top DEGs. Error:", e_heat_specific$message))
      message(paste("Attempted plot dimensions (WxH):", png_width_specific, "x", png_height_specific))
    })
    
  } else {
    message(paste("Not enough Top DEGs (found", length(genes_for_specific_heatmap), ", need at least 2) for the detailed heatmap."))
  }
} else {
  if(is.null(vsd_combined)) message("Detailed heatmap for Top DEGs skipped as VST/rlog transformation failed.")
  if(is.null(top_n_degs_df) || nrow(top_n_degs_df) == 0) message("Detailed heatmap for Top DEGs skipped as no top DEGs were defined (top_n_degs_df is NULL or empty).")
  else if (actual_n_selected_for_analysis == 0) message("Detailed heatmap for Top DEGs skipped as actual_n_selected_for_analysis is 0.")
}


# --- STRINGdb PPI ---
if (!is.null(top_n_degs_df) && nrow(top_n_degs_df) > 0 && "symbol" %in% names(top_n_degs_df)) {
  message(paste("\n--- Attempting to generate PPI network for the selected Top", actual_n_selected_for_analysis, "DEGs ---"))
  
  de_symbols_for_ppi <- na.omit(unique(top_n_degs_df$symbol))
  de_symbols_for_ppi <- de_symbols_for_ppi[de_symbols_for_ppi != "" &
                                             !grepl("^LOC[0-9]+$", de_symbols_for_ppi) &
                                             !grepl("^[0-9]+$", de_symbols_for_ppi) &
                                             !grepl("_dup[0-9]*$", de_symbols_for_ppi) &
                                             !grepl("_metaDup", de_symbols_for_ppi)]
  
  if(length(de_symbols_for_ppi) > 1) {
    message(paste("Found", length(de_symbols_for_ppi), "unique, valid symbols from the top", actual_n_selected_for_analysis, "DEGs for STRINGdb mapping."))
    string_db_ppi_top <- STRINGdb$new(version = "11.5", species = 9606,
                                      score_threshold = 400, input_directory = "")
    
    mapped_genes_ppi_top <- string_db_ppi_top$map(data.frame(gene=de_symbols_for_ppi), "gene", removeUnmappedRows = TRUE)
    message(paste(nrow(mapped_genes_ppi_top), "out of", length(de_symbols_for_ppi), "symbols were mapped by STRINGdb."))
    
    if(nrow(mapped_genes_ppi_top) > 1) {
      ppi_plot_title <- paste("PPI Network of Top", nrow(mapped_genes_ppi_top), "Mapped DEGs (from", actual_n_selected_for_analysis, "selected) -", analysis_name)
      png_filename_ppi <- file.path(output_dir, paste0(analysis_name, "_PPI_network_top", actual_n_selected_for_analysis, "_DEGs.png"))
      
      message(paste("Saving PPI network plot for top DEGs to:", png_filename_ppi))
      png(png_filename_ppi, width=1000, height=1000, bg="white", type="cairo")
      string_db_ppi_top$plot_network(mapped_genes_ppi_top$STRING_id)
      title(main = ppi_plot_title, line = 3) # Use line to adjust title position if needed
      dev.off()
      
      if(interactive()) {
        message("Displaying PPI network for top DEGs in interactive session...")
        string_db_ppi_top$plot_network(mapped_genes_ppi_top$STRING_id)
        title(main = ppi_plot_title, line = 3)
      }
      
      tryCatch({
        enrichment_ppi_top_degs <- string_db_ppi_top$get_enrichment(mapped_genes_ppi_top$STRING_id)
        if(nrow(enrichment_ppi_top_degs) > 0) {
          csv_filename_enrich <- file.path(output_dir, paste0(analysis_name, "_PPI_top", actual_n_selected_for_analysis, "_DEGs_enrichment.csv"))
          message(paste("Saving PPI enrichment data for top DEGs to:", csv_filename_enrich))
          write.csv(enrichment_ppi_top_degs, csv_filename_enrich, row.names=FALSE)
        } else {
          message(paste("No enrichment terms found by STRINGdb for Top", actual_n_selected_for_analysis, "DEGs in", analysis_name))
        }
      }, error = function(e) {
        message(paste("Error getting enrichment from STRINGdb for Top", actual_n_selected_for_analysis, "DEGs in", analysis_name, ":", e$message))
      })
    } else {
      message(paste("Not enough (<=1) Top", actual_n_selected_for_analysis, "DEGs mapped to STRING for PPI network. (Mapped:", nrow(mapped_genes_ppi_top), ", from:",length(de_symbols_for_ppi) ," valid symbols). PPI plot and enrichment skipped."))
    }
  } else {
    message(paste("No suitable gene symbols (or <2 symbols) from Top", actual_n_selected_for_analysis, "DEGs for PPI network. (Valid symbols found:", length(de_symbols_for_ppi), "). PPI plot and enrichment skipped."))
  }
} else {
  message(paste("PPI network for Top", actual_n_selected_for_analysis, "DEGs skipped: No suitable top DEGs selected (top_n_degs_df is NULL, empty, or 'symbol' column missing)."))
}


# --- Gene Ontology and Pathway Analysis ---
# These CSVs contain data for "Benjamini-Hochberg adjusted p-value tables"
if (!is.null(top_n_degs_df) && nrow(top_n_degs_df) > 0) {
  message(paste("\n--- Performing GO and KEGG enrichment analysis for the Top", actual_n_selected_for_analysis, "DEGs ---"))
  entrez_ids_for_enrichment_top_n <- NULL
  entrez_ids_universe_comb <- NULL 
  
  all_tested_ids_universe <- rownames(res_combined) 
  is_entrez_like_universe <- all(grepl("^[0-9]+$", na.omit(all_tested_ids_universe)[1:min(10, length(na.omit(all_tested_ids_universe)))]))
  
  if(is_entrez_like_universe) {
    entrez_ids_universe_comb <- unique(na.omit(all_tested_ids_universe))
    message("Universe for enrichment: Using original row names (assumed Entrez IDs) from all tested genes.")
  } else if ("symbol" %in% names(res_ordered_combined)) { 
    symbols_universe <- na.omit(unique(res_ordered_combined$symbol))
    symbols_universe <- symbols_universe[symbols_universe != "" & !grepl("_dup[0-9]*$", symbols_universe) & !grepl("_metaDup", symbols_universe)]
    if(length(symbols_universe) > 0) {
      message("Universe for enrichment: Mapping symbols from all tested genes to Entrez IDs.")
      map_univ_results <- suppressMessages(bitr(symbols_universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE))
      entrez_ids_universe_comb <- unique(na.omit(map_univ_results$ENTREZID))
      message(paste("Universe: Mapped", length(symbols_universe), "symbols to", length(entrez_ids_universe_comb), "unique Entrez IDs."))
    } else {
      message("Universe for enrichment: No valid symbols found in all_tested_genes for Entrez ID mapping.")
    }
  } else {
    message("Universe for enrichment: Could not determine Entrez IDs for the universe of all tested genes. Enrichment might be affected.")
  }
  entrez_ids_universe_comb <- unique(na.omit(entrez_ids_universe_comb))
  
  
  ids_from_top_n_df <- rownames(top_n_degs_df)
  is_entrez_like_top_n_degs <- all(grepl("^[0-9]+$", na.omit(ids_from_top_n_df)[1:min(10, length(na.omit(ids_from_top_n_df)))]))
  
  if (is_entrez_like_top_n_degs) {
    entrez_ids_for_enrichment_top_n <- unique(na.omit(ids_from_top_n_df))
    message(paste("Gene set for enrichment: Using original row names (assumed Entrez IDs) from top", actual_n_selected_for_analysis, "DEGs."))
  } else if ("symbol" %in% names(top_n_degs_df) && !all(is.na(top_n_degs_df$symbol))) {
    symbols_to_map_top_n_degs <- na.omit(unique(top_n_degs_df$symbol))
    symbols_to_map_top_n_degs <- symbols_to_map_top_n_degs[symbols_to_map_top_n_degs != "" & !grepl("_dup[0-9]*$", symbols_to_map_top_n_degs) & !grepl("_metaDup", symbols_to_map_top_n_degs)]
    if(length(symbols_to_map_top_n_degs) > 0) {
      message(paste("Gene set for enrichment: Mapping symbols from top", actual_n_selected_for_analysis, "DEGs to Entrez IDs."))
      map_results_top_n_degs <- suppressMessages(bitr(symbols_to_map_top_n_degs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE))
      entrez_ids_for_enrichment_top_n <- unique(na.omit(map_results_top_n_degs$ENTREZID))
      message(paste("Gene set: Mapped", length(symbols_to_map_top_n_degs), "symbols to", length(entrez_ids_for_enrichment_top_n), "unique Entrez IDs."))
    } else {
      message(paste("Gene set for enrichment: No valid symbols found in top", actual_n_selected_for_analysis, "DEGs for Entrez ID mapping."))
    }
  } else {
    message(paste("Gene set for enrichment: Could not determine Entrez IDs from top", actual_n_selected_for_analysis,
                  "DEGs (neither rownames are Entrez-like nor 'symbol' column available/valid)."))
  }
  entrez_ids_for_enrichment_top_n <- unique(na.omit(entrez_ids_for_enrichment_top_n))
  
  if (length(entrez_ids_for_enrichment_top_n) > 5 && length(entrez_ids_universe_comb) > 10) {
    tryCatch({
      ego_bp_top_n <- enrichGO(gene = entrez_ids_for_enrichment_top_n, universe = entrez_ids_universe_comb,
                               OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
      if (!is.null(ego_bp_top_n) && nrow(as.data.frame(ego_bp_top_n)) > 0) {
        write.csv(as.data.frame(ego_bp_top_n), file.path(output_dir, paste0(analysis_name, "_GO_BP_enrichment_top", actual_n_selected_for_analysis, "_DEGs.csv")), row.names=FALSE)
        if (nrow(as.data.frame(ego_bp_top_n)) > 1) {
          dp_gobp_top_n <- dotplot(ego_bp_top_n, showCategory=min(15, nrow(as.data.frame(ego_bp_top_n)))) +
            ggtitle(paste("GO BP Enrichment - Top", actual_n_selected_for_analysis, "DEGs -", analysis_name))
          if(interactive()) print(dp_gobp_top_n)
          ggsave(file.path(output_dir, paste0(analysis_name, "_GO_BP_dotplot_top", actual_n_selected_for_analysis, "_DEGs.png")), plot=dp_gobp_top_n, width=10, height=min(15, max(6,nrow(as.data.frame(ego_bp_top_n))*0.25 + 2)))
        }
      } else { message(paste("No significant GO BP terms found for Top", actual_n_selected_for_analysis, "DEGs in", analysis_name)) }
      
      ekegg_top_n <- enrichKEGG(gene = entrez_ids_for_enrichment_top_n, universe = entrez_ids_universe_comb,
                                organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1)
      if (!is.null(ekegg_top_n) && nrow(as.data.frame(ekegg_top_n)) > 0) {
        ekegg_readable_top_n <- setReadable(ekegg_top_n, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        write.csv(as.data.frame(ekegg_readable_top_n), file.path(output_dir, paste0(analysis_name, "_KEGG_enrichment_top", actual_n_selected_for_analysis, "_DEGs.csv")), row.names=FALSE)
        if (nrow(as.data.frame(ekegg_readable_top_n)) > 1) {
          dp_kegg_top_n <- dotplot(ekegg_readable_top_n, showCategory=min(15,nrow(as.data.frame(ekegg_readable_top_n)))) +
            ggtitle(paste("KEGG Pathway Enrichment - Top", actual_n_selected_for_analysis, "DEGs -", analysis_name))
          if(interactive()) print(dp_kegg_top_n)
          ggsave(file.path(output_dir, paste0(analysis_name, "_KEGG_dotplot_top", actual_n_selected_for_analysis, "_DEGs.png")), plot=dp_kegg_top_n, width=10, height=min(15, max(6,nrow(as.data.frame(ekegg_readable_top_n))*0.25 + 2)))
        }
      } else { message(paste("No significant KEGG pathways found for Top", actual_n_selected_for_analysis, "DEGs in", analysis_name)) }
    }, error = function(e) { message(paste("Error during clusterProfiler enrichment for Top", actual_n_selected_for_analysis, "DEGs:", e$message)) })
  } else {
    message(paste("Not enough Entrez IDs for enrichment of Top", actual_n_selected_for_analysis, "DEGs (Gene Set with Entrez:",
                  length(entrez_ids_for_enrichment_top_n), ", Universe with Entrez:", length(entrez_ids_universe_comb),"). Minimum 6 genes in set, 11 in universe required."))
  }
} else {
  message(paste("GO/KEGG enrichment for Top DEGs skipped as no top DEGs were selected (top_n_degs_df is NULL or empty)."))
}


# --- Counts plot for top DEGs ---
if (exists("sig_genes_df_comb") && nrow(sig_genes_df_comb) > 0) { 
  sig_genes_df_comb_ordered_for_counts_plot <- sig_genes_df_comb[order(sig_genes_df_comb$padj), ]
  top_degs_for_plot <- head(rownames(sig_genes_df_comb_ordered_for_counts_plot), min(3, nrow(sig_genes_df_comb_ordered_for_counts_plot)))
  
  if(length(top_degs_for_plot) > 0) {
    for (gene_id_comb in top_degs_for_plot) {
      gene_symbol_for_plot_comb <- gene_id_comb
      if(gene_id_comb %in% rownames(res_ordered_combined) && "symbol" %in% names(res_ordered_combined)){
        symbol_val_comb <- res_ordered_combined[gene_id_comb, "symbol"]
        if(!is.na(symbol_val_comb) && symbol_val_comb != "") gene_symbol_for_plot_comb <- symbol_val_comb
      }
      
      plot_data_counts <- plotCounts(dds_combined, gene = gene_id_comb, intgroup = c("condition", "dataset"),
                                     normalized = TRUE, returnData = TRUE)
      
      p_counts <- ggplot(plot_data_counts, aes(x = condition, y = count, color = dataset)) +
        geom_jitter(width = 0.2, size=2, alpha=0.7) +
        geom_boxplot(aes(fill=condition), outlier.shape=NA, alpha=0.3, width=0.5, position=position_dodge(width=0.8)) +
        labs(title = paste("Normalized Counts:", gene_symbol_for_plot_comb),
             subtitle = analysis_name,
             x = "Condition", y = "Normalized Count") +
        scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      
      if(interactive()) print(p_counts)
      
      safe_gene_symbol <- str_replace_all(gene_symbol_for_plot_comb, "[^A-Za-z0-9_.-]", "_")
      ggsave(file.path(output_dir, paste0(analysis_name, "_counts_plot_", safe_gene_symbol, ".png")), plot=p_counts, width=7, height=6)
    }
  } else {
    message("No DEGs available to plot top counts for general DEG overview.")
  }
} else {
  message("No significant DEGs defined (sig_genes_df_comb) to plot top counts for general DEG overview.")
}

# --- Dispersion Plot ---
if(interactive()) plotDispEsts(dds_combined, main=paste("Dispersion Estimates -", analysis_name))
png(file.path(output_dir, paste0(analysis_name, "_dispersion_plot.png")), width=800, height=600, bg="white")
plotDispEsts(dds_combined, main=paste("Dispersion Estimates -", analysis_name))
dev.off()


message(paste("\n\nCombined differential expression analysis complete."))
message("Results and plots saved in:", output_dir)
message("Top 10 DEGs table data is in '", paste0(analysis_name, "_differential_expression_results.csv"), "'.")
message("Benjamini-Hochberg adjusted p-values for GO/KEGG are in the respective enrichment CSV files ('p.adjust' column).")
message("A detailed heatmap for the specifically selected Top ", actual_n_selected_for_analysis, " DEGs (if any) was generated as '", paste0(analysis_name, "_heatmap_DETAILED_Top", actual_n_selected_for_analysis, "_DEGs.png"),"'.")

gc()
