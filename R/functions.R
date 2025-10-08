# ==============================================================================
# Custom R Functions Library
# Author: Maxwell Konnaris
# Last Updated: 2/12/25
# ==============================================================================


# ==== Data Cleaning ===========================================================

cleandata_runpredictions <- function() {
  source("R/MLP.R")
  library(tidyverse)
  library(dplyr) 
  library(ggplot2)
  library(vegan)
  library(tibble)
  library(here)
  library(readr)
  library(tidyr)
  library(readxl)
  library(data.table)
  library(readODS)
  
  convert_files <- function(file_paths, filepaths, output_format = "csv", delimiter = ",") {
    if (length(file_paths) != length(filepaths)) {
      stop("Error: file_paths and filepaths must have the same length!")
    }
    
    for (i in seq_along(file_paths)) {
      file_path <- file_paths[i]
      output_file <- filepaths[i] 
      
      message("Processing: ", file_path)
      if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
        data <- read.delim(file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
      } else if (grepl("\\.(tsv|txt)$", file_path, ignore.case = TRUE)) {
        data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        message("⚠️ Skipping unsupported file: ", file_path)
        next
      }
      dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
      if (output_format == "csv") {
        write.table(data, output_file, sep = delimiter, row.names = FALSE, quote = FALSE, col.names = TRUE)
      } else if (output_format == "tsv") {
        write.table(data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
      }
      
      message("✅ Saved: ", output_file)
    }
  }
  
  
  find_best_match <- function(loadcol, df) {
    normalize_name <- function(name) {
      name <- gsub("[\\.\\/_]", " ", name)  # Replace ., /, _ with spaces
      name <- gsub("^X(\\d)", "\\1", name)  # Remove 'X' before numbers
      name <- tolower(trimws(name))  # Convert to lowercase and trim whitespace
      return(name)
    }
    generate_variations <- function(name) {
      name_variants <- c(
        name,
        gsub("[\\.\\/_]", " ", name),  # Replace ., /, _ with spaces
        gsub("[\\.\\/_]", ".", name),  # Replace with .
        gsub("[\\.\\/_]", "/", name),  # Replace with /
        gsub("[\\.\\/_]", "_", name),  # Replace with _
        gsub("[\\.\\/_]", "", name),   # Remove special characters
        gsub("^X(\\d)", "\\1", name)   # Remove leading "X" before numbers
      )
      name_variants <- unique(c(name_variants, tolower(name_variants), toupper(name_variants)))
      
      return(name_variants)
    }
    normalized_target <- normalize_name(loadcol)
    possible_names <- generate_variations(loadcol)
    actual_colnames <- colnames(df)
    normalized_colnames <- sapply(actual_colnames, normalize_name)
    exact_match <- actual_colnames[normalized_colnames %in% possible_names]
    if (length(exact_match) > 0) {
      return(exact_match[1]) 
    }
    fuzzy_match <- agrep(normalized_target, normalized_colnames, value = TRUE, ignore.case = TRUE)
    if (length(fuzzy_match) > 0) {
      return(actual_colnames[which(normalized_colnames %in% fuzzy_match)][1])  
    }
    
    return(NA) 
  }
  
  
  process_dataset <- function(
    dataset_name,
    file_path,
    old_file_path,
    model,
    inputformat,
    taxatype,
    trueloadfilepath = NA,
    output_dir = "processedandfinalcomparisons/",
    loadcol = NA,
    sampleid = NA,
    delimiter = NULL) {
    
    message("Processing dataset: ", dataset_name)
    
    if (!file.exists(file_path)) {
      stop("❌ Error: The file does not exist: ", file_path)
    }
    raw_data <- read_and_fix_data(file_path)
    
    if (is.null(raw_data) || nrow(raw_data) == 0 || ncol(raw_data) == 0) {
      message("🔄 Retrying with original file: ", old_file_path)
      
      raw_data <- tryCatch({
        read_and_fix_data(old_file_path)
      }, error = function(e) {
        stop("❌ Error: Unable to read both converted and original files.")
      })
    }
    
    if ("X" %in% colnames(raw_data)) {
      message("🔄 Setting row names from column 'X'...")
      rownames(raw_data) <- raw_data$X
      raw_data$X <- NULL
    } else if ("" %in% colnames(raw_data)) {
      message("🔄 Setting row names from an empty column name...")
      colnames(raw_data)[colnames(raw_data) == ""] <- "sampleID"
      rownames(raw_data) <- raw_data$sampleID
      raw_data$sampleID <- NULL
    } else {
      non_numeric_cols <- which(sapply(raw_data, function(col) !is.numeric(col)))
      if (length(non_numeric_cols) > 0) {
        first_non_numeric_col <- colnames(raw_data)[non_numeric_cols[1]]
        
        if (length(non_numeric_cols) > 1) {
          message("⚠️ Multiple non-numeric columns detected. Using '", first_non_numeric_col, "' as row names.")
        } else {
          message("🔄 Setting row names from non-numeric column: ", first_non_numeric_col)
        }
        rownames(raw_data) <- raw_data[[first_non_numeric_col]]
        raw_data[[first_non_numeric_col]] <- NULL  
      } else {
        warning("⚠️ No suitable column found for row names. Keeping default indices.")
      }
    }
    
    original_colnames <- colnames(raw_data)  
    rownamesIDs = rownames(raw_data)
    message("🔍 Checking for non-numeric columns...")
    non_numeric_cols <- sapply(raw_data, function(x) !is.numeric(x))
    
    if (any(non_numeric_cols)) {
      warning("⚠️ Warning: Non-numeric columns detected! Attempting to convert...")
      print(names(raw_data)[non_numeric_cols])
      
      suppressWarnings({
        raw_data[] <- lapply(raw_data, function(x) {
          if (all(sapply(x, is.character))) {
            as.numeric(x)
          } else {
            x
          }
        })
      })
    }

    if (!all(sapply(raw_data, is.numeric))) {
      stop("❌ Error: The dataset contains non-numeric columns after conversion.")
    }

    if (any(is.na(raw_data))) {
      message("⚠️ Warning: Missing values detected. Imputing with 0...")
      raw_data[is.na(raw_data)] <- 0
    }

    if (nrow(raw_data) == 0 || ncol(raw_data) == 0) {
      stop("❌ Error: Data is empty after removing missing values or filtering.")
    }

    if (!is.null(inputformat)) {
      if (inputformat == 'counts') {
        message("🔢 Normalizing count data...")
        raw_data <- as.data.frame(
          t(apply(raw_data, 1, function(x) {
            if (sum(x, na.rm = TRUE) == 0) {
              return(rep(0, length(x)))
            }
            x / sum(x, na.rm = TRUE)
          }))
        )
        colnames(raw_data) <- original_colnames  
      } else if (inputformat == 'relative') {
        message("🔍 Checking if rows sum to 1 for relative abundance data...")
        
        row_sums <- rowSums(raw_data, na.rm = TRUE)
        rows_to_rescale <- which(abs(row_sums - 1) > 1e-6)
        
        if (length(rows_to_rescale) > 0) {
          warning("⚠️ Warning: Some rows do not sum to 1. Rescaling...")
          
          raw_data[rows_to_rescale, ] <- t(apply(raw_data[rows_to_rescale, , drop = FALSE], 1, function(x) {
            if (sum(x, na.rm = TRUE) == 0) {
              return(rep(0, length(x)))
            }
            x / sum(x, na.rm = TRUE)
          }))
          colnames(raw_data) <- original_colnames  
        } else {
          message("All rows already sum to 1. No rescaling needed.")
        }
      }
    }
    
    rownames(raw_data) = rownamesIDs

    if (ncol(raw_data) > 1) {
      message("🔢 Removing taxa with 0% prevalence across all samples...")
      raw_data <- raw_data[, colSums(raw_data != 0) > 0]
    }
    
    if (ncol(raw_data) == 0) {
      stop("Error: All taxa were removed due to zero prevalence!")
    }

    valid_models <- c("16S", "shotgun metagenomics")
    if (!(model %in% valid_models)) {
      stop("Error: Unknown model type: ", model)
    }

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    if (model == "16S") {
      predictiongalaxy <- NULL
      predictionmetacardis <- NULL
      prediction <- MLP(raw_data, "rdp_train_set_16 (16S rRNA)", "load")
      prediction$load <- log10(prediction$load)
      out_file <- file.path(output_dir, paste0(dataset_name, "_predicted16s.csv"))
      write.csv(prediction, out_file, row.names = FALSE)
      message("✅ Predictions written to: ", out_file)
    } else {
      prediction <- NULL
      out_file_g <- file.path(output_dir, paste0(dataset_name, "_G_predictedshotgunmetagenomics.csv"))
      predictiongalaxy <- MLP(raw_data, taxatype, "galaxy", "load")
      predictiongalaxy$load <- log10(predictiongalaxy$load)
      write.csv(predictiongalaxy, out_file_g, row.names = FALSE)
      message("✅ Galaxy Predictions written to: ", out_file_g)
      
      out_file_mc <- file.path(output_dir, paste0(dataset_name, "_MC_predictedshotgunmetagenomics.csv"))
      predictionmetacardis <- MLP(raw_data, taxatype, "metacardis", "load")
      predictionmetacardis$load <- log10(predictionmetacardis$load)
      write.csv(predictionmetacardis, out_file_mc, row.names = FALSE)
      message("✅ MetaCardis Predictions written to: ", out_file_mc)
    }

    measuredload <- NULL
    
    if (!is.null(trueloadfilepath) && file.exists(trueloadfilepath)) {
      message("📂 Reading measured load file: ", trueloadfilepath)
      measuredload <- read_and_fix_data(trueloadfilepath)
      
      best_loadcol <- if (!is.null(measuredload) && loadcol %in% colnames(measuredload)) {
        loadcol
      } else {
        find_best_match(loadcol, measuredload)
      }
      
      if (!is.null(best_loadcol) && best_loadcol %in% colnames(measuredload)) {
        message("✅ Load column found: ", best_loadcol)
        measuredload[[best_loadcol]] <- suppressWarnings(as.numeric(measuredload[[best_loadcol]]))
        
        if (any(is.na(measuredload[[best_loadcol]]))) {
          measuredload <- measuredload[!is.na(measuredload[[best_loadcol]]), ]
        }
        measuredload[[best_loadcol]] <- log10(measuredload[[best_loadcol]])
        
        out_file <- file.path(output_dir, paste0(dataset_name, "_scale.csv"))
        write.csv(measuredload, out_file, row.names = FALSE)
        message("✅ Predictions written to: ", out_file)
      } else {
        stop("❌ Error: Could not find a valid column matching '", loadcol, "' in measured load.")
      }
    }
    
    return(list(
      study = dataset_name,
      predictedgalaxy = predictiongalaxy,
      predictedmetacardis = predictionmetacardis,
      predictedvandeputte = prediction,
      measuredload = measuredload
    ))
  }
  
  
  read_and_fix_data <- function(file_path) {
    message("📂 Reading file: ", file_path)
    
    if (!file.exists(file_path)) {
      stop("❌ Error: File does not exist: ", file_path)
    }
    delimiters <- c(",", "\t", ";", " ")
    raw_data <- NULL
    for (delim in delimiters) {
      message("🔍 Trying delimiter: '", delim, "'")
      temp_data <- tryCatch(
        read.delim(file_path, sep = delim, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
      if (!is.null(temp_data) && ncol(temp_data) > 1) {
        message("✅ Correct delimiter found: ", delim)
        raw_data <- temp_data
        break
      }
    }

    if (!is.null(raw_data) && ncol(raw_data) == 1) {
      message("⚠️ Still only one column. Retrying without headers...")
      temp_data <- read.delim(file_path, sep = delimiter, header = FALSE, stringsAsFactors = FALSE,  check.names = FALSE)
      
      if (!is.null(temp_data) && nrow(temp_data) > 1) {
        message("✅ Using first row as column names.")
        colnames(raw_data) <- as.character(unlist(temp_data[1, ]))
        raw_data <- temp_data[-1, ] 
      }
    }
    
    return(raw_data)
  }

  my_excel_file <- "processedandfinalcomparisons/inputs.ods"
  datasets_info <- read_ods(my_excel_file, sheet=1)
  convert_files(datasets_info$oldfilepath, datasets_info$filepath, output_format = "csv", delimiter = ",")
  results_list <- lapply(seq_len(nrow(datasets_info)),function(i) {
    
    row_i <- datasets_info[i, ]
    
    process_dataset(
      dataset_name    = row_i$study,
      file_path       = row_i$filepath,
      old_file_path   = row_i$oldfilepath,
      model           = row_i$model,
      inputformat     = row_i$inputformat,
      taxatype        = row_i$taxatype,
      trueloadfilepath= row_i$trueloadfilepath,
      loadcol         = row_i$loadcol,
      sampleid        = row_i$sampleID
    )
  })
}

load_microbialload_dataset <- function(datasets) {
  
  datasets$file_path <- file.path(getwd(), datasets$file_path)
  for (i in seq_len(nrow(datasets))) {
    
    dataset <- datasets[i, "dataset_name"]
    filepath <- datasets[i, "file_path"]
    loadcol <- unique(datasets[i, "loadcol"])
    sourcetype <- datasets[i, "source_type"]
    abundance <- datasets[i, "abundance"]
    trained <- datasets[i, "trainingset"]
    sequencingtype <- datasets[i,"sequencingtype"]
    sampleIDcol <- datasets[i,"sampleIDcol"]
    if (tolower(abundance) %in% c("counts", "relative", "metadata")) {
      next
    }
    
    print(paste0("Loading: ",dataset," ", sourcetype))
    
    
    if (!file.exists(filepath)) {
      warning(paste("File not found for dataset:", dataset, "at path:", filepath))
      next
    }
    
    df <- read_and_fix_data(filepath) %>% data.frame(check.names = FALSE)
    
    if ("X" %in% colnames(df)) {
      df <- df[, !colnames(df) %in% "X"]
    }
    
    if (!loadcol %in% colnames(df)) {
      loadcol <- find_best_match(loadcol, df)
      if (!is.na(sampleIDcol) && loadcol %in% colnames(df)) {
        datasets[i,"loadcol"] = loadcol
      } else {
        message(paste("Approximate load ID not found in predicted dataset:", dataset))
        next
      } 
    }
    if (!is.numeric(df[[loadcol]])) {
      df[[loadcol]] <- as.numeric(df[[loadcol]])
    }
    
    if (!sampleIDcol %in% colnames(df)) {
      sampleIDcol <- find_best_match(sampleIDcol, df)
      if (!is.na(sampleIDcol) && sampleIDcol %in% colnames(df)) {
        datasets[i,"sampleIDcol"] = sampleIDcol
      } else {
        message(paste("Approximate sample ID not found in predicted dataset:", dataset))
        next
      } 
    }
    
    if(sourcetype == 'Measured') {
      df <- subset(df, is.finite(df[[loadcol]]))
    }
    
    if (sequencingtype == "shotgun metagenomics") {
      print(paste("Assigning data to variable:", paste0(tolower(gsub(" ", "_", dataset)), "_",tolower(gsub(" ", "_", trained)),"_", sourcetype)))
      assign(paste0(tolower(gsub(" ", "_", dataset)), "_",tolower(gsub(" ", "_", trained)), "_", sourcetype), df, envir = .GlobalEnv)
    } else {
      print(paste("Assigning data to variable:", paste0(tolower(gsub(" ", "_", dataset)), "_", sourcetype)))
      assign(paste0(tolower(gsub(" ", "_", dataset)), "_", sourcetype), df, envir = .GlobalEnv)
    }
    
  }
}

metadata_loading <- function(datasets) {
  
  datasets$file_path <- file.path(getwd(), datasets$file_path)
  for (i in seq_len(nrow(datasets))) {
    
    dataset <- datasets[i, "dataset_name"]
    filepath <- datasets[i, "file_path"]
    abundance <- datasets[i, "abundance"]
    sampleIDcol <- datasets[i,"sampleIDcol"]
    if (tolower(abundance) %in% c("counts", "relative", "total")) {
      next
    }
    
    print(paste0("Loading: ",dataset))
    
    if (!file.exists(filepath)) {
      warning(paste("File not found for dataset:", dataset, "at path:", filepath))
      next
    }
    
    df <- read_and_fix_data(filepath) %>% data.frame(check.names = FALSE)
    
    if ("X" %in% colnames(df)) {
      df <- df[, !colnames(df) %in% "X"]
    }
    
    if (!sampleIDcol %in% colnames(df)) {
      sampleIDcol <- find_best_match(sampleIDcol, df)
      if (!is.na(sampleIDcol) && sampleIDcol %in% colnames(df)) {
        datasets[i,"sampleIDcol"] = sampleIDcol
      } else {
        message(paste("Approximate sample ID not found in predicted dataset:", dataset))
        next
      } 
    }
    assign(paste0(tolower(gsub(" ", "_", dataset)), "_",abundance), df, envir = .GlobalEnv)
  }
}

convert_datasets <- function(datasets, removal = FALSE, raretaxaremovalthreshold = 0.01, prevalence_threshold = 0.1) {
  datasets$file_path <- file.path(getwd(), datasets$file_path)
  
  for (i in seq_len(nrow(datasets))) {
    
    dataset <- datasets[i, "dataset_name"]
    sequencingtype <- datasets[i, "sequencingtype"]
    loadtype <- datasets[i, "loadtype"]
    filepath <- datasets[i, "file_path"]
    loadcol <- datasets[i, "loadcol"]
    sourcetype <- datasets[i, "source_type"]
    abundance <- datasets[i, "abundance"]
    trained <- datasets[i,"trainingset"]
    
    if (tolower(abundance) %in% c("total", "metadata")) {
      next
    }
    
    
    if (!file.exists(filepath)) {
      warning(paste("File not found for dataset:", dataset, "at path:", filepath))
      next
    }
    
    print(paste("Processing dataset:", dataset))
    data = read_and_fix_data(filepath) %>% data.frame(check.names = FALSE)
    data <- na.omit(data)
    possible_id_cols <- c("X", "Var.1", "V1")
    matched_col <- intersect(possible_id_cols, colnames(data))
    if (length(matched_col) > 0) {
      matched_col <- matched_col[1] 
      rownames(data) <- data[[matched_col]]
      data[[matched_col]] <- NULL
      message(paste0("✅ Assigned row names using column: ", matched_col))
    } else if ("" %in% colnames(data)) {
      message("🔄 Setting row names from an empty column name...")
      colnames(data)[colnames(data) == ""] <- "sampleID"
      rownames(data) <- data$sampleID
      data$sampleID <- NULL
    } else{
      message(paste0("⚠️ No sample ID column found for: ", dataset, " ", abundance))
      non_numeric_cols <- colnames(data)[sapply(data, function(col) !is.numeric(col))]
      if (length(non_numeric_cols) > 0) {
        matched_col <- non_numeric_cols[1]
        rownames(data) <- data[[matched_col]]
        data[[matched_col]] <- NULL
        message(paste0("✅ Assigned row names using first non-numeric column: ", matched_col))
      } else {
        for (col in possible_id_cols) {
          matched_col <- find_best_match(col, data)
          if (!is.na(matched_col) && matched_col %in% colnames(data)) {
            rownames(data) <- data[[matched_col]]
            data[[matched_col]] <- NULL
            message(paste0("✅ Assigned row names using best-matched column: ", matched_col))
            break 
          }
        }
        if (is.na(matched_col) || !(matched_col %in% colnames(data))) {
          message(paste0("❌ No suitable sample ID column found for: ", dataset, " ", abundance))
        }
      }
    }
    
    data[] <- lapply(data, as.numeric)
    data <- data[rowSums(data != 0) > 0, ]
    data<- data[, colSums(data != 0) > 0]
    if (abundance == 'counts') {
      print(paste("Count table available, assigning to variable:", paste0(tolower(gsub(" ", "_", dataset)), "_countdata")))
      assign(paste0(tolower(gsub(" ", "_", dataset)), "_countdata"), data, envir = .GlobalEnv)
      data <- data / ifelse(rowSums(data) == 0, 1, rowSums(data))
    }
    sampleIDcol <- unique(datasets[datasets$dataset_name == dataset & 
                                     datasets$source_type == "Measured" & 
                                     datasets$abundance == "total", "sampleIDcol"])
    
    dataset_key <- tolower(gsub(" ", "_", dataset))
    measured_var <- paste0(dataset_key, if (sequencingtype == "shotgun metagenomics") "__Measured" else "_Measured")
    
    measured <- tryCatch(get(measured_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    
    if (!is.null(measured)) {
      if (sampleIDcol %in% colnames(measured)) {
        aligned_sample_ids <- intersect(measured[[sampleIDcol]], rownames(data))
      } else {
        matched_sampleIDcol <- find_best_match(sampleIDcol, measured)
        if (!is.na(matched_sampleIDcol) && matched_sampleIDcol %in% colnames(measured)) {
          aligned_sample_ids <- intersect(measured[[matched_sampleIDcol]], rownames(data))
          sampleIDcol <- matched_sampleIDcol
          message("✅ Data successfully aligned with measured dataset.")
        } else {
          message(paste("⚠️ SampleIDcol:", sampleIDcol, "not found in measured dataset. Skipping alignment."))
        }
      }
      data <- data[rownames(data) %in% aligned_sample_ids, , drop = FALSE]
    } else {
      message("⚠️ No measured dataset found. Proceeding without alignment.")
    }
    
    sampleIDcol2 <- unique(datasets$sampleIDcol[datasets$dataset_name == dataset & datasets$abundance == "metadata"])
    metadata_var <- paste0(dataset_key, "_metadata") 
    metadata <- tryCatch(get(metadata_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    
    if (!is.null(metadata) & !is.null(measured)) {
      if (sampleIDcol2 %in% colnames(metadata)) {
        aligned_sample_ids <- intersect(measured[[sampleIDcol]], metadata[[sampleIDcol2]])
      } else {
        matched_sampleIDcol2 <- find_best_match(sampleIDcol2, metadata)
        if (!is.na(matched_sampleIDcol2) && matched_sampleIDcol2 %in% colnames(metadata)) {
          aligned_sample_ids <- intersect(measured[[sampleIDcol]], metadata[[matched_sampleIDcol2]])
          sampleIDcol2 <- matched_sampleIDcol2
          message("✅ Data successfully aligned with metadata.")
        } else {
          message(paste("⚠️ SampleIDcol:", sampleIDcol2, "not found in metadata dataset. Skipping alignment."))
        }
      }
      metadata <- metadata[metadata[[sampleIDcol2]] %in% aligned_sample_ids, , drop = FALSE]
      assign(metadata_var, metadata, envir = .GlobalEnv)
    } else if (!is.null(metadata)) {
      assign(metadata_var, metadata, envir = .GlobalEnv)
    } else {
      message(paste("No metadata dataset found for", dataset))
    }
    
    if (removal) {
      min_nonzero_threshold <-prevalence_threshold * nrow(data)  
      low_mean_cols <- colnames(data)[colMeans(data, na.rm = TRUE) < raretaxaremovalthreshold]
      sparse_cols <- colnames(data)[colSums(data > 0, na.rm = TRUE) < min_nonzero_threshold]
      species_to_remove <- union(low_mean_cols, sparse_cols)
      data <- data[, !(colnames(data) %in% species_to_remove)]
    }
    
    print(paste("Assigning processed data to variable:", paste0(tolower(gsub(" ", "_", dataset)), "_relativedata")))
    assign(paste0(tolower(gsub(" ", "_", dataset)), "_relativedata"), data, envir = .GlobalEnv)
  }
}

build_list <- function(datasets) {
  datasets$file_path <- file.path(getwd(), datasets$file_path)
  result_list <- list()
  
  for (dataset in unique(datasets$dataset_name)) {
    dataset_row <- datasets[datasets$dataset_name == dataset, ]
    
    sequencingtype <- unique(dataset_row$sequencingtype)
    loadtype <- unique(dataset_row$loadtype)
    sourcetype <- unique(dataset_row$source_type)
    abundance <- unique(dataset_row$abundance)
    trained <- unique(dataset_row$trainingset)
    
    print(paste("📂 Processing dataset:", dataset))
    
    dataset_key <- tolower(gsub(" ", "_", dataset))
    relativedata_var <- paste0(dataset_key, "_relativedata")
    measured_var <- paste0(dataset_key, "_Measured")
    metadata_var <- paste0(dataset_key, "_metadata")  
    countdata_var <- paste0(dataset_key, "_countdata")
    
    if (sequencingtype == "shotgun metagenomics") {
      measured_var <- paste0(dataset_key, "__Measured")
      galaxypredicted_var <- paste0(dataset_key, "_galaxy_Predicted")
      metacardispredicted_var <- paste0(dataset_key, "_metacardis_Predicted")
    } else {
      predicted_var <- paste0(dataset_key, "_Predicted")
    }
    
    relativedata <- tryCatch(get(relativedata_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    measured <- tryCatch(get(measured_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    metadata <- tryCatch(get(metadata_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL) 
    countdata <- tryCatch(get(countdata_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL) 
    
    if (sequencingtype == "shotgun metagenomics") {
      galaxypredicted <- tryCatch(get(galaxypredicted_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
      metacardispredicted <- tryCatch(get(metacardispredicted_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    } else {
      predicted <- tryCatch(get(predicted_var, envir = .GlobalEnv, inherits = TRUE), error = function(e) NULL)
    }
    
    if (!is.null(relativedata)) {
      if (!is.null(measured)) {
        sampleIDcol <- unique(as.character(datasets$sampleIDcol[datasets$dataset_name == dataset &
                                                                  datasets$source_type == "Measured" &
                                                                  datasets$abundance == "total"]))
        if (length(sampleIDcol) < 1 || !(sampleIDcol %in% colnames(measured))) {
          sampleIDcol <- find_best_match(sampleIDcol, measured)
          if (!is.na(sampleIDcol)) {
            message(paste("✅ Using matched column for Measured:", sampleIDcol))
            datasets$sampleIDcol[datasets$dataset_name == dataset &
                                   datasets$source_type == "Measured" &
                                   datasets$abundance == "total"] <- sampleIDcol
          } else {
            message(paste("⚠️ Sample ID column not found in Measured dataset:", dataset))
          }
        }
        if (!is.na(sampleIDcol)) {
          rownames(measured) <- measured[[sampleIDcol]]
          measured[[sampleIDcol]] <- NULL
        }
      }
      
      if (!is.null(metadata)) {
        sampleIDcol <- unique(as.character(datasets$sampleIDcol[datasets$dataset_name == dataset &
                                                                  datasets$abundance == "metadata"]))
        if (length(sampleIDcol) < 1 || !(sampleIDcol %in% colnames(metadata))) {
          sampleIDcol <- find_best_match(sampleIDcol, metadata)
          if (!is.na(sampleIDcol)) {
            message(paste("✅ Using matched column for Metadata:", sampleIDcol))
            datasets$sampleIDcol[datasets$dataset_name == dataset &
                                   datasets$abundance == "metadata"] <- sampleIDcol
          } else {
            message(paste("⚠️ Sample ID column not found in Metadata dataset:", dataset))
          }
        }
        if (!is.na(sampleIDcol)) {
          rownames(metadata) <- metadata[[sampleIDcol]]
          metadata[[sampleIDcol]] <- NULL
        }
      }
      
      prediction_datasets <- list()
      if (sequencingtype == "shotgun metagenomics") {
        prediction_datasets <- list(metacardispredicted, galaxypredicted)
      } else {
        prediction_datasets <- list(predicted)
      }
      
      sampleIDcol <- unique(as.character(datasets$sampleIDcol[datasets$dataset_name == dataset &
                                                                datasets$source_type == "Predicted" &
                                                                datasets$abundance == "total"]))
      
      for (j in seq_along(prediction_datasets)) {
        pred_data <- prediction_datasets[[j]]
        if (!is.null(pred_data) && (length(sampleIDcol) < 1 || !(sampleIDcol %in% colnames(pred_data)))) {
          sampleIDcol <- find_best_match(sampleIDcol, pred_data)
          if (!is.na(sampleIDcol)) {
            message(paste("✅ Using matched column for Predicted:", sampleIDcol))
            datasets$sampleIDcol[datasets$dataset_name == dataset &
                                   datasets$source_type == "Predicted" &
                                   datasets$abundance == "total"] <- sampleIDcol
          } else {
            message(paste("⚠️ Sample ID column not found in Predicted dataset:", dataset))
          }
        }
        if (!is.na(sampleIDcol)) {
          rownames(pred_data) <- pred_data[[sampleIDcol]]
          pred_data[[sampleIDcol]] <- NULL
          prediction_datasets[[j]] <- pred_data
        }
      }
      
      datasets_list <- c(list(rownames(relativedata)), if (!is.null(countdata)) list(rownames(countdata)) else list(),
                         if (!is.null(measured)) list(rownames(measured)) else list(), 
                         if (!is.null(metadata)) list(rownames(metadata)) else list(),  
                         lapply(prediction_datasets, function(df) if (!is.null(df)) rownames(df) else NULL))
      
      common_rownames <- Reduce(intersect, datasets_list)
      
      relativedata <- relativedata[common_rownames, , drop = FALSE]
      if (!is.null(countdata)) countdata <- countdata[common_rownames, , drop = FALSE] 
      if (!is.null(measured)) measured <- measured[common_rownames, , drop = FALSE]
      if (!is.null(metadata)) metadata <- metadata[common_rownames, , drop = FALSE]  
      prediction_datasets <- lapply(prediction_datasets, function(df) if (!is.null(df)) df[common_rownames, , drop = FALSE])
      
      synchronized_data <- list(
        relabun_df = relativedata,
        counts_df = countdata,
        measured_load_df = measured,
        metadata_df = metadata 
      )
      
      if (sequencingtype == "shotgun metagenomics") {
        synchronized_data$galaxypredicted_load_df <- prediction_datasets[[2]]
        synchronized_data$metacardispredicted_load_df <- prediction_datasets[[1]]
      } else {
        synchronized_data$predicted_load_df <- prediction_datasets[[1]]
      }
      
      result_list[[dataset]] <- synchronized_data
    } else {
      message(paste("⚠️ Skipping dataset:", dataset, "due to missing relativedata."))
    }
  }

  standardize_colnames <- function(df, type = c("measured", "predicted"), target_col) {
    type <- match.arg(type) 
    if (nrow(df) == 0 || is.null(df)) return(df)  
    
    df <- df %>% tibble::rownames_to_column("SampleID")
    
    if (type == "measured") {
      if (length(target_col) == 0) stop("No valid column found for MeasuredLoad")
      colnames(df)[colnames(df) == target_col] <- "MeasuredLoad"
    } else if (type == "predicted") {
      if (!"load" %in% colnames(df)) stop("No 'load' column found in predicted dataframe")
      colnames(df)[colnames(df) == "load"] <- "PredictedLoad"
    }
    
    df <- df %>% dplyr::select(SampleID, everything())
    return(df)
  }

  nomenclature <- names(result_list)
  result_list <- lapply(names(result_list), function(study_name) {
    cat("Processing study:", study_name, "\n")
    study <- result_list[[study_name]]
    
    target_col <- unique(as.character(datasets$loadcol[datasets$dataset_name == study_name &
                                                         datasets$source_type == "Measured" &
                                                         datasets$abundance == "total"]))
    
    if (!is.null(study$measured_load_df)) {
      if (length(target_col) < 1 || !(target_col %in% colnames(study$measured_load_df))) {
        target_col <- find_best_match(target_col, study$measured_load_df)
      }
      study$measured_load_df <- standardize_colnames(study$measured_load_df, type = "measured", target_col)
    }
    if (all(datasets[datasets$dataset_name == study_name, "sequencingtype"] == "shotgun metagenomics")) {
      if (!is.null(study$galaxypredicted_load_df)) {
        study$galaxypredicted_load_df <- standardize_colnames(study$galaxypredicted_load_df, type = "predicted")
      }
      if (!is.null(study$metacardispredicted_load_df)) {
        study$metacardispredicted_load_df <- standardize_colnames(study$metacardispredicted_load_df, type = "predicted")
      }
    } else {
      if (!is.null(study$predicted_load_df)) {
        study$predicted_load_df <- standardize_colnames(study$predicted_load_df, type = "predicted")
      }
    }
    
    return(study)
  })
  
  names(result_list) <- nomenclature
  return(result_list)
}

# ==== Data Visualization ======================================================

prepare_plotting_dataframes <- function(
    cor_results,
    r2_results,
    studies,
    datasets
) {
  
  create_rowlabel <- function(study, reference) {
    study_clean <- gsub("\\s*\\(MLP Training data\\)", "", study)
    dplyr::case_when(
      tolower(reference) == "galaxy" ~ paste0(study_clean, " Model: GALAXY (MLP Training data)"),
      tolower(reference) == "metacardis" ~ paste0(study_clean, " Model: MetaCardis (MLP Training data)"),
      tolower(reference) == "vandeputte2021" ~ paste0(study_clean, " Model: Vandeputte2021 (MLP Training data)"),
      TRUE ~ paste0(study_clean, " Model: ", reference)
    )
  }

  parse_rowlabel <- function(df) {
    df %>%
      dplyr::mutate(
        RowLabel_Study = sub("(.*?) Model:.*", "\\1", RowLabel),
        RowLabel_Model = sub(".*? Model: (.*)", "\\1", RowLabel)
      )
  }

  correlation_df <- cor_results %>%
    dplyr::select(Study, Feature, Pearson_r_residual, Reference) %>%
    tidyr::pivot_longer(
      cols      = c("Pearson_r_residual"),
      names_to  = "LoadType",
      values_to = "Correlation"
    ) %>%
    dplyr::mutate(
      LoadType = sub("Pearson_r_", "", LoadType)
    ) %>%
    dplyr::filter(!is.na(Correlation)) %>%
    dplyr::mutate(
      RowLabel = create_rowlabel(Study, Reference)
    ) %>%
    dplyr::mutate(
      RowLabel = factor(RowLabel, levels = unique(RowLabel[order(Reference)]))
    )

  non_correlation_df <- cor_results %>%
    dplyr::select(Study, Feature, FeatureMean, Reference) %>%
    tidyr::pivot_longer(
      cols      = "FeatureMean",
      names_to  = "LoadType",
      values_to = "Proportion"
    ) %>%
    dplyr::mutate(
      LoadType = sub("N_", "", LoadType),
      RowLabel = create_rowlabel(Study, Reference)
    ) %>%
    dplyr::filter(!is.na(Proportion)) %>%  
    dplyr::distinct()

  reference_values <- non_correlation_df %>%
    dplyr::filter(Study == Reference) %>%
    dplyr::select(Feature, Reference, referenceProportion = Proportion)
  
  non_correlation_df <- non_correlation_df %>%
    dplyr::left_join(reference_values, by = c("Feature", "Reference")) %>%
    dplyr::group_by(Feature) %>%
    dplyr::mutate(
      Distance        = abs(Proportion - referenceProportion),
      Scaled_Distance = (Distance - min(Distance, na.rm = TRUE)) /
        (max(Distance, na.rm = TRUE) - min(Distance, na.rm = TRUE))
    ) %>%
    dplyr::ungroup()
  
  r2_df <- purrr::map_dfr(
    names(r2_results),
    function(study) {
      r2_entry <- r2_results[[study]]
      if (is.list(r2_entry)) {

        tibble::tibble(
          Study  = study,
          Metric = names(r2_entry),
          Value  = as.numeric(r2_entry),
          RowLabel = paste0(study, " Model: ", names(r2_entry))
        )
      } else {

        tibble::tibble(
          Study   = study,
          Metric  = names(r2_entry),
          Value   = r2_entry,
          RowLabel = paste0(study, " Model: Vandeputte2021 (MLP Training data)") 
        )
      }
    }
  ) %>%
    dplyr::mutate(
      Metric = dplyr::case_when(
        Metric == "r2_galaxy"         ~ "GALAXY R²",
        Metric == "r_galaxy"          ~ "GALAXY R",
        Metric == "p_galaxy"          ~ "GALAXY p-value",
        Metric == "r2_metacardis"     ~ "MetaCardis R²",
        Metric == "r_metacardis"      ~ "MetaCardis R",
        Metric == "p_metacardis"      ~ "MetaCardis p-value",
        Metric == "r2_vandeputte2021" ~ "Vandeputte2021 R²",
        Metric == "r_vandeputte2021"  ~ "Vandeputte2021 R",
        Metric == "p_vandeputte2021"  ~ "Vandeputte2021 p-value",
        TRUE ~ Metric
      ),

      RowLabel = dplyr::case_when(
        grepl("galaxy", tolower(Metric)) ~ paste0(
          gsub("\\s*\\(MLP Training data\\)", "", Study),
          " Model: GALAXY (MLP Training data)"
        ),
        grepl("metacardis", tolower(Metric)) ~ paste0(
          gsub("\\s*\\(MLP Training data\\)", "", Study),
          " Model: MetaCardis (MLP Training data)"
        ),
        TRUE ~ paste0(
          gsub("\\s*\\(MLP Training data\\)", "", Study),
          " Model: Vandeputte2021 (MLP Training data)"
        )
      )
    )

  log_stats_df <- purrr::map_dfr(
    .x = names(studies),
    .f = extract_log_stats,
    study_list = studies,
    datasets = datasets
  ) %>%
    dplyr::left_join(datasets, by = c("Study" = "dataset_name"))

  selected_df <- log_stats_df %>%
    dplyr::rename(
      `log(load) mean (measured)`  = log_measured_mean,
      `log(load) sd (measured)`    = log_measured_sd,
      `N (measured)`               = N_measured,
      
      `log(load) mean (predicted)` = log_predicted_mean,
      `log(load) sd (predicted)`   = log_predicted_sd,
      `N (predicted)`              = N_predicted,
      
      `log(load) mean (predicted_galaxy)`    = log_galaxy_mean,
      `log(load) sd (predicted_galaxy)`      = log_galaxy_sd,
      `N (predicted_galaxy)`                 = N_galaxy,
      
      `log(load) mean (predicted_metacardis)` = log_metacardis_mean,
      `log(load) sd (predicted_metacardis)`   = log_metacardis_sd,
      `N (predicted_metacardis)`              = N_metacardis
    ) %>%
    tidyr::pivot_longer(
      cols = matches("^log\\(load\\) (mean|sd) \\("),
      names_to = c("Stat", "Type"),
      names_pattern = "^log\\(load\\) (mean|sd) \\(([^)]+)\\)$",
      values_to = "Value"
    ) %>%
    tidyr::pivot_wider(
      names_from  = Stat,
      values_from = Value,
      names_prefix = "log(load)_"
    ) %>%
    dplyr::mutate(
      RowLabel = dplyr::case_when(
        is.na(trainingset) | trainingset == "" ~ gsub("\\s*\\(MLP Training data\\)", "", Study),
        TRUE ~ paste0(
          gsub("\\s*\\(MLP Training data\\)", "", Study),
          " Model: ",
          trainingset,
          " (MLP Training data)"
        )
      ),
      `log(load) Mean`  = `log(load)_mean`,
      `log(load) SD`    = `log(load)_sd`
    ) %>%
    dplyr::select(-`log(load)_mean`, -`log(load)_sd`) %>%
    dplyr::rename(`N` = `N (measured)`) %>%
    dplyr::select(-`N (predicted)`, -`N (predicted_galaxy)`, -`N (predicted_metacardis)`) %>%
    stats::na.omit()

  selected_df <- selected_df %>%
    dplyr::mutate(
      trainingset = dplyr::if_else(is.na(trainingset), NA_character_, trainingset)
    )
  
  demo_df <- selected_df %>%
    tidyr::pivot_longer(
      cols = c(
        `log(load) Mean`, `log(load) SD`, loadtype, trainingset,
        organismtype, sampletype, N, sequencingtype
      ),
      names_to = "Feature",
      values_to = "Value",
      values_transform = list(Value = as.character)
    ) %>%
    dplyr::select(
      Study, RowLabel, Type, Feature, Value,
      file_path, loadcol, sampleIDcol, abundance, source_type
    ) %>%
    dplyr::distinct(Study, RowLabel, Type, Feature, .keep_all = TRUE) 

  numeric_df <- demo_df %>%
    dplyr::filter(!is.na(suppressWarnings(as.numeric(Value))))

  categorical_df <- demo_df %>%
    dplyr::filter(is.na(suppressWarnings(as.numeric(Value)))) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel))
  
  numeric_df_clean <- numeric_df %>%
    dplyr::mutate(Value = as.numeric(Value)) %>%
    tidyr::pivot_wider(
      id_cols     = c(Study, Feature, RowLabel),
      names_from  = Type,
      values_from = Value,
      values_fn   = list(Value = mean)
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Distance = dplyr::if_else(
        grepl("measured", RowLabel, ignore.case = TRUE),
        0,
        abs(measured - predicted)
      ),
      Distance = dplyr::case_when(
        grepl("Model: MetaCardis", RowLabel, ignore.case=TRUE) ~ abs(measured - predicted_metacardis),
        grepl("Model: GALAXY", RowLabel, ignore.case=TRUE)     ~ abs(measured - predicted_galaxy),
        TRUE ~ Distance
      ),
      Distance = tidyr::replace_na(Distance, 0),
      Annotation = dplyr::case_when(
        grepl("measured", RowLabel, ignore.case=TRUE)           ~ round(measured, 2),
        grepl("Model: MetaCardis", RowLabel, ignore.case=TRUE)  ~ round(predicted_metacardis, 2),
        grepl("Model: GALAXY", RowLabel, ignore.case=TRUE)      ~ round(predicted_galaxy, 2),
        TRUE                                                    ~ round(predicted, 2)
      )
    ) %>%
    dplyr::group_by(Feature) %>%
    dplyr::mutate(
      Scaled_Distance = (Distance - min(Distance, na.rm = TRUE)) /
        (max(Distance, na.rm = TRUE) - min(Distance, na.rm = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel))

  final_levels <- categorical_df %>%
    dplyr::filter(Feature == "sequencingtype") %>%
    dplyr::distinct(RowLabel,Value) %>%
    dplyr::arrange(factor(Value, levels = c("16S rRNA", "shotgun metagenomics"))) %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>% 
    dplyr::mutate(RowLabel = as.character(RowLabel)) %>%
    dplyr::pull(RowLabel) 
  
  correlation_df <- correlation_df %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel)) %>%
    dplyr::arrange(RowLabel)
  
  
  non_correlation_df <- non_correlation_df %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel)) %>%
    dplyr::arrange(RowLabel) 
  
  r2_df <- r2_df %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel)) %>%
    dplyr::arrange(RowLabel) 
  
  numeric_df_clean <- numeric_df_clean %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::arrange(RowLabel) 
  
  demo_df <- demo_df %>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel) & 
                    source_type %in% c('Predicted')) %>%
    dplyr::arrange(RowLabel) 
  
  categorical_df <- categorical_df%>%
    dplyr::mutate(
      RowLabel = gsub("\\s*\\(MLP Training data\\)", "", RowLabel) 
    ) %>%
    dplyr::mutate(RowLabel = factor(RowLabel, levels = final_levels, ordered = TRUE)) %>%
    dplyr::filter(!is.na(RowLabel) & grepl("Model:", RowLabel) & 
                    Type %in% c('predicted', 'predicted_galaxy', 'predicted_metacardis')) %>%
    dplyr::arrange(RowLabel) 
  
  correlation_df     <- parse_rowlabel(correlation_df)
  non_correlation_df <- parse_rowlabel(non_correlation_df)
  r2_df              <- parse_rowlabel(r2_df)
  selected_df        <- parse_rowlabel(selected_df)
  demo_df            <- parse_rowlabel(demo_df)
  numeric_df_clean   <- parse_rowlabel(numeric_df_clean)
  categorical_df     <- parse_rowlabel(categorical_df)
  
  return(list(
    correlation_df     = correlation_df,
    non_correlation_df = non_correlation_df,
    r2_df              = r2_df,
    selected_df        = selected_df,
    demo_df            = demo_df,
    numeric_df_clean   = numeric_df_clean,
    categorical_df     = categorical_df
  ))
}

plot_heatmap <- function(df, type, plotytitle=NULL) {
  if (type == "correlation") {
    p <- ggplot(df, aes(x = Feature, y = RowLabel, fill = Correlation)) +
      geom_tile(color = "grey90", linewidth = 1) +
      geom_text(aes(label = round(Correlation, 2)), size = 7, fontface = "bold") +
      scale_fill_gradient2(
        low = scales::muted("blue"), high = scales::muted("red"), mid = "white",
        midpoint = 0, limits = c(-1, 1), na.value = "grey80",
        name = "Pearson R"
      )+
      labs(title = "Sample Feature Correlation with Residual (Pearson R)")
    
  } else if (type == "non_correlation") {
    p <- ggplot(df, aes(x = Feature, y = RowLabel, fill = Scaled_Distance)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = round(Proportion, 2)), size = 7, fontface = "bold") +
      scale_fill_gradient(
        low = "darkgreen", high = "white",
        limits = c(0, 1),
        name = "Dataset Similarity Distance"
      )  +
      labs(title = "Dataset Taxa Similarity with Training Data")
    
  } else if (type == "r2") {
    df <- df %>%
      mutate(
        Metric_Short = case_when(
          grepl("r2", Metric) ~ "R²",
          grepl("r$", Metric) ~ "R",
          grepl("R²", Metric) ~ "R²",
          grepl("R$", Metric) ~ "R",
          grepl("p", Metric) ~ "p-value",
          grepl("p-value", Metric) ~ "p-value",
          TRUE ~ Metric 
        ),
        ScaleType = case_when(
          Metric_Short == "p-value" ~ "p-value",
          Metric_Short == "R²" ~ "R²",
          Metric_Short == "R" ~ "R",
          TRUE ~ Metric_Short
        ),
        FillColor = case_when(
          ScaleType == "R²" ~ scales::rescale(Value, to = c(0, 1)),
          ScaleType == "R"  ~ scales::rescale(Value, to = c(-1, 1)),
          TRUE ~ NA_real_
        )
      )
    
    p <- ggplot(df, aes(x = Metric_Short, y = RowLabel)) +
      geom_tile(data = subset(df, ScaleType == "R²"), aes(fill = FillColor), color = "black", linewidth = 1) +
      geom_tile(data = subset(df, ScaleType == "R"), aes(fill = FillColor), color = "black", linewidth = 1) +
      geom_tile(data = subset(df, ScaleType == "p-value"), fill = "white", color = "black", linewidth = 1) +
      geom_text(aes(label = round(Value, 5)), size = 7, fontface = "bold") +
      scale_fill_gradientn(
        colors = c("red", "white", "darkgreen"),
        limits = c(-1, 1),
        name = "Value",
        na.value = "grey"
      ) +
      labs(title = "Prediction Performance")
    
  } else if (type == "numeric") {
    subset_feature1 <- df %>% filter(Feature == "log(load) Mean")
    subset_feature2 <- df %>% filter(Feature == "log(load) SD")
    subset_feature3 <- df %>% filter(Feature == "N")
    
    p <- ggplot() +
      geom_tile(data = subset_feature1, aes(x = Feature, y = RowLabel, fill = Distance), color = "grey90", linewidth = 1) +
      scale_fill_gradient(low = "white", high = "red", name = "Distance from Measured (Mean)") +
      ggnewscale::new_scale_fill() +
      geom_tile(data = subset_feature2, aes(x = Feature, y = RowLabel, fill = Distance), color = "grey90", linewidth = 1) +
      scale_fill_gradient(low = "white", high = "red", name = "Distance from Measured (SD)") +
      ggnewscale::new_scale_fill() +
      geom_tile(data = subset_feature3, aes(x = Feature, y = RowLabel), fill = "white", color = "grey90", linewidth = 1) +
      geom_text(data = subset_feature1, aes(x = Feature, y = RowLabel, label = Annotation), size = 7, fontface = "bold") +
      geom_text(data = subset_feature2, aes(x = Feature, y = RowLabel, label = Annotation), size = 7, fontface = "bold") +
      geom_text(data = subset_feature3, aes(x = Feature, y = RowLabel, label = Annotation), size = 7, fontface = "bold") +
      labs(title = "Predicted Summary")
    
  } else if (type == "categorical") {
    p <- ggplot(df, aes(x = Feature, y = RowLabel)) +
      geom_tile(fill = "white", color = "black", linewidth = 1) +
      geom_text(aes(label = Value), size = 7, fontface = "bold") +  # Unified text size+
      labs(title = "Study Characteristics")
  }
  
  # Apply consistent theming to all plots
  p <- p +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 20, face = "bold"),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 18),
      plot.title = element_text(size = 30, hjust = 0.5, face = "bold")
    )
  
  if (!is.null(plotytitle)) {
    p <- p + theme(axis.text.y = element_markdown(size = 22, face = "bold", hjust=1))
  }
  
  return(p)
}

plot_studyfeatures <- function(data, feature, filename = NULL, scalemeasurement=NULL) {
  data |>
    tidyplot(x = Study, y = {{ feature }}, color = {{ scalemeasurement }}, width=65) |>  
    add_violin(draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.5, trim = TRUE) |>
    adjust_x_axis(rotate_labels = 90)|>
    # add_test_asterisks(method = "wilcox_test", p.adjust.method = "BH", ref.group = 1,
    #                    hide.ns = TRUE,
    #                    padding_top = .08,
    #                     label.size = 1.5,
    #                     step.increase = 0.2,
    #                     p.adjust.by = "panel",
    #                     symnum.args = list(
    #                     cutpoints = c(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, Inf),
    #                     symbols = c("****","***", "**", "*", "", "ns")
    #                   ))|>
    adjust_size(width = 25, height = 15, unit = "cm")|>
    adjust_y_axis(title = str_wrap(unique(data$Feature), width = 20))|>
    adjust_legend_title("Scale Measurement")|>
    split_plot(by = sequencingtype)|>
    save_plot(filename, bg="white", dpi=600, width = 10,
              height = 5,
              units = c("in"))
}


plot_studyfeatures_grouped <- function(data, feature, filename = NULL, scalemeasurement=NULL) {
  data |>
    tidyplot(x = Study, y = {{ feature }}, color = {{ scalemeasurement }}, width=65) |>  
    add_violin(draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.5, trim = TRUE) |>
    adjust_x_axis(rotate_labels = 90)|>
    add_test_asterisks(method = "wilcox_test", p.adjust.method = "BH", 
                       hide.ns = TRUE,
                       padding_top = .08,
                        label.size = 1.5,
                        step.increase = 0.2,
                        p.adjust.by = "panel",
                        symnum.args = list(
                        cutpoints = c(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, Inf),
                        symbols = c("****","***", "**", "*", "", "ns")
                      ))|>
    adjust_size(width = 25, height = 15, unit = "cm")|>
    adjust_y_axis(title = "Percent (%)")|>
    adjust_legend_title("Feature")|>
    split_plot(by = Reference)|>
    save_plot(filename, , bg="white", dpi=600, width = 10,
              height = 5,
              units = c("in"))
}

loadcorrelation <- function(study_name, study_data, output_dir = "plots") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  study_name_clean <- gsub("\\s*\\(MLP Training data\\)\\s*", " ", study_name)
  study_name_clean <- trimws(study_name_clean)
  
  measured_df <- study_data$measured_load_df
  p <- NULL
  
  tidyplot_theme <- function(base_size = 16) {  
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                                  margin = margin(b = 5)),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10), hjust = 0.5),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10), hjust = 0.5),
        axis.title = element_text(size = 14, face = "plain"),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),  
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  }
  
  tidyplot_colors <- c("#000000", "#0C58CA", "#FF8C27")
  
  if ("predicted_load_df" %in% names(study_data)) {
    predicted_df <- study_data$predicted_load_df
    df <- merge(measured_df, predicted_df, by = "SampleID")
    
    N <- nrow(df)
    
    model <- lm(PredictedLoad ~ MeasuredLoad, data = df)
    r2 <- round(summary(model)$r.squared, 2)
    pval <- summary(model)$coefficients[2, 4]
    x_annot <- max(df$MeasuredLoad, na.rm = TRUE) - 0.05 * diff(range(df$MeasuredLoad, na.rm = TRUE))
    y_annot <- min(df$PredictedLoad, na.rm = TRUE) + 0.05 * diff(range(df$PredictedLoad, na.rm = TRUE))
    
    p <- ggplot(df, aes(MeasuredLoad, PredictedLoad)) +
      geom_point(size = 1, alpha = 0.4, color = "black", fill = "white", shape = 21) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = tidyplot_colors[1],
                  fill = tidyplot_colors[1]) +
      annotate("text",
               x = x_annot,
               y = y_annot,
               label = sprintf("italic(R)^2 == %.2f", r2),
               size = 7/.pt,
               hjust = 1, 
               vjust = 0,
               parse = TRUE,
               color = tidyplot_colors[1]) +
      labs(title = paste0(study_name_clean, " (N = ", N, ")"),
           x = "Measured",
           y = "Predicted (Nishijima 2024)") +
      tidyplot_theme() +
      coord_cartesian(clip = "off") +
      scale_x_continuous() + 
      scale_y_continuous()
    
  } else if ("galaxypredicted_load_df" %in% names(study_data) &&
             "metacardispredicted_load_df" %in% names(study_data)) {
    galaxy_df <- study_data$galaxypredicted_load_df
    metacardis_df <- study_data$metacardispredicted_load_df
    galaxy_df$Source <- "GALAXY"
    metacardis_df$Source <- "MetaCardis"
    N <- nrow(galaxy_df)
    galaxy_merge <- merge(measured_df, galaxy_df, by = "SampleID")
    metacardis_merge <- merge(measured_df, metacardis_df, by = "SampleID")
    df <- rbind(galaxy_merge, metacardis_merge)
    
    sources <- unique(df$Source)
    colors <- tidyplot_colors[2:3]
    
    stats_df <- df |>
      group_by(Source) |>
      summarise(
        model = list(lm(PredictedLoad ~ MeasuredLoad)),
        .groups = "drop"
      ) |>
      mutate(
        r2 = sapply(model, \(m) round(summary(m)$r.squared, 2)),
        pval = sapply(model, \(m) summary(m)$coefficients[2, 4]),
        x = max(df$MeasuredLoad, na.rm = TRUE) - 0.05 * diff(range(df$MeasuredLoad, na.rm = TRUE)),
        y = min(df$PredictedLoad, na.rm = TRUE) + 0.05 * diff(range(df$PredictedLoad, na.rm = TRUE)) +
          (rev(seq_len(n())) - 1) * 0.05 * diff(range(df$PredictedLoad, na.rm = TRUE)),
        color = colors[match(Source, sources)]
      )
    
    p <- ggplot(df, aes(MeasuredLoad, PredictedLoad, color = Source)) +
      geom_point(size = 1, alpha = 0.4, fill = "white", shape = 21) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1, aes(fill = Source)) +
      scale_color_manual(values = colors) +
      labs(title = paste0(study_name_clean, " (N = ", N, ")"),
           x = "Measured",
           y = "Predicted (Nishijima 2024)") +
      tidyplot_theme() +
      theme(legend.position = "none") +
      coord_cartesian(clip = "off") +
      scale_x_continuous() +  
      scale_y_continuous() +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) 
    
    for (i in seq_len(nrow(stats_df))) {
      p <- p + 
        annotate("text",
                 x = stats_df$x[i],
                 y = stats_df$y[i],
                 label = sprintf("%s:~italic(R)^2 == %.2f", 
                                 stats_df$Source[i], 
                                 stats_df$r2[i]),
                 color = stats_df$color[i],
                 size = 7/.pt,
                 hjust = 1,
                 vjust = 0,
                 parse = TRUE)
    }
  } else {
    warning(paste("Skipping study:", study_name, "— unexpected structure"))
    return(NULL)
  }
  
  ggsave(
    filename = file.path(output_dir, paste0(study_name_clean, "_load_plot.pdf")),
    plot = p,
    dpi = 600,
    width = 5,
    height = 5,
    units = "in",
    bg = "white",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(output_dir, paste0(study_name_clean, "_load_plot.png")),
    plot = p,
    dpi = 600,
    width = 5,
    height = 5,
    units = "in",
    bg = "white"
  )
}


plot_residuals <- function(data_list) {
  rsq <- function(x, y) cor(x, y)^2
  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
  mae <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
  mape <- function(y, yhat) mean(abs((y - yhat) / y), na.rm = TRUE) * 100
  rename_model <- function(x) {
    if (x == "predicted_load_df") {
      "Vandeputte MLP"
    } else if (grepl("^galaxy", x, ignore.case = TRUE)) {
      "GALAXY MLP"
    } else if (grepl("^metacardis", x, ignore.case = TRUE)) {
      "MetaCardis MLP"
    } else {
      x
    }
  }
  plot_list <- list()
  for (study_name in names(data_list)) {
    if (!"measured_load_df" %in% names(data_list[[study_name]])) {
      message("Study '", study_name, "' has no measured_load_df. Skipping.")
      next
    }
    measured_df <- data_list[[study_name]]$measured_load_df
    pred_df_names <- names(data_list[[study_name]])[grepl("predicted_load_df$", 
                                                          names(data_list[[study_name]]))]
    if (length(pred_df_names) == 0) {
      message("Study '", study_name, "' has no predicted_load_df-type data frames. Skipping.")
      next
    }
    residuals_list <- list()
    metrics_list <- list()
    for (pred_name in pred_df_names) {
      pred_df <- data_list[[study_name]][[pred_name]]
      if (nrow(measured_df) != nrow(pred_df)) {
        warning("Study '", study_name, "' has mismatched row counts for '", pred_name, "'.")
      }
      new_model <- rename_model(pred_name)
      r2_value    <- rsq(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      rmse_value  <- rmse(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      mae_value   <- mae(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      mape_value  <- mape(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      metrics_list[[new_model]] <- list(r2 = r2_value, rmse = rmse_value, mae = mae_value, mape = mape_value)
      tmp <- data.frame(
        RowID = seq_len(nrow(measured_df)),
        MeasuredLoad = measured_df$MeasuredLoad,
        PredictedLoad = pred_df$PredictedLoad,
        Residual = measured_df$MeasuredLoad - pred_df$PredictedLoad,
        Model = new_model
      )
      residuals_list[[pred_name]] <- tmp
    }
    residuals_df <- bind_rows(residuals_list)
    p <- ggplot(residuals_df, aes(x = PredictedLoad, y = Residual)) +
      geom_point(aes(color = Model)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~ Model, scales = "free") +
      labs(
        title = paste("Residual Plot for Study:", study_name),
        x = "Predicted Load",
        y = "Residual (Measured - Predicted)"
      ) +
      tidyplot_theme() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none"
      ) +
      scale_color_manual(
        values = c("Vandeputte MLP" = "#000000",
                   "GALAXY MLP" = "#0C58CA",
                   "MetaCardis MLP" = "#FF8C27")
      )
    metrics_str <- paste(
      sapply(names(metrics_list), function(model_name) {
        m <- metrics_list[[model_name]]
        sprintf("\n%s: R²=%.3f, RMSE=%.3f, MAE=%.3f, MAPE=%.1f%%", 
                model_name, m$r2, m$rmse, m$mae, m$mape)
      }),
      collapse = "; "
    )
    p <- p + labs(subtitle = metrics_str)
    ggsave(paste0("plots/residualplot_",study_name,".png"), plot = p, width = 8, height = 5, units = "in", dpi = 600, bg = "white")
    plot_list[[study_name]] <- p
  }
  return(plot_list)
}

plot_measured_vs_predicted <- function(data_list, rank_threshold = 25) {
  rsq <- function(x, y) cor(x, y)^2
  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
  mae <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
  mape <- function(y, yhat) mean(abs((y - yhat) / y), na.rm = TRUE) * 100
  rename_model <- function(x) {
    if (x == "predicted_load_df") {
      "Vandeputte MLP"
    } else if (grepl("^galaxy", x, ignore.case = TRUE)) {
      "GALAXY MLP"
    } else if (grepl("^metacardis", x, ignore.case = TRUE)) {
      "MetaCardis MLP"
    } else {
      x
    }
  }
  plot_list <- list()
  for (study_name in names(data_list)) {
    if (!"measured_load_df" %in% names(data_list[[study_name]])) {
      message("Study '", study_name, "' has no measured_load_df. Skipping.")
      next
    }
    measured_df <- data_list[[study_name]]$measured_load_df
    pred_df_names <- names(data_list[[study_name]])[grepl("predicted_load_df$", 
                                                          names(data_list[[study_name]]))]
    if (length(pred_df_names) == 0) {
      message("Study '", study_name, "' has no predicted_load_df-type data frames. Skipping.")
      next
    }
    combined_list <- list()
    r2_list <- list()
    metrics_list <- list()
    for (pred_name in pred_df_names) {
      pred_df <- data_list[[study_name]][[pred_name]]
      if (nrow(measured_df) != nrow(pred_df)) {
        warning("Study '", study_name, "' has mismatched row counts for '", pred_name, "'.")
      }
      r2_value <- rsq(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      rmse_value <- rmse(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      mae_value  <- mae(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      mape_value <- mape(measured_df$MeasuredLoad, pred_df$PredictedLoad)
      new_model <- rename_model(pred_name)
      r2_list[[new_model]] <- r2_value
      metrics_list[[new_model]] <- list(r2 = r2_value, rmse = rmse_value, mae = mae_value, mape = mape_value)
      tmp <- data.frame(
        RowID = seq_len(nrow(measured_df)),
        MeasuredLoad = measured_df$MeasuredLoad,
        PredictedLoad = pred_df$PredictedLoad,
        Model = new_model
      )
      tmp$RankMeasured <- rank(tmp$MeasuredLoad, ties.method = "average")
      tmp$RankPredicted <- rank(tmp$PredictedLoad, ties.method = "average")
      tmp$RankDiff <- abs(tmp$RankPredicted - tmp$RankMeasured)
      tmp$SignificantChange <- tmp$RankDiff > rank_threshold
      tmp_long <- pivot_longer(tmp, cols = c(MeasuredLoad, PredictedLoad),
                               names_to = "LoadType", values_to = "LoadValue")
      tmp_long$LoadType <- factor(tmp_long$LoadType, levels = c("MeasuredLoad", "PredictedLoad"))
      combined_list[[pred_name]] <- tmp_long
    }
    combined_df <- bind_rows(combined_list)
    lines_data <- combined_df %>%
      group_by(RowID, Model) %>%
      filter(n() == 2) %>%
      summarise(
        MeasuredLoad = first(LoadValue[LoadType == "MeasuredLoad"]),
        PredictedLoad = first(LoadValue[LoadType == "PredictedLoad"]),
        SignificantChange = first(SignificantChange)
      ) %>% ungroup() %>%
      pivot_longer(cols = c(MeasuredLoad, PredictedLoad),
                   names_to = "LoadType", values_to = "LoadValue") %>%
      mutate(LoadType = factor(LoadType, levels = c("MeasuredLoad", "PredictedLoad")))
    metrics_str <- paste(
      sapply(names(metrics_list), function(model_name) {
        m <- metrics_list[[model_name]]
        sprintf("\n%s: R²=%.3f, RMSE=%.3f, MAE=%.3f, MAPE=%.1f%%", model_name, m$r2, m$rmse, m$mae, m$mape)
      }),
      collapse = "; "
    )
    title_text <- paste("Study:", study_name, "\n", metrics_str)
    p <- ggplot() +
      geom_line(data = lines_data %>% filter(SignificantChange == FALSE),
                aes(x = LoadType, y = LoadValue, group = interaction(RowID, Model)),
                color = "gray50", alpha = 0.1) +
      geom_line(data = lines_data %>% filter(SignificantChange == TRUE),
                aes(x = LoadType, y = LoadValue, group = interaction(RowID, Model)),
                color = "red", alpha = 0.5) +
      geom_point(data = combined_df %>% filter(LoadType == "MeasuredLoad"),
                 aes(x = LoadType, y = LoadValue),
                 color = "black", size = 2) +
      geom_point(data = combined_df %>% filter(LoadType == "PredictedLoad"),
                 aes(x = LoadType, y = LoadValue, color = Model),
                 size = 2) +
      labs(
        title = title_text,
        x = NULL,
        y = "Load"
      ) +
      tidyplot_theme() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "top"
      ) +
      scale_color_manual(
        values = c("Measured" = "#929292",
                   "Vandeputte MLP" = "#000000",
                   "GALAXY MLP" = "#0C58CA",
                   "MetaCardis MLP" = "#FF8C27"),
        name = "Model"
      )
    ggsave(paste0("plots/rankplot_",study_name,".png"), plot = p, width = 8, height = 5, units = "in", dpi = 600, bg = "white")
    plot_list[[study_name]] <- p
  }
  return(plot_list)
}



modelingscale_table_plot <- function(nested_df_list, presentation = FALSE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  
  is_presentation <- is.list(presentation)
  factor_filter_list <- if (is_presentation) presentation else NULL
  
  flat_list <- list()
  for (study in names(nested_df_list)) {
    model_dfs <- list()
    for (model_nm in names(nested_df_list[[study]])) {
      df <- as.data.frame(nested_df_list[[study]][[model_nm]])
      if (nrow(df) > 0) {
        df$Study <- study
        df$Model <- model_nm
        model_dfs[[model_nm]] <- df
      }
    }
    if (length(model_dfs) > 0) {
      flat_list[[study]] <- do.call(rbind, model_dfs)
    }
  }
  
  if (length(flat_list) == 0) {
    warning("No data to process.")
    return(list(export_table = NULL, plot = NULL))
  }
  
  big_df <- do.call(rbind, flat_list)
  needed_cols <- c("Study", "Model", "variable", "estimate.mean", "p.val")
  for (colname in needed_cols) {
    if (!colname %in% colnames(big_df)) {
      stop("Missing required column: ", colname)
    }
  }
  
  final_df <- big_df %>%
    dplyr::select(all_of(needed_cols)) %>%
    rename(Beta = estimate.mean, pval = p.val) %>%
    mutate(
      Model    = factor(Model),
      variable = factor(variable)
    )
  
  final_df <- final_df %>%
    group_by(Study) %>%
    filter(any(Model == "measured")) %>%
    ungroup()
  
  final_df <- final_df %>%
    filter(!grepl("intercept", variable, ignore.case = TRUE))
  
  if (is_presentation) {
    final_df <- final_df %>%
      rowwise() %>%
      filter({
        study_list <- factor_filter_list[[Study]]
        if (is.null(study_list)) FALSE else {
          model_vars <- study_list[[as.character(Model)]]
          if (is.null(model_vars)) FALSE else variable %in% model_vars
        }
      }) %>%
      ungroup()
  }
  
  if (is_presentation) {
    pred_model_priority <- c("vandeputte", "GALAXY", "MetaCardis")
    
    predicted_map <- final_df %>%
      filter(Model %in% pred_model_priority) %>%
      group_by(Study) %>%
      summarize(PredictedModel = pred_model_priority[pred_model_priority %in% Model][1], .groups = "drop")
    
    final_df <- final_df %>%
      left_join(predicted_map, by = "Study") %>%
      filter(
        Model == "measured" |
          Model == "Zero Scale" |
          (!is.na(PredictedModel) & Model == PredictedModel)
      ) %>%
      mutate(
        Model = case_when(
          Model == "measured" ~ "Measured",
          Model == "Zero Scale" ~ "ZSM",
          Model == PredictedModel ~ "Predicted (MLP)",
          TRUE ~ Model
        )
      )
  }
  
  get_var_base <- function(varname) {
    varname <- as.character(varname)
    varname <- sub("^X_mat", "", varname)
    varname <- sub("[0-9a-zA-Z]*$", "", varname)
    return(varname)
  }
  
  first_var_df <- final_df %>%
    group_by(Study, Model) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(var_base = get_var_base(variable)) %>%
    dplyr::select(Study, Model, var_base)
  
  final_df <- final_df %>%
    mutate(var_base = get_var_base(variable)) %>%
    inner_join(first_var_df, by = c("Study", "Model", "var_base"))

  final_df <- final_df %>%
    group_by(Study, Model, var_base) %>%
    mutate(factor_index = row_number()) %>%
    ungroup()

  model_levels <- if (is_presentation) {
    c("Measured", "ZSM", "Predicted (MLP)")
  } else {
    c("measured", "Zero Scale", "CLR", "vandeputte", "GALAXY", "MetaCardis")
  }
  
  final_df <- final_df %>%
    mutate(Model = factor(Model, levels = model_levels))

  final_df <- final_df %>%
    arrange(Model, factor_index) %>%
    mutate(
      CoefName = if (is_presentation) {
        "B1"
      } else {
        paste0("B1_", factor_index)
      },
      CoefName = factor(CoefName, levels = unique(CoefName))
    )

  study_order_df <- final_df %>%
    group_by(Study) %>%
    summarize(
      has_vandeputte = any(Model %in% c("vandeputte", "Predicted (MLP)")),
      has_galaxy     = any(Model %in% c("GALAXY", "Predicted (MLP)"))
    ) %>%
    mutate(
      study_group = case_when(
        has_vandeputte ~ "vandeputte",
        has_galaxy     ~ "GALAXY",
        TRUE           ~ "other"
      )
    ) %>%
    arrange(factor(study_group, levels = c("vandeputte", "GALAXY", "other"))) %>%
    mutate(Study = factor(Study, levels = unique(Study)))
  
  final_df <- final_df %>%
    left_join(study_order_df %>% dplyr::select(Study), by = "Study") %>%
    mutate(Study = factor(Study, levels = levels(study_order_df$Study)))

  final_df <- final_df %>%
    mutate(
      Label = paste0("B: ", sprintf("%.3g", Beta), "\n", "p: ", sprintf("%.3g", pval)),
      Beta_signif = ifelse(pval < 0.05, Beta, NA)
    )

  p <- ggplot(final_df, aes(x = CoefName, y = Study, fill = Beta_signif)) +
    geom_tile(color = "gray85", linewidth = 0.4) +  
    geom_text(aes(label = Label), size = 4, color = "black", lineheight = 0.9) +
    scale_fill_gradient2(
      low = "#D73027",
      mid = "white",
      high = "#4575B4",
      midpoint = 0,
      limits = c(min(final_df$Beta, na.rm = TRUE), max(final_df$Beta, na.rm = TRUE)),
      na.value = "white"
    ) +
    labs(
      x = "Condition/Treatment Coefficient",  
      y = NULL,
      title = NULL
    ) +
    facet_grid(~ Model, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_blank(),       
      axis.ticks.x = element_blank(),         
      axis.text.y = element_text(size = 12, hjust = 1, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      strip.text = element_text(size = 13, face = "bold"),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),  
      legend.position = "none",
      plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
    )
  theme_cell <- function(base_size = 11, base_family = "Arial") {
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.text.y = element_text(
          size = rel(0.9),
          color = "black",
          face = "plain",
          margin = margin(r = 2)
        ),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length = unit(0.75, "mm"),
        strip.text = element_text(
          size = rel(0.85),
          face = "bold",
          margin = margin(b = 2, t = 2)
        ),
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.85)), 
        legend.spacing = unit(1.5, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(2.5, "mm"),
        plot.margin = margin(3, 3, 3, 3, "mm"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA)
      )
  }
  model_order <- c("measured", "Zero Scale", "ZSM", "CLR", "vandeputte", "GALAXY", "MetaCardis")
  final_df <- final_df %>%
    mutate(Model = factor(Model, levels = model_order))
  model_levels = final_df %>% distinct(Model) %>% pull(Model)
  offset_range <- 0.3 
  offsets <- seq(from = -offset_range, to = offset_range, length.out = length(model_levels))
  offset_lookup <- setNames(offsets, model_levels)
  plot_df <- final_df %>%
    mutate(
      Study = factor(Study, levels = unique(Study)),
      Model = factor(Model, levels = model_levels),
      model_color = case_when(
        Model == "measured" ~ "#929292",
        Model %in% c("Zero Scale", "CLR", "ZSM") ~ "#781690",
        Model == "vandeputte" ~ "#000000",
        Model == "GALAXY" ~ "#0C58CA",
        Model == "MetaCardis" ~ "#FF8C27",
        TRUE ~ "grey50"
      ),
      arrow_alpha = ifelse(pval < 0.05, 1, 0),
      base_y = 0,
      model_offset = offset_lookup[as.character(Model)]
    ) %>%
    group_by(Study) %>%
    mutate(
      max_abs_beta = max(abs(Beta), na.rm = TRUE),
      max_abs_beta = ifelse(is.finite(max_abs_beta), max_abs_beta, 0)
    ) %>%
    ungroup()
  
  plot_df <- plot_df %>%
    mutate(
      Model = as.character(Model),
      Model = paste0(toupper(substr(Model, 1, 1)), substr(Model, 2, nchar(Model)))
    )
  
  plot_df <- plot_df %>%
    mutate(
      Model = factor(Model, levels = c("Measured", "Zero Scale", "ZSM", "CLR", "Vandeputte", "GALAXY", "MetaCardis"))  
    )
  plot_df$CoefIndex <- as.numeric(factor(plot_df$CoefName))
  x_positions <- sort(unique(plot_df$CoefIndex))
  hline_segments <- data.frame(
    x = x_positions - 0.3,
    xend = x_positions + 0.3,
    y = 0
  )
  plot_df <- plot_df %>%
    mutate(Study = str_replace(Study, "(\\D)(\\d)", "\\1 \\2"))
  p_final <- ggplot(plot_df, aes(x = CoefIndex + model_offset, y = base_y)) +
    geom_segment(
      data = hline_segments,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      color = "grey70",
      linewidth = 0.3
    ) +
    geom_segment(
      aes(
        xend = CoefIndex + model_offset,
        yend = base_y + Beta,
        color = model_color,
        alpha = arrow_alpha,
        linewidth = abs(Beta)
      ),
      arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
      lineend = "round"
    ) +
    geom_point(
      aes(
        fill = model_color
      ),
      size = 4,
      shape = 21,
      color = "white",
      stroke = 0.5
    ) +
    scale_fill_identity(
      guide = "legend",
      name = "Model",
      breaks = c("#929292", "#781690", "#000000", "#0C58CA", "#FF8C27"),
      labels = c("Measured", "Zero Scale/CLR ALDEx3", "Vandeputte 2021", "GALAXY", "MetaCardis")
    ) +
    scale_color_identity(guide = "none") +
    scale_alpha_identity(
      guide = "legend",
      name = "Significance",
      breaks = c(1, 0),
      labels = c("p < 0.05", "p ≥ 0.05")
    ) +
    scale_linewidth_continuous(range = c(0.3, 2), guide = "none") +
    scale_x_continuous(
      breaks = unique(plot_df$CoefIndex),
      labels = unique(plot_df$CoefName),
      expand = expansion(add = c(0.5, 0.5))
    ) +
    facet_wrap(~ Study, ncol = 1, strip.position = "top", scales = "free_y") +
    scale_y_continuous(
      name = "Beta Coefficient",
      breaks = function(limits) {
        max_y <- ceiling(max(abs(limits)) * 1.1 * 2) / 2
        if (max_y == 0) return(c(-1, 0, 1)) 
        seq(-max_y, max_y, length.out = 5)  
      },
      limits = function(limits) {
        max_y <- ceiling(max(abs(limits)) * 1.1 * 2) / 2
        if (max_y == 0) return(c(-1, 1)) 
        c(-max_y, max_y) 
      },
      labels = function(y) sprintf("%.2f", y),
      expand = expansion(mult = c(0.1, 0.1))  
    ) +
    theme_minimal(base_size = 14) +
    labs(x = "Covariate Factors", y = "Beta Coefficient") +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(margin = margin(r=20)),
      axis.title.x = element_text(margin = margin(t=20)),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 16, face = "bold", hjust = 0),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20, 10, 20, 10)
    ) +
    coord_cartesian(
      xlim = c(0.5, length(unique(plot_df$CoefName)) + 0.5),
      clip = "off"
    )
  
  return(list(export_table = final_df, plot = p, betaplot = p_final))
}

modelingscale_arrow_plot <- function(nested_df_list, datasets, model_include = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(grid)
  
  capitalize_first <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  
  flat_list <- list()
  for (study in names(nested_df_list)) {
    model_dfs <- list()
    for (model_nm in names(nested_df_list[[study]])) {
      df <- as.data.frame(nested_df_list[[study]][[model_nm]])
      if (nrow(df) > 0) {
        df$Study <- study
        df$Model <- model_nm
        model_dfs[[model_nm]] <- df
      }
    }
    if (length(model_dfs) > 0) {
      flat_list[[study]] <- do.call(rbind, model_dfs)
    }
  }
  
  if (length(flat_list) == 0) {
    warning("No data to process.")
    return(list(export_table = NULL, plot = NULL))
  }
  
  big_df <- do.call(rbind, flat_list)
  needed_cols <- c("Study", "Model", "variable", "estimate.mean", "p.val")
  for (colname in needed_cols) {
    if (!colname %in% colnames(big_df)) {
      stop("Missing required column: ", colname)
    }
  }
  
  final_df <- big_df %>%
    dplyr::select(all_of(needed_cols)) %>%
    rename(Beta = estimate.mean, pval = p.val) %>%
    mutate(Model = factor(Model), variable = factor(variable)) %>%
    group_by(Study) %>%
    filter(any(Model == "measured")) %>%
    ungroup() %>%
    filter(!grepl("intercept", variable, ignore.case = TRUE))
  
  covariate_cols <- c("covariate1", "covariate2", "covariate3", 
                      "covariate4", "covariate5", "covariate6", "covariate7")
  
  all_studies <- unique(final_df$Study)
  
  study_covariate_list <- purrr::map(
    all_studies,
    function(study_name) {
      subset_df <- datasets[datasets$dataset_name == study_name & datasets$abundance == "metadata", , drop = FALSE]
      if (nrow(subset_df) == 0) return(NA_character_)
      metadatacondition <- unique(subset_df[, intersect(covariate_cols, colnames(subset_df)), drop = FALSE])
      metadatacondition <- metadatacondition[, colSums(!is.na(metadatacondition)) > 0, drop = FALSE]
      metadatacondition <- metadatacondition[, colSums(metadatacondition != "" & metadatacondition != " ") > 0, drop = FALSE]
      if (ncol(metadatacondition) == 0) return(NA_character_)
      return(as.character(unlist(metadatacondition)))
    }
  )
  names(study_covariate_list) <- all_studies
  
  final_df <- final_df %>%
    rowwise() %>%
    mutate(
      var_base = {
        covariates <- study_covariate_list[[Study]]
        matched_base <- covariates[str_detect(as.character(variable), fixed(covariates))]
        if (length(matched_base) > 0) matched_base[1] else NA_character_
      }
    ) %>%
    ungroup()
  
  beta_lookup <- final_df %>%
    distinct(Study, var_base) %>%
    arrange(Study, var_base) %>%
    group_by(Study) %>%
    mutate(BetaIndex = row_number()) %>%
    ungroup()
  
  factor_lookup <- final_df %>%
    distinct(Study, var_base, variable) %>%
    arrange(Study, var_base, variable) %>%
    group_by(Study, var_base) %>%
    mutate(factor_index = row_number()) %>%
    ungroup()
  
  final_df <- final_df %>%
    left_join(beta_lookup, by = c("Study", "var_base")) %>%
    left_join(factor_lookup, by = c("Study", "var_base", "variable")) %>%
    mutate(
      CoefName = paste0("B", BetaIndex, " (", factor_index, ")"),
      CoefNameFull = paste0(Study, "_", CoefName),
      CoefNameFull = factor(CoefNameFull, levels = unique(paste0(Study, "_", CoefName)))
    )
  
  measured_df <- final_df %>%
    filter(Model == "measured") %>%
    dplyr::select(Study, variable, MeasuredBeta = Beta, MeasuredPval = pval)
  
  plot_df <- final_df %>%
    filter(Model != "measured") %>%
    left_join(measured_df, by = c("Study", "variable"))
  
  if (!is.null(model_include)) {
    plot_df <- plot_df %>% filter(Model %in% model_include)
  }
  
  plot_df <- plot_df %>%
    mutate(
      model_color = case_when(
        Model %in% c("Zero Scale", "CLR", "ZSM") ~ "#781690",
        Model == "vandeputte" ~ "#000000",
        Model == "GALAXY" ~ "#0C58CA",
        Model == "MetaCardis" ~ "#FF8C27",
        TRUE ~ "grey50"
      ),
      sig_model = pval < 0.05,
      sig_measured = MeasuredPval < 0.05,
      arrow_alpha = ifelse(sig_model & sig_measured, 1, 0.1)
    )
  
  model_levels <- as.character(unique(plot_df$Model))
  offset_range <- 0.15
  offsets <- seq(from = -offset_range, to = offset_range, length.out = length(model_levels))
  offset_lookup <- setNames(offsets, tolower(model_levels))
  
  ordering_df <- plot_df %>%
    mutate(arrow_direction = MeasuredBeta - Beta) %>%
    group_by(Study, CoefNameFull) %>%
    summarise(avg_arrow = mean(arrow_direction, na.rm = TRUE), .groups = "drop") %>%
    arrange(Study, desc(avg_arrow)) %>%
    group_by(Study) %>%
    mutate(sort_order = row_number()) %>%
    ungroup()
  
  plot_df <- plot_df %>%
    left_join(ordering_df %>% dplyr::select(Study, CoefNameFull, sort_order), by = c("Study", "CoefNameFull")) %>%
    arrange(Study, sort_order) %>%
    mutate(
      CoefNameFull = factor(CoefNameFull, levels = unique(CoefNameFull)),
      Model = capitalize_first(as.character(Model)),
      Model = factor(Model, levels = capitalize_first(model_levels)),
      x_jittered = as.numeric(CoefNameFull) + offset_lookup[tolower(as.character(Model))]
    )
  
  facet_limits <- plot_df %>%
    group_by(Study) %>%
    summarise(ymax = max(abs(c(Beta, MeasuredBeta)), na.rm = TRUE)) %>%
    mutate(ymin = -ymax)
  
  plot_df <- plot_df %>%
    left_join(facet_limits, by = "Study")
  
  p_final <- ggplot(plot_df, aes(x = x_jittered)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    geom_blank(aes(y = ymin)) +
    geom_blank(aes(y = ymax)) +
    
    geom_segment(
      aes(
        y = Beta,
        yend = MeasuredBeta,
        xend = x_jittered,
        color = Model,
        alpha = arrow_alpha
      ),
      linewidth = 1.4,
      lineend = "round",
    ) + 
    
    geom_point(
      aes(
        y = MeasuredBeta,
        color = Model,
        alpha = arrow_alpha
      ),
      size = 2.2,
      shape = 16
    )+
    
    geom_point(
      aes(
        y = Beta,
        fill = Model,
        color = Model,
        shape = sig_model
      ),
      size = 3.8,
      stroke = 1.2
    ) +
    scale_fill_manual(
      name = "Model",
      values = setNames(unique(plot_df$model_color), unique(plot_df$Model))
    ) +
    scale_color_manual(
      name = "Model",
      values = setNames(unique(plot_df$model_color), unique(plot_df$Model))
    ) +
    scale_shape_manual(
      values = c("TRUE" = 21, "FALSE" = 1),
      labels = c("p ≥ 0.05", "p < 0.05"),
      name = "Model Significance"
    ) +
    scale_alpha_identity(
      name = "Agreement",
      breaks = c(1, 0.1),
      labels = c("Both p < 0.05", "Either p ≥ 0.05"),
      guide = "legend"
    ) +
    
    scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      expand = expansion(add = 0.25)
    ) +
    
    facet_wrap(~ Study, ncol = 1, scales = "free", strip.position = "left") +
    labs(y = "Regression Coefficient", x="Covariates") +
    theme_classic(base_size = 18, base_family = "Arial") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x =element_text(size = 18, margin = margin(t=10), hjust = 0.5),
      axis.text.y = element_text(size = 18),
      axis.title.y.placement = "inside",
      axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 1),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20, face = "bold"),
      strip.placement = "outside",
      strip.text.y.left = element_text(size = 20, face = "bold", angle = 90),
      strip.background = element_blank(),
      panel.spacing = unit(1.2, "lines")
    ) + 
    guides(
      shape = guide_legend(
        override.aes = list(
          fill = c("white", "black"), 
          color = "black"
        )
      )
    )
  
  return(list(export_table = plot_df, betaplot = p_final))
}


aldex2_conf_matrix <- function(df_list) {
  compute_conf_matrix <- function(df) {
    true_values <- df %>%
      dplyr::filter(Gamma == "Measured_1e-09") %>%
      dplyr::select(Feature, kw.eBH) %>%
      dplyr::rename(True_Value = kw.eBH)

    df_wide <- df %>%
      dplyr::filter(Gamma != "Measured_1e-09") %>%
      pivot_wider(names_from = Gamma, values_from = kw.eBH) %>%
      left_join(true_values, by = "Feature")

    df_wide <- df_wide %>%
      dplyr::mutate(True_Value = ifelse(True_Value > 0.05, 0, 1))

    models <- setdiff(names(df_wide), c("Feature", "True_Value", "Dataset",
                                        "Externalscalemeasurement", "Uncertainty",
                                        "kw.ep", "glm.ep", "glm.eBH"))

    conf_matrices <- lapply(models, function(model) {
      if (!all(is.na(df_wide[[model]]))) {
        pred <- ifelse(df_wide[[model]] > 0.05, 0, 1)

        true_factor <- factor(df_wide$True_Value, levels = c(0, 1))
        pred_factor <- factor(pred, levels = c(0, 1))

        if (length(unique(true_factor)) < 2 || length(unique(pred_factor)) < 2) {
          return(NULL)
        }

        cm <- as.data.frame(table(Predicted = pred_factor, Actual = true_factor))
        total_samples <- sum(cm$Freq)
        cm <- cm %>%
          mutate(Percent = (Freq / total_samples) * 100)

        return(cm)
      } else {
        return(NULL)
      }
    })

    names(conf_matrices) <- models
    return(conf_matrices)
  }

  results <- lapply(df_list, compute_conf_matrix)
  prepare_cm_for_plot <- function(cm, model_name, dataset_name) {
    if (is.null(cm)) {
      return(NULL)
    }

    cm_df <- as.data.frame(cm)
    colnames(cm_df) <- c("Prediction", "Reference", "Count", "Percent")
    cm_df$Model <- model_name
    cm_df$Dataset <- dataset_name

    return(cm_df)
  }

  all_cm_data <- do.call(rbind, lapply(seq_along(results), function(i) {
    dataset_name <- names(results)[i] 

    do.call(rbind, lapply(names(results[[i]]), function(model) {
      prepare_cm_for_plot(results[[i]][[model]], model, dataset_name)
    }))
  }))

  all_cm_data <- all_cm_data %>% filter(!is.na(Percent))
  all_cm_data$Reference <- as.factor(all_cm_data$Reference)
  all_cm_data$Prediction <- as.factor(all_cm_data$Prediction)

  cm_plot <- ggplot(all_cm_data, aes(x = Reference, y = Prediction, fill = Percent)) +
    geom_tile(color = "black") +
    geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = 0.5, size = 5, fontface = "bold", color = "white") +
    scale_fill_gradient(low = "#f7fbff", high = "#1D4779", name = "Percent") + 
    facet_grid(Dataset ~ Model, scales = "free", space = "free") + 
    labs(
      title = "Confusion Matrices for All Models & Datasets",
      subtitle = "Percentage of Predictions Compared to True Values (Measured_1e-09)",
      x = "Actual (True Label)",
      y = "Predicted (Model Output)"
    ) +
    theme_minimal(base_size = 14) +  
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold", color = "white"),
      strip.background = element_rect(fill = "#2c3e50", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm")
    )

  return(list(conf_matrices = results, plot = cm_plot))
}

aldex3_conf_matrix <- function(
    nested_df_list,
    categorical_df,
    true_model = "Measured_0.5", 
    compare_models = NULL,       
    alpha = 0.05,
    plot_value = "percent",
    split_by = c("custom", "sequencingtype")
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(ggh4x)
  
  split_by <- match.arg(split_by)
  
  # Define default model sets
  galaxy_metacardis_models <- c(
    "GALAXY_0", "GALAXY_0.5", "GALAXY_1", "GALAXY_5",
    "MetaCardis_0", "MetaCardis_0.5", "MetaCardis_1", "MetaCardis_5",
    "ALDEx3_tss_0", "ALDEx3_tss_0.5", "ALDEx3_tss_1", "ALDEx3_tss_5",
    "ALDEx2_clr_0", "ALDEx2_clr_0.5", "ALDEx2_clr_1", "ALDEx2_clr_5",
    "DESeq2_0", "limma_0"
  )
  
  vandeputte_models <- c(
    "vandeputte_0", "vandeputte_0.5", "vandeputte_1", "vandeputte_5",
    "ALDEx3_tss_0", "ALDEx3_tss_0.5", "ALDEx3_tss_1", "ALDEx3_tss_5",
    "ALDEx2_clr_0", "ALDEx2_clr_0.5", "ALDEx2_clr_1", "ALDEx2_clr_5",
    "DESeq2_0", "limma_0"
  )
  
  metric_order <- c("TP", "TN", "FP", "FN", "PPV", "NPV", "Sensitivity", "FDR")
  
  flat_list <- list()
  for (study in names(nested_df_list)) {
    model_dfs <- list()
    for (model in names(nested_df_list[[study]])) {
      df <- as.data.frame(nested_df_list[[study]][[model]])
      if (nrow(df) > 0) {
        required_cols <- c("feature", "covariate", "p.val.adj", "estimate")
        if (!all(required_cols %in% colnames(df))) {
          warning(paste("Model", model, "in study", study, "missing required columns. Skipping."))
          next
        }
        
        model_parts <- unlist(strsplit(model, "_"))
        scale_label <- model_parts[1]
        gamma_value <- if (length(model_parts) > 1) model_parts[2] else "1e09"
        
        df <- df %>%
          mutate(
            scale_label = if ("scale_label" %in% colnames(df)) scale_label else scale_label,
            gamma = if ("gamma" %in% colnames(df)) gamma else gamma_value
          ) %>%
          dplyr::select(feature, covariate, p.val.adj, estimate, scale_label, gamma)
        
        model_dfs[[model]] <- df
      }
    }
    if (length(model_dfs) > 0) {
      flat_list[[study]] <- do.call(rbind, model_dfs)
    }
  }
  
  compute_conf_matrices <- function(df, threshold, true_model) {
    library(dplyr)
    
    df <- df %>%
      mutate(gamma_combined = paste0(scale_label, "_", gsub("-", "", gamma)))
    
    gold <- df %>%
      filter(gamma_combined == true_model) %>% 
      dplyr::select(feature, covariate, gold_pval = p.val.adj, gold_sign = estimate) %>%
      distinct(feature, covariate, .keep_all = TRUE) %>%
      mutate(gold_pval = ifelse(is.na(gold_pval), 1, gold_pval))
    
    if (nrow(gold) == 0) return(NULL)
    
    model_ids <- unique(df$gamma_combined[df$gamma_combined != true_model]) # MODIFIED: Exclude true_model
    conf_matrices <- list()
    
    for (model_id in model_ids) {
      model_df <- df %>%
        filter(gamma_combined == model_id) %>%
        dplyr::select(feature, covariate, model_pval = p.val.adj, model_sign = estimate) %>%
        distinct(feature, covariate, .keep_all = TRUE) %>%
        mutate(model_pval = ifelse(is.na(model_pval), 1, model_pval))
      
      joined <- inner_join(model_df, gold, by = c("feature", "covariate"))
      if (nrow(joined) == 0) next
      
      pred_sig <- joined$model_pval <= threshold
      true_sig <- joined$gold_pval <= threshold
      sign_match <- sign(joined$model_sign) == sign(joined$gold_sign)
      
      TP_vec <- true_sig & pred_sig & sign_match
      TN_vec <- !true_sig & !pred_sig 
      FP_vec <- (pred_sig & !true_sig) | (pred_sig & true_sig & !sign_match)
      FN_vec <- true_sig & !pred_sig
      
      TP <- sum(TP_vec)
      TN <- sum(TN_vec)
      FP <- sum(FP_vec)
      FN <- sum(FN_vec)
      total <- TP + TN + FP + FN
      
      sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
      PPV <- if ((TP + FP) > 0) TP / (TP + FP) else NA
      NPV <- if ((TN + FN) > 0) TN / (TN + FN) else NA
      FDR <- if ((TP + FP) > 0) FP / (TP + FP) else NA
      
      cm_df <- data.frame(
        Predicted = factor(rep(c(0, 1), each = 2), levels = c(0, 1)),
        Actual = factor(rep(c(0, 1), times = 2), levels = c(0, 1)),
        Count = c(TN, FN, FP, TP),
        Percent = c(TN, FN, FP, TP) / total * 100
      )
      
      conf_matrices[[model_id]] <- list(
        plot_df = cm_df,
        metrics = data.frame(
          TP = TP, TN = TN, FP = FP, FN = FN,
          PPV = PPV * 100, NPV = NPV * 100,
          Sensitivity = sensitivity * 100, FDR = FDR * 100
        )
      )
    }
    return(conf_matrices)
  }
  
  results <- lapply(flat_list, function(x) compute_conf_matrices(x, threshold = alpha, true_model = true_model)) # MODIFIED: Pass true_model
  
  all_rows <- list()
  for (study in names(results)) {
    for (model in names(results[[study]])) {
      met <- results[[study]][[model]]$metrics
      if (!is.null(met)) {
        met$Model <- model
        all_rows[[length(all_rows) + 1]] <- data.frame(Study = study, met)
      }
    }
  }
  if (length(all_rows) == 0) {
    warning("No metrics to plot. Returning early.")
    return(list(conf_matrices = results, export_table = NULL, subset_plots = NULL))
  }
  
  export_df <- do.call(rbind, all_rows)
  
  export_table <- export_df %>%
    arrange(Study, Model) %>%
    group_by(Study) %>%
    nest(metrics = c(Model, TP, TN, FP, FN, PPV, NPV, Sensitivity, FDR))
  
  df_value <- export_df %>%
    pivot_longer(cols = c("TP","TN","FP","FN","PPV","NPV","Sensitivity","FDR"),
                 names_to = "Metric", values_to = "Value")
  
  df_percent <- export_df %>%
    rowwise() %>%
    mutate(
      total = TP + TN + FP + FN,
      TP = ifelse(total > 0, (TP/total)*100, NA),
      TN = ifelse(total > 0, (TN/total)*100, NA),
      FP = ifelse(total > 0, (FP/total)*100, NA),
      FN = ifelse(total > 0, (FN/total)*100, NA)
    ) %>%
    ungroup() %>%
    dplyr::select(-total) %>%
    pivot_longer(cols = c("TP","TN","FP","FN","PPV","NPV","Sensitivity","FDR"),
                 names_to = "Metric", values_to = "Percent")
  
  metrics_long <- df_value %>%
    left_join(df_percent, by = c("Study","Model","Metric")) %>%
    filter(!grepl(paste0("^", true_model, "$"), Model)) # MODIFIED: Exclude true_model dynamically
  
  # MODIFIED: Filter metrics_long for plotting based on compare_models
  if (!is.null(compare_models)) {
    metrics_long <- metrics_long %>%
      filter(Model %in% compare_models)
  }
  
  get_studies <- function(filter_expr) {
    categorical_df %>%
      filter({{ filter_expr }}) %>%
      distinct(Study) %>%
      pull(Study)
  }
  
  study_groups <- list()
  
  if (split_by == "custom") {
    study_groups <- list(
      human_fecal_fc_16s = list(
        studies = get_studies(
          organismtype == "Human" &
            sampletype == "Fecal" &
            loadtype == "Flow Cytometry" &
            sequencingtype == "16S rRNA"
        ),
        subtitle = "Human Fecal Samples, Flow Cytometry, 16S rRNA Sequencing"
      ),
      human_fecal_fc_shotgun = list(
        studies = get_studies(
          organismtype == "Human" &
            sampletype == "Fecal" &
            loadtype == "Flow Cytometry" &
            sequencingtype == "shotgun metagenomics"
        ),
        subtitle = "Human Fecal Samples, Flow Cytometry, Shotgun Metagenomics"
      ),
      not_both_human_fecal_shotgun = list(
        studies = get_studies(
          !(organismtype == "Human" & sampletype == "Fecal" & loadtype == "Flow Cytometry") &
            sequencingtype == "shotgun metagenomics"
        ),
        subtitle = "Not (Human, Fecal, Flow Cytometry), Shotgun Metagenomics"
      ),
      not_both_human_fecal_16s = list(
        studies = get_studies(
          !(organismtype == "Human" & sampletype == "Fecal" & loadtype == "Flow Cytometry") &
            sequencingtype == "16S rRNA"
        ),
        subtitle = "Not (Human, Fecal, Flow Cytometry), 16S rRNA Sequencing"
      )
    )
  } else if (split_by == "sequencingtype") {
    seqtypes <- unique(categorical_df$sequencingtype)
    for (stype in seqtypes) {
      studies <- get_studies(sequencingtype == stype)
      if (length(studies) > 0) {
        study_groups[[stype]] <- list(
          studies = studies,
          subtitle = paste("Sequencing Type:", stype)
        )
      }
    }
  }
  
  tidyplot_theme <- function(base_size = 5) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        axis.title.x = element_text(margin = margin(t = 10), hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 10), hjust = 0.5),
        axis.title = element_text(size = 14, face = "plain"),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  }
  
  produce_plot <- function(data, model_order, plot_title, subtitle_info, true_model) { # MODIFIED: Added true_model parameter
    metric_order <- c("TP", "TN", "FP", "FN", "PPV", "NPV", "Sensitivity", "FDR")
    data$Metric <- factor(data$Metric, levels = metric_order)
    model_ranks <- data %>%
      filter(Metric == "FDR") %>%
      group_by(Model) %>%
      summarise(avg_fdr = mean(Value, na.rm = TRUE)) %>%
      arrange(desc(avg_fdr)) %>%
      pull(Model)
    
    data$Model <- factor(data$Model, levels = model_ranks)
    
    study_levels <- unique(data$Study)
    y_scales <- c(
      list(scale_y_discrete()),
      rep(list(scale_y_discrete(labels = NULL, breaks = NULL)), length(study_levels) - 1)
    )
    
    ggplot(data, aes(x = Metric, y = Model, fill = Percent)) +
      geom_tile(color = "black") +
      geom_text(
        aes(label = ifelse(Metric %in% c("TP", "TN", "FP", "FN"),
                           sprintf("%d", as.integer(Value)),
                           sprintf("%.2f", Value))),
        size = 4,
        fontface = "bold"
      ) +
      scale_fill_gradient(low = "#f7fbff", high = "#1F77B4", name = "Percent") +
      scale_x_discrete(position = "top") +
      scale_color_identity() +
      facet_wrap(~ Study, nrow = 1, scales = "free_y", strip.position = "top") +
      facetted_pos_scales(y = y_scales) +
      labs(
        title = plot_title,
        subtitle = paste("Truth =", true_model, ", FDR Threshold <", alpha), # MODIFIED: Dynamic true_model in subtitle
        x = NULL, y = NULL
      ) +
      tidyplot_theme() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray50", margin = margin(b = 20)),
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
        strip.text = element_text(size = 14, face = "bold"),
        strip.placement = "outside",
        panel.spacing = unit(0.3, "cm")
      )
  }
  
  plots <- list()
  for (group_name in names(study_groups)) {
    subset_data <- metrics_long %>% filter(Study %in% study_groups[[group_name]]$studies)
    if (nrow(subset_data) > 0) {
      model_order <- if (any(grepl("^GALAXY_|^MetaCardis_", subset_data$Model))) {
        galaxy_metacardis_models
      } else {
        vandeputte_models
      }

      if (!is.null(compare_models)) {
        model_order <- model_order[model_order %in% compare_models]
      }
      plots[[group_name]] <- produce_plot(
        data = subset_data,
        model_order = model_order,
        plot_title = paste("Differential Abundance:", study_groups[[group_name]]$subtitle),
        subtitle_info = study_groups[[group_name]]$subtitle,
        true_model = true_model
      )
    } else {
      plots[[group_name]] <- NULL
    }
  }
  
  return(list(
    conf_matrices = results,
    export_table = export_table,
    subset_plots = plots,
    metrics_long = metrics_long
  ))
}
aldex_scatterplot <- function(data, split_by_seqtype = FALSE, metrics = c("PPV", "NPV", "Sensitivity", "FDR", "FNR", "FPR", "FOR")) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  tidyplot_theme <- function(base_size = 16) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "#333333"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, color = "#666666", hjust = 0.5),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "top",
        panel.grid.major = element_line(linewidth = 0.4, color = "#DDDDDD"),
        panel.grid.minor = element_blank()
      )
  }
  
  make_plot_grid <- function(df_subset, seqtype_label = NULL) {
    df_wide <- df_subset %>%
      filter(Metric %in% metrics, !is.na(Value)) %>%
      dplyr::select(Study, model_type, Metric, Value) %>%
      pivot_wider(names_from = Metric, values_from = Value, values_fn = list(Value = mean)) %>%
      filter(complete.cases(.))
    
    if (nrow(df_wide) == 0) {
      warning("No complete cases for plotting after filtering NAs.")
      return(NULL)
    }
    
    models <- unique(df_wide$model_type)
    color_values <- setNames(
      c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")[seq_along(models)],
      models
    )
    
    message("Model-to-color mapping:")
    for (model in names(color_values)) {
      message(sprintf("  %s: %s", model, color_values[model]))
    }
    
    metric_pairs <- combn(metrics, 2, simplify = FALSE)
    
    plot_list <- lapply(metric_pairs, function(pair) {
      ggplot(df_wide, aes_string(x = pair[1], y = pair[2], color = "model_type")) +
        geom_point(alpha = 0.8, size = 2.5, position = position_jitter(width = 0.02, height = 0.02)) +
        scale_color_manual(values = color_values, name = "Model") +
        labs(x = pair[1], y = pair[2]) +
        tidyplot_theme()
    })
    
    title_text <- "Pairwise Comparison of Performance Metrics"
    if (!is.null(seqtype_label)) title_text <- paste0(title_text, " (", seqtype_label, ")")
    
    # Compose the grid and add title/subtitle
    (wrap_plots(plotlist = plot_list, ncol = 3, guides = "collect") &
        tidyplot_theme() &
        plot_annotation(
          title = title_text,
          subtitle = expression("Compared to Truth (Measured with "~log[2]~"Standard Deviation = 0.5)"),
          theme = tidyplot_theme()
        )) +
      plot_layout(guides = "collect")
  }
  
  if (split_by_seqtype) {
    seqtypes <- unique(data$sequencingtype)
    plot_list <- lapply(seqtypes, function(type) {
      make_plot_grid(subset(data, sequencingtype == type), seqtype_label = type)
    })
    names(plot_list) <- seqtypes
    return(plot_list)
  } else {
    plot <- make_plot_grid(data)
    if (!is.null(plot)) return(plot) else stop("No plot generated due to insufficient data.")
  }
}

aldex_boxplot <- function(data, split_by_seqtype = FALSE, metrics = c("PPV", "NPV", "Sensitivity", "FDR"), connect_dots_by = NULL) {
  library(scales)
  library(grid)
  library(dplyr)
  library(ggplot2)
  
  if (!"sequencingtype" %in% names(data)) {
    stop("Column `sequencingtype` not found in data.")
  }
  
  if (!"Study" %in% names(data)) {
    stop("Column `Study` not found in data.")
  }
  
  if (!all(metrics %in% unique(data$Metric))) {
    stop("Some specified metrics are not present in the data: ", 
         paste(setdiff(metrics, unique(data$Metric)), collapse = ", "))
  }
  
  if (!is.null(connect_dots_by) && !connect_dots_by %in% names(data)) {
    stop("Column specified in `connect_dots_by` (", connect_dots_by, ") not found in data.")
  }
  
  make_plot <- function(df_subset, seqtype_label = NULL) {
    
    df_subset <- df_subset %>%
      filter(Metric %in% metrics, !is.na(Value))
    
    # Sort models by median Value (averaged across metrics)
    model_order <- df_subset %>%
      group_by(model_type, Metric) %>%
      summarise(median_val = median(Value, na.rm = TRUE), .groups = "drop") %>%
      group_by(model_type) %>%
      summarise(avg_median = mean(median_val, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(avg_median)) %>%
      pull(model_type)
    
    df_subset$model_type <- factor(df_subset$model_type, levels = model_order)
    
    # Define dynamic aesthetic values based on unique studies
    studies <- unique(df_subset$Study[!is.na(df_subset$Study)])
    color_values <- setNames(
      c("black", "#EF9129", "#9932CC", "#004D40", "#489BE4")[seq_along(studies)],
      studies
    )
    
    title_text <- "Differential Abundance Performance Metrics by Model"
    if (!is.null(seqtype_label)) {
      title_text <- paste0(title_text, "\n(", seqtype_label, ")")
    }
    
    tidyplot_theme <- function(base_size = 16) {  
      theme_minimal(base_size = base_size) +
        theme(
          text = element_text(family = "Arial", color = "black"),
          axis.text = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.line = element_line(linewidth = 0.5, color = "black"),
          axis.ticks.length = unit(0.25, "cm"),  
          axis.ticks = element_line(linewidth = 0.5, color = "black"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          plot.subtitle = element_text(size = 20, hjust = 0.5)
        )
    }
    
    # Determine facet layout based on number of metrics
    n_metrics <- length(unique(df_subset$Metric))
    if (n_metrics <= 3) {
      facet_layout <- list(ncol = 3, nrow = 1)
    } else if (n_metrics == 4) {
      facet_layout <- list(ncol = 2, nrow = 2)
    } else {
      facet_layout <- list(ncol = 3, nrow = ceiling(n_metrics / 3))
    }
    
    # Initialize ggplot
    plot <- ggplot(df_subset, aes(x = model_type, y = Value))
    
    # Add lines if connect_dots_by is specified
    if (!is.null(connect_dots_by)) {
      plot <- plot + geom_line(
        aes(group = interaction(.data[[connect_dots_by]], Metric), 
            color = Study),
        linewidth = 0.5,
        alpha = 0.5,
        show.legend = FALSE
      )
    }
    
    # Add boxplots and jittered points
    plot <- plot +
      geom_boxplot(
        color = "black",
        fill = NA,
        size = 1,
        outlier.shape = NA,
        show.legend = FALSE
      ) +
      geom_jitter(
        aes(color = Study),
        width = 0.2,
        size = 4,
        shape = 16,
        stroke = 1,
        show.legend = TRUE
      ) +
      scale_color_manual(
        name = "Study",
        values = color_values,
        guide = guide_legend(
          override.aes = list(
            shape = 16,
            size = 4,
            stroke = 1
          )
        )
      ) +
      scale_x_discrete(
        drop = FALSE,
        expand = expansion(add = c(0.8, 0.8))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      tidyplot_theme(base_size = 16) +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 20, face = "bold", angle = 0, hjust = 1),
        strip.text.x = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        axis.ticks.length.x = unit(0.15, "cm"),
        panel.spacing.y = unit(1.2, "lines"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "vertical"
      ) +
      labs(
        title = title_text,
        subtitle = bquote("Compared to Truth (Measured with " ~ Log[2] ~ "Standard Deviation = 0.5" ~ ")"),
        x = "Model",
        y = "Percentage (%)"
      )
    
    # Apply faceting
    plot <- plot + facet_wrap(
      ~ Metric,
      ncol = facet_layout$ncol,
      nrow = facet_layout$nrow,
      scales = "free_y"
    )
    
    plot
  }
  
  if (split_by_seqtype) {
    seqtypes <- unique(data$sequencingtype)
    plot_list <- lapply(seqtypes, function(type) {
      make_plot(subset(data, sequencingtype == type), seqtype_label = type)
    })
    names(plot_list) <- seqtypes
    return(plot_list)
  } else {
    return(make_plot(data))
  }
}

aldex_lineplot <- function(data, split_by_seqtype = FALSE, metrics = c("PPV", "NPV", "Sensitivity", "FDR")) {
  library(scales)
  library(grid)
  library(dplyr)
  library(ggplot2)
  
  if (!"sequencingtype" %in% names(data) && split_by_seqtype) {
    stop("Column `sequencingtype` not found in data.")
  }
  
  if (!all(metrics %in% unique(data$Metric))) {
    stop("Some specified metrics are not present in the data: ", 
         paste(setdiff(metrics, unique(data$Metric)), collapse = ", "))
  }
  
  make_plot <- function(df_subset, seqtype_label = NULL) {
    
    df_subset <- df_subset %>%
      mutate(model_type = case_when(
        grepl("GALAXY", model_type) ~ "GALAXY/MetaCardis",
        grepl("MetaCardis", model_type) ~ "GALAXY/MetaCardis",
        TRUE ~ model_type
      )) %>%
      filter(Metric %in% metrics)

    df_avg <- df_subset %>%
      filter(!is.na(Value)) %>%
      group_by(gamma, model_type, Metric) %>%
      summarise(
        N = n(),
        meanVal = mean(Value, na.rm = TRUE),
        SD = if (N > 1) sd(Value, na.rm = TRUE) else NA_real_,
        SE = if (N > 1) SD / sqrt(N) else NA_real_,
        CI_lower = if (N > 1) meanVal - qt(0.975, N - 1) * SE else NA_real_,
        CI_upper = if (N > 1) meanVal + qt(0.975, N - 1) * SE else NA_real_,
        .groups = "drop"
      ) %>%
      mutate(Value = meanVal) %>%
      dplyr::select(-meanVal) %>%
      mutate(gamma = ifelse(gamma == "1e09", "<1e-09", as.character(gamma)))
    
    df_avg$gamma <- factor(df_avg$gamma, levels = sort(unique(df_avg$gamma)))
    
    # Define dynamic aesthetic values based on unique model types
    model_types <- unique(df_avg$model_type)
    color_values <- setNames(
      c("black", "#EF9129", "#9932CC", "#004D40", "#489BE4")[seq_along(model_types)],
      model_types
    )
    fill_values <- color_values
    
    # If you want DESeq2 to have dashed lines
    linetype_values <- setNames(rep("solid", length(model_types)), model_types)
    if ("DESeq2" %in% model_types) {
      linetype_values["DESeq2"] <- "dashed"
    }
    
    # Shift DESeq2 and limma horizontally
    nudge_de_limma <- position_nudge(x = 0.1)
    
    title_text <- "Differential Abundance Performance Metrics by Model"
    if (!is.null(seqtype_label)) {
      title_text <- paste0(title_text, "\n(", seqtype_label, ")")
    }
    
    tidyplot_theme <- function(base_size = 14) {  
      theme_minimal(base_size = base_size) +
        theme(
          text = element_text(family = "Arial", color = "black"),
          axis.text = element_text(size = 16, color = "black", face = "bold"),
          axis.line = element_line(linewidth = 0.5, color = "black"),
          axis.ticks.length = unit(0.25, "cm"),  
          axis.ticks = element_line(linewidth = 0.5, color = "black"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          plot.subtitle = element_text(size = 20, hjust = 0.5)
        )
    }
    
    ggplot(df_avg, aes(x = gamma, y = Value)) +
      geom_ribbon(
        aes(ymin = Value - SE, ymax = Value + SE, fill = model_type, group = model_type),
        alpha = 0.1, color = NA,
        show.legend = FALSE
      ) +
      geom_line(
        aes(color = model_type, linetype = model_type, group = model_type),
        linewidth = 2,
        show.legend = FALSE
      ) +
      # Points for models that are NOT DESeq2 or limma
      geom_point(
        data = subset(df_avg, !model_type %in% c("DESeq2", "Limma")),
        aes(fill = model_type),
        size = 5,
        shape = 21,
        stroke = 1,
        color = "black",
        show.legend = TRUE
      ) +
      # Points + errorbars for DESeq2 and limma, nudged right
      geom_point(
        data = subset(df_avg, model_type %in% c("DESeq2", "Limma")),
        aes(fill = model_type),
        size = 5,
        shape = 21,
        stroke = 1,
        color = "black",
        show.legend = TRUE
      ) +
      geom_errorbar(
        data = subset(df_avg, model_type %in% c("DESeq2", "Limma")),
        aes(ymin = CI_lower, ymax = CI_upper, color = model_type),
        position = nudge_de_limma,
        width = 0,
        alpha = 0.2,
        linewidth = 10,
        show.legend = FALSE
      ) +
      scale_color_manual(
        values = color_values,
        guide = "none"
      ) +
      scale_fill_manual(
        name = "Model",
        values = fill_values,
        guide = guide_legend(
          override.aes = list(
            shape = 21, 
            size = 5, 
            color = "black",
            linetype = 0
          )
        )
      ) +
      scale_linetype_manual(values = linetype_values, guide = "none") +
      facet_wrap(~ Metric, ncol = 2, nrow = 2, scales = "free_y") +
      scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.05, 0.05))) +
      tidyplot_theme(base_size = 16) +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 20, face = "bold", angle = 0, hjust = 1),
        strip.text.x = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        axis.ticks.length.x = unit(0.15, "cm"),
        panel.spacing.y = unit(1.2, "lines"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "vertical"
      ) +
      labs(
        title = title_text,
        subtitle = bquote("Compared to Truth (Measured with " ~ Log[2] ~ "Standard Deviation = 0.5" ~ ")"),
        x = bquote("Uncertainty in Scale (" ~ Log[2] ~ " Standard Deviation" ~ ")"),
        y = "Percentage (%)"
      )
  }
  
  if (split_by_seqtype) {
    seqtypes <- unique(data$sequencingtype)
    plot_list <- lapply(seqtypes, function(type) {
      make_plot(subset(data, sequencingtype == type), seqtype_label = type)
    })
    names(plot_list) <- seqtypes
    return(plot_list)
  } else {
    return(make_plot(data))
  }
}

aldex_metriclineplot <- function(data, split_by_seqtype = FALSE, metrics = c("PPV", "NPV", "FDR")) {
  library(scales)
  library(grid)
  library(dplyr)
  library(ggplot2)
  
  if (!"sequencingtype" %in% names(data) && split_by_seqtype) {
    stop("Column `sequencingtype` not found in data.")
  }
  
  if (!all(metrics %in% unique(data$Metric))) {
    stop("Some specified metrics are not present in the data: ", 
         paste(setdiff(metrics, unique(data$Metric)), collapse = ", "))
  }
  
  make_plot <- function(df_subset, seqtype_label = NULL) {
    
    df_subset <- df_subset %>%
      filter(Metric %in% metrics)
    
    df_avg <- df_subset %>%
      filter(!is.na(Value)) %>%
      group_by(Metric, model_type) %>%
      summarise(
        N = n(),
        meanVal = mean(Value, na.rm = TRUE),
        SD = if (N > 1) sd(Value, na.rm = TRUE) else NA_real_,
        SE = if (N > 1) SD / sqrt(N) else NA_real_,
        .groups = "drop"
      ) %>%
      mutate(Value = meanVal) %>%
      dplyr::select(-meanVal)
    
    df_avg$Metric <- factor(df_avg$Metric, levels = metrics)
    
    # Define dynamic aesthetic values based on unique model types
    model_types <- unique(df_avg$model_type)
    color_values <- setNames(
      c("black", "#EF9129", "#9932CC", "#004D40", "#489BE4")[seq_along(model_types)],
      model_types
    )
    fill_values <- color_values
    
    # DESeq2 with dashed lines
    linetype_values <- setNames(rep("solid", length(model_types)), model_types)
    if ("DESeq2" %in% model_types) {
      linetype_values["DESeq2"] <- "dashed"
    }
    
    title_text <- "Differential Abundance Performance Metrics by Model"
    if (!is.null(seqtype_label)) {
      title_text <- paste0(title_text, "\n(", seqtype_label, ")")
    }
    
    tidyplot_theme <- function(base_size = 14) {  
      theme_minimal(base_size = base_size) +
        theme(
          text = element_text(family = "Arial", color = "black"),
          axis.text = element_text(size = 16, color = "black", face = "bold"),
          axis.line = element_line(linewidth = 0.5, color = "black"),
          axis.ticks.length = unit(0.25, "cm"),  
          axis.ticks = element_line(linewidth = 0.5, color = "black"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          plot.subtitle = element_text(size = 20, hjust = 0.5)
        )
    }
    
    ggplot(df_avg, aes(x = Metric, y = Value)) +
      geom_ribbon(
        aes(ymin = Value - SE, ymax = Value + SE, fill = model_type, group = model_type),
        alpha = 0.1, color = NA,
        show.legend = FALSE
      ) +
      geom_line(
        aes(color = model_type, linetype = model_type, group = model_type),
        linewidth = 2,
        show.legend = FALSE
      ) +
      geom_point(
        aes(fill = model_type),
        size = 5,
        shape = 21,
        stroke = 1,
        color = "black",
        show.legend = TRUE
      ) +
      scale_color_manual(
        values = color_values,
        guide = "none"
      ) +
      scale_fill_manual(
        name = "Model",
        values = fill_values,
        guide = guide_legend(
          override.aes = list(
            shape = 21, 
            size = 5, 
            color = "black",
            linetype = 0
          )
        )
      ) +
      scale_linetype_manual(values = linetype_values, guide = "none") +
      scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.05, 0.05))) +
      tidyplot_theme(base_size = 16) +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 20, face = "bold", angle = 0, hjust = 1),
        strip.text.x = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        axis.ticks.length.x = unit(0.15, "cm"),
        panel.spacing.y = unit(1.2, "lines"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.box = "vertical"
      ) +
      labs(
        title = title_text,
        subtitle = bquote("Compared to Truth (Measured with " ~ Log[2] ~ "Standard Deviation = 0.5" ~ ")"),
        x = "Metric",
        y = "Percentage (%)"
      )
  }
  
  if (split_by_seqtype) {
    seqtypes <- unique(data$sequencingtype)
    plot_list <- lapply(seqtypes, function(type) {
      make_plot(subset(data, sequencingtype == type), seqtype_label = type)
    })
    names(plot_list) <- seqtypes
    return(plot_list)
  } else {
    return(make_plot(data))
  }
}


plot_jaccard_heatmap <- function(
  subset_df,
  height = 20,
  width = 20,
  prefix_label = "16S rRNA Amplicon",
  use_consortium = TRUE
  ) {
  
  tidyplot_theme <- function(base_size = 8) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "black"),
        plot.title = element_text(
          size = base_size * 1.4,        
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.subtitle = element_text(
          size = base_size * 1.25,       
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        axis.title = element_text(
          size = base_size * 1.2,        
          face = "plain"
        ),
        axis.title.x = element_text(
          margin = margin(t = 8), 
          hjust = 0.5
        ),
        axis.title.y = element_text(
          margin = margin(r = 8), 
          hjust = 0.5
        ),
        axis.text = element_text(
          size = base_size,
          color = "black"
        ),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = base_size * 1.2),
        legend.text  = element_text(size = base_size)
      )
  }
  
  num_studies <- nrow(subset_df)
  similarity_matrix   <- matrix(0, nrow = num_studies, ncol = num_studies)
  shared_authors_mat  <- matrix(0, nrow = num_studies, ncol = num_studies)
  rownames(similarity_matrix)  <- subset_df$Study.Identifier
  colnames(similarity_matrix)  <- subset_df$Study.Identifier
  rownames(shared_authors_mat) <- subset_df$Study.Identifier
  colnames(shared_authors_mat) <- subset_df$Study.Identifier
  combined_authors_list <- mapply(function(auth_str, cons_str) {
    authors <- str_split(auth_str, ",\\s*")[[1]]
    consortium <- if (use_consortium) str_split(cons_str, ",\\s*")[[1]] else character(0)
    unique(c(authors, consortium))
  }, subset_df$Authors, subset_df$Consortium, SIMPLIFY = FALSE)
  for (i in seq_len(num_studies)) {
    for (j in i:num_studies) {
      A <- combined_authors_list[[i]]
      B <- combined_authors_list[[j]]
      shared <- length(intersect(A, B))
      total_unique <- length(unique(c(A, B)))
      similarity_percent <- ifelse(total_unique == 0, 0, (shared / total_unique) * 100)
      similarity_matrix[i, j]  <- similarity_percent
      similarity_matrix[j, i]  <- similarity_percent
      shared_authors_mat[i, j] <- shared
      shared_authors_mat[j, i] <- shared
    }
  }
  dist_mat <- as.dist(1 - (similarity_matrix / 100))
  hc <- hclust(dist_mat, method = "ward.D2")
  new_order <- hc$labels[hc$order]
  similarity_matrix   <- similarity_matrix[new_order, new_order]
  shared_authors_mat  <- shared_authors_mat[new_order, new_order]
  sim_data    <- melt(similarity_matrix)
  colnames(sim_data) <- c("Study_A", "Study_B", "Similarity")
  shared_data <- melt(shared_authors_mat)
  colnames(shared_data) <- c("Study_A", "Study_B", "Shared_Authors")
  heatmap_data <- merge(sim_data, shared_data, by = c("Study_A", "Study_B"))
  subtitle_text <- paste0(
    "Jaccard Index (%) with Shared Author Counts (Authors",
    if (use_consortium) " + Consortium" else "",
    ", Flow Cytometry Studies)"
  )
  heatmap_plot <- ggplot(heatmap_data, aes(x = Study_A, y = Study_B, fill = Similarity)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white",
      high = tidyplot_colors[2],
      name = "Jaccard\nSimilarity (%)"  
    ) +
    geom_text(aes(label = Shared_Authors), color = tidyplot_colors[1], size = 14) +
    tidyplot_theme() +
    labs(
      title = paste0("Author Similarity Heatmap (", prefix_label, ")"),
      subtitle = subtitle_text,
      x = NULL,
      y = NULL
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 26, margin = margin(b = 10)),
      legend.text = element_text(size = 24),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 30),
      axis.text.y = element_text(hjust = 1, size = 30),
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 28, hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_fixed()
  out_name <- paste0("plots/Author_similarity_heatmap_", gsub("\\s+", "_", prefix_label), ".png")
  ggsave(out_name, heatmap_plot, width = width, height = height, dpi = 600, bg = "white")
  print(heatmap_plot)
}

plot_count_keyauthor_presence <- function(subset_df, key_authors, width=30, height=32,prefix_label = "Subset") {
  
  
  tidyplot_theme <- function(base_size = 8) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "black"),
        plot.title = element_text(
          size = base_size * 1.4,        
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.subtitle = element_text(
          size = base_size * 1.25,       
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        axis.title = element_text(
          size = base_size * 1.2,        
          face = "plain"
        ),
        axis.title.x = element_text(
          margin = margin(t = 8), 
          hjust = 0.5
        ),
        axis.title.y = element_text(
          margin = margin(r = 8), 
          hjust = 0.5
        ),
        axis.text = element_text(
          size = base_size,
          color = "black"
        ),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = base_size * 1.2),
        legend.text  = element_text(size = base_size)
      )
  }
  
  num_studies <- nrow(subset_df)
  study_ids <- subset_df$Study.Identifier
  presence_mat <- matrix(0, nrow = num_studies, ncol = num_studies)
  color_mat <- matrix(NA, nrow = num_studies, ncol = num_studies)
  
  rownames(presence_mat) <- study_ids
  colnames(presence_mat) <- study_ids
  rownames(color_mat) <- study_ids
  colnames(color_mat) <- study_ids
  okabe_ito_colors <- c(
    "#0C58CA", "#009E73", "#CC79A7","#FF8C27", "#F0E442")
  if (length(key_authors) > length(okabe_ito_colors)) {
    stop("Too many key authors for available color palette.")
  }
  author_colors <- setNames(okabe_ito_colors[seq_along(key_authors)], key_authors)
  study_authors_list <- mapply(function(auth_str, cons_str) {
    authors <- str_split(auth_str, ",\\s*")[[1]]
    consortium <- str_split(cons_str, ",\\s*")[[1]]
    unique(c(authors, consortium))
  }, subset_df$Authors, subset_df$Consortium, SIMPLIFY = FALSE)
  for (i in seq_len(num_studies)) {
    for (j in i:num_studies) {
      shared <- intersect(study_authors_list[[j]], key_authors)
      presence_mat[i, j] <- length(shared)
      presence_mat[j, i] <- length(shared)
      
      if (length(shared) > 0) {
        shared_colors <- col2rgb(author_colors[shared])
        avg_rgb <- rowMeans(shared_colors)
        hex_color <- rgb(avg_rgb[1]/255, avg_rgb[2]/255, avg_rgb[3]/255)
        color_mat[i, j] <- hex_color
        color_mat[j, i] <- hex_color
      } else {
        color_mat[i, j] <- "#FFFFFF"
        color_mat[j, i] <- "#FFFFFF"
      }
    }
  }
  heatmap_data <- melt(presence_mat)
  colnames(heatmap_data) <- c("Study_A", "Study_B", "Author_Count")
  color_data <- melt(color_mat)
  colnames(color_data) <- c("Study_A", "Study_B", "Fill_Color")
  heatmap_data <- left_join(heatmap_data, color_data, by = c("Study_A", "Study_B"))
  author_label <- paste(key_authors, collapse = ", ")
  file_suffix <- gsub("\\s+", "_", paste(key_authors, collapse = "_"))
  legend_df <- data.frame(
    Author = names(author_colors),
    Color = unname(author_colors)
  )
  heatmap_plot <- ggplot(heatmap_data, aes(x = Study_A, y = Study_B)) +
    geom_tile(aes(fill = Fill_Color), color = "white") +
    scale_fill_identity(guide = "none") +
    geom_text(aes(label = Author_Count), color = tidyplot_colors[1], size = 14) +
    tidyplot_theme() +
    labs(
      title = paste0("Key Author Presence: ", author_label),
      subtitle = "Tile color = blended shared authors | Number = overlap count",
      x = NULL, y = NULL
    ) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(size = 20, hjust = 0.5)
    ) +
    coord_fixed()
  out_name <- paste0("plots/Author_blended_heatmap_", gsub("\\s+", "_", prefix_label), "_", file_suffix, ".png")
  
  library(ggforce)
  
  plot_venn_color_legend <- function(author_colors) {
    stopifnot(length(author_colors) <= 3)
    venn_layout <- data.frame(
      Author = names(author_colors),
      x = c(0, 1, 0.5)[1:length(author_colors)],
      y = c(0, 0, 0.8)[1:length(author_colors)],
      radius = 1
    )
    col_rgb <- lapply(author_colors, function(clr) col2rgb(clr))
    blend_rgb <- function(rgb_list) {
      avg_rgb <- rowMeans(do.call(cbind, rgb_list))
      rgb(avg_rgb[1]/255, avg_rgb[2]/255, avg_rgb[3]/255)
    }
    blends <- list()
    blends[[1]] <- blend_rgb(list(col_rgb[[1]]))
    blends[[2]] <- blend_rgb(list(col_rgb[[2]]))
    if (length(author_colors) == 3) blends[[3]] <- blend_rgb(list(col_rgb[[3]]))
    blends[["12"]] <- blend_rgb(col_rgb[1:2])
    if (length(author_colors) == 3) {
      blends[["13"]] <- blend_rgb(col_rgb[c(1,3)])
      blends[["23"]] <- blend_rgb(col_rgb[c(2,3)])
      blends[["123"]] <- blend_rgb(col_rgb[1:3])
    }
    venn_plot <- ggplot() +
      geom_circle(data = venn_layout, aes(x0 = x, y0 = y, r = radius, fill = Author), 
                  alpha = 0.5, color = "black") +
      scale_fill_manual(values = author_colors) +
      coord_fixed() +
      theme_void() +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 24)
      ) +
      labs(title = " ")
    return(venn_plot)
  }
  author_colors <- setNames(okabe_ito_colors[1:length(key_authors)], key_authors)
  venn_plot <- plot_venn_color_legend(author_colors)
  combined_plot <- venn_plot / heatmap_plot +
    plot_layout(heights = c(0.2, 1))
  ggsave(out_name, combined_plot, width = width, height = height, dpi = 600, bg = 'white')
  print(combined_plot)
}


plot_key_author_presence <- function(subset_df, key_authors, width = 20, height = 20, prefix_label = "Subset") {
  
  tidyplot_theme <- function(base_size = 8) {
    theme_minimal(base_size = base_size) +
      theme(
        text = element_text(family = "Arial", color = "black"),
        plot.title = element_text(
          size = base_size * 1.4,        
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.subtitle = element_text(
          size = base_size * 1.25,       
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        axis.title = element_text(
          size = base_size * 1.2,        
          face = "plain"
        ),
        axis.title.x = element_text(
          margin = margin(t = 8), 
          hjust = 0.5
        ),
        axis.title.y = element_text(
          margin = margin(r = 8), 
          hjust = 0.5
        ),
        axis.text = element_text(
          size = base_size,
          color = "black"
        ),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = base_size * 1.2),
        legend.text  = element_text(size = base_size)
      )
  }
  
  num_studies <- nrow(subset_df)
  presence_mat <- matrix(0, nrow = num_studies, ncol = num_studies)
  rownames(presence_mat) <- subset_df$Study.Identifier
  colnames(presence_mat) <- subset_df$Study.Identifier
  
  for (i in seq_len(num_studies)) {
    for (j in i:num_studies) {
      authors_j    <- str_split(subset_df$Authors[j], ",\\s*")[[1]]
      consortium_j <- str_split(subset_df$Consortium[j], ",\\s*")[[1]]
      combined_j   <- unique(c(authors_j, consortium_j))
      
      is_present <- any(combined_j %in% key_authors)
      val <- ifelse(is_present, 100, 0)
      
      presence_mat[i, j] <- val
      presence_mat[j, i] <- val
    }
  }
  
  heatmap_data <- melt(presence_mat)
  colnames(heatmap_data) <- c("Study_A", "Study_B", "Presence")
  author_label <- paste(key_authors, collapse = ", ")
  file_suffix  <- gsub("\\s+", "_", paste(key_authors, collapse = "_"))
  heatmap_plot <- ggplot(heatmap_data, aes(x = Study_A, y = Study_B, fill = Presence)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = tidyplot_colors[2], name = "Presence (%)") +
    geom_text(aes(label = ifelse(Presence == 100, "✔", "")), 
              color = tidyplot_colors[1], size = 10) +
    tidyplot_theme() +
    labs(
      title    = paste0("Key Author Presence: ", author_label),
      subtitle = "16S rRNA Amplicon",
      x = NULL, 
      y = NULL
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(size = 28, hjust = 0.5)
    ) +
    coord_fixed()
  if (!dir.exists("plots")) dir.create("plots")
  out_name <- paste0("plots/Author_presence_heatmap_", gsub("\\s+", "_", prefix_label), "_", file_suffix, ".png")
  ggsave(out_name, heatmap_plot, width = width, height = height, dpi = 600, bg = 'white')
  print(heatmap_plot)
}

# ==== Analysis ================================================================

remove_unwanted_features <- function(features) {
  setdiff(features, c(".outcome", "Shannon diversity"))
}

get_predictions <- function(study_name, study_list, datasets) {
  
  sdata <- study_list[[study_name]]
  
  if ("galaxypredicted_load_df" %in% names(sdata) &&
      "metacardispredicted_load_df" %in% names(sdata)) {
    galaxy_preds <- sdata$galaxypredicted_load_df %>%
      dplyr::select(PredictedLoad) %>%
      dplyr::mutate(Model = "GALAXY")
    metacardis_preds <- sdata$metacardispredicted_load_df %>%
      dplyr::select('PredictedLoad') %>%
      dplyr::mutate(Model = "MetaCardis")
    
    preds <- bind_rows(galaxy_preds, metacardis_preds)
  } else if ("predicted_load_df" %in% names(sdata)) {
    preds <- sdata$predicted_load_df %>%
      dplyr::select('PredictedLoad') %>%
      dplyr::mutate(Model = "Vandeputte")
  } 
  preds %>%
    dplyr::mutate(Study = study_name)
}

get_sample_size <- function(study_name, study_list, datasets) {
  
  sdata <- study_list[[study_name]]
  
  if ("galaxypredicted_load_df" %in% names(sdata) &&
      "metacardispredicted_load_df" %in% names(sdata)) {
    galaxy_preds <- sdata$galaxypredicted_load_df %>%
      dplyr::select(PredictedLoad) %>%
      dplyr::mutate(Model = "GALAXY")
    metacardis_preds <- sdata$metacardispredicted_load_df %>%
      dplyr::select('PredictedLoad') %>%
      dplyr::mutate(Model = "MetaCardis")
    preds <- bind_rows(galaxy_preds, metacardis_preds)
  } else if ("predicted_load_df" %in% names(sdata)) {
    preds <- sdata$predicted_load_df %>%
      dplyr::select('PredictedLoad') %>%
      dplyr::mutate(Model = "Vandeputte")
  } 
  preds %>%
    dplyr::mutate(Study = study_name)
}


extract_log_stats <- function(study_name, study_list, datasets) {
  is_shotgun <- all(datasets[datasets$dataset_name == study_name, "sequencingtype"] == "shotgun metagenomics")
  sdata <- study_list[[study_name]]
  if (is_shotgun) {
    if (!("measured_load_df" %in% names(sdata)) ||
        !("galaxypredicted_load_df" %in% names(sdata)) ||
        !("metacardispredicted_load_df" %in% names(sdata))) {
      return(
        tibble(
          Study = study_name,
          log_measured_mean       = NA_real_,
          log_measured_sd         = NA_real_,
          N_measured              = NA_integer_,
          
          log_galaxy_mean         = NA_real_,
          log_galaxy_sd           = NA_real_,
          N_galaxy                = NA_integer_,
          
          log_metacardis_mean     = NA_real_,
          log_metacardis_sd       = NA_real_,
          N_metacardis            = NA_integer_
        )
      )
    }
    
    measured_df       <- sdata$measured_load_df
    galaxy_df         <- sdata$galaxypredicted_load_df
    metacardis_df     <- sdata$metacardispredicted_load_df
    
    if (!("MeasuredLoad" %in% colnames(measured_df))) {
      stop("For study: ", study_name, ", measured_load_df must have a 'MeasuredLoad' column.")
    }
    if (!("PredictedLoad" %in% colnames(galaxy_df))) {
      stop("For study: ", study_name, ", galaxypredicted_load_df must have a 'PredictedLoad' column.")
    }
    if (!("PredictedLoad" %in% colnames(metacardis_df))) {
      stop("For study: ", study_name, ", metacardispredicted_load_df must have a 'PredictedLoad' column.")
    }
    
    measured_vals    <- measured_df$MeasuredLoad
    galaxy_vals      <- galaxy_df$PredictedLoad
    metacardis_vals  <- metacardis_df$PredictedLoad
    
    tibble(
      Study               = study_name,
      log_measured_mean   = mean(measured_vals, na.rm = TRUE),
      log_measured_sd     = sd(measured_vals, na.rm = TRUE),
      N_measured          = sum(!is.na(measured_vals)),
      
      log_galaxy_mean     = mean(galaxy_vals, na.rm = TRUE),
      log_galaxy_sd       = sd(galaxy_vals, na.rm = TRUE),
      N_galaxy            = sum(!is.na(galaxy_vals)),
      
      log_metacardis_mean = mean(metacardis_vals, na.rm = TRUE),
      log_metacardis_sd   = sd(metacardis_vals, na.rm = TRUE),
      N_metacardis        = sum(!is.na(metacardis_vals))
    )
    
  } else {
    if (!("measured_load_df" %in% names(sdata)) || !("predicted_load_df" %in% names(sdata))) {
      return(
        tibble(
          Study = study_name,
          log_measured_mean  = NA_real_,
          log_measured_sd    = NA_real_,
          log_predicted_mean = NA_real_,
          log_predicted_sd   = NA_real_,
          N_measured         = NA_integer_,
          N_predicted        = NA_integer_
        )
      )
    }
    
    measured_df  <- sdata$measured_load_df
    predicted_df <- sdata$predicted_load_df
    if (!("MeasuredLoad" %in% colnames(measured_df))) {
      stop("For study: ", study_name, ", measured_load_df must have a 'MeasuredLoad' column.")
    }
    if (!("PredictedLoad" %in% colnames(predicted_df))) {
      stop("For study: ", study_name, ", predicted_load_df must have a 'PredictedLoad' column.")
    }
    
    measured_vals  <- measured_df$MeasuredLoad
    predicted_vals <- predicted_df$PredictedLoad
    
    tibble(
      Study              = study_name,
      log_measured_mean  = mean(measured_vals, na.rm = TRUE),
      log_measured_sd    = sd(measured_vals, na.rm = TRUE),
      N_measured         = sum(!is.na(measured_vals)),
      
      log_predicted_mean = mean(predicted_vals, na.rm = TRUE),
      log_predicted_sd   = sd(predicted_vals, na.rm = TRUE),
      N_predicted        = sum(!is.na(predicted_vals))
    )
  }
}


# AA. Compute R^2
compute_r2 <- function(study_data, study_name, datasets) {
  is_shotgun <- all(datasets[datasets$dataset_name == study_name, "sequencingtype"] == "shotgun metagenomics")
  
  if (is_shotgun) {
    if (!all(c("measured_load_df", "galaxypredicted_load_df", "metacardispredicted_load_df") %in% names(study_data))) {
      warning("Missing one or more required data frames for shotgun study. Returning NA.")
      return(list(r2_galaxy = NA, r_galaxy = NA, p_galaxy = NA, 
                  r2_metacardis = NA, r_metacardis = NA, p_metacardis = NA))
    }
    
    measured_df   <- study_data$measured_load_df
    galaxy_df     <- study_data$galaxypredicted_load_df
    metacardis_df <- study_data$metacardispredicted_load_df
    
    compute_metrics <- function(measured, predicted) {
      if (nrow(measured) == 0 || nrow(predicted) == 0) {
        return(list(r2 = NA, r = NA, p = NA))
      }
      
      merged <- merge(measured, predicted, by = "SampleID", all = FALSE)
      if (!all(c("MeasuredLoad", "PredictedLoad") %in% colnames(merged))) {
        warning("Missing MeasuredLoad or PredictedLoad columns.")
        return(list(r2 = NA, r = NA, p = NA))
      }
      
      r2_value <- cor(merged$MeasuredLoad, merged$PredictedLoad, use = "complete.obs")^2
      cor_test <- cor.test(merged$MeasuredLoad, merged$PredictedLoad, method = "pearson")
      r_value <- round(cor_test$estimate, 2)
      p_value <- format(cor_test$p.value, scientific = TRUE, digits = 2)
      
      return(list(r2 = r2_value, r = r_value, p = p_value))
    }
    
    galaxy_metrics <- compute_metrics(measured_df, galaxy_df)
    metacardis_metrics <- compute_metrics(measured_df, metacardis_df)
    
    return(list(
      r2_galaxy = galaxy_metrics$r2, r_galaxy = galaxy_metrics$r, p_galaxy = galaxy_metrics$p,
      r2_metacardis = metacardis_metrics$r2, r_metacardis = metacardis_metrics$r, p_metacardis = metacardis_metrics$p
    ))
    
  } else {
    if (!all(c("measured_load_df", "predicted_load_df") %in% names(study_data))) {
      warning("Missing measured_load_df or predicted_load_df for non-shotgun study. Returning NA.")
      return(list(r2 = NA, r = NA, p = NA))
    }
    
    if (nrow(study_data$measured_load_df) == 0 || nrow(study_data$predicted_load_df) == 0) {
      warning("Empty dataframes in study. Skipping R^2 computation.")
      return(list(r2 = NA, r = NA, p = NA))
    }
    
    merged_df <- merge(
      study_data$measured_load_df,
      study_data$predicted_load_df,
      by = "SampleID",   
      all = FALSE  
    )
    
    if (!all(c("MeasuredLoad", "PredictedLoad") %in% colnames(merged_df))) {
      warning("MeasuredLoad or PredictedLoad columns missing. Returning NA.")
      return(list(r2 = NA, r = NA, p = NA))
    }
    
    r2_value <- cor(merged_df$MeasuredLoad, merged_df$PredictedLoad, use = "complete.obs")^2
    cor_test <- cor.test(merged_df$MeasuredLoad, merged_df$PredictedLoad, method = "pearson")
    r_value <- round(cor_test$estimate, 2)
    p_value <- format(cor_test$p.value, scientific = TRUE, digits = 2)
    
    return(list(r2 = r2_value, r = r_value, p = p_value))
  }
}


# A. Alpha Diversity
compute_alpha_diversity <- function(relabun_df) {
  # Observed = number of taxa with abundance > 0
  Observed <- rowSums(relabun_df > 0)
  
  # Shannon, Simpson from vegan
  Shannon <- diversity(relabun_df, index = "shannon")
  Simpson <- diversity(relabun_df, index = "simpson")
  
  # Evenness (Pielou's) = Shannon / log(Observed)
  Evenness <- Shannon / log(Observed)
  
  alpha_df <- data.frame(
    Observed = Observed,
    Shannon  = Shannon,
    Simpson  = Simpson,
    Evenness = Evenness*100,
    stringsAsFactors = FALSE
  )
  return(alpha_df)
}

# B. Sparsity of each sample
compute_sparsity_per_sample <- function(relabun_df) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a dataframe or matrix.")
  }
  
  sparsity <- apply(relabun_df, 1, function(row) {
    mean(row == 0) 
  })
  return(sparsity)
}

# C. Amount of rare taxa per sample below threshold (same threshold used in MLP study)
evaluate_rare_taxa <- function(relabun_df, threshold = 0.1) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a dataframe or matrix.")
  }
  
  rare_taxa_count <- apply(relabun_df, 1, function(row) {
    sum(row < threshold & row > 0)  
  })
  
  return(rare_taxa_count)
}

# rare taxa compared to reference taxa
evaluate_rare_taxa_reference <- function(relabun_df, referencetaxa, threshold = 0.1) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a dataframe or matrix.")
  }
  common_taxa <- intersect(colnames(relabun_df), referencetaxa)
  
  if (length(common_taxa) == 0) {
    stop("No overlapping taxa between input data and reference taxa.")
  }
  relabun_subset <- relabun_df[, common_taxa, drop = FALSE]
  rare_taxa_percent <- apply(relabun_subset, 1, function(row) {
    present_taxa <- sum(row > 0)
    rare_taxa <- sum(row > threshold & row > 0)
    if (present_taxa == 0) {
      return(NA_real_) 
    } else {
      return((rare_taxa / present_taxa) * 100)
    }
  })
  return(rare_taxa_percent)
}


# D. Obtain presence absence, but we dont really use this feature
presence_absence <- function(relabun_df) {
  pa_mat <- (relabun_df > 0) * 1
  colnames(pa_mat) <- paste0("PA_", colnames(pa_mat))
  return(pa_mat)
}

# E. Compute the weighted proportion of overlap per sample. Essentially the intersection between dataset taxa and the referencetaxa are taken and then of the present taxa for each sample calculate the sum of the present taxa and divide by the number of total taxa in the dataset
compute_proportion_overlap_per_sample <- function(relabun_df, referencetaxa) {
  overlapping_taxa <- intersect(colnames(relabun_df), referencetaxa)
  prop_vec <- (rowSums(relabun_df[, overlapping_taxa, drop = FALSE] > 0, na.rm = TRUE) /
    length(relabun_df))*100
  
  return(prop_vec)
}

# F. Compute the number of taxa in the dataset that overlap with the reference dataset and divide by the total taxa in the reference dataset
compute_proportion_overlap_total <- function(relabun_df, referencetaxa) {
  
  overlapping_taxa <- intersect(colnames(relabun_df), referencetaxa)
  total_overlap <- (length(overlapping_taxa) / length(referencetaxa))*100

  return(rep(total_overlap, nrow(relabun_df)))
}

# G. Compute the weighted proportion of non overlap per sample. Essentially the different between dataset taxa and the referencetaxa are taken and then of the present taxa for each sample calculate the sum of the present taxa (that dont overlap with reference) and divide by the number of total taxa in the dataset
compute_proportion_nonoverlap_per_sample <- function(relabun_df, referencetaxa) {
  
  nonoverlapping_taxa <- setdiff(colnames(relabun_df), referencetaxa)
  if (length(nonoverlapping_taxa) == 0) {
    return(rep(0, nrow(relabun_df)))
  }
  prop_vec <- (rowSums(relabun_df[, nonoverlapping_taxa, drop = FALSE] > 0, na.rm = TRUE) /
    length(nonoverlapping_taxa))*100
  
  return(prop_vec)
}

# H. Compute the number of taxa in the dataset that do not overlap with the reference dataset and divide by the total taxa in the reference dataset
compute_proportion_nonoverlap_total <- function(relabun_df, referencetaxa) {
  
  nonoverlapping_training_taxa <- setdiff(referencetaxa, colnames(relabun_df))
  total_nonoverlap <- (length(nonoverlapping_training_taxa) / length(referencetaxa))*100
  
  return(rep(total_nonoverlap, nrow(relabun_df)))
}

# I. Compute the number of overlapping taxa in the dataset with the reference dataset and this time, divide by total taxa in reference dataset
compute_proportion_overlap_totalreference <- function(relabun_df, referencetaxa) {
  
  overlapping_taxa <- intersect(referencetaxa, colnames(relabun_df))
  total_overlapping_taxa <- (length(overlapping_taxa) / length(referencetaxa))*100
  
  return(rep(total_overlapping_taxa, nrow(relabun_df)))
}

# J. Compute the beta diversity - bray curtis between the dataset and reference dataset.
compute_beta_diversity <- function(dataset, reference_taxa, method = "bray") {
  if (!is.data.frame(dataset) && !is.matrix(dataset)) {
    stop("`dataset` must be a data frame or matrix.")
  }
  if (!is.data.frame(reference_taxa) && !is.matrix(reference_taxa)) {
    stop("`reference_taxa` must be a data frame or matrix.")
  }
  if (is.null(colnames(dataset))) {
    stop("`dataset` must have column names.")
  }
  if (is.null(colnames(reference_taxa))) {
    stop("`reference_taxa` must have column names.")
  }
  all_taxa <- union(colnames(dataset), colnames(reference_taxa))
  if (length(all_taxa) == 0) {
    stop("No overlapping taxa between `dataset` and `reference_taxa`.")
  }
  dataset_aligned <- matrix(0, nrow = nrow(dataset), ncol = length(all_taxa),
                            dimnames = list(rownames(dataset), all_taxa))
  reference_aligned <- matrix(0, nrow = nrow(reference_taxa), ncol = length(all_taxa),
                              dimnames = list(rownames(reference_taxa), all_taxa))
  dataset_aligned[, colnames(dataset)] <- as.matrix(dataset)
  reference_aligned[, colnames(reference_taxa)] <- as.matrix(reference_taxa)
  dataset_aligned[is.na(dataset_aligned)] <- 0
  reference_aligned[is.na(reference_aligned)] <- 0
  dataset_sum <- colSums(dataset_aligned)
  reference_sum <- colSums(reference_aligned)
  combined <- rbind(dataset_sum, reference_sum)
  combined[is.na(combined)] <- 0
  beta_div <- vegdist(combined, method = method)
  if (length(beta_div) == 1) {
    beta_div_vector <- rep(as.numeric(beta_div), nrow(dataset))
  } else {
    stop("Beta diversity computation resulted in multiple values. Ensure your input data is correct.")
  }
  
  return(beta_div_vector)
}

numberoftaxausedformodeltesting <- function(dataset, reference_taxa) {
  if (!is.data.frame(dataset) && !is.matrix(dataset)) {
    stop("`dataset` must be a data frame or matrix.")
  }
  if (is.null(colnames(dataset))) {
    stop("`dataset` must have column names.")
  }
  all_taxa <- intersect(colnames(dataset), reference_taxa)
  return(length(all_taxa))
}

numberoftaxaaddedformodel <- function(dataset, reference_taxa) {
  if (!is.data.frame(dataset) && !is.matrix(dataset)) {
    stop("`dataset` must be a data frame or matrix.")
  }
  if (is.null(colnames(dataset))) {
    stop("`dataset` must have column names.")
  }
  all_taxa <- setdiff(reference_taxa, colnames(dataset))
  return(length(all_taxa))
}

build_features_for_study <- function(study_name,
                                     relabun_df,        
                                     measured_load_df,  
                                     predicted_load_df, 
                                     galaxypredicted_load_df = NULL, 
                                     metacardispredicted_load_df = NULL,  
                                     reference_taxa,
                                     reference_name,
                                     datasets) {
  
  print(paste("Processing:", study_name))
  
  if ("SampleID" %in% colnames(relabun_df)) {
    rownames(relabun_df) <- relabun_df$SampleID
    relabun_df$SampleID  <- NULL
  }
  is_shotgun <- any(datasets[datasets$dataset_name == study_name, "sequencingtype"] == "shotgun metagenomics")
  if (is_shotgun) {
    if (reference_taxa == "motus25" & reference_name == 'GALAXY (MLP Training data)') {
      reference_taxa = galaxymotus25features
    } else if (reference_taxa == "motus25" & reference_name == 'MetaCardis (MLP Training data)') {
      reference_taxa = metacardismotus25features
    } else if (reference_taxa == "motus3" & reference_name == 'GALAXY (MLP Training data)') {
      reference_taxa = galaxymotus3features
      reference_dataset = NA
    } else if (reference_taxa == "motus3" & reference_name == 'MetaCardis (MLP Training data)') {
      reference_taxa = metacardismotus3features
      reference_dataset = NA
    } else if (reference_taxa == "metaphlan3" & reference_name == 'GALAXY (MLP Training data)') {
      reference_taxa = galaxymetaphlan3features
      reference_dataset = NA
    } else if (reference_taxa == "metaphlan3" & reference_name == 'MetaCardis (MLP Training data)') {
      reference_taxa = metacardismetaphlan3features
      reference_dataset = NA
    } else if (reference_taxa == "metaphlan4" & reference_name == 'GALAXY (MLP Training data)') {
      reference_taxa = galaxymetaphlan4features
      reference_dataset = NA
    } else if (reference_taxa == "metaphlan4" & reference_name == 'MetaCardis (MLP Training data)') {
      reference_taxa = metacardismetaphlan4features
      reference_dataset = NA
    } else {
      stop(paste("dataset:",study_name,"does not have the classifier column filled out in datasets.csv or the type of classifier is incorrectly spelled: motus25, motus3, metaphlan3, metaphlan4"))
    }
  }
  
  if (reference_name %in% names(studies)) {
    if (!is.null(studies[[reference_name]]$relabun_df)) {
      possible_dataset <- tryCatch(
        {
          rel_df <- studies[[reference_name]]$relabun_df
          common_taxa <- intersect(reference_taxa, colnames(rel_df))
          
          if (length(common_taxa) == 0) {
            message(
              "No requested reference_taxa found in '", reference_name,
              "'. Setting beta diversity = NA."
            )
            NA
          } else {
            rel_df[, common_taxa, drop = FALSE]
          }
        },
        error = function(e) {
          message("Error extracting reference dataset: ", e$message)
          return(NA)
        }
      )
      if (is.data.frame(possible_dataset)) {
        reference_dataset <- possible_dataset
      }
    } else {
      message("No 'relabun_df' found in '", reference_name, "' - setting beta diversity = NA.")
    }
  } else {
    message("Reference name '", reference_name, "' not found in studies - setting beta diversity = NA.")
  }
  
  if (is.null(rownames(relabun_df))) {
    stop(paste("relabun_df for", study_name, 
               "must have rownames = Sample IDs or a 'SampleID' column."))
  }
  relabun_df[] <- lapply(relabun_df, function(col) {
    if (is.character(col) || is.factor(col)) {
      as.numeric(as.character(col))
    } else {
      col
    }
  })
  relabun_df <- relabun_df[, colSums(is.na(relabun_df)) < nrow(relabun_df), drop = FALSE]
  sample_ids <- rownames(relabun_df)
  
  # Proportions of overlap
  prop_overlap_vec <- compute_proportion_overlap_per_sample(relabun_df, reference_taxa)
  prop_overlap_total_vec <- compute_proportion_overlap_total(relabun_df, reference_taxa)
  prop_nonoverlap_vec <- compute_proportion_nonoverlap_per_sample(relabun_df, reference_taxa)
  prop_nonoverlap_total_vec <- compute_proportion_nonoverlap_total(relabun_df, reference_taxa)
  prop_overlap_complete <- compute_proportion_overlap_totalreference(relabun_df, reference_taxa)

  # Taxa used and added into dataset for model
  numberoftaxausedformodeltesting = numberoftaxausedformodeltesting(relabun_df, reference_taxa)
  numberoftaxaaddedformodel = numberoftaxaaddedformodel(relabun_df, reference_taxa)
  
  # Compute alpha diversity
  alpha_df <- compute_alpha_diversity(relabun_df)
  
  # Sparsity per sample
  sparsitysamples <- compute_sparsity_per_sample(relabun_df)

  # Proportion of rare taxa
  raretaxa <- evaluate_rare_taxa(relabun_df)
  raretaxamlp <- evaluate_rare_taxa_reference(relabun_df, reference_taxa)
  
  # Compute beta diversity
  if (is.data.frame(reference_dataset)) {
    betadiversity <- compute_beta_diversity(relabun_df, reference_dataset)
  } else {
    betadiversity <- NA
  }
  relabun_df[setdiff(reference_taxa, colnames(relabun_df))] <- 0
  relabun_df <- relabun_df[, reference_taxa, drop = FALSE]
  compute_alpha_diversity_after_model <- function(relabun_df) {
    # Observed = number of taxa with abundance > 0
    `Observed After Model Predictions` <- rowSums(relabun_df > 0)
    
    # Shannon, Simpson from vegan
    `Shannon After Model Predictions` <- diversity(relabun_df, index = "shannon")
    `Simpson After Model Predictions` <- diversity(relabun_df, index = "simpson")
    
    # Evenness (Pielou's) = Shannon / log(Observed)
    `Evenness After Model Predictions` <- `Shannon After Model Predictions` / log(`Observed After Model Predictions`)
    
    alpha_df <- data.frame(
      `Observed of Taxa Used in MLP` = `Observed After Model Predictions`,
      `Shannon of Taxa Used in MLP`  = `Shannon After Model Predictions`,
      `Simpson of Taxa Used in MLP`  = `Simpson After Model Predictions`,
      `Evenness of Taxa Used in MLP` = `Evenness After Model Predictions`,
      stringsAsFactors = FALSE
    )
    return(alpha_df)
  }
  #alpha_dftrainingtaxaonly <- compute_alpha_diversity_after_model(relabun_df)
  raretaxatrainingtaxaonly <- evaluate_rare_taxa(relabun_df)
  sparsitysamplestrainingtaxaonly <- compute_sparsity_per_sample(relabun_df)
  prop_overlap_vectrainingtaxaonly <- compute_proportion_overlap_per_sample(relabun_df, reference_taxa)
  prop_nonoverlap_vectrainingtaxaonly <- compute_proportion_nonoverlap_per_sample(relabun_df, reference_taxa)
  
  feature_df <- data.frame(
    SampleID            = sample_ids,
    alpha_df,           
    PropOverlapRef      = prop_overlap_vec, 
    PropOverlapReftrainingtaxaonly  = prop_overlap_vectrainingtaxaonly, 
    PropOverlapRefTotal = prop_overlap_total_vec, 
    PropNonoverlapRef   = prop_nonoverlap_vec, 
    PropNonoverlapReftrainingtaxaonly  = prop_nonoverlap_vectrainingtaxaonly, 
    PropNonoverlapRefTotal  = prop_nonoverlap_total_vec, 
    PropOverlapComplete = prop_overlap_complete,
    datasetbetadiversity = betadiversity,
    sparsityofsample = sparsitysamples,
    raretaxaofsample = raretaxa,
    raretaxamlp = raretaxamlp,
    #alpha_dftrainingtaxaonly,
    raretaxatrainingtaxaonly = raretaxatrainingtaxaonly,
    sparsitysamplestrainingtaxaonly = sparsitysamplestrainingtaxaonly,
    numberoftaxausedformodeltesting = numberoftaxausedformodeltesting,
    numberoftaxaaddedformodel = numberoftaxaaddedformodel,
    stringsAsFactors    = FALSE
  )
  feature_df$Study <- study_name
  if (!is.null(measured_load_df) && nrow(measured_load_df) > 0) {
    feature_df <- feature_df %>%
      left_join(measured_load_df, by = "SampleID")
  } else {
    feature_df$MeasuredLoad <- NA
  }
  
  if (is_shotgun) {
    if (!is.null(galaxypredicted_load_df) && nrow(galaxypredicted_load_df) > 0) {
      feature_df <- feature_df %>%
        left_join(galaxypredicted_load_df, by = "SampleID", suffix = c("", "_GALAXY"))
      colnames(feature_df)[colnames(feature_df) == "PredictedLoad"] <- "GALAXYPredictedLoad"
    } else {
      feature_df$GALAXYPredictedLoad <- NA
    }
    if (!is.null(metacardispredicted_load_df) && nrow(metacardispredicted_load_df) > 0) {
      feature_df <- feature_df %>%
        left_join(metacardispredicted_load_df, by = "SampleID", suffix = c("", "_MetaCardis"))
      colnames(feature_df)[colnames(feature_df) == "PredictedLoad"] <- "MetaCardisPredictedLoad"
    } else {
      feature_df$MetaCardisPredictedLoad <- NA
    }
    if (reference_name == "GALAXY (MLP Training data)") {
      feature_df <- feature_df %>%
        mutate(PredictedLoad = GALAXYPredictedLoad)
    } else if (reference_name == "MetaCardis (MLP Training data)") {
      feature_df <- feature_df %>%
        mutate(PredictedLoad = MetaCardisPredictedLoad)
    }
  } else {
    if (!is.null(predicted_load_df)) {
      feature_df <- feature_df %>%
        left_join(predicted_load_df, by = "SampleID")
    } else {
      feature_df$PredictedLoad <- NA
    }
  } 
  return(feature_df)
}


build_and_process_features <- function(datasets, studies) {
  dataset_names <- names(studies)
  all_features <- map_dfr(dataset_names, function(dn) {
    stype <- datasets %>%
      filter(dataset_name == dn) %>%
      pull(sequencingtype) %>%
      unique()
    sdata <- studies[[dn]]
    if (length(stype) != 1) {
      stop("Dataset '", dn, "' has multiple or zero sequencing types in 'datasets'.")
    }
    if (stype == "16S rRNA") {
      ref_name <- "Vandeputte2021 (MLP Training data)"
      out_df <- build_features_for_study(
        study_name                     = dn,
        relabun_df                     = sdata$relabun_df,
        measured_load_df               = sdata$measured_load_df,
        predicted_load_df              = sdata$predicted_load_df,
        galaxypredicted_load_df        = NULL,
        metacardispredicted_load_df    = NULL,
        reference_taxa                 = rdpclassifierfeatures,
        reference_name                 = ref_name,
        datasets                       = datasets 
      )
      out_df <- out_df %>%
        mutate(Reference = ref_name)
      return(out_df)
    } else if (stype == "shotgun metagenomics") {
      ref_names <- c("GALAXY (MLP Training data)", "MetaCardis (MLP Training data)")
      out_list <- lapply(ref_names, function(ref_n) {
        build_features_for_study(
          study_name                     = dn,
          relabun_df                     = sdata$relabun_df,
          measured_load_df               = sdata$measured_load_df,
          predicted_load_df              = NULL,  
          galaxypredicted_load_df        = sdata$galaxypredicted_load_df,
          metacardispredicted_load_df    = sdata$metacardispredicted_load_df,
          reference_taxa                 = unique(datasets[datasets$dataset_name == dn & datasets$source_type == "Predicted", "classifier"])[1],
          reference_name                 = ref_n,
          datasets                       = datasets
        ) %>%
          mutate(Reference = ref_n)
      })
      
      out_df <- bind_rows(out_list)
      return(out_df)
    } else {
      message("Skipping dataset ", dn, " due to unknown sequencingtype: ", stype)
      return(NULL)
    }
  })
  cols_to_drop <- c("Subject_ID", "Read.number", "Sample.Coverage",
                    "replicate", "time..hours.", "condition", "medium", "index",
                    "time (hours)","Read number","Sample Coverage","Var.1",
                    "Sample_ID","Avg_Ct","Std_Ct","Std_Calibrated_rConc",
                    "cohort","logCopyNumber")
  all_features <- all_features %>%
    dplyr::select(-any_of(cols_to_drop))
  all_features <- all_features %>%
    dplyr::rename(
      'Beta Diversity: Study vs Training Study' = datasetbetadiversity,
      'Shared Taxa (Study ∩ Training Study)' = PropOverlapComplete,
      'Model Taxa (MLP Features ∉ Study)' = PropNonoverlapRefTotal,
      'Taxa >0 (Sample ∉ Training Study)' = PropNonoverlapRef,
      'Taxa >0 (Sample ∉ MLP Features)' = PropNonoverlapReftrainingtaxaonly,
      'Taxa >0 (Sample ∩ Training Study)' = PropOverlapRef,
      'Taxa >0 (Sample ∩ MLP Features)' = PropOverlapReftrainingtaxaonly,
      'Sample Taxa Sparsity' = sparsityofsample,
      'Sample Rare Taxa' = raretaxaofsample,
      'Sample ∩ MLP Feature Abundant Taxa' = raretaxamlp,
      'Sample ∩ MLP Feature Sparsity' = sparsitysamplestrainingtaxaonly,
      '# Taxa ∩ MLP Features' = numberoftaxausedformodeltesting,
      '# Taxa ∉ MLP Features' = numberoftaxaaddedformodel
    )
  feature_cols <- c(
    "Observed", 
    "Shannon", 
    "Simpson", 
    "Evenness",
    "Shared Taxa (Study ∩ Training Study)",
    "Shared Taxa (Study ∩ MLP Features)",
    "Taxa >0 (Sample ∩ Training Study)",
    "Taxa >0 (Sample ∩ MLP Features)",
    "Taxa >0 (Sample ∉ Training Study)",
    "Taxa >0 (Sample ∉ MLP Features)",
    "# Taxa ∩ MLP Features",
    "# Taxa ∉ MLP Features",
    "Beta Diversity: Study vs Training Study",
    "Sample Taxa Sparsity",
    "Sample Rare Taxa",
    "Sample ∩ MLP Feature Abundant Taxa",
    "Sample ∩ MLP Feature Sparsity"
  )
  long_df <- all_features %>%
    pivot_longer(
      cols      = any_of(feature_cols),
      names_to  = "Feature",
      values_to = "FeatureValue"
    )
  cor_results <- long_df %>%
    mutate(Residual = PredictedLoad - MeasuredLoad) %>%  
    group_by(Study, Feature, Reference) %>%
    summarize(
      FeatureVariance = var(FeatureValue, na.rm = TRUE),
      FeatureMean     = ifelse(FeatureVariance == 0, mean(FeatureValue, na.rm = TRUE), NA_real_),
      Pearson_r_measured = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(MeasuredLoad)) < 2,
        NA,
        cor(FeatureValue, MeasuredLoad, method = "pearson", use = "complete.obs")
      ),
      Spearman_r_measured = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(MeasuredLoad)) < 2,
        NA,
        cor(FeatureValue, MeasuredLoad, method = "spearman", use = "complete.obs")
      ),
      N_measured = sum(!is.na(FeatureValue) & !is.na(MeasuredLoad)),
      Pearson_r_predicted = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(PredictedLoad)) < 2,
        NA,
        cor(FeatureValue, PredictedLoad, method = "pearson", use = "complete.obs")
      ),
      Spearman_r_predicted = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(PredictedLoad)) < 2,
        NA,
        cor(FeatureValue, PredictedLoad, method = "spearman", use = "complete.obs")
      ),
      N_predicted = sum(!is.na(FeatureValue) & !is.na(PredictedLoad)),
      Pearson_r_residual = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(Residual)) < 2,
        NA,
        cor(FeatureValue, Residual, method = "pearson", use = "complete.obs")
      ),
      Spearman_r_residual = ifelse(
        FeatureVariance == 0 || sum(!is.na(FeatureValue) & !is.na(Residual)) < 2,
        NA,
        cor(FeatureValue, Residual, method = "spearman", use = "complete.obs")
      ),
      N_residual = sum(!is.na(FeatureValue) & !is.na(Residual)),
      
      .groups = "drop"
    )
  
  list(
    all_features = all_features,
    long_df      = long_df,
    cor_results  = cor_results
  )
}

# <pre><code class="markdown"> 
#   ### Feature Summary Table | **Feature Name** | **Description** | **Function** |
# |---------------------------------------------|---------------------------------------------------------------------------------------------------|---------------------------------------------------| 
# | **Observed** | Number of taxa with nonzero abundance per sample | `compute_alpha_diversity()` | 
# | **Shannon** | Shannon diversity index per sample | `compute_alpha_diversity()` | 
# | **Simpson** | Simpson diversity index per sample | `compute_alpha_diversity()` | 
# | **Evenness** | Pielou’s evenness (Shannon / log(Observed)) | `compute_alpha_diversity()` | 
# | **Observed of Taxa Used in MLP** | Observed diversity restricted to MLP model features | `compute_alpha_diversity_after_model()` | 
# | **Shannon of Taxa Used in MLP** | Shannon index over MLP model features | `compute_alpha_diversity_after_model()` | 
# | **Simpson of Taxa Used in MLP** | Simpson index over MLP model features | `compute_alpha_diversity_after_model()` | 
# | **Evenness of Taxa Used in MLP** | Evenness over MLP model features | `compute_alpha_diversity_after_model()` | 
# | **Sample Taxa Sparsity** | Fraction of zero values across all taxa per sample | `compute_sparsity_per_sample()` | 
# | **Sample ∩ MLP Feature Sparsity** | Sparsity computed only across taxa used by the MLP model | `compute_sparsity_per_sample()` | 
# | **Sample Rare Taxa** | Count of taxa < 0.1 abundance (but > 0) per sample | `evaluate_rare_taxa()` | 
# | **Sample ∩ MLP Feature Rare Taxa** | Count of rare taxa < 0.1 within the MLP model feature set | `evaluate_rare_taxa()` | 
# | **Taxa >0 (Sample ∩ Training Study)** | Proportion of taxa present in the sample and also present in the **training study** | `compute_proportion_overlap_per_sample()` | 
# | **Taxa >0 (Sample ∩ MLP Features)** | Proportion of taxa present in the sample and also used in the **MLP model** | `compute_proportion_overlap_per_sample()` | 
# | **Taxa >0 (Sample ∉ Training Study)** | Proportion of present taxa in the sample **not found** in the training study | `compute_proportion_nonoverlap_per_sample()` | 
# | **Taxa >0 (Sample ∉ MLP Features)** | Proportion of present taxa in the sample **not included** in the MLP model | `compute_proportion_nonoverlap_per_sample()` | 
# | **Shared Taxa (Study ∩ Training Study)** | Proportion of taxa in study overlapping with training study taxa (global, not per sample) | `compute_proportion_overlap_total()` | 
# | **Shared Taxa (Study ∩ MLP Features)** | Proportion of MLP model features found in the study (global) | `compute_proportion_overlap_totalreference()` | 
# | **Model Taxa (MLP Features ∉ Study)** | Proportion of MLP features **not found** in the study (global) | `compute_proportion_nonoverlap_total()` | 
# | **# Taxa ∩ MLP Features** | Count of model taxa present in the study | `numberoftaxausedformodeltesting()` | 
# | **# Taxa ∉ MLP Features** | Count of model taxa **missing** from the study | `numberoftaxaaddedformodel()` | 
# | **Beta Diversity: Study vs Training Study** | Bray-Curtis beta diversity between aggregated sample and training reference | `compute_beta_diversity()` | 
# </code></pre>

extracttrainingfeatures <- function(studies, datasets) {
  results <- list()  
  for (dataset_name in names(studies)) {
    dataset <- studies[[dataset_name]]
    if (!"relabun_df" %in% names(dataset)) {
      warning(paste("Skipping", dataset_name, "- no relabun_df found."))
      next
    }
    relabun_df <- dataset$relabun_df
    if ("SampleID" %in% colnames(relabun_df)) {
      rownames(relabun_df) <- relabun_df$SampleID
      relabun_df$SampleID <- NULL
    }
    shannondiversity = diversity(relabun_df)
    dataset_info <- datasets[datasets$dataset_name == dataset_name & datasets$source_type == "Predicted", ]
    if (nrow(dataset_info) == 0) {
      warning(paste("Skipping", dataset_name, "- no 'Predicted' source_type found in datasets."))
      next
    }
    for (i in 1:nrow(dataset_info)) {
      sequencing_type <- dataset_info$sequencingtype[i]
      reference_name <- dataset_info$trainingset[i]
      classifier <- dataset_info$classifier[i]
      if (sequencing_type == "16S rRNA") {
        reference_taxa <- rdpclassifierfeatures
      } else if (sequencing_type == "shotgun metagenomics") {
        reference_mappings <- list(
          "motus25" = list(
            "GALAXY" = galaxymotus25features,
            "MetaCardis" = metacardismotus25features
          ),
          "motus3" = list(
            "GALAXY" = galaxymotus3features,
            "MetaCardis" = metacardismotus3features
          ),
          "metaphlan3" = list(
            "GALAXY" = galaxymetaphlan3features,
            "MetaCardis" = metacardismetaphlan3features
          ),
          "metaphlan4" = list(
            "GALAXY" = galaxymetaphlan4features,
            "MetaCardis" = metacardismetaphlan4features
          )
        )
        if (!classifier %in% names(reference_mappings) || 
            !reference_name %in% names(reference_mappings[[classifier]])) {
          warning(paste("Skipping", dataset_name, "- invalid classifier or reference name."))
          next
        }
        reference_taxa <- reference_mappings[[classifier]][[reference_name]]
      } else {
        warning(paste("Skipping", dataset_name, "- unknown sequencing type."))
        next
      }
      intersected_taxa <- intersect(colnames(relabun_df), reference_taxa)
      filtered_df <- relabun_df[, intersected_taxa, drop = FALSE]
      missing_taxa <- setdiff(reference_taxa, colnames(filtered_df))
      if (length(missing_taxa) > 0) {
        missing_df <- matrix(0, nrow = nrow(filtered_df), ncol = length(missing_taxa))
        colnames(missing_df) <- missing_taxa
        rownames(missing_df) <- rownames(filtered_df)
        filtered_df <- cbind(filtered_df, as.data.frame(missing_df))
      }
      if (sequencing_type == "shotgun metagenomics") {
        filtered_df$`Shannon diversity` <- as.numeric(shannondiversity)
      }
      results[[dataset_name]] <- filtered_df
    }
  }
  
  return(results)
}


flatten_and_separate <- function(results, datasets) {
  list_16S <- list()
  list_shotgun <- list()
  dataset_info <- datasets[datasets$dataset_name %in% names(results), c("dataset_name", "sequencingtype")]
  for (dataset_name in names(results)) {
    df <- results[[dataset_name]]
    sequencing_type <- unique(dataset_info$sequencingtype[dataset_info$dataset_name == dataset_name])
    df$dataset_name <- dataset_name  
    if (sequencing_type == "16S rRNA") {
      list_16S[[dataset_name]] <- df
    } else if (sequencing_type == "shotgun metagenomics") {
      list_shotgun[[dataset_name]] <- df
    }
  }
  df_16S <- if (length(list_16S) > 0) do.call(rbind, list_16S) else data.frame()
  if (length(list_shotgun) > 0) {
    all_features <- unique(unlist(lapply(list_shotgun, colnames)))
    list_shotgun <- lapply(list_shotgun, function(df) {
      missing_cols <- setdiff(all_features, colnames(df))
      if (length(missing_cols) > 0) {
        df[, missing_cols] <- 0 
      }
      return(df[, all_features])  
    })
    df_shotgun <- do.call(rbind, list_shotgun)
  } else {
    df_shotgun <- data.frame()
  }
  return(list("16S" = df_16S, "shotgun" = df_shotgun))
}

###### AUTHOR SIMILARITY ######

# author_similarity <- function(authors_A, authors_B) {
#   authors_A <- str_split(authors_A, ",\\s*")[[1]]  
#   authors_B <- str_split(authors_B, ",\\s*")[[1]] 
#   shared_authors <- length(intersect(authors_A, authors_B))
#   total_unique_authors <- length(unique(c(authors_A, authors_B))) 
#   similarity_percent <- ifelse(total_unique_authors == 0, 0, (shared_authors / total_unique_authors) * 100)
#   return(list(similarity_percent, shared_authors))
# }

corresponding_author_similarity <- function(authors_B) {
  authors_B_list <- str_split(authors_B, ",\\s*")[[1]] 
  is_present <- corresponding_author %in% authors_B_list
  return(ifelse(is_present, 100, 0))  
}

author_presence <- function(authors_B) {
  authors_B_list <- str_split(authors_B, ",\\s*")[[1]]  
  is_present <- first_author %in% authors_B_list | corresponding_author %in% authors_B_list
  return(ifelse(is_present, 100, 0))  
}
combine_authors <- function(author_str, consortium_str) {
  
  if (is.na(author_str)) author_str <- ""
  if (is.na(consortium_str)) consortium_str <- ""
  
  authors_list <- str_split(author_str, ",\\s*")[[1]]
  consortium_list <- str_split(consortium_str, ",\\s*")[[1]]
  combined <- unique(c(authors_list, consortium_list))
  combined <- combined[combined != ""]
  return(combined)
}

combine_names <- function(name_strings) {
  name_strings <- name_strings[!is.na(name_strings) & name_strings != ""]
  if (length(name_strings) == 0) return("")
  
  names <- unlist(str_split(name_strings, ",\\s*"))
  names <- unique(trimws(names))
  names <- names[names != ""]
  paste(names, collapse = ", ")
}

author_presence_absence <- function(authors_A_combined, authors_B_combined) {
  shared_count <- length(intersect(authors_A_combined, authors_B_combined))
  return(ifelse(shared_count > 0, 1, 0))
}

author_similarity <- function(authors_A, authors_B) {
  A_list <- str_split(authors_A, ",\\s*")[[1]]
  B_list <- str_split(authors_B, ",\\s*")[[1]]
  
  shared_authors <- length(intersect(A_list, B_list))
  total_unique   <- length(unique(c(A_list, B_list)))
  
  similarity_percent <- ifelse(
    total_unique == 0, 
    0, 
    (shared_authors / total_unique) * 100
  )
  return(list(similarity_percent, shared_authors))
}

has_overlapping_author <- function(authors_str, consortium_str, training_list) {
  people <- c(
    if (!is.na(authors_str)) str_split(authors_str, ",\\s*")[[1]] else character(0),
    if (!is.na(consortium_str) && consortium_str != "") str_split(consortium_str, ",\\s*")[[1]] else character(0)
  )
  any(trimws(people) %in% training_list) %>% as.integer()
}


###########ALDEX2##############
run_aldex2 <- function(countdata, conds, scaledata = NULL, feature_names = NULL,
                      scale_label = NULL, uncertainty = NULL) {
  if ("SampleID" %in% colnames(conds)) {
    conds <- conds[, !(colnames(conds) %in% "SampleID")]
  }
  conds <- as.data.frame(conds)
  covariates <- data.frame(lapply(conds, as.factor))
  num_covariates <- ncol(covariates)
  if (num_covariates == 1) {
    num_levels <- length(unique(covariates[[1]]))
    
    if (num_levels == 2) {
      test <- "t"  
    } else {
      test <- "kw" 
    }
  } else if (num_covariates > 1) {
    test <- "glm"  
    formula_str <- paste("~", paste(names(covariates), collapse = " + "))
    model_formula <- as.formula(formula_str)
    conditions <- model.matrix(model_formula, data = covariates)
  } else {
    stop("No valid covariates available.")
  }
  message("Test selected: ", test)
  if (!is.null(scaledata)) {
    if (is.null(uncertainty)) {
      uncertainty <- rep(1e-9, length(scaledata))  
    }
    scale_samps <- matrix(
      rnorm(length(scaledata) * 128, mean = rep(scaledata, each = 128), 
            sd = rep(uncertainty, each = 128)), 
      nrow = length(scaledata), ncol = 128
    )
    clr <- aldex.clr(countdata,
                     conditions,
                     denom = "all",
                     mc.samples = 128,
                     gamma = scale_samps)
    mod <- aldex.glm(clr)
  } else {
    gamma_to_test <- c(1e-3, 0.1, 0.25, 0.5)
    clr <- aldex.clr(countdata,
                     conditions,
                     denom = "all",
                     mc.samples = 128)
    
    mod <- aldex.senAnalysis(clr, test = test, gamma = gamma_to_test)
  }
  if (!is.null(clr)) {
    clr_values <- as.data.frame(do.call(cbind, clr@analysisData))
    rownames(clr_values) <- rownames(clr@reads) 
    if (is.matrix(conditions)) {
      unique_groups <- colnames(conditions) 
    } else {
      unique_groups <- unique(conds[[1]]) 
    }
    pairwise_comparisons <- combn(unique_groups, 2, simplify=FALSE)
    lfc_list <- lapply(pairwise_comparisons, function(pair) {
      group1 <- pair[1]
      group2 <- pair[2]
      if (is.matrix(conditions)) {
        if (!(group1 %in% colnames(conditions)) || !(group2 %in% colnames(conditions))) {
          return(NULL)  
        }
        group1_samples <- which(conditions[, group1] == 1)
        group2_samples <- which(conditions[, group2] == 1)
      } else {
        group1_samples <- which(as.character(conds[[1]]) == group1)
        group2_samples <- which(as.character(conds[[1]]) == group2)
      }
      if (length(group1_samples) == 0 || length(group2_samples) == 0) {
        message("Skipping comparison: ", group1, " vs ", group2, " (group missing)")
        return(NULL)
      }
      clr_group1 <- clr_values[, group1_samples, drop=FALSE]
      clr_group2 <- clr_values[, group2_samples, drop=FALSE]
      mean_lfc <- rowMeans(clr_group2, na.rm=TRUE) - rowMeans(clr_group1, na.rm=TRUE)
      return(data.frame(
        Feature = rownames(clr_values),
        Group1 = group1,
        Group2 = group2,
        LogFoldChange = mean_lfc
      ))
    })
  } else {
    return(NULL)
  }
  
  ### **Check if ALDEx2 returned results**
  if (is.null(mod) || length(mod) == 0) {
    message("Warning: ALDEx2 returned no valid results")
    return(NULL)
  }
  
  if (!is.data.frame(mod)) {
    result_list <- lapply(names(mod), function(gamma_val) {
      res <- mod[[gamma_val]]
      
      if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
        return(NULL)
      }
      
      if (!is.null(rownames(res))) {
        feature_count <- length(rownames(res))
        if (!is.null(feature_names) && feature_count == length(feature_names)) {
          res$Feature <- feature_names  
        } else {
          message("Warning: Feature count mismatch! Using original rownames.")
          res$Feature <- rownames(res)
        }
      }
      
      data.frame(
        Dataset = rep(study_name, nrow(res)),
        Externalscalemeasurement = rep(scale_label, nrow(res)),
        Uncertainty = rep(if (is.null(uncertainty)) NA else uncertainty, nrow(res)),
        Gamma = rep(gamma_val, nrow(res)),  
        res,
        stringsAsFactors = FALSE
      )
    })
  } else {
    
    res <- mod 
    
    if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
      message("Warning: ALDEx2 returned no valid results")
      return(NULL)
    }
    
    if (!is.null(rownames(res))) {
      feature_count <- length(rownames(res))
      if (!is.null(feature_names) && feature_count == length(feature_names)) {
        res$Feature <- feature_names 
      } else {
        message("Warning: Feature count mismatch! Using original rownames.")
        res$Feature <- rownames(res)
      }
    }
    
    result_list <- list(data.frame(
      Dataset = rep(study_name, nrow(res)),
      Externalscalemeasurement = rep(scale_label, nrow(res)),
      Uncertainty = rep(if (is.null(uncertainty)) NA else uncertainty, nrow(res)),
      Gamma = rep(paste0(scale_label,"_",uncertainty), nrow(res)), 
      res,
      stringsAsFactors = FALSE
    ))
  }
  
  result_list <- Filter(Negate(is.null), result_list)
  if (length(result_list) == 0) {
    return(NULL)
  }
  
  return(list(ALDEx_results = do.call(rbind, result_list), LogFoldChanges = lfc_list))
}

process_study_with_aldex2 <- function(study_name, study_data, sequencing_type) {
  
  if (!is.null(study_data$counts_df)) {
    feature_names <- colnames(study_data$counts_df)
    counts_df <- t(study_data$counts_df)
    counts_df <- matrix(as.integer(counts_df),
                        nrow = nrow(counts_df),
                        ncol = ncol(counts_df))
  } else {
    return(NULL)
  }
  
  message("Differential Abundance Analysis for: ", study_name)
  
  measured_df <- study_data$measured_load_df
  metadatacondition <- unique(
    datasets[datasets$dataset_name == study_name & datasets$abundance == "metadata", 
             c("covariate1", "covariate2", "covariate3", "covariate4", "covariate5", "covariate6", "covariate7")]
  )
  
  metadatacondition <- as.data.frame(metadatacondition)
  metadatacondition <- metadatacondition[, !(colnames(metadatacondition) %in% c("", " "))]
  metadatacondition <- metadatacondition[, colSums(!is.na(metadatacondition)) > 0, drop = FALSE]
  metadatacondition <- metadatacondition[, colSums(metadatacondition != "" & metadatacondition != " ") > 0, drop = FALSE]
  
  if (ncol(metadatacondition) == 0) {
    message("Skipping study: ", study_name, " - No valid covariates available")
    return(NULL)
  }
  values_from_metadatacondition <- as.character(unlist(metadatacondition))
  conds <- study_data$metadata_df[, values_from_metadatacondition, drop = FALSE]
  valid_samples <- rowSums(!is.na(conds)) > 0
  counts_df <- counts_df[, valid_samples, drop = FALSE]
  conds <- conds[valid_samples, , drop = FALSE]
  conds <- data.frame(lapply(conds, function(x) ifelse(is.na(x), "Missing", as.character(x))), 
                      stringsAsFactors = TRUE)
  
  if (!is.null(measured_df)) {
    measured_df <- measured_df[valid_samples, , drop = FALSE]
  }
  
  ALDEx <- NULL
  measured <- NULL
  measured_with_added_uncertainty <- NULL
  galaxy <- NULL
  galaxy_with_added_uncertainty <- NULL
  metacardis <- NULL
  metacardis_with_added_uncertainty <- NULL
  vandeputte <- NULL
  vandeputte_with_added_uncertainty <- NULL
  
  message("Run Sensitivity Analysis without scale measurements")
  ALDEx <- run_aldex2(counts_df, conds, NULL,
                               feature_names,
                               "None", uncertainty = NULL)
  
  if (!is.null(measured_df)) {
    measured_load <- measured_df$MeasuredLoad
    message("Run measured scale without uncertainty")
    measured <- run_aldex2(counts_df, conds, measured_load,
                          feature_names,
                          "Measured", uncertainty = 1e-9)
    
    message("Run measured scale with uncertainty")
    measured_with_added_uncertainty <- run_aldex2(counts_df, conds, measured_load,
                                                 feature_names,
                                                 "Measured", uncertainty = 0.1)
  }
  
  if (sequencing_type == "shotgun metagenomics") {
    predictions <- list(
      galaxy     = study_data$galaxypredicted_load_df$PredictedLoad,
      metacardis = study_data$metacardispredicted_load_df$PredictedLoad,
      vandeputte = NULL
    )
  } else {
    predictions <- list(
      galaxy     = NULL,
      metacardis = NULL,
      vandeputte = study_data$predicted_load_df$PredictedLoad
    )
  }
  
  for (pred in names(predictions)) {
    if (!is.null(predictions[[pred]])) {
      predictions[[pred]] <- predictions[[pred]][valid_samples]
    }
  }
  
  for (pred_name in names(predictions)) {
    pred_data <- predictions[[pred_name]]
    if (!is.null(pred_data)) {
      message("Run predicted scale without uncertainty")
      assign(pred_name, run_aldex2(counts_df, conds, pred_data, 
                                  feature_names,
                                  pred_name, uncertainty = 1e-9))
      
      message("Run predicted scale with uncertainty")
      assign(paste0(pred_name, "_with_added_uncertainty"), 
             run_aldex2(counts_df, conds, pred_data,
                       feature_names,
                       pred_name, uncertainty = 0.1))
    }
  }
  
  final_results_list <- list(
    ALDEx = ALDEx,
    measured = measured,
    measured_with_added_uncertainty = measured_with_added_uncertainty,
    galaxy = galaxy,
    galaxy_with_added_uncertainty = galaxy_with_added_uncertainty,
    metacardis = metacardis,
    metacardis_with_added_uncertainty = metacardis_with_added_uncertainty,
    vandeputte = vandeputte,
    vandeputte_with_added_uncertainty = vandeputte_with_added_uncertainty
  )
  
  aldex_results_list <- lapply(final_results_list, function(x) if (!is.null(x$ALDEx_results)) x$ALDEx_results else data.frame())
  
  if (all(sapply(aldex_results_list, nrow) == 0)) {
    return(NULL)
  }
  
  lfc_results_list <- lapply(names(final_results_list), function(study_name) {
    x <- final_results_list[[study_name]]
    
    tryCatch({
      if (!is.null(x$LogFoldChanges) && is.data.frame(x$LogFoldChanges)) {
        x$LogFoldChanges$Study <- study_name  # Add study name to each row
        return(x$LogFoldChanges)
      } else {
        return(data.frame())  # Return an empty dataframe
      }
    }, error = function(e) {
      message("Error processing LogFoldChanges for study: ", study_name, " - ", e$message)
      return(data.frame())  # Return an empty dataframe on failure
    })
  })
  
  return(list(
    ALDEx_results = data.table::rbindlist(aldex_results_list, fill = TRUE),
    LogFoldChanges = data.table::rbindlist(lfc_results_list, fill = TRUE),
    Total = final_results_list
  ))
}

########## ALDEX3 ##############################################################

# run_aldex3 <- function(countdata, conds, scaledata = NULL,
#                       scale_label = NULL, study_name = NULL) {
#   
#   if ("SampleID" %in% colnames(conds)) {
#     conds <- conds[, !(colnames(conds) %in% "SampleID")]
#   }
#   
#   covariates <- as.data.frame(conds)
#   num_covariates <- ncol(covariates)
#   
#   if (num_covariates == 0) {
#     stop("No valid covariates available.")
#   }
#   
#   formula_str <- paste("~", paste(names(covariates), collapse = " + "))
#   modelformula <- as.formula(formula_str)
#   message(modelformula)
#   aldex_out <- NULL
# 
#   if (!is.null(scaledata)) {
#     aldex_out <- aldex.lm(Y=countdata, X=modelformula, data=covariates,stream=3000,
#                           scale=function(logWpara, X, Y) {
#                             scalesensitivityexternal(logWpara, X, externalscale = scaledata, gamma = c(1e-9,0.5, 1, 5))
#                           })
#   } else {
#     aldex_out <- aldex.lm(Y=countdata, X=modelformula, data=covariates,stream=3000,
#                           scale=function(logWpara, X, Y) {
#                             scalesensitivity(logWpara, X, gamma = c(1e-9,0.5, 1, 5))
#                           })
#   }
# 
#   if (is.null(aldex_out) || length(aldex_out) == 0) {
#     message("Warning: No valid results from ALDEx analysis (placeholder).")
#     return(NULL)
#   }
# 
#   flatten_aldex_results <- function(aldex_out, study_name, scale_label = "ZSM", use = c("mean", "median")) {
#     use <- match.arg(use)
#     
#     results <- lapply(names(aldex_out), function(gamma) {
#       res <- aldex_out[[gamma]]
#       
#       # Summarize posterior estimates
#       est <- switch(
#         use,
#         mean = apply(res$estimate, c(1, 2), mean),
#         median = apply(res$estimate, c(1, 2), median)
#       )
#       
#       # Always drop the first row (intercept)
#       est <- est[-1, , drop = FALSE]
#       p.val <- res$p.val[-1, , drop = FALSE]
#       p.val.adj <- res$p.val.adj[-1, , drop = FALSE]
#       
#       # Assign generic covariate and feature names if needed
#       rownames(est) <- paste0("covariate_", seq_len(nrow(est)))
#       colnames(est) <- paste0("feature_", seq_len(ncol(est)))
#       rownames(p.val) <- rownames(p.val.adj) <- rownames(est)
#       colnames(p.val) <- colnames(p.val.adj) <- colnames(est)
#       
#       # Tidy format
#       estimate_df <- as.data.frame(est) %>%
#         rownames_to_column("covariate") %>%
#         pivot_longer(-covariate, names_to = "feature", values_to = "estimate")
#       
#       pval_df <- as.data.frame(p.val) %>%
#         rownames_to_column("covariate") %>%
#         pivot_longer(-covariate, names_to = "feature", values_to = "p.val")
#       
#       padj_df <- as.data.frame(p.val.adj) %>%
#         rownames_to_column("covariate") %>%
#         pivot_longer(-covariate, names_to = "feature", values_to = "p.val.adj")
#       
#       # Merge and label
#       final_df <- estimate_df %>%
#         left_join(pval_df, by = c("covariate", "feature")) %>%
#         left_join(padj_df, by = c("covariate", "feature")) %>%
#         mutate(
#           gamma = as.numeric(gamma),
#           study_name = study_name,
#           scale_label = scale_label
#         )
#       
#       return(final_df)
#     })
#     
#     bind_rows(results)
#   }
#   
#   aldex_results = flatten_aldex_results(aldex_out, study_name, scale_label, use="mean")
#   return(ALDEx3_results = aldex_results)
# }


run_aldex3 <- function(countdata, conds, scaledata = NULL,
                       scale_label = NULL, study_name = NULL,
                       gammas = c(1e-9, 0.5, 1, 5)) {
  
  if ("SampleID" %in% colnames(conds)) {
    conds <- conds[, !(colnames(conds) %in% "SampleID")]
  }
  
  covariates <- as.data.frame(conds)
  num_covariates <- ncol(covariates)
  
  if (num_covariates == 0) {
    stop("No valid covariates available.")
  }
  
  formula_str <- paste("~", paste(names(covariates), collapse = " + "))
  modelformula <- as.formula(formula_str)
  message(modelformula)
  
  aldex_out <- list()
  
  for (g in gammas) {
    message(sprintf("Running ALDEx with gamma = %s", g))
    if (!is.null(scaledata)) {
      externalscale <- function(X, logWpara, gamma = 0.5, scaledata) {
        N <- length(scaledata)
        nsample <- dim(logWpara)[3]
        gamma_vec <- if (length(gamma) == 1) rep(gamma, N) else gamma
        logWperp <- t(sapply(seq_len(N), function(i) {
          rnorm(nsample, mean = scaledata[i], sd = gamma_vec[i])
        }))
        return(logWperp)
      }
      
      result <- aldex(
        Y = countdata, 
        X = modelformula, 
        data = covariates,
        stream = 3000,
        scale = externalscale,
        gamma = g, 
        scaledata = scaledata
      )
    } else {
      tss <- function(X, logWpara, gamma=0.5) {
         P <- nrow(X)
         nsample <- dim(logWpara)[3]
         Lambdaperp <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
         logWperp <- t(X)%*% Lambdaperp
         return(logWperp)
      }
      result <- aldex(
        Y = countdata, 
        X = modelformula, 
        data = covariates,
        stream = 3000,
        scale = tss,
        gamma = g
      )
    }
    
    aldex_out[[as.character(g)]] <- result
  }
  
  if (length(aldex_out) == 0) {
    message("Warning: No valid results from ALDEx analysis.")
    return(NULL)
  }
  
  aldex_results <- flatten_aldex_results(aldex_out, countdata, study_name, scale_label, use = "mean")
  return(ALDEx3_results = aldex_results)
}

process_study_with_aldex3 <- function(study_name, study_data, sequencing_type) {
  if (!is.null(study_data$counts_df)) {
    feature_names <- colnames(study_data$counts_df)
    counts_df <- t(study_data$counts_df)
    counts_df <- matrix(as.integer(counts_df),
                        nrow = nrow(counts_df),
                        ncol = ncol(counts_df))
  } else {
    return(NULL)
  }
  message("Differential Abundance Analysis for: ", study_name)
  measured_df <- study_data$measured_load_df
  metadatacondition <- unique(
    datasets[datasets$dataset_name == study_name & datasets$abundance == "metadata", 
             c("covariate1", "covariate2", "covariate3", "covariate4", "covariate5", "covariate6", "covariate7")]
  )
  metadatacondition <- as.data.frame(metadatacondition)
  metadatacondition <- metadatacondition[, !(colnames(metadatacondition) %in% c("", " "))]
  metadatacondition <- metadatacondition[, colSums(!is.na(metadatacondition)) > 0, drop = FALSE]
  metadatacondition <- metadatacondition[, colSums(metadatacondition != "" & metadatacondition != " ") > 0, drop = FALSE]
  if (ncol(metadatacondition) == 0) {
    message("Skipping study: ", study_name, " - No valid covariates available")
    return(NULL)
  }
  values_from_metadatacondition <- as.character(unlist(metadatacondition))
  conds <- study_data$metadata_df[, values_from_metadatacondition, drop = FALSE]
  valid_samples <- rowSums(!is.na(conds)) > 0
  counts_df <- counts_df[, valid_samples, drop = FALSE]
  conds <- conds[valid_samples, , drop = FALSE]
  message(str(conds))
  
  if (!is.null(measured_df)) {
    measured_df <- measured_df[valid_samples, , drop = FALSE]
  }
  ALDEx <- NULL
  measured <- NULL
  GALAXY <- NULL
  MetaCardis <- NULL
  vandeputte <- NULL
  message("Run Sensitivity Analysis without scale measurements")
  ALDEx <- run_aldex3(counts_df, conds, scaledata=NULL, scale_label = "ALDEx3", 
                                study_name = study_name)
  if (!is.null(measured_df)) {
    measured_load <- measured_df$MeasuredLoad
    message("Run measured scale")
    measured <- run_aldex3(counts_df, conds, scaledata=measured_load, 
                           scale_label = "Measured", study_name = study_name)
  }
  if (sequencing_type == "shotgun metagenomics") {
    predictions <- list(
      GALAXY     = study_data$galaxypredicted_load_df$PredictedLoad,
      MetaCardis = study_data$metacardispredicted_load_df$PredictedLoad,
      vandeputte = NULL
    )
  } else {
    predictions <- list(
      GALAXY     = NULL,
      MetaCardis = NULL,
      vandeputte = study_data$predicted_load_df$PredictedLoad
    )
  }
  for (pred in names(predictions)) {
    if (!is.null(predictions[[pred]])) {
      predictions[[pred]] <- predictions[[pred]][valid_samples]
    }
  }
  for (pred_name in names(predictions)) {
    pred_data <- predictions[[pred_name]]
    if (!is.null(pred_data)) {
      message("Run predicted scale")
      assign(pred_name, run_aldex3(counts_df, conds, scaledata=pred_data, 
                                   scale_label = pred_name, study_name = study_name))
    }
  }
  final_results_list <- list(
    ALDEx = ALDEx,
    measured = measured,
    GALAXY = GALAXY,
    MetaCardis = MetaCardis,
    vandeputte = vandeputte
  )
  aldex_results_list <- lapply(final_results_list, function(x) if (!is.null(x)) x else data.frame())
  if (all(sapply(aldex_results_list, nrow) == 0)) {
    return(NULL)
  }
  return(aldex_results_list)
}

run_aldex3_with_benchmarks <- function(countdata, conds, scaledata = NULL,
                                       scale_label = NULL, study_name = NULL,
                                       gammas = c(0, 0.5, 1, 5)) {
  
  if ("SampleID" %in% base::colnames(conds)) {
    conds <- conds[, !(base::colnames(conds) %in% "SampleID")]
  }
  
  covariates <- base::as.data.frame(conds)
  num_covariates <- base::ncol(covariates)
  
  if (num_covariates == 0) {
    base::stop("No valid covariates available.")
  }
  
  base::colnames(covariates) <- base::make.names(base::colnames(covariates), unique = TRUE)
  for (col in base::colnames(covariates)) {
    if (base::is.character(covariates[[col]]) || base::is.factor(covariates[[col]])) {
      covariates[[col]] <- base::factor(covariates[[col]])
    } else if (base::is.numeric(covariates[[col]])) {
      covariates[[col]] <- base::scale(covariates[[col]], center = TRUE, scale = TRUE)
    }
  }
  
  formula_str <- base::paste("~", base::paste(base::names(covariates), collapse = " + "))
  modelformula <- as.formula(formula_str)
  base::message("Design formula: ", formula_str)
  
  zero_count_rows <- rowSums(countdata) == 0
  if (sum(zero_count_rows) > 0) {
    base::message(paste("Found", sum(zero_count_rows), 
                        "features with zero counts across all samples. Filtering them out..."))
    countdata <- countdata[!zero_count_rows, ]
  } else {
    base::message("No features with zero counts across all samples found. Proceeding without filtering.")
  }
  
  results <- base::list()
  aldex_out <- base::list()
  
  # ALDEx3 analysis
  if (base::is.null(scaledata)) {
    base::message("Running ALDEx3 without scale data")
    aldex_out_tss <- list()
    aldex_out_clr <- list()
    for (g in gammas) {
      base::message(base::sprintf("Running ALDEx3 with tss scale, gamma = %s", g))
      tss <- function(X, logWpara, gamma = 0.5) {
        P <- base::nrow(X)
        nsample <- base::dim(logWpara)[3]
        if (gamma == 0) {
          Lambdaperp <- base::matrix(0, P, nsample)
        } else {
          Lambdaperp <- base::matrix(stats::rnorm(P * nsample, 0, gamma), P, nsample)
        }
        logWperp <- base::t(X) %*% Lambdaperp
        return(logWperp)
      }
      result_tss <- aldex(
        Y = countdata, 
        X = modelformula, 
        data = covariates,
        stream = 3000,
        scale = tss,
        gamma = g
      )
      aldex_out_tss[[base::as.character(g)]] <- result_tss
      base::message(base::sprintf("Running ALDEx3 with clr scale, gamma = %s", g))
      clr <- function(X, logWpara, gamma = 0.5) {
        P <- base::nrow(X)
        nsample <- base::dim(logWpara)[3]
        logWperp <- -base::colMeans(logWpara, dims = 1)
        if (gamma == 0) {
          Lambdaperp <- base::matrix(0, P, nsample)
        } else {
          tmp <- P * nsample
          Lambdaperp <- base::matrix(stats::rnorm(tmp, 0, gamma), P, nsample)
        }
        logWperp <- logWperp + base::t(X) %*% Lambdaperp
        return(logWperp)
      }
      result_clr <- aldex(
        Y = countdata, 
        X = modelformula, 
        data = covariates,
        stream = 3000,
        scale = clr,
        gamma = g
      )
      aldex_out_clr[[base::as.character(g)]] <- result_clr
    }
    if (base::length(aldex_out_tss) > 0) {
      results$ALDEx3_tss <- flatten_aldex_results(
        aldex_out_tss, 
        countdata, 
        study_name, 
        scale_label = "ALDEx3_tss", 
        use = "mean"
      )
    } else {
      base::message("Warning: No valid results from ALDEx3 tss analysis.")
      results$ALDEx3_tss <- base::data.frame()
    }
    if (base::length(aldex_out_clr) > 0) {
      results$ALDEx2_clr <- flatten_aldex_results(
        aldex_out_clr, 
        countdata, 
        study_name, 
        scale_label = "ALDEx2_clr", 
        use = "mean"
      )
    } else {
      base::message("Warning: No valid results from ALDEx3 clr analysis.")
      results$ALDEx3_clr <- base::data.frame()
    }
  } else {
    base::message("Running ALDEx3 with scale data")
    externalscale <- function(X, logWpara, gamma = 0.5, scaledata) {
      N <- base::length(scaledata)
      nsample <- base::dim(logWpara)[3]
      gamma_vec <- if (base::length(gamma) == 1) base::rep(gamma, N) else gamma
      
      logWperp <- base::t(base::sapply(base::seq_len(N), function(i) {
        if (gamma == 0) {
          base::rep(scaledata[i], nsample)
        } else {
          stats::rnorm(nsample, mean = scaledata[i], sd = gamma_vec[i])
        }
      }))
      return(logWperp)
    }
    for (g in gammas) {
      base::message(base::sprintf("Running ALDEx with gamma = %s", g))
      result <- aldex(
        Y = countdata, 
        X = modelformula, 
        data = covariates,
        stream = 3000,
        scale = externalscale,
        gamma = g,
        scaledata = scaledata
      )
      aldex_out[[base::as.character(g)]] <- result
    }
    if (base::length(aldex_out) > 0) {
      results$ALDEx3 <- flatten_aldex_results(aldex_out, countdata, study_name, 
                                              scale_label = scale_label, use = "mean")
    } else {
      base::message("Warning: No valid results from ALDEx3 analysis.")
      results$ALDEx3 <- base::data.frame()
    }
  }
  
  # DESeq2 analysis
  if (base::is.null(scaledata)) {
    base::message("Running DESeq2 without scale data")
    for (col in base::colnames(covariates)) {
      if (base::is.ordered(covariates[[col]])) {
        covariates[[col]] <- base::factor(covariates[[col]], ordered = FALSE)
      }
    }
    design <- stats::model.matrix(modelformula, data = covariates)
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = countdata,  
      colData = covariates,
      design = design
    )
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::DESeq(dds)
    coef_names <- DESeq2::resultsNames(dds)
    design_names <- base::colnames(design)[-1]  
    
    if (length(design_names) != length(coef_names) - 1) {
      base::stop("Mismatch between design matrix columns and DESeq2 coefficients.")
    }
    coef_map <- c("Intercept", design_names)
    names(coef_map) <- coef_names
    deseq_results <- base::lapply(coef_names[-1], function(cov) {
      res <- DESeq2::results(dds, name = cov)
      res_df <- base::as.data.frame(res) %>%
        tibble::rownames_to_column(var = "feature") %>%
        dplyr::mutate(feature = paste0("V", feature)) %>%  
        dplyr::select(feature, log2FoldChange, pvalue, padj) %>%
        dplyr::rename(estimate = log2FoldChange, p.val = pvalue, p.val.adj = padj) %>%
        dplyr::mutate(covariate = unname(as.character(coef_map[cov])), 
                      study_name = study_name, 
                      scale_label = "DESeq2", 
                      gamma = "0")
      return(res_df)
    })
    results$DESeq2 <- dplyr::bind_rows(deseq_results)
  } else {
    results$DESeq2 <- base::data.frame()  
  }
  
  # limma analysis
  if (base::is.null(scaledata)) {
    base::message("Running limma with voom without scale data")
    design <- stats::model.matrix(modelformula, data = covariates)
    v <- limma::voom(countdata, design, plot = FALSE) 
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    
    limma_results <- base::lapply(base::colnames(design)[-1], function(cov) {
      res <- limma::topTable(fit, coef = cov, number = Inf, adjust.method = "BH")
      res_df <- base::as.data.frame(res) %>%
        dplyr::mutate(feature = paste0("V",base::rownames(res))) %>%
        dplyr::select(feature, logFC, P.Value, adj.P.Val) %>%
        dplyr::rename(estimate = logFC, p.val = P.Value, p.val.adj = adj.P.Val) %>%
        dplyr::mutate(covariate = cov, study_name = study_name, scale_label = "limma",
                      gamma="0")
      return(res_df)
    })
    results$limma <- dplyr::bind_rows(limma_results)
  } else {
    results$limma <- base::data.frame()  
  }
  
  return(results)
}

process_study_with_benchmarks <- function(study_name, study_data, sequencing_type) {
  if (!base::is.null(study_data$counts_df)) {
    feature_names <- base::colnames(study_data$counts_df)
    counts_df <- base::t(study_data$counts_df)
    counts_df <- base::matrix(base::as.integer(counts_df),
                              nrow = base::nrow(counts_df),
                              ncol = base::ncol(counts_df))
  } else {
    return(NULL)
  }
  
  base::message("Differential Abundance Analysis for: ", study_name)
  measured_df <- study_data$measured_load_df
  metadatacondition <- base::unique(
    datasets[datasets$dataset_name == study_name & datasets$abundance == "metadata", 
             c("covariate1", "covariate2", "covariate3", "covariate4", "covariate5", "covariate6", "covariate7")]
  )
  metadatacondition <- base::as.data.frame(metadatacondition)
  metadatacondition <- metadatacondition[, !(base::colnames(metadatacondition) %in% c("", " "))]
  metadatacondition <- metadatacondition[, base::colSums(!base::is.na(metadatacondition)) > 0, drop = FALSE]
  metadatacondition <- metadatacondition[, base::colSums(metadatacondition != "" & metadatacondition != " ") > 0, drop = FALSE]
  
  if (base::ncol(metadatacondition) == 0) {
    base::message("Skipping study: ", study_name, " - No valid covariates available")
    return(NULL)
  }
  
  values_from_metadatacondition <- base::as.character(base::unlist(metadatacondition))
  conds <- study_data$metadata_df[, values_from_metadatacondition, drop = FALSE]
  valid_samples <- base::rowSums(!base::is.na(conds)) > 0
  counts_df <- counts_df[, valid_samples, drop = FALSE]
  conds <- conds[valid_samples, , drop = FALSE]
  base::message(str(conds))
  
  if (!base::is.null(measured_df)) {
    measured_df <- measured_df[valid_samples, , drop = FALSE]
  }
  
  final_results_list <- base::list(
    ALDEx3_tss = NULL,
    ALDEx2_clr = NULL,
    DESeq2 = NULL,
    limma = NULL,
    measured = NULL,
    GALAXY = NULL,
    MetaCardis = NULL,
    vandeputte = NULL
  )
  
  base::message("Run Sensitivity Analysis without scale measurements")
  no_scale_results <- run_aldex3_with_benchmarks(counts_df, conds, scaledata = NULL, 
                                                 scale_label = "ALDEx3", study_name = study_name)
  final_results_list$ALDEx3 <- no_scale_results$ALDEx3
  final_results_list$DESeq2 <- no_scale_results$DESeq2
  final_results_list$limma <- no_scale_results$limma
  
  if (!base::is.null(measured_df)) {
    measured_load <- measured_df$MeasuredLoad
    base::message("Run measured scale")
    measured_results <- run_aldex3_with_benchmarks(counts_df, conds, scaledata = measured_load, 
                                                   scale_label = "Measured", study_name = study_name)
    final_results_list$measured <- measured_results$ALDEx3
  }
  
  if (sequencing_type == "shotgun metagenomics") {
    predictions <- base::list(
      GALAXY     = study_data$galaxypredicted_load_df$PredictedLoad,
      MetaCardis = study_data$metacardispredicted_load_df$PredictedLoad,
      vandeputte = NULL
    )
  } else {
    predictions <- base::list(
      GALAXY     = NULL,
      MetaCardis = NULL,
      vandeputte = study_data$predicted_load_df$PredictedLoad
    )
  }
  
  for (pred in base::names(predictions)) {
    if (!base::is.null(predictions[[pred]])) {
      predictions[[pred]] <- predictions[[pred]][valid_samples]
    }
  }
  
  for (pred_name in base::names(predictions)) {
    pred_data <- predictions[[pred_name]]
    if (!base::is.null(pred_data)) {
      base::message("Run predicted scale: ", pred_name)
      pred_results <- run_aldex3_with_benchmarks(counts_df, conds, scaledata = pred_data, 
                                                 scale_label = pred_name, study_name = study_name)
      final_results_list[[pred_name]] <- pred_results$ALDEx3
    }
  }
  
  aldex_results_list <- base::lapply(final_results_list, function(x) if (!base::is.null(x)) x else base::data.frame())
  
  if (base::all(base::sapply(aldex_results_list, base::nrow) == 0)) {
    return(NULL)
  }
  
  return(aldex_results_list)
}

flatten_aldex_results <- function(aldex_out, countdata, study_name, scale_label = "ZSM", use = c("mean", "median")) {
  use <- base::match.arg(use)
  results <- base::lapply(base::names(aldex_out), function(gamma) {
    res <- aldex_out[[gamma]]
    est <- switch(
      use,
      mean = base::apply(res$estimate, c(1, 2), base::mean),
      median = base::apply(res$estimate, c(1, 2), stats::median)
    )
    est <- est[-1, , drop = FALSE]  # Remove intercept
    p.val <- res$p.val[-1, , drop = FALSE]
    p.val.adj <- res$p.val.adj[-1, , drop = FALSE]
    covariate_names <- base::rownames(res$X)[-1]
    feature_names <- base::colnames(countdata)
    base::rownames(est) <- base::rownames(p.val) <- base::rownames(p.val.adj) <- covariate_names
    base::colnames(est) <- base::colnames(p.val) <- base::colnames(p.val.adj) <- feature_names
    estimate_df <- base::as.data.frame(est) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "estimate")
    pval_df <- base::as.data.frame(p.val) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "p.val")
    padj_df <- base::as.data.frame(p.val.adj) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "p.val.adj")
    out_df <- estimate_df %>%
      dplyr::left_join(pval_df, by = c("covariate", "feature")) %>%
      dplyr::left_join(padj_df, by = c("covariate", "feature"))
    metadata <- res$data %>%
      dplyr::mutate(feature = base::rownames(res$data)) %>%
      dplyr::relocate(feature) 
    
    out_df <- out_df %>%
      dplyr::left_join(metadata, by = "feature") %>%
      dplyr::mutate(
        gamma = base::as.numeric(gamma),
        study_name = study_name,
        scale_label = scale_label
      )
    
    return(out_df)
  })
  dplyr::bind_rows(results)
}



lm.scale <- function(Y, X, data = NULL, scale = NULL, p.adjust.method = "BH") {
  N <- ncol(Y)
  D <- nrow(Y) 

  if (is.null(scale)) {
    # CLR-like approach: 
    #  1) add 0.5 pseudo-count
    #  2) compute proportions
    #  3) log2 transform
    #  4) average across features => colMeans(log2(proportions))
    #  5) multiply by -1
    closure <- t(t(Y + 0.5) / colSums(Y + 0.5))  # (D x N)
    log_closure <- log2(closure)
    Z <- -colMeans(log_closure)                 # length N
  } else {
    if (length(scale) != N) {
      stop("Length of 'scale' vector must match number of samples (ncol(Y)).")
    }
    Z <- scale
  }
  if (inherits(X, "formula")) {
    if (is.null(data)) {
      stop("Data frame must be supplied if X is a formula.")
    }
    X_mat <- model.matrix(X, data) 
  } else {
    X_mat <- t(X)
    if (nrow(X_mat) != N) {
      stop("X dimensions must match the number of samples.")
    }
  }
  
  fit <- lm(Z ~ X_mat - 1)
  coefs  <- coef(fit)
  ses    <- sqrt(diag(vcov(fit)))
  pvals  <- summary(fit)$coefficients[, 4]
  padj   <- p.adjust(pvals, method = p.adjust.method)

  return(list(
    estimate  = coefs,
    std.error = ses,
    p.val     = pvals,
    p.val.adj = padj
  ))
}

run_lm_scale <- function(countdata, conds, scaledata = NULL,
                         scale_label = NULL, study_name = NULL) {

  if ("SampleID" %in% colnames(conds)) {
    conds <- conds[, !(colnames(conds) %in% "SampleID"), drop = FALSE]
  }

  covariate_names <- colnames(conds)
  if (length(covariate_names) == 0) {
    stop("No valid covariates in 'conds'.")
  }
  formula_str    <- paste("~", paste(covariate_names, collapse = " + "))
  modelformula   <- as.formula(formula_str)
  message(modelformula)

  fit_out <- lm.scale(
    Y     = countdata,
    X     = modelformula,
    data  = conds,
    scale = scaledata  
  )

  result_df <- data.frame(
    variable       = names(fit_out$estimate),
    estimate.mean  = fit_out$estimate,
    std.error      = fit_out$std.error,
    p.val          = fit_out$p.val,
    p.val.adj      = fit_out$p.val.adj,
    study_name     = study_name,
    scale_label    = scale_label,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

process_study_with_lm_scale <- function(study_name, study_data, sequencing_type) {
  
  # 3A) Read the counts
  if (!is.null(study_data$counts_df)) {
    counts_df <- t(study_data$counts_df)
    counts_df <- matrix(as.integer(counts_df),
                        nrow = nrow(counts_df),
                        ncol = ncol(counts_df))
  } else {
    message("No counts found, returning NULL")
    return(NULL)
  }

  measured_df <- study_data$measured_load_df
  metadatacondition <- unique(
    datasets[datasets$dataset_name == study_name & datasets$abundance == "metadata",
             c("covariate1", "covariate2", "covariate3", 
               "covariate4", "covariate5", "covariate6", "covariate7")]
  )
  metadatacondition <- as.data.frame(metadatacondition)
  metadatacondition <- metadatacondition[, !(colnames(metadatacondition) %in% c("", " ")), drop = FALSE]
  metadatacondition <- metadatacondition[, colSums(!is.na(metadatacondition)) > 0, drop = FALSE]
  metadatacondition <- metadatacondition[, colSums(metadatacondition != "" & metadatacondition != " ") > 0, drop = FALSE]
  if (ncol(metadatacondition) == 0) {
    message("Skipping study: ", study_name, " - no valid covariates available")
    return(NULL)
  }
  values_from_metadatacondition <- as.character(unlist(metadatacondition))
  conds <- study_data$metadata_df[, values_from_metadatacondition, drop = FALSE]
  valid_samples <- rowSums(!is.na(conds)) > 0
  counts_df <- counts_df[, valid_samples, drop = FALSE]
  conds <- conds[valid_samples, , drop = FALSE]

  if (!is.null(measured_df)) {
    measured_df <- measured_df[valid_samples, , drop = FALSE]
  }
  
  # 1) CLR scale (no external scale -> uses closure approach)
  # message("Analyzing CLR-based scale (no external scale)")
  # CLR_result <- run_lm_scale(
  #   countdata  = counts_df,
  #   conds      = conds,
  #   scaledata  = NULL,
  #   scale_label= "CLR",
  #   study_name = study_name
  # )
  
  message("Analyzing Zero scale model (no external scale)")
  Zero_result <- run_lm_scale(
    countdata  = counts_df,
    conds      = conds,
    scaledata  = NULL,
    scale_label= "Zero Scale",
    study_name = study_name
  )
  
  measured_result <- data.frame()
  if (!is.null(measured_df) && "MeasuredLoad" %in% colnames(measured_df)) {
    message("Analyzing measured scale")
    measured_load  <- measured_df$MeasuredLoad
    measured_result <- run_lm_scale(
      countdata   = counts_df,
      conds       = conds,
      scaledata   = measured_load,
      scale_label = "Measured",
      study_name  = study_name
    )
  }

  GALAXY     <- data.frame()
  MetaCardis <- data.frame()
  vandeputte <- data.frame()
  
  if (sequencing_type == "shotgun metagenomics") {
    if (!is.null(study_data$galaxypredicted_load_df)) {
      galaxy_load <- study_data$galaxypredicted_load_df$PredictedLoad[valid_samples]
      message("Analyzing GALAXY predicted scale")
      GALAXY <- run_lm_scale(
        countdata   = counts_df,
        conds       = conds,
        scaledata   = galaxy_load,
        scale_label = "GALAXY",
        study_name  = study_name
      )
    }
    if (!is.null(study_data$metacardispredicted_load_df)) {
      metacardis_load <- study_data$metacardispredicted_load_df$PredictedLoad[valid_samples]
      message("Analyzing MetaCardis predicted scale")
      MetaCardis <- run_lm_scale(
        countdata   = counts_df,
        conds       = conds,
        scaledata   = metacardis_load,
        scale_label = "MetaCardis",
        study_name  = study_name
      )
    }
  } else {
    if (!is.null(study_data$predicted_load_df)) {
      vandeputte_load <- study_data$predicted_load_df$PredictedLoad[valid_samples]
      message("Analyzing vandeputte predicted scale")
      vandeputte <- run_lm_scale(
        countdata   = counts_df,
        conds       = conds,
        scaledata   = vandeputte_load,
        scale_label = "vandeputte",
        study_name  = study_name
      )
    }
  }

  final_results_list <- list(
    `Zero Scale` = Zero_result,
    measured   = measured_result,
    GALAXY     = GALAXY,
    MetaCardis = MetaCardis,
    vandeputte = vandeputte
  )

  if (all(sapply(final_results_list, nrow) == 0)) {
    message("All results are empty for study: ", study_name)
    return(NULL)
  }
  
  return(final_results_list)
}



# ==== Utility Functions =======================================================

read_and_fix_data <- function(file_path) {
  message("📂 Reading file: ", file_path)
  
  if (!file.exists(file_path)) {
    stop("❌ Error: File does not exist: ", file_path)
  }
  
  delimiters <- c(",", "\t", ";", " ")
  raw_data <- NULL
  
  for (delim in delimiters) {
    message("🔍 Trying delimiter: '", delim, "'")
    temp_data <- tryCatch(
      read.delim(file_path, sep = delim, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>% data.frame(check.names = F),
      error = function(e) NULL
    )
    
    if (!is.null(temp_data) && ncol(temp_data) > 1) {
      message("✅ Correct delimiter found: ", delim)
      raw_data <- temp_data
      break
    }
  }
  
  if (!is.null(raw_data) && ncol(raw_data) == 1) {
    message("⚠️ Still only one column. Retrying without headers...")
    temp_data <- read.delim(file_path, sep = delimiter, header = FALSE, stringsAsFactors = FALSE) %>% data.frame(check.names = F)
    
    if (!is.null(temp_data) && nrow(temp_data) > 1) {
      message("✅ Using first row as column names.")
      colnames(raw_data) <- as.character(unlist(temp_data[1, ]))
      raw_data <- temp_data[-1, ]
    }
  }
  
  return(raw_data)
}


find_best_match <- function(loadcol, df) {
  normalize_name <- function(name) {
    name <- gsub("[\\.\\/_]", " ", name)  # Replace ., /, _ with spaces
    name <- gsub("^X(\\d)", "\\1", name)  # Remove 'X' before numbers
    name <- tolower(trimws(name))  # Convert to lowercase and trim whitespace
    return(name)
  }

  generate_variations <- function(name) {
    name_variants <- c(
      name,
      gsub("[\\.\\/_]", " ", name),  # Replace ., /, _ with spaces
      gsub("[\\.\\/_]", ".", name),  # Replace with .
      gsub("[\\.\\/_]", "/", name),  # Replace with /
      gsub("[\\.\\/_]", "_", name),  # Replace with _
      gsub("[\\.\\/_]", "", name),   # Remove special characters
      gsub("^X(\\d)", "\\1", name)   # Remove leading "X" before numbers
    )

    name_variants <- unique(c(name_variants, tolower(name_variants), toupper(name_variants)))
    
    return(name_variants)
  }

  normalized_target <- normalize_name(loadcol)
  possible_names <- generate_variations(loadcol)

  actual_colnames <- colnames(df)
  normalized_colnames <- sapply(actual_colnames, normalize_name)

  exact_match <- actual_colnames[normalized_colnames %in% possible_names]
  if (length(exact_match) > 0) {
    return(exact_match[1])  
  }

  fuzzy_match <- agrep(normalized_target, normalized_colnames, value = TRUE, ignore.case = TRUE)
  if (length(fuzzy_match) > 0) {
    return(actual_colnames[which(normalized_colnames %in% fuzzy_match)][1])  
  }
  
  return(NA)  
}

log10_to_log2 <- function(x) {
  log2(10^x)
}

# ==================================================
# End of Custom Functions
# ==================================================
