# Analysis functions for mlpanalysis.Rmd

unzip_all_in_dir <- function(dir) {
  zip_files <- list.files(path = dir, pattern = "\\.zip$", full.names = TRUE)
  for (zip in zip_files) {
    message("Unzipping: ", zip)
    unzip(zip, exdir = dir)
  }
}

filter_studies <- function(repository, specs, demographics) {
  studies <- imap(repository, function(study, name) {
    filtered <- filter_study_by_specs(study, specs[[name]])
    filtered
  }) %>% compact()
  imap(studies, function(study, name) {
    if (name %in% names(demographics)) {
      study$studydemographics <- demographics[[name]]
    }
    study
  })
}

process_mlp_predictions <- function(datasets, repo) {
  old_warn <- options(warn = 1)
  on.exit(options(old_warn))
  run_16s <- function(df, path) {
    if (is.null(df) || length(df) < 1 || nrow(df) < 1) return(NULL)
    message("Processing ", paste(path, collapse = " → "))
    df <- df %>%
      mutate(across(
        everything(),
        ~ replace_na(.x, 0)
      ))
    df <- df[rowSums(df) > 0, colSums(df) > 0]
    if (nrow(df) < 1 || ncol(df) < 1) {
      message("No data to process")
      return(NULL)
    }
    total    <- length(as.matrix(df))
    nonzero  <- sum(df != 0, na.rm = TRUE)
    sparsity <- (total - nonzero) / total * 100
    row_sparsity <- apply(df, 1, function(x) mean(x == 0, na.rm = TRUE) * 100)
    rdpclassifier <- readRDS("inst/extdata/16S_rRNA/model.16S_rRNA.rds")
    referencetaxa <- names(rdpclassifier$trainingData)
    if (any(c("Shannon diversity", ".outcome") %in% referencetaxa)) {
      referencetaxa <- setdiff(referencetaxa, c("Shannon diversity", ".outcome"))
    }
    samplesparsity <- compute_sparsity_per_sample(df, referencetaxa)
    if (mean(samplesparsity == 1, na.rm = TRUE) > 0.95) {
      print(paste("Model features for dataset are almost completely sparse:", mean(samplesparsity == 1, na.rm = TRUE)))
    }
    overlap_prop    <- prop_presenttaxa_overlap_per_sample(df, referencetaxa)
    nonoverlap_prop <- prop_taxa_notusedinmodel_per_sample(df, referencetaxa)
    psuedocount <- compute_psuedocount(df, referencetaxa)
    observed <- compute_observed(df)
    shannon <- compute_shannon(df)
    simpson <- compute_simpson(df)
    evenness <- compute_evenness(df, observed)
    df2 <- model_relative_abundances(df, referencetaxa)
    alpha_df <- compute_alpha_diversity(df2)
    rare_taxa <- evaluate_rare_taxa(df2)
    #ace <- compute_ace(df2)
    berger_parker <- compute_berger_parker(df2)
    inv_simpson <- compute_inv_simpson(df2)
    #chao1  <- compute_chao1(df2)
    rank_abundance_slope <- compute_rank_abundance_slope(df2)
    #pareto_alpha <- compute_pareto_alpha(df2)
    kurtosis <- compute_kurtosis(df2)
    skewness <- compute_skewness(df2)
    top_taxa <- get_top_taxa_per_sample(df2, top_n = 1)
    msg <- capture.output(
      tmp_res <- MLP(df,"rdp_train_set_16 (16S rRNA)", "load"),
      type = "message"
    )
    res <- tmp_res
    if (!is.null(res)) {
      full_msg <- paste(msg, collapse = " ")
      m        <- regexec("([0-9]+)\\D*([0-9]+\\.?[0-9]*)\\s*%", full_msg)
      mm       <- regmatches(full_msg, m)[[1]]
      if (length(mm) >= 3) {
        res$model_says_sharedtaxa         <- as.integer(mm[2])
        res$model_says_sharedtaxa_percent <- as.numeric(mm[3])
      } else {
        res$model_says_sharedtaxa         <- NA_integer_
        res$model_says_sharedtaxa_percent <- NA_real_
      }
      res$dataset_sparsity_percent <- sparsity
      res$sample_sparsity_percent_alltaxa <- row_sparsity
      freqs <- table(res$load)
      mode_value <- as.numeric(names(freqs)[which.max(freqs)])
      mode_freq  <- max(freqs)
      mode_pct   <- mode_freq / length(res$load) * 100
      res$mode_of_load               <- mode_value
      res$modefrequencyofload        <- mode_freq
      res$modefrequencyofload_percent <- mode_pct
      res$profile      <- "rdp"
      res$training_set <- "Vandeputte2021"
      res$log10_mlp    <- log10(res$load)
      res$log2_mlp     <- log2(res$load)
      res <- res %>% rename(mlp = load)
      res <- cbind(
        res,
        sample_sparsity_percent = samplesparsity*100,
        observed = observed,
        shannon = shannon,
        simpson = simpson,
        evenness = evenness,
        psuedocount = psuedocount,
        alpha_df,
        rare_taxa_percent = rare_taxa/alpha_df$Observed_modelfeatures * 100,
        taxa_used_in_model_percent    = overlap_prop*100,
        taxa_not_used_in_model_percent = nonoverlap_prop*100,
        #ace = ace,
        berger_parker = berger_parker,
        inv_simpson = inv_simpson,
        #chao1 = chao1,
        rank_abundance_slope = rank_abundance_slope,
        #pareto_alpha = pareto_alpha,
        kurtosis = kurtosis,
        skewness = skewness,
        top_taxa = top_taxa
      )
    }
    res
  }
  run_shotgun <- function(df, path) {
    if (is.null(df) || length(df) < 1 || nrow(df) < 1) return(NULL)
    message("Processing ", paste(path, collapse = " → "))
    df <- df %>%
      mutate(across(
        everything(),
        ~ replace_na(.x, 0)
      ))
    df <- df[rowSums(df) > 0, colSums(df) > 0]
    if (nrow(df) < 1 || ncol(df) < 1) {
      message("No data to process")
      return(NULL)
    }
    total    <- length(as.matrix(df))
    nonzero  <- sum(df != 0, na.rm = TRUE)
    sparsity <- (total - nonzero) / total * 100
    row_sparsity <- apply(df, 1, function(x) mean(x == 0, na.rm = TRUE) * 100)
    profiles      <- c(
      "motus25", 
      "metaphlan4.mpa_vJun23_CHOCOPhlAnSGB_202403"
    )
    training_sets <- c("metacardis", "galaxy")
    preds <- list()
    for (prof in profiles) {
      for (ts in training_sets) {
        message("Processing ", paste(c(path, prof, ts), collapse = " → "))
        model_path <- file.path(
          "inst/extdata/",
          ts,
          paste0("model.", prof, ".rds")
        )
        if (file.exists(model_path)) {
          model <- readRDS(model_path)
          referencetaxa <- names(model$trainingData)
          if (any(c("Shannon diversity", ".outcome") %in% referencetaxa)) {
            referencetaxa <- setdiff(referencetaxa, c("Shannon diversity", ".outcome"))
          }
        } else {
          warning("Missing model file: ", model_path)
          referencetaxa <- character(0)
        }
        colnames(df) <- colnames(df) %>% str_replace("unassigned", "-1")
        colnames(df) <- colnames(df) %>% str_replace("UNCLASSIFIED", "-1")
        samplesparsity <- compute_sparsity_per_sample(df, referencetaxa)
        if (mean(samplesparsity == 1, na.rm = TRUE) > 0.95) {
          print(paste("Model features for dataset are almost completely sparse:", mean(samplesparsity == 1, na.rm = TRUE)))
        }
        overlap_prop    <- prop_presenttaxa_overlap_per_sample(df, referencetaxa)
        nonoverlap_prop <- prop_taxa_notusedinmodel_per_sample(df, referencetaxa)
        psuedocount <- compute_psuedocount(df, referencetaxa)
        observed <- compute_observed(df)
        shannon <- compute_shannon(df)
        simpson <- compute_simpson(df)
        evenness <- compute_evenness(df, observed)
        df2 <- model_relative_abundances(df, referencetaxa)
        alpha_df <- compute_alpha_diversity(df2)
        rare_taxa <- evaluate_rare_taxa(df2)
        #ace <- compute_ace(df2)
        berger_parker <- compute_berger_parker(df2)
        inv_simpson <- compute_inv_simpson(df2)
        #chao1  <- compute_chao1(df2)
        rank_abundance_slope <- compute_rank_abundance_slope(df2)
        #pareto_alpha <- compute_pareto_alpha(df2)
        kurtosis <- compute_kurtosis(df2)
        skewness <- compute_skewness(df2)
        top_taxa <- get_top_taxa_per_sample(df2, top_n = 1)
        msg <- capture.output(
          tmp_res <- MLP(df, ts, "load", profiler = prof),
          type = "message"
        )
        res <- tmp_res
        if (!is.null(res)) {
          full_msg <- paste(msg, collapse = " ")
          m        <- regexec("([0-9]+)\\D*([0-9]+\\.?[0-9]*)\\s*%", full_msg)
          mm       <- regmatches(full_msg, m)[[1]]
          if (length(mm) >= 3) {
            res$model_says_sharedtaxa         <- as.integer(mm[2])
            res$model_says_sharedtaxa_percent <- as.numeric(mm[3])
          } else {
            res$model_says_sharedtaxa         <- NA_integer_
            res$model_says_sharedtaxa_percent <- NA_real_
          }
          res$dataset_sparsity_percent <- sparsity
          res$sample_sparsity_percent_alltaxa <- row_sparsity
          freqs <- table(res$load)
          mode_value <- as.numeric(names(freqs)[which.max(freqs)])
          mode_freq  <- max(freqs)
          mode_pct   <- mode_freq / length(res$load) * 100
          res$mode_of_load               <- mode_value
          res$modefrequencyofload        <- mode_freq
          res$modefrequencyofload_percent <- mode_pct
          res$profile      <- prof
          res$training_set <- ts
          res$log10_mlp    <- log10(res$load)
          res$log2_mlp     <- log2(res$load)
          res <- res %>% rename(mlp = load)
          res <- cbind(
            res,
            sample_sparsity_percent = samplesparsity*100,
            observed = observed,
            shannon = shannon,
            simpson = simpson,
            evenness = evenness,
            psuedocount = psuedocount,
            alpha_df,
            rare_taxa_percent = rare_taxa/alpha_df$Observed_modelfeatures * 100,
            taxa_used_in_model_percent    = overlap_prop*100,
            taxa_not_used_in_model_percent = nonoverlap_prop*100,
            #ace = ace,
            berger_parker = berger_parker,
            inv_simpson = inv_simpson,
            #chao1 = chao1,
            rank_abundance_slope = rank_abundance_slope,
            #pareto_alpha = pareto_alpha,
            kurtosis = kurtosis,
            skewness = skewness,
            top_taxa = top_taxa
          )
          preds[[length(preds) + 1L]] <- res
        }
      }
    }
    if (length(preds)) bind_rows(preds) else NULL
  }
  recurse <- function(node, path, runner) {
    if (is.data.frame(node)) return(runner(node, path))
    if (is.list(node)) {
      out <- lapply(names(node), function(nm)
        recurse(node[[nm]], c(path, nm), runner))
      names(out) <- names(node)
      out[lengths(out) == 0] <- NULL
      return(out)
    }
    NULL
  }
  results <- list()
  for (dn in names(datasets)) {
    seqtype <- tolower(datasets[[dn]]$sequencingtype %||% "")
    runner <- switch(
      seqtype,
      "amplicon"             = run_16s,
      "16s rrna"             = run_16s,
      "16s rna"              = run_16s,
      "16srna"               = run_16s,
      "16s"                  = run_16s,
      "metagenome"           = run_shotgun,
      "metagenomic"          = run_shotgun,
      "metagenomics"         = run_shotgun,
      "shotgun"              = run_shotgun,
      "shotgun metagenomics" = run_shotgun,
      stop("Unknown sequencing type: ", seqtype)
    )
    prop_root <- repo[[dn]]$proportions
    if (is.null(prop_root)) next
    preds <- list()
    for (ptype in names(prop_root)) {
      prop_data <- prop_root[[ptype]]
      if (is.null(prop_data)) next
      pr <- recurse(prop_data, c(dn, ptype), runner)
      if (!is.null(pr)) preds[[ptype]] <- pr
    }
    if (length(preds) == 0) next
    study <- repo[[dn]]
    study$studydemographics <- datasets[[dn]]
    study$predictions       <- preds
    results[[dn]]           <- study
  }
  return(results)
}

get_by_path <- function(study_list, path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  reduce(parts, function(obj, nm) {
    if (is.list(obj) && nm %in% names(obj)) obj[[nm]] else NULL
  }, .init = study_list)
}

pick_and_filter <- function(study_list, specs) {
  map(specs, function(spec) {
    df <- get_by_path(study_list, spec$path)
    if (is.null(df)) return(NULL)
    df_filt <- df %>%
      filter(profile == spec$profile,
             training_set == spec$training_set)
    if (nrow(df_filt) == 0) return(NULL)
    df_filt
  }) %>%
  discard(is.null) %>%
  set_names(map_chr(specs, ~ tail(strsplit(.x$path, "/", TRUE)[[1]], 1)))
}

pair_scale_to_predictions <- function(study, studyname) {
  get_dfs <- function(x) {
    if (is.data.frame(x)) {
      list(x)
    } else if (is.list(x)) {
      unlist(lapply(x, get_dfs), recursive = FALSE)
    } else {
      NULL
    }
  }
  scale_dfs <- get_dfs(study$scale)
  merge_one <- function(pred_df, path) {
    id_col <- "sample ID"
    if (!id_col %in% names(pred_df) || nrow(pred_df) < 1) return(NULL)
    pid <- as.character(pred_df[[id_col]])
    for (df in scale_dfs) {
      char_df <- df %>% mutate(across(everything(), as.character))
      potential_id_cols <- names(char_df)[grepl("sample|id|accession", tolower(names(char_df)))]
      for (col in potential_id_cols) {
        candidate_ids <- as.character(char_df[[col]])
        if (all(pid %in% candidate_ids)) {
          matched <- df %>% rename(!!id_col := all_of(col)) %>%
            mutate(!!id_col := as.character(!!sym(id_col)))
          rename_log_cols <- names(matched)[grepl("^(log10_|log2_)", names(matched)) & 
                                            !grepl("^(log10_mlp|log2_mlp)$", names(matched))]
          if (length(rename_log_cols) > 0) {
            new_names <- setNames(
              rename_log_cols,
              sapply(rename_log_cols, function(col) {
                prefix <- sub("^(log[0-9]+_).*", "\\1", col)
                rest <- sub("^log[0-9]+_", "", col)
                paste0(prefix, "truth_", rest)
              })
            )
            matched <- matched %>% rename(!!!new_names)
          }
          joined <- withCallingHandlers(
            left_join(pred_df, matched, by = id_col),
            warning = function(w) {
              message("Warning in process_mlp_predictions():\n  ", w$message)
            }
          )
          char_cols <- c(id_col)
          if ("Accession" %in% names(pred_df)) {
            char_cols <- c(char_cols, "Accession", "cohort")
          }
          char_cols <- unique(c(char_cols, 
            names(joined)[grepl("sample|id|accession|profile|training_set|cohort", tolower(names(joined)))]
          ))
          numeric_candidates <- setdiff(names(joined), char_cols)
          numeric_cols <- numeric_candidates[
            map_lgl(joined[numeric_candidates], function(x) all(grepl("^-?\\d*\\.?\\d*$", x)))
          ]
          joined <- joined %>%
            mutate(across(.cols = all_of(numeric_cols), ~ suppressWarnings(as.numeric(.))))
          return(joined)
        }
      }
    }
    message(sprintf("No matching scale columns found for predictions in study '%s' at path '%s'", 
                    studyname, paste(path, collapse = " → ")))
    NULL
  }
  merge_recursive <- function(x, path = character()) {
    if (is.data.frame(x)) {
      merge_one(x, path)
    } else if (is.list(x)) {
      out <- map2(x, names(x), function(val, nm) {
        merge_recursive(val, c(path, nm))
      })
      keep <- !map_lgl(out, is.null)
      out[keep]
    } else {
      x
    }
  }
  study$predictions <- merge_recursive(study$predictions)
  study
}

subsetNestedByCohort <- function(studyList, cohortValue) {
  if (is.null(studyList$metadata) || !"Sample" %in% colnames(studyList$metadata)) {
    stop("Expected a top‐level element named 'metadata' with a column 'Sample'.")
  }
  studyList$metadata <- studyList$metadata %>%
    filter(Cohort == cohortValue)
  sample_ids <- studyList$metadata$Sample
  recurse <- function(node, inCounts = FALSE, inProps = FALSE, inPreds = FALSE) {
    if (is.data.frame(node)) {
      if (inCounts || inProps) {
        keep_rows <- intersect(rownames(node), sample_ids)
        return(node[keep_rows, , drop = FALSE])
      }
      if (inPreds) {
        possible_cols <- intersect(c("sample ID", "sample.ID", "Sample"),
                                   colnames(node))
        if (length(possible_cols) == 0) {
          warning("Found a data.frame under 'predictions' that lacks any recognized sample‐ID column; leaving it unchanged.")
          return(node)
        }
        col_pred <- possible_cols[[1]]
        return(node %>% filter(.data[[col_pred]] %in% sample_ids))
      }
      return(node)
    }
    if (!is.list(node)) {
      return(node)
    }
    for (nm in names(node)) {
      if (identical(nm, "counts")) {
        node[[nm]] <- recurse(node[[nm]],
                              inCounts = TRUE,
                              inProps  = FALSE,
                              inPreds  = FALSE)
      } else if (identical(nm, "proportions")) {
        node[[nm]] <- recurse(node[[nm]],
                              inCounts = FALSE,
                              inProps  = TRUE,
                              inPreds  = FALSE)
      } else if (identical(nm, "predictions")) {
        node[[nm]] <- recurse(node[[nm]],
                              inCounts = FALSE,
                              inProps  = FALSE,
                              inPreds  = TRUE)
      } else if (nm %in% c("tax", "studydemographics")) {
        next
      } else {
        node[[nm]] <- recurse(node[[nm]],
                              inCounts = inCounts,
                              inProps  = inProps,
                              inPreds  = inPreds)
      }
    }
    return(node)
  }
  studyList <- recurse(studyList,
                       inCounts = FALSE,
                       inProps  = FALSE,
                       inPreds  = FALSE)
  return(studyList)
}

filter_study_by_specs <- function(study, specs_for_study) {
  if (is.character(specs_for_study)) {
    specs_for_study <- lapply(specs_for_study, function(p) list(path = p))
  }
  spec_paths <- vapply(specs_for_study, `[[`, "", "path")
  specs_tbl  <- data.frame(
    path         = spec_paths,
    profile      = vapply(specs_for_study, function(s) if (is.null(s$profile)) NA_character_ else s$profile, character(1)),
    training_set = I(lapply(specs_for_study, function(s) if (is.null(s$training_set)) character() else s$training_set)),
    stringsAsFactors = FALSE
  )
  cat("Normalized specs_tbl:\n"); print(specs_tbl)
  leaf_paths <- function(x, path = character()) {
    if (is.list(x) && !is.data.frame(x)) {
      unlist(lapply(names(x), function(nm) leaf_paths(x[[nm]], c(path, nm))))
    } else if (is.data.frame(x)) {
      paste(path, collapse = "/")
    } else {
      NULL
    }
  }
  is_subseq <- function(full_parts, spec_parts) {
    start <- 1L
    for (p in spec_parts) {
      idx <- which(full_parts[start:length(full_parts)] == p)
      if (length(idx) == 0) return(FALSE)
      start <- start + idx[1]
    }
    TRUE
  }
  prune_by_paths <- function(x, spec_paths, path = character()) {
    for (sp in spec_paths) {
      if (is_subseq(path, strsplit(sp, "/", fixed = TRUE)[[1]])) {
        return(x)
      }
    }
    if (is.data.frame(x)) return(NULL)
    if (is.list(x)) {
      out <- lapply(names(x), function(nm)
        prune_by_paths(x[[nm]], spec_paths, c(path, nm))
      )
      names(out) <- names(x)
      out <- out[!vapply(out, is.null, logical(1))]
      if (length(out)) return(out)
    }
    NULL
  }
  filter_predictions <- function(x, specs_tbl, path = character()) {
    if (is.list(x) && !is.data.frame(x)) {
      out <- lapply(names(x), function(nm)
        filter_predictions(x[[nm]], specs_tbl, c(path, nm))
      )
      names(out) <- names(x)
      out <- out[!vapply(out, is.null, logical(1))]
      if (length(out)) return(out)
      return(NULL)
    }
    if (is.data.frame(x)) {
      this_path <- paste(path, collapse = "/")
      row_i <- which(specs_tbl$path == this_path)
      if (length(row_i) != 1) {
        return(NULL)
      }
      df <- x
      prof <- specs_tbl$profile[row_i]
      if (!is.na(prof)) {
        df <- df[df$profile == prof, , drop = FALSE]
      }
      ts <- specs_tbl$training_set[[row_i]]
      if (length(ts) > 0) {
        df <- df[df$training_set %in% ts, , drop = FALSE]
      }
      if (nrow(df) == 0) return(NULL)
      return(df)
    }
    NULL
  }
  for (slot in c("counts", "proportions", "tax", "predictions")) {
    if (is.null(study[[slot]])) next
    cat(sprintf("\n-- SLOT: %s before prune:\n", slot))
    print(leaf_paths(study[[slot]]))
    pruned <- prune_by_paths(study[[slot]], spec_paths)
    if (slot == "predictions") {
      pruned <- filter_predictions(pruned, specs_tbl)
    }
    cat(sprintf("-- SLOT: %s after prune+filter:\n", slot))
    print(if (!is.null(pruned)) leaf_paths(pruned) else NULL)
    study[[slot]] <- pruned
  }
  empty_slot <- function(x) is.null(x) || (is.list(x) && length(x) == 0L)
  if (all(vapply(study[c("counts","proportions","tax","predictions")],
                 empty_slot, logical(1)))) {
    return(NULL)
  }
  study
}

filter_study_by_specs_global <- function(study, specs_for_study) {
  if (is.character(specs_for_study)) {
    specs_for_study <- lapply(specs_for_study, function(p) list(path = p))
  }
  spec_paths <- vapply(specs_for_study, `[[`, "", "path")
  specs_tbl  <- data.frame(
    path         = spec_paths,
    profile      = vapply(specs_for_study, function(s) if (is.null(s$profile)) NA_character_ else s$profile, character(1)),
    training_set = I(lapply(specs_for_study, function(s) if (is.null(s$training_set)) character() else s$training_set)),
    stringsAsFactors = FALSE
  )
  specs_tbl <- specs_tbl[!duplicated(specs_tbl$path), ]
  cat("Normalized specs_tbl:\n"); print(specs_tbl)
  leaf_paths <- function(x, path = character()) {
    if (is.list(x) && !is.data.frame(x)) {
      unlist(lapply(names(x), function(nm) leaf_paths(x[[nm]], c(path, nm))))
    } else if (is.data.frame(x)) {
      paste(path, collapse = "/")
    } else {
      NULL
    }
  }
  is_subseq <- function(full_parts, spec_parts) {
    start <- 1L
    for (p in spec_parts) {
      idx <- which(full_parts[start:length(full_parts)] == p)
      if (length(idx) == 0) return(FALSE)
      start <- start + idx[1]
    }
    TRUE
  }
  prune_by_paths <- function(x, spec_paths, path = character()) {
    for (sp in spec_paths) {
      if (is_subseq(path, strsplit(sp, "/", fixed = TRUE)[[1]])) {
        return(x)
      }
    }
    if (is.data.frame(x)) return(NULL)
    if (is.list(x)) {
      out <- lapply(names(x), function(nm)
        prune_by_paths(x[[nm]], spec_paths, c(path, nm))
      )
      names(out) <- names(x)
      out <- out[!vapply(out, is.null, logical(1))]
      if (length(out)) return(out)
    }
    NULL
  }
  filter_predictions <- function(x, specs_tbl, path = character()) {
    if (is.list(x) && !is.data.frame(x)) {
      out <- lapply(names(x), function(nm)
        filter_predictions(x[[nm]], specs_tbl, c(path, nm))
      )
      names(out) <- names(x)
      out <- out[!vapply(out, is.null, logical(1))]
      if (length(out)) return(out)
      return(NULL)
    }
    if (is.data.frame(x)) {
      this_path <- paste(path, collapse = "/")
      row_i <- which(specs_tbl$path == this_path)
      if (length(row_i) != 1) {
        return(NULL)
      }
      df <- x
      prof <- specs_tbl$profile[row_i]
      if (!is.na(prof)) {
        df <- df[df$profile == prof, , drop = FALSE]
      }
      ts <- specs_tbl$training_set[[row_i]]
      if (length(ts) > 0) {
        df <- df[df$training_set %in% ts, , drop = FALSE]
      }
      if (nrow(df) == 0) return(NULL)
      return(df)
    }
    NULL
  }
  for (slot in c("predictions")) {
    if (is.null(study[[slot]])) next

    cat(sprintf("\n-- SLOT: %s before prune:\n", slot))
    print(leaf_paths(study[[slot]]))

    pruned <- prune_by_paths(study[[slot]], spec_paths)

    if (slot == "predictions") {
      pruned <- filter_predictions(pruned, specs_tbl)
    }
    cat(sprintf("-- SLOT: %s after prune+filter:\n", slot))
    print(if (!is.null(pruned)) leaf_paths(pruned) else NULL)
    study[[slot]] <- pruned
  }
  empty_slot <- function(x) is.null(x) || (is.list(x) && length(x) == 0L)
  if (all(vapply(study[c("counts","proportions","tax","predictions")],
                 empty_slot, logical(1)))) {
    return(NULL)
  }
  study
}


compute_nested_r2<-function(predictionandtruthresults,exclude_loadtypes=NULL,exclude_truth_columns=NULL){

compute_r2<-function(df){
df<-df[,!grepl("sd|stdev",tolower(names(df)))]
log10_truth_cols<-grep("^log10_truth_",names(df),value=TRUE)
if(!is.null(exclude_truth_columns)){drop10<-paste0("log10_truth_",exclude_truth_columns);log10_truth_cols<-setdiff(log10_truth_cols,drop10)}
log2_truth_cols<-grep("^log2_truth_",names(df),value=TRUE)
if(!is.null(exclude_truth_columns)){drop2<-paste0("log2_truth_",exclude_truth_columns);log2_truth_cols<-setdiff(log2_truth_cols,drop2)}
do_r2<-function(pred_col,truth_col,data){
if(!(pred_col%in%names(data)&&truth_col%in%names(data)))return(NA_real_)
v<-data[,c(pred_col,truth_col),drop=FALSE]
ok<-complete.cases(v)&apply(v,1,function(z)all(is.finite(z)))
v<-v[ok,,drop=FALSE]
if(nrow(v)<2L)return(NA_real_)
pr<-v[[pred_col]]-mean(v[[pred_col]],na.rm=TRUE)
tr<-v[[truth_col]]-mean(v[[truth_col]],na.rm=TRUE)
sse<-sum((pr-tr)^2);sst<-sum(tr^2);if(sst==0)return(NA_real_)
1-sse/sst
}
do_r<-function(pred_col,truth_col,data){
if(!(pred_col%in%names(data)&&truth_col%in%names(data)))return(NA_real_)
v<-data[,c(pred_col,truth_col),drop=FALSE]
ok<-complete.cases(v)&apply(v,1,function(z)all(is.finite(z)))
v<-v[ok,,drop=FALSE]
if(nrow(v)<2L)return(NA_real_)
pr<-v[[pred_col]]-mean(v[[pred_col]],na.rm=TRUE)
tr<-v[[truth_col]]-mean(v[[truth_col]],na.rm=TRUE)
cor(pr,tr,use="complete.obs")
}
handle_one_group<-function(sub){
log10_r2<-setNames(vapply(log10_truth_cols,function(tr)do_r2("log10_mlp",tr,sub),numeric(1)),log10_truth_cols)
log2_r2<-setNames(vapply(log2_truth_cols,function(tr)do_r2("log2_mlp",tr,sub),numeric(1)),log2_truth_cols)
log10_r<-setNames(vapply(log10_truth_cols,function(tr)do_r("log10_mlp",tr,sub),numeric(1)),log10_truth_cols)
log2_r<-setNames(vapply(log2_truth_cols,function(tr)do_r("log2_mlp",tr,sub),numeric(1)),log2_truth_cols)
constant_cols<-names(sub)[vapply(sub,function(col)length(unique(col[!is.na(col)]))==1,logical(1))]
scalar_meta<-setNames(lapply(constant_cols,function(col)sub[[col]][1]),constant_cols)
c(list(log10_r2=log10_r2,log2_r2=log2_r2,log10_r=log10_r,log2_r=log2_r,profile=unique(sub$profile),training_set=unique(sub$training_set)),scalar_meta)
}
if(all(c("profile","training_set")%in%names(df))){
groups<-unique(df[,c("profile","training_set")],MARGIN=1)
results<-vector("list",nrow(groups))
names(results)<-apply(groups,1,paste,collapse="_")
for(i in seq_len(nrow(groups))){
sub<-df[((is.na(groups$profile[i])&is.na(df$profile))|df$profile==groups$profile[i])&((is.na(groups$training_set[i])&is.na(df$training_set))|df$training_set==groups$training_set[i]),,drop=FALSE]
results[[i]]<-handle_one_group(sub)
}
}else results<-list(overall=handle_one_group(df))
results
}
extract_dfs<-function(x,path=character()){
out<-list()
if(is.data.frame(x))out[[paste(path,collapse=">")]]<-x
else if(is.list(x))for(nm in names(x))out<-c(out,extract_dfs(x[[nm]],c(path,nm)))
out
}
safe_r2<-function(df){if(!is.data.frame(df))return(NULL);tryCatch(compute_r2(df),error=function(e)NULL)}
tidy_r2_rows<-function(r2_res,study,df_path,loadtypes){
rows<-list()
for(grp in names(r2_res)){
g<-r2_res[[grp]]
for(truth_col in names(g$log10_r2)){
truth_column<-sub("^log10_truth_","",truth_col)
loadtype_val<-if(length(loadtypes)==1L||all(is.na(loadtypes)))loadtypes[1]else{lc<-tolower(truth_column);lt_val<-NA_character_;for(lt in loadtypes){lt_l<-tolower(lt);if((lt_l=="flow cytometry"&&grepl("fc|flow|facs",lc))||(lt_l=="qpcr"&&grepl("qpcr|copy number|ct",lc))||(lt_l=="ddpcr"&&grepl("ddpcr",lc))||(lt_l=="mk_spike"&&grepl("mk_spike",lc))||(lt_l=="cfu"&&grepl("cfu",lc))){lt_val<-lt;break}};lt_val}
base<-list(study=study,dataframe=df_path,profile=g$profile,training_set=g$training_set,truth_column=truth_column,loadtype=loadtype_val,log10_r2=g$log10_r2[[truth_col]],log2_r2=g$log2_r2[[sub("^log10_","log2_",truth_col)]],log10_r=g$log10_r[[truth_col]],log2_r=g$log2_r[[sub("^log10_","log2_",truth_col)]])
meta<-g[setdiff(names(g),names(base))]
rows[[length(rows)+1]]<-as.data.frame(c(base,meta),stringsAsFactors=FALSE)
}
}
if(length(loadtypes)>1L&&length(rows)>0L){
used<-vapply(rows,`[[`,FUN.VALUE=character(1),"loadtype")
extra<-setdiff(loadtypes,unique(used))
if(length(extra)==1L)for(i in seq_along(rows))if(is.na(rows[[i]]$loadtype))rows[[i]]$loadtype<-extra
}
rows
}
study_results<-list()
for(study_name in names(predictionandtruthresults)){
study_obj<-predictionandtruthresults[[study_name]]
preds<-study_obj$predictions
loadtypes<-study_obj$studydemographics$loadtype%||%NA_character_
all_rows<-list()
dfs<-extract_dfs(preds)
for(df_path in names(dfs)){
r2<-safe_r2(dfs[[df_path]])
if(!is.null(r2))all_rows<-c(all_rows,tidy_r2_rows(r2,study_name,df_path,loadtypes))
}
if(length(all_rows)>0L)study_results[[study_name]]<-bind_rows(all_rows)
}
combined<-if(length(study_results))bind_rows(study_results)else NULL
if(!is.null(combined)){
if(!is.null(exclude_truth_columns))combined<-combined%>%filter(!truth_column%in%exclude_truth_columns)
if(!is.null(exclude_loadtypes))combined<-combined%>%filter(!loadtype%in%exclude_loadtypes)
}
list(per_study=study_results,combined=combined)
}

mode_plot <- function(var, thresh, ylab) {
  #r2results = r2results %>% filter(!(study == "GALAXY"),!(study == "MetaCardis"),!(study == "Vandeputte2017 Disease Cohort"), !(study == "Vandeputte2017 Study Cohort"))
  df2 <- r2results %>%
    distinct(study, dataframe, training_set, .keep_all = TRUE) %>%
    mutate(
      training_set = recode(
        training_set,
        galaxy     = "GALAXY",
        metacardis = "MetaCardis"
      )
    ) %>%
    group_by(training_set) %>%
    mutate(study = fct_reorder(study, .data[[var]])) %>%
    ungroup()
  theme_cell <- function(base_size = 18, base_family = "Arial") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        axis.line       = element_blank(),
        axis.line.x     = element_line(size = 1.2, color = "black"),
        axis.line.y     = element_line(size = 1.2, color = "black"),
        axis.ticks      = element_line(size = 1.2, color = "black"),
        axis.title.x    = element_blank(),
        axis.title.y    = element_text(size = base_size, face = "bold", margin=margin(0,10,0,0)),
        axis.text.x       = element_text(size = base_size, color = "black", angle=90, hjust=1, vjust=1),
        axis.text.y       = element_text(size = base_size, color = "black"),
        plot.title      = element_blank(),
        strip.text       = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        plot.margin        = margin(10, 10, 10, 10)
      )
  }
  if (var == "modefrequencyofload") {
    p <- ggplot(df2, aes(x = study, y = .data[[var]], fill = training_set)) +
      geom_col(color    = "black",
              width    = 0.7,
              position = position_identity()) +
      geom_text(
        data     = df2 %>% filter(.data[[var]] >= 5),
        aes(label = .data[[var]]),
        vjust    = -0.3,
        size     = 5,
        color    = "black",
        position = position_identity()
      ) +

      facet_grid(
        cols   = vars(training_set),
        scales = "free_x",
        space  = "free_x",
        switch = "x"              
      ) +
      scale_y_continuous(limits = c(0, 1000),          
                        expand = expansion(add = c(0, 20)),              
                        breaks = seq(0, 1000, by = 100)) +
      scale_fill_manual(values = c(GALAXY     = "#0C58CA",
                                  MetaCardis = "#FF8C27")) +
      theme(
        strip.placement = "outside",     
        strip.background = element_blank()
      )
  } else {
    p <- ggplot(df2, aes(x = training_set, y = .data[[var]])) +
      geom_boxplot(
        size  = 1.2,
        fill  = "white",
        color = "black",
        width = 0.6,
        staplewidth = 0.2,
        outlier.shape = NA
      ) +
      geom_line(
        aes(group = study),
        colour   = "grey70",
        linewidth = 0.8,
        alpha    = 0.6
      ) +
      geom_point(
        data = subset(df2, training_set != "Vandeputte2021"), 
        aes(fill = training_set),
        shape = 21,
        color = "black",
        size  = 3,
        alpha = 0.6
      ) +
      geom_jitter(
        data = subset(df2, training_set == "Vandeputte2021"), 
        aes(fill = training_set),
        shape = 21,
        color = "black",
        size  = 3,
        alpha = 0.6,
        width = 0.25,    
        height = 0       
      ) +
      scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, by = 10),
        expand = expansion(add = c(0, 2))
      ) +
      scale_x_discrete(
        drop   = FALSE,
        labels = c(
          "GALAXY"         = '<span style="color:#0C58CA">GALAXY</span>',
          "MetaCardis"     = '<span style="color:#FF8C27">MetaCardis</span>',
          "Vandeputte2021" = '<span style="color:#1B1B1B">Vandeputte2021</span>'
        )      
      ) + 
      scale_fill_manual(
        values = c(
          GALAXY         = "#0C58CA",
          MetaCardis     = "#FF8C27",
          Vandeputte2021 = "#1B1B1B"
        ),
        guide  = FALSE  
      )
  }
  p +
    labs(
      y     = ylab,
      title = str_to_title(gsub("_", " ", var))
    ) +
    theme_cell() +
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, face = "bold"))
}

compute_sample_residuals <- function(results,
                                     exclude_truth_columns = NULL,
                                     exclude_loadtypes     = NULL,
                                     pred_col              = "log10_mlp",
                                     truth_prefix          = "^log10_truth_",
                                     meta_vars             = c("sample_sparsity_percent","observed","shannon","simpson","evenness",
                                                               "Shannon_modelfeatures","Simpson_modelfeatures","Evenness_modelfeatures",
                                                               "Observed_modelfeatures","rare_taxa_percent","taxa_used_in_model_percent",
                                                               "taxa_not_used_in_model_percent","berger_parker","inv_simpson",
                                                               "rank_abundance_slope","kurtosis","skewness","top_taxa")) {
  extract_dfs <- function(x) {
    if (is.data.frame(x)) return(list(x))
    if (is.list(x))   return(flatten(map(x, extract_dfs)))
    list()
  }
  clean_training_set <- function(ts_vec) {
    ts_low <- stringr::str_to_lower(ts_vec)
    dplyr::case_when(
      str_detect(ts_low, "galaxy")     ~ "GALAXY",
      str_detect(ts_low, "metacardis") ~ "MetaCardis",
      TRUE                              ~ ts_vec
    )
  }
  full_df <- map_dfr(names(results), function(study_name) {
    study     <- results[[study_name]]
    dfs       <- extract_dfs(study$predictions)
    demos     <- study$studydemographics %||% list(loadtype = NA_character_)
    loadtypes <- demos$loadtype %||% NA_character_
    map_dfr(seq_along(dfs), function(i) {
      df_path <- names(dfs)[i] %||% paste0("df", i)
      df      <- as_tibble(dfs[[i]])
      truth_cols <- grep(truth_prefix, names(df), value = TRUE)
      map_dfr(truth_cols, function(truth_col) {
        if (!pred_col %in% names(df)) return(tibble())
        truth_column <- sub(truth_prefix, "", truth_col)
        lc           <- tolower(truth_column)
        loadtype_val <- NA_character_
        if (length(loadtypes) == 1L || all(is.na(loadtypes))) {
          loadtype_val <- loadtypes[1]
        } else {
          for (lt in loadtypes) {
            lt_l <- tolower(lt)
            if (lt_l == "flow cytometry" && grepl("fc|flow|facs", lc)) {
              loadtype_val <- lt; break
            }
            if (lt_l == "qpcr"            && grepl("qpcr|copy number|ct", lc)) {
              loadtype_val <- lt; break
            }
            if (lt_l == "ddpcr"           && grepl("ddpcr", lc)) {
              loadtype_val <- lt; break
            }
            if (lt_l == "mk_spike"        && grepl("mk_spike", lc)) {
              loadtype_val <- lt; break
            }
            if (lt_l == "cfu"             && grepl("cfu", lc)) {
              loadtype_val <- lt; break
            }
          }
        }
        subdf <- df %>%
          select(any_of(c("sample ID", pred_col, truth_col, meta_vars, "training_set", "profile"))) %>%
          filter(
            !is.na(.data[[pred_col]]),
            !is.na(.data[[truth_col]]),
            is.finite(.data[[pred_col]]),
            is.finite(.data[[truth_col]])
          )
        if (nrow(subdf) == 0) return(tibble())
        pred_s  <- mean_center(subdf[[pred_col]])
        truth_s <- mean_center(subdf[[truth_col]])
        ts_clean <- clean_training_set(subdf$training_set)
        tibble(
          study         = study_name,
          ID            = subdf$`sample ID`,
          dataframe     = df_path,
          truth_column  = truth_column,
          loadtype      = loadtype_val,
          profile       = subdf$profile,
          training_set  = ts_clean,
          pred_scaled   = pred_s,
          truth_scaled  = truth_s,
          residual      = abs(pred_s - truth_s),
          !!!subdf[meta_vars]
        )
      })
    })
  })
  if (!is.null(exclude_truth_columns)) {
    full_df <- full_df %>%
      filter(!truth_column %in% exclude_truth_columns)
  }
  if (!is.null(exclude_loadtypes)) {
    full_df <- full_df %>%
      filter(!loadtype %in% exclude_loadtypes)
  }
  full_df
}


subset_r2 <- function(r2_out, specs) {
  per_study <- map(specs, function(study_specs, study_name) {
    study_df <- r2_out$per_study[[study_name]]
    if (is.null(study_df)) return(NULL)
    map(study_specs, function(spec) {
      df <- study_df %>%
        filter(
          dataframe    == spec$path,
          profile      == spec$profile,
          (is.null(spec$training_set) | training_set %in% spec$training_set)
        )
      if (!is.null(spec$loadtype)) {
        df <- df %>% filter(loadtype %in% spec$loadtype)
      }
      if (!is.null(spec$load_columns)) {
        df <- df %>% select(all_of(c("study","dataframe","profile",
                                     "training_set","truth_column",
                                     spec$load_columns)))
      }
      df
    })
  }, .id = "study")
  combined <- map2_dfr(specs, names(specs), function(study_specs, study_name) {
    map_dfr(study_specs, function(spec) {
      df <- r2_out$combined %>%
        filter(
          study        == study_name,
          dataframe    == spec$path,
          profile      == spec$profile,
          (is.null(spec$training_set) | training_set %in% spec$training_set)
        )
      if (!is.null(spec$loadtype)) {
        df <- df %>% filter(loadtype %in% spec$loadtype)
      }
      df
    })
  })
  list(per_study = per_study, combined = combined)
}

remove_empty_named_lists <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    x <- lapply(x, remove_empty_named_lists)
    x <- x[!(vapply(x, function(e) is.list(e) && length(e) == 0, logical(1)))]
  }
  x
}

prop_presenttaxa_overlap_per_sample <- function(relabun_df, referencetaxa) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a data.frame or matrix.")
  }
  relabun_df <- as.data.frame(relabun_df, stringsAsFactors = FALSE)
  missing_taxa <- setdiff(referencetaxa, colnames(relabun_df))
  if (length(missing_taxa) > 0) {
    relabun_df[missing_taxa] <- 0
  }
  df_ref <- relabun_df[, referencetaxa, drop = FALSE]
  n_ref <- ncol(df_ref)
  if (n_ref == 0) {
    warning("No reference taxa found—returning NA for every sample.")
    return(rep(NA_real_, nrow(df_ref)))
  }
  present_counts <- rowSums(df_ref > 0, na.rm = TRUE)
  prop_present   <- present_counts / n_ref
  return(prop_present)
}

prop_taxa_notusedinmodel_per_sample <- function(relabun_df, referencetaxa) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a data.frame or matrix.")
  }
  relabun_df <- as.data.frame(relabun_df, stringsAsFactors = FALSE)
  total_cols <- ncol(relabun_df)
  if (total_cols == 0) {
    warning("Input has zero columns—returning NA for every sample.")
    return(rep(NA_real_, nrow(relabun_df)))
  }
  non_model_taxa <- setdiff(colnames(relabun_df), referencetaxa)
  if (length(non_model_taxa) == 0) {
    return(rep(0, nrow(relabun_df)))
  }
  non_counts <- rowSums(relabun_df[, non_model_taxa, drop = FALSE] > 0, na.rm = TRUE)
  prop_non   <- non_counts / total_cols
  return(prop_non)
}

compute_psuedocount <- function(relabun_df, referencetaxa) {
  overlapping_taxa <- intersect(colnames(relabun_df), referencetaxa)
  relabun_df <- relabun_df[, overlapping_taxa, drop = FALSE]
  nonzero_vals <- unlist(relabun_df)
  nonzero_vals <- nonzero_vals[nonzero_vals > 0]
  if (length(nonzero_vals) == 0) {
    return(0)
  }
  input.min <- min(nonzero_vals) / 2
  return(input.min)
}

model_relative_abundances <- function(relabun_df, reference_taxa) {
  overlapping_taxa <- intersect(colnames(relabun_df), reference_taxa)
  missing_taxa <- setdiff(reference_taxa, colnames(relabun_df))
  relabun_df <- relabun_df[, overlapping_taxa, drop = FALSE]
  if (length(missing_taxa) > 0) {
    zeros_mat <- matrix(0, nrow = nrow(relabun_df), ncol = length(missing_taxa))
    colnames(zeros_mat) <- missing_taxa
    rownames(zeros_mat) <- rownames(relabun_df)
    relabun_df <- cbind(relabun_df, zeros_mat)
  }
  relabun_df <- relabun_df[, reference_taxa, drop = FALSE]
  nonzero_vals <- unlist(relabun_df)
  nonzero_vals <- nonzero_vals[nonzero_vals > 0]
  if (length(nonzero_vals) == 0) {
    relabun_df <- as.data.frame(matrix(
      0,
      nrow = nrow(relabun_df),
      ncol = ncol(relabun_df),
      dimnames = list(rownames(relabun_df), colnames(relabun_df))
    ))
  } else {
    input.min <- min(nonzero_vals) / 2
    relabun_df[relabun_df == 0] <- input.min
  }
  return(relabun_df)
}

compute_observed <- function(relabun_df) {
  rowSums(relabun_df > 0)
}

compute_shannon <- function(relabun_df) {
  vegan::diversity(relabun_df, index = "shannon")
}

compute_simpson <- function(relabun_df) {
  vegan::diversity(relabun_df, index = "simpson")
}

compute_evenness <- function(relabun_df, observed) {
  vegan::diversity(relabun_df, index = "shannon") / log(observed)
}

compute_alpha_diversity <- function(relabun_df) {
  Observed <- rowSums(relabun_df > 0)
  Shannon <- diversity(relabun_df, index = "shannon")
  Simpson <- diversity(relabun_df, index = "simpson")
  Evenness <- Shannon / log(Observed)
  alpha_df <- data.frame(
    Observed_modelfeatures = Observed,
    Shannon_modelfeatures  = Shannon,
    Simpson_modelfeatures  = Simpson,
    Evenness_modelfeatures = Evenness*100,
    stringsAsFactors = FALSE
  )
  return(alpha_df)
}

get_top_taxa_per_sample <- function(relabun_df, top_n = 1) {
  out <- apply(relabun_df, 1, function(x) {
    names(sort(x, decreasing = TRUE))[1:top_n]
  })
  setNames(out, rownames(relabun_df))
}

compute_chao1 <- function(relabun_df) {
  out <- sapply(1:nrow(relabun_df), function(i) {
    est <- estimateR(relabun_df[i, ] * sum(relabun_df[i, ]))
    est["S.chao1"]
  })
  setNames(out, rownames(relabun_df))
}

compute_ace <- function(relabun_df) {
  out <- sapply(1:nrow(relabun_df), function(i) {
    est <- estimateR(relabun_df[i, ] * sum(relabun_df[i, ]))
    est["ace"]
  })
  setNames(out, rownames(relabun_df))
}

compute_berger_parker <- function(relabun_df) {
  out <- apply(relabun_df, 1, max)
  setNames(out, rownames(relabun_df))
}

compute_inv_simpson <- function(relabun_df) {
  out <- apply(relabun_df, 1, function(x) 1 / sum(x^2))
  setNames(out, rownames(relabun_df))
}

compute_skewness <- function(relabun_df) {
  out <- apply(relabun_df, 1, e1071::skewness)
  setNames(out, rownames(relabun_df))
}

compute_kurtosis <- function(relabun_df) {
  out <- apply(relabun_df, 1, e1071::kurtosis)
  setNames(out, rownames(relabun_df))
}

compute_rank_abundance_slope <- function(relabun_df) {
  out <- sapply(1:nrow(relabun_df), function(i) {
    abund <- unlist(relabun_df[i, ])
    abund <- sort(abund[abund > 0], decreasing = TRUE)
    if (length(abund) < 2) return(NA_real_)
    ranks <- seq_along(abund)
    coef(lm(log10(abund) ~ log10(ranks)))[2] 
  })
  setNames(out, rownames(relabun_df))
}

compute_pareto_alpha <- function(relabun_df, tail_prop = 0.1) {
  library(poweRlaw)
  if (!requireNamespace("poweRlaw", quietly = TRUE)) {
    stop("The 'poweRlaw' package is required. Install it via install.packages('poweRlaw').")
  }
  out <- sapply(1:nrow(relabun_df), function(i) {
    x <- unlist(relabun_df[i, ])
    x <- sort(x[x > 0], decreasing = TRUE)
    if (length(x) < 5) return(NA_real_)
    cutoff <- quantile(x, 1 - tail_prop)
    tail_x <- x[x >= cutoff]
    if (length(tail_x) < 2) return(NA_real_)
    pl <- poweRlaw::conpl$new(tail_x) 
    est_xmin <- poweRlaw::estimate_xmin(pl)
    pl$setXmin(est_xmin)
    est_pars <- poweRlaw::estimate_pars(pl)
    return(est_pars$pars) 
  })
  setNames(out, rownames(relabun_df))
}

evaluate_rare_taxa <- function(relabun_df, threshold = 0.01) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a dataframe or matrix.")
  }
  rare_taxa_count <- apply(relabun_df, 1, function(row) {
    sum(row < threshold & row > 0)  
  })
  return(rare_taxa_count)
}

compute_sparsity_per_sample <- function(relabun_df, referencetaxa) {
  if (!is.data.frame(relabun_df) && !is.matrix(relabun_df)) {
    stop("Input must be a dataframe or matrix.")
  }
  relabun_df <- as.data.frame(relabun_df, stringsAsFactors = FALSE)
  missing_taxa <- setdiff(referencetaxa, colnames(relabun_df))
  if (length(missing_taxa) > 0) {
    relabun_df[missing_taxa] <- 0
  }
  relabun_df <- relabun_df[, referencetaxa, drop = FALSE]
  sparsity <- apply(relabun_df, 1, function(row) {
    mean(row == 0, na.rm = TRUE)
  })
  return(sparsity)
}

title_model <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, regex("galaxy",     ignore_case = TRUE)) ~ "GALAXY",
    stringr::str_detect(x, regex("metacardis", ignore_case = TRUE)) ~ "MetaCardis",
    stringr::str_detect(x, regex("vandeputte", ignore_case = TRUE)) ~ "Vandeputte2021",
    TRUE                                                             ~ stringr::str_to_title(x)
  )
}

plot_features <- function(data,x_var,y_var, colour_var = NULL, label_thresh = 0, save_as = NULL,
                          width  = 6,height = 4, dpi = 600) {
  if (is.character(x_var)) x_sym <- sym(x_var)     else x_sym <- ensym(x_var)
  if (is.character(y_var)) y_sym <- sym(y_var)     else y_sym <- ensym(y_var)
  if (!is.null(colour_var)) {
    if (is.character(colour_var)) col_sym <- sym(colour_var)
    else                           col_sym <- ensym(colour_var)
    col_nm <- as_name(col_sym)
  } else {
    col_sym <- NULL; col_nm <- NULL
  }
  x_nm <- as_name(x_sym)
  y_nm <- as_name(y_sym)
  if ("training_set" %in% names(data) && is.character(data$training_set)) {
    data <- data %>% mutate(training_set = title_model(training_set))
  }
  df <- data %>%
    filter(
      !is.na(!!x_sym), !is.na(!!y_sym),
      is.finite(!!x_sym), is.finite(!!y_sym)
    )
  if (is.null(col_sym)) {
    ann_df <- df %>%
      summarise(
        fit   = list(lm(!!y_sym ~ !!x_sym, data = .)),
        x_pos = max(!!x_sym, na.rm = TRUE),
        y_pos = min(!!y_sym, na.rm = TRUE)
      ) %>%
      mutate(
        R2    = map_dbl(fit, ~ glance(.x)$r.squared),
        label = paste0("R^2 == ", formatC(R2, digits = 2, format = "f"))
      )
  } else {
    ann_df <- df %>%
      group_by(!!col_sym) %>%
      summarise(
        fit   = list(lm(!!y_sym ~ !!x_sym, data = .)),
        x_pos = max(!!x_sym, na.rm = TRUE),
        y_pos = min(!!y_sym, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        R2    = map_dbl(fit, ~ glance(.x)$r.squared),
        label = paste0("R^2 == ", formatC(R2, digits = 2, format = "f"))
      )
  }
  study_sym <- sym("study")
  lbl_df <- df %>% filter(!!y_sym > label_thresh)
  aes_map <- if (is.null(col_sym)) {
    aes(x = !!x_sym, y = !!y_sym)
  } else {
    aes(x = !!x_sym, y = !!y_sym, colour = !!col_sym)
  }
  p <- ggplot(df, aes_map) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(
      method    = "lm",
      formula   = y ~ x,
      se        = TRUE,
      level     = 0.95,
      linewidth = 0.8
    ) +
    geom_text(
      data       = ann_df,
      aes(x = x_pos, y = y_pos, label = label,
          colour = if (!is.null(col_sym)) !!col_sym else NULL),
      parse      = TRUE,
      hjust      = 1, 
      vjust      = 0,
      size       = 4,
      show.legend= FALSE
    ) +
    geom_text_repel(
      data       = lbl_df,
      aes(x = !!x_sym, y = !!y_sym, label = !!study_sym,
          colour = if (!is.null(col_sym)) !!col_sym else NULL),
      size       = 3,
      min.segment.length = 0.1,
      show.legend= FALSE
    ) +
    labs(
      x      = str_to_title(gsub("_", " ", x_nm)),
      y      = str_to_title(gsub("_", " ", y_nm)),
      colour = if (!is.null(col_sym))
                 str_to_title(gsub("_", " ", col_nm))
               else NULL
    ) +
    { if (!is.null(col_sym) && col_nm == "training_set")
        scale_colour_manual(values = c(
          Vandeputte2021 = "#1B1B1B",
          GALAXY         = "#0C58CA",
          MetaCardis     = "#FF8C27"
        ))
    } +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 14),
      axis.text        = element_text(size = 12),
      legend.position  = "right"
    )
  if (!is.null(save_as)) {
    device_to_use <- if (grepl("\\.png$", save_as, ignore.case = TRUE))
                       "png" else cairo_pdf
    ggsave(
      filename = save_as,
      plot     = p,
      width    = width,
      height   = height,
      units    = "in",
      dpi      = dpi,
      device   = device_to_use
    )
  }
  invisible(p)
}

choose_truth_col <- function(df) {
  truth_cols <- grep("^log10_truth_", names(df), value = TRUE)
  valid     <- truth_cols[vapply(truth_cols, function(col) any(is.finite(df[[col]])), logical(1))]
  if (length(valid) == 0) return(NA_character_)
  valid[[1]]
}

flatten_repo_to_long <- function(repo) {
  rows <- list(); i <- 1L
  recurse <- function(node, path = character()) {
    if (inherits(node, "data.frame")) {
      df           <- node
      study        <- path[1]
      dataframe_nm <- path[length(path)]
      prof         <- df$profile
      trn          <- df$training_set
      if ("log10_mlp" %in% names(df)) {
        rows[[i]] <<- tibble(
          Dataset     = study,
          DataFrame   = dataframe_nm,
          Profile     = prof,
          TrainingSet = trn,
          Source      = "Predicted",
          Label       = "Predicted",
          Load        = df$log10_mlp
        )
        i <<- i + 1L
      }
      truth_cols <- grep("^log10_truth_", names(df), value = TRUE)
      if (length(truth_cols)) {
        rows[[i]] <<- df %>%
          select(all_of(truth_cols)) %>%
          pivot_longer(everything(),
                       names_to  = "TruthCol",
                       values_to = "Load") %>%
          filter(is.finite(Load)) %>%
          mutate(
            Dataset     = study,
            DataFrame   = dataframe_nm,
            Profile     = "Measured",
            TrainingSet = study,
            Source      = "Measured",
            Label       = "Measured"
          ) %>%
          select(Dataset, DataFrame, Profile, TrainingSet, Source, Label, Load)
        i <<- i + 1L
      }
    } else if (is.list(node)) {
      nm <- names(node)
      for (j in seq_along(node)) {
        name_j <- if (!is.null(nm) && nzchar(nm[j])) nm[j] else paste0("[[", j, "]]")
        recurse(node[[j]], c(path, name_j))
      }
    }
  }
  recurse(repo)
  bind_rows(rows)
}

output_structure <- function(pred_and_demo_list) {
  inspect_one <- function(output, study_name = NULL) {
    if (!is.list(output)) {
      stop(sprintf(
        "%s: Output must be a list",
        ifelse(is.null(study_name), "", study_name)
      ))
    }
    get_structure <- function(obj) {
      tryCatch({
        if (is.null(obj)) {
          "NULL"
        } else if (is.matrix(obj)) {
          sprintf("MATRIX(%d x %d)", nrow(obj), ncol(obj))
        } else if (is.data.frame(obj)) {
          sprintf("DATAFRAME(%d x %d)", nrow(obj), ncol(obj))
        } else if (is.list(obj)) {
          if (length(obj) == 0) return("EMPTY_LIST")
          nested <- sapply(names(obj), function(nm) {
            paste(nm, "=", get_structure(obj[[nm]]))
          })
          sprintf("LIST={%s}", paste(nested, collapse = ", "))
        } else if (is.logical(obj) && all(is.na(obj))) {
          "NULL"
        } else if (is.vector(obj)) {
          sprintf("VECTOR(length=%d)", length(obj))
        } else {
          class(obj)[1]
        }
      }, error = function(e) {
        sprintf("ERROR: %s", e$message)
      })
    }
    result <- map(
      names(output),
      function(elem_name) {
        get_structure(output[[elem_name]])
      }
    )
    names(result) <- names(output)
    return(result)
  }
  all_structures <- map(
    names(pred_and_demo_list),
    function(study_name) {
      inspect_one(pred_and_demo_list[[study_name]], study_name = study_name)
    }
  )
  names(all_structures) <- names(pred_and_demo_list)
  return(print(all_structures))
}

get_label <- function(ts) {
  if (str_detect(ts, regex("vandeputte",   TRUE))) {
    "Vandeputte2021"
  } else if (str_detect(ts, regex("galaxy",    TRUE))) {
    "GALAXY"
  } else if (str_detect(ts, regex("metacardis",TRUE))) {
    "MetaCardis"
  } else {
    ts
  }
}

z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.finite(s) && s > 0) (x - mean(x, na.rm = TRUE)) / s else rep(0, length(x))
}

mean_center <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.finite(m)) x - m else rep(0, length(x))
}


write_all_load_plots <- function(predictions_only,
                                 output_dir = "plots",
                                 source_palette = c(
                                   Vandeputte2021 = "#000000",
                                   GALAXY         = "#0C58CA",
                                   MetaCardis     = "#FF8C27"
                                 )) {

  tidyplot_theme <- function(base_size = 18) {
    theme_minimal(base_size = base_size) +
      theme(
        text              = element_text(family = "Helvetica", color = "black"),
        plot.title        = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        axis.title.x      = element_text(size = 16, face = "bold", margin = margin(t = 10), hjust = 0.5),
        axis.title.y      = element_markdown(size = 16, face = "bold", margin = margin(r = 10), hjust = 0.5),
        axis.title        = element_text(size = 14, face = "plain"),
        axis.text         = element_text(size = 12, color = "black"),
        axis.line         = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks        = element_line(linewidth = 0.5, color = "black"),
        panel.grid.major  = element_line(color = "grey90", size = 0.1),
        panel.grid.minor  = element_line(color = "grey95", size = 0.05),
        panel.border      = element_blank(),
        legend.position   = "none",
        plot.background   = element_rect(fill = "white", color = NA),
        panel.background  = element_rect(fill = "white", color = NA)
      )
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  flatten <- function(x, path = "") {
    dfs <- list()
    if (is.list(x)) {
      for (nm in names(x)) {
        el <- x[[nm]]
        new_path <- if (path == "") nm else paste(path, nm, sep = "/")
        if (is.data.frame(el)) {
          dfs[[new_path]] <- el
        } else if (is.list(el)) {
          dfs <- c(dfs, flatten(el, new_path))
        }
      }
    }
    dfs
  }

  for (study_name in names(predictions_only)) {
    study_entry <- predictions_only[[study_name]]
    id_col      <- "sample ID"
    clean_study_title <- gsub("\\s*\\(.*?\\)\\s*", " ", study_name) %>% str_trim()
    safe_study        <- gsub("[^[:alnum:]]+", "_", clean_study_title)
    dfs_list <- flatten(study_entry)

    for (df_path in names(dfs_list)) {
      df <- dfs_list[[df_path]]
      if (!all(c(id_col, "log10_mlp", "profile", "training_set") %in% names(df))) next
      truth_cols <- grep("^log10_truth_", names(df), value = TRUE)
      if (length(truth_cols) == 0) next

      for (truth_col in truth_cols) {
        df_pair <- df %>%
          select(
            ID        = !!sym(id_col),
            Predicted = log10_mlp,
            Truth     = !!sym(truth_col),
            profile,
            training_set
          ) %>%
          filter(is.finite(Truth), is.finite(Predicted)) %>%
          mutate(Source = vapply(training_set, get_label, FUN.VALUE = character(1)))

        sample_n <- nrow(df_pair)
        if (sample_n < 3 || var(df_pair$Truth) == 0) next

        df_pair <- df_pair %>%
          group_by(Source) %>%
          mutate(
            Pred_z  = mean_center(Predicted),
            Truth_z = mean_center(Truth)
          ) %>%
          ungroup()

        unique_sources <- unique(df_pair$Source)

        if (length(unique_sources) == 1) {

          src     <- unique_sources
          pearson <- cor(df_pair$Pred_z, df_pair$Truth_z, use = "complete.obs")
          r2 <- {
            sse <- sum((df_pair$Pred_z - df_pair$Truth_z)^2)
            sst <- sum((df_pair$Truth_z)^2)
            round(1 - sse / sst, 2)
          }
          if (grepl("Vandeputte2021", src, ignore.case = TRUE)) {
            y_label <- "**<span style=\"color:#1B1B1B\">Vandeputte2021</span>**"
          } else if (grepl("GALAXY", src, ignore.case = TRUE)) {
            y_label <- "**<span style=\"color:#0C58CA\">GALAXY</span>**"
           } else if (grepl("MetaCardis", src, ignore.case = TRUE)) {
            y_label <- "**<span style=\"color:#FF8C27\">MetaCardis</span>**"
          }

          p <- ggplot(df_pair, aes(x = Truth_z, y = Pred_z, color = Source, fill = Source)) +

            geom_point(size = 1, alpha = 0.7, shape = 21) +

            geom_smooth(aes(fill = Source),
                        method   = "lm",
                        formula  = y ~ x,
                        se       = TRUE,
                        colour   = NA) +

            geom_smooth(aes(x = Truth_z, y = Pred_z),
                        inherit.aes = FALSE,
                        method      = "lm",
                        formula     = y ~ x,
                        se          = FALSE,
                        colour      = "white",
                        linewidth   = 1.8) +

            geom_smooth(aes(color = Source),
                        method    = "lm",
                        formula   = y ~ x,
                        se        = FALSE,
                        linewidth = 1) +
            labs(
              title = sprintf("%s (n = %d)", clean_study_title, sample_n),
              x     = "Measured Total Microbial Load",
              y     = y_label
            ) +
            scale_color_manual(values = source_palette) +
            scale_fill_manual(values = source_palette) +
            tidyplot_theme() +
            coord_cartesian(clip = "off") +
            scale_x_continuous(minor_breaks = waiver()) +
            scale_y_continuous(minor_breaks = waiver()) +
            geom_richtext(
              data = data.frame(),
              aes(
                x = Inf, y = -Inf,
                label = paste0(
                  "<span style='color:", source_palette[src], "'>", src, "</span> ",
                  "<span style='color:black'>ρ: ", sprintf("%.2f", pearson), ", R²: ", sprintf("%.2f", r2), "</span>"
                )
              ),
              inherit.aes = FALSE,
              hjust = 1, vjust = 0,
              size = 4,
              fill = NA, label.color = NA
            )

        } else {
          y_label <- "**<span style=\"color:#0C58CA\">GALAXY</span> / <span style=\"color:#FF8C27\">MetaCardis</span>**"

          stats_df <- df_pair %>%
            group_by(Source) %>%
            summarize(
              pearson = cor(Pred_z, Truth_z, use = "complete.obs"),
              r2 = {
                sse <- sum((Pred_z - Truth_z)^2)
                sst <- sum((Truth_z)^2)
                round(1 - sse / sst, 2)
              },
              .groups = "drop"
            ) %>%
            arrange(match(Source, names(source_palette))) %>%        
            mutate(label = paste0(
              "<span style='color:", source_palette[Source], "'>", Source, "</span> ",
              "<span style='color:black'>ρ: ", sprintf('%.2f', pearson),
              ", R²: ", sprintf('%.2f', r2), "</span>"
            ))

          combined_label <- paste(stats_df$label, collapse = "<br>")

          p <- ggplot(df_pair, aes(x = Truth_z, y = Pred_z, color = Source, fill = Source)) +
            geom_point(size = 1, alpha = 0.7, shape = 21) +
            geom_smooth(aes(fill = Source),
                        method   = "lm", formula = y ~ x,
                        se       = TRUE,  colour  = NA) +
            geom_smooth(aes(x = Truth_z, y = Pred_z),
                        inherit.aes = TRUE,
                        method      = "lm", formula = y ~ x,
                        se          = FALSE, colour = "white", linewidth = 1.8) +
            geom_smooth(aes(color = Source),
                        method    = "lm", formula = y ~ x,
                        se        = FALSE, linewidth = 1) +
            labs(
              title = sprintf("%s (n = %d)", clean_study_title, sample_n),
              x     = "Measured Total Microbial Load",
              y     = y_label
            ) +
            scale_color_manual(values = source_palette) +
            scale_fill_manual(values = source_palette) +
            tidyplot_theme() +
            coord_cartesian(clip = "off") +
            scale_x_continuous(minor_breaks = waiver()) +
            scale_y_continuous(minor_breaks = waiver()) +
            ggtext::geom_richtext(
              data        = data.frame(x = Inf, y = -Inf, label = combined_label),
              aes(x, y, label = label),
              inherit.aes = FALSE,
              hjust       = 1,
              vjust       = 0,
              size        = 4,
              fill        = NA,
              label.color = NA
            )
        }

        safe_truth <- sub("^log10_truth_", "", truth_col)
        safe_df    <- gsub("[^[:alnum:]]+", "_", df_path)
        fname_base <- paste(safe_study, safe_df, safe_truth, sep = "_")

        ggsave(
          file.path(output_dir, paste0(fname_base, ".pdf")),
          plot   = p,
          device = cairo_pdf,
          width  = 8, height = 8, units = "in", dpi = 600, bg = "white"
        )
        ggsave(
          file.path(output_dir, paste0(fname_base, ".png")),
          plot  = p,
          width = 5, height = 5, units = "in", dpi = 600, bg = "white"
        )
      }
    }
  }
}

clean_nested_dataframes <- function(obj, exact_names, substr_patterns) {
  if (is.data.frame(obj)) {
    if (!is.null(exact_names)) {
      obj <- obj[, !colnames(obj) %in% exact_names, drop = FALSE]
    }
    if (!is.null(substr_patterns)) {
      matches <- sapply(substr_patterns, function(p)
        grepl(p, colnames(obj), ignore.case = TRUE))
      obj <- obj[, !rowSums(matches) > 0, drop = FALSE]
    }
    return(obj)
  } else if (is.list(obj)) {
    return(lapply(obj, clean_nested_dataframes, exact_names, substr_patterns))
  } else {
    return(obj)
  }
}

is_empty_counts <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.atomic(x)) return(all(is.na(x)))
  if (is.data.frame(x)) {
    return(nrow(x) == 0 || all(sapply(x, function(col) all(is.na(col)))))
  }
  if (is.list(x)) {
    return(length(x) == 0 || all(sapply(x, is_empty_counts)))
  }
  return(FALSE)
}
prune_empty_lists <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.list(x)) {
    return(x)
  }
  x <- lapply(x, prune_empty_lists)
  keep_idx <- !sapply(x, function(el) {
    is.null(el) ||
    (is.list(el) && length(el) == 0) ||
    (is.atomic(el) && all(is.na(el))) ||
    (is.list(el) && all(sapply(el, function(y) {
       is.null(y) ||
       (is.atomic(y) && all(is.na(y))) ||
       (is.list(y) && length(y) == 0)
    })))
  })
  x[keep_idx]
}

extract_leaf_df <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && length(x) > 0 && all(sapply(x, is.atomic))) {
    return(as.data.frame(x))
  }
  if (is.list(x)) {
    for (el in x) {
      df <- extract_leaf_df(el)
      if (!is.null(df)) return(df)
    }
  }
  return(NULL)
}
# compute all metrics (Pearson r, R², SMAPE, RMSE, MAE, Bias, MAPE, MDAPE, Spearman's p, CCC, CV) 
compute_all_metrics <- function(pred, truth) {
  idx <- which(is.finite(pred) & is.finite(truth))
  pred <- mean_center(pred[idx])
  truth <- mean_center(truth[idx])
  n <- length(pred)
  if (n < 3 || var(truth) == 0) return(NULL)
  resid <- pred - truth
  pearson_r <- cor(pred, truth, method = "pearson", use = "complete.obs")
  r2 <- {
    sse <- sum(resid^2)
    sst <- sum((truth)^2)
    1 - sse / sst
  }
  explained_variance <- if (var(truth) > 0) 1 - var(resid) / var(truth) else NA_real_
  rmse <- sqrt(mean(resid^2))
  mae <- mean(abs(resid))
  medae <- median(abs(resid))
  bias <- mean(resid)
  nz <- which(truth != 0)
  if (length(nz) > 0) {
    smape <- 200 * mean(abs(pred[idx] - truth[idx]) / (abs(pred[idx]) + abs(truth[idx])))
    mape  <- 100 * mean(abs((pred[nz] - truth[nz]) / truth[nz]))
    mdape <- 100 * median(abs((pred[nz] - truth[nz]) / truth[nz]))
  } else {
    mape  <- NA_real_
    mdape <- NA_real_
  }
  spearman_p <- cor(pred, truth, method = "spearman", use = "complete.obs")
  mu_p   <- mean(pred)
  mu_t   <- mean(truth)
  var_p  <- var(pred)
  var_t  <- var(truth)
  cov_pt <- cov(pred, truth)
  denom_ccc <- var_p + var_t + (mu_p - mu_t)^2
  ccc_val <- if (denom_ccc > 0) 2 * cov_pt / denom_ccc else NA_real_
  cv_val <- if (mean(truth) != 0) 100 * sd(resid) / mean(truth) else NA_real_
  tibble(
    Pearson_r  = pearson_r,
    R2         = r2,
    Explained_Variance = explained_variance,
    SMAPE      = smape,
    RMSE       = rmse,
    MAE        = mae,
    MedAE      = medae,
    N          = sum(!is.na(truth)),
    #Bias       = bias,
    #MAPE       = mape,
    #MDAPE      = mdape,
    #Spearman_p = spearman_p,
    #CCC        = ccc_val,
    #CV_percent = cv_val
  )
}

flatten_dfs <- function(x, path = "") {
  out <- list()
  if (is.data.frame(x)) {
    out[[path]] <- x
  } else if (is.list(x)) {
    for (nm in names(x)) {
      new_path <- if (path == "") nm else paste(path, nm, sep = " > ")
      out <- c(out, flatten_dfs(x[[nm]], new_path))
    }
  }
  out
}

infer_loadtype <- function(truth_column, loadtypes) {
  if (length(loadtypes) == 1 || all(is.na(loadtypes))) {
    return(loadtypes[1])
  }
  lc_truth <- tolower(truth_column)
  chosen   <- NA_character_
  for (lt in loadtypes) {
    lt_l <- tolower(lt)
    if (lt_l == "flow cytometry" && grepl("fc|flow|facs", lc_truth)) {
      chosen <- lt; break
    }
    if (lt_l == "qpcr" && grepl("qpcr|copy number|ct", lc_truth)) {
      chosen <- lt; break
    }
    if (lt_l == "ddpcr" && grepl("ddpcr", lc_truth)) {
      chosen <- lt; break
    }
    if (lt_l == "mk_spike" && grepl("mk_spike", lc_truth)) {
      chosen <- lt; break
    }
    if (lt_l == "cfu" && grepl("cfu", lc_truth)) {
      chosen <- lt; break
    }
  }
  chosen
}

build_metrics_df <- function(
  predictions_and_demographics,
  exclude_truth_columns = NULL,  # e.g. c("fc_cells_per_ul_r1", "copies_ul", ...)
  exclude_studies       = NULL,  # e.g. c("Vandeputte2021", "Vandeputte2017")
  exclude_dataframes    = NULL,  # e.g. c("reprocessed > rdp16", "original > rdp16")
  exclude_loadtypes     = NULL,  # e.g. c("CFU", "MK_spike")
  verbose               = TRUE
) {
  all_metrics <- list()
  all_samples <- list()
  for (study_name in names(predictions_and_demographics)) {
    study_obj <- predictions_and_demographics[[study_name]]
    preds     <- study_obj$predictions
    loadtypes <- study_obj$studydemographics$loadtype %||% NA_character_
    dfs_list <- flatten_dfs(preds)
    if (length(dfs_list) == 0) next
    for (df_path in names(dfs_list)) {
      df <- dfs_list[[df_path]]
      req_cols <- c("sample ID", "log10_mlp", "profile", "training_set")
      if (!all(req_cols %in% names(df))) next
      truth_cols <- grep("^log10_truth_", names(df), value = TRUE)
      if (verbose) {
      message("-- ",paste(study_name), " -- has available truth cols before filtering: ", paste(truth_cols, collapse=", "))
      }
      if (length(truth_cols) == 0) next
      for (truth_col in truth_cols) {
        df_pair <- df %>%
          select(
            ID        = `sample ID`,
            Predicted = log10_mlp,
            Truth     = !!sym(truth_col),
            profile,
            training_set
          ) %>%
          filter(is.finite(Truth), is.finite(Predicted))
        if (nrow(df_pair) < 3 || var(df_pair$Truth) == 0) next
        df_pair <- df_pair %>%
          mutate(Source = vapply(training_set, get_label, FUN.VALUE = character(1))) %>%
          group_by(Source) %>%
          mutate(
            Pred_z  = mean_center(Predicted),
            Truth_z = mean_center(Truth)
          ) %>%
          ungroup()
        df_errors <- df_pair %>%
          mutate(
            residual      = Pred_z - Truth_z,
            study         = study_name,
            dataframe     = df_path,
            truth_column  = str_replace(truth_col, "^log10_truth_", ""),
            loadtype      = infer_loadtype(truth_col, loadtypes),
            predicted     = Predicted,
            truth         = Truth,
            Pred_z        = Pred_z,
            Truth_z       = Truth_z,
            Source        = Source,
            training_set  = training_set,
            profile       = profile
          ) %>%
          select(ID, study, dataframe, truth_column, Source, training_set, profile, loadtype,
                 predicted, truth, Pred_z, Truth_z, residual)
        all_samples[[length(all_samples) + 1]] <- df_errors
        srcs <- unique(df_pair$Source)
        if (length(srcs) == 1) {
          met <- compute_all_metrics(df_pair$Pred_z, df_pair$Truth_z)
          if (!is.null(met)) {
            all_metrics[[length(all_metrics) + 1]] <-
              met %>%
              mutate(
                study        = study_name,
                dataframe    = df_path,
                truth_column = str_replace(truth_col, "^log10_truth_", ""),
                Source       = srcs,
                loadtype     = infer_loadtype(truth_col, loadtypes)
              )
          }
        } else {
            by_src <- df_pair %>%
            group_by(Source) %>%
            summarize(
              compute_all_metrics(Pred_z, Truth_z),
              .groups = "drop"
            ) %>%
            mutate(
              study        = study_name,
              dataframe    = df_path,
              truth_column = str_replace(truth_col, "^log10_truth_", ""),
              loadtype     = infer_loadtype(truth_col, loadtypes)
            )
            all_metrics[[length(all_metrics) + 1]] <- by_src
        }
      }
    }
  }
  sample_df <- if (length(all_samples)) bind_rows(all_samples) else NULL
  study_df  <- if (length(all_metrics)) bind_rows(all_metrics)  else NULL
  if (!is.null(sample_df)) {
    sample_df <- sample_df %>%
      filter(
        !is.na(residual),
        if (!is.null(exclude_loadtypes))      !loadtype      %in% exclude_loadtypes      else TRUE,
        if (!is.null(exclude_truth_columns))  !truth_column  %in% exclude_truth_columns  else TRUE,
        if (!is.null(exclude_studies))        !study         %in% exclude_studies        else TRUE,
        if (!is.null(exclude_dataframes))     !dataframe     %in% exclude_dataframes     else TRUE
      )
  }
  if (!is.null(study_df)) {
    study_df <- study_df %>%
      filter(
        if (!is.null(exclude_loadtypes))      !loadtype      %in% exclude_loadtypes      else TRUE,
        if (!is.null(exclude_truth_columns))  !truth_column  %in% exclude_truth_columns  else TRUE,
        if (!is.null(exclude_studies))        !study         %in% exclude_studies        else TRUE,
        if (!is.null(exclude_dataframes))     !dataframe     %in% exclude_dataframes     else TRUE
      )
  }
  list(
    study_level  = study_df,
    sample_level = sample_df
  )
}

analyze_R2_by_author_presence <- function(
  predictions_and_demographics,
  seqtype,
  authors = c("Doris Vandeputte", "Jeroen Raes"),
  fill_by_source = FALSE,
  exclude_truth_columns = NULL,
  exclude_loadtypes = NULL,
  facetbyloadtype = FALSE
) {
  r2_list <- list()

  for (study_name in names(predictions_and_demographics)) {
    study_obj   <- predictions_and_demographics[[study_name]]
    dfs_list    <- flatten_dfs(study_obj$predictions)
    seqtype_all <- study_obj$studydemographics$sequencingtype %||% NA_character_
    loadtypes   <- study_obj$studydemographics$loadtype      %||% NA_character_

    for (df_path in names(dfs_list)) {
      df <- dfs_list[[df_path]]

      truth_cols <- grep("^log10_truth_", names(df), value = TRUE)
      if (!is.null(exclude_truth_columns)) {
        drop_names  <- paste0("log10_truth_", exclude_truth_columns)
        truth_cols  <- setdiff(truth_cols, drop_names)
      }
      if (!"log10_mlp" %in% names(df) || length(truth_cols) == 0) next

      for (truth_col in truth_cols) {
        df_pair <- df %>%
          select(
            Predicted = log10_mlp,
            Truth     = !!sym(truth_col),
            training_set
          ) %>%
          filter(is.finite(Predicted), is.finite(Truth))

        if (nrow(df_pair) < 3 || var(df_pair$Truth) == 0) next

        df_pair <- df_pair %>%
          mutate(
            Source = vapply(training_set, get_label, FUN.VALUE = character(1))
          ) %>%
          group_by(Source) %>%
          mutate(
            Pred_centered      = Predicted - mean(Predicted),
            Truth_centered     = Truth     - mean(Truth),
            One_minus_SSE_SST  = 1 - sum((Predicted - Truth)^2) / sum(Truth_centered^2)
          ) %>%
          ungroup()

        this_load <- infer_loadtype(truth_col, loadtypes)

        r2_df <- df_pair %>%
          group_by(Source) %>%
          summarize(
            sse = sum((Pred_centered - Truth_centered)^2),
            sst = sum((Truth_centered - mean(Truth_centered))^2),
            R2  = 1 - sse/sst,
            rho = cor(Pred_centered, Truth_centered, method="pearson"),
            .groups = "drop"
          ) %>%
          mutate(
            study        = study_name,
            sequencingtype = seqtype_all,
            loadtype     = this_load,
            dataframe    = df_path,
            truth_column = str_replace(truth_col, "^log10_truth_", "")
          )

        r2_list[[length(r2_list) + 1]] <- r2_df
      }
    }
  }

  tidyplot_theme <- function(base_size = 18) {
    theme_minimal(base_size = base_size) +
      theme(
        text              = element_text(family = "Helvetica", color = "black"),
        axis.text         = element_text(face = "bold", color = "black"),
        axis.line         = element_line(size = 0.75, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks        = element_line(size = 0.5, color = "black"),
        panel.grid        = element_blank(),
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        legend.title      = element_text(face = "bold"),
        legend.text       = element_text(),
        plot.subtitle     = element_text(hjust = 0.5)
      )
  }

  study_r2 <- bind_rows(r2_list)

  if (!is.null(exclude_truth_columns)) study_r2 <- study_r2 %>% filter(!truth_column %in% exclude_truth_columns)
  if (!is.null(exclude_loadtypes))    study_r2 <- study_r2 %>% filter(!loadtype     %in% exclude_loadtypes)

  studyinfo_df <- imap_dfr(
    predictions_and_demographics,
    ~{
      demo <- .x$studydemographics
      info <- .x$studyinfo
      tibble(
        study         = .y,
        sequencingtype = demo$sequencingtype %||% NA_character_,
        Authors       = unique(info$Authors) %>% paste(collapse = ";")
      )
    }
  ) %>%
    filter(!is.na(Authors), Authors != "") %>%
    mutate(authors_list = str_split(Authors, "\\s*;\\s*"))

  info   <- studyinfo_df %>% filter(sequencingtype == seqtype)
  r2_df  <- study_r2     %>% filter(study %in% info$study)

  pres_tbl <- info %>% select(study, authors_list)
  for (auth in authors) {
    pres_tbl <- pres_tbl %>% mutate(!!auth := map_lgl(authors_list, ~ auth %in% .x))
  }
  pres_tbl <- pres_tbl %>% select(-authors_list)

  df <- r2_df %>%
    filter(!is.na(R2)) %>%
    left_join(pres_tbl, by = "study") %>%
    pivot_longer(
      cols      = all_of(authors),
      names_to  = "Author",
      values_to = "Present"
    ) %>%
    filter(!is.na(Present))

  df <- df %>%
    filter(
      !(Source == "Vandeputte2021" & study == "Vandeputte2021"),
      !(Source == "Vandeputte2021" & study == "Vandeputte2017" &
        (dataframe == "original > rdp16" | dataframe == "reprocessed > rdp16")),
      !(Source == "Vandeputte2021" & study == "Vandeputte2017"),
      !(Source == "Vandeputte2021" & study == "Pereira2023" & dataframe == "reprocessed > rdp16"),
      !(Source == "Vandeputte2021" & study == "Krawczyk2022" & dataframe == "original > rdp16"),
      !(Source == "Vandeputte2021" & study == "Marotz2021" & dataframe == "reprocessed > rdp16" & truth_column == "all_flow_cellsperul_avg"),
      !(Source == "Vandeputte2021" & study == "Morton2019 16S" & dataframe == "reprocessed > amplicon > rdp16" & truth_column == "all_flow_cellsperul_avg")
    )

  status_tbl <- df %>%
    group_by(study, loadtype) %>%
    summarize(
      status   = if_else(any(Present), "Present", "Absent"),
      .groups  = "drop"
    )

  df <- df %>% left_join(status_tbl, by = c("study", "loadtype"))
  df2 = df %>%
    group_by(status, Source) %>%
    distinct(study, R2) %>%
    mutate(n = n()) %>%
    ungroup()

  valid_sources <- df2 %>%
  distinct(Source, status, n) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  filter(Absent > 3, Present > 3) %>%
  pull(Source)

  if (facetbyloadtype) {
    df2 = df %>%
    group_by(status, Source, loadtype) %>%
    distinct(study, R2) %>%
    mutate(n = n()) %>%
    ungroup()
    valid_sources <- df2 %>%
    distinct(Source, status, n, loadtype) %>%
    pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    filter(Absent > 3, Present > 3) %>%
    pull(Source)
  }
  comparisons_list <- list()
  for (src in valid_sources) {
    comparisons_list[[length(comparisons_list) + 1]] <- c("Absent", "Present")
  }
  highlight_sources <- c(
    "Vandeputte2017",
    "Vandeputte2017 Disease Cohort",
    "Vandeputte2017 Study Cohort",
    "GALAXY",
    "MetaCardis"
  )

  high_df <- df2 %>% 
    filter(study %in% highlight_sources)  %>%
    filter(
      !(study == "Vandeputte2017 Study Cohort"))

  df2 <- df2 %>%
    filter(
      !(study == "Vandeputte2017 Study Cohort"))
  y_limits <- c(0.6, 0.6)

  y_label <- switch(
      seqtype,
      "16S rRNA" = "<span style='color:#1B1B1B;'>Vandeputte2021</span>",
      "Shotgun Metagenomics" = paste0(
        "<span style='color:#0C58CA;'>GALAXY</span> / ",
        "<span style='color:#FF8C27;'>MetaCardis</span>"
      ),
      "R<sup>2</sup>"
    )

  g <- ggplot(df2, aes(
    x     = status,
    y     = R2,
    fill  = Source,
    color = Source
  )) +
    geom_violin(
      data     = filter(df2, n > 3),
      aes(fill = Source),
      alpha    = 0.5,
      width    = 0.8,
      size     = 0.8,
      position = position_dodge(width = 0.8),
      color    = "black",
      trim     = FALSE
    ) +
    geom_boxplot(
      data     = filter(df2, n > 3),
      if (seqtype == "Shotgun Metagenomics") { 
        aes(group = Source) 
      },
      fill         = "white",
      alpha        = 0.8,
      outlier.shape= NA,
      width        = 0.25,
      color        = "black", 
      position     = position_dodge(width = 0.8),
      size         = 0.8
    ) +
    geom_jitter(
      data     = filter(df2,
                    n <= 3#,
                    #!Source %in% highlight_sources
                    ),
      aes(color = Source),
      position = position_jitterdodge(dodge.width = 0.8),
      size     = 3,
      alpha    = 1,
      shape    = 21,
      stroke   = 1.5,
      fill     = "white"
    ) +
    # geom_jitter(
    #   data     = high_df,
    #   aes(x = status, y = R2, color = Source),
    #   position = position_jitterdodge(dodge.width = 0.8),
    #   size     = 3,
    #   alpha    = 1,
    #   shape    = 21,
    #   stroke   = 1.5,
    #   fill     = "white"
    # ) +
    # ggrepel::geom_text_repel(
    #   data       = high_df,
    #   aes(x = status, y = R2, fill = Source, label = study),
    #   position = position_jitterdodge(dodge.width = 0.8),
    #   size               = 3,
    #   max.overlaps       = Inf,
    #   min.segment.length = 0,
    #   box.padding        = 1,
    #   point.padding      = 0.5,
    #   segment.color      = "black",
    #   segment.angle      = 45,   
    #   force_pull         = 0.5,  
    #   direction          = "y",  
    #   seed               = 42    
    # ) +
    scale_fill_manual(
      name   = "Source",
      values = c(
        Vandeputte2021 = "#000000",
        GALAXY         = "#0C58CA",
        MetaCardis     = "#FF8C27"
      )
    ) +
    scale_color_manual(
      name   = "Source",
      values = c(
        Vandeputte2021 = "#000000",
        GALAXY         = "#0C58CA",
        MetaCardis     = "#FF8C27"
      )
    ) +
    scale_x_discrete(labels = c("Absent", "Present")) +
    scale_y_continuous(
      limits      = c(-1, 1),
      expand      = c(0, 0),
      breaks      = seq(-1, 1, by = 0.25),
      minor_breaks = seq(-1, 1, by = 0.125)
    ) +
    labs(x = NULL, y =' R<sup>2</sup>', title = y_label) +
    tidyplot_theme(base_size = 18) +
    theme(
      plot.title = element_markdown(face = "bold", hjust = 0.5),
      axis.title.x = element_text(face = "bold", margin = margin(t = 15),size = 18),
      axis.title.y = element_markdown(                    
                           face = "bold",
                           margin = margin(r = 15),
                           size = 18
                         ),
      legend.title = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.ticks = element_line(size = 0.5, color = "black"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      legend.position = "none",
      axis.text = element_text(size = 16, color = "black"),
      strip.text   = element_text(face = "bold", size = 20),
      strip.background = element_blank(),
      strip.placement = "outside"

  )  +
    stat_compare_means(
      data        = filter(df, Source %in% valid_sources),
      method      = "wilcox.test",
      comparisons = comparisons_list,  
      aes(group    = Source),          
      label       = "p.format",
      hide.ns     = FALSE,
      size        = 5,
      tip.length  = 0.02,
      vjust       = 0,
      label.y = 0.723,
      step.increase = 0.1,
    )

  if (facetbyloadtype) {
    g <- g + facet_wrap(~ loadtype)
  }

  return(list(plot = g, data = df))
}

plot_author_heatmap <- function(sim_df) {
  library(scales)
  library(ggtext)
  sim_clean <- sim_df %>%
    filter(
      !study1 %in% c("Vandeputte2017 Disease Cohort", "Vandeputte2017 Study Cohort", "Alessandri2024_vaginal"),
      !study2 %in% c("Vandeputte2017 Disease Cohort", "Vandeputte2017 Study Cohort", "Alessandri2024_vaginal")
    ) %>%
    mutate(
      study1 = ifelse(study1 == "Alessandri2024_fecal", "Alessandri2024", study1),
      study2 = ifelse(study2 == "Alessandri2024_fecal", "Alessandri2024", study2),
    )                                                                                                      
  mat_df <- sim_clean %>%
    pivot_wider(
      id_cols    = study1,
      names_from = study2,
      values_from = similarity
    )
  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$study1
  dist_mat <- as.dist(1 - mat)
  hc <- hclust(dist_mat, method = "average")
  ord <- hc$labels[hc$order]
  sim_ord <- sim_clean %>%
    mutate(
      study1 = factor(study1, levels = ord),
      study2 = factor(study2, levels = ord)
    )
  highlight <- c("Vandeputte2021", "GALAXY", "MetaCardis")
  lab_text <- sapply(ord, function(x) {
    if (x %in% highlight) paste0("<b>", x, "</b>")
    else x
  })
  y_i <- which(levels(sim_ord$study2) %in% highlight)
  p <- ggplot(sim_ord, aes(x = study1, y = study2, fill = similarity)) +
    geom_tile(color = "white") +
    annotate(
      "rect",
      xmin = 0.5, 
      xmax = length(levels(sim_ord$study1)) + 0.5,
      ymin = y_i   - 0.5,
      ymax = y_i   + 0.5,
      fill = NA,
      color = "gold",
      size  = 1,
      inherit.aes = FALSE
    ) +
    scale_fill_gradient(
      name   = "Author Similarity (%)",
      low    = "white",
      high   = "#0707a9",
      limits = c(0, 1),
      labels = percent_format(accuracy = 1)
    ) +
    scale_x_discrete(labels = lab_text) +
    scale_y_discrete(labels = lab_text) +
    coord_fixed() +
    theme_minimal(base_size = 18,  base_family = "Helvetica") +
    theme(
      axis.text.x = element_markdown(angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_markdown(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      axis.line   = element_line(size = 0.5, color = "black"),
      legend.title = element_text(size=16, margin=margin(b=15)),
      axis.ticks   = element_line(size = 0.5, color = "black")
    ) + geom_text(
    data = sim_ord %>% filter(study2 %in% highlight),
    aes(
      x     = study1,
      y     = study2,
      label = sprintf("%d", intersect),
      color = ifelse(similarity * 100 > 85, "white", "black")
    ),
    size       = 4,
    family     = "Helvetica",
    fontface   = "italic",
    show.legend = FALSE
  ) +
  scale_color_identity()

  return(p)
}

compute_author_similarity <- function(studyinfo_df, seqtype) {
  df <- studyinfo_df %>%
    filter(sequencingtype == seqtype)
  studs <- df$study
  lists <- df$authors_list
  mat <- tibble(expand.grid(i = seq_along(studs), j = seq_along(studs))) %>%
    mutate(
      study1 = studs[i],
      study2 = studs[j],
      intersect = map2_int(lists[i], lists[j], ~length(intersect(.x, .y))),
      union     = map2_int(lists[i], lists[j], ~length(union(.x, .y))),
      similarity = if_else(union > 0, intersect/union, 0)
    ) %>%
    select(study1, study2, similarity, intersect)
  
  mat
}

extract_studyinfo <- function(preds_and_demo) {
  imap_dfr(preds_and_demo, ~{
    info <- .x$studyinfo
    demo <- .x$studydemographics
    seqtype <- demo$sequencingtype %||% NA_character_
    authors <- info$Authors %>% unique() %>% paste(collapse = ";")
    tibble(
      study          = .y,
      sequencingtype = seqtype,
      Authors        = authors
    )
  }) %>%
  filter(!is.na(Authors), Authors != "") %>%
  mutate(
    authors_list = str_split(Authors, "\\s*;\\s*")
  )
}


run_aldex3_with_benchmarks <- function(countdata,conds,scaledata = NULL,scale_label = NULL,
                                       study_name = NULL,gammas = c(0, 0.5, 1, 5),
                                       checkpoint_dir = NULL,truthsamplesubset = NULL) {
  results <- list()
  if (!is.null(checkpoint_dir)) {
    if (is.null(study_name) || is.null(scale_label) || is.null(truthsamplesubset)) {
      stop("If you pass checkpoint_dir, you must supply 'study_name', 'scale_label', and 'truthsamplesubset'.")
    }
    if (!dir.exists(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE)
    }
    subset_label <- paste0("truthsub_", paste0(make.names(truthsamplesubset), collapse = "_"))
    ckpt_file <- file.path(
      checkpoint_dir,
      paste0(make.names(study_name), "_", subset_label, "_differentialabundance.rds")
    )
    if (file.exists(ckpt_file)) {
      message("► Loading existing checkpoint for ", study_name,
              " / ", scale_label, " / subset: ", paste(truthsamplesubset, collapse = ", "))
      results <- readRDS(ckpt_file)
    }
  }
  if ("SampleID" %in% colnames(conds)) {
    conds <- conds[, !(colnames(conds) %in% "SampleID"), drop = FALSE]
  }
  covariates <- as.data.frame(conds)
  drop_cols <- sapply(covariates, function(x) length(unique(x[!is.na(x)])) < 2)
  if (any(drop_cols)) {
    covariates <- covariates[, !drop_cols, drop = FALSE]
  }
  num_covariates <- ncol(covariates)
  if (num_covariates == 0) {
    message("No valid covariates available after dropping empty/constant columns.")
    return(results)
  }
  colnames(covariates) <- make.names(colnames(covariates), unique = TRUE)
  for (col in colnames(covariates)) {
    if (is.character(covariates[[col]]) || is.factor(covariates[[col]])) {
      covariates[[col]] <- factor(covariates[[col]])
    } else if (is.numeric(covariates[[col]])) {
      covariates[[col]] <- scale(covariates[[col]], center = TRUE, scale = TRUE)
    }
  }
  formula_str <- paste("~", paste(colnames(covariates), collapse = " + "))
  modelformula <- as.formula(formula_str)
  message("Design formula: ", formula_str)
  zero_count_rows <- rowSums(countdata) == 0
  if (sum(zero_count_rows) > 0) {
    message(sprintf("  • Filtering out %d features with zero counts", sum(zero_count_rows)))
    countdata <- countdata[!zero_count_rows, , drop = FALSE]
  } else {
    message("  • No zero-sum features found.")
  }
  if ("unclassified" %in% colnames(countdata)) {
    message("  • Filtering out unclassified features")
    countdata <- countdata[, !colnames(countdata) %in% "unclassified", drop = FALSE]
  }
  zero_sample_cols <- colSums(countdata) == 0
  if (sum(zero_sample_cols) > 0) {
    message(sprintf("  • Filtering out %d samples with zero counts", sum(zero_sample_cols)))
    keep_cols <- !zero_sample_cols
    countdata  <- countdata[, keep_cols, drop = FALSE]
    covariates <- covariates[keep_cols, , drop = FALSE]
    if (!is.null(scaledata)) {
      scaledata <- scaledata[keep_cols]
      names(scaledata) <- NULL
    }
  } else {
    message("  • No zero-sum samples found.")
  }
  if (!is.null(scaledata)) {
    scaledata <- as.numeric(scaledata)
  }
  if (nrow(countdata) == 0 || ncol(countdata) == 0) {
    message("No samples or features remain after filtering; skipping analysis for this subset.")
    return(results)
  }
  if (!is.null(scaledata)) {
    stopifnot(
      length(scaledata) == ncol(countdata),   
      length(scaledata) == nrow(covariates), 
      !any(is.na(scaledata))                  
    )
  }
  stopifnot(
    ncol(countdata) == nrow(covariates)
  )
  # ------------ ALDEx3 (TSS & CLR) ---------------------------------------------------------
  .append_results <- function(old_df, new_df) {
    if (is.null(old_df)) return(new_df)
    dplyr::bind_rows(old_df, new_df) |>
      dplyr::distinct(covariate, feature, gamma, scale_label, .keep_all = TRUE)
  }
  if (is.null(scaledata)) {
    message("Running ALDEx3 (no external scale)")
    existing_tss <- if ("ALDEx3_tss" %in% names(results)) results$ALDEx3_tss else NULL
    existing_clr <- if ("ALDEx3_clr" %in% names(results)) results$ALDEx3_clr else NULL
    for (g in gammas) {
      if (!is.null(existing_tss) && any(existing_tss$gamma == g)) {
        message(sprintf("  • Skipping ALDEx3[TSS], γ = %s (already in checkpoint)", g))
      } else {
        message(sprintf("  • ALDEx3[TSS], γ = %s", g))
        tss <- function(X, logWpara, gamma = 0.5) {
          P  <- nrow(X)
          ns <- dim(logWpara)[3]
          L  <- if (gamma == 0) matrix(0, P, ns) else
                matrix(rnorm(P * ns, 0, gamma), P, ns)
          t(X) %*% L
        }
        res_tss <- aldex(
          Y      = countdata,
          X      = modelformula,
          data   = covariates,
          stream = 3000,
          scale  = tss,
          gamma  = g
        )
        flat_tss <- flatten_aldex_results(
          setNames(list(res_tss), g),
          countdata,
          study_name,
          scale_label = "ALDEx3_tss",
          use         = "mean"
        )
        existing_tss       <- .append_results(existing_tss, flat_tss)
        results$ALDEx3_tss <- existing_tss
        if (!is.null(checkpoint_dir)) {
          saveRDS(results, ckpt_file)
          message("    ⇢ Saved checkpoint (ALDEx3_tss up to γ=", g, ") to: ", ckpt_file)
        }
      }
      if (!is.null(existing_clr) && any(existing_clr$gamma == g)) {
        message(sprintf("  • Skipping ALDEx3[CLR], γ = %s (already in checkpoint)", g))
      } else {
        message(sprintf("  • ALDEx3[CLR], γ = %s", g))
        clr <- function(X, logWpara, gamma = 0.5) {
          P  <- nrow(X)
          ns <- dim(logWpara)[3]
          mu <- -colMeans(logWpara, dims = 1)
          L  <- if (gamma == 0) matrix(0, P, ns) else
                matrix(rnorm(P * ns, 0, gamma), P, ns)
          mu + t(X) %*% L
        }
        res_clr <- aldex(
          Y      = countdata,
          X      = modelformula,
          data   = covariates,
          stream = 1000,
          scale  = clr,
          gamma  = g
        )
        flat_clr <- flatten_aldex_results(
          setNames(list(res_clr), g),
          countdata,
          study_name,
          scale_label = "ALDEx3_clr",
          use         = "mean"
        )
        existing_clr       <- .append_results(existing_clr, flat_clr)
        results$ALDEx3_clr <- existing_clr
        if (!is.null(checkpoint_dir)) {
          saveRDS(results, ckpt_file)
          message("    ⇢ Saved checkpoint (ALDEx3_clr up to γ=", g, ") to: ", ckpt_file)
        }
      }
    }
  } else {
    message("Running ALDEx3 (external scale)")
    result_name  <- paste0("ALDEx3_", scale_label)
    existing_ext <- if (result_name %in% names(results)) results[[result_name]] else NULL
    for (g in gammas) {
      if (!is.null(existing_ext) && any(existing_ext$gamma == g)) {
        message(sprintf("  • Skipping ALDEx3[%s], γ = %s (already in checkpoint)", scale_label, g))
      } else {
        message(sprintf("  • ALDEx3[%s], γ = %s", scale_label, g))
        externalscale <- function(X, logWpara, gamma = 0.5, scaledata) {
          N         <- length(scaledata)
          nsample   <- dim(logWpara)[3]
          gamma_vec <- if (length(gamma) == 1) rep(gamma, N) else gamma
          t(sapply(seq_len(N), function(i) {
            if (gamma == 0) {
              rep(scaledata[i], nsample)
            } else {
              rnorm(nsample, mean = scaledata[i], sd = gamma_vec[i])
            }
          }))
        }
        res_ext <- aldex(
          Y         = countdata,
          X         = modelformula,
          data      = covariates,
          stream    = 3000,
          scale     = externalscale,
          gamma     = g,
          scaledata = scaledata
        )
        flat_ext            <- flatten_aldex_results(
          setNames(list(res_ext), g),
          countdata,
          study_name,
          scale_label = scale_label,
          use         = "mean"
        )
        existing_ext        <- .append_results(existing_ext, flat_ext)
        results[[result_name]] <- existing_ext
        if (!is.null(checkpoint_dir)) {
          saveRDS(results, ckpt_file)
          message("    ⇢ Saved checkpoint (", result_name, " up to γ=", g, ") to: ", ckpt_file)
        }
      }
    }
  }
  # ------------ DESeq2 ---------------------------------------------------------
  if (is.null(scaledata)) {
    if (!("DESeq2" %in% names(results)) ||
        (is.data.frame(results$DESeq2) && nrow(results$DESeq2) == 0)) {
      message("Running DESeq2")
      for (col in colnames(covariates)) {
        if (is.ordered(covariates[[col]])) {
          covariates[[col]] <- factor(covariates[[col]], ordered = FALSE)
        }
      }
      design <- model.matrix(modelformula, data = covariates)
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countdata,
        colData   = covariates,
        design    = design
      )
      dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
      dds <- DESeq2::DESeq(dds)
      coef_names   <- DESeq2::resultsNames(dds)
      design_names <- colnames(design)[-1]
      if (length(design_names) != length(coef_names) - 1) {
        stop("Mismatch between design matrix columns and DESeq2 coefficients.")
      }
      coef_map <- c("Intercept", design_names)
      names(coef_map) <- coef_names
      deseq_results <- lapply(coef_names[-1], function(cov) {
        res <- DESeq2::results(dds, name = cov)
        res_df <- as.data.frame(res) %>%
          tibble::rownames_to_column(var = "feature") %>%
          dplyr::select(feature, log2FoldChange, pvalue, padj) %>%
          dplyr::rename(estimate = log2FoldChange,
                        p.val    = pvalue,
                        p.val.adj= padj) %>%
          dplyr::mutate(
            covariate   = unname(as.character(coef_map[cov])),
            study_name  = study_name,
            scale_label = "DESeq2",
            gamma       = "0"
          )
        res_df
      })
      results$DESeq2 <- dplyr::bind_rows(deseq_results)
      if (!is.null(checkpoint_dir)) {
        saveRDS(results, ckpt_file)
        message("    ⇢ Saved checkpoint (DESeq2) to: ", ckpt_file)
      }
    } else {
      message("Skipping DESeq2 (already in checkpoint).")
    }
  }
  # ------------ limma-voom ---------------------------------------------------------
  if (is.null(scaledata)) {
    if (!("limma" %in% names(results)) ||
        (is.data.frame(results$limma) && nrow(results$limma) == 0)) {
      message("Running limma-voom")
      design <- model.matrix(modelformula, data = covariates)
      v <- limma::voom(countdata, design, plot = FALSE)
      fit <- limma::lmFit(v, design)
      fit <- limma::eBayes(fit)
      limma_results <- lapply(colnames(design)[-1], function(cov) {
        res <- limma::topTable(fit, coef = cov, number = Inf, adjust.method = "BH")
        res_df <- as.data.frame(res) %>%
          dplyr::mutate(feature = paste0(rownames(res))) %>%
          dplyr::select(feature, logFC, P.Value, adj.P.Val) %>%
          dplyr::rename(estimate = logFC,
                        p.val    = P.Value,
                        p.val.adj= adj.P.Val) %>%
          dplyr::mutate(
            covariate   = cov,
            study_name  = study_name,
            scale_label = "limma voom",
            gamma       = "0"
          )
        res_df
      })
      limma_results <- dplyr::bind_rows(limma_results)
      rownames(limma_results) <- NULL
      results$limma <- limma_results
      if (!is.null(checkpoint_dir)) {
        saveRDS(results, ckpt_file)
        message("    ⇢ Saved checkpoint (limma-voom) to: ", ckpt_file)
      }
    } else {
      message("Skipping limma-voom (already in checkpoint).")
    }
  }

    
  # ------------ edgeR ---------------------------------------------------------
  if (is.null(scaledata)) {
    ## ---------------- edgeR (GLM) ----------------
    if (!("edgeR" %in% names(results)) ||
        (is.data.frame(results$edgeR) && nrow(results$edgeR) == 0)) {
      message("Running edgeR (GLM)")
      design <- model.matrix(modelformula, data = covariates)
      y <- edgeR::DGEList(counts = countdata)
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, design = design)
      fit <- edgeR::glmFit(y, design)
      coef_names <- colnames(design)[-1] 
      edger_results <- lapply(seq_along(coef_names), function(i) {
        lrt <- edgeR::glmLRT(fit, coef = i + 1) 
        tt  <- edgeR::topTags(lrt, n = Inf)$table
        df  <- tt |>
          tibble::rownames_to_column(var = "feature") |>
          dplyr::transmute(
            feature,
            estimate   = logFC,
            p.val      = PValue,
            p.val.adj  = FDR,
            covariate  = coef_names[i],
            study_name = study_name,
            scale_label= "edgeR",
            gamma      = "0"
          )
        df
      })
      results$edgeR <- dplyr::bind_rows(edger_results)
      if (!is.null(checkpoint_dir)) {
        saveRDS(results, ckpt_file)
        message("    ⇢ Saved checkpoint (edgeR) to: ", ckpt_file)
      }
    } else {
      message("Skipping edgeR (already in checkpoint).")
    }
  }

  # ------------ metagenomeSeq (CSS + fitZig) -----------------------------------
  if (is.null(scaledata)) {
    if (!("metagenomeSeq" %in% names(results)) ||
        (is.data.frame(results$metagenomeSeq) && nrow(results$metagenomeSeq) == 0)) {

      if (!requireNamespace("metagenomeSeq", quietly = TRUE) ||
          !requireNamespace("Biobase", quietly = TRUE)) {
        message("Skipping metagenomeSeq: package(s) not installed.")
        results$metagenomeSeq <- data.frame()

      } else {
        message("Running metagenomeSeq (CSS + fitZig)")

        ## --- align: rownames(covariates) must match colnames(countdata)
        if (is.null(rownames(covariates))) rownames(covariates) <- colnames(countdata)
        if (!identical(rownames(covariates), colnames(countdata))) {
          common <- intersect(rownames(covariates), colnames(countdata))
          if (!length(common)) stop("No common sample IDs between covariates and countdata.")
          covariates <- covariates[common, , drop = FALSE]
          countdata  <- countdata[,  common, drop = FALSE]
        }
        stopifnot(identical(rownames(covariates), colnames(countdata)))

        ## --- parameters for pruning
        min_nonzero <- 3L
        min_prop    <- 0.10

        ## --- initial sample prune: ≤1 non-zero feature
        nz_per_sample <- colSums(countdata > 0, na.rm = TRUE)
        drop_samp <- names(nz_per_sample[nz_per_sample <= 1L])
        if (length(drop_samp)) {
          message("  • Dropping ", length(drop_samp),
                  " sample(s) with ≤1 non-zero feature: ",
                  paste(drop_samp, collapse = ", "))
          keep <- setdiff(colnames(countdata), drop_samp)
          countdata  <- countdata[, keep, drop = FALSE]
          covariates <- covariates[keep, , drop = FALSE]
        }

        ## --- iterative prune: features then samples until stable
        sweeps  <- 0L
        changed <- TRUE
        while (changed && sweeps < 6L) {
          sweeps  <- sweeps + 1L
          changed <- FALSE

          # feature filter: ≥ min_nonzero AND present in ≥ ceil(min_prop * n_samples)
          nz_feat  <- rowSums(countdata > 0, na.rm = TRUE)
          thresh   <- ceiling(min_prop * ncol(countdata))
          keep_feat <- (nz_feat >= min_nonzero) & (nz_feat >= thresh)
          if (any(!keep_feat)) {
            message("  • Filtering out ", sum(!keep_feat),
                    " features failing nonzero thresholds (≥", min_nonzero,
                    " and ≥", thresh, " samples).")
            countdata <- countdata[keep_feat, , drop = FALSE]
            changed <- TRUE
          }

          # sample filter again (after feature drops)
          nz_samp   <- colSums(countdata > 0, na.rm = TRUE)
          drop_samp <- names(nz_samp[nz_samp <= 1L])
          if (length(drop_samp)) {
            message("  • Dropping ", length(drop_samp),
                    " sample(s) now ≤1 non-zero after feature filter: ",
                    paste(drop_samp, collapse = ", "))
            keep <- setdiff(colnames(countdata), drop_samp)
            countdata  <- countdata[, keep, drop = FALSE]
            covariates <- covariates[keep, , drop = FALSE]
            changed <- TRUE
          }
        }

        ## --- final guard
        message("  • Final after filtering: ", nrow(countdata), " features × ",
                ncol(countdata), " samples.")
        if (nrow(countdata) < 2L || ncol(countdata) < 3L) {
          message("  • Not enough features/samples left; storing empty frame.")
          results$metagenomeSeq <- data.frame()
          if (!is.null(checkpoint_dir)) { saveRDS(results, ckpt_file) }
        } else {
          ## --- build design (ensure full rank)
          design <- model.matrix(modelformula, data = covariates)
          qrX <- qr(design)
          if (qrX$rank < ncol(design)) {
            keep_idx <- qrX$pivot[seq_len(qrX$rank)]
            dropped  <- setdiff(seq_len(ncol(design)), keep_idx)
            message("  • Design not full rank; dropping columns: ",
                    paste(colnames(design)[dropped], collapse = ", "))
            design <- design[, keep_idx, drop = FALSE]
          }
          if (ncol(design) <= 1L) {
            message("  • Design has no non-intercept coefficients; storing empty frame.")
            results$metagenomeSeq <- data.frame()
            if (!is.null(checkpoint_dir)) { saveRDS(results, ckpt_file) }
          } else {
            ## --- try CSS; if it fails, fall back to TMM+limma
            css_ok <- TRUE
            pheno  <- Biobase::AnnotatedDataFrame(as.data.frame(covariates))
            obj    <- metagenomeSeq::newMRexperiment(counts = as.matrix(countdata), phenoData = pheno)

            p_css <- tryCatch(
              metagenomeSeq::cumNormStatFast(obj),
              error = function(e) { message("  • cumNormStatFast error: ", e$message); css_ok <<- FALSE; NA_real_ }
            )
            if (css_ok) {
              obj <- tryCatch(
                metagenomeSeq::cumNorm(obj, p = p_css),
                error = function(e) { message("  • cumNorm error: ", e$message); css_ok <<- FALSE; obj }
              )
            }

            if (!css_ok) {
              message("  • Falling back to edgeR TMM + limma (CSS unavailable).")
              y <- edgeR::DGEList(counts = as.matrix(countdata))
              y <- edgeR::calcNormFactors(y, method = "TMM")
              v <- limma::voom(y, design, plot = FALSE)
              fit <- limma::lmFit(v, design)
              fit <- limma::eBayes(fit)

              coef_names <- colnames(design)[-1]
              lim_res <- lapply(seq_along(coef_names), function(i) {
                tt <- limma::topTable(fit, coef = i + 1, number = Inf, adjust.method = "BH")
                data.frame(
                  feature     = rownames(tt),
                  estimate    = tt$logFC,
                  p.val       = tt$P.Value,
                  p.val.adj   = tt$adj.P.Val,
                  covariate   = coef_names[i],
                  study_name  = study_name,
                  scale_label = "metagenomeSeq (TMM+limma fallback)",
                  gamma       = "0",
                  stringsAsFactors = FALSE
                )
              })
              results$metagenomeSeq <- dplyr::bind_rows(lim_res)

            } else {
              ## --- CSS succeeded; try fitZig, then eBayes/topTable
              zig <- tryCatch(
                metagenomeSeq::fitZig(obj = obj, mod = design),
                error = function(e) { message("  • fitZig error: ", conditionMessage(e)); NULL }
              )
              if (is.null(zig)) {
                message("  • Falling back to CSS-normalized limma (fitZig error).")
                css <- metagenomeSeq::MRcounts(obj, norm = TRUE)
                v   <- limma::voom(css, design, plot = FALSE)
                fit <- limma::lmFit(v, design)
                fit <- limma::eBayes(fit)

                coef_names <- colnames(design)[-1]
                lim_res <- lapply(seq_along(coef_names), function(i) {
                  tt <- limma::topTable(fit, coef = i + 1, number = Inf, adjust.method = "BH")
                  data.frame(
                    feature     = rownames(tt),
                    estimate    = tt$logFC,
                    p.val       = tt$P.Value,
                    p.val.adj   = tt$adj.P.Val,
                    covariate   = coef_names[i],
                    study_name  = study_name,
                    scale_label = "metagenomeSeq (CSS+limma fallback)",
                    gamma       = "0",
                    stringsAsFactors = FALSE
                  )
                })
                results$metagenomeSeq <- dplyr::bind_rows(lim_res)

              } else {
                zigFit <- NULL
                if (methods::is(zig, "fitZigResults")) {
                  if ("eb" %in% slotNames(zig) && methods::is(zig@eb, "MArrayLM")) zigFit <- zig@eb
                  else if ("fit" %in% slotNames(zig)) zigFit <- zig@fit
                } else if (is.list(zig)) {
                  if (!is.null(zig$eb) && inherits(zig$eb, "MArrayLM")) zigFit <- zig$eb
                  else if (!is.null(zig$fit)) zigFit <- zig$fit
                }
                if (is.null(zigFit)) stop("Unexpected object returned by fitZig: ", paste(class(zig), collapse = ", "))
                if (is.null(zigFit$s2.post)) zigFit <- limma::eBayes(zigFit)

                coef_names <- colnames(design)[-1]
                mgs_results <- lapply(seq_along(coef_names), function(i) {
                  tt <- limma::topTable(zigFit, coef = i + 1, number = Inf, adjust.method = "BH")
                  data.frame(
                    feature     = rownames(tt),
                    estimate    = tt$logFC,
                    p.val       = tt$P.Value,
                    p.val.adj   = tt$adj.P.Val,
                    covariate   = coef_names[i],
                    study_name  = study_name,
                    scale_label = "metagenomeSeq (CSS+fitZig)",
                    gamma       = "0",
                    stringsAsFactors = FALSE
                  )
                })
                results$metagenomeSeq <- dplyr::bind_rows(mgs_results)
              }
            }

            if (!is.null(checkpoint_dir)) {
              saveRDS(results, ckpt_file)
              message("    ⇢ Saved checkpoint (metagenomeSeq) to: ", ckpt_file)
            }
          }
        }
      }
    } else {
      message("Skipping metagenomeSeq (already in checkpoint).")
    }
  }

  # ------------ LinDA (MicrobiomeStat) -------------------------------------------------
  if (is.null(scaledata)) {
    if (!("LinDA" %in% names(results)) ||
        (is.data.frame(results$LinDA) && nrow(results$LinDA) == 0)) {

      if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
        message("Skipping LinDA: MicrobiomeStat package not installed.")
        results$LinDA <- data.frame()

      } else {
        message("Running LinDA via MicrobiomeStat::linda (tidy output; auto-handle variants)")

        # Ensure features = rows, samples = cols
        feature_table <- as.matrix(countdata)
        meta_data     <- as.data.frame(covariates)

        # Make sure samples align by name (rownames(meta) must equal colnames(counts))
        if (is.null(rownames(meta_data)) || !identical(rownames(meta_data), colnames(feature_table))) {
          if (!is.null(colnames(feature_table))) {
            rownames(meta_data) <- colnames(feature_table)
          } else {
            colnames(feature_table) <- rownames(meta_data)
          }
        }

        # MicrobiomeStat::linda expects a character formula (e.g., "~ Group + Age")
        form_str <- deparse(modelformula)

        run_linda_once <- function(is_winsor = TRUE, outlier_pct = 0.03) {
          tryCatch(
            MicrobiomeStat::linda(
              feature.dat      = feature_table,
              meta.dat         = meta_data,
              formula          = form_str,
              feature.dat.type = "count",
              prev.filter      = 0,          # don't drop by prevalence here
              is.winsor        = is_winsor,  # winsorization toggle
              outlier.pct      = outlier_pct,
              p.adj.method     = "BH",
              alpha            = 0.05,
              verbose          = TRUE
            ),
            error = function(e) {
              message("  • LinDA failed: ", conditionMessage(e))
              NULL
            }
          )
        }

        parse_linda_output <- function(linda_out) {
          if (is.null(linda_out)) return(data.frame())

          norm_one <- function(df, cov_name) {
            df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
            feat <- rownames(df)
            if (is.null(feat) || anyNA(feat) || any(feat == "")) {
              feat <- if ("taxon" %in% names(df)) df$taxon else
                      if ("feature" %in% names(df)) df$feature else
                      as.character(seq_len(nrow(df)))
            }
            need <- c("log2FoldChange","pvalue","padj")
            if (!all(need %in% names(df))) return(data.frame())
            data.frame(
              feature     = as.character(feat),
              estimate    = as.numeric(df[["log2FoldChange"]]),
              p.val       = as.numeric(df[["pvalue"]]),
              p.val.adj   = as.numeric(df[["padj"]]),
              covariate   = as.character(cov_name),
              study_name  = study_name,
              scale_label = "LinDA",
              gamma       = "0",
              stringsAsFactors = FALSE
            )
          }

          if (!is.null(linda_out$output)) {
            if (is.list(linda_out$output)) {
              parts <- lapply(names(linda_out$output),
                              function(v) norm_one(linda_out$output[[v]], v))
              out <- dplyr::bind_rows(parts)
              return(out)
            } else if (is.data.frame(linda_out$output)) {
              cov_name <- if (!is.null(linda_out$variables) && length(linda_out$variables) >= 1)
                            as.character(linda_out$variables[[1]]) else "(covariate)"
              return(norm_one(linda_out$output, cov_name))
            }
          }
          if (!is.null(linda_out$res) && is.data.frame(linda_out$res)) {
            cov_name <- if (!is.null(linda_out$variables) && length(linda_out$variables) >= 1)
                          as.character(linda_out$variables[[1]]) else "(covariate)"
            return(norm_one(linda_out$res, cov_name))
          }
          data.frame()
        }

        # First attempt (winsorization on; no lib.size arg!)
        linda_out  <- run_linda_once(is_winsor = TRUE, outlier_pct = 0.03)
        linda_res  <- parse_linda_output(linda_out)

        # If empty, try once more with winsorization OFF (some datasets behave better)
        if (!is.data.frame(linda_res) || nrow(linda_res) == 0) {
          message("  • LinDA produced empty results; retrying with winsorization OFF...")
          linda_out2 <- run_linda_once(is_winsor = FALSE, outlier_pct = 0.03)
          linda_res  <- parse_linda_output(linda_out2)
        }

        if (!is.data.frame(linda_res)) linda_res <- data.frame()
        results$LinDA <- linda_res

        if (!nrow(results$LinDA)) {
          message("  • LinDA returned no significant/parsable results; storing empty frame.")
        }

        if (!is.null(checkpoint_dir)) {
          saveRDS(results, ckpt_file)
          message("    ⇢ Saved checkpoint (LinDA) to: ", ckpt_file)
        }
      }
    } else {
      message("Skipping LinDA (already in checkpoint).")
    }
  }

  # ------------ ANCOM-BC2 ---------------------------------------------------------
  if (is.null(scaledata)) {
    if (!("ANCOM_BC" %in% names(results)) ||
        (is.data.frame(results$ANCOM_BC) && nrow(results$ANCOM_BC) == 0)) {
      message("Running ANCOM-BC2")
      feature_table <- as.matrix(countdata)
      meta_data <- covariates
      if (!("sample_id" %in% colnames(meta_data))) {
        meta_data <- meta_data %>%
          tibble::rownames_to_column(var = "sample_id")
      }
      if (!all(meta_data$sample_id == colnames(feature_table))) {
        meta_data$sample_id <- colnames(feature_table)
      }
      rownames(meta_data) <- meta_data$sample_id
      formula_str2 <- sub("^\\s*~\\s*", "", formula_str)
      group_var <- colnames(covariates)[1]
      ancombc2_out <- ANCOMBC::ancombc2(
        data         = feature_table,
        meta_data    = meta_data,
        fix_formula  = formula_str2,
        group        = group_var,
        p_adj_method = "BH",
        prv_cut      = 0,
        lib_cut      = 0,
        struc_zero   = FALSE,
        neg_lb       = FALSE,
        alpha        = 0.05,
        verbose      = TRUE
      )
      res_df      <- ancombc2_out$res
      taxon_names <- res_df$taxon
      lfc_cols    <- grep("^lfc_", colnames(res_df), value = TRUE)
      lfc_cols    <- lfc_cols[!grepl("\\(Intercept\\)$", lfc_cols)]
      ancombc2_results <- lapply(lfc_cols, function(lfc_col) {
        covariate_name <- sub("^lfc_", "", lfc_col)
        data.frame(
          feature       = taxon_names,
          estimate      = res_df[[lfc_col]],
          p.val         = res_df[[paste0("p_", covariate_name)]],
          p.val.adj     = res_df[[paste0("q_", covariate_name)]],
          covariate     = covariate_name,
          study_name    = study_name,
          scale_label   = "ANCOM-BC2",
          gamma         = "0",
          stringsAsFactors = FALSE
        )
      })
      results$ANCOM_BC <- dplyr::bind_rows(ancombc2_results)
      if (!is.null(checkpoint_dir)) {
        saveRDS(results, ckpt_file)
        message("    ⇢ Saved checkpoint (ANCOM-BC2) to: ", ckpt_file)
      }
    } else {
      message("Skipping ANCOM-BC2 (already in checkpoint).")
    }
  }
  if (!is.null(checkpoint_dir)) {
    saveRDS(results, ckpt_file)
    message("✔︎ Finished all models. Final checkpoint written to: ", ckpt_file)
  }
  return(results)
}

process_study_with_benchmarks <- function(study_name, study_data, checkpoint_dir = NULL) {
  message("Differential Abundance Analysis for: ", study_name)

  demo   <- study_data$studydemographics
  id_col <- demo$ID
  counts_df <- study_data$counts
  df_meta   <- study_data$metadata
  df_pred   <- study_data$predictions

  if (is.null(counts_df)) {
    message("No counts matrix for ", study_name, ". Skipping.")
    return(NULL)
  }
  if (!(id_col %in% colnames(df_meta))) {
    stop("ID column '", id_col, "' not found in metadata for ", study_name, ".")
  }

  sample_id_col <- grep("^sample[._ ]?ID$", colnames(df_pred), ignore.case = TRUE, value = TRUE)
  if (length(sample_id_col) != 1) {
    stop("`predictions` must contain a column named `sample ID` (or `sample.ID`) for ", study_name, ".")
  }

  # 1) Get sample IDs for counts (rows are samples here)
  if (!is.null(rownames(counts_df)) && !any(rownames(counts_df) == "")) {
    sample_ids_counts <- rownames(counts_df)
  } else if (nrow(counts_df) == nrow(df_meta)) {
    sample_ids_counts <- as.character(df_meta[[id_col]])
    rownames(counts_df) <- sample_ids_counts
  } else {
    stop("Cannot infer sample IDs for counts in ", study_name, ".")
  }

  # 2) Intersect sample IDs across data sources
  meta_ids   <- as.character(df_meta[[id_col]])
  pred_ids   <- as.character(df_pred[[sample_id_col]])
  common_ids <- Reduce(intersect, list(sample_ids_counts, meta_ids, pred_ids))
  if (!length(common_ids)) {
    stop("No common sample IDs among counts, metadata, and predictions for ", study_name, ".")
  }

  # 3) Subset to common IDs (keep metadata order flexible for now)
  counts_df <- counts_df[common_ids, , drop = FALSE]
  df_meta   <- df_meta[match(common_ids, df_meta[[id_col]]), , drop = FALSE]
  df_pred   <- df_pred[df_pred[[sample_id_col]] %in% common_ids, , drop = FALSE]

  # 4) Covariates and complete cases
  covariate_cols <- demo$covariates
  if (is.null(covariate_cols) || !length(covariate_cols)) {
    message("Skipping study: ", study_name, " – no covariates defined.")
    return(NULL)
  }
  conds_full <- df_meta[, covariate_cols, drop = FALSE]
  rownames(conds_full) <- df_meta[[id_col]]  # <<< IMPORTANT: rows carry sample IDs

  keep_idx <- complete.cases(conds_full)
  if (!any(keep_idx)) {
    message("Skipping study: ", study_name, " – zero samples have complete covariate data.")
    return(NULL)
  }

  # Keep only complete rows in ALL tables (preserve rownames = sample IDs)
  counts_df  <- counts_df[keep_idx, , drop = FALSE]
  df_meta    <- df_meta[keep_idx, , drop = FALSE]
  df_pred    <- df_pred[df_pred[[sample_id_col]] %in% df_meta[[id_col]], , drop = FALSE]
  conds_full <- conds_full[keep_idx, , drop = FALSE]

  # 5) Drop rare factor levels and keep sample IDs as rownames
  for (v in covariate_cols) {
    if (!is.factor(conds_full[[v]]) && !is.character(conds_full[[v]])) next
    freq <- table(conds_full[[v]])
    rare <- names(freq[freq < 2])
    if (length(rare)) {
      drop_idx <- conds_full[[v]] %in% rare
      message("Dropping ", sum(drop_idx), " sample(s) from rare level(s) in ", v, ": ", paste(rare, collapse = ", "))
      keep2        <- !drop_idx
      conds_full   <- conds_full[ keep2, , drop = FALSE]
      counts_df    <- counts_df[ keep2, , drop = FALSE]
      df_meta      <- df_meta[   keep2, , drop = FALSE]
      # df_pred filtered later by ID map
    }
  }

  # 6) Build count matrix with samples as COLUMNS and align strictly to conds_full row order
  counts_mat <- t(as.matrix(counts_df))
  # Ensure colnames(counts_mat) are sample IDs
  colnames(counts_mat) <- rownames(counts_df)

  # <<< ALIGNMENT CORE: use covariate rownames to index count columns >>>
  # Reorder counts columns to exactly match the row order of conds_full
  sample_order <- rownames(conds_full)
  counts_mat   <- counts_mat[, sample_order, drop = FALSE]

  storage.mode(counts_mat) <- "integer"

  # 7) Prepare outputs
  final_results_list <- list()
  truth_cols <- grep("^log2_truth_", colnames(df_pred), value = TRUE)
  message("DEBUG: truth_cols = ", if (length(truth_cols)) paste(truth_cols, collapse = ", ") else "<none>")

  # For convenient lookups later
  # Named vectors keyed by sample ID
  pred_id_vec <- as.character(df_pred[[sample_id_col]])
  names(pred_id_vec) <- pred_id_vec

  if (!length(truth_cols)) {
    # -------- No truth: just run with all samples, already aligned by rownames(conds_full)
    n_samples <- nrow(conds_full)
    message(" >> no_scale on all ", n_samples, " samples <<")

    no_scale_res <- run_aldex3_with_benchmarks(
      countdata        = counts_mat,                  # columns = rownames(conds_full)
      conds            = conds_full,                  # rows = sample IDs
      scaledata        = NULL,
      scale_label      = "no_scale",
      study_name       = study_name,
      checkpoint_dir   = checkpoint_dir,
      truthsamplesubset = "no_truth_subset"
    )
    if (length(no_scale_res)) {
      no_scale_res <- purrr::imap(no_scale_res, \(df, .) {
        if (is.data.frame(df) && nrow(df))
          dplyr::mutate(df, n_samples = n_samples)
      })
    }
    final_results_list[["no_truth_subset_no_scale"]] <-
      no_scale_res[names(no_scale_res) %in% c("ALDEx3_tss", "ALDEx3_clr", "DESeq2", "limma", "ANCOM_BC", "LinDA","edgeR","metagenomeSeq")]

  } else {
    # -------- With truths: build named vectors by sample ID, then index by covariate rows
    for (truth_col in truth_cols) {

      # Map truth to sample IDs (named vector), then pull in covariate row order
      truth_map <- setNames(df_pred[[truth_col]], pred_id_vec)
      truth_all <- truth_map[rownames(conds_full)]   # <<< rows index columns here

      non_na_idx <- which(!is.na(truth_all))
      n_samples  <- length(non_na_idx)
      if (!n_samples) {
        message("  Skipping ", truth_col, " (0 samples).")
        next
      }

      message(" >> truth = ", truth_col, " | ", n_samples, " samples <<")

      sample_ids_sub <- rownames(conds_full)[non_na_idx]
      counts_sub     <- counts_mat[, sample_ids_sub, drop = FALSE]  # columns = these sample IDs
      conds_sub      <- conds_full[sample_ids_sub, , drop = FALSE]  # rows   = these sample IDs
      truth_sub      <- truth_all[non_na_idx]

      # Sanity check dimensions
      if (length(colnames(counts_sub)) != nrow(conds_sub)) {
        message("Dimensions mismatch for ", study_name, " (", truth_col, ").")
        next
      }

      # No-scale for this subset
      no_scale_res <- run_aldex3_with_benchmarks(
        countdata        = counts_sub,
        conds            = conds_sub,
        scaledata        = NULL,
        scale_label      = "no_scale",
        study_name       = study_name,
        checkpoint_dir   = checkpoint_dir,
        truthsamplesubset = truth_col
      )
      if (length(no_scale_res)) {
        no_scale_res <- purrr::imap(no_scale_res, \(df, .) {
          if (is.data.frame(df) && nrow(df))
            dplyr::mutate(df, n_samples = n_samples)
        })
      }
      final_results_list[[truth_col]][["no_scale"]] <-
        no_scale_res[names(no_scale_res) %in% c("ALDEx3_tss", "ALDEx3_clr", "DESeq2", "limma", "ANCOM_BC", "LinDA","edgeR","metagenomeSeq")]

      # Measured scale = truth_sub (same sample ordering)
      if (length(truth_sub) != length(colnames(counts_sub)) ||
          length(truth_sub) != nrow(conds_sub)) {
        message("Truth have the wrong number of samples/mismatched dimensions for ",
                study_name, " (", truth_col, ").")
        next
      }

      meas_res <- run_aldex3_with_benchmarks(
        countdata        = counts_sub,
        conds            = conds_sub,
        scaledata        = truth_sub,
        scale_label      = truth_col,
        study_name       = study_name,
        checkpoint_dir   = checkpoint_dir,
        truthsamplesubset = truth_col
      )
      if (length(meas_res)) {
        meas_res <- purrr::imap(meas_res, \(df, .) {
          if (is.data.frame(df) && nrow(df))
            dplyr::mutate(df, n_samples = n_samples)
        })
      }
      final_results_list[[truth_col]][["measured"]] <-
        meas_res[names(meas_res) %in% c(paste0("ALDEx3_", truth_col))]

      # Optional: MLP blocks
      final_results_list[[truth_col]][["mlp"]] <- list()
      if (all(c("log2_mlp", "training_set", "profile") %in% colnames(df_pred))) {
        # Restrict to this truth's non-NA rows (by original df_pred rows), then keep only those with non-NA MLP
        pred_df <- df_pred[df_pred[[sample_id_col]] %in% rownames(conds_full), , drop = FALSE]
        pred_df <- pred_df[!is.na(pred_df$log2_mlp), , drop = FALSE]
        if (!nrow(pred_df)) {
          message("  Skipping MLP – no non-NA predictions in this truth subset.")
        } else {
          pred_df$sample_id <- pred_df[[sample_id_col]]
          message("DEBUG MLP: raw combos → ", paste0(sort(unique(pred_df$training_set)), collapse = ", "))
          pred_df %>%
            dplyr::group_by(training_set, profile) %>%
            dplyr::group_walk(~{
              combo <- .y; grp <- .x
              training_set <- combo$training_set; profile <- combo$profile

              # Named vector of MLP by sample ID, then index to covariate row order
              mlp_map <- setNames(grp$log2_mlp, grp$sample_id)
              mlp_all <- mlp_map[rownames(conds_full)]
              mlp_sub <- mlp_all[sample_ids_sub]   # exactly the same subset/order as counts_sub/conds_sub

              n_samples_combo <- length(mlp_sub)
              message(" >> MLP (", training_set, " / ", profile, ") | ", n_samples_combo, " samples <<")

              scale_lab <- paste0("Nishijima2024_", make.names(training_set), "_", make.names(profile))

              if (length(mlp_sub) != length(colnames(counts_sub)) ||
                  length(mlp_sub) != nrow(conds_sub)) {
                message("MLP predictions have the wrong number of samples for ",
                        study_name, " (", training_set, " / ", profile, ").")
                return(invisible(NULL))
              }

              mlp_res <- run_aldex3_with_benchmarks(
                countdata        = counts_sub,
                conds            = conds_sub,
                scaledata        = mlp_sub,
                scale_label      = scale_lab,
                study_name       = study_name,
                checkpoint_dir   = checkpoint_dir,
                truthsamplesubset = truth_col
              )
              if (length(mlp_res)) {
                mlp_res <- purrr::imap(mlp_res, \(df, nm) {
                  if (is.data.frame(df) && nrow(df)) df %>% dplyr::mutate(n_samples = n_samples_combo) else df
                })
              }
              matched_keys <- grep(scale_lab, names(mlp_res), value = TRUE)
              if (length(matched_keys) == 1) {
                savingmlpresultname <- paste0("ALDEx3_", scale_lab)
                final_results_list[[truth_col]][["mlp"]][[savingmlpresultname]] <<- mlp_res[[matched_keys]]
              } else if (length(matched_keys) > 1) {
                stop("Multiple matches found for scale_lab = ", scale_lab, ":\n", paste(matched_keys, collapse = ", "))
              } else {
                stop("No match found for scale_lab = ", scale_lab)
              }
            })
        }
      } else {
        message("DEBUG: MLP block skipped – columns ‘log2_mlp’, ‘training_set’, or ‘profile’ missing.")
      }
    }
  }

  all_empty <- all(purrr::map_lgl(final_results_list, \(lst) {
    all(purrr::map_lgl(lst, \(df) is.data.frame(df) && !nrow(df)))
  }))
  if (all_empty) {
    message("No non-empty benchmarking results for ", study_name, ". Returning NULL.")
    return(NULL)
  }
  return(final_results_list)
}


flatten_aldex_results <- function(aldex_out, countdata, study_name,
                                  scale_label = "ZSM", use = c("mean", "median")) {
  use <- match.arg(use)
  results <- lapply(names(aldex_out), function(gamma) {
    res <- aldex_out[[gamma]]
    est_matrix <- switch(
      use,
      mean   = apply(res$estimate, c(1, 2), mean),
      median = apply(res$estimate, c(1, 2), median)
    )
    est_matrix <- est_matrix[-1, , drop = FALSE]
    pval_matrix <- res$p.val[-1, , drop = FALSE]
    padj_matrix <- res$p.val.adj[-1, , drop = FALSE]
    covariate_names <- rownames(res$X)[-1]
    feature_names <- rownames(countdata)
    rownames(est_matrix)   <- covariate_names
    rownames(pval_matrix)  <- covariate_names
    rownames(padj_matrix)  <- covariate_names
    colnames(est_matrix)   <- feature_names
    colnames(pval_matrix)  <- feature_names
    colnames(padj_matrix)  <- feature_names
    estimate_df <- as.data.frame(est_matrix) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "estimate")
    pval_df <- as.data.frame(pval_matrix) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "p.val")
    padj_df <- as.data.frame(padj_matrix) %>%
      dplyr::mutate(covariate = covariate_names) %>%
      tidyr::pivot_longer(-covariate, names_to = "feature", values_to = "p.val.adj")
    out_df <- estimate_df %>%
      dplyr::left_join(pval_df, by = c("covariate", "feature")) %>%
      dplyr::left_join(padj_df, by = c("covariate", "feature"))
    metadata_df <- res$data %>%
      dplyr::mutate(feature = rownames(res$data)) %>%
      dplyr::relocate(feature)
    out_df <- out_df %>%
      dplyr::left_join(metadata_df, by = "feature") %>%
      dplyr::mutate(
        gamma       = as.numeric(gamma),
        study_name  = study_name,
        scale_label = scale_label
      )
    return(out_df)
  })
  dplyr::bind_rows(results)
}

aldex3_conf_matrix_by_study <- function(nested_df_list, true_gamma = 0.5, compare_models = NULL, alpha = 0.05) {
  results <- list()
  for (study in names(nested_df_list)) {
    study_list <- nested_df_list[[study]]
    truth_columns <- names(study_list)
    results[[study]] <- list()
    for (truth in truth_columns) {
      methods_group <- study_list[[truth]]
      rows <- list()
      for (group in names(methods_group)) {
        for (method_name in names(methods_group[[group]])) {
          df0 <- as.data.frame(methods_group[[group]][[method_name]])
          if (nrow(df0) == 0) next
          df0 <- df0 %>% mutate(gamma = as.numeric(as.character(gamma)))
          if (group == "measured") {
            df0 <- df0 %>% filter(gamma == true_gamma)
            if (nrow(df0) == 0) next
          }
          req <- c("feature","covariate","p.val.adj","estimate","scale_label","gamma")
          if (!all(req %in% colnames(df0))) next
          df1 <- df0 %>%
            mutate(
              Study = study,
              CombinedModel = paste0(scale_label, "_", gamma)
            ) %>%
            select(Study, covariate, feature, estimate, p.val.adj,
                   scale_label, gamma, CombinedModel)
          rows[[paste(group, method_name, sep="__")]] <- df1
        }
      }
      if (length(rows) == 0) next
      df_flat <- bind_rows(rows)
      true_model_name <- paste0(truth, "_", true_gamma)
      gold <- df_flat %>%
        filter(CombinedModel == true_model_name) %>%
        select(feature, covariate,
               gold_pval = p.val.adj, gold_sign = estimate) %>%
        distinct(feature, covariate, .keep_all = TRUE) %>%
        mutate(gold_pval = ifelse(is.na(gold_pval), 1, gold_pval))
      conf_cats <- list()
      if (nrow(gold) > 0) {
        model_ids <- setdiff(unique(df_flat$CombinedModel), true_model_name)
        for (m in model_ids) {
          md <- df_flat %>%
            filter(CombinedModel == m) %>%
            select(feature, covariate,
                   model_pval = p.val.adj, model_sign = estimate) %>%
            distinct(feature, covariate, .keep_all = TRUE) %>%
            mutate(model_pval = ifelse(is.na(model_pval), 1, model_pval))
          joined <- inner_join(md, gold, by = c("feature","covariate"))
          if (nrow(joined) == 0) next
          sig_pred  <- joined$model_pval <= alpha
          sig_true  <- joined$gold_pval  <= alpha
          same_sign <- sign(joined$model_sign) == sign(joined$gold_sign)
          TP <- sum(sig_true & sig_pred & same_sign)
          TN <- sum(!sig_true & !sig_pred)
          FP <- sum((sig_pred & !sig_true) | (sig_pred & sig_true & !same_sign))
          FN <- sum(sig_true & !sig_pred)
          metrics <- tibble::tibble(
            TP          = TP,
            TN          = TN,
            FP          = FP,
            FN          = FN,
            PPV         = if_else((TP + FP) > 0, TP  / (TP + FP) * 100, NA_real_),
            NPV         = if_else((TN + FN) > 0, TN  / (TN + FN) * 100, NA_real_),
            Sensitivity = if_else((TP + FN) > 0, TP  / (TP + FN) * 100, NA_real_),
            FDR         = if_else((TP + FP) > 0, FP  / (TP + FP) * 100, 0)
          )
          conf_cats[[m]] <- metrics
        }
      }
      export_rows <- purrr::imap_dfr(
        conf_cats,
        ~ data.frame(Study = study, Model = .y, conf_cats[[.y]], stringsAsFactors = FALSE)
      )
      export_table <- tibble::tibble(
        Study   = study,
        metrics = list(export_rows)
      )
      df_val <- export_rows %>%
        pivot_longer(
          cols = c(TP, TN, FP, FN, PPV, NPV, Sensitivity, FDR),
          names_to  = "Metric",
          values_to = "Value"
        )
      df_pct <- export_rows %>%
        rowwise() %>%
        mutate(
          total = TP + TN + FP + FN,
          TP    = if (total > 0) TP   / total * 100 else NA,
          TN    = if (total > 0) TN   / total * 100 else NA,
          FP    = if (total > 0) FP   / total * 100 else NA,
          FN    = if (total > 0) FN   / total * 100 else NA
        ) %>%
        ungroup() %>%
        select(-total) %>%
        pivot_longer(
          cols = c(TP, TN, FP, FN, PPV, NPV, Sensitivity, FDR),
          names_to  = "Metric",
          values_to = "Percent"
        )
      metrics_long <- df_val %>%
        left_join(df_pct, by = c("Study", "Model", "Metric"))

      if (!is.null(compare_models)) {
        metrics_long <- metrics_long %>% filter(Model %in% compare_models)
      }
      results[[study]][[truth]] <- list(
        conf_matrices = conf_cats,
        export_table  = export_table,
        metrics_long  = metrics_long
      )
    }
  }
  combined_metrics_long <- purrr::imap_dfr(
    results,
    function(study_list, study_name) {
      purrr::imap_dfr(
        study_list,
        function(truth_list, truth_name) {
          truth_list$metrics_long %>%
            mutate(Study = study_name, Truth = truth_name)
        }
      )
    }
  )
  combined_metrics_long <- combined_metrics_long %>%
  extract(
    col = Model,
    into = c("model_type", "gamma"),
    regex = "^(.*?)(?:_([0-9.]+))?$",  
    remove = FALSE
  ) 
  return(list(studyresults = results, longresults = combined_metrics_long))
}

aldex_boxplot <- function(data, gamma_level = NULL, stat = FALSE, output_dir = ".", split_by_seqtype = FALSE, split_by_loadtype = FALSE,
                          metrics = c("PPV", "NPV", "FDR"),  highlight_studies = NULL, color_by_study = TRUE, 
                          color_by_loadtype = FALSE, add_boxplot = TRUE) {
  
  require(ggsignif)
  require(purrr)
  require(glue)
  require(scales)
  
  data <- data %>%
    mutate(model_type = str_replace_all(model_type, "_", " "))
  
  data <- if (is.null(gamma_level)) {
    data %>% filter(
      str_detect(model_type, "^ALDEx3") |
      (!str_detect(model_type, "^ALDEx3") & (is.na(gamma) | gamma == 0))
    )
  } else {
    data %>% filter(
      (str_detect(model_type, "^ALDEx3") & gamma == gamma_level) |
      (!str_detect(model_type, "^ALDEx3") & (is.na(gamma) | gamma == 0))
    )
  }
  
  tidyplot_theme <- function(base_size = 18) {
    theme_minimal(base_size = base_size) +
      theme(
        text              = element_text(family = "Arial", color = "black"),
        axis.text         = element_text(size = 16, color = "black"),
        axis.line         = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks        = element_line(size = 0.5, color = "black"),
        panel.grid        = element_blank(),
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        legend.title      = element_text(size = 18, face = "bold"),
        legend.text       = element_text(size = 16),
        plot.subtitle     = element_blank(),
        plot.title = element_blank()
      )
  }

  make_plot <- function(df_subset, panel_label = NULL, add_boxplot = TRUE) {
    df_subset <- df_subset %>%
      filter(Metric %in% metrics, !is.na(Value)) %>%
      mutate(
        model_label = case_when(
          str_detect(model_type, regex("ALDEx3 tss", ignore_case=TRUE)) ~ glue("ALDEx3 TSS (γ={gamma})"),
          str_detect(model_type, regex("^limma voom",     ignore_case=TRUE)) ~ "Limma Voom",
          str_detect(model_type, regex("^metagenomeSeq \\(CSS\\+fitZig\\)", ignore_case=TRUE)) ~ "metagenomeSeq",
          str_detect(model_type, regex("Nishijima2024.*galaxy",     ignore_case=TRUE)) ~ "Nishijima2024 GALAXY",
          str_detect(model_type, regex("Nishijima2024.*metacardis", ignore_case=TRUE)) ~ "Nishijima2024 MetaCardis",
          str_detect(model_type, regex("Nishijima2024.*Vandeputte2021", ignore_case=TRUE)) ~ "Nishijima2024 Vandeputte2021",
          TRUE ~ model_type
        )
      ) %>%
      filter(
        !(sequencingtype == "shotgun metagenomics" & model_label == "Nishijima2024 Vandeputte2021"),
        !(sequencingtype == "16s rrna"              & model_label %in% c("Nishijima2024 GALAXY","Nishijima2024 MetaCardis"))
      ) 
    
    nishi <- if (is.null(panel_label)) {
      c("Nishijima2024 GALAXY","Nishijima2024 MetaCardis","Nishijima2024 Vandeputte2021")
    } else if (panel_label %in% unique(data$sequencingtype)) {
      if (panel_label == "shotgun metagenomics") c("Nishijima2024 GALAXY","Nishijima2024 MetaCardis")
      else if (panel_label == "16s rrna")        "Nishijima2024 Vandeputte2021"
      else                                       character(0)
    } else {
      c("Nishijima2024 GALAXY","Nishijima2024 MetaCardis","Nishijima2024 Vandeputte2021")
    }
    
    base_models <- c("ANCOM-BC2","DESeq2","Limma Voom", "edgeR","metagenomeSeq", "LinDA")
    tss_names   <- unique(df_subset$model_label[grepl("^ALDEx3 TSS", df_subset$model_label)])
    extract_gamma <- function(lbl) as.numeric(str_extract(lbl, "(?<=γ=)[0-9.]+"))
    tss_names <- tss_names[order(sapply(tss_names, extract_gamma))]
    desired <- c(base_models, tss_names, nishi)
    
    cols <- c(GALAXY = "#0C58CA", MetaCardis = "#FF8C27", Vandeputte2021 = "#1B1B1B")
    levels_html <- vapply(desired, function(lbl) {
      if (lbl == "ALDEx3 TSS (γ=0)") {
        "<strong>ALDEx3 TSS (γ=0)</strong>"
      } else if (grepl("GALAXY", lbl, fixed = TRUE)) {
        paste0("Nishijima2025 <span style='color:", cols["GALAXY"], 
               "'><strong><em>GALAXY</em></strong></span>")
      } else if (grepl("MetaCardis", lbl, fixed = TRUE)) {
        paste0("Nishijima2025 <span style='color:", cols["MetaCardis"], 
               "'><strong><em>MetaCardis</em></strong></span>")
      } else if (grepl("Vandeputte2021", lbl, fixed = TRUE)) {
        paste0("Nishijima2025 <span style='color:", cols["Vandeputte2021"], 
               "'><strong><em>Vandeputte2021</em></strong></span>")
      } else {
        lbl
      }
    }, character(1))
    names(levels_html) <- desired
    
    df_subset <- df_subset %>%
      mutate(label_html = factor(levels_html[model_label], levels = levels_html))
    
    present_labels <- unique(df_subset$model_label)
    all_others    <- c("ALDEx3 TSS (γ=0.5)", "ALDEx3 TSS (γ=1)", nishi)
    want_labels   <- intersect(all_others, present_labels)
    want_html     <- levels_html[want_labels]
    baseline_plain <- "ALDEx3 TSS (γ=0)"
    baseline_html  <- levels_html[baseline_plain]
    my_comparisons <- lapply(want_html, function(x_html) c(baseline_html, x_html))
    my_comparisons <- lapply(my_comparisons, unname)
    
    studies     <- sort(unique(df_subset$Study))
    study_cols  <- setNames(hue_pal()(length(studies)), studies)
    if (!is.null(highlight_studies)) study_cols[setdiff(studies, highlight_studies)] <- "grey80"
    study_breaks <- if (!is.null(highlight_studies)) highlight_studies else studies
    
    lts       <- sort(unique(df_subset$loadtype))
    load_cols <- setNames(hue_pal()(length(lts)), lts)
    load_breaks <- lts
    n_met      <- length(unique(df_subset$Metric))
    facet_dims <- if (n_met <= 3) c(3,1) else if (n_met == 4) c(2,2) else c(3, ceiling(n_met/3))

    df_subset <- df_subset %>%
      group_by(label_html, Metric) %>%
      mutate(
        Q1        = quantile(Value, .25),
        Q3        = quantile(Value, .75),
        IQR       = Q3 - Q1,
        lower     = Q1 - 1.5 * IQR,
        upper     = Q3 + 1.5 * IQR,
        is_outlier = Value < lower | Value > upper
      ) %>%
      ungroup()

    if (stat) {
    sig_data <- df_subset %>%
      group_by(Metric) %>%
      group_modify(~ {
        y_max <- max(.x$Value, na.rm = TRUE)
        y_range <- diff(range(.x$Value, na.rm = TRUE))
        map_dfr(seq_along(my_comparisons), function(i) {
          comp <- my_comparisons[[i]]
          group1 <- comp[1]
          group2 <- comp[2]
          paired_data <- .x %>%
            filter(label_html %in% c(group1, group2)) %>%
            select(Study, loadtype, label_html, Value) %>%
            pivot_wider(
              names_from = label_html, 
              values_from = Value,
              values_fn = mean
            )
          if (!all(c(group1, group2) %in% colnames(paired_data))) {
            return(tibble())
          }
          test <- tryCatch(
            wilcox.test(
              paired_data[[group1]], 
              paired_data[[group2]], 
              paired = TRUE,
              exact = FALSE
            ),
            error = function(e) list(p.value = NA_real_)
          )
          tibble(
            group1 = group1,
            group2 = group2,
            p_value = test$p.value,
            y_pos = y_max + (0.07 * y_range * i) 
          )
        })
      }) %>%
      ungroup() %>%
      mutate(
        label = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(label)) 
    }

    p <- ggplot(df_subset, aes(x = label_html, y = Value)) +
      { if (color_by_loadtype) {
          geom_jitter(aes(color = loadtype), width = 0.2, size = 1.5, shape = 16, stroke = 1)
        } else if (color_by_study) {
          geom_jitter(aes(color = Study), width = 0.2, size = 1.5, shape = 16, stroke = 1)
        } else {
          list(
            geom_jitter(
              data = filter(df_subset, !is_outlier),
              color = "grey80", width = 0.2, size = 1.5, shape = 16
            ),
            geom_jitter(
              data = filter(df_subset, is_outlier),
              color = "black", size = 1.5, width = 0.1, shape = 23
            )
          )
        }
      } +
      { if (add_boxplot) {
          list(
            geom_boxplot(color = "black", fill = "white", size = 0.8, 
                        width = 0.5, alpha = 0.5, outlier.shape = NA),
            stat_boxplot(
              geom = "errorbar", width = 0.5, size = 0.7,
              aes(ymin = after_stat(ymax), ymax = after_stat(ymax))
            ),
            stat_boxplot(
              geom = "errorbar", width = 0.5, size = 0.7,
              aes(ymin = after_stat(ymin), ymax = after_stat(ymin))
            )
          )
        }
      } +
      { if (stat && nrow(sig_data) > 0) {
          ggsignif::geom_signif(  
            data = sig_data,
            aes(xmin = group1, xmax = group2, annotations = label, y_position = y_pos),
            manual = TRUE,
            textsize = 4, 
            tip_length = 0.01, 
            vjust = 0.5  
          )
        }
      } +
      { if (color_by_loadtype) {
          scale_color_manual(name = "Loadtype", values = load_cols, breaks = load_breaks,
                             guide = guide_legend(override.aes = list(size = 4)))
        } else if (color_by_study) {
          scale_color_manual(name = "Study", values = study_cols, breaks = study_breaks,
                             guide = guide_legend(override.aes = list(size = 4)))
        } 
      } +
      scale_x_discrete(drop = FALSE, expand = expansion(add = c(0.8,0.8))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + 
      tidyplot_theme() +
      theme(
        
        plot.margin = margin(t = 50, r = 5, b = 5, l = 5),
        axis.text.x      = element_markdown(size = 16, angle = 45, hjust = 1, vjust = 1, color = "black"),
        strip.placement  = "outside",
        panel.grid.major.y   = element_blank(),
        panel.grid.major.x   = element_blank(),
        strip.text.y     = element_text(size = 20, face = "bold"),
        axis.title.x     = element_text(size = 20, face = "bold", margin = margin(t = 10)),
        axis.title.y     = element_text(size = 20, face = "bold", margin = margin(r = 10)),
        panel.spacing.y  = unit(4, "lines"),
        plot.subtitle    = element_blank(),
        plot.title       = element_blank(),
        legend.title        = element_blank(),
        legend.key.size     = unit(1, "lines"),
        legend.text         = element_text(size = 6, color = "black", face = "bold"),
        legend.spacing.x    = unit(0.2, "cm"),
        legend.spacing.y    = unit(0.15, "cm"),
        legend.box.margin   = margin(5, 5, 5, 5),
        legend.box.spacing  = unit(0.3, "cm"),
        legend.position  = ifelse(color_by_loadtype | color_by_study, "top", "none"),
        legend.box       = "vertical",
        legend.justification = "left",
        legend.margin      = margin(l = -0.2, unit = "npc")
      ) +
      labs(x = NULL, y = "Percentage (%)") +
      facet_wrap(~ Metric, ncol = facet_dims[2], nrow = facet_dims[1], scales = "free_y", switch = "y")
    
    p <- p +
      facetted_pos_scales(
        y = list(
          Metric == "NPV" ~ scale_y_continuous(
            limits = if(any(df_subset$sequencingtype == "16s rrna")) c(65, 100) else c(92.5, 100),
            expand = expansion(mult = c(0.05, 0.05))
          ),
          Metric == "PPV" ~ scale_y_continuous(
            limits = c(NA, 100),
            expand = expansion(mult = c(0.05, 0.05))
          ),
          Metric == "FDR" ~ scale_y_continuous(
            limits = c(NA, 100),
            expand = expansion(mult = c(0.05, 0.05))
          )
        )
      ) +
      coord_cartesian(clip = "off")
    
    fname <- ifelse(
      is.null(panel_label),
      "aldex_boxplot_all",
      paste0("aldex_boxplot_", str_replace_all(tolower(panel_label), " ", "_"))
    )
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 6, device = cairo_pdf)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 6, dpi = 300)
    
    return(p)
  }
  
  if (split_by_seqtype && split_by_loadtype) {
    seqs <- unique(data$sequencingtype)
    lts  <- unique(data$loadtype)
    plots <- list()
    for (st in seqs) {
      for (lt in lts) {
        df_sub <- filter(data, sequencingtype == st, loadtype == lt)
        panel_label <- paste0(st, " | ", lt)
        key <- paste0(st, "__", lt)
        plots[[key]] <- make_plot(df_sub, panel_label, add_boxplot)
      }
    }
    return(plots)
  }
  
  if (split_by_loadtype) {
    lts   <- unique(data$loadtype)
    plots <- lapply(lts, function(lt) {
      make_plot(filter(data, loadtype == lt), panel_label = lt, add_boxplot)
    })
    names(plots) <- lts
    return(plots)
  }
  
  if (split_by_seqtype) {
    sts   <- unique(data$sequencingtype)
    plots <- lapply(sts, function(st) {
      make_plot(filter(data, sequencingtype == st), panel_label = st, add_boxplot)
    })
    names(plots) <- sts
    return(plots)
  }
  
  make_plot(data, panel_label = NULL, add_boxplot)
}

compute_clr <- function(counts) {
  closure <- t(t(counts + 0.5) / colSums(counts + 0.5))
  log_closure <- log2(closure)
  -colMeans(log_closure)
}

lm.scale <- function(Y, X, data = NULL, scale = NULL, p.adjust.method = "BH") {
  N <- ncol(Y)
  if (is.null(scale)) {
    closure <- t(t(Y + 0.5) / colSums(Y + 0.5))
    log_closure <- log2(closure)
    Z <- -colMeans(log_closure)
  } else {
    if (length(scale) != N) stop("Length of 'scale' vector must match number of samples (ncol(Y)).")
    Z <- scale
  }
  Z <- mean_center(Z)
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("Data frame required when X is a formula.")
    X_mat <- model.matrix(X, data)
  } else if (is.matrix(X) || is.data.frame(X)) {
    X_mat <- as.matrix(X)
    if (nrow(X_mat) != N) stop("Rows of X must equal number of samples.")
  } else {
    stop("X must be a formula or a matrix/data.frame.")
  }
  df_lm <- as.data.frame(X_mat, check.names = TRUE)
  df_lm$Z <- Z
  fit <- stats::lm(Z ~ . - 1, data = df_lm)
  s <- summary(fit)
  ss    <- s$coefficients
  coefs <- setNames(ss[, "Estimate"],    rownames(ss))
  ses   <- setNames(ss[, "Std. Error"], rownames(ss))
  tvals <- setNames(ss[, "t value"],    rownames(ss))
  pvals <- setNames(ss[, "Pr(>|t|)"],   rownames(ss))
  ci <- tryCatch({
    ci_mat <- confint(fit)
    ci_mat
  }, error = function(e) {
    matrix(NA, nrow = length(coefs), ncol = 2)
  })
  ci_low  <- setNames(ci[,1], rownames(ss))
  ci_high <- setNames(ci[,2], rownames(ss))
  adj_r2  <- s$adj.r.squared
  aic_val <- AIC(fit)
  bic_val <- BIC(fit)
  fstat   <- s$fstatistic
  f_pval  <- if (!is.null(fstat)) pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) else NA
  list(
    estimate    = coefs,
    std.error   = ses,
    statistic   = tvals,
    p.val       = pvals,
    p.val.adj   = p.adjust(pvals, method = p.adjust.method),
    ci.low      = ci_low,
    ci.high     = ci_high,
    adj.r.squared = adj_r2,
    aic           = aic_val,
    bic           = bic_val,
    f.statistic   = if (!is.null(fstat)) fstat[1] else NA,
    f.df1         = if (!is.null(fstat)) fstat[2] else NA,
    f.df2         = if (!is.null(fstat)) fstat[3] else NA,
    f.pvalue      = f_pval,
    n_coef        = length(coefs)
  )
}

evaluate_predictions <- function(truth_vec, pred_vec, conds_df, p.adjust.method = "BH",conf.level = 0.95) {
  z_truth <- mean_center(truth_vec)
  z_pred  <- mean_center(pred_vec)
  resid   <- z_truth - z_pred
  df <- as.data.frame(conds_df, stringsAsFactors = FALSE)
  df$resid <- resid
  fit_res <- lm(resid ~ ., data = df)
  s_res   <- summary(fit_res)
  coef_mat <- s_res$coefficients
  vars     <- rownames(coef_mat)
  ci_mat <- tryCatch(
    confint(fit_res, level = conf.level),
    error = function(e) {
      matrix(NA, nrow = length(vars), ncol = 2,
             dimnames = list(vars, c("2.5 %", "97.5 %")))
    }
  )
  coef_df <- data.frame(
    variable   = vars,
    estimate   = coef_mat[, "Estimate"],
    std.error  = coef_mat[, "Std. Error"],
    statistic  = coef_mat[, "t value"],
    p.value    = coef_mat[, "Pr(>|t|)"],
    p.adj      = p.adjust(coef_mat[, "Pr(>|t|)"], method = p.adjust.method),
    ci.low     = ci_mat[, 1],
    ci.high    = ci_mat[, 2],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  fstat <- s_res$fstatistic
  metrics_df <- data.frame(
    r.squared     = s_res$r.squared,
    adj.r.squared = s_res$adj.r.squared,
    aic           = AIC(fit_res),
    bic           = BIC(fit_res),
    f.statistic   = if (!is.null(fstat)) fstat[1] else NA_real_,
    f.df1         = if (!is.null(fstat)) fstat[2] else NA_real_,
    f.df2         = if (!is.null(fstat)) fstat[3] else NA_real_,
    f.pvalue      = if (!is.null(fstat))
                       pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
                     else NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  list(
    coefficients = coef_df,
    metrics      = metrics_df
  )
}

run_lm_scale <- function(countdata, conds, scaledata = NULL, scale_label = NULL, study_name = NULL) {
  conds <- conds[, setdiff(names(conds), "SampleID"), drop = FALSE]
  covariates <- make.names(names(conds))
  colnames(conds) <- covariates
  modelformula <- as.formula(paste("~", paste(covariates, collapse=" + ")))
  fit_out <- tryCatch({
    lm.scale(
      Y     = countdata,
      X     = modelformula,
      data  = conds,
      scale = scaledata
    )
  }, error = function(e) {
    message("Error in lm.scale: ", e$message)
    return(NULL)
  })
  if (is.null(fit_out)) {
    return(list(coefficients = NULL, metrics = NULL))
  }
  vars   <- names(fit_out$estimate)
  est    <- fit_out$estimate[vars]
  ses    <- fit_out$std.error[vars]
  tval   <- fit_out$statistic[vars]
  pval   <- fit_out$p.val[vars]
  padj   <- fit_out$p.val.adj[vars]
  ci_low  <- fit_out$ci.low[vars]
  ci_high <- fit_out$ci.high[vars]
  coef_df <- data.frame(
    variable      = vars,
    estimate.mean = est,
    std.error     = ses,
    statistic     = tval,
    p.val         = pval,
    p.val.adj     = padj,
    ci.low        = ci_low,
    ci.high       = ci_high,
    study_name    = study_name,
    scale_label   = scale_label,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  var_clean <- make.names( gsub("`", "", coef_df$variable) )
  covs_sorted <- covariates[order(-nchar(covariates))]
  coef_df$BaseCov <- vapply(
    var_clean,
    FUN.VALUE = character(1),
    function(v) {
      hits <- covs_sorted[ startsWith(v, covs_sorted) ]
      if (length(hits)) hits[1] else NA_character_
    }
  )
  coef_df$BetaIndex <- match(coef_df$BaseCov, covariates)
  coef_df$BaseCov <- NULL
  metrics_df <- data.frame(
    adj.r.squared = fit_out$adj.r.squared,
    aic           = fit_out$aic,
    bic           = fit_out$bic,
    f.statistic   = fit_out$f.statistic,
    f.df1         = fit_out$f.df1,
    f.df2         = fit_out$f.df2,
    f.pvalue      = fit_out$f.pvalue,
    study_name    = study_name,
    scale_label   = scale_label,
    stringsAsFactors = FALSE
  )
  list(
    coefficients = coef_df,
    metrics      = metrics_df
  )
}

process_study_with_lm_scale <- function(study_name, study_data) {
  message("Analyzing: ", study_name)
  id_col    <- study_data$studydemographics$ID
  df_meta   <- study_data$metadata
  df_pred   <- study_data$predictions
  counts_df <- study_data$counts
  sid_col <- grep("^sample[._ ]?ID$", colnames(df_pred),
                  ignore.case = TRUE, value = TRUE)
  if (length(sid_col) != 1) {
    stop("Need exactly one sample-ID column in predictions")
  }
  if (is.null(rownames(counts_df))) {
    rownames(counts_df) <- as.character(df_meta[[id_col]])
  } else if (any(rownames(counts_df) == "")) {
    miss <- which(rownames(counts_df) == "")
    rownames(counts_df)[miss] <- as.character(df_meta[[id_col]])[miss]
  }
  common_ids <- Reduce(intersect, list(
    rownames(counts_df),
    as.character(df_meta[[id_col]]),
    as.character(df_pred[[sid_col]])
  ))
  if (!length(common_ids)) {
    message("No overlapping samples; skipping ", study_name)
    return(NULL)
  }
  counts_df <- counts_df[common_ids, , drop = FALSE]
  df_meta   <- df_meta[match(common_ids, df_meta[[id_col]]), , drop = FALSE]
  df_pred   <- df_pred[df_pred[[sid_col]] %in% common_ids, , drop = FALSE]
  counts_mat <- t(as.matrix(counts_df))
  storage.mode(counts_mat) <- "integer"
  cov_cols <- study_data$studydemographics$covariates
  conds_df  <- df_meta[, cov_cols, drop = FALSE]
  structural_results <- list()
  evaluation_results <- list()
  truth_cols <- grep("^log2_truth_", colnames(df_pred), value = TRUE)
  for (tc in truth_cols) {
    has_truth  <- !is.na(df_pred[[tc]])
    truth_sids <- intersect(unique(df_pred[[sid_col]][has_truth]), common_ids)
    if (!length(truth_sids)) next
    sub_counts <- counts_mat[, truth_sids, drop = FALSE]
    sub_meta   <- df_meta[df_meta[[id_col]] %in% truth_sids, , drop = FALSE]
    sub_meta   <- sub_meta[match(truth_sids, sub_meta[[id_col]]), , drop = FALSE]
    sub_conds  <- sub_meta[, cov_cols, drop = FALSE]
    sub_pred   <- df_pred[df_pred[[sid_col]] %in% truth_sids, , drop = FALSE]
    truth_vec  <- sub_pred[[tc]][match(truth_sids, sub_pred[[sid_col]])]
    if (length(truth_vec) != ncol(sub_counts)) {
      stop("Truth length (", length(truth_vec),
           ") != #samples (", ncol(sub_counts), ") for ", tc)
    }
    keep <- complete.cases(sub_conds)
    sub_counts <- sub_counts[, keep, drop = FALSE]
    sub_conds  <- sub_conds[keep, , drop = FALSE]
    sub_pred   <- sub_pred[keep, , drop = FALSE]
    valid_ids  <- colnames(sub_counts)
    sub_pred  <- sub_pred[sub_pred[[sid_col]] %in% valid_ids, , drop = FALSE]
    truth_vec <- sub_pred[[tc]][match(valid_ids, sub_pred[[sid_col]])]
    keep_cov <- vapply(sub_conds, function(x) length(unique(x)) >= 2, logical(1))
    sub_conds <- sub_conds[, keep_cov, drop = FALSE]
    if (!ncol(sub_conds)) next
    clr_vec <- compute_clr(sub_counts)
    tss_vec <- colSums(sub_counts)
    clr_fit   <- run_lm_scale(sub_counts, sub_conds, clr_vec,
                              scale_label = paste0(tc, "_CLR"),
                              study_name  = study_name)
    tss_fit   <- run_lm_scale(sub_counts, sub_conds, tss_vec,
                              scale_label = paste0(tc, "_TSS"),
                              study_name  = study_name)
    truth_fit <- run_lm_scale(sub_counts, sub_conds, truth_vec,
                              scale_label = tc,
                              study_name  = study_name)
    eval_clr <- evaluate_predictions(truth_vec, clr_vec, sub_conds)
    eval_tss <- evaluate_predictions(truth_vec, tss_vec, sub_conds)
    mlp_struct   <- list()
    mlp_evaluate <- list()
    if (all(c("log2_mlp","training_set","profile") %in% colnames(sub_pred))) {
      for (tr in unique(sub_pred$training_set)) {
        for (pr in unique(sub_pred$profile[sub_pred$training_set == tr])) {
          label <- paste0(tr,"_",pr)
          idx   <- which(sub_pred$training_set==tr & sub_pred$profile==pr)
          mlp_vec <- sub_pred$log2_mlp[idx]
          mlp_struct[[label]]   <- run_lm_scale(sub_counts, sub_conds, mlp_vec,
                                               scale_label = label,
                                               study_name  = study_name)
          mlp_evaluate[[label]] <- evaluate_predictions(
            truth_vec, mlp_vec, sub_conds
          )
        }
      }
    }
    structural_results[[tc]] <- list(
      CLR   = clr_fit,
      TSS   = tss_fit,
      Truth = truth_fit,
      MLP   = mlp_struct
    )
    evaluation_results[[tc]] <- c(
      list(CLR   = eval_clr,
           TSS   = eval_tss),
      mlp_evaluate
    )
  }
  list(
    structural_results = structural_results,
    residual_regression_results = evaluation_results
  )
}

modelingscale <- function(nested_truth_list, aldexrepo = NULL, model_include = NULL) {  
  capitalize_first <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  flat_coef_list <- list()
  flat_metrics_list <- list()
  for(study in names(nested_truth_list)) {
    struct_results <- nested_truth_list[[study]]$structural_results 
    for(truth in names(struct_results)) {
      models <- struct_results[[truth]]
      for(model_nm in setdiff(names(models), "MLP")) {
        model_data <- models[[model_nm]]
        
        if(!is.null(model_data$coefficients)) {
          df_coef <- model_data$coefficients
          df_coef$Study <- study
          df_coef$Truth <- truth
          df_coef$Model <- model_nm
          df_coef$ModelType <- model_nm
          flat_coef_list[[length(flat_coef_list) + 1]] <- df_coef
        }
        if(!is.null(model_data$metrics)) {
          df_metrics <- model_data$metrics
          df_metrics$Study <- study
          df_metrics$Truth <- truth
          df_metrics$Model <- model_nm
          df_metrics$ModelType <- model_nm
          flat_metrics_list[[length(flat_metrics_list) + 1]] <- df_metrics
        }
      }
      if("MLP" %in% names(models)) {
        mlp_models <- models[["MLP"]]
        for(mlp_name in names(mlp_models)) {
          mlp_model <- mlp_models[[mlp_name]]
          if(!is.null(mlp_model$coefficients)) {
            df_coef <- mlp_model$coefficients
            df_coef$Study <- study
            df_coef$Truth <- truth
            df_coef$Model <- mlp_name
            df_coef$ModelType <- "MLP"
            flat_coef_list[[length(flat_coef_list) + 1]] <- df_coef
          }
          if(!is.null(mlp_model$metrics)) {
            df_metrics <- mlp_model$metrics
            df_metrics$Study <- study
            df_metrics$Truth <- truth
            df_metrics$Model <- mlp_name
            df_metrics$ModelType <- "MLP"
            flat_metrics_list[[length(flat_metrics_list) + 1]] <- df_metrics
          }
        }
      }
    }
  }
  if(length(flat_coef_list) == 0) {
    warning("No data to process.")
    return(list(export_table = NULL))
  }
  big_coef_df <- bind_rows(flat_coef_list)
  big_metrics_df <- if(length(flat_metrics_list) > 0) bind_rows(flat_metrics_list) else NULL
  required_coef <- c("Study", "Truth", "Model", "variable", "estimate.mean", "p.val")
  missing_coef <- setdiff(required_coef, colnames(big_coef_df))
  if(length(missing_coef)) {
    stop("Missing required column(s) in coefficients: ", 
         paste(missing_coef, collapse = ", "))
  }
  final_coef_df <- big_coef_df %>%
    dplyr::rename(Beta = estimate.mean, pval = p.val) %>%
    filter(!str_detect(variable, regex("intercept", ignore_case = TRUE))) 
  study_covs <- map(
    names(nested_truth_list),
    ~ {
      sd <- aldexrepo[[.x]]$studydemographics
      if(!is.null(sd$covariates)) sd$covariates else character(0)
    }
  ) %>% set_names(names(nested_truth_list))
  final_coef_df <- final_coef_df %>%
    rowwise() %>%
    mutate(
      var_base = {
        covs <- study_covs[[Study]]
        hits <- covs[str_detect(variable, fixed(covs))]
        if(length(hits) > 0) hits[1] else NA_character_
      }
    ) %>%
    ungroup()
  factor_lookup <- final_coef_df %>%
    distinct(Study, Truth, var_base, variable) %>%
    group_by(Study, Truth, var_base) %>%
    arrange(variable) %>%
    mutate(factor_index = row_number()) %>%
    ungroup()
  final_coef_df <- final_coef_df %>%
    left_join(factor_lookup, by = c("Study", "Truth", "var_base", "variable")) %>%
    mutate(
      CoefName = paste0("B", BetaIndex, " (", factor_index, ")"),
      CoefNameFull = factor(
        paste0(Study, "__", Truth, "__", CoefName),
        levels = unique(paste0(Study, "__", Truth, "__", CoefName))
      )
    )
  measured_df <- final_coef_df %>%
    filter(ModelType == "Truth") %>%
    select(Study, Truth, variable, MeasuredBeta = Beta, MeasuredPval = pval)
  plot_coef_df <- final_coef_df %>%
    left_join(measured_df, by = c("Study", "Truth", "variable"))
  if(!is.null(model_include)) {
    plot_coef_df <- filter(plot_coef_df, Model %in% model_include)
  }
  combined_df <- plot_coef_df
  if(!is.null(big_metrics_df)) {
    plot_coef_df$model_id <- paste(plot_coef_df$Study, plot_coef_df$Truth, 
                                   plot_coef_df$Model, sep = "::")
    big_metrics_df$model_id <- paste(big_metrics_df$Study, big_metrics_df$Truth, 
                                     big_metrics_df$Model, sep = "::")
    combined_df <- left_join(
      plot_coef_df,
      big_metrics_df %>% select(-Study, -Truth, -Model, -scale_label),
      by = "model_id"
    ) %>% select(-model_id)
  }
  combined_df <- combined_df %>%
    mutate(
      model_color = case_when(
        Model %in% c("CLR", "ZSM", "Zero Scale") ~ "#781690",
        str_detect(Model, "vandeputte") ~ "black",
        str_detect(Model, "GALAXY") ~ "#0C58CA",
        str_detect(Model, "MetaCardis") ~ "#FF8C27",
        TRUE ~ "grey50"
      ),
      sig_model = pval < 0.05,
      sig_meas = MeasuredPval < 0.05,
      alpha_val = ifelse(sig_model & sig_meas, 1, 0.1)
    )
  models <- unique(combined_df$Model)
  offsets <- seq(-0.15, 0.15, length.out = length(models))
  off_map <- set_names(offsets, models)
  ordering <- combined_df %>%
    mutate(dir = MeasuredBeta - Beta) %>%
    group_by(Study, Truth, CoefNameFull) %>%
    summarise(mean_dir = mean(dir, na.rm = TRUE), .groups = "drop") %>%
    arrange(Study, Truth, desc(mean_dir)) %>%
    group_by(Study, Truth) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  combined_df <- combined_df %>%
    left_join(ordering %>% select(Study, Truth, CoefNameFull, rank),
              by = c("Study", "Truth", "CoefNameFull")) %>%
    arrange(Study, Truth, rank) %>%
    mutate(
      CoefNameFull = factor(CoefNameFull, levels = unique(CoefNameFull)),
      Model = capitalize_first(Model),
      x_jit = as.numeric(CoefNameFull) + off_map[Model]
    )
  facet_limits <- combined_df %>%
    group_by(Study, Truth) %>%
    summarise(
      y_max = max(abs(c(Beta, MeasuredBeta)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_min = -y_max)
  combined_df <- left_join(combined_df, facet_limits, by = c("Study", "Truth"))
  if (!is.null(aldexrepo)) {
    combined_df <- combined_df %>%
      rowwise() %>%
      mutate(
        loadtype = infer_loadtype(
          truth_column = Truth,
          loadtypes    = aldexrepo[[Study]]$studydemographics$loadtype
        )
      ) %>%
      ungroup()
  }
  list(export_table = combined_df)
}
