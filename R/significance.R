# =============================================================================
# R/significance.R
# Functions for identifying significant moderators after model fitting.
# All functions return a list with at minimum a `$table` data frame and a
# `$summary` character string.
# =============================================================================

#' Identify significant moderators of the latent mean
#'
#' Extracts regression coefficients for latent mean regressions
#' (`LV ~ m_LV_01*mod1 + ...`) from the parameter data frame, applies an
#' optional p-value correction, and returns a table of results together with
#' a plain-text summary.
#'
#' @param mxsem_result A list returned by [run_mxsem()].
#' @param alpha Significance level.  Default is `0.05`.
#' @param correction_method Method passed to [stats::p.adjust()].  Common
#'   choices: `"none"`, `"BH"`, `"BY"`, `"bonferroni"`.  Default is `"none"`.
#'
#' @return A named list:
#' \describe{
#'   \item{`table`}{Data frame with one row per moderator–latent combination.}
#'   \item{`summary`}{Character string describing the results.}
#' }
#'
#' @seealso [run_mxsem()], [combine_sig_moderators()]
#' @export
significant_mean_moderators <- function(mxsem_result, alpha = 0.05,
                                        correction_method = "none") {
  df         <- mxsem_result$param_df
  df$matrix  <- as.character(df$matrix)

  mean_rows <- df[df$type == "latent_mean" & df$matrix == "A", , drop = FALSE]
  if (nrow(mean_rows) == 0) {
    return(list(table = data.frame(),
                summary = "No latent mean regressions found."))
  }

  latents       <- unique(mean_rows$row)
  combined      <- list()
  summary_lines <- c()

  for (latent in latents) {
    regressors <- mean_rows[mean_rows$row == latent &
                              mean_rows$col != "one", , drop = FALSE]
    if (nrow(regressors) == 0) {
      summary_lines <- c(summary_lines,
        paste0("No moderators for latent mean of ", latent, "."))
      next
    }
    regressors$significant           <- regressors$p <= alpha
    regressors$p_corrected           <- stats::p.adjust(regressors$p,
                                                        method = correction_method)
    regressors$significant_corrected <- regressors$p_corrected <= alpha
    regressors$param_name            <- regressors$name
    regressors$moderator             <- regressors$col
    regressors$latent                <- latent

    combined[[latent]] <- regressors[, c(
      "param_name", "latent", "moderator", "Estimate", "Std.Error", "z", "p",
      "p_corrected", "significant", "significant_corrected", "facet"
    )]

    sig_mods <- na.omit(regressors$moderator[regressors$significant_corrected])
    summary_lines <- c(summary_lines,
      if (length(sig_mods) == 0)
        paste0("No significant moderators for latent mean of ", latent,
               " at alpha = ", alpha, " after ", correction_method, " correction.")
      else
        paste0("Latent mean of ", latent, ": significant moderators (",
               correction_method, "): ", paste(unique(sig_mods), collapse = ", "), ".")
    )
  }

  if (length(combined) == 0)
    return(list(table   = data.frame(),
                summary = paste(summary_lines, collapse = "\n")))

  return(list(
    table   = do.call(rbind, combined),
    summary = paste(summary_lines, collapse = "\n")
  ))
}


#' Identify significant moderators of the latent variance
#'
#' Extracts `v_*` definition-variable parameters from the parameter data frame
#' and tests their significance.
#'
#' @inheritParams significant_mean_moderators
#'
#' @return A named list with `$table` and `$summary`.
#'
#' @seealso [add_latent_variance_moderation()], [combine_sig_moderators()]
#' @export
significant_var_moderators <- function(mxsem_result, alpha = 0.05,
                                       correction_method = "none") {
  df         <- mxsem_result$param_df
  moderators <- .get_moderators(mxsem_result)

  variance_mods <- df[df$matrix == "new_parameters" & grepl("^v_", df$name), ]
  if (nrow(variance_mods) == 0) {
    return(list(table   = data.frame(),
                summary = "No latent variance moderation terms found."))
  }

  mod_index <- suppressWarnings(
    as.integer(gsub(".*_(\\d+)$", "\\1", variance_mods$name)))
  variance_mods$moderator <- if (!is.null(moderators) &&
                                 length(moderators) >= max(mod_index, na.rm = TRUE))
    moderators[mod_index] else NA

  variance_mods$significant           <- variance_mods$p <= alpha
  variance_mods$p_corrected           <- stats::p.adjust(variance_mods$p,
                                                         method = correction_method)
  variance_mods$significant_corrected <- variance_mods$p_corrected <= alpha
  variance_mods$param_name            <- variance_mods$name

  mod_table <- variance_mods[, c("param_name", "moderator", "Estimate",
                                  "Std.Error", "z", "p", "p_corrected",
                                  "significant", "significant_corrected")]
  sig_mods     <- mod_table$moderator[mod_table$significant_corrected]
  summary_text <- if (length(sig_mods) == 0)
    sprintf("No significant latent variance moderators at alpha = %.3f after %s correction.",
            alpha, correction_method)
  else
    sprintf("Significant latent variance moderators at alpha = %.3f after %s correction: %s.",
            alpha, correction_method, paste(unique(sig_mods), collapse = ", "))

  return(list(table = mod_table, summary = summary_text))
}


#' Identify significant moderators of latent covariances
#'
#' Extracts `r{k}_cov_LV1_LV2` parameters from the parameter data frame and
#' tests their significance.
#'
#' @inheritParams significant_mean_moderators
#'
#' @return A named list with `$table` and `$summary`.
#'
#' @seealso [add_covariance_moderation()]
#' @export
significant_covar_moderators <- function(mxsem_result, alpha = 0.05,
                                         correction_method = "none") {
  if ("param_df" %in% names(mxsem_result[[1]])) {
    mxsem_result <- mxsem_result[[1]]
  }
  df <- mxsem_result$param_df
  if (is.null(df) || !"name" %in% names(df)) {
    return(list(table = data.frame(), summary = "No parameter data frame found."))
  }

  mod_rows        <- grepl("^r[1-9]\\d*_cov_[A-Z]_\\d+_[A-Z]_\\d+$", df$name)
  covariance_mods <- df[mod_rows, ]
  if (nrow(covariance_mods) == 0) {
    return(list(table   = data.frame(),
                summary = "No moderated covariance parameters found."))
  }

  valid_mods <- covariance_mods[!is.na(covariance_mods$Std.Error) &
                                  covariance_mods$Std.Error != 0, ]
  valid_mods$param_name            <- valid_mods$name
  valid_mods$p_corrected           <- stats::p.adjust(valid_mods$p,
                                                      method = correction_method)
  valid_mods$significant           <- valid_mods$p <= alpha
  valid_mods$significant_corrected <- valid_mods$p_corrected <= alpha

  covariance_mods$param_name <- covariance_mods$name
  covariance_mods <- merge(
    covariance_mods,
    valid_mods[, c("param_name", "p_corrected", "significant", "significant_corrected")],
    by = "param_name", all.x = TRUE
  )

  mod_index  <- as.integer(sub("^r(\\d+)_.*", "\\1", covariance_mods$param_name))
  moderators <- .get_moderators(mxsem_result)

  if (!is.null(moderators) && max(mod_index, na.rm = TRUE) <= length(moderators)) {
    covariance_mods$moderator <- moderators[mod_index]
  } else {
    covariance_mods$moderator <- covariance_mods$param_name
    warning("Could not assign moderator names. Using param_name instead.")
  }

  mod_table <- covariance_mods[, c("param_name", "type", "moderator", "Estimate",
                                    "Std.Error", "z", "p", "p_corrected",
                                    "significant", "significant_corrected",
                                    "facet", "domain")]
  sig_mods     <- na.omit(mod_table$moderator[mod_table$significant_corrected])
  summary_text <- if (length(sig_mods) == 0)
    paste0("No significant latent covariance moderators at alpha = ", alpha,
           " after ", correction_method, " correction.")
  else
    paste0("Significant latent covariance moderators at alpha = ", alpha,
           " after ", correction_method, " correction: ",
           paste(unique(sig_mods), collapse = ", "), ".")

  return(list(table = mod_table, summary = summary_text))
}


#' Identify significant item-level moderators (screening step)
#'
#' Extracts loading (`l_*`) and intercept (`int_*`) moderation parameters
#' from a single-item screening model result and tests their significance.
#'
#' @inheritParams significant_mean_moderators
#'
#' @return A named list with `$table` and `$summary`.
#'
#' @seealso [moderate_loadings_and_intercepts()], [run_mxsem()]
#' @export
significant_item_moderators <- function(mxsem_result, alpha = 0.05,
                                        correction_method = "none") {
  df         <- mxsem_result$param_df
  moderators <- .get_moderators(mxsem_result)

  loading_mods   <- subset(df, matrix == "new_parameters" &
                             grepl("^l_", name) & !grepl("_0_", name))
  intercept_mods <- subset(df, matrix == "A" & grepl("^int_", name))

  if (nrow(loading_mods)   > 0) loading_mods$type   <- "loading"
  if (nrow(intercept_mods) > 0) intercept_mods$type <- "intercept"

  item_mods <- rbind(loading_mods, intercept_mods)
  if (nrow(item_mods) == 0) {
    return(list(table   = data.frame(),
                summary = "No item-level moderation terms found."))
  }

  mod_index <- suppressWarnings(
    as.integer(gsub(".*_(\\d+)$", "\\1", item_mods$name)))
  item_mods$moderator <- if (!is.null(moderators) &&
                             length(moderators) >= max(mod_index, na.rm = TRUE))
    moderators[mod_index] else NA

  item_mods$significant           <- item_mods$p <= alpha
  item_mods$p_corrected           <- stats::p.adjust(item_mods$p,
                                                     method = correction_method)
  item_mods$significant_corrected <- item_mods$p_corrected <= alpha
  item_mods$param_name            <- item_mods$name

  mod_table <- item_mods[, c("param_name", "type", "moderator", "Estimate",
                             "Std.Error", "z", "p", "p_corrected",
                             "significant", "significant_corrected")]
  sig_mods     <- mod_table$moderator[mod_table$significant_corrected]
  summary_text <- if (length(sig_mods) == 0)
    sprintf("No significant item-level moderators at alpha = %.3f after %s correction.",
            alpha, correction_method)
  else
    sprintf("Significant item-level moderators at alpha = %.3f after %s correction: %s (%d).",
            alpha, correction_method,
            paste(unique(sig_mods), collapse = ", "), length(sig_mods))

  return(list(table = mod_table, summary = summary_text))
}


#' Extract fit indices from a list of item-level screening models
#'
#' Collates the `fit_df` component from each per-item model result into a
#' single data frame with `facet` and `item` columns added.
#'
#' @param facet_item_fits Named list (level 1: facets; level 2: items) of
#'   [run_mxsem()] results.
#'
#' @return A data frame with one row per facet–item combination.
#'
#' @export
extract_item_fit_indices <- function(facet_item_fits) {
  results <- list()
  for (facet_name in names(facet_item_fits)) {
    facet_models <- facet_item_fits[[facet_name]]
    for (model_name in names(facet_models)) {
      fit_df <- facet_models[[model_name]]$fit_df
      if (!is.null(fit_df)) {
        fit_row       <- fit_df[1, ]
        fit_row$facet <- facet_name
        fit_row$item  <- model_name
        results[[paste(facet_name, model_name, sep = ".")]] <- fit_row
      }
    }
  }
  final_df <- do.call(rbind, results)
  rownames(final_df) <- NULL
  return(final_df)
}


#' Combine screening-step significance results into a single data frame
#'
#' Merges the outputs of [significant_mean_moderators()],
#' [significant_var_moderators()], and per-item calls to
#' [significant_item_moderators()] into a standardised data frame that can
#' be passed to [build_full_domain_models()].
#'
#' @param sig.means Named list of latent mean significance results, one element
#'   per facet.  Each element is the direct output of
#'   [significant_mean_moderators()] (a list with `$table` and `$summary`).
#'   List names are the facet names (e.g. `list(P_1 = ..., P_2 = ...)`).
#' @param sig.variances Named list of latent variance significance results,
#'   same structure as `sig.means`, from [significant_var_moderators()].
#' @param item.mods Named list of item-level significance results.  Keys must
#'   use the composite `"facet_name.item_name"` format
#'   (e.g. `"P_1.p1_i5"`).  Each value is the output of
#'   [significant_item_moderators()].  In the pipeline this list is produced by
#'   `unlist(sig_items, recursive = FALSE)` after the per-facet screening loop.
#'
#' @return A data frame with columns: `facet`, `item`, `type`, `moderator`,
#'   `Estimate`, `Std.Error`, `z`, `p`, `p_corrected`, `significant`,
#'   `significant_corrected`, `param_name`, `source`.  The `type` column takes
#'   values `"mean"`, `"variance"`, `"loading"`, or `"intercept"`.
#'
#' @seealso [build_full_domain_models()], [significant_item_moderators()]
#' @export
combine_sig_moderators <- function(sig.means = NULL, sig.variances = NULL,
                                   item.mods = NULL) {
  extract_df <- function(list_obj, type, facet_from_name = TRUE) {
    if (is.null(list_obj)) return(NULL)
    do.call(rbind, lapply(names(list_obj), function(name) {
      obj <- list_obj[[name]]
      if (is.null(obj) || !("table" %in% names(obj))) return(NULL)
      tbl <- obj$table
      if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) return(NULL)
      df        <- tbl
      df$facet  <- if (facet_from_name) sub("\\..*", "", name) else name
      df$item   <- if (type == "item") sub(".*\\.", "", name) else NA
      if (type != "item") df$type <- type
      df$source <- name
      df
    }))
  }

  means_df     <- extract_df(sig.means,     type = "mean",     facet_from_name = FALSE)
  variances_df <- extract_df(sig.variances, type = "variance", facet_from_name = FALSE)
  items_df     <- extract_df(item.mods,     type = "item",     facet_from_name = TRUE)

  combined <- dplyr::bind_rows(means_df, variances_df, items_df)
  combined[, c("facet", "item", "type", "moderator",
               "Estimate", "Std.Error", "z", "p", "p_corrected",
               "significant", "significant_corrected", "param_name", "source")]
}


#' Identify significant item moderators from a fitted domain model
#'
#' Used in the re-screening step after fitting the full domain model.
#' Extracts loading and intercept moderation parameters and maps them back to
#' item names using the matrix label information stored in the fitted model.
#'
#' @inheritParams significant_mean_moderators
#'
#' @return A named list with `$table` (columns: `param_name`, `type`,
#'   `moderator`, `Estimate`, `Std.Error`, `z`, `p`, `p_corrected`,
#'   `significant`, `significant_corrected`, `facet`, `item`).
#'
#' @seealso [extract_all_significant_moderators()]
#' @export
significant_item_moderators_from_domain <- function(mxsem_result,
                                                     alpha = 0.05,
                                                     correction_method = "none") {
  if (!"param_df" %in% names(mxsem_result)) {
    if (is.list(mxsem_result) && length(mxsem_result) > 0) {
      mxsem_result <- mxsem_result[[1]]
    } else {
      stop("Input must be a fitted mxsem model or a list of such models.")
    }
  }

  df             <- mxsem_result$param_df
  moderators     <- .get_moderators(mxsem_result)
  manifest_items <- mxsem_result$script@manifestVars
  loading_labels <- mxsem_result$script@matrices$A@labels

  loading_mods   <- subset(df, matrix == "new_parameters" & grepl("^l_", name))
  loading_mods   <- loading_mods[!grepl("_0_\\d+$", loading_mods$name), ]
  intercept_mods <- subset(df, matrix == "A" & grepl("^int_", name))

  loading_mods$type   <- "loading"
  intercept_mods$type <- "intercept"

  loading_mods$facet   <- gsub("^l_([^_]+_[^_]+)_.*",   "\\1", loading_mods$name)
  intercept_mods$facet <- gsub("^int_([^_]+_[^_]+)_.*", "\\1", intercept_mods$name)

  mod_index <- suppressWarnings(
    as.integer(gsub(".*_(\\d+)$", "\\1",
                    c(loading_mods$name, intercept_mods$name)))
  )
  moderator_names <- if (!is.null(moderators) &&
                         length(moderators) >= max(mod_index, na.rm = TRUE))
    moderators[mod_index] else NA

  loading_mods$moderator   <- moderator_names[seq_len(nrow(loading_mods))]
  intercept_mods$moderator <- moderator_names[
    (nrow(loading_mods) + 1):length(moderator_names)]
  intercept_mods$item      <- intercept_mods$row

  get_item_for_loading <- function(param_name) {
    parts <- unlist(strsplit(param_name, "_"))
    if (length(parts) < 4) return(NA_character_)
    facet         <- paste(parts[2:3], collapse = "_")
    item_idx      <- sprintf("%02d", as.integer(parts[4]))
    label_pattern <- paste0(facet, "_load_", item_idx)
    match_row     <- rownames(loading_labels)[
      apply(loading_labels, 1,
            function(row) any(grepl(label_pattern, row, fixed = TRUE)))]
    if (length(match_row) > 0) match_row[1] else NA_character_
  }

  loading_mods$item <- vapply(loading_mods$name, get_item_for_loading, character(1))

  item_mods <- rbind(loading_mods, intercept_mods)
  item_mods$significant           <- item_mods$p <= alpha
  item_mods$p_corrected           <- stats::p.adjust(item_mods$p,
                                                     method = correction_method)
  item_mods$significant_corrected <- item_mods$p_corrected <= alpha
  item_mods$param_name            <- item_mods$name

  mod_table <- item_mods[, c("param_name", "type", "moderator", "Estimate",
                             "Std.Error", "z", "p", "p_corrected",
                             "significant", "significant_corrected",
                             "facet", "item")]
  return(list(table = mod_table))
}


#' Extract all significant moderators from a list of domain model results
#'
#' Iterates over a named list of [run_mxsem()] results (one per domain) and
#' calls the four significance functions (mean, variance, item, covariance)
#' with `tryCatch()` so that errors in one domain do not abort the whole
#' extraction.
#'
#' @param mxsem_results A named list of [run_mxsem()] results, one per domain.
#'   A single (unwrapped) result is also accepted.
#' @param alpha Significance level.  Default is `0.05`.
#' @param correction_method P-value correction method.  Default is `"none"`.
#'
#' @return A named list (one element per domain) of data frames, each with the
#'   columns required by [add_all_moderation_single_facet()].
#'
#' @seealso [significant_mean_moderators()], [significant_var_moderators()],
#'   [significant_item_moderators_from_domain()],
#'   [significant_covar_moderators()]
#' @export
extract_all_significant_moderators <- function(mxsem_results, alpha = 0.05,
                                               correction_method = "none") {
  if (!is.list(mxsem_results) ||
      all(c("param_df", "script") %in% names(mxsem_results))) {
    mxsem_results <- list(SINGLE = mxsem_results)
  }

  required_cols <- c("param_name", "type", "moderator", "Estimate", "Std.Error",
                     "z", "p", "p_corrected", "significant", "significant_corrected",
                     "facet", "item", "domain")
  ensure_columns <- function(df) {
    missing <- setdiff(required_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[, required_cols, drop = FALSE]
  }

  all_results <- list()

  for (domain_name in names(mxsem_results)) {
    message("Processing domain: ", domain_name)
    model_result <- mxsem_results[[domain_name]]

    if (is.list(model_result) && length(model_result) == 1 &&
        inherits(model_result[[1]], "list")) {
      model_result <- model_result[[1]]
    }

    combined <- list()

    tryCatch({
      message(" -> Checking latent mean moderators")
      m <- significant_mean_moderators(model_result, alpha, correction_method)$table
      if (nrow(m) > 0) {
        m$type   <- "mean"
        m$item   <- NA
        m$domain <- domain_name
        combined[["mean"]] <- ensure_columns(m)
      }
    }, error = function(e) message("   [mean] ", conditionMessage(e)))

    tryCatch({
      message(" -> Checking latent variance moderators")
      v <- significant_var_moderators(model_result, alpha, correction_method)$table
      if (nrow(v) > 0) {
        v$type   <- "variance"
        v$item   <- NA
        v$facet  <- sub("^v_([^_]+_[^_]+).*", "\\1", v$param_name)
        v$domain <- domain_name
        combined[["variance"]] <- ensure_columns(v)
      }
    }, error = function(e) message("   [variance] ", conditionMessage(e)))

    tryCatch({
      message(" -> Checking item-level moderators")
      i <- significant_item_moderators_from_domain(model_result, alpha,
                                                   correction_method)$table
      message("    -> Found ", nrow(i), " item rows")
      if (!is.null(i) && nrow(i) > 0) {
        i$domain <- domain_name
        combined[["item"]] <- ensure_columns(i)
      }
    }, error = function(e) message("   [item] ", conditionMessage(e)))

    tryCatch({
      message(" -> Checking latent covariance moderators")
      cv <- significant_covar_moderators(model_result, alpha, correction_method)$table
      if (nrow(cv) > 0) {
        cv$type   <- "covariance"
        cv$item   <- NA
        cv$domain <- domain_name
        combined[["covariance"]] <- ensure_columns(cv)
      }
    }, error = function(e) message("   [covariance] ", conditionMessage(e)))

    all_results[[domain_name]] <- if (length(combined) > 0)
      do.call(rbind, combined) else data.frame()
  }

  return(all_results)
}
