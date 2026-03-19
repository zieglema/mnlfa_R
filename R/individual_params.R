# =============================================================================
# R/individual_params.R
# Functions for computing person-specific parameter values from a fitted
# moderated SEM.  Each function takes the parameter data frame (`param_df`)
# from a [run_mxsem()] result and the original data matrix, and returns a
# data frame with one row per person.
# =============================================================================

#' Compute person-specific item intercepts
#'
#' For each item with an intercept moderation in `param_df`, calculates the
#' individual-level intercept as:
#' `intercept_base + sum(beta_k * moderator_k)`.
#'
#' @param param_df Data frame returned in the `$param_df` slot of a
#'   [run_mxsem()] result.
#' @param data Data frame with the original observations and moderator columns.
#' @param significance_level Optional numeric. If supplied, only moderators
#'   with `p <= significance_level` are included.
#' @param p_adjust_method Optional character. P-value adjustment method passed
#'   to [stats::p.adjust()] before applying `significance_level`.
#' @param verbose Logical. If `TRUE`, additional diagnostic columns are added.
#'   Default is `FALSE`.
#'
#' @return A data frame with `nrow(data)` rows and one prediction column per
#'   item (named `<item>_pred`) plus intercept columns (`<item>_intercept`).
#'
#' @seealso [build_all_individual_parameters()]
#' @export
build_individual_intercepts <- function(param_df, data,
                                        significance_level = NULL,
                                        p_adjust_method    = NULL,
                                        verbose            = FALSE) {
  int_params <- param_df[param_df$type == "intercept", ]
  items      <- unique(int_params$row)
  results    <- list()

  for (item in items) {
    item_params   <- int_params[int_params$row == item, ]
    intercept_row <- item_params[grepl(paste0("^int_.*_", item, "$"),
                                       item_params$name), ]
    intercept     <- if (nrow(intercept_row) != 1) 0 else intercept_row$Estimate

    mod_params <- item_params[!item_params$name %in% intercept_row$name, ]
    if (!is.null(significance_level)) {
      if (!is.null(p_adjust_method)) {
        adj_p      <- stats::p.adjust(mod_params$p, method = p_adjust_method)
        mod_params <- mod_params[adj_p <= significance_level | is.na(adj_p), ]
      } else {
        mod_params <- mod_params[mod_params$p <= significance_level |
                                   is.na(mod_params$p), ]
      }
    }

    all_mods <- unique(mod_params$col)
    item_df  <- data.frame(row_id = seq_len(nrow(data)))
    item_df[[paste0(item, "_intercept")]] <- intercept
    item_df[[paste0(item, "_pred")]]      <- intercept

    for (mod in all_mods) {
      item_df[[mod]] <- data[[mod]]
      if (verbose) {
        item_df[[paste0(item, "_", mod, "_used")]]   <- mod %in% mod_params$col
        item_df[[paste0(item, "_", mod, "_effect")]] <-
          if (mod %in% mod_params$col)
            mod_params$Estimate[mod_params$col == mod] * data[[mod]]
          else 0
      }
    }

    for (j in seq_len(nrow(mod_params))) {
      mod <- mod_params$col[j]
      est <- mod_params$Estimate[j]
      item_df[[paste0(item, "_pred")]] <-
        item_df[[paste0(item, "_pred")]] + est * data[[mod]]
    }
    results[[item]] <- item_df
  }

  final_df <- Reduce(function(x, y) {
    common_cols <- setdiff(intersect(names(x), names(y)), "row_id")
    merge(x, y[, !names(y) %in% common_cols, drop = FALSE], by = "row_id", all = TRUE)
  }, results)
  final_df$row_id <- NULL
  return(final_df)
}


#' Compute person-specific factor loadings
#'
#' For each item with a loading moderation in `param_df`, calculates the
#' individual-level loading as:
#' `baseline_loading + sum(beta_k * moderator_k)`.
#'
#' @inheritParams build_individual_intercepts
#'
#' @return A data frame with one prediction column per item
#'   (`<item>_pred`).
#'
#' @seealso [build_all_individual_parameters()]
#' @export
build_individual_loadings <- function(param_df, data,
                                      significance_level = NULL,
                                      p_adjust_method    = NULL,
                                      verbose            = FALSE) {
  load_params <- param_df[param_df$type == "loading", ]
  facets      <- unique(load_params$facet)
  results     <- list()

  for (facet in facets) {
    facet_params <- load_params[load_params$facet == facet, ]
    items        <- unique(facet_params$row)
    item_results <- list()

    for (item in items) {
      item_params <- facet_params[facet_params$row == item, ]
      if (nrow(item_params) == 0) next
      base_row    <- item_params[grepl("_0_", item_params$name), ]
      if (nrow(base_row) != 1) {
        warning(paste("No/multiple base loading for", item)); next
      }

      intercept  <- base_row$Estimate
      mod_params <- item_params[!grepl("_0_", item_params$name), ]

      if (!is.null(significance_level)) {
        if (!is.null(p_adjust_method)) {
          adj_p      <- stats::p.adjust(mod_params$p, method = p_adjust_method)
          mod_params <- mod_params[adj_p <= significance_level | is.na(adj_p), ]
        } else {
          mod_params <- mod_params[mod_params$p <= significance_level |
                                     is.na(mod_params$p), ]
        }
      }

      item_df <- data.frame(row_id = seq_len(nrow(data)))
      item_df[[paste0(item, "_intercept")]] <- intercept
      item_df[[paste0(item, "_pred")]]      <- intercept

      for (j in seq_len(nrow(mod_params))) {
        mod    <- mod_params$col[j]
        est    <- mod_params$Estimate[j]
        if (mod %in% names(data)) {
          effect <- est * data[[mod]]
          item_df[[paste0(item, "_", mod, "_effect")]] <- effect
          item_df[[paste0(item, "_pred")]] <- item_df[[paste0(item, "_pred")]] + effect
          item_df[[paste0(item, "_", mod, "_used")]] <- TRUE
          if (verbose) message("Used ", mod, " for loading of ", item)
        } else {
          item_df[[paste0(item, "_", mod, "_effect")]] <- 0
          item_df[[paste0(item, "_", mod, "_used")]]   <- FALSE
          warning(paste("Moderator", mod, "not in data for item", item))
        }
      }
      item_results[[item]] <- item_df
    }

    if (length(item_results) > 0) {
      facet_df <- do.call(cbind, lapply(item_results,
                                        function(df) df[, -1, drop = FALSE]))
      facet_df <- cbind(row_id = seq_len(nrow(facet_df)), facet_df)
      facet_df$row_id <- NULL
      results[[facet]] <- facet_df
    }
  }

  if (length(results) == 0) return(NULL)
  final_df <- do.call(cbind, results)
  names(final_df) <- gsub(".*\\.", "", names(final_df))
  return(final_df)
}


#' Compute person-specific latent variances
#'
#' For each facet with a variance moderation in `param_df`, calculates the
#' individual-level latent variance as:
#' `exp(sum(gamma_k * moderator_k))`.
#'
#' @inheritParams build_individual_intercepts
#'
#' @return A data frame with prediction columns `<facet>_pred` and, for
#'   `verbose = TRUE`, additional detail columns.
#'
#' @seealso [build_all_individual_parameters()]
#' @export
build_individual_latent_variances <- function(param_df, data,
                                              significance_level = NULL,
                                              p_adjust_method    = NULL,
                                              verbose            = FALSE) {
  var_params      <- param_df[param_df$type == "variance", ]
  facets          <- unique(var_params$facet)
  results         <- list()
  used_moderators <- c()

  for (facet in facets) {
    facet_params <- var_params[var_params$facet == facet, ]

    if (!is.null(significance_level)) {
      if (!is.null(p_adjust_method)) {
        adj_p        <- stats::p.adjust(facet_params$p, method = p_adjust_method)
        facet_params <- facet_params[adj_p <= significance_level | is.na(adj_p), ]
      } else {
        facet_params <- facet_params[facet_params$p <= significance_level |
                                       is.na(facet_params$p), ]
      }
    }
    if (nrow(facet_params) == 0) next

    valid_params <- facet_params[facet_params$col %in% names(data), ]
    if (nrow(valid_params) == 0) next

    if (verbose && nrow(valid_params) < nrow(facet_params)) {
      warning(sprintf("Moderator(s) %s not in data.",
                      paste(setdiff(facet_params$col, names(data)), collapse = ", ")))
    }

    mods            <- valid_params$col
    effects         <- valid_params$Estimate
    used_moderators <- union(used_moderators, mods)

    # Look up the free log-variance intercept v0_<facet>
    v0_name  <- paste0("v0_", facet)
    v0_row   <- param_df[param_df$name == v0_name & !is.na(param_df$name), ]
    v0_value <- if (nrow(v0_row) == 1) v0_row$Estimate else 0

    mod_matrix <- sapply(mods, function(mod) data[[mod]])
    if (is.null(dim(mod_matrix))) mod_matrix <- matrix(mod_matrix, ncol = 1)

    pred_log_var <- v0_value + mod_matrix %*% effects
    pred_var     <- exp(pred_log_var)

    df <- data.frame(row_id = seq_len(nrow(data)))
    df[[paste0(facet, "_intercept")]] <- 0
    df[[paste0(facet, "_pred")]]      <- as.vector(pred_var)

    if (verbose) {
      for (i in seq_along(mods)) {
        df[[mods[i]]]                                <- data[[mods[i]]]
        df[[paste0(facet, "_", mods[i], "_used")]]   <- TRUE
        df[[paste0(facet, "_", mods[i], "_effect")]] <- data[[mods[i]]] * effects[i]
      }
    }
    results[[facet]] <- df
  }

  if (length(results) == 0) return(data.frame())
  final_df <- Reduce(function(x, y) merge(x, y, by = "row_id", all = TRUE), results)
  final_df$row_id <- NULL
  for (mod in used_moderators) {
    if (!mod %in% names(final_df)) final_df[[mod]] <- data[[mod]]
  }
  return(final_df)
}


#' Compute person-specific latent covariance for one covariance parameter
#'
#' Uses the Fisher-z formula to compute individual-level covariances from a
#' `rho` parameter that is a linear function of moderators.
#'
#' @param param_df Data frame from `$param_df` of a [run_mxsem()] result.
#' @param data Data frame with the original observations and moderators.
#' @param cov_name Character. Name of the covariance (e.g. `"cov_E_1_E_2"`).
#' @inheritParams build_individual_intercepts
#'
#' @return A data frame with columns `<cov_name>_rho` and
#'   `<cov_name>_pred`.
#'
#' @seealso [build_all_individual_parameters()]
#' @export
build_individual_latent_covariance <- function(param_df, data,
                                               cov_name,
                                               significance_level = NULL,
                                               p_adjust_method    = NULL,
                                               verbose            = FALSE) {
  rho_rows   <- grepl(paste0("^r\\d+_", cov_name, "$"), param_df$name)
  rho_params <- param_df[rho_rows, ]
  if (nrow(rho_params) == 0) stop(paste("No parameters found for", cov_name))

  rho_params$mod_index <- as.integer(gsub("^r(\\d+)_.*", "\\1", rho_params$name))

  if (!is.null(attr(param_df, "moderators"))) {
    moderators           <- attr(param_df, "moderators")
    rho_params$moderator <- moderators[rho_params$mod_index]
  } else if ("col" %in% names(param_df)) {
    rho_params <- merge(rho_params, param_df[, c("name", "col")],
                        by = "name", all.x = TRUE)
    names(rho_params)[names(rho_params) == "col.y"] <- "moderator"
  } else {
    rho_params$moderator <- NA
    warning("Moderator names could not be matched.")
  }

  rho_params$param_name <- rho_params$name
  intercept  <- rho_params$Estimate[rho_params$mod_index == 0]
  mod_params <- rho_params[rho_params$mod_index > 0, ]

  if (!is.null(significance_level)) {
    if (!is.null(p_adjust_method)) {
      mod_params$p_corrected <- stats::p.adjust(mod_params$p,
                                                method = p_adjust_method)
      mod_params <- mod_params[mod_params$p_corrected <= significance_level |
                                 is.na(mod_params$p_corrected), ]
    } else {
      mod_params <- mod_params[mod_params$p <= significance_level |
                                 is.na(mod_params$p), ]
    }
  }

  rho_vec   <- rep(intercept, nrow(data))
  used_mods <- character()

  if (nrow(mod_params) > 0) {
    for (i in seq_len(nrow(mod_params))) {
      mod_name <- mod_params$moderator[i]
      est      <- mod_params$Estimate[i]
      if (is.na(mod_name) || !mod_name %in% names(data)) {
        if (verbose) warning(paste("Moderator", mod_name, "not in data."))
        next
      }
      used_mods <- c(used_mods, mod_name)
      rho_vec   <- rho_vec + est * data[[mod_name]]
    }
  }

  vars      <- param_df[grepl("^var_", param_df$name), ]
  latents   <- regmatches(cov_name, gregexpr("V_\\d+", cov_name))[[1]]
  var_names <- paste0("var_", latents)
  var_values <- sapply(var_names, function(vn) {
    vrow <- vars[vars$name == vn, ]
    if (nrow(vrow) == 0) stop(paste("Missing variance parameter:", vn))
    vrow$Estimate[1]
  })

  cov_pred <- sqrt(prod(var_values)) * (exp(2 * rho_vec) - 1) / (exp(2 * rho_vec) + 1)

  out <- data.frame(rho = rho_vec, pred = cov_pred)
  names(out) <- c(paste0(cov_name, "_rho"), paste0(cov_name, "_pred"))

  if (verbose && length(used_mods) > 0) {
    for (mod in used_mods) {
      out[[mod]] <- data[[mod]]
      out[[paste0(cov_name, "_", mod, "_effect")]] <-
        data[[mod]] * mod_params$Estimate[mod_params$moderator == mod]
    }
  }
  return(out)
}


#' Compute person-specific latent means
#'
#' For each facet with a latent mean moderation in `param_df`, calculates the
#' individual-level latent mean as:
#' `mean_intercept + sum(alpha_k * moderator_k)`.
#'
#' @inheritParams build_individual_intercepts
#'
#' @return A data frame with `<facet>_pred` and `<facet>_intercept` columns
#'   for each facet with latent mean parameters.
#'
#' @seealso [build_all_individual_parameters()]
#' @export
build_individual_latent_means <- function(param_df, data,
                                          verbose            = FALSE,
                                          significance_level = NULL,
                                          p_adjust_method    = NULL) {
  lm_params       <- param_df[param_df$type == "latent_mean", ]
  facets          <- unique(lm_params$facet)
  results         <- list()
  used_moderators <- c()

  for (facet in facets) {
    facet_params  <- lm_params[lm_params$facet == facet, ]

    intercept_row <- facet_params[facet_params$matrix == "M" &
                                    facet_params$row == facet &
                                    is.na(facet_params$col), ]
    if (nrow(intercept_row) == 0) {
      if (verbose) warning(paste("No intercept for", facet, "- assuming 0"))
      intercept <- 0
    } else if (nrow(intercept_row) == 1) {
      intercept <- intercept_row$Estimate
    } else {
      warning(paste("Multiple intercepts for", facet, "- skipping"))
      next
    }

    mod_params <- facet_params[facet_params$name != intercept_row$name, ]
    if (!is.null(significance_level)) {
      if (!is.null(p_adjust_method)) {
        adj_p      <- stats::p.adjust(mod_params$p, method = p_adjust_method)
        mod_params <- mod_params[adj_p <= significance_level | is.na(adj_p), ]
      } else {
        mod_params <- mod_params[mod_params$p <= significance_level |
                                   is.na(mod_params$p), ]
      }
    }

    all_mods        <- unique(mod_params$col)
    used_moderators <- union(used_moderators, all_mods)

    facet_df <- data.frame(row_id = seq_len(nrow(data)))
    facet_df[[paste0(facet, "_pred")]]      <- intercept
    facet_df[[paste0(facet, "_intercept")]] <- intercept

    if (nrow(mod_params) == 0) {
      if (verbose) message("No significant moderators for ", facet, " - intercept only")
      results[[facet]] <- facet_df
      next
    }

    for (mod in all_mods) {
      if (verbose) {
        facet_df[[paste0(facet, "_", mod, "_used")]] <- mod %in% mod_params$col
        facet_df[[paste0(facet, "_", mod, "_effect")]] <-
          if (mod %in% mod_params$col)
            mod_params$Estimate[mod_params$col == mod] * data[[mod]]
          else 0
      }
    }
    for (j in seq_len(nrow(mod_params))) {
      mod  <- mod_params$col[j]
      est  <- mod_params$Estimate[j]
      facet_df[[paste0(facet, "_pred")]] <-
        facet_df[[paste0(facet, "_pred")]] + est * data[[mod]]
    }
    results[[facet]] <- facet_df
  }

  if (length(results) == 0) return(NULL)
  final_df <- Reduce(function(x, y) merge(x, y, by = "row_id", all = TRUE), results)
  final_df$row_id <- NULL
  for (mod in used_moderators) {
    if (!mod %in% names(final_df)) final_df[[mod]] <- data[[mod]]
  }
  return(final_df)
}


#' Compute all person-specific parameters at once
#'
#' Convenience wrapper that calls all four individual-parameter builders
#' ([build_individual_latent_means()], [build_individual_latent_variances()],
#' [build_individual_loadings()], [build_individual_intercepts()]) and returns
#' their results in a named list.
#'
#' @param param_df Data frame from the `$param_df` slot of a [run_mxsem()]
#'   result.
#' @param data Data frame with the original observations and moderators.
#' @param significance_level Optional numeric significance threshold.
#' @param p_adjust_method Optional p-value adjustment method.
#' @param verbose Logical.  Default is `FALSE`.
#'
#' @return A named list with four data frames: `latent_means`,
#'   `latent_variances`, `loadings`, `intercepts`.
#'
#' @seealso [build_individual_latent_means()],
#'   [build_individual_latent_variances()], [build_individual_loadings()],
#'   [build_individual_intercepts()]
#' @export
build_all_individual_parameters <- function(param_df, data,
                                            significance_level = NULL,
                                            p_adjust_method    = NULL,
                                            verbose            = FALSE) {
  message("build_all_individual_parameters: ", nrow(param_df), " parameter rows")
  list(
    latent_means     = build_individual_latent_means(
      param_df, data, verbose = verbose,
      significance_level = significance_level,
      p_adjust_method    = p_adjust_method),
    latent_variances = build_individual_latent_variances(
      param_df, data,
      significance_level = significance_level,
      p_adjust_method    = p_adjust_method,
      verbose = verbose),
    loadings         = build_individual_loadings(
      param_df, data,
      significance_level = significance_level,
      p_adjust_method    = p_adjust_method,
      verbose = verbose),
    intercepts       = build_individual_intercepts(
      param_df, data,
      significance_level = significance_level,
      p_adjust_method    = p_adjust_method,
      verbose = verbose)
  )
}
