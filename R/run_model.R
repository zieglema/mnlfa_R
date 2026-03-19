# =============================================================================
# R/run_model.R
# Model estimation and summary extraction
# =============================================================================

#' Fit a moderated SEM via mxsem / OpenMx
#'
#' Converts a lavaan-style model string to an OpenMx model object using
#' [mxsem::mxsem()], fits it with [OpenMx::mxTryHard()], and returns a
#' richly annotated list that is consumed by the downstream significance- and
#' parameter-extraction functions.
#'
#' @param model_string A character string containing the lavaan-style model
#'   syntax (single string with `\n` separators, or a character vector of
#'   lines).
#' @param data A `data.frame` (or coercible object) containing the observed
#'   variables and moderators.
#' @param extract_fun Function used to extract fit statistics and parameters
#'   from the OpenMx summary.  Defaults to [extract_mx_info()].
#' @param scale_loadings Logical. If `TRUE` (default), the first factor
#'   loading is fixed to 1 for identification.  Set to `FALSE` when
#'   the loading constraint is encoded in `model_string` (e.g. `1*item`).
#' @param scale_latent_variances Logical. If `TRUE`, the latent variance is
#'   fixed to 1.  Use `TRUE` for screening-step models, `FALSE` for full
#'   domain models.  Default is `FALSE`.
#' @param start_values Optional named numeric vector of starting values passed
#'   to [mxsem::set_starting_values()].
#' @param add_ref_models Logical. Whether to compute reference models needed
#'   for CFI/TLI/RMSEA.  Default is `TRUE`.  **For all MNLFA pipeline steps,
#'   set this to `FALSE`**: MNLFA models use definition variables (person-
#'   specific parameters), which makes the standard independence reference model
#'   inapplicable.  The SABIC (returned in `fit_df$SABIC`) is used for all
#'   model comparisons instead.  Setting `add_ref_models = FALSE` also
#'   substantially speeds up model fitting.
#' @param moderators Optional character vector of moderator variable names.
#'   When supplied, the names are stored as an attribute on `param_df` so that
#'   downstream significance functions can recover the original moderator
#'   labels from numeric indices in parameter names.  The order must match the
#'   numeric suffixes used in parameter labels (e.g. `m_E_1_01` corresponds to
#'   `moderators[1]`).
#' @param seed Optional integer seed for reproducibility.
#' @param extraTries Number of additional optimisation attempts passed to
#'   [OpenMx::mxTryHard()].  Default is `10`.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`model_string`}{The model syntax as supplied.}
#'   \item{`script`}{The OpenMx model object (before fitting).}
#'   \item{`fit`}{The fitted OpenMx model object.}
#'   \item{`ref_models`}{Reference model list (or `NULL` if
#'     `add_ref_models = FALSE`).}
#'   \item{`summary`}{The OpenMx model summary.}
#'   \item{`extracted`}{Raw output of `extract_fun`.}
#'   \item{`fit_df`}{Data frame of fit statistics.}
#'   \item{`param_df`}{Data frame of parameter estimates, standard errors,
#'     z-values, and p-values, enriched with columns `type`, `facet`, and
#'     `domain`.}
#' }
#'
#' @seealso [extract_mx_info()], [significant_mean_moderators()],
#'   [significant_item_moderators()]
#' @export
run_mxsem <- function(model_string,
                      data,
                      extract_fun            = extract_mx_info,
                      scale_loadings         = TRUE,
                      scale_latent_variances = FALSE,
                      start_values           = NULL,
                      add_ref_models         = TRUE,
                      moderators             = NULL,
                      seed                   = NULL,
                      extraTries             = 10) {

  if (!is.null(seed)) set.seed(seed)

  script <- mxsem::mxsem(model_string,
                          data                   = as.data.frame(data),
                          scale_loadings         = scale_loadings,
                          scale_latent_variances = scale_latent_variances)

  if (!is.null(start_values)) {
    script <- mxsem::set_starting_values(script, values = start_values)
  }

  fit <- OpenMx::mxTryHard(script, extraTries = extraTries, greenOK = TRUE)

  if (add_ref_models) {
    ref_models  <- OpenMx::mxRefModels(fit, run = TRUE)
    summary_fit <- summary(fit, refModels = ref_models)
  } else {
    ref_models  <- NULL
    summary_fit <- summary(fit)
  }

  extracted <- extract_fun(summary_fit)

  if (!is.null(moderators)) {
    attr(fit,       "moderators") <- moderators
    attr(extracted, "moderators") <- moderators
  }
  attr(fit,       "raw_data") <- data
  attr(extracted, "raw_data") <- data

  if (!is.null(extracted$parameters)) {
    param_df <- extracted$parameters

    # Assign parameter types
    param_df$type <- NA_character_
    param_df$type[grepl("^m_",   param_df$name)] <- "latent_mean"
    param_df$type[grepl("^v_",   param_df$name)] <- "variance"
    param_df$type[grepl("^l_",   param_df$name)] <- "loading"
    param_df$type[grepl("^int_", param_df$name)] <- "intercept"
    param_df$type[is.na(param_df$type)]           <- "other"

    # Assign facet labels
    param_df$facet <- NA_character_
    m_rows    <- grepl("^m_",   param_df$name)
    int_rows  <- grepl("^int_", param_df$name)
    var_rows  <- which(param_df$type == "variance")
    load_rows <- which(param_df$type == "loading")

    param_df$facet[m_rows]    <- sub("^m_([A-Z]_[0-9]+).*",    "\\1", param_df$name[m_rows])
    param_df$facet[int_rows]  <- sub("^int_([A-Z]_[0-9]+)_.*", "\\1", param_df$name[int_rows])
    param_df$facet[var_rows]  <- sub("^v_([A-Z]_[0-9]+)_.*",   "\\1", param_df$name[var_rows])
    param_df$facet[load_rows] <- sub("^l_([A-Z]_[0-9]+)_.*",   "\\1", param_df$name[load_rows])
    param_df$facet[grepl("^m_acq_", param_df$name)] <- "acq"
    param_df$facet[param_df$row == "acq"]            <- "acq"

    # Parse model syntax to build a factor -> items map (uses stringr::)
    model_lines <- stringr::str_split(model_string, "\n")[[1]]
    factor_map  <- list()
    for (line in model_lines) {
      line <- stringr::str_trim(line)
      if (stringr::str_detect(line, "=~")) {
        parts       <- stringr::str_split(line, "=~")[[1]]
        factor_name <- stringr::str_trim(parts[1])
        items_part  <- stringr::str_trim(parts[2])
        items       <- stringr::str_split(items_part, "\\+")[[1]]
        items       <- stringr::str_trim(items)
        items       <- gsub(".*\\*", "", items)
        factor_map[[factor_name]] <- items
      }
    }

    # Variance parameters: annotate row/col
    for (i in var_rows) {
      param_df$row[i] <- param_df$facet[i]
      index <- as.numeric(sub(".*_([0-9]+)$", "\\1", param_df$name[i]))
      if (!is.na(index) && !is.null(moderators) && index <= length(moderators)) {
        param_df$col[i] <- moderators[index]
      }
    }

    # Latent mean intercepts: row/col
    base_mean_rows <- param_df$type == "latent_mean" & param_df$matrix == "M"
    param_df$row[base_mean_rows] <- param_df$col[base_mean_rows]
    param_df$col[base_mean_rows] <- NA

    # Manifest intercept base parameters: row/col
    base_int_rows <- param_df$type == "intercept" & param_df$matrix == "M"
    param_df$row[base_int_rows] <- param_df$col[base_int_rows]
    param_df$col[base_int_rows] <- NA

    # Covariance moderation parameters
    covar_rows <- grepl("^r\\d+_cov_[A-Z]_\\d+_[A-Z]_\\d+$", param_df$name)
    param_df$type[covar_rows] <- "covariance"
    for (i in which(covar_rows)) {
      index <- as.numeric(sub("^r(\\d+)_.*", "\\1", param_df$name[i]))
      if (!is.na(index) && !is.null(moderators) && index >= 1 && index <= length(moderators)) {
        param_df$col[i] <- moderators[index]
      }
    }

    # Loading parameters: annotate row/col
    if (!is.null(moderators)) {
      base_item_map <- list()
      for (i in load_rows) {
        label <- param_df$name[i]
        parts <- unlist(strsplit(label, "_"))
        if (length(parts) == 5) {
          facet        <- paste0(parts[2], "_", parts[3])
          base_or_idx  <- parts[4]
          index_last   <- parts[5]
          facet_items  <- factor_map[[facet]]
          if (!is.null(facet_items)) {
            if (base_or_idx == "0") {
              item_index <- as.numeric(index_last)
              if (!is.na(item_index) && item_index <= length(facet_items)) {
                item_name <- facet_items[item_index]
                param_df$row[i] <- item_name
                param_df$col[i] <- NA
                base_item_map[[paste0(facet, "_", item_index)]] <- item_name
              }
            } else {
              item_index      <- as.numeric(base_or_idx)
              moderator_index <- as.numeric(index_last)
              item_name       <- base_item_map[[paste0(facet, "_", item_index)]]
              param_df$row[i] <- if (!is.null(item_name)) item_name else NA
              if (!is.na(moderator_index) && moderator_index <= length(moderators)) {
                param_df$col[i] <- moderators[moderator_index]
              }
            }
          }
        }
      }
    }

    param_df$domain      <- sub("_.*", "", param_df$facet)
    extracted$parameters <- param_df

    if (!is.null(moderators)) {
      attr(extracted$parameters, "moderators") <- moderators
    }
  }

  return(list(
    model_string = model_string,
    script       = script,
    fit          = fit,
    ref_models   = ref_models,
    summary      = summary_fit,
    extracted    = extracted,
    fit_df       = extracted$fit,
    param_df     = extracted$parameters
  ))
}


#' Extract fit statistics and parameters from an OpenMx model summary
#'
#' Computes SABIC from the model log-likelihood and number of parameters,
#' and organises the parameter data frame with z-values and two-sided p-values.
#'
#' @param model An object returned by [base::summary()] on a fitted OpenMx
#'   model (optionally including reference models for CFI/TLI).
#'
#' @return A named list with two components:
#' \describe{
#'   \item{`fit`}{Single-row data frame of fit statistics: `minus2LL`, `AIC`,
#'     `BIC`, `SABIC`, `CFI`, `TLI`, `RMSEA`, `RMSEA_lower`, `RMSEA_upper`,
#'     `Chi2`, `df`, `p_value`, `n`.}
#'   \item{`parameters`}{Data frame with columns `model`, `name`, `matrix`,
#'     `row`, `col`, `Estimate`, `Std.Error`, `z`, `p`.}
#' }
#'
#' @seealso [run_mxsem()]
#' @export
extract_mx_info <- function(model) {
  model_name <- model$modelName
  LL    <- -0.5 * model$Minus2LogLikelihood
  k     <- model$estimatedParameters
  n     <- model$numObs
  SABIC <- -2 * LL + k * log((n + 2) / 24)

  fit_stats <- data.frame(
    minus2LL    = model$Minus2LogLikelihood,
    AIC         = model$AIC.Mx,
    BIC         = model$BIC.Mx,
    SABIC       = SABIC,
    CFI         = model$CFI,
    TLI         = model$TLI,
    RMSEA       = model$RMSEA,
    RMSEA_lower = model$RMSEACI["lower"],
    RMSEA_upper = model$RMSEACI["upper"],
    Chi2        = model$Chi,
    df          = model$ChiDoF,
    p_value     = round(model$p, 5),
    n           = n,
    stringsAsFactors = FALSE
  )

  if (!is.null(model$parameters) && nrow(model$parameters) > 0) {
    param_df       <- model$parameters
    param_df$model <- model_name
    if (!"row" %in% names(param_df)) param_df$row <- NA
    if (!"col" %in% names(param_df)) param_df$col <- NA
    param_df <- param_df[, c("model", "name", "matrix", "row", "col",
                             "Estimate", "Std.Error")]
    param_df$z <- param_df$Estimate / param_df$Std.Error
    param_df$p <- 2 * (1 - stats::pnorm(abs(param_df$z)))
  } else {
    param_df <- data.frame()
  }

  out        <- list(fit = fit_stats, parameters = param_df)
  moderators <- attr(model, "moderators")
  if (!is.null(moderators)) attr(out, "moderators") <- moderators
  return(out)
}
