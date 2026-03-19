# =============================================================================
# R/model_building.R
# Functions that construct or modify lavaan-style model syntax strings.
# All functions accept model_string either as a single \n-delimited string or
# as a character vector of lines, and return the same format (character vector
# of lines unless stated otherwise).
# =============================================================================

#' Add latent mean intercept and moderator regressions to model syntax
#'
#' Appends a labelled intercept line (`latent ~ m_latent*1`) and one or more
#' moderator regression lines (`latent ~ m_latent_01*mod1 + ...`) to the
#' model string.  Any pre-existing mean/intercept lines for the latent variable
#' are removed first to avoid duplication.
#'
#' @param model_string Character. Lavaan-style model syntax (single string or
#'   character vector of lines).
#' @param moderators Character vector of moderator variable names.
#' @param identification Character. Scale identification strategy used when
#'   fitting the resulting model with [run_mxsem()].  `"variance"` (default,
#'   Bauer 2017) keeps the `LV ~~ 1*LV` constraint so the latent variance
#'   is fixed to 1; pass `scale_latent_variances = TRUE` to [run_mxsem()].
#'   `"loading"` removes the variance constraint so that scale identification
#'   is via the first loading being fixed to 1; pass
#'   `scale_loadings = TRUE, scale_latent_variances = FALSE` to [run_mxsem()].
#'
#' @return A single character string (lines joined by `\n`).
#'
#' @examples
#' m <- "F =~ x1 + x2 + x3\nF ~~ 1*F\nF ~ 0*1"
#' add_latent_mean_regression(m, c("age", "sex"))
#' add_latent_mean_regression(m, c("age", "sex"), identification = "loading")
#'
#' @export
add_latent_mean_regression <- function(model_string, moderators,
                                       identification = c("variance", "loading")) {
  identification <- match.arg(identification)
  model_lines <- unlist(strsplit(model_string, "\n"))

  lv_index <- grep("=~", model_lines)[1]
  lv_line  <- gsub("\\s+", "", model_lines[lv_index])
  latent   <- sub("=~.*", "", lv_line)
  rhs      <- sub(".*=~", "", lv_line)
  terms    <- unlist(strsplit(rhs, "\\+"))

  first_indicator       <- sub(".*\\*", "", terms[1])
  terms[1]              <- first_indicator
  new_loading_line      <- paste0(latent, " =~ ", paste(terms, collapse = " + "))
  model_lines[lv_index] <- new_loading_line

  model_lines <- model_lines[!grepl(paste0("^", latent, " ~ mean_"),              model_lines)]
  model_lines <- model_lines[!grepl(paste0("^", latent, "~~var_"),                model_lines)]
  model_lines <- model_lines[!grepl(paste0("^", latent, " ~ m_", latent, "\\*1$"), model_lines)]
  model_lines <- model_lines[!grepl(" ~ int_.*\\*1$",                             model_lines)]

  # Loading identification: remove variance = 1 constraint so that run_mxsem()
  # can identify the scale via the first loading instead.
  if (identification == "loading") {
    model_lines <- model_lines[
      !grepl(paste0("^", latent, "\\s*~~\\s*1\\*", latent), model_lines)]
  }

  intercept_label <- paste0("m_", latent)
  intercept_line  <- paste0(latent, " ~ ", intercept_label, "*1")

  preds           <- paste0("m_", latent, "_", sprintf("%02d", seq_along(moderators)),
                            "*", moderators)
  regression_line <- paste0(latent, " ~ ", paste(preds, collapse = " + "))

  updated_lines <- c(model_lines, intercept_line, regression_line)
  return(paste(updated_lines, collapse = "\n"))
}


#' Add latent variance moderation to model syntax (screening step)
#'
#' Replaces the raw variance statement for the latent variable with a
#' definition-variable formula `var_LV := exp(v0_LV + v_LV_01*data.mod1 + ...)`,
#' which implements a log-linear model for the conditional variance.
#'
#' @param model_string Character. Lavaan-style model syntax.
#' @param moderators Character vector of moderator names.
#' @param latent_name Optional character string. Name of the latent variable.
#'   If `NULL` (default), the first `=~` line is used.
#'
#' @return A single character string.
#'
#' @export
add_latent_variance_moderation <- function(model_string, moderators,
                                           latent_name = NULL) {
  lines <- unlist(strsplit(model_string, "\n"))
  if (is.null(latent_name)) {
    lv_line     <- lines[grepl("=~", lines)][1]
    latent_name <- sub("=~.*", "", gsub("\\s+", "", lv_line))
  }

  lines <- lines[!grepl(paste0("^", latent_name, "~~"), lines)]

  var_label       <- paste0("var_", latent_name)
  intercept_param <- paste0("v0_", latent_name)
  mod_params      <- paste0("v_", latent_name, "_", sprintf("%02d", seq_along(moderators)))

  declarations  <- c(paste0("!", intercept_param, ";"), paste0("!", mod_params, ";"))
  mod_terms     <- paste0(mod_params, "*data.", moderators)
  formula       <- paste0(var_label, " := exp(", intercept_param, " + ",
                          paste(mod_terms, collapse = " + "), ")")
  variance_line <- paste0(latent_name, " ~~ ", var_label, "*", latent_name)

  updated_lines <- c(lines, variance_line, declarations, formula)
  return(paste(updated_lines, collapse = "\n"))
}


#' Build per-item screening models (moderated loadings and intercepts)
#'
#' For each item in the measurement model, creates a separate model string in
#' which that item's loading and intercept are each expressed as a linear
#' function of the supplied moderators.
#'
#' ## Identification strategies
#'
#' Two identification strategies are supported via `scale_loadings`:
#'
#' * **`scale_loadings = FALSE` (default):** The latent variance is fixed to 1
#'   (`LV ~~ 1*LV`) in every item model.  This is the recommended default
#'   because it allows every item — including the first — to have a freely
#'   moderated loading without any identification problem.  The returned models
#'   should be passed to [run_mxsem()] with `scale_loadings = FALSE` and
#'   `scale_latent_variances = TRUE` (or the variance constraint already
#'   present in the model string is sufficient).
#'
#' * **`scale_loadings = TRUE`:** The latent variance is free; scale
#'   identification relies on fixing the **first loading in the `=~` line**
#'   to 1 (handled automatically by [run_mxsem()] via `scale_loadings = TRUE`).
#'   **Special case — marker item:** when the item being tested for DIF *is*
#'   the first item (the current marker), fixing its loading to 1 and
#'   simultaneously allowing it to vary with the moderators is a
#'   contradiction.  The function therefore **swaps items 1 and 2** in the
#'   `=~` line for that specific model so that [run_mxsem()] fixes item 2's
#'   loading instead, leaving item 1's loading free for the DIF test.  The
#'   returned models should all be passed to [run_mxsem()] with
#'   `scale_loadings = TRUE` and `scale_latent_variances = FALSE`.
#'
#' **Parameter naming convention** used internally:
#' \itemize{
#'   \item `[facet]_load_[ii]` — definition-variable label for the target
#'     item's loading (zero-padded two-digit item index `ii`).
#'   \item `l_[facet]_0_[i]` — baseline (intercept) of the loading formula.
#'   \item `l_[facet]_[i]_[k]` — slope of the loading on moderator `k`.
#'   \item `int_[facet]_[item]` — baseline intercept of the target item.
#'   \item `int_[facet]_[item]_[kk]` — slope of the item intercept on
#'     moderator `k` (zero-padded two-digit moderator index `kk`).
#' }
#'
#' @param model_string Character. Baseline model syntax (one latent variable).
#' @param moderators Character vector of moderator names.  The names must
#'   match column names in the data frame that will be passed to [run_mxsem()].
#' @param scale_loadings Logical.  If `FALSE` (default) the latent variance is
#'   fixed to 1 in every item model (variance-based identification).  If
#'   `TRUE` identification is via the first loading; a marker-swap is applied
#'   automatically when the first item is the DIF target.
#'
#' @return A named list of character strings, one per item in the measurement
#'   model.  List names are the item names as they appear in `model_string`.
#'   When `scale_loadings = TRUE`, fit each returned model with
#'   `run_mxsem(..., scale_loadings = TRUE, scale_latent_variances = FALSE)`;
#'   when `scale_loadings = FALSE`, use
#'   `run_mxsem(..., scale_loadings = FALSE, scale_latent_variances = TRUE)`.
#'
#' @seealso [significant_item_moderators()], [run_mxsem()]
#' @export
moderate_loadings_and_intercepts <- function(model_string, moderators,
                                             scale_loadings = FALSE) {
  lines      <- unlist(strsplit(model_string, "\n"))
  lv_line    <- lines[grepl("=~", lines)][1]
  latent     <- sub("=~.*", "", gsub("\\s+", "", lv_line))
  rhs        <- sub(".*=~", "", lv_line)
  indicators <- gsub("\\s+", "", unlist(strsplit(rhs, "\\+")))

  # The first indicator is the marker item when scale_loadings = TRUE
  marker_item <- indicators[1]

  models <- list()
  for (i in seq_along(indicators)) {
    target_item    <- indicators[i]
    load_label     <- paste0(latent, "_load_", sprintf("%02d", i))
    load_intercept <- paste0("l_", latent, "_0_", i)
    load_params    <- paste0("l_", latent, "_", i, "_", seq_along(moderators))

    # ── Build the =~ line ────────────────────────────────────────────────────
    if (scale_loadings && i == 1L && length(indicators) >= 2L) {
      # Marker item is the DIF target: swap items 1 & 2 so run_mxsem() fixes
      # item 2 as the temporary marker, leaving item 1's loading free.
      new_loading_parts <- c(
        indicators[2L],                           # item 2 becomes temp marker
        paste0(load_label, "*", indicators[1L]),  # item 1 gets moderated label
        if (length(indicators) > 2L) indicators[seq(3L, length(indicators))]
      )
    } else {
      new_loading_parts <- ifelse(
        seq_along(indicators) == i,
        paste0(load_label, "*", indicators),
        indicators
      )
    }
    new_loading_line <- paste0(latent, " =~ ",
                               paste(new_loading_parts, collapse = " + "))

    loading_formula <- paste0(
      load_label, " := ", load_intercept, " + ",
      paste0("data.", moderators, " * ", load_params, collapse = " + ")
    )

    int_base       <- paste0("int_", latent, "_", target_item)
    int_mod_params <- paste0(int_base, "_", sprintf("%02d", seq_along(moderators)))
    intercept_line <- paste0(target_item, " ~ ", int_base, "*1")
    moderator_line <- paste0(target_item, " ~ ",
                             paste0(int_mod_params, "*", moderators, collapse = " + "))

    declarations <- c(paste0("!", load_intercept, ";"), paste0("!", load_params, ";"))

    new_lines <- lines
    new_lines[grepl("=~", lines)] <- new_loading_line
    new_lines <- new_lines[!grepl(paste0("^", target_item, " ~"), new_lines)]

    if (scale_loadings) {
      # Remove any hard-coded variance = 1 constraint from the baseline string;
      # identification is now via the first loading in the =~ line.
      new_lines  <- new_lines[!grepl(paste0("^", latent, "\\s*~~\\s*1\\*", latent),
                                     new_lines)]
      model_code <- c(new_lines, intercept_line, moderator_line,
                      declarations, loading_formula)
    } else {
      # Default: append an explicit variance = 1 constraint (safe for all items,
      # including the first).
      latent_variance_line <- paste0(latent, " ~~ 1*", latent)
      model_code <- c(new_lines, latent_variance_line, intercept_line,
                      moderator_line, declarations, loading_formula)
    }

    models[[target_item]] <- paste(model_code, collapse = "\n")
  }
  return(models)
}


#' Add covariance moderation to a multi-facet domain model string
#'
#' Locates existing `LV1 ~~ cov_LV1_LV2*LV2` lines and replaces each raw
#' covariance with a Fisher-z-transformed formula driven by the supplied
#' moderators.  Variance terms needed for the formula are added if not yet
#' present.
#'
#' @param model_string Character. Lavaan-style domain model syntax.
#' @param moderators Character vector of moderator names.
#'
#' @return A single character string with covariance moderation appended.
#'
#' @seealso [build_full_domain_models()], [significant_covar_moderators()]
#' @export
add_covariance_moderation <- function(model_string, moderators) {
  model_lines     <- strsplit(model_string, "\n")[[1]]
  cov_lines       <- grep("~~\\s*cov_", model_lines, value = TRUE)
  if (length(cov_lines) == 0) return(model_string)

  cov_mod_lines   <- c()
  added_variances <- character()

  for (line in cov_lines) {
    parts     <- strsplit(line, "~~")[[1]]
    lv1       <- trimws(parts[1])
    rhs_parts <- strsplit(trimws(parts[2]), "\\*")[[1]]
    cov_name  <- trimws(rhs_parts[1])
    lv2       <- trimws(rhs_parts[2])

    for (lv in c(lv1, lv2)) {
      var_label       <- paste0("var_", lv)
      already_defined <- any(grepl(paste0("^", lv, "\\s*~~\\s*", var_label, "\\*", lv),
                                   model_lines))
      if (!already_defined && !(var_label %in% added_variances)) {
        cov_mod_lines   <- c(cov_mod_lines, paste0(lv, " ~~ ", var_label, "*", lv))
        added_variances <- c(added_variances, var_label)
      }
    }

    rho_label <- paste0("rho_", cov_name)
    r_params  <- paste0("r", seq_along(moderators), "_", cov_name)
    r_defs    <- c(paste0("!r0_", cov_name, ";"), paste0("!", r_params, ";"))
    rho_terms <- paste0(r_params, "*data.", moderators)
    rho_expr  <- paste(c(paste0("r0_", cov_name), rho_terms), collapse = " + ")

    cov_mod_lines <- c(
      cov_mod_lines,
      r_defs,
      paste0(cov_name, " := sqrt(var_", lv1, " * var_", lv2, ") * ",
             "(exp(2 * ", rho_label, ") - 1)/(exp(2 * ", rho_label, ") + 1)"),
      paste0(rho_label, " := ", rho_expr)
    )
  }

  paste(c(model_lines, cov_mod_lines), collapse = "\n")
}


#' Add latent variance moderation from a significance data frame (full-model step)
#'
#' Reads the significant variance moderators for `facet_name` from
#' `all_sig_mods` and appends the corresponding definition-variable formula
#' to `model_string`.  Unlike [add_latent_variance_moderation()], this
#' function takes its moderator list from the screening results rather than
#' from a manually supplied vector.
#'
#' @param model_string Character. Lavaan-style model syntax.
#' @param facet_name Character. Name of the facet (latent variable).
#' @param all_sig_mods Data frame returned by [combine_sig_moderators()] or
#'   [extract_all_significant_moderators()].
#'
#' @return A single character string, or the original `model_string` unchanged
#'   if no significant variance moderators are found.
#'
#' @seealso [build_full_domain_models()]
#' @export
add_latent_variance_moderation_from_df <- function(model_string, facet_name,
                                                    all_sig_mods) {
  var_mods <- all_sig_mods[
    all_sig_mods$facet == facet_name &
    all_sig_mods$type  == "variance"  &
    all_sig_mods$significant_corrected == TRUE, ]

  if (nrow(var_mods) == 0) return(model_string)

  moderators <- unique(var_mods$moderator)
  lines      <- unlist(strsplit(model_string, "\n"))

  lv_line <- lines[grepl("=~", lines)][1]
  latent  <- sub("=~.*", "", gsub("\\s+", "", lv_line))

  variance_label    <- paste0("var_", facet_name)
  mod_params        <- paste0("v_", facet_name, "_", seq_along(moderators))

  lines <- lines[!grepl(paste0("^", latent, "~~"), lines)]
  lines <- lines[!grepl("^!", lines)]
  lines <- lines[!grepl(paste0("^", variance_label, " :="), lines)]

  v0_param          <- paste0("v0_", facet_name)
  variance_link     <- paste0(latent, " ~~ ", variance_label, "*", latent)
  declaration_lines <- c(paste0("!", v0_param, ";"), paste0("!", mod_params, ";"))
  linear_formula    <- paste0(c(v0_param, paste0(mod_params, "*data.", moderators)),
                              collapse = " + ")
  variance_formula  <- paste0(variance_label, " := exp(", linear_formula, ")")

  updated_lines <- c(lines, variance_link, declaration_lines, variance_formula)
  return(paste(updated_lines, collapse = "\n"))
}


#' Add item-level moderation for one facet (full-model step)
#'
#' Modifies the measurement line for `facet_name` so that items with
#' significant loading or intercept moderation receive the appropriate labels
#' and definition formulas.  Anchor items (in `anchor_items`) are left
#' unchanged.
#'
#' @param model_string Character. Lavaan-style model syntax.
#' @param facet_name Character. Name of the facet.
#' @param all_sig_mods Data frame of all significant moderations (see
#'   [combine_sig_moderators()]).
#' @param anchor_items Character vector of anchor item identifiers in the form
#'   `"facet_name.item_name"`.
#'
#' @return A single character string.
#'
#' @seealso [add_all_moderation_single_facet()], [build_full_domain_models()]
#' @export
add_item_moderation_single_facet <- function(model_string, facet_name,
                                             all_sig_mods, anchor_items) {
  lines      <- unlist(strsplit(model_string, "\n"))
  lv_line    <- lines[grepl("=~", lines)][1]
  latent     <- sub("=~.*", "", gsub("\\s+", "", lv_line))
  rhs        <- sub(".*=~", "", lv_line)
  indicators <- trimws(unlist(strsplit(rhs, "\\+")))

  loading_parts <- indicators
  extra_lines   <- c()

  for (i in seq_along(indicators)) {
    item      <- indicators[i]
    full_item <- paste0(facet_name, ".", item)
    if (full_item %in% anchor_items) next

    load_mods <- all_sig_mods[
      all_sig_mods$facet == facet_name &
      all_sig_mods$item  == item        &
      all_sig_mods$type  == "item"      &
      grepl("^l", all_sig_mods$param_name) &
      all_sig_mods$significant_corrected == TRUE, ]

    if (nrow(load_mods) > 0) {
      moderators   <- unique(load_mods$moderator)
      load_label   <- paste0(facet_name, "_load_", sprintf("%02d", i))
      param_names  <- paste0("l_", facet_name, "_", i, "_", seq_along(moderators))
      loading_parts[i] <- paste0(load_label, "*", item)
      declaration_lines <- paste0("!", param_names, ";")
      formula_line <- paste0(load_label, " := ",
                             paste0(param_names, "*data.", moderators, collapse = " + "))
      extra_lines  <- c(extra_lines, declaration_lines, formula_line)
    }

    int_mods <- all_sig_mods[
      all_sig_mods$facet == facet_name &
      all_sig_mods$item  == item        &
      all_sig_mods$type  == "item"      &
      grepl("^int_", all_sig_mods$param_name) &
      all_sig_mods$significant_corrected == TRUE, ]

    if (nrow(int_mods) > 0) {
      moderators   <- unique(int_mods$moderator)
      param_labels <- paste0("int_", facet_name, "_", item, "_",
                             sprintf("%02d", seq_along(moderators)))
      intercept_line <- paste0(item, " ~ ",
                               paste0(param_labels, "*", moderators, collapse = " + "))
      extra_lines  <- c(extra_lines, intercept_line)
    }
  }

  lines[grepl("=~", lines)] <- paste0(latent, " =~ ",
                                      paste(loading_parts, collapse = " + "))
  full_model <- c(lines, extra_lines)
  return(paste(full_model, collapse = "\n"))
}


#' Add all significant moderations for one facet (primary full-model builder)
#'
#' The main workhorse for building per-facet model syntax in the full domain
#' model.  It handles the identification strategy switch (screening: variance
#' fixed to 1; full model: first anchor loading fixed to 1), excludes all
#' anchor items from moderation labels, and builds definition formulas for
#' loading and intercept moderations.
#'
#' This is the function that should be used in production pipelines; it is
#' called internally by [build_full_domain_models()].
#'
#' @param model_string Character. Single-facet baseline model syntax
#'   (character string or vector of lines).
#' @param facet_name Character. Name of the facet (must match the latent
#'   variable name in `model_string`).
#' @param all_sig_mods Data frame from [combine_sig_moderators()].
#' @param anchor_items Character vector of anchor identifiers
#'   (`"facet_name.item_name"` format).
#' @param moderators Optional character vector of all moderator names in their
#'   original order.  Used to compute correct numeric indices for parameter
#'   names when `all_sig_mods` contains a subset of moderators.
#'
#' @return A character vector of model syntax lines (ready to be
#'   `paste(collapse = "\n\n")`d together with other facets).
#'
#' @seealso [build_full_domain_models()]
#' @export
add_all_moderation_single_facet <- function(model_string, facet_name,
                                             all_sig_mods, anchor_items,
                                             moderators = NULL) {
  # Normalise: single \n-string → line vector
  model_string <- .normalise_model_string(model_string)

  facet_mods <- subset(all_sig_mods, facet == facet_name & significant_corrected)
  if (nrow(facet_mods) == 0) return(model_string)

  measurement_line <- grep(paste0("^", facet_name, " =~"), model_string, value = TRUE)
  items <- unlist(strsplit(gsub(paste0("^", facet_name, " =~ "), "",
                                measurement_line), " \\+ "))
  model_string <- model_string[!grepl(paste0("^", facet_name, " =~"), model_string)]

  # Determine first anchor item for scale identification
  anchor_item_for_scaling <- NULL
  for (item in items) {
    if (paste0(facet_name, ".", item) %in% anchor_items) {
      anchor_item_for_scaling <- item
      break
    }
  }
  if (is.null(anchor_item_for_scaling)) {
    stop(paste0("No anchor item found for facet ", facet_name))
  }

  # Build measurement line
  # BUG FIX: All anchor items (not just the first) must remain label-free.
  new_items <- sapply(seq_along(items), function(i) {
    item      <- items[i]
    full_item <- paste0(facet_name, ".", item)
    if (item == anchor_item_for_scaling) {
      # First anchor: loading fixed to 1 (scale identification)
      paste0("1*", item)
    } else if (full_item %in% anchor_items) {
      # Additional anchors: free but unmoderated
      item
    } else {
      # Non-anchor items: label only when actually moderated
      is_moderated <- any(facet_mods$type == "loading" & facet_mods$item == item)
      if (is_moderated)
        paste0(paste0(facet_name, "_load_", sprintf("%02d", i)), "*", item)
      else
        item
    }
  })

  measurement_line_new <- paste0(facet_name, " =~ ", paste(new_items, collapse = " + "))
  lines <- c(measurement_line_new, paste0(facet_name, " ~ m_", facet_name, "*1"))

  # Latent mean moderation
  latent_mean <- subset(facet_mods, type == "mean")
  if (nrow(latent_mean) > 0) {
    mean_parts <- character()
    for (i in seq_len(nrow(latent_mean))) {
      row        <- latent_mean[i, ]
      mod_index  <- if (!is.null(moderators)) which(moderators == row$moderator) else i
      param_name <- paste0("m_", facet_name, "_", sprintf("%02d", mod_index))
      mean_parts <- c(mean_parts, paste0(param_name, "*", row$moderator))
    }
    lines <- c(lines, paste0(facet_name, " ~ ", paste(mean_parts, collapse = " + ")))
  }

  # Latent variance moderation
  latent_var <- subset(facet_mods, type == "variance")
  if (nrow(latent_var) > 0) {
    lines <- c(lines, paste0(facet_name, " ~~ var_", facet_name, "*", facet_name))
    v0_param    <- paste0("v0_", facet_name)
    var_labels  <- paste0("!", v0_param, ";")
    var_formula <- v0_param   # free intercept: log-variance at moderator = 0
    for (i in seq_len(nrow(latent_var))) {
      row        <- latent_var[i, ]
      mod_index  <- if (!is.null(moderators)) which(moderators == row$moderator) else i
      param_name <- paste0("v_", facet_name, "_", mod_index)
      mod_label  <- paste0("data.", row$moderator)
      var_labels  <- c(var_labels,  paste0("!", param_name, ";"))
      var_formula <- c(var_formula, paste0(param_name, "*", mod_label))
    }
    var_line <- paste0("var_", facet_name, " := exp(",
                       paste(var_formula, collapse = " + "), ")")
    lines <- c(lines, var_labels, var_line)
  }

  # Loading moderation (non-anchor items only)
  loading_mods <- subset(facet_mods,
    type == "loading" & !(paste0(facet, ".", item) %in% anchor_items))
  if (nrow(loading_mods) > 0) {
    for (item_name in unique(loading_mods$item)) {
      item_index     <- which(items == item_name)
      label          <- paste0(facet_name, "_load_", sprintf("%02d", item_index))
      baseline_param <- paste0("l_", facet_name, "_0_", item_index)
      lines          <- c(lines, paste0("!", baseline_param, ";"))
      item_mods      <- subset(loading_mods, item == item_name)
      mod_terms      <- character()
      for (i in seq_len(nrow(item_mods))) {
        row        <- item_mods[i, ]
        mod_index  <- if (!is.null(moderators)) which(moderators == row$moderator) else i
        param_name <- paste0("l_", facet_name, "_", item_index, "_", mod_index)
        lines      <- c(lines, paste0("!", param_name, ";"))
        mod_terms  <- c(mod_terms, paste0("data.", row$moderator, " * ", param_name))
      }
      lines <- c(lines, paste0(label, " := ",
                               paste(c(baseline_param, mod_terms), collapse = " + ")))
    }
  }

  # Intercept moderation (non-anchor items only)
  intercept_mods <- subset(facet_mods,
    type == "intercept" & !(paste0(facet, ".", item) %in% anchor_items))
  if (nrow(intercept_mods) > 0) {
    for (item_name in unique(intercept_mods$item)) {
      item_rows  <- subset(intercept_mods, item == item_name)
      base_param <- paste0("int_", facet_name, "_", item_name)
      lines      <- c(lines, paste0(item_name, " ~ ", base_param, "*1"))
      mod_terms  <- character()
      for (i in seq_len(nrow(item_rows))) {
        row        <- item_rows[i, ]
        mod_index  <- if (!is.null(moderators)) which(moderators == row$moderator) else i
        param_name <- paste0("int_", facet_name, "_", item_name, "_",
                             sprintf("%02d", mod_index))
        mod_terms  <- c(mod_terms, paste0(param_name, "*", row$moderator))
      }
      if (length(mod_terms) > 0)
        lines <- c(lines, paste0(item_name, " ~ ", paste(mod_terms, collapse = " + ")))
    }
  }

  return(lines)
}


#' Build a complete single-facet moderated SEM model string
#'
#' Alternative full-model builder for single-latent-variable models.  Unlike
#' [add_all_moderation_single_facet()], this function works directly with
#' item names (not `facet.item` composite keys) in `anchor_items`, and
#' returns a character vector of lines.
#'
#' @param model_string Character. Baseline model syntax.
#' @param facet_name Character. Latent variable name.
#' @param all_sig_mods Data frame of significant moderations.
#' @param anchor_items Character vector of bare item names (without
#'   `facet_name.` prefix) that serve as anchors.
#' @param moderators Optional character vector of all moderators in order.
#'
#' @return A character vector of model syntax lines.
#'
#' @seealso [add_all_moderation_single_facet()]
#' @export
build_moderated_sem_one_latent <- function(model_string, facet_name,
                                           all_sig_mods, anchor_items,
                                           moderators = NULL) {
  # Normalise: single \n-string → line vector
  model_string <- .normalise_model_string(model_string)

  meas_line <- grep(paste0("^", facet_name, " =~"), model_string, value = TRUE)
  if (length(meas_line) == 0) stop("Measurement line not found.")

  items <- unlist(strsplit(gsub(paste0("^", facet_name, " =~\\s*"), "",
                                meas_line), "\\s*\\+\\s*"))
  anchors       <- anchor_items[anchor_items %in% items]
  anchor1       <- if (length(anchors) > 0) anchors[1] else NULL
  anchor_others <- setdiff(anchors, anchor1)

  facet_mods <- subset(all_sig_mods, item %in% items | is.na(item))
  load_mods  <- subset(facet_mods, type == "loading")
  int_mods   <- subset(facet_mods, type == "intercept")
  mean_mods  <- subset(facet_mods, type == "mean")

  load_items       <- unique(load_mods$item)
  int_items        <- unique(int_mods$item)
  all_mod_items    <- union(load_items, int_items)
  non_anchor_mod   <- intersect(setdiff(items, anchors), all_mod_items)
  non_anchor_unmod <- setdiff(setdiff(items, anchors), all_mod_items)

  item_index      <- setNames(seq_along(items), items)
  meas_terms      <- character()
  label_lines     <- character()
  def_lines       <- character()
  intercept_lines <- character()

  if (!is.null(anchor1))         meas_terms <- c(meas_terms, paste0("1*", anchor1))
  if (length(anchor_others) > 0) meas_terms <- c(meas_terms, anchor_others)
  meas_terms <- c(meas_terms, non_anchor_unmod)

  for (item in non_anchor_mod) {
    idx <- item_index[item]

    item_load_mods <- load_mods[load_mods$item == item, ]
    if (nrow(item_load_mods) > 0) {
      base_label  <- sprintf("l_%s_0_%02d", facet_name, idx)
      param_label <- sprintf("%s_load_%02d", facet_name, idx)
      label_lines <- c(label_lines, paste0("!", base_label, ";"))
      mod_terms   <- character()
      seen_mods   <- character()
      for (i in seq_len(nrow(item_load_mods))) {
        mod_name <- item_load_mods$moderator[i]
        if (mod_name %in% seen_mods) next
        seen_mods   <- c(seen_mods, mod_name)
        mod_idx     <- if (!is.null(moderators)) which(moderators == mod_name) else i
        param       <- sprintf("l_%s_%02d_%d", facet_name, idx, mod_idx)
        label_lines <- c(label_lines, paste0("!", param, ";"))
        mod_terms   <- c(mod_terms, paste0("data.", mod_name, " * ", param))
      }
      def_lines  <- c(def_lines, paste0(param_label, " := ",
                                        paste(c(base_label, mod_terms), collapse = " + ")))
      meas_terms <- c(meas_terms, paste0(param_label, "*", item))
    } else {
      meas_terms <- c(meas_terms, item)
    }

    item_int_mods <- int_mods[int_mods$item == item, ]
    if (nrow(item_int_mods) > 0) {
      base_param      <- sprintf("int_%s_%s", facet_name, item)
      intercept_lines <- c(intercept_lines, paste0(item, " ~ ", base_param, "*1"))
      mod_terms <- character()
      seen_mods <- character()
      for (i in seq_len(nrow(item_int_mods))) {
        mod_name <- item_int_mods$moderator[i]
        if (mod_name %in% seen_mods) next
        seen_mods <- c(seen_mods, mod_name)
        mod_idx   <- if (!is.null(moderators)) which(moderators == mod_name) else i
        param     <- sprintf("int_%s_%s_%02d", facet_name, item, mod_idx)
        mod_terms <- c(mod_terms, paste0(param, "*", mod_name))
      }
      intercept_lines <- c(intercept_lines,
                           paste0(item, " ~ ", paste(mod_terms, collapse = " + ")))
    }
  }

  # Latent mean moderation
  mean_lines <- character()
  if (nrow(mean_mods) > 0) {
    mean_lines <- c(mean_lines, paste0(facet_name, " ~ m_", facet_name, "*1"))
    mod_terms  <- character()
    seen_mods  <- character()
    for (i in seq_len(nrow(mean_mods))) {
      mod_name <- mean_mods$moderator[i]
      if (mod_name %in% seen_mods) next
      seen_mods <- c(seen_mods, mod_name)
      mod_idx   <- if (!is.null(moderators)) which(moderators == mod_name) else i
      param     <- sprintf("m_%s_%02d", facet_name, mod_idx)
      mod_terms <- c(mod_terms, paste0(param, "*", mod_name))
    }
    if (length(mod_terms) > 0)
      mean_lines <- c(mean_lines,
                      paste0(facet_name, " ~ ", paste(mod_terms, collapse = " + ")))
  }

  meas_line_new <- paste0(facet_name, " =~ ", paste(meas_terms, collapse = " + "))
  return(c(meas_line_new, mean_lines, label_lines, def_lines, intercept_lines))
}
