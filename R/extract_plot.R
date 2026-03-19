# =============================================================================
# R/extract_plot.R
# Functions for extracting and plotting conditional moderation effects.
# =============================================================================


#' Extract a summary table of all moderation effects
#'
#' Collects every moderation-slope parameter from `param_df` (the `$param_df`
#' slot of a [run_mxsem()] result) into a single, human-readable table.  Each
#' row represents the partial effect of one moderator on one target parameter
#' (loading, intercept, latent mean, latent variance, or covariance).
#'
#' Baseline parameters (intercepts of the formulas, e.g. `l_P_1_0_3`,
#' `v0_P_1`, `m_P_1`) are excluded; only the slope coefficients are returned.
#'
#' @param param_df Data frame from the `$param_df` slot of [run_mxsem()].
#' @param moderators Optional character vector of moderator names in the order
#'   used during model fitting.  Recovered from `attr(param_df, "moderators")`
#'   if `NULL`.
#' @param alpha Significance level.  Default `0.05`.
#' @param correction_method P-value correction method passed to
#'   [stats::p.adjust()].  Default `"none"`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`param_type`}{One of `"latent_mean"`, `"latent_variance"`,
#'     `"loading"`, `"intercept"`, `"covariance"`.}
#'   \item{`facet`}{Latent variable name (e.g. `"P_1"`).}
#'   \item{`item`}{Item name for `"loading"` / `"intercept"` rows; `NA`
#'     otherwise.}
#'   \item{`moderator`}{Name of the moderator variable.}
#'   \item{`Estimate`}{Regression coefficient (unstandardised).}
#'   \item{`Std.Error`}{Standard error.}
#'   \item{`z`}{z-statistic.}
#'   \item{`p`}{Two-sided p-value.}
#'   \item{`p_corrected`}{Corrected p-value (via `correction_method`).}
#'   \item{`significant`}{Logical: `p_corrected <= alpha`.}
#' }
#'
#' @seealso [compute_conditional_parameters()], [run_mxsem()]
#' @export
extract_moderation_table <- function(param_df, moderators = NULL,
                                     alpha             = 0.05,
                                     correction_method = "none") {
  if (is.null(moderators)) moderators <- attr(param_df, "moderators")

  rows_list <- list()

  # ── 1. Latent mean moderations ─────────────────────────────────────────────
  # Params: m_Facet_kk stored in the A matrix (regression on moderator)
  lm_rows <- param_df[!is.na(param_df$type)   & param_df$type   == "latent_mean" &
                      !is.na(param_df$matrix) & param_df$matrix == "A", ]
  if (nrow(lm_rows) > 0) {
    lm_rows$param_type <- "latent_mean"
    lm_rows$item       <- NA_character_
    lm_rows$moderator  <- lm_rows$col
    rows_list[["mean"]] <- lm_rows[, c("param_type", "facet", "item",
                                       "moderator", "Estimate", "Std.Error",
                                       "z", "p")]
  }

  # ── 2. Latent variance moderations ────────────────────────────────────────
  # Params: v_Facet_k in new_parameters; exclude the v0_ baseline intercept
  v_rows <- param_df[!is.na(param_df$type) & param_df$type == "variance" &
                     !is.na(param_df$name) & !grepl("^v0_", param_df$name), ]
  if (nrow(v_rows) > 0) {
    v_rows$param_type <- "latent_variance"
    v_rows$item       <- NA_character_
    v_rows$moderator  <- v_rows$col
    rows_list[["variance"]] <- v_rows[, c("param_type", "facet", "item",
                                          "moderator", "Estimate", "Std.Error",
                                          "z", "p")]
  }

  # ── 3. Item loading moderations ────────────────────────────────────────────
  # Params: l_Facet_i_k in new_parameters; exclude l_Facet_0_i (baseline)
  l_rows <- param_df[!is.na(param_df$type)   & param_df$type   == "loading" &
                     !is.na(param_df$matrix) & param_df$matrix == "new_parameters" &
                     !grepl("_0_", param_df$name), ]
  if (nrow(l_rows) > 0) {
    l_rows$param_type <- "loading"
    l_rows$item       <- l_rows$row
    l_rows$moderator  <- l_rows$col
    rows_list[["loading"]] <- l_rows[, c("param_type", "facet", "item",
                                         "moderator", "Estimate", "Std.Error",
                                         "z", "p")]
  }

  # ── 4. Item intercept moderations ──────────────────────────────────────────
  # Params: int_Facet_item_kk in the A matrix
  int_rows <- param_df[!is.na(param_df$type)   & param_df$type   == "intercept" &
                       !is.na(param_df$matrix) & param_df$matrix == "A", ]
  if (nrow(int_rows) > 0) {
    int_rows$param_type <- "intercept"
    int_rows$item       <- int_rows$row
    int_rows$moderator  <- int_rows$col
    rows_list[["intercept"]] <- int_rows[, c("param_type", "facet", "item",
                                              "moderator", "Estimate",
                                              "Std.Error", "z", "p")]
  }

  # ── 5. Latent covariance moderations ──────────────────────────────────────
  # Params: rk_cov_Facet1_Facet2 (type "covariance")
  cov_rows <- param_df[!is.na(param_df$type) & param_df$type == "covariance", ]
  if (nrow(cov_rows) > 0) {
    cov_rows$param_type <- "covariance"
    cov_rows$item       <- NA_character_
    cov_rows$moderator  <- cov_rows$col
    # Facet label = the cov_* part of the parameter name
    cov_rows$facet <- sub("^r[0-9]+_", "", cov_rows$name)
    rows_list[["covariance"]] <- cov_rows[, c("param_type", "facet", "item",
                                               "moderator", "Estimate",
                                               "Std.Error", "z", "p")]
  }

  if (length(rows_list) == 0) {
    return(data.frame(
      param_type = character(), facet = character(), item = character(),
      moderator  = character(), Estimate = numeric(), Std.Error = numeric(),
      z = numeric(), p = numeric(), p_corrected = numeric(),
      significant = logical(), stringsAsFactors = FALSE
    ))
  }

  combined          <- do.call(rbind, rows_list)
  rownames(combined) <- NULL
  combined$p_corrected <- stats::p.adjust(combined$p, method = correction_method)
  combined$significant <- !is.na(combined$p_corrected) &
    combined$p_corrected <= alpha

  ord <- order(combined$param_type,
               combined$facet,
               ifelse(is.na(combined$item), "", combined$item))
  combined[ord, ]
}


#' Compute parameter values conditional on a moderator grid
#'
#' For each moderated parameter in `param_df`, evaluates the parameter value
#' at every row of `moderator_grid`.  This produces a long-format data frame
#' that is directly suitable for plotting conditional (DIF) curves.
#'
#' **Interaction plots.** Passing an `expand.grid()` of a continuous and a
#' binary moderator (e.g. `expand.grid(age_z = seq(-2, 2, .5), sex = c(0, 1))`)
#' produces separate trajectories per binary level, which can then be rendered
#' as an interaction graph via [plot_conditional_parameters()] or `ggplot2`.
#'
#' **Formula per parameter type:**
#' \itemize{
#'   \item Latent mean: \eqn{m_0 + \sum_k m_k \cdot \text{mod}_k}
#'   \item Loading:     \eqn{l_0 + \sum_k l_k \cdot \text{mod}_k}
#'   \item Intercept:   \eqn{\tau_0 + \sum_k \tau_k \cdot \text{mod}_k}
#'   \item Variance:    \eqn{\exp(v_0 + \sum_k v_k \cdot \text{mod}_k)}
#' }
#'
#' @param param_df Data frame from the `$param_df` slot of [run_mxsem()].
#' @param moderator_grid A `data.frame` whose columns are the moderator
#'   variables.  Typically produced with [base::expand.grid()], e.g.
#'   `expand.grid(age_z = seq(-2, 2, 0.5), sex = c(0, 1))`.
#' @param moderators Optional character vector.  Recovered from
#'   `attr(param_df, "moderators")` if `NULL`.
#' @param param_types Character vector restricting which parameter types are
#'   included.  Any subset of
#'   `c("latent_mean", "latent_variance", "loading", "intercept")`.
#'   Default: all four.
#'
#' @return A long-format `data.frame` combining all columns of
#'   `moderator_grid` with:
#' \describe{
#'   \item{`parameter`}{Unique label, e.g. `"P_1_mean"`, `"p1_i5_intercept"`,
#'     `"p3_i7_loading"`.}
#'   \item{`param_type`}{Parameter type.}
#'   \item{`facet`}{Latent variable name.}
#'   \item{`item`}{Item name or `NA`.}
#'   \item{`value`}{Conditional parameter value at the given moderator values.}
#' }
#'
#' @examples
#' \dontrun{
#' grid <- expand.grid(age_z = seq(-2, 2, 0.5), sex = c(0, 1))
#' cond <- compute_conditional_parameters(fit_domain$param_df, grid)
#'
#' # Interaction plot for the loading of p3_i7
#' plot_conditional_parameters(cond, x_var = "age_z", group_var = "sex",
#'                              param_types = "loading")
#' }
#'
#' @seealso [extract_moderation_table()], [plot_conditional_parameters()]
#' @export
compute_conditional_parameters <- function(param_df,
                                           moderator_grid,
                                           moderators   = NULL,
                                           param_types  = c("latent_mean",
                                                            "latent_variance",
                                                            "loading",
                                                            "intercept")) {
  if (is.null(moderators)) moderators <- attr(param_df, "moderators")

  # Helper: evaluate intercept + slopes at all grid rows
  eval_at_grid <- function(intercept, slope_df, grid) {
    # slope_df: data.frame with columns "moderator" and "Estimate"
    pred <- rep(as.numeric(intercept), nrow(grid))
    if (nrow(slope_df) == 0) return(pred)
    for (i in seq_len(nrow(slope_df))) {
      mod <- slope_df$moderator[i]
      if (!is.na(mod) && mod %in% names(grid)) {
        pred <- pred + slope_df$Estimate[i] * grid[[mod]]
      }
    }
    pred
  }

  results <- list()

  # ── 1. Latent means ─────────────────────────────────────────────────────────
  if ("latent_mean" %in% param_types) {
    facets <- unique(param_df$facet[!is.na(param_df$type) &
                                      param_df$type == "latent_mean"])
    for (facet in na.omit(facets)) {
      fp <- param_df[!is.na(param_df$facet) & param_df$facet == facet &
                     !is.na(param_df$type)  & param_df$type  == "latent_mean", ]
      base_row   <- fp[!is.na(fp$matrix) & fp$matrix == "M", ]
      base_val   <- if (nrow(base_row) == 1) base_row$Estimate else 0
      slope_rows <- fp[!is.na(fp$matrix) & fp$matrix == "A", ]
      slope_df   <- data.frame(moderator = slope_rows$col,
                               Estimate  = slope_rows$Estimate,
                               stringsAsFactors = FALSE)
      if (nrow(slope_df) == 0) next
      out            <- moderator_grid
      out$parameter  <- paste0(facet, "_mean")
      out$param_type <- "latent_mean"
      out$facet      <- facet
      out$item       <- NA_character_
      out$value      <- eval_at_grid(base_val, slope_df, moderator_grid)
      results[[paste0(facet, "_mean")]] <- out
    }
  }

  # ── 2. Latent variances ─────────────────────────────────────────────────────
  if ("latent_variance" %in% param_types) {
    v_params <- param_df[!is.na(param_df$type) & param_df$type == "variance" &
                         !is.na(param_df$name) & !grepl("^v0_", param_df$name), ]
    facets   <- unique(v_params$facet[!is.na(v_params$facet)])
    for (facet in facets) {
      fp      <- v_params[!is.na(v_params$facet) & v_params$facet == facet, ]
      v0_name <- paste0("v0_", facet)
      v0_row  <- param_df[!is.na(param_df$name) & param_df$name == v0_name, ]
      v0_val  <- if (nrow(v0_row) == 1) v0_row$Estimate else 0
      slope_df <- data.frame(moderator = fp$col,
                              Estimate  = fp$Estimate,
                              stringsAsFactors = FALSE)
      if (nrow(slope_df) == 0) next
      out            <- moderator_grid
      out$parameter  <- paste0(facet, "_variance")
      out$param_type <- "latent_variance"
      out$facet      <- facet
      out$item       <- NA_character_
      out$value      <- exp(eval_at_grid(v0_val, slope_df, moderator_grid))
      results[[paste0(facet, "_variance")]] <- out
    }
  }

  # ── 3. Item loadings ────────────────────────────────────────────────────────
  if ("loading" %in% param_types) {
    l_params <- param_df[!is.na(param_df$type) & param_df$type == "loading", ]
    items    <- unique(l_params$row[!is.na(l_params$row)])
    for (item in items) {
      ip       <- l_params[!is.na(l_params$row) & l_params$row == item, ]
      base_row <- ip[grepl("_0_", ip$name), ]
      if (nrow(base_row) != 1) next
      base_val  <- base_row$Estimate
      facet_val <- if (!is.na(base_row$facet[1])) base_row$facet[1] else NA_character_
      slope_rows <- ip[!grepl("_0_", ip$name), ]
      slope_df   <- data.frame(moderator = slope_rows$col,
                               Estimate  = slope_rows$Estimate,
                               stringsAsFactors = FALSE)
      if (nrow(slope_df) == 0) next
      out            <- moderator_grid
      out$parameter  <- paste0(item, "_loading")
      out$param_type <- "loading"
      out$facet      <- facet_val
      out$item       <- item
      out$value      <- eval_at_grid(base_val, slope_df, moderator_grid)
      results[[paste0(item, "_loading")]] <- out
    }
  }

  # ── 4. Item intercepts ──────────────────────────────────────────────────────
  if ("intercept" %in% param_types) {
    int_params <- param_df[!is.na(param_df$type) & param_df$type == "intercept", ]
    items      <- unique(int_params$row[!is.na(int_params$row)])
    for (item in items) {
      ip        <- int_params[!is.na(int_params$row) & int_params$row == item, ]
      # Baseline: M matrix (item mean intercept)
      base_row  <- ip[!is.na(ip$matrix) & ip$matrix == "M", ]
      base_val  <- if (nrow(base_row) == 1) base_row$Estimate else 0
      facet_val <- unique(ip$facet[!is.na(ip$facet)])
      # Slopes: A matrix
      slope_rows <- ip[!is.na(ip$matrix) & ip$matrix == "A", ]
      slope_df   <- data.frame(moderator = slope_rows$col,
                               Estimate  = slope_rows$Estimate,
                               stringsAsFactors = FALSE)
      if (nrow(slope_df) == 0) next
      out            <- moderator_grid
      out$parameter  <- paste0(item, "_intercept")
      out$param_type <- "intercept"
      out$facet      <- if (length(facet_val) > 0) facet_val[1] else NA_character_
      out$item       <- item
      out$value      <- eval_at_grid(base_val, slope_df, moderator_grid)
      results[[paste0(item, "_intercept")]] <- out
    }
  }

  if (length(results) == 0) return(data.frame())
  out_df            <- do.call(rbind, results)
  rownames(out_df)  <- NULL
  out_df
}


#' Plot conditional parameter values across a moderator grid
#'
#' Renders the output of [compute_conditional_parameters()] as line plots —
#' one panel per moderated parameter.  If **ggplot2** is installed, a list of
#' `ggplot` objects is returned; otherwise base-R graphics are drawn.
#'
#' **Interaction / Zweifachinteraktions-plots.**  Specify `x_var` as the
#' continuous moderator and `group_var` as the binary (or categorical) one.
#' Each level of `group_var` is drawn as a separate line, producing a
#' classic two-way interaction graph.
#'
#' @param cond_params Data frame returned by [compute_conditional_parameters()].
#' @param x_var Character. Column in `cond_params` to use as x-axis (typically
#'   a continuous moderator, e.g. `"age_z"`).
#' @param group_var Optional character. Column to use for separate lines
#'   (e.g. `"sex"`).  `NULL` draws a single line per parameter.
#' @param param_types Optional character vector to restrict which parameter
#'   types are plotted.  Defaults to all types present in `cond_params`.
#' @param xlab Character. X-axis label.  Defaults to `x_var`.
#' @param ylab Character. Y-axis label.  Default `"Parameter value"`.
#' @param main_prefix Character. Prefix prepended to each panel title.
#'   Default `""`.
#' @param ncol Integer. Number of plot columns in the base-R layout.
#'   Default `3`.
#' @param ... Further arguments passed to `plot()` (base-R mode only).
#'
#' @return In ggplot2 mode: a named list of `ggplot` objects (one per
#'   parameter), returned invisibly.  Each plot is also printed immediately
#'   so it appears in RStudio / a graphics window without an explicit
#'   `print()` call.  In base-R mode: `NULL` invisibly; plots are drawn as
#'   a side-effect.
#'
#' @seealso [compute_conditional_parameters()], [extract_moderation_table()]
#' @export
plot_conditional_parameters <- function(cond_params,
                                        x_var,
                                        group_var    = NULL,
                                        param_types  = NULL,
                                        xlab         = x_var,
                                        ylab         = "Parameter value",
                                        main_prefix  = "",
                                        ncol         = 3L,
                                        ...) {
  if (is.null(param_types)) param_types <- unique(cond_params$param_type)
  cond_params <- cond_params[cond_params$param_type %in% param_types, ]
  if (nrow(cond_params) == 0) {
    message("No conditional parameters to plot for the selected param_types.")
    return(invisible(NULL))
  }

  parameters  <- unique(cond_params$parameter)
  use_ggplot  <- requireNamespace("ggplot2", quietly = TRUE)

  # ── ggplot2 path ────────────────────────────────────────────────────────────
  if (use_ggplot) {
    make_one <- function(param) {
      df    <- cond_params[cond_params$parameter == param, ]
      title <- paste0(main_prefix, param)
      if (!is.null(group_var) && group_var %in% names(df)) {
        df[[group_var]] <- factor(df[[group_var]])
        aes_mapping <- ggplot2::aes(
          x        = .data[[x_var]],
          y        = .data[["value"]],
          colour   = .data[[group_var]],
          linetype = .data[[group_var]],
          group    = .data[[group_var]]
        )
      } else {
        aes_mapping <- ggplot2::aes(
          x = .data[[x_var]],
          y = .data[["value"]]
        )
      }
      ggplot2::ggplot(df, aes_mapping) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::labs(title  = title, x = xlab, y = ylab,
                      colour = group_var, linetype = group_var) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    }
    plot_list        <- lapply(parameters, make_one)
    names(plot_list) <- parameters
    for (p in plot_list) print(p)   # auto-print each plot (required outside interactive sessions)
    return(invisible(plot_list))
  }

  # ── Base-R path ─────────────────────────────────────────────────────────────
  n_plots <- length(parameters)
  n_col   <- min(as.integer(ncol), n_plots)
  n_row   <- ceiling(n_plots / n_col)
  old_par <- graphics::par(mfrow = c(n_row, n_col),
                           mar   = c(4, 4, 2.5, 1),
                           oma   = c(0, 0, 0, 0))
  on.exit(graphics::par(old_par), add = TRUE)

  palette_cols <- c("#1b7837", "#762a83", "#e08214", "#2166ac",
                    "#d6604d", "#4d9221", "#8073ac", "#bf812d")

  for (param in parameters) {
    df    <- cond_params[cond_params$parameter == param, ]
    title <- paste0(main_prefix, param)

    if (!is.null(group_var) && group_var %in% names(df)) {
      groups  <- sort(unique(df[[group_var]]))
      y_range <- range(df$value, na.rm = TRUE)
      x_range <- range(df[[x_var]], na.rm = TRUE)
      # Pad y-range slightly for legend
      y_pad   <- diff(y_range) * 0.12
      y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)

      graphics::plot(NA, xlim = x_range, ylim = y_range,
                     xlab = xlab, ylab = ylab, main = title, ...)
      for (g_idx in seq_along(groups)) {
        g_df <- df[df[[group_var]] == groups[g_idx], ]
        g_df <- g_df[order(g_df[[x_var]]), ]
        graphics::lines(g_df[[x_var]], g_df$value,
                        col = palette_cols[(g_idx - 1L) %% length(palette_cols) + 1L],
                        lwd = 1.8, lty = g_idx)
      }
      graphics::legend("topright",
                       legend  = paste0(group_var, " = ", groups),
                       col     = palette_cols[seq_along(groups)],
                       lty     = seq_along(groups),
                       lwd     = 1.5,
                       cex     = 0.78,
                       bty     = "n")
    } else {
      df <- df[order(df[[x_var]]), ]
      graphics::plot(df[[x_var]], df$value,
                     type = "l", lwd = 1.8,
                     xlab = xlab, ylab = ylab, main = title,
                     col  = palette_cols[1L], ...)
    }
  }
  invisible(NULL)
}


# =============================================================================
# Anchor-candidate summary and visualisation
# =============================================================================

#' Summarise item-screening fit indices for anchor selection
#'
#' Computes Δ-AIC, Δ-BIC, and Δ-SABIC for every item relative to the
#' corresponding facet's baseline CFA, and annotates each item with its DIF
#' type (none / intercept / loading / loading + intercept) based on the
#' significance results from the screening step.
#'
#' The returned data frame is intended to be inspected visually (via
#' [plot_anchor_candidates()]) before the user decides how many items to use
#' as anchors (`N_ANCHORS`) and which information criterion to rank by
#' (`ANCHOR_CRITERION`).
#'
#' @param item_fit_results Named list (level 1: facets; level 2: items) of
#'   `fit_df` data frames as returned by `run_mxsem()$fit_df`.
#' @param baseline_fits Named list (one element per facet) of `run_mxsem()`
#'   results for the baseline CFAs.  The `$fit_df` component is used.
#' @param sig_items Optional nested list (level 1: facets; level 2: items) of
#'   `significant_item_moderators()` results.  When supplied, each item is
#'   annotated with its DIF type.
#' @param criterion Character.  The information criterion used for sorting
#'   within facets.  One of `"SABIC"` (default), `"AIC"`, or `"BIC"`.
#'
#' @return A data frame with columns `facet`, `item`, `delta_AIC`,
#'   `delta_BIC`, `delta_SABIC`, `dif_type` (character: `"none"`,
#'   `"intercept"`, `"loading"`, or `"loading + intercept"`).  Rows are
#'   sorted descending by the chosen criterion within each facet (highest
#'   Δ-IC first = best anchor candidate).  The chosen `criterion` is stored
#'   as an attribute.
#'
#' @seealso [plot_anchor_candidates()], [moderate_loadings_and_intercepts()],
#'   [significant_item_moderators()]
#' @export
summarise_anchor_candidates <- function(item_fit_results,
                                        baseline_fits,
                                        sig_items  = NULL,
                                        criterion  = c("SABIC", "AIC", "BIC")) {
  criterion <- match.arg(criterion)

  rows <- list()
  for (facet_name in names(item_fit_results)) {
    base_fi    <- baseline_fits[[facet_name]]$fit_df
    base_aic   <- base_fi$AIC
    base_bic   <- base_fi$BIC
    base_sabic <- base_fi$SABIC

    for (item_name in names(item_fit_results[[facet_name]])) {
      fi <- item_fit_results[[facet_name]][[item_name]]
      if (is.null(fi)) next

      # Determine DIF annotation
      dif_type <- "none"
      if (!is.null(sig_items)) {
        tbl <- sig_items[[facet_name]][[item_name]]$table
        if (!is.null(tbl) && nrow(tbl) > 0) {
          has_load <- any(tbl$type == "loading"   & tbl$significant_corrected, na.rm = TRUE)
          has_int  <- any(tbl$type == "intercept" & tbl$significant_corrected, na.rm = TRUE)
          dif_type <- if (has_load && has_int) "loading + intercept"
                      else if (has_load)        "loading"
                      else if (has_int)         "intercept"
                      else                      "none"
        }
      }

      rows[[paste0(facet_name, ".", item_name)]] <- data.frame(
        facet       = facet_name,
        item        = item_name,
        delta_AIC   = round(fi$AIC   - base_aic,   2),
        delta_BIC   = round(fi$BIC   - base_bic,   2),
        delta_SABIC = round(fi$SABIC - base_sabic,  2),
        dif_type    = dif_type,
        stringsAsFactors = FALSE
      )
    }
  }

  df       <- do.call(rbind, rows)
  rownames(df) <- NULL
  crit_col <- paste0("delta_", criterion)
  # Sort descending within facet (highest = best anchor candidate)
  df <- do.call(rbind, lapply(split(df, df$facet), function(d)
    d[order(-d[[crit_col]]), ]))
  rownames(df) <- NULL
  attr(df, "criterion") <- criterion
  df
}


#' Plot anchor candidate information criteria after item screening
#'
#' Produces a dot/lollipop plot of Δ-AIC, Δ-BIC, or Δ-SABIC per item for
#' each facet.  Items are coloured by DIF type and the proposed anchor items
#' (top `n_anchors` clean candidates per facet) are highlighted.  A dashed
#' vertical line at 0 separates items that improve fit when freed (left) from
#' those that do not (right).
#'
#' The plot is printed immediately so it appears in RStudio / a graphics
#' device without an explicit `print()` call.
#'
#' @param anchor_df Data frame as returned by [summarise_anchor_candidates()].
#' @param criterion Character.  Information criterion to plot.  Defaults to
#'   the `"criterion"` attribute of `anchor_df`, or `"SABIC"` if absent.
#' @param n_anchors Integer.  Number of anchor items per facet to highlight.
#'   Should match the `N_ANCHORS` constant used in the pipeline.  Default 2.
#'
#' @return In ggplot2 mode: the `ggplot` object, invisibly.
#'   In base-R mode: `NULL` invisibly.  Plots are drawn as side-effects in
#'   both cases.
#'
#' @seealso [summarise_anchor_candidates()]
#' @export
plot_anchor_candidates <- function(anchor_df,
                                   criterion = NULL,
                                   n_anchors = 2L) {
  if (is.null(criterion)) criterion <- attr(anchor_df, "criterion")
  if (is.null(criterion)) criterion <- "SABIC"
  crit_col  <- paste0("delta_", criterion)
  n_anchors <- as.integer(n_anchors)

  # Mark proposed anchor items (top n_anchors clean candidates per facet)
  anchor_df <- do.call(rbind, lapply(split(anchor_df, anchor_df$facet), function(d) {
    clean     <- d[d$dif_type %in% c("none", "intercept"), ]
    pool      <- if (nrow(clean) >= n_anchors) clean else d
    top_items <- pool[order(-pool[[crit_col]]), "item"][seq_len(
                   min(n_anchors, nrow(pool)))]
    d$proposed_anchor <- d$item %in% top_items
    d
  }))
  rownames(anchor_df) <- NULL

  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  # ── ggplot2 path ─────────────────────────────────────────────────────────────
  if (use_ggplot) {
    dif_colours <- c("none"                = "#2166ac",
                     "intercept"           = "#f4a582",
                     "loading"             = "#d6604d",
                     "loading + intercept" = "#b2182b")

    p <- ggplot2::ggplot(
      anchor_df,
      ggplot2::aes(
        x      = .data[[crit_col]],
        y      = stats::reorder(.data[["item"]], .data[[crit_col]]),
        colour = .data[["dif_type"]],
        shape  = .data[["proposed_anchor"]]
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          colour = "grey50", linewidth = 0.5) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x    = 0,
          xend = .data[[crit_col]],
          yend = stats::reorder(.data[["item"]], .data[[crit_col]])
        ),
        colour = "grey75", linewidth = 0.4
      ) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::scale_colour_manual(
        values = dif_colours,
        name   = "DIF type"
      ) +
      ggplot2::scale_shape_manual(
        values = c("TRUE" = 16L, "FALSE" = 1L),
        name   = paste0("Proposed anchor\n(top ", n_anchors, ")")
      ) +
      ggplot2::facet_wrap(~ facet, scales = "free_y") +
      ggplot2::labs(
        title = paste0("Anchor candidate selection  \u2013  \u0394-", criterion,
                       "  (higher = better anchor)"),
        x     = paste0("\u0394-", criterion, " vs. baseline model"),
        y     = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "right",
                     strip.text      = ggplot2::element_text(face = "bold"))

    print(p)
    return(invisible(p))
  }

  # ── Base-R path ───────────────────────────────────────────────────────────────
  facets  <- unique(anchor_df$facet)
  n_col   <- length(facets)
  old_par <- graphics::par(mfrow = c(1L, n_col), mar = c(4, 7, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  dif_col_map <- c("none"                = "#2166ac",
                   "intercept"           = "#f4a582",
                   "loading"             = "#d6604d",
                   "loading + intercept" = "#b2182b")

  for (fn in facets) {
    d   <- anchor_df[anchor_df$facet == fn, ]
    d   <- d[order(d[[crit_col]]), ]
    col <- dif_col_map[d$dif_type]
    pch <- ifelse(d$proposed_anchor, 16L, 1L)
    graphics::dotchart(d[[crit_col]], labels = d$item,
                       main  = fn,
                       xlab  = paste0("\u0394-", criterion),
                       color = col, pch = pch, pt.cex = 1.4)
    graphics::abline(v = 0, lty = 2L, col = "grey50")
  }
  invisible(NULL)
}
