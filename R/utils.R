# =============================================================================
# R/utils.R
# Internal helper functions (not exported)
# =============================================================================

# Suppress R CMD check "no visible binding for global variable" notes that
# arise from subset() and similar NSE-based functions using bare column names.
# These are all data-frame column names, not actual global variables.
utils::globalVariables(c(
  # all_sig_mods / param_df column names used with subset()
  "facet", "significant_corrected", "type", "item",
  "name", "matrix", "moderator",
  # param enrichment columns
  "Estimate", "Std.Error", "z", "p", "p_corrected",
  "significant", "param_name", "domain", "source",
  "row", "col", "model"
))

#' Extract moderator names from a run_mxsem result object
#'
#' Checks first on `param_df`, then falls back to `extracted` for backward
#' compatibility with result objects created before the attribute was added to
#' `param_df`.
#'
#' @param mxsem_result A list returned by [run_mxsem()].
#' @return A character vector of moderator names, or `NULL`.
#' @keywords internal
.get_moderators <- function(mxsem_result) {
  mods <- attr(mxsem_result$param_df,    "moderators")
  if (is.null(mods)) mods <- attr(mxsem_result$extracted, "moderators")
  mods
}

#' Normalise a model string to a character vector of lines
#'
#' Several functions receive `model_string` either as a single `\n`-delimited
#' string or already as a character vector of individual lines.  This helper
#' ensures a consistent line-vector representation in both cases.
#'
#' @param model_string Character vector (length 1 or > 1) representing a
#'   lavaan-style model syntax.
#' @return A character vector with one element per non-empty line.
#' @keywords internal
.normalise_model_string <- function(model_string) {
  if (length(model_string) == 1 && grepl("\n", model_string)) {
    model_string <- trimws(unlist(strsplit(model_string, "\n")))
    model_string <- model_string[nchar(model_string) > 0]
  }
  model_string
}
