# =============================================================================
# R/domain_models.R
# Functions for assembling multi-facet domain models and the acquiescence factor
# =============================================================================

#' Build full domain models with all significant moderations
#'
#' For each domain (defined as the prefix before the first `_` in the facet
#' name), this function:
#' \enumerate{
#'   \item Calls [add_all_moderation_single_facet()] for each facet.
#'   \item Adds labelled covariance lines between all pairs of facets within
#'     the domain.
#'   \item Optionally appends extra syntax lines (e.g. acquiescence factor
#'     lines from [create_acquiescence_lines()]).
#' }
#'
#' The result is a complete domain model string ready for [run_mxsem()] with
#' `scale_loadings = FALSE` and `scale_latent_variances = FALSE`.  Scale
#' identification is achieved by fixing the first anchor item's loading to 1
#' (handled automatically by [add_all_moderation_single_facet()]).
#'
#' **Important**: if a facet has *no* significant moderations in `all_sig_mods`,
#' [add_all_moderation_single_facet()] returns the original `model_string`
#' unchanged (including its `LV ~~ 1*LV` baseline variance constraint).  This
#' situation arises only when all moderation tests are non-significant and is
#' typically not a problem in practice; however users should verify anchor
#' selection is still applied in such cases.
#'
#' @param facet_models Named list of per-facet baseline model strings.  Names
#'   must follow the convention `"[LETTER]_[DIGIT(S)]"` (e.g. `"E_1"`,
#'   `"P_3"`).  All facets with the same letter prefix (e.g. `P_1`, `P_2`,
#'   `P_3`) are assembled into a single domain model.
#' @param all_sig_mods Data frame from [combine_sig_moderators()].
#' @param anchor_items Character vector of anchor identifiers in the form
#'   `"facet_name.item_name"`.
#' @param moderators Optional character vector of all moderator names (used to
#'   compute correct numeric parameter indices).
#' @param extra_lines_per_domain Optional named list (one element per domain)
#'   of additional syntax lines appended after the facet models and covariance
#'   lines.  Typically the output of [create_acquiescence_lines()].
#'
#' @return A named list of character strings, one per domain.
#'
#' @seealso [add_all_moderation_single_facet()], [run_mxsem()],
#'   [create_acquiescence_lines()]
#' @export
build_full_domain_models <- function(facet_models, all_sig_mods, anchor_items,
                                     moderators = NULL,
                                     extra_lines_per_domain = NULL) {

  # ── Level-3 safeguard: warn if any anchor item has significant loading DIF ──
  # An anchor item with loading DIF compromises scale identification because the
  # latent variable's unit is defined by that item's (non-invariant) loading.
  # This situation should have been prevented by anchor selection (Level 2), but
  # may occur when screening power is low.
  if (!is.null(all_sig_mods) && nrow(all_sig_mods) > 0) {
    loading_dif_rows <- all_sig_mods[
      !is.na(all_sig_mods$type) &
      all_sig_mods$type == "loading" &
      !is.na(all_sig_mods$significant_corrected) &
      all_sig_mods$significant_corrected == TRUE, ]
    if (nrow(loading_dif_rows) > 0) {
      dif_keys        <- paste0(loading_dif_rows$facet, ".", loading_dif_rows$item)
      anchors_with_dif <- anchor_items[anchor_items %in% dif_keys]
      if (length(anchors_with_dif) > 0) {
        warning(
          "Scale identification warning: the following anchor item(s) have ",
          "significant loading DIF and may compromise identification:\n  ",
          paste(unique(anchors_with_dif), collapse = ", "), "\n",
          "Consider rerunning anchor selection after excluding these items ",
          "(see ?moderate_loadings_and_intercepts, section 'Identification').",
          call. = FALSE
        )
      }
    }
  }

  facet_names  <- names(facet_models)
  domain_names <- unique(sub("_.*", "", facet_names))
  domain_models <- list()

  for (domain in domain_names) {
    domain_facet_names <- grep(paste0("^", domain, "_"), facet_names, value = TRUE)

    domain_facet_models <- lapply(domain_facet_names, function(facet_name) {
      add_all_moderation_single_facet(
        model_string = facet_models[[facet_name]],
        facet_name   = facet_name,
        all_sig_mods = all_sig_mods,
        anchor_items = anchor_items,
        moderators   = moderators
      )
    })

    cov_lines <- c()
    if (length(domain_facet_names) > 1) {
      for (i in 1:(length(domain_facet_names) - 1)) {
        for (j in (i + 1):length(domain_facet_names)) {
          lv1       <- domain_facet_names[i]
          lv2       <- domain_facet_names[j]
          cov_label <- paste0("cov_", lv1, "_", lv2)
          cov_lines <- c(cov_lines, paste0(lv1, " ~~ ", cov_label, "*", lv2))
        }
      }
    }

    extra_lines <- extra_lines_per_domain[[domain]]

    domain_models[[domain]] <- paste(
      c(unlist(domain_facet_models), cov_lines, extra_lines),
      collapse = "\n\n"
    )
  }
  return(domain_models)
}


#' Generate acquiescence bifactor syntax lines for a domain
#'
#' Creates the syntax for an orthogonal acquiescence factor (`acq`) with fixed
#' loadings of `+1` (regular items) or `-1` (reversed items, identified by
#' the substring `_rec_` in the item name).  Optionally adds latent mean
#' regressions on the acquiescence factor.
#'
#' The returned list is intended to be passed as `extra_lines_per_domain` to
#' [build_full_domain_models()].
#'
#' **Reversed-item convention.**  An item is treated as reversed if and only if
#' its name contains the exact substring `_rec_` (case-sensitive).  This must
#' be present in both the data column name and the model syntax string.
#' Example: rename `e1_i8` → `e1_rec_i8` in your data frame and use the same
#' name in the measurement line.  Items without `_rec_` receive a `+1` loading
#' on the acquiescence factor.
#'
#' The acquiescence factor is constrained to be orthogonal to all content
#' factors (`acq ~~ 0 * LV`) and its mean is fixed to zero (`acq ~ 0*1`)
#' because acquiescence is assumed to be independent of and unconfounded with
#' the content factors.
#'
#' @param facet_models Named list of per-facet model strings (same as the
#'   first argument of [build_full_domain_models()]).  Item names are parsed
#'   from the `=~` measurement lines.
#' @param moderators Optional character vector.  If supplied, acquiescence
#'   mean regressions (`acq ~ m_acq_k * mod_k`) are added for each moderator.
#'
#' @return A named list (one element per domain) of character vectors with
#'   acquiescence factor syntax lines.
#'
#' @seealso [build_full_domain_models()]
#' @export
create_acquiescence_lines <- function(facet_models, moderators = NULL) {
  domain_names <- unique(sub("_.*", "", names(facet_models)))
  extra_lines  <- list()

  for (domain in domain_names) {
    facet_names <- grep(paste0("^", domain, "_"), names(facet_models), value = TRUE)

    all_items <- unlist(lapply(facet_names, function(facet) {
      lines         <- unlist(strsplit(facet_models[[facet]], "\n"))
      loading_lines <- grep("=~", lines, value = TRUE)
      item_parts    <- unlist(strsplit(loading_lines, "\\+"))
      item_names    <- trimws(gsub(".*=~", "", item_parts))
      item_names[item_names != ""]
    }))

    item_clean <- gsub("^[^*]*\\*", "", all_items)
    weights    <- ifelse(grepl("_rec_", item_clean), "-1", "1")
    loadings   <- paste0(weights, "*", item_clean)
    acq_line   <- paste0("acq=~", paste(loadings, collapse = " + "))
    orthogonal <- paste0("acq~~0*", facet_names)
    mean_line  <- "acq ~ 0*1"
    mod_lines  <- if (!is.null(moderators))
      paste0("acq ~ m_acq_", seq_along(moderators), "*", moderators)
    else
      character(0)

    extra_lines[[domain]] <- c(acq_line, orthogonal, mean_line, mod_lines)
  }
  return(extra_lines)
}
