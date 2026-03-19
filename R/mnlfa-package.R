#' mnlfa: Moderated Nonlinear Factor Analysis
#'
#' @description
#' The **mnlfa** package implements the stepwise Moderated Nonlinear Factor
#' Analysis (MNLFA) procedure for assessing differential item functioning (DIF)
#' and measurement invariance in the presence of one or more observed moderator
#' variables.
#'
#' Unlike existing R implementations, **mnlfa** fully supports
#' *multidimensional* measurement models with multiple correlated latent
#' variables (facets) per domain — the typical case when working with
#' hierarchical personality or ability constructs.
#'
#' @section Main workflow:
#' A typical MNLFA analysis proceeds in six steps:
#' \enumerate{
#'   \item **Baseline CFAs**: fit one CFA per facet with
#'     [run_mxsem()] (`scale_latent_variances = TRUE`).
#'   \item **Screen latent mean moderation**: call
#'     [add_latent_mean_regression()] + [run_mxsem()] +
#'     [significant_mean_moderators()] for each facet.
#'   \item **Screen item-level moderation**: call
#'     [moderate_loadings_and_intercepts()] + [run_mxsem()] +
#'     [significant_item_moderators()] for each item within each facet.
#'   \item **Select anchor items**: rank items by Δ-SABIC; items with the
#'     smallest improvement are best anchors.  At least two per facet are
#'     required.
#'   \item **Build full domain model**: assemble the moderated syntax with
#'     [build_full_domain_models()] and re-fit with [run_mxsem()]
#'     (`scale_loadings = FALSE, scale_latent_variances = FALSE`).
#'     Optionally add covariance moderation ([add_covariance_moderation()])
#'     and an acquiescence factor ([create_acquiescence_lines()]).
#'   \item **Compute person-specific parameters**: use
#'     [build_all_individual_parameters()] to obtain individual factor loadings,
#'     intercepts, latent means, and latent variances.
#' }
#'
#' @section Anchor item format:
#' Anchor items are identified by the composite key `"facet_name.item_name"`,
#' e.g. `"E_1.e1_i1"`.  This convention arises naturally from
#' `unlist(facet.items, recursive = FALSE)` in the pipeline.  The only
#' exception is [build_moderated_sem_one_latent()], which uses bare item names.
#'
#' @section Data requirements and naming conventions:
#'
#' **Data format.**  The data must be a `data.frame` in wide format with one
#' row per person.  All item columns and all moderator columns must be present.
#' Moderators should be numeric; binary variables (e.g. sex coded 0/1) are
#' supported.  Centring continuous moderators (e.g. z-standardisation) is
#' strongly recommended so that the log-variance intercept `v0_*` and the
#' latent-mean intercept `m_*` have a meaningful reference point.
#'
#' **Facet names.**  Latent variable names must follow the pattern
#' `[LETTER]_[DIGIT(S)]`, e.g. `E_1`, `P_3`, `N_12`.  Several internal
#' regex patterns rely on this convention to extract domain names and parameter
#' labels.  The domain name is the part before the first underscore
#' (`sub("_.*", "", facet_name)`), so `P_1`, `P_2`, and `P_3` all belong to
#' domain `P` and will be assembled into a single domain model by
#' [build_full_domain_models()].
#'
#' **Item names.**  Item names must not contain dots (`.`), because the package
#' uses `.` as a separator in the composite anchor key `"facet.item"`.
#' Underscores are fine.  Item names must match exactly between the model
#' syntax string and the column names in the data frame.
#'
#' **Reversed items and the acquiescence factor.**  [create_acquiescence_lines()]
#' assigns a loading of `+1` to regular items and `-1` to reversed items.
#' An item is treated as reversed if and only if its name contains the
#' substring `_rec_` (case-sensitive), e.g. `e1_rec_i8`.  This substring must
#' appear in both the data column name and the model syntax string; the
#' simulation script demonstrates the required column renaming step:
#' \preformatted{
#'   names(data)[names(data) == "e1_i8"] <- "e1_rec_i8"
#' }
#'
#' **Fit indices.**  MNLFA models rely on definition variables (person-specific
#' parameters), which makes the standard reference model approach for computing
#' CFI, TLI, and RMSEA inapplicable.  All calls to [run_mxsem()] in the
#' pipeline should therefore use `add_ref_models = FALSE`.  The sample-size
#' adjusted BIC (SABIC) is used instead for model comparison (Δ-SABIC for
#' anchor selection; SABIC for variance and mean moderation screening).
#'
#' @section References:
#' Bauer, D. J. (2017). A more general model for testing measurement
#' invariance and differential item functioning.
#' *Psychological Methods*, *22*(3), 507–526.
#' \doi{10.1037/met0000077}
#'
#' Gottfredson, N. C., Cole, V. T., Giordano, M. L., Grimm, K. J.,
#' Bauer, D. J., & Hussong, A. M. (2019). Simplifying the implementation of
#' modern scale scoring methods with an automated R package: Automated Moderated
#' Nonlinear Factor Analysis (AutoMNLFA).
#' *Psychological Methods*, *24*(5), 669–682.
#' \doi{10.1037/met0000193}
#'
#' Kolbe, L., Molenaar, D., Tučková, T., Kořínek, D., & Mayer, A. (2022).
#' Testing measurement invariance with many groups using MNLFA.
#' *Structural Equation Modeling*, *29*(3), 384–401.
#' \doi{10.1080/10705511.2021.1984772}
#'
#' @keywords internal
#'
#' @importFrom OpenMx mxTryHard mxRefModels
#' @importFrom mxsem mxsem set_starting_values
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split str_trim str_detect
#' @importFrom stats na.omit p.adjust pnorm setNames
#' @importFrom graphics par plot lines legend
"_PACKAGE"
