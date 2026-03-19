# =============================================================================
# tests/testthat/helper-mnlfa.R
# Shared test fixtures loaded automatically by testthat before any test file
# =============================================================================

# Zwei-Facetten-Domain: D_1 (4 Items) und D_2 (4 Items)
make_domain_models <- function() {
  list(
    D_1 = paste0("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4\n",
                 "D_1 ~~ 1*D_1\n",
                 "D_1 ~ 0*1"),
    D_2 = paste0("D_2 =~ d2_i1 + d2_i2 + d2_i3 + d2_i4\n",
                 "D_2 ~~ 1*D_2\n",
                 "D_2 ~ 0*1")
  )
}

# Minimales all_sig_mods data.frame für Tests
make_sig_mods <- function(facet, type, item = NA_character_,
                          moderator = "mod_z", significant_corrected = TRUE) {
  data.frame(
    facet                 = facet,
    item                  = item,
    type                  = type,
    moderator             = moderator,
    significant_corrected = significant_corrected,
    significant           = significant_corrected,
    param_name            = paste0("p_", facet, "_01"),
    Estimate              = 0.25,
    Std.Error             = 0.05,
    z                     = 5.0,
    p                     = 0.001,
    p_corrected           = 0.003,
    stringsAsFactors      = FALSE
  )
}

moderators_t <- c("mod_z", "sex")
