# =============================================================================
# tests/testthat/test-integration.R
# Teil 3: Integrations-Tests mit synthetischen Daten
# Erfordert: mxsem, OpenMx – werden mit skip_if_not_installed() geprüft
# =============================================================================

# Synthetischer Datensatz erzeugen (einmal für alle Integration-Tests)
local({
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")

  set.seed(2025)
  N     <- 500
  mod_z <- as.numeric(scale(rnorm(N)))
  sex   <- rbinom(N, 1, 0.5)

  eta1  <- 0.30 * mod_z + rnorm(N, sd = sqrt(1 - 0.09))
  eta2  <- 0.25 * sex   + rnorm(N, sd = sqrt(1 - 0.0625))

  gen_items <- function(eta, item_effects = NULL, n = 4) {
    sapply(seq_len(n), function(i) {
      extra <- if (!is.null(item_effects[[i]])) item_effects[[i]] else 0
      round(0.7 * eta + extra + rnorm(N, sd = sqrt(0.51)), 3)
    })
  }

  items_d1 <- gen_items(eta1,
    item_effects = list(NULL, NULL, 0.20 * mod_z, NULL))
  items_d2 <- gen_items(eta2)

  colnames(items_d1) <- paste0("d1_i", 1:4)
  colnames(items_d2) <- paste0("d2_i", 1:4)

  synth <<- data.frame(items_d1, items_d2, mod_z, sex)

  facet_models_synth <<- list(
    D_1 = "D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4\nD_1 ~~ 1*D_1\nD_1 ~ 0*1",
    D_2 = "D_2 =~ d2_i1 + d2_i2 + d2_i3 + d2_i4\nD_2 ~~ 1*D_2\nD_2 ~ 0*1"
  )
})

# ---------------------------------------------------------------------------
# Test 3.1: Basis-CFA
# ---------------------------------------------------------------------------

test_that("[Integration] Basis-CFA laeuft fuer beide Facetten", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  fit_d1 <- run_mxsem(facet_models_synth$D_1, data = synth,
                      scale_loadings = FALSE, scale_latent_variances = TRUE,
                      add_ref_models = FALSE, moderators = moderators_t, seed = 42)
  fit_d2 <- run_mxsem(facet_models_synth$D_2, data = synth,
                      scale_loadings = FALSE, scale_latent_variances = TRUE,
                      add_ref_models = FALSE, moderators = moderators_t, seed = 42)

  for (fit in list(fit_d1, fit_d2)) {
    expect_named(fit, c("model_string", "script", "fit", "ref_models",
                        "summary", "extracted", "fit_df", "param_df"))
    expect_true(nrow(fit$param_df) > 0, label = "param_df leer")
    expect_false(is.null(attr(fit$param_df, "moderators")),
                 label = "Moderatoren-Attribut fehlt auf param_df")
  }
})

# ---------------------------------------------------------------------------
# Test 3.2: Screening – Latente Mittelwert-Moderation
# ---------------------------------------------------------------------------

test_that("[Integration] significant_mean_moderators erkennt mod_z-Effekt auf D_1", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  model_lm <- add_latent_mean_regression(facet_models_synth$D_1, moderators_t)
  fit <- run_mxsem(model_lm, data = synth, scale_loadings = FALSE,
                   scale_latent_variances = TRUE, add_ref_models = FALSE,
                   moderators = moderators_t, seed = 42)

  sig <- significant_mean_moderators(fit, alpha = 0.05, correction_method = "none")
  expect_s3_class(sig$table, "data.frame")

  expect_true(
    any(sig$table$significant & sig$table$moderator == "mod_z"),
    label = "mod_z sollte signifikant auf D_1 latenten Mittelwert wirken"
  )

  sex_row <- sig$table[sig$table$moderator == "sex", ]
  if (nrow(sex_row) > 0) {
    expect_false(sex_row$significant[1],
                 label = "sex sollte nicht signifikant auf D_1 Mittelwert wirken")
  }
})

# ---------------------------------------------------------------------------
# Test 3.3: Screening – Item-Level (anchor selection)
# ---------------------------------------------------------------------------

test_that("[Integration] Anker-Kandidat hat kleineres delta SABIC als moderiertes Item", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  item_models <- moderate_loadings_and_intercepts(facet_models_synth$D_1, moderators_t)

  fit_anchor <- run_mxsem(item_models$d1_i1, data = synth,
                          scale_loadings = FALSE, add_ref_models = FALSE,
                          moderators = moderators_t, seed = 42)
  fit_mod    <- run_mxsem(item_models$d1_i3, data = synth,
                          scale_loadings = FALSE, add_ref_models = FALSE,
                          moderators = moderators_t, seed = 42)
  fit_base   <- run_mxsem(facet_models_synth$D_1, data = synth,
                          scale_loadings = FALSE, scale_latent_variances = TRUE,
                          add_ref_models = FALSE, moderators = moderators_t, seed = 42)

  delta_anchor <- fit_anchor$fit_df$SABIC - fit_base$fit_df$SABIC
  delta_mod    <- fit_mod$fit_df$SABIC    - fit_base$fit_df$SABIC

  expect_true(delta_anchor > delta_mod,
              label = paste0(
                "d1_i1 (Anker-Kandidat) sollte kleineres delta SABIC haben als d1_i3. ",
                "delta_anchor=", round(delta_anchor, 1),
                " delta_mod=",   round(delta_mod, 1)))
})

# ---------------------------------------------------------------------------
# Test 3.4: Domain-Modell – Multi-Faktor, mit Anker-Items
# ---------------------------------------------------------------------------

test_that("[Integration] Domain-Modell mit 2 Facetten laeuft durch", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  all_sig_mods_synth <- rbind(
    make_sig_mods("D_1", "mean",      moderator = "mod_z"),
    make_sig_mods("D_2", "mean",      moderator = "sex"),
    make_sig_mods("D_1", "intercept", item = "d1_i3", moderator = "mod_z")
  )
  all_sig_mods_synth$significant_corrected <- TRUE
  anchor_items_synth <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  domain_models_synth <- build_full_domain_models(
    facet_models = facet_models_synth,
    all_sig_mods = all_sig_mods_synth,
    anchor_items = anchor_items_synth,
    moderators   = moderators_t
  )

  lines <- unlist(strsplit(domain_models_synth$D, "\n"))
  expect_true(any(grepl("D_1 ~~ cov_D_1_D_2\\*D_2", lines)),
              label = "Kovarianz zwischen D_1 und D_2 fehlt")

  fit_domain <- run_mxsem(
    domain_models_synth$D, data = synth,
    scale_loadings = FALSE, scale_latent_variances = FALSE,
    add_ref_models = FALSE, moderators = moderators_t,
    seed = 42, extraTries = 20
  )

  expect_s3_class(fit_domain$param_df, "data.frame")
  expect_true(nrow(fit_domain$param_df) > 0)

  mean_mods <- fit_domain$param_df[
    fit_domain$param_df$type == "latent_mean" &
    !is.na(fit_domain$param_df$col) &
    fit_domain$param_df$col == "mod_z", ]
  expect_true(nrow(mean_mods) > 0,
              label = "Latente Mittelwert-Moderation durch mod_z fehlt in param_df")
})

# ---------------------------------------------------------------------------
# Test 3.5: extract_all_significant_moderators
# ---------------------------------------------------------------------------

test_that("[Integration] extract_all_significant_moderators liefert strukturiertes Ergebnis", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  all_sig_mods_synth <- rbind(
    make_sig_mods("D_1", "mean", moderator = "mod_z"),
    make_sig_mods("D_2", "mean", moderator = "sex")
  )
  all_sig_mods_synth$significant_corrected <- TRUE
  anchor_items_synth <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  domain_models_synth <- build_full_domain_models(
    facet_models = facet_models_synth,
    all_sig_mods = all_sig_mods_synth,
    anchor_items = anchor_items_synth,
    moderators   = moderators_t
  )

  fit_d <- run_mxsem(domain_models_synth$D, data = synth,
                     scale_loadings = FALSE, scale_latent_variances = FALSE,
                     add_ref_models = FALSE, moderators = moderators_t,
                     seed = 42, extraTries = 20)

  dom_fits <- list(D = list(D = fit_d))
  suppressMessages({
    all_mods <- extract_all_significant_moderators(dom_fits, alpha = 0.10,
                                                    correction_method = "none")
  })

  expect_type(all_mods, "list")
  expect_true("D" %in% names(all_mods), label = "Domain D fehlt im Ergebnis")
  expect_s3_class(all_mods$D, "data.frame")
  expect_true("type"      %in% names(all_mods$D))
  expect_true("moderator" %in% names(all_mods$D))
})

# ---------------------------------------------------------------------------
# Test 3.6: build_all_individual_parameters
# ---------------------------------------------------------------------------

test_that("[Integration] build_all_individual_parameters erzeugt individuelle Vorhersagen", {
  skip_if_not_installed("mxsem")
  skip_if_not_installed("OpenMx")
  skip_on_cran()

  N <- nrow(synth)

  all_sig_mods_synth <- rbind(
    make_sig_mods("D_1", "mean", moderator = "mod_z"),
    make_sig_mods("D_2", "mean", moderator = "sex")
  )
  all_sig_mods_synth$significant_corrected <- TRUE
  anchor_items_synth <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  dom_model <- build_full_domain_models(
    facet_models = facet_models_synth,
    all_sig_mods = all_sig_mods_synth,
    anchor_items = anchor_items_synth,
    moderators   = moderators_t
  )

  fit_d <- run_mxsem(dom_model$D, data = synth, scale_loadings = FALSE,
                     scale_latent_variances = FALSE, add_ref_models = FALSE,
                     moderators = moderators_t, seed = 42, extraTries = 20)

  suppressMessages({
    ind_params <- build_all_individual_parameters(
      param_df           = fit_d$param_df,
      data               = synth,
      p_adjust_method    = "none",
      significance_level = 0.10,
      verbose            = FALSE
    )
  })

  expect_named(ind_params,
               c("latent_means", "latent_variances", "loadings", "intercepts"))

  lm <- ind_params$latent_means
  expect_false(is.null(lm), label = "latent_means ist NULL")
  expect_equal(nrow(lm), N,
               label = paste("latent_means hat", nrow(lm), "statt", N, "Zeilen"))
  expect_true(any(grepl("D_1_pred", names(lm))),
              label = "Spalte D_1_pred fehlt in latent_means")
})
