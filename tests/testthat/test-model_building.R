# =============================================================================
# tests/testthat/test-model_building.R
# Teil 1: String-Manipulations-Tests (kein OpenMx erforderlich)
# Teil 2: Anker-Item-Logik (kein OpenMx erforderlich)
# =============================================================================

# ============================================================================
# TEIL 1: STRING-MANIPULATIONS-TESTS
# ============================================================================

# ---------------------------------------------------------------------------
# Test 1.1: add_latent_mean_regression
# ---------------------------------------------------------------------------

test_that("add_latent_mean_regression fuegt Intercept + Regressions-Zeilen hinzu", {
  model  <- make_domain_models()$D_1
  result <- add_latent_mean_regression(model, moderators_t)
  lines  <- unlist(strsplit(result, "\n"))

  expect_true(any(grepl("^D_1 ~ m_D_1\\*1$",     lines)), label = "Intercept-Zeile fehlt")
  expect_true(any(grepl("m_D_1_01.*mod_z",        lines)), label = "mod_z fehlt")
  expect_true(any(grepl("m_D_1_02.*sex",          lines)), label = "sex fehlt")
  expect_true(any(grepl("D_1 =~",                 lines)), label = "Messgleichung fehlt")
})

test_that("add_latent_mean_regression ist idempotent", {
  model   <- make_domain_models()$D_1
  result1 <- add_latent_mean_regression(model, moderators_t)
  result2 <- add_latent_mean_regression(result1, moderators_t)
  lines   <- unlist(strsplit(result2, "\n"))

  expect_equal(sum(grepl("^D_1 ~ m_D_1\\*1$", lines)), 1,
               label = "Intercept-Zeile nach zweitem Aufruf doppelt vorhanden")
})

# ---------------------------------------------------------------------------
# Test 1.2: add_latent_variance_moderation
# ---------------------------------------------------------------------------

test_that("add_latent_variance_moderation erzeugt exp()-Formel mit Intercept", {
  model  <- make_domain_models()$D_1
  result <- add_latent_variance_moderation(model, moderators_t)
  lines  <- unlist(strsplit(result, "\n"))

  formula_line <- grep(":= exp\\(", lines, value = TRUE)
  expect_length(formula_line, 1)
  expect_true(grepl("v0_D_1",       formula_line), label = "Intercept v0_ fehlt")
  expect_true(grepl("v_D_1_01.*mod_z", formula_line), label = "mod_z fehlt in Varianz")
  expect_true(any(grepl("D_1 ~~ var_D_1\\*D_1", lines)), label = "Varianz-Link fehlt")
  expect_true(any(grepl("^!v0_D_1;",            lines)), label = "!v0_ Deklaration fehlt")
})

# ---------------------------------------------------------------------------
# Test 1.3: moderate_loadings_and_intercepts
# ---------------------------------------------------------------------------

test_that("moderate_loadings_and_intercepts gibt Liste mit einem Modell pro Item zurueck", {
  model   <- make_domain_models()$D_1
  results <- moderate_loadings_and_intercepts(model, moderators_t)

  expect_type(results, "list")
  expect_length(results, 4)
  expect_setequal(names(results), c("d1_i1", "d1_i2", "d1_i3", "d1_i4"))
})

test_that("moderate_loadings_and_intercepts: Ziel-Item hat Load-Label, andere nicht", {
  model   <- make_domain_models()$D_1
  results <- moderate_loadings_and_intercepts(model, moderators_t)

  lines_i1  <- unlist(strsplit(results$d1_i1, "\n"))
  load_line <- grep("D_1 =~", lines_i1, value = TRUE)

  expect_true(grepl("D_1_load_01\\*d1_i1", load_line),
              label = "d1_i1 hat kein Load-Label")
  expect_false(grepl("D_1_load_02\\*d1_i2", load_line),
               label = "d1_i2 hat faelschlich ein Load-Label")
})

test_that("moderate_loadings_and_intercepts: Varianz-Fixierung fuer Identifikation", {
  model   <- make_domain_models()$D_1
  results <- moderate_loadings_and_intercepts(model, moderators_t)
  lines   <- unlist(strsplit(results$d1_i1, "\n"))

  expect_true(any(grepl("D_1 ~~ 1\\*D_1", lines)),
              label = "Latente Varianz sollte auf 1 fixiert sein (Screening-Schritt)")
})

# ---------------------------------------------------------------------------
# Test 1.4: add_covariance_moderation – Zwei-Faktor-Modell
# ---------------------------------------------------------------------------

test_that("add_covariance_moderation fuegt rho und Varianz-Links korrekt ein", {
  model  <- paste0("D_1 =~ d1_i1 + d1_i2 + d1_i3\n",
                   "D_2 =~ d2_i1 + d2_i2 + d2_i3\n",
                   "D_1 ~~ cov_D_1_D_2*D_2")
  result <- add_covariance_moderation(model, moderators_t)
  lines  <- unlist(strsplit(result, "\n"))

  expect_true(any(grepl("rho_cov_D_1_D_2 :=", lines)), label = "rho-Definition fehlt")
  expect_true(any(grepl("cov_D_1_D_2 := sqrt",lines)), label = "sqrt-Formel fehlt")
  expect_true(any(grepl("D_1 ~~ var_D_1\\*D_1",lines)), label = "Varianz-Link D_1 fehlt")
  expect_true(any(grepl("D_2 ~~ var_D_2\\*D_2",lines)), label = "Varianz-Link D_2 fehlt")
  expect_equal(sum(grepl("D_1 ~~ var_D_1\\*D_1", lines)), 1,
               label = "Varianz-Link D_1 doppelt vorhanden")
})


# ============================================================================
# TEIL 2: ANKER-ITEM-LOGIK
# ============================================================================

# ---------------------------------------------------------------------------
# Test 2.1: Format der Anchor Items ("facet.item" Konvention)
# ---------------------------------------------------------------------------

test_that("add_all_moderation_single_facet: Anker im Format 'facet.item' wird erkannt", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4",
                    "D_1 ~~ 1*D_1", "D_1 ~ 0*1")
  all_sig_mods <- make_sig_mods("D_1", "mean")
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2")

  result    <- add_all_moderation_single_facet(model_lines, "D_1", all_sig_mods,
                                               anchor_items, moderators_t)
  meas_line <- result[grepl("D_1 =~", result)]

  expect_true(grepl("1\\*d1_i1", meas_line),
              label = "Erstes Anker-Item sollte Ladung 1 haben")
  expect_true(grepl("d1_i2", meas_line),
              label = "Zweites Anker-Item fehlt in Messgleichung")
  expect_false(grepl("D_1_load_.*\\*d1_i2", meas_line),
               label = "Zweites Anker-Item darf kein Moderations-Label haben")
})

test_that("add_all_moderation_single_facet: Anker ohne facet-Praefix fuehrt zu Fehler", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4")
  all_sig_mods <- make_sig_mods("D_1", "mean")
  wrong_format <- c("d1_i1", "d1_i2")

  expect_error(
    add_all_moderation_single_facet(model_lines, "D_1", all_sig_mods,
                                    wrong_format, moderators_t),
    regexp = "No anchor item found",
    label  = "Fehler bei falschem Anker-Format erwartet"
  )
})

# ---------------------------------------------------------------------------
# Test 2.2: Anker-Items werden korrekt aus Moderation ausgeschlossen
# ---------------------------------------------------------------------------

test_that("Anker-Items erhalten keine Moderations-Labels fuer Ladungen", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4",
                    "D_1 ~~ 1*D_1", "D_1 ~ 0*1")
  all_sig_mods <- make_sig_mods("D_1", "loading", item = "d1_i2")
  all_sig_mods$significant_corrected <- TRUE
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2")

  result    <- add_all_moderation_single_facet(model_lines, "D_1", all_sig_mods,
                                               anchor_items, moderators_t)
  meas_line <- result[grepl("D_1 =~", result)]

  expect_false(grepl("D_1_load_.*\\*d1_i2", meas_line),
               label = "Anker-Item d1_i2 hat faelschlich ein Moderations-Label")
})

test_that("Anker-Items erhalten keine Intercept-Moderation", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4",
                    "D_1 ~~ 1*D_1", "D_1 ~ 0*1")
  all_sig_mods <- make_sig_mods("D_1", "intercept", item = "d1_i1")
  all_sig_mods$significant_corrected <- TRUE
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2")

  result <- add_all_moderation_single_facet(model_lines, "D_1", all_sig_mods,
                                             anchor_items, moderators_t)

  expect_false(any(grepl("int_D_1_d1_i1", result)),
               label = "Anker-Item d1_i1 hat faelschlich Intercept-Moderation")
})

test_that("Nicht-Anker-Items koennen Moderations-Labels erhalten", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4",
                    "D_1 ~~ 1*D_1", "D_1 ~ 0*1")
  all_sig_mods <- make_sig_mods("D_1", "loading", item = "d1_i3")
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2")

  result    <- add_all_moderation_single_facet(model_lines, "D_1", all_sig_mods,
                                               anchor_items, moderators_t)
  meas_line <- result[grepl("D_1 =~", result)]

  expect_true(grepl("D_1_load_03\\*d1_i3", meas_line),
              label = "Nicht-Anker d1_i3 sollte Moderations-Label haben")
})

# ---------------------------------------------------------------------------
# Test 2.3: Identifikationsstrategie
# ---------------------------------------------------------------------------

test_that("Screening-Modelle haben Varianz = 1", {
  model <- make_domain_models()$D_1
  expect_true(grepl("D_1 ~~ 1\\*D_1", model),
              label = "Screening-Modell sollte latente Varianz = 1 haben")
})

test_that("Domain-Modell nach build_full_domain_models hat Ladung = 1 (nicht Varianz)", {
  facet_models <- make_domain_models()
  all_sig_mods <- make_sig_mods("D_1", "mean")
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  domain_model <- build_full_domain_models(
    facet_models = facet_models,
    all_sig_mods = all_sig_mods,
    anchor_items = anchor_items,
    moderators   = moderators_t
  )

  lines <- unlist(strsplit(domain_model$D, "\n"))

  expect_true(any(grepl("D_1 =~ 1\\*d1_i1", lines)),
              label = "Anker-Ladung 1*d1_i1 fehlt im Domain-Modell")
  expect_false(any(grepl("D_1 ~~ 1\\*D_1", lines)),
               label = "Domain-Modell darf Varianz D_1 nicht mehr auf 1 fixieren")
})

# ---------------------------------------------------------------------------
# Test 2.4: Multi-Faktor-Domain – Kovarianz-Struktur
# ---------------------------------------------------------------------------

test_that("build_full_domain_models fuegt Kovarianz zwischen allen Facetten-Paaren ein", {
  facet_models <- make_domain_models()
  all_sig_mods <- make_sig_mods("D_1", "mean")
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  domain_model <- build_full_domain_models(
    facet_models = facet_models,
    all_sig_mods = all_sig_mods,
    anchor_items = anchor_items,
    moderators   = moderators_t
  )

  lines <- unlist(strsplit(domain_model$D, "\n"))

  expect_true(any(grepl("D_1 ~~ cov_D_1_D_2\\*D_2", lines)),
              label = "Kovarianz zwischen D_1 und D_2 fehlt")
})

test_that("build_full_domain_models: extra_lines_per_domain werden angehaengt", {
  facet_models <- make_domain_models()
  all_sig_mods <- make_sig_mods("D_1", "mean")
  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")
  extra        <- list(D = c("d1_i1 ~~ d1_i2"))

  domain_model <- build_full_domain_models(
    facet_models           = facet_models,
    all_sig_mods           = all_sig_mods,
    anchor_items           = anchor_items,
    moderators             = moderators_t,
    extra_lines_per_domain = extra
  )

  lines <- unlist(strsplit(domain_model$D, "\n"))
  expect_true(any(grepl("d1_i1 ~~ d1_i2", lines)),
              label = "extra_lines wurden nicht angehaengt")
})

test_that("build_full_domain_models: Drei-Facetten-Domain hat 3 Kovarianz-Paare", {
  facet_models_3 <- list(
    E_1 = "E_1 =~ e1_i1 + e1_i2 + e1_i3\nE_1 ~~ 1*E_1\nE_1 ~ 0*1",
    E_2 = "E_2 =~ e2_i1 + e2_i2 + e2_i3\nE_2 ~~ 1*E_2\nE_2 ~ 0*1",
    E_3 = "E_3 =~ e3_i1 + e3_i2 + e3_i3\nE_3 ~~ 1*E_3\nE_3 ~ 0*1"
  )
  all_sig_mods <- make_sig_mods("E_1", "mean")
  anchor_items <- c("E_1.e1_i1", "E_1.e1_i2",
                    "E_2.e2_i1", "E_2.e2_i2",
                    "E_3.e3_i1", "E_3.e3_i2")

  domain_model <- build_full_domain_models(
    facet_models = facet_models_3,
    all_sig_mods = all_sig_mods,
    anchor_items = anchor_items,
    moderators   = moderators_t
  )

  lines <- unlist(strsplit(domain_model$E, "\n"))
  n_cov <- sum(grepl("~~ cov_", lines))
  expect_equal(n_cov, 3,
               label = "Drei-Facetten-Domain sollte genau 3 Kovarianz-Zeilen haben")
})

# ---------------------------------------------------------------------------
# Test 2.5: Format-Inkonsistenz zwischen den beiden Hilfsfunktionen (dokumentiert)
# ---------------------------------------------------------------------------

test_that("[DOKUMENTIERT] build_moderated_sem_one_latent nutzt anderes Anker-Format", {
  model_lines  <- c("D_1 =~ d1_i1 + d1_i2 + d1_i3 + d1_i4")
  all_sig_mods <- make_sig_mods("D_1", "mean")

  result_ok <- build_moderated_sem_one_latent(
    model_string = model_lines,
    facet_name   = "D_1",
    all_sig_mods = all_sig_mods,
    anchor_items = c("d1_i1", "d1_i2"),
    moderators   = moderators_t
  )
  expect_true(any(grepl("1\\*d1_i1", result_ok)),
              label = "Anker-Ladung mit reinem item-Format nicht erkannt")

  skip("Inkonsistenz ist bekannt und dokumentiert – kein Fehler, sondern Design-Entscheidung")
})

# ---------------------------------------------------------------------------
# Test 2.6: Moderationen in mehreren Facetten
# ---------------------------------------------------------------------------

test_that("build_full_domain_models verarbeitet Moderationen in mehreren Facetten", {
  facet_models <- make_domain_models()

  all_sig_mods <- rbind(
    make_sig_mods("D_1", "mean",      moderator = "mod_z"),
    make_sig_mods("D_1", "intercept", item = "d1_i3", moderator = "sex"),
    make_sig_mods("D_2", "mean",      moderator = "sex")
  )
  all_sig_mods$significant_corrected <- TRUE

  anchor_items <- c("D_1.d1_i1", "D_1.d1_i2", "D_2.d2_i1", "D_2.d2_i2")

  domain_model <- build_full_domain_models(
    facet_models = facet_models,
    all_sig_mods = all_sig_mods,
    anchor_items = anchor_items,
    moderators   = moderators_t
  )

  lines <- unlist(strsplit(domain_model$D, "\n"))

  expect_true(any(grepl("D_1 ~ m_D_1_01\\*mod_z", lines)),
              label = "D_1 Mittelwert-Moderation durch mod_z fehlt")
  expect_true(any(grepl("D_2 ~ m_D_2_02\\*sex", lines)),
              label = "D_2 Mittelwert-Moderation durch sex fehlt")
  expect_true(any(grepl("int_D_1_d1_i3", lines)),
              label = "Intercept-Moderation fuer d1_i3 fehlt")
})
