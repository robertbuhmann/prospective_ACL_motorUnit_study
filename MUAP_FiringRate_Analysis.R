# MUAP_FiringRate_Analysis.R
# Secondary exploratory motor-unit action potential (MUAP) and firing-rate analysis
# Prospective ACL motor unit study
#
# This script is the R implementation intended for the public reproducibility
# repository. It uses the de-identified MUAP_FiringRate_Data.csv file.
#
# IMPORTANT:
# 1. Surface-detected MUAP amplitude is treated as an electrophysiological
#    characteristic, not a direct anatomical measure of motor-unit size.
# 2. MUAP amplitude is log-transformed and standardized within each contraction.
# 3. The same original Delsys MU can occur in multiple holds from one recording.
#    Crossed random intercepts therefore account for participant, recording-MU,
#    and contraction bout.
# 4. Recruitment threshold is not reconstructed in this analysis.
#
# The analysis is intentionally aligned with the final Stage 2F v3 analysis:
#   Primary outcome: plateau_train_rate_hz
#   Predictor: within-bout z(log MUAP amplitude)
#   Covariates: intensity, muscle, limb, timepoint
#   Random intercepts: participant, original recording-MU, contraction bout
#   Secondary: MUAP x limb
#   Sensitivities: >=5 firings; A+B; A-clean; mean instantaneous firing rate
#
# The original Stage 2F v3 analysis used a custom Gibbs sampler. This public
# R implementation uses brms/Stan for a conventional and transparent Bayesian
# multilevel implementation. Numerical posterior summaries should be very close
# but will not be bit-for-bit identical because the sampler and prior
# parameterisation differ.
#
# Packages ---------------------------------------------------------------

library(brms)
library(tidyverse)
library(posterior)

# Reproducibility settings -----------------------------------------------

options(mc.cores = parallel::detectCores())

seed_set <- 260817
chains_set <- 4
iter_set <- 8000
warmup_set <- 3000

# Input/output ------------------------------------------------------------

data_file <- "MUAP_FiringRate_Data.csv"
reference_file <- "MUAP_FiringRate_ReferenceResults.csv"

output_dir <- "MUAP_FiringRate_outputs"
dir.create(output_dir, showWarnings = FALSE)

# Data -------------------------------------------------------------------

d_raw <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  as_tibble()

required_columns <- c(
  "participant_id",
  "recording_id",
  "bout_id",
  "original_mu_number",
  "timepoint",
  "limb",
  "muscle",
  "contraction_intensity_percent",
  "plateau_n_firings",
  "plateau_train_rate_hz",
  "plateau_mean_instantaneous_rate_hz",
  "muap_peak_to_peak_mean_raw",
  "primary_eligible",
  "sensitivity_5plus_eligible",
  "sensitivity_Aclean_eligible",
  "sensitivity_A_or_B_eligible"
)

missing_columns <- setdiff(required_columns, names(d_raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# Helper: prepare one inferential dataset --------------------------------
#
# Within-bout standardisation is recalculated AFTER the dataset-specific
# eligibility restriction. Bouts with <2 usable MUs or zero within-bout
# log-MUAP SD cannot contribute to a within-bout MUAP slope.

prepare_muap_dataset <- function(
    data,
    eligibility_column,
    outcome_column
) {

  data %>%
    filter(
      .data[[eligibility_column]] == "Yes",
      is.finite(.data[[outcome_column]]),
      is.finite(muap_peak_to_peak_mean_raw),
      muap_peak_to_peak_mean_raw > 0
    ) %>%
    mutate(
      original_mu_number = as.integer(original_mu_number),
      recording_mu_id = interaction(
        recording_id,
        original_mu_number,
        drop = TRUE,
        lex.order = TRUE
      ),
      log_muap = log(muap_peak_to_peak_mean_raw)
    ) %>%
    group_by(bout_id) %>%
    mutate(
      model_bout_n_mus = n(),
      bout_log_muap_sd = sd(log_muap),
      z_log_muap_within_bout =
        (log_muap - mean(log_muap)) / bout_log_muap_sd
    ) %>%
    ungroup() %>%
    filter(
      model_bout_n_mus >= 2,
      is.finite(bout_log_muap_sd),
      bout_log_muap_sd > 0,
      is.finite(z_log_muap_within_bout)
    ) %>%
    mutate(
      participant_id = factor(participant_id),
      recording_mu_id = factor(recording_mu_id),
      bout_id = factor(bout_id),

      # Protocol categories
      intensity = factor(
        as.character(contraction_intensity_percent),
        levels = c("20", "50", "75")
      ),
      muscle = factor(muscle, levels = c("VM", "VL")),
      limb = factor(limb, levels = c("Opp", "ACL")),
      timepoint = factor(
        timepoint,
        levels = c("Pre-op", "2wks", "6wks", "3months", "6months")
      ),

      # Explicit dummy variables, then centered. This mirrors the final
      # Stage 2F v3 parameterisation and avoids relying on contrast defaults.
      i50 = as.numeric(intensity == "50"),
      i75 = as.numeric(intensity == "75"),
      vl = as.numeric(muscle == "VL"),
      acl = as.numeric(limb == "ACL"),
      t2 = as.numeric(timepoint == "2wks"),
      t6w = as.numeric(timepoint == "6wks"),
      t3m = as.numeric(timepoint == "3months"),
      t6m = as.numeric(timepoint == "6months"),

      i50_c = i50 - mean(i50),
      i75_c = i75 - mean(i75),
      vl_c = vl - mean(vl),
      acl_c = acl - mean(acl),
      t2_c = t2 - mean(t2),
      t6w_c = t6w - mean(t6w),
      t3m_c = t3m - mean(t3m),
      t6m_c = t6m - mean(t6m),

      # Center rather than scale the raw-Hz outcome. Slopes remain in Hz.
      y = .data[[outcome_column]],
      y_centered = y - mean(y),

      # Raw ACL indicator is retained for the interaction so that the
      # z_log_muap coefficient in the limb model is the contralateral slope.
      z_muap_x_acl = z_log_muap_within_bout * acl
    )
}

# Prepare analysis datasets ----------------------------------------------

d_primary <- prepare_muap_dataset(
  d_raw,
  "primary_eligible",
  "plateau_train_rate_hz"
)

d_5plus <- prepare_muap_dataset(
  d_raw,
  "sensitivity_5plus_eligible",
  "plateau_train_rate_hz"
)

d_A_or_B <- prepare_muap_dataset(
  d_raw,
  "sensitivity_A_or_B_eligible",
  "plateau_train_rate_hz"
)

d_A_clean <- prepare_muap_dataset(
  d_raw,
  "sensitivity_Aclean_eligible",
  "plateau_train_rate_hz"
)

d_mean_instantaneous <- prepare_muap_dataset(
  d_raw,
  "primary_eligible",
  "plateau_mean_instantaneous_rate_hz"
)

# Coverage checks ---------------------------------------------------------

coverage <- bind_rows(
  primary = d_primary,
  sensitivity_5plus = d_5plus,
  sensitivity_A_or_B = d_A_or_B,
  sensitivity_A_clean = d_A_clean,
  sensitivity_mean_instantaneous_rate = d_mean_instantaneous,
  .id = "dataset"
) %>%
  group_by(dataset) %>%
  summarise(
    n_mu_observations = n(),
    n_unique_recording_mus = n_distinct(recording_mu_id),
    n_bouts = n_distinct(bout_id),
    n_participants = n_distinct(participant_id),
    .groups = "drop"
  )

write.csv(
  coverage,
  file.path(output_dir, "R_model_dataset_coverage.csv"),
  row.names = FALSE
)

print(coverage)

# Expected primary structure from final Stage 2F v3:
# 4241 observations, 3806 unique recording-MUs, 510 bouts, 25 participants.
expected_primary <- c(
  n_mu_observations = 4241,
  n_unique_recording_mus = 3806,
  n_bouts = 510,
  n_participants = 25
)

observed_primary <- coverage %>%
  filter(dataset == "primary") %>%
  select(
    n_mu_observations,
    n_unique_recording_mus,
    n_bouts,
    n_participants
  ) %>%
  unlist(use.names = TRUE)

if (!all(as.numeric(observed_primary) == as.numeric(expected_primary))) {
  warning(
    "Primary dataset structure does not match the final Stage 2F v3 ",
    "reference counts. Inspect eligibility and factor coding before ",
    "interpreting the model."
  )
}

# Model formulas ----------------------------------------------------------
#
# No population-level intercept is fitted. The outcome has been centered and
# nuisance dummy variables are centered. This parallels the final Stage 2F v3
# computational parameterisation.

formula_primary <- bf(
  y_centered ~ 0 +
    z_log_muap_within_bout +
    i50_c + i75_c +
    vl_c +
    acl_c +
    t2_c + t6w_c + t3m_c + t6m_c +
    (1 | participant_id) +
    (1 | recording_mu_id) +
    (1 | bout_id)
)

formula_limb <- bf(
  y_centered ~ 0 +
    z_log_muap_within_bout +
    i50_c + i75_c +
    vl_c +
    acl_c +
    t2_c + t6w_c + t3m_c + t6m_c +
    z_muap_x_acl +
    (1 | participant_id) +
    (1 | recording_mu_id) +
    (1 | bout_id)
)

# Priors -----------------------------------------------------------------
#
# Coefficients are in Hz because only the outcome mean was removed.
# These are deliberately weakly informative relative to plausible motor-unit
# firing-rate differences.

priors_primary <- c(
  prior(normal(0, 5), class = "b"),
  prior(student_t(3, 0, 5), class = "sd"),
  prior(student_t(3, 0, 5), class = "sigma")
)

# Fit helper --------------------------------------------------------------

fit_model <- function(formula, data, file_stem) {

  brm(
    formula = formula,
    data = data,
    family = gaussian(),
    prior = priors_primary,
    chains = chains_set,
    iter = iter_set,
    warmup = warmup_set,
    seed = seed_set,
    backend = "rstan",
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    save_pars = save_pars(all = TRUE),
    file = file.path(output_dir, file_stem),
    file_refit = "on_change"
  )
}

# Fit final model set -----------------------------------------------------

fit_primary <- fit_model(
  formula_primary,
  d_primary,
  "fit_primary"
)

fit_limb <- fit_model(
  formula_limb,
  d_primary,
  "fit_limb_modifier"
)

fit_5plus <- fit_model(
  formula_primary,
  d_5plus,
  "fit_sensitivity_5plus"
)

fit_A_or_B <- fit_model(
  formula_primary,
  d_A_or_B,
  "fit_sensitivity_A_or_B"
)

fit_A_clean <- fit_model(
  formula_primary,
  d_A_clean,
  "fit_sensitivity_A_clean"
)

fit_mean_instantaneous <- fit_model(
  formula_primary,
  d_mean_instantaneous,
  "fit_sensitivity_mean_instantaneous"
)

# Posterior summary helper ------------------------------------------------

summarise_parameter <- function(
    fit,
    parameter,
    analysis,
    data
) {

  draws <- as_draws_df(fit)[[parameter]]

  tibble(
    analysis = analysis,
    parameter = parameter,
    posterior_mean_Hz = mean(draws),
    posterior_sd_Hz = sd(draws),
    lower_95_CrI_Hz = quantile(draws, 0.025),
    upper_95_CrI_Hz = quantile(draws, 0.975),
    Pr_gt_0 = mean(draws > 0),
    Pr_lt_0 = mean(draws < 0),
    R_hat = rhat(draws),
    n_MU_observations = nrow(data),
    n_unique_recording_MUs = n_distinct(data$recording_mu_id),
    n_bouts = n_distinct(data$bout_id),
    n_participants = n_distinct(data$participant_id)
  )
}

# brms names population-level coefficients with b_
primary_parameter <- "b_z_log_muap_within_bout"
limb_interaction_parameter <- "b_z_muap_x_acl"

publication_results <- bind_rows(
  summarise_parameter(
    fit_primary,
    primary_parameter,
    "Primary",
    d_primary
  ),
  summarise_parameter(
    fit_5plus,
    primary_parameter,
    "Sensitivity: >=5 plateau firings",
    d_5plus
  ),
  summarise_parameter(
    fit_A_or_B,
    primary_parameter,
    "Sensitivity: A-clean + B-segmentation review",
    d_A_or_B
  ),
  summarise_parameter(
    fit_A_clean,
    primary_parameter,
    "Sensitivity: A-clean only",
    d_A_clean
  ),
  summarise_parameter(
    fit_mean_instantaneous,
    primary_parameter,
    "Sensitivity: mean instantaneous firing-rate outcome",
    d_mean_instantaneous
  ),
  summarise_parameter(
    fit_limb,
    limb_interaction_parameter,
    "Secondary: limb effect modification",
    d_primary
  )
)

write.csv(
  publication_results,
  file.path(output_dir, "R_publication_results.csv"),
  row.names = FALSE
)

print(publication_results)

# Limb-specific MUAP slopes ----------------------------------------------

limb_draws <- as_draws_df(fit_limb)

opp_slope <- limb_draws$b_z_log_muap_within_bout
acl_difference <- limb_draws$b_z_muap_x_acl
acl_slope <- opp_slope + acl_difference

summarise_draw_vector <- function(draws, label) {
  tibble(
    slope = label,
    posterior_mean_Hz_per_1SD_logMUAP = mean(draws),
    posterior_sd = sd(draws),
    lower_95_CrI = quantile(draws, 0.025),
    upper_95_CrI = quantile(draws, 0.975),
    Pr_gt_0 = mean(draws > 0),
    Pr_lt_0 = mean(draws < 0)
  )
}

limb_slopes <- bind_rows(
  summarise_draw_vector(opp_slope, "Opp"),
  summarise_draw_vector(acl_slope, "ACL"),
  summarise_draw_vector(acl_difference, "ACL_minus_Opp")
)

write.csv(
  limb_slopes,
  file.path(output_dir, "R_limb_modifier_slopes.csv"),
  row.names = FALSE
)

print(limb_slopes)

# Diagnostics -------------------------------------------------------------

diagnostics <- bind_rows(
  Primary = as.data.frame(summary(fit_primary)$fixed) %>%
    rownames_to_column("parameter"),
  Limb_modifier = as.data.frame(summary(fit_limb)$fixed) %>%
    rownames_to_column("parameter"),
  Sensitivity_5plus = as.data.frame(summary(fit_5plus)$fixed) %>%
    rownames_to_column("parameter"),
  Sensitivity_A_or_B = as.data.frame(summary(fit_A_or_B)$fixed) %>%
    rownames_to_column("parameter"),
  Sensitivity_A_clean = as.data.frame(summary(fit_A_clean)$fixed) %>%
    rownames_to_column("parameter"),
  Sensitivity_mean_instantaneous =
    as.data.frame(summary(fit_mean_instantaneous)$fixed) %>%
    rownames_to_column("parameter"),
  .id = "model"
)

write.csv(
  diagnostics,
  file.path(output_dir, "R_fixed_effect_diagnostics.csv"),
  row.names = FALSE
)

# Compare with the archived Stage 2F v3 reference estimates ---------------
#
# The reference file contains the final custom-Gibbs results used to verify
# that this independent R/brms implementation gives the same substantive
# conclusions. Exact equality is not expected.

if (file.exists(reference_file)) {

  reference <- read.csv(reference_file, stringsAsFactors = FALSE)

  comparison <- publication_results %>%
    select(
      analysis,
      R_posterior_mean_Hz = posterior_mean_Hz,
      R_lower_95_CrI_Hz = lower_95_CrI_Hz,
      R_upper_95_CrI_Hz = upper_95_CrI_Hz
    ) %>%
    left_join(
      reference %>%
        select(
          analysis,
          reference_posterior_mean_Hz = posterior_mean_Hz,
          reference_lower_95_CrI_Hz = lower_95_CrI_Hz,
          reference_upper_95_CrI_Hz = upper_95_CrI_Hz
        ),
      by = "analysis"
    ) %>%
    mutate(
      mean_difference_Hz =
        R_posterior_mean_Hz - reference_posterior_mean_Hz
    )

  write.csv(
    comparison,
    file.path(output_dir, "R_vs_Stage2Fv3_reference.csv"),
    row.names = FALSE
  )

  print(comparison)
}

# Posterior predictive checks --------------------------------------------

pdf(file.path(output_dir, "posterior_predictive_checks.pdf"))

print(pp_check(fit_primary, ndraws = 100))
print(pp_check(fit_limb, ndraws = 100))

dev.off()

# Save textual model summaries -------------------------------------------

sink(file.path(output_dir, "R_model_summaries.txt"))

cat("PRIMARY MODEL\n")
print(summary(fit_primary))

cat("\n\nLIMB MODIFIER MODEL\n")
print(summary(fit_limb))

cat("\n\nSENSITIVITY >=5 FIRINGS\n")
print(summary(fit_5plus))

cat("\n\nSENSITIVITY A+B\n")
print(summary(fit_A_or_B))

cat("\n\nSENSITIVITY A-CLEAN\n")
print(summary(fit_A_clean))

cat("\n\nSENSITIVITY MEAN INSTANTANEOUS RATE\n")
print(summary(fit_mean_instantaneous))

sink()

# Record software versions ------------------------------------------------

sink(file.path(output_dir, "sessionInfo_MUAP_analysis.txt"))
sessionInfo()
sink()

# Final console message ---------------------------------------------------

cat("\n================================================\n")
cat("MUAP FIRING-RATE ANALYSIS COMPLETE\n")
cat("================================================\n")
cat("Primary model rows:", nrow(d_primary), "\n")
cat(
  "Primary unique recording-MUs:",
  n_distinct(d_primary$recording_mu_id),
  "\n"
)
cat("Primary bouts:", n_distinct(d_primary$bout_id), "\n")
cat(
  "Primary participants:",
  n_distinct(d_primary$participant_id),
  "\n"
)
cat("\nResults written to:", output_dir, "\n")
cat("================================================\n")
