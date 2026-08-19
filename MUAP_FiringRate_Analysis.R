# Secondary MUAP–firing-rate analysis
# Prospective ACL motor-unit study
#
# This script reproduces the final crossed-hierarchical analysis used
# to examine the relationship between surface-detected MUAP amplitude
# and motor-unit firing rate.
#
# Required packages:
#   rstan
#   tidyverse
#
# Input:
#   MUAP_FiringRate_Data.csv
#
# Model structure:
#   fixed effects:
#     within-bout z(log MUAP), intensity, muscle, limb, timepoint
#   random intercepts:
#     participant, original recording-MU, contraction bout
#
# Secondary model:
#   adds within-bout MUAP x limb interaction
#
# Sensitivity analyses:
#   >=5 plateau firings
#   A-clean + B-segmentation-review
#   A-clean only
#   mean instantaneous firing-rate outcome
#
# MUAP amplitude is interpreted as a surface-detected
# electrophysiological characteristic, not as a direct anatomical
# measure of motor-unit size.

library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

chains_set <- 4
iter_set <- 8000
warmup_set <- 3000
thin_set <- 5
seed_set <- 260817

coefficient_prior_sd <- 2.5
variance_prior_shape <- 2.5
variance_prior_rate <- 0.5

data_file <- "MUAP_FiringRate_Data.csv"

# ------------------------------------------------------------
# Load and type data
# ------------------------------------------------------------

d <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  mutate(
    participant_id = factor(participant_id),
    recording_id = factor(recording_id),
    original_mu_number = as.integer(original_mu_number),
    recording_mu_id = factor(
      paste(recording_id, sprintf("MU%03d", original_mu_number), sep = "_")
    ),
    bout_id = factor(bout_id),
    contraction_intensity_percent = as.character(contraction_intensity_percent),
    timepoint = factor(
      timepoint,
      levels = c("Pre-op", "2wks", "6wks", "3months", "6months")
    ),
    limb = factor(limb, levels = c("Opp", "ACL")),
    muscle = factor(muscle, levels = c("VM", "VL")),
    plateau_n_firings = as.numeric(plateau_n_firings),
    plateau_train_rate_hz = as.numeric(plateau_train_rate_hz),
    plateau_mean_instantaneous_rate_hz =
      as.numeric(plateau_mean_instantaneous_rate_hz),
    muap_peak_to_peak_mean_raw = as.numeric(muap_peak_to_peak_mean_raw)
  )

# ------------------------------------------------------------
# Prepare an inferential dataset
# ------------------------------------------------------------

prepare_dataset <- function(data, eligibility_variable, outcome_variable) {

  data %>%
    filter(
      .data[[eligibility_variable]] == "Yes",
      is.finite(.data[[outcome_variable]]),
      is.finite(muap_peak_to_peak_mean_raw),
      muap_peak_to_peak_mean_raw > 0
    ) %>%
    mutate(
      log_muap = log(muap_peak_to_peak_mean_raw)
    ) %>%
    group_by(bout_id) %>%
    mutate(
      n_usable_mus_in_bout = n(),
      bout_log_muap_sd = sd(log_muap)
    ) %>%
    filter(
      n_usable_mus_in_bout >= 2,
      is.finite(bout_log_muap_sd),
      bout_log_muap_sd > 0
    ) %>%
    mutate(
      z_log_muap_within_bout =
        (log_muap - mean(log_muap)) / sd(log_muap)
    ) %>%
    ungroup() %>%
    mutate(
      outcome = .data[[outcome_variable]]
    )
}

primary <- prepare_dataset(
  d,
  "primary_eligible",
  "plateau_train_rate_hz"
)

sens_5plus <- prepare_dataset(
  d,
  "sensitivity_5plus_eligible",
  "plateau_train_rate_hz"
)

sens_A_or_B <- prepare_dataset(
  d,
  "sensitivity_A_or_B_eligible",
  "plateau_train_rate_hz"
)

sens_A_clean <- prepare_dataset(
  d,
  "sensitivity_A_clean_eligible",
  "plateau_train_rate_hz"
)

sens_mean_inst <- prepare_dataset(
  d,
  "primary_eligible",
  "plateau_mean_instantaneous_rate_hz"
)

# Expected primary coverage from the final analysis:
# 4241 MU observations
# 3806 unique original recording-MUs
# 510 bouts
# 25 participants

coverage <- function(x, label) {
  tibble(
    dataset = label,
    n_mu_observations = nrow(x),
    n_unique_recording_mus = n_distinct(x$recording_mu_id),
    n_bouts = n_distinct(x$bout_id),
    n_participants = n_distinct(x$participant_id)
  )
}

coverage_table <- bind_rows(
  coverage(primary, "primary"),
  coverage(sens_5plus, "sensitivity_5plus"),
  coverage(sens_A_or_B, "sensitivity_A_or_B"),
  coverage(sens_A_clean, "sensitivity_A_clean"),
  coverage(sens_mean_inst, "sensitivity_mean_instantaneous_rate")
)

write.csv(
  coverage_table,
  "MUAP_FiringRate_R_model_coverage.csv",
  row.names = FALSE
)

print(coverage_table)

# ------------------------------------------------------------
# Build the fixed-effect design matrix
# ------------------------------------------------------------
#
# There is deliberately no global intercept in the Stan model.
# The outcome is standardized internally, and nuisance dummy
# predictors are mean-centered. The interaction uses the raw ACL
# indicator so the z-MUAP coefficient in the limb-modifier model
# represents the contralateral-limb slope.

build_design <- function(data, limb_modifier = FALSE) {

  temp <- data %>%
    mutate(
      i50 = as.numeric(contraction_intensity_percent == "50"),
      i75 = as.numeric(contraction_intensity_percent == "75"),
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
      t6m_c = t6m - mean(t6m)
    )

  X <- temp %>%
    transmute(
      z_log_MUAP_within_bout = z_log_muap_within_bout,
      Intensity_50_vs_20 = i50_c,
      Intensity_75_vs_20 = i75_c,
      Muscle_VL_vs_VM = vl_c,
      Limb_ACL_vs_Opp = acl_c,
      Time_2wks_vs_Preop = t2_c,
      Time_6wks_vs_Preop = t6w_c,
      Time_3months_vs_Preop = t3m_c,
      Time_6months_vs_Preop = t6m_c
    )

  if (limb_modifier) {
    X$zMUAP_x_LimbACL <-
      temp$z_log_muap_within_bout * temp$acl
  }

  as.matrix(X)
}

# ------------------------------------------------------------
# Stan model
# ------------------------------------------------------------
#
# This encodes the final Stage 2F v3 model directly.
# Variance components have Inv-Gamma(2.5, 0.5) priors on
# standardized outcome variance, matching the final analysis plan.

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;

  int<lower=1> J_participant;
  array[N] int<lower=1, upper=J_participant> participant;

  int<lower=1> J_mu;
  array[N] int<lower=1, upper=J_mu> recording_mu;

  int<lower=1> J_bout;
  array[N] int<lower=1, upper=J_bout> bout;

  real<lower=0> coefficient_prior_sd;
  real<lower=0> variance_prior_shape;
  real<lower=0> variance_prior_rate;
}
parameters {
  vector[K] beta;

  vector[J_participant] participant_re;
  vector[J_mu] recording_mu_re;
  vector[J_bout] bout_re;

  real<lower=0> sigma2_participant;
  real<lower=0> sigma2_recording_mu;
  real<lower=0> sigma2_bout;
  real<lower=0> sigma2_error;
}
model {
  beta ~ normal(0, coefficient_prior_sd);

  sigma2_participant ~ inv_gamma(
    variance_prior_shape,
    variance_prior_rate
  );
  sigma2_recording_mu ~ inv_gamma(
    variance_prior_shape,
    variance_prior_rate
  );
  sigma2_bout ~ inv_gamma(
    variance_prior_shape,
    variance_prior_rate
  );
  sigma2_error ~ inv_gamma(
    variance_prior_shape,
    variance_prior_rate
  );

  participant_re ~ normal(0, sqrt(sigma2_participant));
  recording_mu_re ~ normal(0, sqrt(sigma2_recording_mu));
  bout_re ~ normal(0, sqrt(sigma2_bout));

  y ~ normal(
    X * beta
      + participant_re[participant]
      + recording_mu_re[recording_mu]
      + bout_re[bout],
    sqrt(sigma2_error)
  );
}
"

stan_mod <- stan_model(model_code = stan_code)

# ------------------------------------------------------------
# Model fitting helper
# ------------------------------------------------------------

fit_muap_model <- function(
  data,
  model_name,
  limb_modifier = FALSE,
  seed = seed_set
) {

  X <- build_design(
    data,
    limb_modifier = limb_modifier
  )

  y_mean <- mean(data$outcome)
  y_sd <- sd(data$outcome)
  y_std <- (data$outcome - y_mean) / y_sd

  participant_index <- as.integer(factor(data$participant_id))
  mu_index <- as.integer(factor(data$recording_mu_id))
  bout_index <- as.integer(factor(data$bout_id))

  stan_data <- list(
    N = nrow(data),
    K = ncol(X),
    X = X,
    y = as.vector(y_std),

    J_participant = max(participant_index),
    participant = participant_index,

    J_mu = max(mu_index),
    recording_mu = mu_index,

    J_bout = max(bout_index),
    bout = bout_index,

    coefficient_prior_sd = coefficient_prior_sd,
    variance_prior_shape = variance_prior_shape,
    variance_prior_rate = variance_prior_rate
  )

  fit <- sampling(
    stan_mod,
    data = stan_data,
    chains = chains_set,
    iter = iter_set,
    warmup = warmup_set,
    thin = thin_set,
    seed = seed,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 15
    ),
    refresh = 250
  )

  beta_draws_std <- as.matrix(fit, pars = "beta")
  beta_draws_hz <- beta_draws_std * y_sd
  colnames(beta_draws_hz) <- colnames(X)

  stan_summary <- summary(
    fit,
    pars = c(
      "beta",
      "sigma2_participant",
      "sigma2_recording_mu",
      "sigma2_bout",
      "sigma2_error"
    )
  )$summary

  beta_rhat <- stan_summary[
    paste0("beta[", seq_len(ncol(X)), "]"),
    "Rhat"
  ]

  fixed_summary <- map_dfr(
    seq_len(ncol(beta_draws_hz)),
    function(j) {

      draws <- beta_draws_hz[, j]

      tibble(
        model = model_name,
        parameter_type = "fixed_effect",
        parameter = colnames(beta_draws_hz)[j],
        posterior_mean = mean(draws),
        posterior_sd = sd(draws),
        posterior_median = median(draws),
        lower_95_CrI = quantile(draws, 0.025),
        upper_95_CrI = quantile(draws, 0.975),
        Pr_gt_0 = mean(draws > 0),
        Pr_lt_0 = mean(draws < 0),
        R_hat = beta_rhat[j],
        n_rows = nrow(data),
        n_participants = n_distinct(data$participant_id),
        n_recording_MUs = n_distinct(data$recording_mu_id),
        n_bouts = n_distinct(data$bout_id)
      )
    }
  )

  variance_names <- c(
    sigma2_participant = "SD_participant_intercept",
    sigma2_recording_mu = "SD_recording_MU_intercept",
    sigma2_bout = "SD_bout_intercept",
    sigma2_error = "SD_residual"
  )

  variance_summary <- map_dfr(
    names(variance_names),
    function(par) {

      variance_draws <- as.vector(
        as.matrix(fit, pars = par)
      )

      sd_draws_hz <- sqrt(variance_draws) * y_sd

      tibble(
        model = model_name,
        parameter_type = "standard_deviation",
        parameter = variance_names[[par]],
        posterior_mean = mean(sd_draws_hz),
        posterior_sd = sd(sd_draws_hz),
        posterior_median = median(sd_draws_hz),
        lower_95_CrI = quantile(sd_draws_hz, 0.025),
        upper_95_CrI = quantile(sd_draws_hz, 0.975),
        Pr_gt_0 = mean(sd_draws_hz > 0),
        Pr_lt_0 = mean(sd_draws_hz < 0),
        R_hat = stan_summary[par, "Rhat"],
        n_rows = nrow(data),
        n_participants = n_distinct(data$participant_id),
        n_recording_MUs = n_distinct(data$recording_mu_id),
        n_bouts = n_distinct(data$bout_id)
      )
    }
  )

  list(
    fit = fit,
    fixed_draws_hz = beta_draws_hz,
    summary = bind_rows(fixed_summary, variance_summary),
    data = data,
    X = X
  )
}

# ------------------------------------------------------------
# Fit final model set
# ------------------------------------------------------------

fit_primary <- fit_muap_model(
  primary,
  "primary_overall",
  limb_modifier = FALSE,
  seed = seed_set + 1
)

fit_limb <- fit_muap_model(
  primary,
  "secondary_limb_modifier",
  limb_modifier = TRUE,
  seed = seed_set + 2
)

fit_5plus <- fit_muap_model(
  sens_5plus,
  "sensitivity_5plus_firings",
  limb_modifier = FALSE,
  seed = seed_set + 3
)

fit_A_or_B <- fit_muap_model(
  sens_A_or_B,
  "sensitivity_A_or_B",
  limb_modifier = FALSE,
  seed = seed_set + 4
)

fit_A_clean <- fit_muap_model(
  sens_A_clean,
  "sensitivity_A_clean",
  limb_modifier = FALSE,
  seed = seed_set + 5
)

fit_mean_inst <- fit_muap_model(
  sens_mean_inst,
  "sensitivity_mean_instantaneous_rate",
  limb_modifier = FALSE,
  seed = seed_set + 6
)

all_model_summary <- bind_rows(
  fit_primary$summary,
  fit_limb$summary,
  fit_5plus$summary,
  fit_A_or_B$summary,
  fit_A_clean$summary,
  fit_mean_inst$summary
)

write.csv(
  all_model_summary,
  "MUAP_FiringRate_R_model_posterior_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Derived limb-specific slopes
# ------------------------------------------------------------

opp_slope <- fit_limb$fixed_draws_hz[
  ,
  "z_log_MUAP_within_bout"
]

limb_difference <- fit_limb$fixed_draws_hz[
  ,
  "zMUAP_x_LimbACL"
]

acl_slope <- opp_slope + limb_difference

summarise_draws <- function(draws, label) {
  tibble(
    slope = label,
    posterior_mean_Hz_per_1SD_logMUAP = mean(draws),
    posterior_sd = sd(draws),
    posterior_median = median(draws),
    lower_95_CrI = quantile(draws, 0.025),
    upper_95_CrI = quantile(draws, 0.975),
    Pr_gt_0 = mean(draws > 0),
    Pr_lt_0 = mean(draws < 0)
  )
}

limb_slopes <- bind_rows(
  summarise_draws(opp_slope, "Opp"),
  summarise_draws(acl_slope, "ACL"),
  summarise_draws(limb_difference, "ACL_minus_Opp")
)

write.csv(
  limb_slopes,
  "MUAP_FiringRate_R_limb_modifier_slopes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Compact publication/reviewer results
# ------------------------------------------------------------

get_fixed <- function(fit_object, parameter) {
  fit_object$summary %>%
    filter(
      parameter_type == "fixed_effect",
      .data$parameter == parameter
    )
}

publication_results <- bind_rows(
  get_fixed(
    fit_primary,
    "z_log_MUAP_within_bout"
  ) %>%
    mutate(analysis = "Primary"),

  get_fixed(
    fit_5plus,
    "z_log_MUAP_within_bout"
  ) %>%
    mutate(analysis = "Sensitivity: >=5 plateau firings"),

  get_fixed(
    fit_A_or_B,
    "z_log_MUAP_within_bout"
  ) %>%
    mutate(
      analysis =
        "Sensitivity: A-clean + B-segmentation review"
    ),

  get_fixed(
    fit_A_clean,
    "z_log_MUAP_within_bout"
  ) %>%
    mutate(analysis = "Sensitivity: A-clean only"),

  get_fixed(
    fit_mean_inst,
    "z_log_MUAP_within_bout"
  ) %>%
    mutate(
      analysis =
        "Sensitivity: mean instantaneous firing-rate outcome"
    ),

  get_fixed(
    fit_limb,
    "zMUAP_x_LimbACL"
  ) %>%
    mutate(analysis = "Secondary: limb effect modification")
) %>%
  select(
    analysis,
    parameter,
    posterior_mean,
    lower_95_CrI,
    upper_95_CrI,
    Pr_gt_0,
    Pr_lt_0,
    R_hat,
    n_rows,
    n_recording_MUs,
    n_bouts,
    n_participants
  )

write.csv(
  publication_results,
  "MUAP_FiringRate_R_publication_results.csv",
  row.names = FALSE
)

print(publication_results)
print(limb_slopes)

# ------------------------------------------------------------
# Convergence checks
# ------------------------------------------------------------

convergence <- all_model_summary %>%
  select(model, parameter, R_hat) %>%
  mutate(
    warning = if_else(
      is.na(R_hat) | R_hat <= 1.05,
      "",
      "R_hat_above_1.05"
    )
  )

write.csv(
  convergence,
  "MUAP_FiringRate_R_model_diagnostics.csv",
  row.names = FALSE
)

if (any(convergence$warning != "")) {
  warning(
    "One or more monitored parameters have R-hat > 1.05. ",
    "Inspect MUAP_FiringRate_R_model_diagnostics.csv before ",
    "using the results."
  )
}

# ------------------------------------------------------------
# Save fitted objects and R session information
# ------------------------------------------------------------

saveRDS(
  list(
    primary = fit_primary$fit,
    limb_modifier = fit_limb$fit,
    sensitivity_5plus = fit_5plus$fit,
    sensitivity_A_or_B = fit_A_or_B$fit,
    sensitivity_A_clean = fit_A_clean$fit,
    sensitivity_mean_instantaneous_rate = fit_mean_inst$fit
  ),
  "MUAP_FiringRate_R_fits.rds"
)

capture.output(
  sessionInfo(),
  file = "sessionInfo_MUAP_analysis.txt"
)

# ------------------------------------------------------------
# Reference values from the original Stage 2F v3 implementation
# ------------------------------------------------------------
#
# These values provide an independent check after running this R/Stan
# implementation. Small Monte Carlo differences are expected.
#
# Primary:
#   -3.7128 Hz, 95% CrI -3.7839 to -3.6392
#
# >=5 firings:
#   -3.7383 Hz, 95% CrI -3.8089 to -3.6677
#
# A-clean + B:
#   -3.7010 Hz, 95% CrI -3.7731 to -3.6287
#
# A-clean:
#   -3.7450 Hz, 95% CrI -3.8172 to -3.6717
#
# Mean instantaneous rate:
#   -3.8013 Hz, 95% CrI -3.8781 to -3.7221
#
# ACL minus contralateral slope difference:
#    0.4680 Hz, 95% CrI 0.3250 to 0.6078
#
# Contralateral slope:
#   -3.9298 Hz, 95% CrI -4.0291 to -3.8318
#
# ACL slope:
#   -3.4618 Hz, 95% CrI -3.5647 to -3.3586
