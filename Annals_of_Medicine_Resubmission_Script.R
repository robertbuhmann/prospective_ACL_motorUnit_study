# ================================================================
# PRIMARY BAYESIAN ANALYSIS
# Prospective ACL motor-unit study
#
# Corrected analysis:
# - excludes 100% MVC motor-unit observations
# - analyses intended submaximal contractions only:
#       20%, 50%, and 75% MVC
# - uses CART multiple imputation across 5 datasets
# - six-month outcome is NOT imputed
# - model is fitted on recorded kg-force scale
# - model summaries are converted to Newtons for reporting
# - ROPE = +/- 15 N = +/- 1.5296 kg-force
#
# ================================================================


# ================================================================
# 1. PACKAGES
# ================================================================

library(mice)
library(rstan)
library(brms)
library(tidyverse)
library(janitor)
library(naniar)
library(posterior)
library(purrr)


# ================================================================
# 2. SETTINGS
# ================================================================

cores_set  <- 8
iter_set   <- 4000
chains_set <- 4
warmup_set <- 2000
seed_set   <- 123

# Force conversion
kgf_to_N <- 9.80665

# Region of practical equivalence
rope_N   <- 15
rope_kgf <- rope_N / kgf_to_N

rope_lower <- -rope_kgf
rope_upper <-  rope_kgf

cat(
  "ROPE used for posterior calculations:",
  round(rope_lower, 4), "to",
  round(rope_upper, 4), "kg-force\n"
)

cat(
  "Equivalent to:",
  -rope_N, "to", rope_N, "N\n"
)


# ================================================================
# 3. LOAD DATA
# ================================================================

dat_raw <- read.csv("ModellingData.csv") %>%
  clean_names()


# ------------------------------------------------
# Check required variables are present
# ------------------------------------------------

required_variables <- c(
  "force",
  "force_six_months",
  "limb",
  "muscle",
  "mean_fr",
  "time",
  "sex",
  "id",
  "contraction_intensity_s",
  "mass"
)

missing_variables <- setdiff(
  required_variables,
  names(dat_raw)
)

if (length(missing_variables) > 0) {
  stop(
    paste(
      "Missing required variables:",
      paste(missing_variables, collapse = ", ")
    )
  )
}


# ================================================================
# 4. PREPARE INTENDED SUBMAXIMAL ANALYSIS DATA
# ================================================================

analysis_raw <- dat_raw %>%
  
  mutate(
    
    time =
      as.numeric(as.character(time)),
    
    contraction_intensity_s =
      as.numeric(as.character(contraction_intensity_s)),
    
    force =
      as.numeric(as.character(force)),
    
    force_six_months =
      as.numeric(as.character(force_six_months)),
    
    mean_fr =
      as.numeric(as.character(mean_fr)),
    
    peak_fr =
      as.numeric(as.character(peak_fr)),
    
    mass =
      as.numeric(as.character(mass))
    
  ) %>%
  
  # ------------------------------------------------
# IMPORTANT:
# Predictor measurements are pre-op through
# 3 months only.
# Six-month motor-unit observations are excluded.
# ------------------------------------------------

filter(time %in% c(1, 2, 3, 4)) %>%
  
  # ------------------------------------------------
# IMPORTANT CORRECTION:
# Keep SUBMAXIMAL contractions only.
#
# Excludes 100% MVC observations.
# ------------------------------------------------

filter(
  contraction_intensity_s %in%
    c(0.20, 0.50, 0.75)
)


# ================================================================
# 5. VERIFY ANALYSIS DATA
# ================================================================

cat("\nContraction intensities retained:\n")
print(
  sort(
    unique(
      analysis_raw$contraction_intensity_s
    )
  )
)

cat("\nTime points retained:\n")
print(
  sort(
    unique(
      analysis_raw$time
    )
  )
)

cat("\nNumber of participants:\n")
print(
  n_distinct(analysis_raw$id)
)

cat("\nNumber of analysis rows:\n")
print(
  nrow(analysis_raw)
)


# Hard stop if 100% MVC somehow remains
if (
  any(
    analysis_raw$contraction_intensity_s == 1,
    na.rm = TRUE
  )
) {
  stop(
    "ERROR: 100% MVC observations remain in analysis data."
  )
}


# ------------------------------------------------
# Save counts for audit/reproducibility
# ------------------------------------------------

analysis_counts <- analysis_raw %>%
  
  count(
    time,
    contraction_intensity_s,
    limb,
    muscle,
    name = "n_rows"
  )

write.csv(
  analysis_counts,
  "Primary_analysis_data_counts.csv",
  row.names = FALSE
)


# ================================================================
# 6. STANDARDISE FIRING RATE AND BODY MASS
# ================================================================

# IMPORTANT:
# These statistics are calculated AFTER 100% MVC observations
# have been excluded.

mean_fr_mu <- mean(
  analysis_raw$mean_fr,
  na.rm = TRUE
)

mean_fr_sd <- sd(
  analysis_raw$mean_fr,
  na.rm = TRUE
)

mass_mu <- mean(
  analysis_raw$mass,
  na.rm = TRUE
)

mass_sd <- sd(
  analysis_raw$mass,
  na.rm = TRUE
)


cat(
  "\nMean firing rate scaling:\nMean =",
  mean_fr_mu,
  "\nSD =",
  mean_fr_sd,
  "\n"
)

cat(
  "\nBody mass scaling:\nMean =",
  mass_mu,
  "\nSD =",
  mass_sd,
  "\n"
)


# ================================================================
# 7. CREATE MODELLING DATASET
# ================================================================

dsub <- analysis_raw %>%
  
  mutate(
    
    # Standardised variables
    mean_fr_s =
      (mean_fr - mean_fr_mu) /
      mean_fr_sd,
    
    mass_s =
      (mass - mass_mu) /
      mass_sd,
    
    # Explicit factor/reference levels
    id = factor(id),
    
    time = factor(
      time,
      levels = c(1, 2, 3, 4)
    ),
    
    limb = factor(
      limb,
      levels = c("ACL", "Opp")
    ),
    
    muscle = factor(
      muscle,
      levels = c("VL", "VM")
    ),
    
    sex = factor(
      sex,
      levels = c("F", "M")
    ),
    
    contraction_intensity_s = factor(
      contraction_intensity_s,
      levels = c(
        0.20,
        0.50,
        0.75
      )
    )
    
  ) %>%
  
  # Original analysis uses mean_fr_s rather than raw mean_fr
  select(
    -any_of(
      c(
        "peak_fr",
        "mean_fr",
        "force_six_months_norm"
      )
    )
  )


# ================================================================
# 8. MISSING-DATA CHECKS
# ================================================================

cat("\nMissing observations by variable:\n")

print(
  colSums(
    is.na(dsub)
  )
)


# Optional visual checks
vis_miss(dsub)

gg_miss_fct(
  x = dsub,
  fct = time
)


# ================================================================
# 9. CART MULTIPLE IMPUTATION
# ================================================================

# Create method and predictor matrices without first
# running an unnecessary complete imputation.

meth <- make.method(dsub)

pred <- make.predictorMatrix(dsub)

diag(pred) <- 0


# ------------------------------------------------
# Use CART for variables requiring imputation
# ------------------------------------------------

meth[meth != ""] <- "cart"


# ------------------------------------------------
# IMPORTANT:
# Six-month strength is the outcome and is NOT
# imputed.
# ------------------------------------------------

meth["force_six_months"] <- ""


cat("\nImputation methods:\n")
print(meth)


# Run 5 imputations

imp_datasets <- mice(
  dsub,
  m = 5,
  method = meth,
  predictorMatrix = pred,
  seed = seed_set,
  printFlag = TRUE
)


cat("\nMICE logged events:\n")
print(
  imp_datasets$loggedEvents
)


# ================================================================
# 10. PRIMARY BAYESIAN MODEL
# ================================================================

mod_formula <-
  
  force_six_months ~
  
  time +
  
  limb +
  
  muscle +
  
  sex +
  
  contraction_intensity_s +
  
  mass_s +
  
  limb * time * mean_fr_s +
  
  (1 | id)


cat("\nModel formula:\n")
print(mod_formula)


# ------------------------------------------------
# Priors
#
# Preserve original primary analysis:
# beta coefficients = Normal(0, 10)
#
# Remaining parameters use brms defaults.
# ------------------------------------------------

fit <- brm_multiple(
  
  formula = mod_formula,
  
  family = gaussian(),
  
  data = imp_datasets,
  
  prior =
    prior(
      normal(0, 10),
      class = b
    ),
  
  warmup = warmup_set,
  
  cores = cores_set,
  
  iter = iter_set,
  
  chains = chains_set,
  
  seed = seed_set,
  
  control =
    list(
      max_treedepth = 20
    )
)


# ================================================================
# 11. SAVE FINAL MODEL
# ================================================================

saveRDS(
  fit,
  "model_submax_corrected.rds"
)


# ================================================================
# 12. MODEL SUMMARY AND DIAGNOSTICS
# ================================================================

print(
  summary(fit)
)


capture.output(
  summary(fit),
  file =
    "Primary_model_summary.txt"
)


# ------------------------------------------------
# Save priors actually used
# ------------------------------------------------

capture.output(
  prior_summary(fit),
  file =
    "Primary_model_priors.txt"
)


# ------------------------------------------------
# Posterior predictive checks
# ------------------------------------------------

pp1 <- pp_check(
  fit,
  ndraws = 100
)

ggsave(
  "Posterior_predictive_density.png",
  pp1,
  width = 6,
  height = 4,
  dpi = 600
)


pp2 <- pp_check(
  fit,
  type = "error_hist",
  ndraws = 20
)

ggsave(
  "Posterior_predictive_error_hist.png",
  pp2,
  width = 6,
  height = 4,
  dpi = 600
)


# ------------------------------------------------
# Bayesian R-squared
# ------------------------------------------------

r2_results <- bayes_R2(
  fit,
  probs = c(0.025, 0.975)
)

print(
  r2_results
)

write.csv(
  as.data.frame(r2_results),
  "Bayesian_R2.csv",
  row.names = FALSE
)


# ================================================================
# 13. EXTRACT POSTERIOR DRAWS
# ================================================================

post <- as_draws_df(
  fit
)


# ================================================================
# 14. SUPPLEMENTARY TABLE S2
# ALL FIXED EFFECTS
# ================================================================

fixed_effect_columns <-
  
  grep(
    "^b_",
    names(post),
    value = TRUE
  )


fixed_effects <- map_dfr(
  
  fixed_effect_columns,
  
  function(parameter_name) {
    
    draws <-
      post[[parameter_name]]
    
    tibble(
      
      parameter =
        parameter_name,
      
      # Posterior summaries converted to N
      beta_N =
        mean(draws) *
        kgf_to_N,
      
      posterior_SD_N =
        sd(draws) *
        kgf_to_N,
      
      CrI_lower_N =
        as.numeric(
          quantile(
            draws,
            0.025
          )
        ) *
        kgf_to_N,
      
      CrI_upper_N =
        as.numeric(
          quantile(
            draws,
            0.975
          )
        ) *
        kgf_to_N,
      
      # ROPE probabilities calculated on the
      # matching kg-force scale
      Pr_below_ROPE =
        mean(
          draws <
            rope_lower
        ),
      
      Pr_within_ROPE =
        mean(
          draws >=
            rope_lower &
            draws <=
            rope_upper
        ),
      
      Pr_above_ROPE =
        mean(
          draws >
            rope_upper
        ),
      
      Pr_outside_ROPE =
        mean(
          abs(draws) >
            rope_kgf
        )
      
    )
  }
)


write.csv(
  fixed_effects,
  "Supplementary_Table_S2_fixed_effects.csv",
  row.names = FALSE
)


print(
  fixed_effects
)


# ================================================================
# 15. MODEL VARIANCE / SCALE PARAMETERS
# ================================================================

variance_columns <-
  
  intersect(
    
    c(
      "sd_id__Intercept",
      "sigma"
    ),
    
    names(post)
    
  )


variance_parameters <- map_dfr(
  
  variance_columns,
  
  function(parameter_name) {
    
    draws <-
      post[[parameter_name]]
    
    tibble(
      
      parameter =
        parameter_name,
      
      estimate_N =
        mean(draws) *
        kgf_to_N,
      
      posterior_SD_N =
        sd(draws) *
        kgf_to_N,
      
      CrI_lower_N =
        as.numeric(
          quantile(
            draws,
            0.025
          )
        ) *
        kgf_to_N,
      
      CrI_upper_N =
        as.numeric(
          quantile(
            draws,
            0.975
          )
        ) *
        kgf_to_N
      
    )
  }
)


write.csv(
  variance_parameters,
  "Supplementary_Table_S2_variance_parameters.csv",
  row.names = FALSE
)


# ================================================================
# 16. CONDITIONAL FIRING-RATE SLOPES
#
# Association between 1-SD greater mean firing rate
# and six-month strength at each time x limb combination.
# ================================================================


# ------------------------------------------------
# Helper to handle possible brms interaction
# name ordering
# ------------------------------------------------

get_draw <- function(candidates) {
  
  hit <-
    candidates[
      candidates %in%
        names(post)
    ]
  
  if (length(hit) == 0) {
    
    stop(
      paste(
        "Could not find coefficient. Tried:",
        paste(
          candidates,
          collapse = ", "
        )
      )
    )
    
  }
  
  post[[hit[1]]]
  
}


# Main firing-rate coefficient

b_fr <- get_draw(
  c(
    "b_mean_fr_s"
  )
)


# Time x firing rate

b_t2_fr <- get_draw(
  c(
    "b_time2:mean_fr_s"
  )
)

b_t3_fr <- get_draw(
  c(
    "b_time3:mean_fr_s"
  )
)

b_t4_fr <- get_draw(
  c(
    "b_time4:mean_fr_s"
  )
)


# Limb x firing rate

b_opp_fr <- get_draw(
  c(
    "b_limbOpp:mean_fr_s",
    "b_mean_fr_s:limbOpp"
  )
)


# Three-way interactions

b_t2_opp_fr <- get_draw(
  c(
    "b_time2:limbOpp:mean_fr_s",
    "b_limbOpp:time2:mean_fr_s"
  )
)

b_t3_opp_fr <- get_draw(
  c(
    "b_time3:limbOpp:mean_fr_s",
    "b_limbOpp:time3:mean_fr_s"
  )
)

b_t4_opp_fr <- get_draw(
  c(
    "b_time4:limbOpp:mean_fr_s",
    "b_limbOpp:time4:mean_fr_s"
  )
)


# ------------------------------------------------
# Construct conditional slopes
# ------------------------------------------------

slopes <- list(
  
  "ACL_Pre-op" =
    b_fr,
  
  "ACL_2wk" =
    b_fr +
    b_t2_fr,
  
  "ACL_6wk" =
    b_fr +
    b_t3_fr,
  
  "ACL_3mo" =
    b_fr +
    b_t4_fr,
  
  "Contralateral_Pre-op" =
    b_fr +
    b_opp_fr,
  
  "Contralateral_2wk" =
    b_fr +
    b_t2_fr +
    b_opp_fr +
    b_t2_opp_fr,
  
  "Contralateral_6wk" =
    b_fr +
    b_t3_fr +
    b_opp_fr +
    b_t3_opp_fr,
  
  "Contralateral_3mo" =
    b_fr +
    b_t4_fr +
    b_opp_fr +
    b_t4_opp_fr
  
)


conditional_slopes <- bind_rows(
  
  lapply(
    
    names(slopes),
    
    function(condition_name) {
      
      draws <-
        slopes[[condition_name]]
      
      tibble(
        
        condition =
          condition_name,
        
        beta_N =
          mean(draws) *
          kgf_to_N,
        
        posterior_SD_N =
          sd(draws) *
          kgf_to_N,
        
        CrI_lower_N =
          as.numeric(
            quantile(
              draws,
              0.025
            )
          ) *
          kgf_to_N,
        
        CrI_upper_N =
          as.numeric(
            quantile(
              draws,
              0.975
            )
          ) *
          kgf_to_N,
        
        Pr_below_ROPE =
          mean(
            draws <
              rope_lower
          ),
        
        Pr_within_ROPE =
          mean(
            abs(draws) <=
              rope_kgf
          ),
        
        Pr_above_ROPE =
          mean(
            draws >
              rope_upper
          ),
        
        Pr_outside_ROPE =
          mean(
            abs(draws) >
              rope_kgf
          ),
        
        Pr_above_zero =
          mean(
            draws > 0
          ),
        
        Pr_below_zero =
          mean(
            draws < 0
          )
        
      )
      
    }
    
  )
  
)


write.csv(
  conditional_slopes,
  "Conditional_firing_rate_slopes.csv",
  row.names = FALSE
)


print(
  conditional_slopes
)


# ================================================================
# 17. FIGURE 1
# GROUPED POSTERIOR FIXED-EFFECT COEFFICIENTS
# ================================================================

fig1_dat <- fixed_effects %>%
  
  mutate(
    
    parameter =
      sub(
        "^b_",
        "",
        parameter
      ),
    
    label =
      case_when(
        
        # ----------------------------------------
        # Time
        # ----------------------------------------
        
        parameter == "time2" ~
          "2 weeks",
        
        parameter == "time3" ~
          "6 weeks",
        
        parameter == "time4" ~
          "3 months",
        
        
        # ----------------------------------------
        # Time x firing rate
        # ----------------------------------------
        
        parameter ==
          "time2:mean_fr_s" ~
          "2 weeks × Mean firing rate",
        
        parameter ==
          "time3:mean_fr_s" ~
          "6 weeks × Mean firing rate",
        
        parameter ==
          "time4:mean_fr_s" ~
          "3 months × Mean firing rate",
        
        
        # ----------------------------------------
        # Limb
        # ----------------------------------------
        
        parameter ==
          "limbOpp" ~
          "Contralateral limb",
        
        
        # ----------------------------------------
        # Time x limb
        # ----------------------------------------
        
        parameter %in%
          c(
            "time2:limbOpp",
            "limbOpp:time2"
          ) ~
          "Contralateral × 2 weeks",
        
        parameter %in%
          c(
            "time3:limbOpp",
            "limbOpp:time3"
          ) ~
          "Contralateral × 6 weeks",
        
        parameter %in%
          c(
            "time4:limbOpp",
            "limbOpp:time4"
          ) ~
          "Contralateral × 3 months",
        
        
        # ----------------------------------------
        # Limb x firing rate
        # ----------------------------------------
        
        parameter %in%
          c(
            "limbOpp:mean_fr_s",
            "mean_fr_s:limbOpp"
          ) ~
          "Contralateral × Mean firing rate",
        
        
        # ----------------------------------------
        # Three-way interactions
        # ----------------------------------------
        
        parameter %in%
          c(
            "time2:limbOpp:mean_fr_s",
            "limbOpp:time2:mean_fr_s"
          ) ~
          "Contralateral × 2 weeks × Mean firing rate",
        
        parameter %in%
          c(
            "time3:limbOpp:mean_fr_s",
            "limbOpp:time3:mean_fr_s"
          ) ~
          "Contralateral × 6 weeks × Mean firing rate",
        
        parameter %in%
          c(
            "time4:limbOpp:mean_fr_s",
            "limbOpp:time4:mean_fr_s"
          ) ~
          "Contralateral × 3 months × Mean firing rate",
        
        
        # ----------------------------------------
        # Firing rate
        # ----------------------------------------
        
        parameter ==
          "mean_fr_s" ~
          "Mean firing rate (1 SD)",
        
        
        # ----------------------------------------
        # Contraction intensity
        # 20% MVC is the reference category.
        # ----------------------------------------
        
        parameter ==
          "contraction_intensity_s0.5" ~
          "Contraction intensity (50%)",
        
        parameter ==
          "contraction_intensity_s0.75" ~
          "Contraction intensity (75%)",
        
        
        # ----------------------------------------
        # Covariates
        # ----------------------------------------
        
        parameter ==
          "mass_s" ~
          "Body mass (1 SD)",
        
        parameter ==
          "sexM" ~
          "Male sex",
        
        parameter ==
          "muscleVM" ~
          "Muscle: VM (vs VL)",
        
        
        TRUE ~
          parameter
        
      )
    
  ) %>%
  
  # Do not plot intercept
  filter(
    parameter != "Intercept"
  )


# ------------------------------------------------
# Order coefficient groups
# ------------------------------------------------

plot_order <- c(
  
  # Time
  "2 weeks",
  "6 weeks",
  "3 months",
  
  "SPACE1",
  
  # Time x firing rate
  "2 weeks × Mean firing rate",
  "6 weeks × Mean firing rate",
  "3 months × Mean firing rate",
  
  "SPACE2",
  
  # Limb
  "Contralateral limb",
  
  "SPACE3",
  
  # Time x limb
  "Contralateral × 2 weeks",
  "Contralateral × 6 weeks",
  "Contralateral × 3 months",
  
  "SPACE4",
  
  # Limb x firing rate
  "Contralateral × Mean firing rate",
  
  "SPACE5",
  
  # Three-way interactions
  "Contralateral × 2 weeks × Mean firing rate",
  "Contralateral × 6 weeks × Mean firing rate",
  "Contralateral × 3 months × Mean firing rate",
  
  "SPACE6",
  
  # Firing rate
  "Mean firing rate (1 SD)",
  
  "SPACE7",
  
  # Intensity
  "Contraction intensity (50%)",
  "Contraction intensity (75%)",
  
  "SPACE8",
  
  # Covariates
  "Body mass (1 SD)",
  "Male sex",
  "Muscle: VM (vs VL)"
  
)


# ------------------------------------------------
# Dummy spacer rows
# ------------------------------------------------

spacer_dat <- tibble(
  
  parameter =
    NA_character_,
  
  beta_N =
    NA_real_,
  
  posterior_SD_N =
    NA_real_,
  
  CrI_lower_N =
    NA_real_,
  
  CrI_upper_N =
    NA_real_,
  
  Pr_below_ROPE =
    NA_real_,
  
  Pr_within_ROPE =
    NA_real_,
  
  Pr_above_ROPE =
    NA_real_,
  
  Pr_outside_ROPE =
    NA_real_,
  
  label =
    paste0(
      "SPACE",
      1:8
    )
  
)


fig1_plot_dat <- bind_rows(
  
  fig1_dat,
  spacer_dat
  
) %>%
  
  mutate(
    
    label =
      factor(
        label,
        levels =
          rev(plot_order)
      )
    
  )


# ------------------------------------------------
# Save figure data
# ------------------------------------------------

write.csv(
  fig1_dat,
  "Figure1_plot_data.csv",
  row.names = FALSE
)


# ------------------------------------------------
# Plot
# ------------------------------------------------

figure1 <- ggplot(
  
  fig1_plot_dat,
  
  aes(
    x = beta_N,
    y = label
  )
  
) +
  
  # ROPE
  annotate(
    "rect",
    xmin = -rope_N,
    xmax = rope_N,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.15
  ) +
  
  # Null effect
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  
  # 95% CrI
  geom_errorbar(
    aes(
      xmin = CrI_lower_N,
      xmax = CrI_upper_N
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  
  # Posterior mean
  geom_point(
    size = 2.2,
    na.rm = TRUE
  ) +
  
  # Hide spacer names
  scale_y_discrete(
    drop = FALSE,
    labels = function(x) {
      ifelse(
        grepl(
          "^SPACE",
          x
        ),
        "",
        x
      )
    }
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid =
      element_blank(),
    
    axis.title.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank()
  ) +
  
  labs(
    x =
      "Posterior coefficient estimate (N)",
    y = NULL
  )


figure1


ggsave(
  "Figure1_updated_submax.png",
  figure1,
  width = 8,
  height = 7,
  dpi = 600
)


# ================================================================
# 18. FIGURE 2
# CONDITIONAL FIRING-RATE EFFECTS BY TIME AND LIMB
# ================================================================

ce <- conditional_effects(
  
  fit,
  
  effects =
    "mean_fr_s",
  
  conditions =
    make_conditions(
      fit,
      vars =
        c(
          "time",
          "limb"
        )
    ),
  
  prob =
    0.95
  
)


ce_df <-
  ce$mean_fr_s


# ------------------------------------------------
# Convert:
#
# x axis:
# standardised firing rate -> Hz
#
# y axis:
# kg-force -> Newtons
# ------------------------------------------------

ce_plot <- ce_df %>%
  
  mutate(
    
    mean_fr_hz =
      mean_fr_s *
      mean_fr_sd +
      mean_fr_mu,
    
    estimate_N =
      estimate__ *
      kgf_to_N,
    
    lower_N =
      lower__ *
      kgf_to_N,
    
    upper_N =
      upper__ *
      kgf_to_N,
    
    time =
      recode_factor(
        
        as.character(time),
        
        "1" =
          "Pre-surgery",
        
        "2" =
          "2-weeks post-surgery",
        
        "3" =
          "6-weeks post-surgery",
        
        "4" =
          "3-months post-surgery"
        
      ),
    
    limb =
      recode_factor(
        
        limb,
        
        "ACL" =
          "ACL",
        
        "Opp" =
          "Contralateral"
        
      )
    
  )


write.csv(
  ce_plot,
  "Figure2_plot_data.csv",
  row.names = FALSE
)


figure2 <- ggplot(
  
  ce_plot,
  
  aes(
    
    x =
      mean_fr_hz,
    
    y =
      estimate_N,
    
    colour =
      limb,
    
    fill =
      limb,
    
    group =
      limb
    
  )
  
) +
  
  geom_ribbon(
    
    aes(
      ymin =
        lower_N,
      
      ymax =
        upper_N
    ),
    
    alpha =
      0.25,
    
    colour =
      NA
    
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  facet_wrap(
    ~time,
    nrow = 1
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid =
      element_blank(),
    
    legend.position =
      "bottom"
  ) +
  
  labs(
    
    x =
      "Mean firing rate (Hz)",
    
    y =
      "Knee extension strength (N)",
    
    colour =
      "Limb",
    
    fill =
      "Limb"
    
  ) +
  
  scale_colour_viridis_d(
    end = 0.6,
    option = "D"
  ) +
  
  scale_fill_viridis_d(
    end = 0.6,
    option = "D"
  )


figure2


ggsave(
  "Figure2_updated_submax.png",
  figure2,
  width = 7,
  height = 2.5,
  dpi = 600
)
