# Analysis
# February 2024
# RB

#libraries
library(mice)
library(rstan)
library(ranger)
library(brms)
library(tidyverse)
library(janitor)
library(naniar)
library(ggplot2)

# Load data
dsub <- read.csv("ModellingData.csv") %>%
  clean_names() %>%
  mutate(across(
    c(age_s,
      time,
      contraction_intensity_s,
      force,
      peak_fr,
      mean_fr),
    ~ as.numeric(as.character(.))),
    id = factor(id),
    mean_fr_s = as.numeric(scale(mean_fr, center = TRUE, scale = TRUE)),
    peak_fr_s = as.numeric(scale(peak_fr, center = TRUE, scale = TRUE)),
    time = as.factor(time))

#separate out by different contraction intensities- run separate models for each contraction intensity,
#contraction intensities aren't possible at all time points so also select out appropriate timepoints
contraction_intensity <- 0.20 
time_points <- c(1,2,3,4)

#subset data
dsub_contraction_intensity <- dsub %>%
  filter(contraction_intensity_s == contraction_intensity, 
         muscle == "VL", #change to VM or VL
         time %in% time_points)

dsub_contraction_intensity <- dsub_contraction_intensity %>%
  select(-matches("peak_fr")) %>% #remove mean or peak firing rate
  select(-force)

# Impute missing values
# Imputation done by random forests (DF)
# RF don't assume linearity and make no assumptions about interactions in the data

imp_datasets <- mice(dsub_contraction_intensity, m = 5, method = "rf", seed = 123)
imp_datasets$loggedEvents

stripplot(imp_datasets, mean_fr, #change to mean or peak
          pch = 19, xlab = "Imputation number") # look at where data is imputed to

# Bayesian model
# Settings
cores_set = 8
iter_set = 4000
chains_set = 2
seed_set = 123
warmup = 2000

#Set priors
priors <- c(
  # Intercept: centred on observed mean
  set_prior("student_t(3, 1.1, 0.4)", class = "Intercept"),
  
  # All fixed effects: weakly informative around 0
  set_prior("student_t(3, 0, 0.3)", class = "b"),
  
  # Random intercept SD for id
  set_prior("student_t(3, 0, 0.25)", class = "sd", group = "id", lb = 0),
  
  # Residual SD
  set_prior("student_t(3, 0, 0.3)", class = "sigma", lb = 0),
  
  # Student-t df
  set_prior("gamma(2, 0.1)", class = "nu")
)

# Model
fit <- brm_multiple(force_six_months_norm ~
                      time +
                      limb +
                      mean_fr+ #change to mean or peak
                      sex+
                      time:mean_fr + 
                      (1|id)
                    ,
                    family = student,
                    data = imp_datasets,
                    prior = priors,
                    warmup = warmup,
                    cores = cores_set,
                    iter = iter_set,
                    chains = chains_set,
                    seed = seed_set,
                    control = list(max_treedepth = 10), 
                    sample_prior = "yes") 

# Check Rhats
summary<- summary(fit)
round(fit$rhats, 3)

#posterior predicted values check
pp_check(fit, type = "dens_overlay", ndraws = 100) +
  coord_cartesian(xlim = c(-5, 5))

#predicted vs residuals
point_preds <- fitted(fit)[, 1]
point_errs <- residuals(fit)[, 1]
qplot(point_preds, point_errs)

# Save model
muscle = "VL" #"VM" or "VL"
model_type = "mean_FR" #"Mean_FR" or "Peak_FR"
intensity = "twenty" #"twenty", "fifty", "seventyFive"
filename <- paste(muscle, "_", model_type, "_", intensity, ".csv", sep = "")

#traceplot
traceplot(fit$fit, pars = c("b_Intercept",
                            "b_time2",
                            "b_time3",
                            "b_time4", 
                            "b_mean_fr",
                            "b_limbOpp", 
                            "b_sexM", 
                            "b_time2:mean_fr", 
                            "b_time3:mean_fr", 
                            "b_time4:mean_fr",
                            "sigma"))

# Conditional effects
conditional_effects(fit)

# extract posterior draws
post <- as_draws_df(fit)

# Probability that the slope for mean_fr is greater than 0
# compute posterior probabilities > 0
prob_table <- tibble(
  parameter = c("limbOpp", "sexM", "mean_fr"),
  prob_gt_zero = c(
    mean(post$b_limbOpp > 0),
    mean(post$b_sexM    > 0),
    mean(post$b_mean_fr    > 0)
  )
)

prob_table