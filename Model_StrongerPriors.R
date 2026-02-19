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
library(purrr)

# Set WD
here::here()

# Load data
dsub <- read.csv("ModellingData.csv") %>%
  filter(!(time %in% c("5"))) %>%
  clean_names() %>%
  mutate(
    across(
      c(age_s, time, contraction_intensity_s, force, peak_fr, mean_fr, mass),
      ~ as.numeric(as.character(.))
    ),
    id = factor(id),
    mean_fr_s = as.numeric(scale(mean_fr)),
    mass_s    = as.numeric(scale(mass)),
    time = as.factor(time),
    contraction_intensity_s = as.factor(contraction_intensity_s)
  ) %>%
  select(-any_of(c("peak_fr", "mean_fr", "force_six_months_norm")))

#### EDA
# force_six_months_norm ~ time + limb + mean_fr + sex + time:mean_fr +  (1|id) # Original model

dsub %>%
  ggplot(aes(x = mean_fr_s,
             y = force_six_months,
             colour = sex))+
  geom_point(aes(colour = sex,
                 shape = as.factor(contraction_intensity_s)))+
  facet_grid(muscle+limb~time)+
  stat_smooth(method = 'lm', se = FALSE)

dsub %>%
  ggplot(aes(x = mean_fr_s,
             y = force_six_months))+
  geom_point(aes(colour = contraction_intensity_s))+
  facet_grid(contraction_intensity_s~time)+
  stat_smooth(method = 'lm', se = FALSE) # nothing doing here

# Check some individuals
dsub %>%
  filter(id == 'P004' & contraction_intensity_s == '0.2') %>%
  ggplot(aes(x = mean_fr_s,
             y = force_six_months))+
  geom_point(aes(colour = muscle,
                 shape = limb))+
  facet_wrap(~time)

# Relationship over time for a single intensity
dsub %>%
  filter(contraction_intensity_s == '0.2') %>%
  ggplot(aes(x = mean_fr_s,
             y = force_six_months))+
  geom_point(aes(colour = limb))+
  facet_wrap(~time)

# Distributions
hist(dsub$force_six_months)
hist(dsub$mean_fr_s)

# Check need for random slope term for time
dsub %>%
  filter(id %in% c('P001','P002','P003','P004','P005','P006','P007','P008','P009')) %>%
  ggplot(aes(x = mean_fr_s,
             y = force_six_months))+
  geom_point(aes(colour = limb))+
  facet_wrap(id~time)



# Look at missing data, and over time
library(naniar)
vis_miss(dsub)
gg_miss_fct(x = dsub, fct = time)
gg_miss_upset(dsub, 
              nsets = 10,
              nintersects = NA)


#### Bayesian models
# Settings
cores_set = 8
iter_set = 4000
chains_set = 4
seed_set = 123
warmup = 2000

# Model
mod_formula <- force_six_months ~ time + limb + muscle + sex + contraction_intensity_s + mass_s + limb*time*mean_fr_s + (1|id) #
# Random slope for time not needed
# There's no limb by mean_fr_s effect, so left out
# sex by mean_fr_s was 91% chance of being positive, but wondering if this is just noise, RB?

# Notes:
# Checked limb by muscle; not needed.

fit0 <- brm(formula = mod_formula,
           #family = student,
           prior = prior(normal(0, 10), class = b),
           data = dsub,
           warmup = warmup,
           cores = cores_set,
           iter = iter_set,
           chains = chains_set,
           seed = seed_set,
           control = list(max_treedepth = 20))

# Checks
pp_check(fit0, ndraws = 100)
pp_check(fit0, type = "error_hist", ndraws = 20)
pp_check(fit0, type = "error_scatter_avg")+
  coord_flip()

# Summary
summary(fit0)
conditional_effects(fit0)

bayes_R2(fit0)
loo_R2(fit0) # similar to above

# Posterior probs
post <- as_draws_df(fit0)
mean(post$b_mean_fr_s < 0)
mean(post$`b_time2:mean_fr_s` > 0)
mean(post$`b_time3:mean_fr_s` > 0)
mean(post$`b_time4:mean_fr_s` > 0)

mean(post$`b_limbOpp:mean_fr_s` > 0)



# Now do with imputed data set for missings
# Impute data
imp_datasets <- mice(dsub,
                     m = 5,
                     method = "cart",
                     seed = 123, 
                     predictorMatrix = (1 - diag(ncol(dsub))))

# Set y to not be imputed
pred <- imp_datasets$predictorMatrix
meth <- imp_datasets$method
meth["force_six_months"] <- ""

imp_datasets <- mice(dsub,
                     m = 5,
                     method = meth,
                     seed = 123, 
                     predictorMatrix = pred)

imp_datasets$loggedEvents

# Checks
stripplot(imp_datasets, force, pch = 19, xlab = "Imputation number") # checks
stripplot(imp_datasets, mean_fr_s, pch = 19, xlab = "Imputation number") # checks

#Set priors
priors <- c(
  set_prior("student_t(3, 0, 0.2)", class = "b"),
  set_prior("student_t(3, 0, 0.5)", class = "Intercept"),  # often keep intercept a bit wider
  set_prior("student_t(3, 0, 0.3)", class = "sd", group = "id", lb = 0),
  set_prior("student_t(3, 0, 0.3)", class = "sigma", lb = 0),
  set_prior("gamma(2, 0.1)", class = "nu")
)

# Now fit model to imputed data set
fit <- brm_multiple(formula = mod_formula,
                    #family = student,
                    data = imp_datasets,
                    prior = priors,
                    warmup = warmup,
                    #combine = FALSE,
                    cores = cores_set,
                    iter = iter_set,
                    chains = chains_set,
                    seed = seed_set,
                    control = list(max_treedepth = 20))

# Save model
saveRDS(fit, "model.rds")

# Checks
pp_check(fit, ndraws = 100)
pp_check(fit, type = "error_hist", ndraws = 11)
pp_check(fit, type = "error_scatter_avg")+
  coord_flip()

# Convergence checks
# These are all < 1.01; have to check individual models. Have to use "combine = FALSE" in the model argument to do the checks
# length(fit)  # Should be 5
# rhat(fit[[1]]); max(rhat(fit[[1]]))
# rhat(fit[[2]]); max(rhat(fit[[2]]))
# rhat(fit[[3]]); max(rhat(fit[[3]]))
# rhat(fit[[4]]); max(rhat(fit[[4]]))
# rhat(fit[[5]]); max(rhat(fit[[5]]))

# Posterior probs
post <- as_draws_df(fit)
mean(post$b_mean_fr_s < 0)
mean(post$`b_time2:mean_fr_s` > 0)
mean(post$`b_time3:mean_fr_s` > 0)
mean(post$`b_time4:mean_fr_s` > 0)

mean(post$`b_limbOpp:mean_fr_s` > 0)

mean(post$b_muscleVM>0)

# For ROPE
rope_lower <- -15 # half a SD change
rope_upper <- 15

# Calculate proportion INSIDE ROPE (to comment on equivalence)
mean(post$`b_mean_fr_s` > rope_lower & post$`b_mean_fr_s` < rope_upper)

# Or proportion OUTSIDE ROPE (more common for decision making)
mean(post$`b_mean_fr_s` < rope_lower | post$`b_mean_fr_s` > rope_upper)
mean(post$`b_time2:mean_fr_s` < rope_lower | post$`b_time2:mean_fr_s` > rope_upper)
mean(post$`b_time3:mean_fr_s` < rope_lower | post$`b_time3:mean_fr_s` > rope_upper)
mean(post$`b_time4:mean_fr_s` < rope_lower | post$`b_time4:mean_fr_s` > rope_upper)
mean(post$`b_limbOpp` < rope_lower | post$`b_limbOpp` > rope_upper)
mean(post$`b_sexM` < rope_lower | post$`b_sexM` > rope_upper)
mean(post$`b_mass_s` < rope_lower | post$`b_mass_s` > rope_upper)


# R-squared
bayes_R2(fit)



# For interaction or multiple predictors
ce <- conditional_effects(fit, 
                          effects = "mean_fr_s",
                          conditions = make_conditions(fit, vars = c("time","limb")),
                          prob = 0.95)

# Extract data
ce_df <- ce$mean_fr_s

# Plot with limb as groups (different colors/lines)
d_trans <- read.csv("ModellingData.csv") %>%filter(!(time %in% c('5'))) %>% clean_names() # to transform firing rate

mu <- mean(d_trans$mean_fr, na.rm = T)
sigma <- sd(d_trans$mean_fr, na.rm = T)

ggplot(ce_df %>%
         mutate(time = recode_factor(time,
                                     '1' = 'Pre-surgery',
                                     '2' = '2-weeks post-surgery',
                                     '3' = '6-weeks post-surgery',
                                     '4' = '3-months post-surgery'),
                mean_fr_s = mean_fr_s*sigma+mu),
       aes(x = mean_fr_s, y = estimate__, 
                  color = limb, fill = limb, group = limb)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha = 0.25, color = NA) +
  geom_line(size = 1) +
  facet_wrap(~time, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  labs(x = "Mean firing rate (Hz)",
       y = "Knee extension strength (N)",
       color = "Limb",
       fill = "Limb")+
  scale_colour_viridis_d(end = 0.6, option = 'D')+
  scale_fill_viridis_d(end = 0.6, option = 'D')

ggsave(filename = 'figure3.png',
       width = 7,
       height = 2,
       dpi = 600)


#### Here down is you RB

# Check Rhats


summary(fit)

summary <- summary(fit)
round(fit$rhats, 3)

#posterior predicted values check
pp_check(fit, type = "dens_overlay", ndraws = 100) #+
  #coord_cartesian(xlim = c(-5, 5))

#predicted vs residuals
point_preds <- fitted(fit)[, 1]
point_errs <- residuals(fit)[, 1]
qplot(point_preds, point_errs)

# Save model
muscle = "VM" #"VM" or "VL"
model_type = "mean_FR" #"Mean_FR" or "Peak_FR"
intensity = "twenty" #"twenty", "fifty", "seventyFive"
filename <- paste(muscle, "_", model_type, "_", intensity, ".csv", sep = "")

#traceplot
traceplot(fit$fit, pars = c("b_Intercept",
                            "b_time2",
                            "b_time3",
                            "b_time4", 
                            "b_mean_fr_s",
                            "b_limbOpp", 
                            "b_sexM", 
                            "b_time2:mean_fr_s", 
                            "b_time3:mean_fr_s", 
                            "b_time4:mean_fr_s",
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
    mean(post$b_mean_fr_s    > 0)
  )
)

prob_table

# ROPE for slope (N/kg per Hz)
rope_low  <- -0.025
rope_high <-  0.025

# % of posterior draws inside ROPE for the baseline slope
p_in_rope_mean_fr <- mean(post$b_mean_fr_s >= rope_low & post$b_mean_fr_s <= rope_high)

tibble(
  parameter = "b_mean_fr (baseline)",
  p_in_rope = p_in_rope_mean_fr,
  pct_in_rope = 100 * p_in_rope_mean_fr,
  pct_outside_rope = 100 * (1 - p_in_rope_mean_fr)
)

# Identify the interaction coefficient columns
coef_cols <- names(post)[
  grepl("^b_", names(post)) &
    grepl("time", names(post)) &
    grepl("mean_fr", names(post))
]

# Summarise ROPE proportions for each interaction coefficient
rope_summary_time_mean_fr <- map_dfr(coef_cols, function(coef_name) {
  draws <- post[[coef_name]]
  p_in  <- mean(draws >= rope_low & draws <= rope_high)
  
  tibble(
    parameter = coef_name,
    p_in_rope = p_in,
    pct_in_rope = 100 * p_in,
    pct_outside_rope = 100 * (1 - p_in)
  )
})

rope_summary_time_mean_fr
