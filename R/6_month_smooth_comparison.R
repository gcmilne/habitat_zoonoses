#------------------------------------------------------------------------------#
# 6_month_smooth_comparison.R
#
# Purpose: compares two smooths for the seasonal (month) effect in the final
#          model: the rw2 used in the final model and a cyclic rw1, in which
#          December joins back to January.
#
# Inputs:  data/df_modNDVI2.csv, via prepare_data() in R/1_prepare_data.R
#
# Methods: 1. fit the final model (same priors and INLA settings as
#             R/5_final_model.R) with the month smooth as rw2 and as cyclic
#             rw1 (PC prior on precision, P(SD > 0.5) = 0.01)
#          2. compare fixed effects, and model fit by DIC, WAIC and mean
#             log CPO
#          3. plot each month effect (log-odds) with its 95% credible interval
#
# Outputs: fixed-effect and fit tables printed to the console; two-panel
#          month-effect plot drawn but not saved
#
# Usage:   source("R/6_month_smooth_comparison.R")   # runtime ~40 s (single-threaded INLA)
#------------------------------------------------------------------------------#

# 1. LOAD PACKAGES ----
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(knitr)
})

# 2. DATA AND PLOT THEME ----

## a. Load data ----
if (!exists("prepare_data")) source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

## b. Plot theme ----
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# 3. MODEL SPECIFICATION ----

## a. Priors ----
# Same as the final model
prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),  # N(0, 1.5) intercept
                    mean = 0, prec = 1 / (1^2))                        # N(0, 1) slopes
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))          # PC prior on RE SDs

## b. RW2 month smooth ----
# Same as the final model
formula_rw2 <- number_positive ~ 1 +
  prob_occur_within  +
  prob_occur_between +
  log_tested  +
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

## c. Cyclic RW1 month smooth ----
formula_rw1 <- number_positive ~ 1 +
  prob_occur_within  +
  prob_occur_between +
  log_tested         +
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw1", cyclic = TRUE,   # cyclic: December links to January
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))

# 4. FIT MODELS ----
# Same INLA settings as the final model, so only the month smooth differs.
# Three settings make each fit give identical results on every run; without
# them INLA chooses some settings by timing them, which varies between runs:
#   num.threads = "1:1"    single-threaded fitting
#   internal.opt = FALSE   no timing-based choice of internal strategies
#   reordering = "metis"   fixed sparse-matrix reordering

## a. Cyclic RW1 model ----
m_rw1 <- inla(
  formula_rw1, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
                         internal.opt = FALSE),
  control.inla    = list(strategy = "adaptive", reordering = "metis"),
  num.threads     = "1:1", verbose = FALSE)

## b. RW2 model ----
m_rw2 <- inla(
  formula_rw2, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
                         internal.opt = FALSE),
  control.inla    = list(strategy = "adaptive", reordering = "metis"),
  num.threads     = "1:1", verbose = FALSE)

# 5. COMPARE MODELS ----

## a. Fixed effects ----
fx_cols <- c("mean", "sd", "0.025quant", "0.975quant")
print(kable(round(m_rw1$summary.fixed[, fx_cols], 3),
      caption = "Fixed-effect posterior summaries, cyclic RW1 month smooth (log-odds scale)."))
print(kable(round(m_rw2$summary.fixed[, fx_cols], 3),
      caption = "Fixed-effect posterior summaries, RW2 month smooth (log-odds scale)."))

## b. Model fit ----
# Lower DIC/WAIC and higher mean log CPO indicate better fit
comparison_fit <- data.frame(
  Model = c("RW1 (Cyclic)", "RW2 (Original)"),
  DIC   = c(m_rw1$dic$dic, m_rw2$dic$dic),
  WAIC  = c(m_rw1$waic$waic, m_rw2$waic$waic),
  Mean_Log_CPO = c(mean(log(m_rw1$cpo$cpo), na.rm = TRUE),
                   mean(log(m_rw2$cpo$cpo), na.rm = TRUE))
)

print(kable(comparison_fit, digits = 3, caption = "Model fit by month smooth."))

# 6. PLOT MONTH EFFECTS ----

## a. Cyclic RW1 ----
rw1_trend <- m_rw1$summary.random$month

p1 <- ggplot(rw1_trend, aes(x = ID, y = mean)) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2, fill = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Month",
    y = "Log-odds deviation (linear predictor scale)"
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_pub

## b. RW2 ----
rw2_trend <- m_rw2$summary.random$month

p2 <- ggplot(rw2_trend, aes(x = ID, y = mean)) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`),
              alpha = 0.2, fill = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Month",
    y = "Log-odds deviation (linear predictor scale)"
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_pub

## c. Combine panels ----
print(wrap_plots(p1, p2, nrow = 2))   # rw1 on top, rw2 below
