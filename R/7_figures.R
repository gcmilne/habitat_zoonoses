#--------------------------------------#
# Create multi-panel figure of results #
#--------------------------------------#

# 1. Load packages ----
pacman::p_load(dplyr, ggplot2, patchwork, ggridges)

# 2. Set ggplot theme & colour palettes ----
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

col_ref <- "#B2182B"; col_pt <- "#2166AC"

# 3. Load & prepare data & model outputs ----
source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE) # to get observed prevalence for plots
obs_prev <- mean(df$number_positive/df$number_tested)

ref  <- readRDS("outputs/reference_coefficient.rds") # reference coefficient
logo <- read.csv("outputs/logo_cv/logo_cv.csv") # Leave one group out
fam  <- read.csv("outputs/logo_cv/family_strata.csv") # family level stratification
sel  <- read.csv("outputs/comparison/model_selection.csv") # model selection
cv   <- read.csv("outputs/comparison/cv_logscore.csv") # cross-validation on model selection
null_summary <- read.csv("outputs/null_test/null_test_summary.csv") # summary table of null permutation test
null_coefs   <- read.csv("outputs/null_test/perm_coefs.csv") # null permutation test coefficients 
suit_vs_prev <- read.csv("outputs/nichepos_vs_seroprev.csv") # relationship btwn habitat suitability (standardised) and prevalence (on response scale)
fixed_df <- read.csv("outputs/fullmodel_fixed_marginal_posteriors.csv") # fixed effect marginal posterior distributions
month_resp <- read.csv("outputs/month_effect_response_scale.csv")
sampling_effort <- read.csv("outputs/samplingeffort_vs_seroprev.csv")

# 3. Plots ----

## i. LOGO ----
p_logo <- logo %>%
  group_by(grouping) %>% 
  mutate(ord = rank(beta_within, ties.method = "first")) %>% 
  ungroup() %>%
  mutate(overlaps_zero = if_else(q975 >= 0, "Yes", "No"),
         overlaps_zero = as.factor(overlaps_zero)) %>% 
  ggplot(aes(x = beta_within, y = ord, col = overlaps_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref["mean"], colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = q025, xmax = q975), width = .1, linewidth = .1) +
  geom_point(size = .2, alpha = .7) +
  facet_wrap(~grouping, scales = "free_y") +
  scale_color_manual(values = c("black", "lightgrey")) +
  labs(x = bquote(beta[within] ~ "(one group held out)"), y = "Held-out level (ranked)", col = "Overlaps zero") +
  theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ii. Family-level stratification ----
p_fam <- fam %>% 
  mutate(overlaps_zero = if_else(q975 >= 0, "Yes", "No"),
         overlaps_zero = as.factor(overlaps_zero)) %>%
  ggplot(aes(x = beta_within, y = family, col = overlaps_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref["mean"], colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = q025, xmax = q975), width = .1) +
  geom_point(size = 2.6) +
  labs(x = "Within-species coefficient (family-stratified)", y = NULL, col = "Overlaps zero") +
  scale_color_manual(values = c("black", "lightgrey")) +
  theme_pub

## iii. Model comparison ----
REF_MODEL <- "Core + REALM"
ref_beta <- sel$beta_within[sel$model == REF_MODEL]

p_beta <- sel %>% 
  mutate(overlaps_zero = if_else(beta_q975 >= 0, "Yes", "No"),
         overlaps_zero = as.factor(overlaps_zero)) %>%
  ggplot(aes(x = beta_within, y = reorder(model, beta_within), col = overlaps_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref_beta, colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = beta_q025, xmax = beta_q975), width = .1) +
  geom_point(size = 2.4) +
  labs(x = bquote(beta[within] ~ "(log-odds)"), y = NULL, col = "Overlaps zero") + 
  scale_color_manual(values = c("black", "lightgrey")) +
  theme_pub

## iV. Model comparison cross-validation ----
cvp <- cv %>% filter(model != REF_MODEL) %>%
  mutate(model = reorder(model, delta_elpd_vs_ref),
         lo = delta_elpd_vs_ref - 2 * se_diff, hi = delta_elpd_vs_ref + 2 * se_diff)

p_cv <- ggplot(cvp, aes(x = delta_elpd_vs_ref, y = model)) +
  geom_vline(xintercept = 0, colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = .1) +
  geom_point(size = 2.6) +
  labs(x = "Out-of-sample log-score difference vs Core + REALM (+/- 2 SE)", y = NULL, col = "Overlaps zero") +
  theme_pub

## v. Coefficient null distribution
null_summary <- read.csv("outputs/null_test/null_test_summary.csv") # summary table of null permutation test
null_coefs   <- read.csv("outputs/null_test/perm_coefs.csv") # null permutation test coefficients 
null_coefs <- null_coefs$perm_coef
obs_coef <- -0.1076454  # !!!! check why different between different figures !!!!

p_null <- ggplot(data.frame(perm = null_coefs), aes(perm)) +
  geom_histogram(bins = 30, fill = "darkgrey", colour = "white") +
  geom_vline(xintercept = obs_coef, colour = col_ref, linewidth = 1.1) +
  geom_vline(xintercept = mean(null_coefs), colour = "black", linetype = "dashed") +
  labs(x = bquote(beta[within] ~ "null"), y = "Count") +
  theme_pub

## vi. Relationship btwn habitat suitability & seroprevalence ----
p_suit <- suit_vs_prev %>% 
  ggplot(aes(x = prob_occur_stnd, y = med_seroprev)) +
  geom_ribbon(aes(ymin = lo_seroprev, ymax = hi_seroprev), fill = RColorBrewer::brewer.pal(3, "PuOr")[c(1,3)][1], alpha = 0.6) +
  geom_line(colour = "black", linewidth = 1.1) +
  geom_hline(yintercept = mean(obs_prev), linetype = "dashed", colour = "#B2182B") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Within-species suitability",
       y = "Prevalence") +
  theme_pub

## Fixed effect marginal posteriors ----

# Define custom parameter labels
par_labs <- c(
  "(Intercept)" = expression(alpha),
  "log_tested" = expression(beta[sampling]),
  "prob_occur_within" = expression(beta[within]),
  "prob_occur_between" = expression(beta[between])
)

# Prepare data
plot_df <- fixed_df %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    parameter = factor(parameter, levels = c("log_tested", "prob_occur_between", "prob_occur_within")),
    param_group = case_when(
      parameter == "log_tested" ~ "Sampling bias",
      parameter %in% c("prob_occur_within", "prob_occur_between") ~ "Rodent ecology"),
    param_group = factor(param_group, levels = c("Sampling bias", "Rodent ecology"))
  )

cols <- c("Rodent ecology" = RColorBrewer::brewer.pal(3, "PuOr")[1],
          "Sampling bias"  = RColorBrewer::brewer.pal(3, "PuOr")[3])

make_inset <- function(par, pad = 0.15, n_breaks = 4) {
  d <- filter(plot_df, parameter == par)
  rng <- c(unique(d$lower95), unique(d$upper95))
  w   <- diff(range(rng))
  xlim <- c(min(rng) - pad * w, max(rng) + pad * w)
  
  ggplot(d, aes(x = x, y = parameter)) +
    geom_density_ridges(aes(fill = param_group), scale = 1, alpha = .8, col = NA) +
    geom_segment(aes(x = lower95, xend = upper95, y = parameter, yend = parameter), linewidth = 0.4) +
    geom_segment(aes(x = lower80, xend = upper80, y = parameter, yend = parameter), linewidth = 0.8) +
    geom_point(aes(x = mean), size = 1) +
    geom_vline(xintercept = 0, linetype = 2) +
    coord_cartesian(xlim = xlim) +
    scale_x_continuous(breaks = scales::breaks_extended(n = n_breaks),
                       guide = guide_axis(check.overlap = TRUE)) +
    scale_fill_manual(values = cols, guide = "none") +
    theme_bw(base_size = 7) +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", colour = "grey30"),
          plot.margin = margin(3, 3, 3, 3))   # even on all sides
}

# Inset placement: y range is defined relative to each row's baseline,
# since ridges extend upward from the row position
x_max <- max(plot_df$x)
ix <- c(0.35, 0.85) * x_max
inset_y <- function(row, offset = -0.1, height = 1.0) c(row + offset, row + offset + height)

p_posteriors <- ggplot(plot_df, aes(x = x, y = parameter)) +
  geom_density_ridges(aes(fill = param_group), scale = 1, alpha = .8, col = NA) +
  geom_segment(aes(x = lower95, xend = upper95, y = parameter, yend = parameter),
               linewidth = 0.4, inherit.aes = FALSE) +
  geom_segment(aes(x = lower80, xend = upper80, y = parameter, yend = parameter),
               linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(x = mean, y = parameter), size = 1.5, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_y_discrete(labels = par_labs, limits = rev) +
  scale_fill_manual(values = cols) +
  coord_cartesian(clip = "off") +
  labs(x = "Marginal value", y = "", fill = "") +
  theme_pub +
  theme(legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal") +
  annotation_custom(ggplotGrob(make_inset("log_tested")),
                    xmin = ix[1], xmax = ix[2],
                    ymin = inset_y(3)[1], ymax = inset_y(3)[2]) +
  annotation_custom(ggplotGrob(make_inset("prob_occur_within")),
                    xmin = ix[1], xmax = ix[2],
                    ymin = inset_y(1)[1], ymax = inset_y(1)[2])


## Random month effect ----
p_month <- ggplot(month_resp, aes(month, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.6, fill = "lightgrey") +
  geom_line(colour = "black", linewidth = 1.1) +
  labs(x = "Month", y = "Prevalence") +
  scale_x_continuous(breaks = seq(1,12,by=3), labels = month.abb[seq(1,12,by=3)]) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_pub

## Sampling effort effect ----
n_range <- quantile(df$number_tested, c(0.01, 0.99))
df_subset <- df %>%
  filter(number_tested >= min (n_range) & number_tested <= n_range)

p_sampling <- ggplot(sampling_effort, aes(n_tested, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = RColorBrewer::brewer.pal(3, "PuOr")[3], alpha = 0.6) +
  geom_line(colour = "black", linewidth = 1.1) +
  geom_rug(data = df_subset, aes(x = number_tested), inherit.aes = FALSE,
           sides = "b", alpha = 0.05) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "N tested (log scale)",
       y = "Prevalence") +
  theme_pub

## Multi-panel ----
design <- "AAABB
           AAACC
           AAADD
           EEFFF
           EEFFF"
wrap_plots(A = p_posteriors, 
           B = p_sampling, 
           C = p_month,
           D = p_suit,
           E = p_beta + theme(legend.position = "none"), 
           F = p_logo + theme(legend.position = "none"), 
           design = design) + 
  plot_annotation(tag_levels = 'a', tag_suffix = "")

# Nat Eco Evo double column width figure
ggsave("plots/results_multipanel.pdf", device = cairo_pdf, width = 18.3, height = 18.3, units = "cm")
ggsave("plots/results_multipanel.png", bg = "white", dpi = 1000, width = 18.3, height = 18.3, units = "cm")


### Supplementary

# Null distribution figure
p_null
ggsave("plots/null_coefficient.png", bg = "white", dpi = 600, height = 4, width = 4, units = "in")


