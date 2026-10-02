#------------------------------------------------------------------------------#
# 5_final_model.R
#
# Purpose: fits the final model, checks it, and produces the effect estimates
#          and curves used in the results figures.
#
# Inputs:  data/df_modNDVI2.csv, via prepare_data() in R/1_prepare_data.R
#
# Methods: 1. fit the final beta-binomial INLA model (within- and
#             between-species probability of occurrence, log number tested;
#             random intercepts for host, pathogen, study and realm; rw2
#             smooth on month)
#          2. model checks: posterior predictive density check, Moran's I on
#             residuals (k = 8 nearest neighbours), random-effect estimates,
#             and sensitivity of the within-species coefficient to the
#             fixed-effect prior SD (0.5, 1, 2)
#          3. from 500 joint posterior samples, predict seroprevalence across
#             the within-species gradient, across number tested, and by
#             month, with other covariates at their means and random effects
#             set to zero
#
# Outputs: outputs/nichepos_vs_seroprev.csv                (prediction curve)
#          outputs/fullmodel_fixed_marginal_posteriors.csv (fixed-effect densities)
#          outputs/fullmodel_random_month_rw2_summary.csv  (month effect, log-odds)
#          outputs/month_effect_response_scale.csv         (month effect, prevalence)
#          outputs/samplingeffort_vs_seroprev.csv          (prevalence vs number tested)
#          Check plots are drawn but not saved (ggsave lines commented out).
#
# Usage:   source("R/5_final_model.R")   # runtime ~1 min (single-threaded INLA)
#------------------------------------------------------------------------------#

# 1. LOAD PACKAGES ----
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(knitr)
})

# 2. DATA AND PLOT THEME ----

## a. Load data ----
source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

## b. Plot theme ----
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# 3. MODEL SPECIFICATION ----

## a. Priors ----
# Weakly informative, matched to the Stan reference model
prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),  # N(0, 1.5) intercept
                    mean = 0, prec = 1 / (1^2))                        # N(0, 1) slopes
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))          # PC prior on RE SDs

# Default priors for the overdispersion and rw2 hyperparameters can be viewed with:
# inla.models()$likelihood$betabinomial$hyper
# inla.models()$latent$rw2$hyper

## b. Model formula ----
final_formula <- number_positive ~ 1 +
  prob_occur_within  +          # deviation from host-species mean
  prob_occur_between +          # host-species mean
  log_tested +                  # sampling-effort control
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

# 4. FIT MODEL ----
# config = TRUE is needed for inla.posterior.sample() below.
# num.threads = "1:1", internal.opt = FALSE and reordering = "metis" make the
# fit identical on every run; otherwise INLA picks some settings by timing them.
m_final <- inla(
  final_formula, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
                         internal.opt = FALSE),
  control.inla    = list(strategy = "adaptive", reordering = "metis"),
  num.threads     = "1:1", verbose = FALSE)

fx <- m_final$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
kable(round(fx, 3), caption = "Fixed-effect posterior summaries (log-odds scale).")

# 5. MODEL CHECKS ----

## a. Posterior predictive check ----
S <- 500
# Joint posterior samples, reused in later sections. Both seeds are needed:
# set.seed() fixes R's random numbers and seed = 42L fixes INLA's internal sampler
# (a non-zero INLA seed requires single-threaded sampling, num.threads = "1:1")
set.seed(42)
ps <- inla.posterior.sample(S, m_final, seed = 42L, num.threads = "1:1")
pr <- grep("^Predictor", rownames(ps[[1]]$latent))
hyp <- grep("overdispersion", names(ps[[1]]$hyperpar), ignore.case = TRUE, value = TRUE)[1]
nt <- df$number_tested; obs_prev <- df$number_positive / nt
# Simulate one replicate dataset of prevalences from a posterior sample
draw <- function(s) {
  p <- plogis(s$latent[pr]); rho <- min(max(s$hyperpar[[hyp]], 1e-8), 1 - 1e-8)
  sab <- (1 - rho) / rho
  rbinom(length(p), nt, rbeta(length(p), p * sab, (1 - p) * sab)) / nt
}
set.seed(42)   # fixes which samples are used and the simulated counts
reps <- sapply(ps[sample(S, 60)], draw)   # 60 replicate datasets
rep_long <- data.frame(prev = as.vector(reps), draw = rep(seq_len(ncol(reps)), each = nrow(reps)))

# Replicated (blue) vs observed (red) prevalence densities
ggplot() +
  geom_line(data = rep_long, aes(prev, group = draw), stat = "density",
            colour = "#2166AC", alpha = 0.15, linewidth = .4) +
  geom_line(data = data.frame(prev = obs_prev), aes(prev), stat = "density",
            colour = "#B2182B", linewidth = .6) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Seroprevalence", y = "Density") +
  theme_pub

# ggsave("plots/posterior_predictions.png", bg = "white", dpi=600, width = 4, height = 4, units = "in")

## b. Residual spatial autocorrelation ----
resid_vec <- obs_prev - m_final$summary.fitted.values$mean[seq_len(nrow(df))]
if (requireNamespace("spdep", quietly = TRUE)) {
  coords <- as.matrix(df[, c("longitude", "latitude")])
  # Row-standardised weights from 8 nearest neighbours (great-circle distance)
  lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords, k = 8, longlat = TRUE)),
                        style = "W")
  mc <- spdep::moran.mc(resid_vec, lw, nsim = 499)
  cat(sprintf("Moran's I on residuals (kNN=8): I = %.3f, p = %.3f\n",
              mc$statistic, mc$p.value))
} else cat("Install 'spdep' for the spatial test.\n")

## c. Random effect estimates ----
re_specs <- c(host_species_int = "Host species", pathogen_id_f = "Pathogen",
              study_id_f = "Study", REALM_f = "Realm")
re_df <- do.call(rbind, lapply(names(re_specs), function(rn) {
  s <- m_final$summary.random[[rn]]
  data.frame(group = re_specs[[rn]], mean = s$mean,
             lo = s$`0.025quant`, hi = s$`0.975quant`)
})) %>% group_by(group) %>% mutate(ord = rank(mean, ties.method = "first")) %>% ungroup()

ggplot(re_df, aes(mean, ord)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.1, colour = "grey60") +
  geom_point(size = 1, colour = "#2166AC") +
  facet_wrap(~group, scales = "free", ncol = 2) +
  labs(x = "Random effect (log-odds)", y = "Ranked level") +
  theme_pub + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ggsave("plots/random_effects.png", bg = "white", dpi=600, width = 6, height = 6, units = "in")

## d. Prior sensitivity ----
# Refit with wider and narrower priors on the fixed-effect slopes
variants <- list(
  "Baseline N(0,1)"   = prior_fixed,
  "Wider N(0,2)"      = modifyList(prior_fixed, list(prec = 1 / (2^2))),
  "Narrower N(0,0.5)" = modifyList(prior_fixed, list(prec = 1 / (0.5^2))))
ps_tab <- do.call(rbind, lapply(names(variants), function(nm) {
  f <- if (nm == "Baseline N(0,1)") m_final else
    inla(final_formula, family = "betabinomial", Ntrials = df$number_tested,
         data = df, control.fixed = variants[[nm]],
         control.compute = list(internal.opt = FALSE),
         control.inla = list(strategy = "adaptive", reordering = "metis"),
         num.threads = "1:1")   # same repeatable settings as m_final
  r <- f$summary.fixed["prob_occur_within", ]
  data.frame(prior = nm, mean = r[["mean"]], q025 = r[["0.025quant"]], q975 = r[["0.975quant"]])
}))
kable(round(ps_tab[, -1], 3), row.names = FALSE,
      caption = paste("Within-species coefficient under alternative priors:",
                      paste(ps_tab$prior, collapse = ", ")))

ps_tab %>%
  mutate(prior = case_when(
    prior == "Baseline N(0,1)" ~ "1",
    prior == "Wider N(0,2)" ~ "2",
    prior == "Narrower N(0,0.5)" ~ "0.5"
  )) %>%
  ggplot(aes(x = mean, y = prior)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = q025, xmax = q975), width =.1) +
  geom_point(size = 1.5) +
  labs(x = "Coefficient effect size", y = "Prior distribution standard deviation") +
  theme_pub

# ggsave("plots/prior_width.png", bg = "white", dpi=600, width = 4, height = 4, units = "in")

# 6. WITHIN-SPECIES EFFECT ON SEROPREVALENCE SCALE ----

## a. Prediction curve ----
# Fixed-effect draws: one row per coefficient, one column per posterior sample
B <- sapply(ps, function(s) s$latent[paste0(m_final$names.fixed, ":1"), 1])
rownames(B) <- m_final$names.fixed
# Linear predictor excluding the within-species term, other covariates at their means
const <- B["(Intercept)", ] + B["prob_occur_between", ] * mean(df$prob_occur_between) +
  B["log_tested", ] * mean(df$log_tested)
bw <- B["prob_occur_within", ]
xg <- seq(quantile(df$prob_occur_within, .01), quantile(df$prob_occur_within, .99), length = 100)
eta <- outer(xg, bw) + matrix(const, length(xg), length(const), byrow = TRUE)
prev <- plogis(eta)
ribbon <- data.frame(prob_occur_stnd = xg, med_seroprev = apply(prev, 1, median),
                     lo_seroprev = apply(prev, 1, quantile, .025), hi_seroprev = apply(prev, 1, quantile, .975))

write.csv(ribbon, "outputs/nichepos_vs_seroprev.csv", row.names = FALSE)

## b. Change across the niche gradient ----
# Predicted prevalence at the 1st and 99th percentiles of prob_occur_within
prev_lo <- plogis(const + bw * quantile(df$prob_occur_within, .01))
prev_hi <- plogis(const + bw * quantile(df$prob_occur_within, .99))
cat(sprintf("Predicted prevalence: low niche %.1f%% vs high niche %.1f%%\n",
            100 * median(prev_lo), 100 * median(prev_hi)))

round(quantile(prev_lo, c(0.5, 0.025, 0.975))*100, 1)  # lowest-quality habitat
round(quantile(prev_hi, c(0.5, 0.025, 0.975))*100, 1)  # highest-quality habitat

cat(sprintf("Relative change (high/low): %.2f [%.2f, %.2f]\n",
            median(prev_hi / prev_lo), quantile(prev_hi / prev_lo, .025),
            quantile(prev_hi / prev_lo, .975)))

# 7. SAMPLING EFFORT EFFECT ----

## a. Helper functions ----
# Linear predictor excluding log_tested, other covariates at their means
const_n <- B["(Intercept)", ] +
  B["prob_occur_within", ]  * mean(df$prob_occur_within) +
  B["prob_occur_between", ] * mean(df$prob_occur_between)
bn <- B["log_tested", ]

p_at <- function(n) plogis(const_n + bn * log(n))   # prevalence per draw at n tested
summ <- function(x) c(median = median(x), lo = quantile(x, .025), hi = quantile(x, .975))

n_ref <- median(df$number_tested)

## b. Change per doubling of number tested ----
# Average change in prevalence (percentage points) across observed records
X   <- cbind(1, df$prob_occur_within, df$prob_occur_between, df$log_tested)
Bf  <- B[c("(Intercept)", "prob_occur_within", "prob_occur_between", "log_tested"), ]
eta <- X %*% Bf                                   # N_obs x 500, random effects = 0
d_pp <- 100 * colMeans(plogis(eta + matrix(bn * log(2), nrow(X), ncol(Bf), byrow = TRUE)) - plogis(eta))
summ(d_pp)

summ(exp(bn * log(2)))                 # odds ratio per doubling

mean(p_at(2 * n_ref) < p_at(n_ref))    # posterior probability of a negative effect; same as mean(bn < 0)

# 8. MARGINAL POSTERIORS ----

## a. Extraction function ----
# Density, mean, and 80% and 95% intervals for one INLA marginal
extract_marginal <- function(marginal, parameter, effect) {

  density_df <- as.data.frame(marginal)
  names(density_df) <- c("x", "density")

  mean_val <- inla.emarginal(
    fun = identity,
    marginal = marginal
  )

  q <- inla.qmarginal(
    c(0.025, 0.10, 0.90, 0.975),
    marginal = marginal
  )

  tibble(
    parameter = parameter,
    effect = effect,
    x = density_df$x,
    density = density_df$density,
    mean = mean_val,
    lower80 = q[2],
    upper80 = q[3],
    lower95 = q[1],
    upper95 = q[4]
  )
}

## b. Fixed effects ----
fixed_df <- map_dfr(
  names(m_final$marginals.fixed),
  ~ extract_marginal(
    marginal = m_final$marginals.fixed[[.x]],
    parameter = .x,
    effect = "fixed"
  )
)

write.csv(fixed_df, "outputs/fullmodel_fixed_marginal_posteriors.csv", row.names = FALSE)

## c. Month effects ----
month_effects <- m_final$summary.random$month
write.csv(month_effects, "outputs/fullmodel_random_month_rw2_summary.csv")

# 9. MONTH EFFECT ON SEROPREVALENCE SCALE ----
mi        <- grep("^month:", rownames(ps[[1]]$latent))
month_ids <- m_final$summary.random$month$ID
M         <- sapply(ps, function(s) s$latent[mi, 1])  # n_months x 500

# Fixed part with all covariates at their means
const_m <- B["(Intercept)", ] +
  B["prob_occur_within", ]  * mean(df$prob_occur_within) +
  B["prob_occur_between", ] * mean(df$prob_occur_between) +
  B["log_tested", ]         * mean(df$log_tested)

# Prevalence per month per draw; other random effects set to zero
prev_m <- plogis(M + matrix(const_m, nrow(M), ncol(M), byrow = TRUE))

month_resp <- data.frame(
  month = month_ids,
  med   = apply(prev_m, 1, median),
  lo    = apply(prev_m, 1, quantile, .025),
  hi    = apply(prev_m, 1, quantile, .975))

write.csv(month_resp, "outputs/month_effect_response_scale.csv", row.names = FALSE)

# 10. SEROPREVALENCE VS NUMBER TESTED ----
X0  <- cbind(1, df$prob_occur_within, df$prob_occur_between)
eta0 <- X0 %*% B[c("(Intercept)", "prob_occur_within", "prob_occur_between"), ]  # N_obs x 500
bn   <- B["log_tested", ]

# Log-spaced grid from the 1st to 99th percentile of number tested
n_lo <- quantile(df$number_tested, .01); n_hi <- quantile(df$number_tested, .99)
n_grid <- unique(round(exp(seq(log(n_lo), log(n_hi), length = 60))))

# For each n, mean prevalence across records per posterior draw (500 x length(n_grid))
prev_n <- sapply(n_grid, function(n)
  colMeans(plogis(eta0 + matrix(bn * log(n), nrow(eta0), ncol(eta0), byrow = TRUE))))

curve_n <- data.frame(n_tested = n_grid,
                      med = apply(prev_n, 2, median),
                      lo  = apply(prev_n, 2, quantile, .025),
                      hi  = apply(prev_n, 2, quantile, .975))

write.csv(curve_n, "outputs/samplingeffort_vs_seroprev.csv", row.names = FALSE)
