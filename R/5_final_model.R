#-------------#
# Final model #
#-------------#

# 1. Load packages ----
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

# 2. Prepare data ----
source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)
imp <- attr(df, "imputation_summary")

# 2. Set ggplot theme & colour palettes ----
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# 3. Set priors & model formula ----

# Weakly-informative priors (matched to the validated Stan reference model)
prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),  # N(0, 1.5) intercept
                    mean = 0, prec = 1 / (1^2))                        # N(0, 1) slopes
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))          # PC prior on RE SDs

# inla.models()$likelihood$betabinomial$hyper  # default prior for beta-binomial overdispersion precision hyperparameter
# inla.models()$latent$rw2$hyper  # default prior for random walk hyperparameter

# Model formula
final_formula <- number_positive ~ 1 +
  prob_occur_within  +          # focal effect (population deviation from species mean)
  prob_occur_between +          # species-level companion term
  log_tested +                  # sampling-effort control
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

# 4. Fit model ----
m_final <- inla(
  final_formula, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.inla    = list(strategy = "adaptive"), verbose = FALSE)

# Summarise
fx <- m_final$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
# kable(round(fx, 3), caption = "Fixed-effect posterior summaries (log-odds scale).")

# 5. Posterior predictive checks ----
S <- 500
ps <- inla.posterior.sample(S, m_final)
pr <- grep("^Predictor", rownames(ps[[1]]$latent))
hyp <- grep("overdispersion", names(ps[[1]]$hyperpar), ignore.case = TRUE, value = TRUE)[1]
nt <- df$number_tested; obs_prev <- df$number_positive / nt
draw <- function(s) {
  p <- plogis(s$latent[pr]); rho <- min(max(s$hyperpar[[hyp]], 1e-8), 1 - 1e-8)
  sab <- (1 - rho) / rho
  rbinom(length(p), nt, rbeta(length(p), p * sab, (1 - p) * sab)) / nt
}
reps <- sapply(ps[sample(S, 60)], draw)
rep_long <- data.frame(prev = as.vector(reps), draw = rep(seq_len(ncol(reps)), each = nrow(reps)))

# Plot
ggplot() +
  geom_line(data = rep_long, aes(prev, group = draw), stat = "density",
            colour = "#2166AC", alpha = 0.15, linewidth = .4) +
  geom_line(data = data.frame(prev = obs_prev), aes(prev), stat = "density",
            colour = "#B2182B", linewidth = .6) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Seroprevalence", y = "Density") + 
  theme_pub

# ggsave("plots/posterior_predictions.png", bg = "white", dpi=600, width = 4, height = 4, units = "in")

# 6. Residual spatial autocorrelation ----
resid_vec <- obs_prev - m_final$summary.fitted.values$mean[seq_len(nrow(df))]
if (requireNamespace("spdep", quietly = TRUE)) {
  coords <- as.matrix(df[, c("longitude", "latitude")])
  lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords, k = 8, longlat = TRUE)),
                        style = "W")
  mc <- spdep::moran.mc(resid_vec, lw, nsim = 499)
  cat(sprintf("Moran's I on residuals (kNN=8): I = %.3f, p = %.3f\n",
              mc$statistic, mc$p.value))
} else cat("Install 'spdep' for the spatial test.\n")

# 7. Random effect estimates ----
re_specs <- c(host_species_int = "Host species", pathogen_id_f = "Pathogen",
              study_id_f = "Study", REALM_f = "Realm")
re_df <- do.call(rbind, lapply(names(re_specs), function(rn) {
  s <- m_final$summary.random[[rn]]
  data.frame(group = re_specs[[rn]], mean = s$mean,
             lo = s$`0.025quant`, hi = s$`0.975quant`)
})) %>% group_by(group) %>% mutate(ord = rank(mean, ties.method = "first")) %>% ungroup()

ggplot(re_df, aes(mean, ord)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "grey60") +
  geom_point(size = 1, colour = "#2166AC") +
  facet_wrap(~group, scales = "free", ncol = 2) +
  labs(x = "Random effect (log-odds)", y = "Ranked level") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme_pub

# ggsave("plots/random_effects.png", bg = "white", dpi=600, width = 6, height = 6, units = "in")


# 8. Prior sensitivity ----
variants <- list(
  "Baseline N(0,1)"   = prior_fixed,
  "Wider N(0,2)"      = modifyList(prior_fixed, list(prec = 1 / (2^2))),
  "Narrower N(0,0.5)" = modifyList(prior_fixed, list(prec = 1 / (0.5^2))))
ps_tab <- do.call(rbind, lapply(names(variants), function(nm) {
  f <- if (nm == "Baseline N(0,1)") m_final else
    inla(final_formula, family = "betabinomial", Ntrials = df$number_tested,
         data = df, control.fixed = variants[[nm]],
         control.inla = list(strategy = "adaptive"))
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


# 9. Effect size on seroprevalence scale ----
B <- sapply(ps, function(s) s$latent[paste0(m_final$names.fixed, ":1"), 1])
rownames(B) <- m_final$names.fixed
const <- B["(Intercept)", ] + B["prob_occur_between", ] * mean(df$prob_occur_between) +
  B["log_tested", ] * mean(df$log_tested)
bw <- B["prob_occur_within", ]
xg <- seq(quantile(df$prob_occur_within, .01), quantile(df$prob_occur_within, .99), length = 100)
eta <- outer(xg, bw) + matrix(const, length(xg), length(const), byrow = TRUE)
prev <- plogis(eta)
ribbon <- data.frame(prob_occur_stnd = xg, med_seroprev = apply(prev, 1, median),
                     lo_seroprev = apply(prev, 1, quantile, .025), hi_seroprev = apply(prev, 1, quantile, .975))

write.csv(ribbon, "outputs/nichepos_vs_seroprev.csv")

# ggplot(ribbon, aes(x = prob_occur_stnd, y = med_seroprev)) +
#   geom_ribbon(aes(ymin = lo_seroprev, ymax = hi_seroprev), fill = "lightgrey", alpha = 0.6) +
#   geom_line(colour = "black", linewidth = 1.1) +
#   geom_hline(yintercept = mean(obs_prev), linetype = "dashed", colour = "#B2182B") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Within-species probability of occurrence (standardised)",
#        y = "Predicted seroprevalence") +
#   theme_pub

# Calculate change across the niche gradient
prev_lo <- plogis(const + bw * quantile(df$prob_occur_within, .01))
prev_hi <- plogis(const + bw * quantile(df$prob_occur_within, .99))
cat(sprintf("Predicted prevalence: low niche %.1f%% vs high niche %.1f%%\n",
            100 * median(prev_lo), 100 * median(prev_hi)))

round(quantile(prev_lo, c(0.5, 0.025, 0.975))*100, 1)  # lowest-quality habitat
round(quantile(prev_hi, c(0.5, 0.025, 0.975))*100, 1)  # highest-quality habitat

cat(sprintf("Relative change (high/low): %.2f [%.2f, %.2f]\n",
            median(prev_hi / prev_lo), quantile(prev_hi / prev_lo, .025),
            quantile(prev_hi / prev_lo, .975)))


# 10. Effect of sampling effort on seroprevalence scale ----
const_n <- B["(Intercept)", ] +
  B["prob_occur_within", ]  * mean(df$prob_occur_within) +
  B["prob_occur_between", ] * mean(df$prob_occur_between)
bn <- B["log_tested", ]

# Helper: prevalence for each posterior draw at a given N
p_at <- function(n) plogis(const_n + bn * log(n))
summ <- function(x) c(median = median(x), lo = quantile(x, .025), hi = quantile(x, .975))

# median number tested
n_ref <- median(df$number_tested)

# Calculate average change in prevalence for every doubling of no. tested individuals
X   <- cbind(1, df$prob_occur_within, df$prob_occur_between, df$log_tested)
Bf  <- B[c("(Intercept)", "prob_occur_within", "prob_occur_between", "log_tested"), ]
eta <- X %*% Bf                                   # N_obs x 500, random effects = 0
d_pp <- 100 * colMeans(plogis(eta + matrix(bn * log(2), nrow(X), ncol(Bf), byrow = TRUE)) - plogis(eta))
summ(d_pp)                                        # average pp change per doubling

# Odds ratio per doubling no. tested
summ(exp(bn * log(2))) 

# Posterior probability that the effect is negative
mean(p_at(2 * n_ref) < p_at(n_ref))  # same as mean(bn < 0)


# 11. Extract marginal posterior densities ----

# Function to extract one marginal distribution

extract_marginal <- function(marginal, parameter, effect) {
  
  # Convert INLA marginal to a data frame
  density_df <- as.data.frame(marginal)
  names(density_df) <- c("x", "density")
  
  # Posterior summaries
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

# Extract fixed effects
fixed_df <- map_dfr(
  names(m_final$marginals.fixed),
  ~ extract_marginal(
    marginal = m_final$marginals.fixed[[.x]],
    parameter = .x,
    effect = "fixed"
  )
)

write.csv(fixed_df, "outputs/fullmodel_fixed_marginal_posteriors.csv")

# Extract month effects
month_effects <- m_final$summary.random$month
write.csv(month_effects, "outputs/fullmodel_random_month_rw2_summary.csv")


# 12. Put month effect on the response scale ----
# Month latent field from the same 500 joint posterior samples
mi        <- grep("^month:", rownames(ps[[1]]$latent))
month_ids <- m_final$summary.random$month$ID
M         <- sapply(ps, function(s) s$latent[mi, 1])  # n_months x 500

# Fixed part with the other covariates at their means
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

# 13. Expected seroprevalence vs number tested ----
X0  <- cbind(1, df$prob_occur_within, df$prob_occur_between)
eta0 <- X0 %*% B[c("(Intercept)", "prob_occur_within", "prob_occur_between"), ]  # N_obs x 500
bn   <- B["log_tested", ]

# Grid spans 1st-99th percentile of observed n, log-spaced
n_lo <- quantile(df$number_tested, .01); n_hi <- quantile(df$number_tested, .99)
n_grid <- unique(round(exp(seq(log(n_lo), log(n_hi), length = 60))))

# For each n: average prevalence over observations, per posterior draw
prev_n <- sapply(n_grid, function(n)
  colMeans(plogis(eta0 + matrix(bn * log(n), nrow(eta0), ncol(eta0), byrow = TRUE))))
# prev_n is 500 x length(n_grid)

curve_n <- data.frame(n_tested = n_grid,
                      med = apply(prev_n, 2, median),
                      lo  = apply(prev_n, 2, quantile, .025),
                      hi  = apply(prev_n, 2, quantile, .975))

write.csv(curve_n, "outputs/samplingeffort_vs_seroprev.csv", row.names = FALSE)
