###############################################################################
## 2_run_model_comparison.R
## -----------------------------------------------------------------------------
## Justification for the final model structure. This is the SLOW script
## (k-fold refitting); run it once to (re)generate the comparison tables that
## the lab book reads. Typical runtime ~20-30 min.
##
##   source("R/1_prepare_data.R")
##   source("R/2_run_model_comparison.R")        # runs on source
##
## It answers: does adding any optional term (REALM, climate smooths, human
## footprint, a host random slope, travel time) actually improve the model?
##
## Two complementary comparisons:
##   (A) In-sample information criteria (DIC, WAIC). NOTE: for this
##       beta-binomial + rw2-smooth setup these are NUMERICALLY UNRELIABLE
##       (effective-parameter estimates blow up into the thousands). The script
##       computes them but FLAGS the unreliability -- do not select on them.
##   (B) Proper out-of-sample 10-fold cross-validation log-score (the log
##       posterior predictive density of held-out records, scored correctly for
##       the beta-binomial via posterior sampling). Models are compared to the
##       Core + REALM reference with a paired difference and its standard error
##       (elpd_diff +/- 2 SE; the loo-style rule).
##
## Outputs land in outputs/comparison/.
###############################################################################

suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(42)
K_FOLDS   <- 10     # CV folds
N_SAMPLE  <- 500    # posterior draws per fold for the predictive density
REF_MODEL <- "Core + REALM"

if (!exists("prepare_data")) source("R/1_prepare_data.R")
df_model <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)
n_obs    <- nrow(df_model)
ntrials  <- df_model$number_tested
y_obs    <- df_model$number_positive

dir_cmp <- file.path("outputs", "comparison")
dir.create(dir_cmp, showWarnings = FALSE, recursive = TRUE)

theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"))
col_ref <- "#B2182B"; col_pt <- "#2166AC"

# ---- Priors -----------------------------------------------------------------
prior_fixed  <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),
                     mean = 0, prec = 1 / (1^2))
re           <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
prior_smooth <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# ---- Term library + candidate models ---------------------------------------
core_terms <- c(
  "1", "prob_occur_within", "prob_occur_between", "log_tested",
  "f(host_species_int, model = 'iid', hyper = re)",
  "f(pathogen_id_f,    model = 'iid', hyper = re)",
  "f(study_id_f,       model = 'iid', hyper = re)",
  "f(month,            model = 'rw2')"
)
opt_terms <- list(
  REALM   = "f(REALM_f, model = 'iid', hyper = re)",
  slope   = "f(host_species_slope, prob_occur_within, model = 'iid', hyper = re)",
  travel  = "travel_z",
  footprint = "human_footprint_z",
  NDVI    = "f(NDVI_z_disc,         model = 'rw2', hyper = prior_smooth, scale.model = TRUE)",
  precip  = "f(precip_ratio_z_disc, model = 'rw2', hyper = prior_smooth, scale.model = TRUE)",
  temp    = "f(temp_anom_z_disc,    model = 'rw2', hyper = prior_smooth, scale.model = TRUE)"
)
specs <- list(
  "Core (no REALM)"     = character(0),
  "Core + REALM"        = c("REALM"),
  "+ climate"           = c("REALM", "NDVI", "precip", "temp"),
  "+ human footprint"   = c("REALM", "footprint"),
  "+ host random slope" = c("REALM", "slope"),
  "+ travel time"       = c("REALM", "travel"),
  "Full (all terms)"    = c("REALM", "slope", "travel", "footprint", "NDVI", "precip", "temp")
)

build_formula <- function(opt_keys) {
  terms <- c(core_terms, unlist(opt_terms[opt_keys], use.names = FALSE))
  f <- as.formula(paste("number_positive ~", paste(terms, collapse = " + ")))
  environment(f) <- environment(); f
}

# =============================================================================
# (A) In-sample information criteria
# =============================================================================
message("Fitting candidate models for DIC/WAIC ...")
sel_rows <- lapply(names(specs), function(nm) {
  fit <- inla(build_formula(specs[[nm]]), family = "betabinomial",
              Ntrials = ntrials, data = df_model, control.fixed = prior_fixed,
              control.compute = list(dic = TRUE, waic = TRUE),
              control.inla = list(strategy = "adaptive"), verbose = FALSE)
  fx <- fit$summary.fixed["prob_occur_within", ]
  data.frame(model = nm, has_REALM = "REALM" %in% specs[[nm]],
             DIC = fit$dic$dic, WAIC = fit$waic$waic,
             DIC_p_eff = fit$dic$p.eff, WAIC_p_eff = fit$waic$p.eff,
             beta_within = fx[["mean"]], beta_q025 = fx[["0.025quant"]],
             beta_q975 = fx[["0.975quant"]], stringsAsFactors = FALSE)
})
sel <- do.call(rbind, sel_rows)
peff_max <- max(c(sel$WAIC_p_eff, sel$DIC_p_eff))
sel$IC_reliable <- peff_max <= 600
write.csv(sel, file.path(dir_cmp, "model_selection.csv"), row.names = FALSE)
message(sprintf("DIC/WAIC done. Max effective parameters = %.0f (%s).",
                peff_max, if (peff_max > 600) "IC UNRELIABLE" else "IC ok"))

# =============================================================================
# (B) 10-fold out-of-sample log-score
# =============================================================================
logsumexp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

score_heldout <- function(fit, test_idx, y_true, n_t, S) {
  samp <- inla.posterior.sample(S, fit)
  pred_rows <- grep("^Predictor", rownames(samp[[1]]$latent))
  hyp_nm <- grep("overdispersion", names(samp[[1]]$hyperpar),
                 ignore.case = TRUE, value = TRUE)[1]
  eta <- vapply(samp, function(s) s$latent[pred_rows[test_idx]], numeric(length(test_idx)))
  if (length(test_idx) == 1) eta <- matrix(eta, nrow = 1)
  rho <- vapply(samp, function(s) s$hyperpar[[hyp_nm]], numeric(1))
  rho <- pmin(pmax(rho, 1e-8), 1 - 1e-8); sab <- (1 - rho) / rho
  vapply(seq_along(test_idx), function(j) {
    p <- plogis(eta[j, ]); a <- p * sab; b <- (1 - p) * sab
    y <- y_true[j]; n <- n_t[j]
    lpmf <- lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)
    logsumexp(lpmf) - log(S)
  }, numeric(1))
}

folds <- sample(rep_len(seq_len(K_FOLDS), n_obs))   # fixed across all models
pointwise <- list()
for (nm in names(specs)) {
  message("CV: ", nm)
  form <- build_formula(specs[[nm]])
  ll <- numeric(n_obs)
  for (k in seq_len(K_FOLDS)) {
    test_idx <- which(folds == k)
    dtrain <- df_model; dtrain$number_positive[test_idx] <- NA
    fit_k <- inla(form, family = "betabinomial", Ntrials = ntrials, data = dtrain,
                  control.fixed = prior_fixed, control.compute = list(config = TRUE),
                  control.inla = list(strategy = "adaptive"), verbose = FALSE)
    ll[test_idx] <- score_heldout(fit_k, test_idx, y_obs[test_idx], ntrials[test_idx], N_SAMPLE)
  }
  pointwise[[nm]] <- ll
}

ref_pw <- pointwise[[REF_MODEL]]
cv <- do.call(rbind, lapply(names(specs), function(nm) {
  pw <- pointwise[[nm]]; d <- pw - ref_pw
  delta <- sum(d)
  se    <- if (nm == REF_MODEL) 0 else sqrt(length(d)) * sd(d)
  verdict <- if (nm == REF_MODEL) "reference"
             else if (delta >  2 * se) "favours this model"
             else if (delta < -2 * se) "favours reference"
             else "indistinguishable"
  data.frame(model = nm, elpd = sum(pw), mean_logscore = mean(pw),
             delta_elpd_vs_ref = delta, se_diff = se, verdict = verdict,
             stringsAsFactors = FALSE)
}))
cv <- merge(cv, sel[, c("model", "beta_within", "beta_q025", "beta_q975")], by = "model")
cv <- cv[order(-cv$elpd), ]
write.csv(cv, file.path(dir_cmp, "cv_logscore.csv"), row.names = FALSE)

# =============================================================================
# Figures
# =============================================================================
sel <- read.csv("outputs/comparison/model_selection.csv")
cv <- read.csv("outputs/comparison/cv_logscore.csv")

ref_beta <- sel$beta_within[sel$model == REF_MODEL]
p_beta <- sel %>% 
  mutate(overlaps_zero = if_else(beta_q975 >= 0, "Yes", "No"),
         overlaps_zero = as.factor(overlaps_zero)) %>% 
  ggplot(aes(x = beta_within, y = reorder(model, beta_within), col = overlaps_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref_beta, colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = beta_q025, xmax = beta_q975), width = 0.1) +
  geom_point(size = 1.5) +
  labs(x = bquote(beta[within] ~ "(log-odds)"), y = NULL, col = "Overlaps zero") + 
  scale_color_manual(values = c("black", "lightgrey")) +
  theme_pub + 
  theme(legend.position = "bottom")

ggsave(file.path(dir_cmp, "fig_beta_stability.png"), p_beta, width = 4, height = 4,
       units = "in", dpi = 600, bg = "white")

cvp <- cv %>% filter(model != REF_MODEL) %>%
  mutate(model = reorder(model, delta_elpd_vs_ref),
         lo = delta_elpd_vs_ref - 2 * se_diff, hi = delta_elpd_vs_ref + 2 * se_diff)

p_cv <- cvp %>% 
  mutate(overlaps_zero = if_else(hi >= 0 & lo <=0, "Yes", "No"),
         overlaps_zero = as.factor(overlaps_zero)) %>% 
  ggplot(aes(x = delta_elpd_vs_ref, y = model, col = overlaps_zero)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.1) +
  geom_point(size = 1.5) +
  labs(x = bquote(Delta ~ "ELPD vs. main model (+/- 2 SE)"), y = NULL, col = "Overlaps zero") +
  scale_color_manual(values = c("lightgrey", "black")) +
  theme_pub + 
  theme(legend.position = "bottom")

p_cv

ggsave(file.path(dir_cmp, "fig_cv_logscore.png"), p_cv, width = 4, height = 4, units = "in",
       dpi = 600, bg = "white")

message("\nDONE. Comparison outputs under ", normalizePath(dir_cmp))
print(cv[, c("model", "elpd", "delta_elpd_vs_ref", "se_diff", "verdict")], row.names = FALSE)
