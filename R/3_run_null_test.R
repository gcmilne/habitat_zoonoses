#------------------------------------------------------------------------------#
# 3_run_null_test.R
#
# Purpose: tests whether the observed within-species coefficient could arise
#          by chance from the data structure, using a within-species
#          permutation null.
#
# Inputs:  data/df_modNDVI2.csv, via prepare_data() in R/1_prepare_data.R
#
# Methods: 1. fit the final model to the observed data and record the
#             posterior mean of the within-species coefficient
#          2. for each of 500 permutations, shuffle prob_occur_within within
#             each host species (species means are unchanged), refit the
#             final model and record the coefficient
#          3. compare the observed coefficient to the permuted (null)
#             distribution for one- and two-sided randomisation p-values
#          Progress is saved every 10 permutations; with RESUME = TRUE an
#          interrupted run continues from the saved coefficients.
#
# Outputs: outputs/null_test/perm_coefs.csv        (permuted coefficients)
#          outputs/null_test/null_test_summary.csv (observed value, null
#                                                   quantiles, p-values)
#          outputs/null_test/fig_null_test.png / .pdf
#
# Usage:   source("R/3_run_null_test.R")   # runtime ~1.5-2 h (single-threaded INLA)
#------------------------------------------------------------------------------#

# 1. LOAD PACKAGES ----
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
})

# 2. SETTINGS AND DATA ----

## a. Settings ----
SAVE_EVERY <- 10       # write progress to disk every this many permutations
RESUME <- TRUE         # continue from an existing perm_coefs.csv if present

## b. Load data ----
if (!exists("prepare_data")) source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

## c. Output paths ----
dir_nt <- file.path("outputs", "null_test")
dir.create(dir_nt, showWarnings = FALSE, recursive = TRUE)
coef_path <- file.path(dir_nt, "perm_coefs.csv")

# 3. MODEL SPECIFICATION ----

## a. Priors ----
prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),
                    mean = 0, prec = 1 / (1^2))
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))

## b. Final model formula ----
final_formula <- number_positive ~ 1 +
  prob_occur_within + prob_occur_between + log_tested +
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

## c. Fitting function ----
# Returns the posterior mean of prob_occur_within. Retries failed INLA fits and
# returns NA if all tries fail, so one crash does not stop the run.
# Ntrials must be the bare column name: INLA evaluates it inside `data`, and a
# vector passed from within a function is not found.
# num.threads = "1:1", internal.opt = FALSE and reordering = "metis" make each
# fit identical on every run; otherwise INLA picks some settings by timing them.
fit_coef <- function(data, tries = 3) {
  for (t in seq_len(tries)) {
    fit <- tryCatch(
      inla(final_formula, family = "betabinomial", Ntrials = number_tested,
           data = data, control.fixed = prior_fixed,
           control.compute = list(dic = FALSE, waic = FALSE, internal.opt = FALSE),
           control.inla = list(strategy = "adaptive", reordering = "metis"),
           num.threads = "1:1", verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) return(fit$summary.fixed["prob_occur_within", "mean"])
  }
  NA_real_
}

# 4. OBSERVED COEFFICIENT ----
obs_coef <- fit_coef(df)
message(sprintf("Observed within-species coefficient: %.4f", obs_coef))

# 5. SUMMARY AND FIGURE FUNCTION ----
# Writes the summary table and histogram for the permutations completed so far;
# called during the loop and at the end
write_outputs <- function(perm_coefs, obs_coef) {
  perm_coefs <- perm_coefs[is.finite(perm_coefs)]            # drop failed fits
  summary_tab <- data.frame(
    observed_coef = obs_coef,
    n_perm        = length(perm_coefs),
    null_mean     = mean(perm_coefs),
    null_sd       = sd(perm_coefs),
    null_q025     = quantile(perm_coefs, 0.025),
    null_q975     = quantile(perm_coefs, 0.975),
    p_one_sided   = mean(perm_coefs <= obs_coef),          # as or more negative
    p_two_sided   = mean(abs(perm_coefs) >= abs(obs_coef)) # as or more extreme
  )
  write.csv(summary_tab, file.path(dir_nt, "null_test_summary.csv"), row.names = FALSE)
  p <- ggplot(data.frame(perm = perm_coefs), aes(perm)) +
    geom_histogram(bins = 30, fill = "#92C5DE", colour = "white") +
    geom_vline(xintercept = obs_coef, colour = "#B2182B", linewidth = 1.1) +
    geom_vline(xintercept = mean(perm_coefs), colour = "grey30", linetype = "dashed") +
    labs(x = "Within-species coefficient under the null (permuted)", y = "Count",
         title = "Within-species permutation null distribution",
         subtitle = sprintf("Observed = %.3f (red); null mean = %.3f; one-sided p = %.3f (n = %d perms)",
                            obs_coef, mean(perm_coefs), summary_tab$p_one_sided, length(perm_coefs))) +
    theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(dir_nt, "fig_null_test.png"), p, width = 7.5, height = 4.8,
         dpi = 320, bg = "white")
  ggsave(file.path(dir_nt, "fig_null_test.pdf"), p, width = 7.5, height = 4.8)
  summary_tab
}

# 6. PERMUTATION TEST ----

## a. Resume from saved progress ----
perm_coefs <- numeric(0)
if (RESUME && file.exists(coef_path)) {
  perm_coefs <- read.csv(coef_path)$perm_coef
  message(sprintf("Resuming: %d permutations already on disk.", length(perm_coefs)))
}

## b. Permutation loop ----
start <- length(perm_coefs) + 1
if (start <= N_PERM) {
  for (i in start:N_PERM) {
    set.seed(42+i) # repeatable but random for each i
    df_perm <- df %>%
      group_by(host_species) %>%
      mutate(prob_occur_within = sample(prob_occur_within)) %>%  # shuffle within species
      ungroup()
    perm_coefs[i] <- fit_coef(df_perm)
    if (i %% SAVE_EVERY == 0 || i == N_PERM) {
      write.csv(data.frame(perm_coef = perm_coefs[seq_len(i)]), coef_path, row.names = FALSE)
      if (i >= 2) write_outputs(perm_coefs[seq_len(i)], obs_coef)
      cat(sprintf("  permutation %d / %d (null mean so far %.4f)\n",
                  i, N_PERM, mean(perm_coefs[seq_len(i)])))
    }
  }
}

# 7. FINAL SUMMARY ----
summary_tab <- write_outputs(perm_coefs, obs_coef)
cat("\n==== Permutation null test ====\n"); print(summary_tab, row.names = FALSE)
message("\nDONE. Null-test outputs under ", normalizePath(dir_nt))
