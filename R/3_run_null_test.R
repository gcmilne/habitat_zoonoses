###############################################################################
## 3_run_null_test.R
## -----------------------------------------------------------------------------
## PRIMARY INFERENTIAL TEST: within-species permutation (randomisation) null.
##
## The within/between decomposition removes between-species confounding, but we
## still need to know whether the observed within-species coefficient could
## arise from the data structure by chance. This test breaks ONLY the
## within-species link, holding everything else fixed:
##
##   For each permutation, shuffle prob_occur_within *within each host species*
##   (so each species keeps its own set of within-values and its between-species
##   mean is untouched), then refit the final model and record the coefficient.
##
## The distribution of those permuted coefficients is the null. The observed
## coefficient is compared to it for a randomisation p-value. This is the test
## that licenses a causal-flavoured reading of the within-species effect.
##
## SLOW: ~500 refits, typically a few hours. Results are saved INCREMENTALLY to
## outputs/null_test/perm_coefs.csv so a crash never loses progress, and the
## summary + figure are written at the end. Set RESUME = TRUE to continue an
## interrupted run.
##
##   source("R/1_prepare_data.R")
##   source("R/3_run_null_test.R")
###############################################################################

suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)
N_PERM <- 500          # number of permutations
SAVE_EVERY <- 10       # write progress to disk every this many permutations
RESUME <- TRUE         # continue from an existing perm_coefs.csv if present

if (!exists("prepare_data")) source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

dir_nt <- file.path("outputs", "null_test")
dir.create(dir_nt, showWarnings = FALSE, recursive = TRUE)
coef_path <- file.path(dir_nt, "perm_coefs.csv")

prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),
                    mean = 0, prec = 1 / (1^2))
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))

## Final model (identical to the lab book) -- only the data changes per permutation
final_formula <- number_positive ~ 1 +
  prob_occur_within + prob_occur_between + log_tested +
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

# Resilient single fit: INLA occasionally crashes transiently; retry a couple of
# times, and return NA rather than aborting a multi-hour run if it persists.
# NOTE: pass Ntrials as the BARE column name `number_tested` so INLA resolves it
# inside `data`. Passing `data$number_tested` (or a local vector) fails INLA's
# non-standard evaluation when called from inside a function ("Ntrials expanded
# to NULL"), because INLA evaluates Ntrials in the data-frame environment.
fit_coef <- function(data, tries = 3) {
  for (t in seq_len(tries)) {
    fit <- tryCatch(
      inla(final_formula, family = "betabinomial", Ntrials = number_tested,
           data = data, control.fixed = prior_fixed,
           control.compute = list(dic = FALSE, waic = FALSE),
           control.inla = list(strategy = "adaptive"), verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) return(fit$summary.fixed["prob_occur_within", "mean"])
  }
  NA_real_
}

## Observed coefficient (fit once on the real data)
obs_coef <- fit_coef(df)
message(sprintf("Observed within-species coefficient: %.4f", obs_coef))

## Write summary table + figure for whatever permutations are done so far.
## Called inside the loop AND at the end, so an interrupted run still leaves a
## complete, usable summary (this is what was missing before).
write_outputs <- function(perm_coefs, obs_coef) {
  perm_coefs <- perm_coefs[is.finite(perm_coefs)]            # drop any failed fits
  summary_tab <- data.frame(
    observed_coef = obs_coef,
    n_perm        = length(perm_coefs),
    null_mean     = mean(perm_coefs),
    null_sd       = sd(perm_coefs),
    null_q025     = quantile(perm_coefs, 0.025),
    null_q975     = quantile(perm_coefs, 0.975),
    p_one_sided   = mean(perm_coefs <= obs_coef),          # as/more negative
    p_two_sided   = mean(abs(perm_coefs) >= abs(obs_coef)) # as/more extreme
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

## Resume support
perm_coefs <- numeric(0)
if (RESUME && file.exists(coef_path)) {
  perm_coefs <- read.csv(coef_path)$perm_coef
  message(sprintf("Resuming: %d permutations already on disk.", length(perm_coefs)))
}

start <- length(perm_coefs) + 1
if (start <= N_PERM) {
  for (i in start:N_PERM) {
    df_perm <- df %>%
      group_by(host_species) %>%
      mutate(prob_occur_within = sample(prob_occur_within)) %>%  # shuffle WITHIN species
      ungroup()
    perm_coefs[i] <- fit_coef(df_perm)
    if (i %% SAVE_EVERY == 0 || i == N_PERM) {
      write.csv(data.frame(perm_coef = perm_coefs[seq_len(i)]), coef_path, row.names = FALSE)
      if (i >= 2) write_outputs(perm_coefs[seq_len(i)], obs_coef)   # incremental summary
      cat(sprintf("  permutation %d / %d (null mean so far %.4f)\n",
                  i, N_PERM, mean(perm_coefs[seq_len(i)])))
    }
  }
}

summary_tab <- write_outputs(perm_coefs, obs_coef)
cat("\n==== Permutation null test ====\n"); print(summary_tab, row.names = FALSE)
message("\nDONE. Null-test outputs under ", normalizePath(dir_nt))
