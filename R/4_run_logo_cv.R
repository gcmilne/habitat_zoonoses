#------------------------------------------------------------------------------#
# 4_run_logo_cv.R
#
# Purpose: checks that the within-species coefficient is not driven by any
#          single study, pathogen, realm or pathogen family.
#
# Inputs:  data/df_modNDVI2.csv, via prepare_data() in R/1_prepare_data.R
#
# Methods: 1. fit the final model to the full data (reference coefficient)
#          2. leave-one-group-out (LOGO): for each study, pathogen and realm,
#             drop that level, refit the final model and record the
#             within-species coefficient and 95% credible interval
#          3. pathogen-family stratification: refit the final model
#             separately for each family (Hantaviridae, Arenaviridae); two
#             families are too few for a family random effect
#
# Outputs: outputs/reference_coefficient.rds (full-data coefficient)
#          outputs/logo_cv/logo_cv.csv        (one row per held-out level)
#          outputs/logo_cv/logo_summary.csv   (range and sign changes per grouping)
#          outputs/logo_cv/family_strata.csv  (coefficient per family)
#          outputs/logo_cv/fig_logo_stability.png / .pdf
#          outputs/logo_cv/fig_family_strata.png / .pdf
#
# Usage:   source("R/4_run_logo_cv.R")   # runtime ~35 min (single-threaded INLA)
#------------------------------------------------------------------------------#

# 1. LOAD PACKAGES ----
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
})

# 2. SETTINGS AND DATA ----

## a. Load data ----
if (!exists("prepare_data")) source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

## b. Output directory ----
dir_lg <- file.path("outputs", "logo_cv")
dir.create(dir_lg, showWarnings = FALSE, recursive = TRUE)

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
# Returns the mean and 95% interval of prob_occur_within; retries failed fits
# and returns NAs if all tries fail. Ntrials is the bare column name so INLA
# finds it inside `data` when called from a function.
# num.threads = "1:1", internal.opt = FALSE and reordering = "metis" make each
# fit identical on every run; otherwise INLA picks some settings by timing them.
fit_beta <- function(data, tries = 3) {
  for (t in seq_len(tries)) {
    fit <- tryCatch(
      inla(final_formula, family = "betabinomial", Ntrials = number_tested,
           data = data, control.fixed = prior_fixed,
           control.compute = list(internal.opt = FALSE),
           control.inla = list(strategy = "adaptive", reordering = "metis"),
           num.threads = "1:1", verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) {
      r <- fit$summary.fixed["prob_occur_within", ]
      return(c(mean = r[["mean"]], q025 = r[["0.025quant"]], q975 = r[["0.975quant"]]))
    }
  }
  c(mean = NA, q025 = NA, q975 = NA)
}

# 4. REFERENCE COEFFICIENT ----
ref <- fit_beta(df)
message(sprintf("Full-data within coefficient: %.4f [%.4f, %.4f]",
                ref["mean"], ref["q025"], ref["q975"]))
saveRDS(ref, file.path("outputs", "reference_coefficient.rds"))

# 5. LEAVE-ONE-GROUP-OUT ----

## a. Refit with each level held out ----
groupings <- c(study_id_f = "Study", pathogen_id_f = "Pathogen", REALM_f = "Realm")
logo_rows <- list()
for (g in names(groupings)) {
  levs <- levels(df[[g]])
  message(sprintf("LOGO over %s (%d levels) ...", groupings[[g]], length(levs)))
  for (lv in levs) {
    b <- fit_beta(df[df[[g]] != lv, ])
    logo_rows[[paste(g, lv)]] <- data.frame(
      grouping = unname(groupings[[g]]), held_out = lv,
      beta_within = b["mean"], q025 = b["q025"], q975 = b["q975"],
      stringsAsFactors = FALSE)
  }
  # Save after each grouping so progress is kept if the run stops
  write.csv(do.call(rbind, logo_rows), file.path(dir_lg, "logo_cv.csv"), row.names = FALSE)
}
logo <- do.call(rbind, logo_rows); rownames(logo) <- NULL
write.csv(logo, file.path(dir_lg, "logo_cv.csv"), row.names = FALSE)

## b. Summarise per grouping ----
logo_summary <- logo %>%
  group_by(grouping) %>%
  summarise(n_levels = n(),
            beta_min = min(beta_within, na.rm = TRUE),
            beta_max = max(beta_within, na.rm = TRUE),
            n_sign_flip = sum(beta_within > 0, na.rm = TRUE),                # fits with a positive estimate
            n_ci_excludes_0 = sum(sign(q025) == sign(q975), na.rm = TRUE),   # 95% interval excludes zero
            .groups = "drop")
write.csv(logo_summary, file.path(dir_lg, "logo_summary.csv"), row.names = FALSE)

# 6. PATHOGEN-FAMILY STRATIFICATION ----
fam_rows <- lapply(sort(unique(df$pathogen_family)), function(fam) {
  sub <- df[df$pathogen_family == fam, ]
  b <- fit_beta(sub)
  data.frame(family = fam, n_records = nrow(sub),
             beta_within = b["mean"], q025 = b["q025"], q975 = b["q975"])
})
fam <- do.call(rbind, fam_rows); rownames(fam) <- NULL
write.csv(fam, file.path(dir_lg, "family_strata.csv"), row.names = FALSE)

# 7. FIGURES ----

## a. Plot theme ----
col_ref <- "#B2182B"; col_pt <- "#2166AC"
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

## b. Leave-one-group-out stability ----
p_logo <- logo %>%
  group_by(grouping) %>% mutate(ord = rank(beta_within, ties.method = "first")) %>% ungroup() %>%   # rank for y-axis order
  ggplot(aes(beta_within, ord)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref["mean"], colour = col_ref, linewidth = 0.6) +
  geom_point(size = 1, colour = col_pt, alpha = 0.7) +
  facet_wrap(~grouping, scales = "free_y") +
  labs(x = "Within-species coefficient (one group held out)", y = "Held-out level (ranked)",
       title = "Leave-one-group-out stability of the within-species effect",
       subtitle = "Red line = full-data estimate; dashed = zero") +
  theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(dir_lg, "fig_logo_stability.png"), p_logo, width = 9, height = 4,
       dpi = 320, bg = "white")
ggsave(file.path(dir_lg, "fig_logo_stability.pdf"), p_logo, width = 9, height = 4)

## c. Pathogen-family stratification ----
p_fam <- ggplot(fam, aes(beta_within, family)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref["mean"], colour = col_ref, linewidth = 0.6) +
  geom_errorbar(aes(xmin = q025, xmax = q975), width = 0.1, colour = col_pt) +
  geom_point(size = 2.6, colour = col_pt) +
  labs(x = "Within-species coefficient (family-stratified)", y = NULL,
       title = "Pathogen-family stratification",
       subtitle = sprintf("Full-data estimate (red) = %.3f", ref["mean"])) +
  theme_pub
ggsave(file.path(dir_lg, "fig_family_strata.png"), p_fam, width = 7.5, height = 3.2,
       dpi = 320, bg = "white")
ggsave(file.path(dir_lg, "fig_family_strata.pdf"), p_fam, width = 7.5, height = 3.2)

# 8. PRINT SUMMARY ----
cat("\n==== LOGO summary ====\n"); print(logo_summary, row.names = FALSE)
cat("\n==== Family stratification ====\n"); print(fam, row.names = FALSE)
message("\nDONE. LOGO outputs under ", normalizePath(dir_lg))
