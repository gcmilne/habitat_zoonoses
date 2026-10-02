###############################################################################
## 4_run_logo_cv.R
## -----------------------------------------------------------------------------
## Leave-one-group-out (LOGO) coefficient-stability checks + family stratification.
## SLOW (~30 min): refits the final model many times. Run once; the lab book
## reads the saved tables.
##
##   source("R/1_prepare_data.R")
##   source("R/4_run_logo_cv.R")
##
## (A) LOGO stability. For each grouping in {study, pathogen, realm}, drop one
##     level at a time, refit the final model on the rest, and record the
##     within-species coefficient. If the coefficient stays negative and similar
##     across every leave-one-out fit, the result is not driven by any single
##     study, pathogen or biogeographic realm.
##
## (B) Pathogen-family stratification. There are only two viral families
##     (Hantaviridae, Arenaviridae), so a family "random effect" is degenerate;
##     the informative test is to refit the model WITHIN each family and check
##     the within-species effect holds in both (equivalently, leave-one-family-out).
##
## Outputs land in outputs/logo_cv/.
###############################################################################

suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)
if (!exists("prepare_data")) source("R/1_prepare_data.R")
df <- prepare_data("data/df_modNDVI2.csv", verbose = FALSE)

dir_lg <- file.path("outputs", "logo_cv")
dir.create(dir_lg, showWarnings = FALSE, recursive = TRUE)

prior_fixed <- list(mean.intercept = 0, prec.intercept = 1 / (1.5^2),
                    mean = 0, prec = 1 / (1^2))
re <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))

final_formula <- number_positive ~ 1 +
  prob_occur_within + prob_occur_between + log_tested +
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

# Bare `number_tested` so INLA resolves Ntrials inside `data` (works in a function)
fit_beta <- function(data, tries = 3) {
  for (t in seq_len(tries)) {
    fit <- tryCatch(
      inla(final_formula, family = "betabinomial", Ntrials = number_tested,
           data = data, control.fixed = prior_fixed,
           control.inla = list(strategy = "adaptive"), verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fit)) {
      r <- fit$summary.fixed["prob_occur_within", ]
      return(c(mean = r[["mean"]], q025 = r[["0.025quant"]], q975 = r[["0.975quant"]]))
    }
  }
  c(mean = NA, q025 = NA, q975 = NA)
}

## Full-data reference coefficient
ref <- fit_beta(df)
message(sprintf("Full-data within coefficient: %.4f [%.4f, %.4f]",
                ref["mean"], ref["q025"], ref["q975"]))
saveRDS(ref, file.path("outputs", "reference_coefficient.rds"))

# =============================================================================
# (A) Leave-one-group-out over study / pathogen / realm
# =============================================================================
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
  # incremental save
  write.csv(do.call(rbind, logo_rows), file.path(dir_lg, "logo_cv.csv"), row.names = FALSE)
}
logo <- do.call(rbind, logo_rows); rownames(logo) <- NULL
write.csv(logo, file.path(dir_lg, "logo_cv.csv"), row.names = FALSE)

logo_summary <- logo %>%
  group_by(grouping) %>%
  summarise(n_levels = n(),
            beta_min = min(beta_within, na.rm = TRUE),
            beta_max = max(beta_within, na.rm = TRUE),
            n_sign_flip = sum(beta_within > 0, na.rm = TRUE),
            n_ci_excludes_0 = sum(sign(q025) == sign(q975), na.rm = TRUE),
            .groups = "drop")
write.csv(logo_summary, file.path(dir_lg, "logo_summary.csv"), row.names = FALSE)

# =============================================================================
# (B) Pathogen-family stratification
# =============================================================================
fam_rows <- lapply(sort(unique(df$pathogen_family)), function(fam) {
  sub <- df[df$pathogen_family == fam, ]
  b <- fit_beta(sub)
  data.frame(family = fam, n_records = nrow(sub),
             beta_within = b["mean"], q025 = b["q025"], q975 = b["q975"])
})
fam <- do.call(rbind, fam_rows); rownames(fam) <- NULL
write.csv(fam, file.path(dir_lg, "family_strata.csv"), row.names = FALSE)

# =============================================================================
# Figures
# =============================================================================
col_ref <- "#B2182B"; col_pt <- "#2166AC"
theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

p_logo <- logo %>%
  group_by(grouping) %>% mutate(ord = rank(beta_within, ties.method = "first")) %>% ungroup() %>%
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

p_fam <- ggplot(fam, aes(beta_within, family)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = ref["mean"], colour = col_ref, linewidth = 0.6) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.15, colour = col_pt) +
  geom_point(size = 2.6, colour = col_pt) +
  labs(x = "Within-species coefficient (family-stratified)", y = NULL,
       title = "Pathogen-family stratification",
       subtitle = sprintf("Full-data estimate (red) = %.3f", ref["mean"])) +
  theme_pub
ggsave(file.path(dir_lg, "fig_family_strata.png"), p_fam, width = 7.5, height = 3.2,
       dpi = 320, bg = "white")
ggsave(file.path(dir_lg, "fig_family_strata.pdf"), p_fam, width = 7.5, height = 3.2)

cat("\n==== LOGO summary ====\n"); print(logo_summary, row.names = FALSE)
cat("\n==== Family stratification ====\n"); print(fam, row.names = FALSE)
message("\nDONE. LOGO outputs under ", normalizePath(dir_lg))
