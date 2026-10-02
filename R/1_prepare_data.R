###############################################################################
## 1_prepare_data.R
## -----------------------------------------------------------------------------
## Shared data-preparation step for the niche-position / seroprevalence analysis.
## Source this file, then call prepare_data() to obtain the modelling frame.
##
##   source("R/1_prepare_data.R")
##   df <- prepare_data("data/df_modNDVI2.csv")
##
## What it does:
##   1. Loads the survey table and applies the inclusion filters.
##   2. Builds the WITHIN- and BETWEEN-species decomposition of the focal niche
##      measure (probability of occurrence). This is the heart of the design:
##        within  = a record's value minus its host-species mean  (deviation)
##        between = the host-species mean                          (species level)
##      Both are then standardised (mean 0, SD 1).
##   3. Imputes the few missing values in the *nuisance* covariates
##      (climate, travel time, human footprint) with the biogeographic-realm
##      mean (global-mean fallback), and keeps a "*_imputed" flag for each.
##      The focal measure, the response and all grouping variables are complete.
##   4. Creates the factor/index columns and the month variable that INLA needs.
##
## Only base R + dplyr are required here.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

## Impute NAs with the group (realm) mean, falling back to the global mean.
.impute_realm_mean <- function(x, grp) {
  gm <- ave(x, grp, FUN = function(v) mean(v, na.rm = TRUE))
  gl <- mean(x, na.rm = TRUE)
  filled <- ifelse(is.na(gm), gl, gm)
  ifelse(is.na(x), filled, x)
}

prepare_data <- function(path = "data/df_modNDVI2.csv", verbose = TRUE) {

  ## ---- 1. Load + inclusion filters ----------------------------------------
  raw <- read.csv(path) %>%
    filter(number_tested >= 1,                 # at least one animal tested
           number_positive <= number_tested,   # sane counts
           !is.na(host_species),
           !is.na(pathogen_family),
           !is.na(longitude), !is.na(latitude),
           !is.na(Suitability_orderNorm),       # niche layers available
           !is.na(assay_group))

  ## ---- 2. Realm-mean imputation of nuisance covariates --------------------
  ## These are controls only; none is the focal predictor. human_footprint is
  ## the most incomplete (~18% missing) and is used only in the model
  ## comparison, never in the final model.
  impute_vars <- c("NDVI", "precip_ratio_3mo", "temp_anom_3mo",
                   "travel_time_health", "human_footprint")
  impute_vars <- impute_vars[impute_vars %in% names(raw)]
  imp_summary <- data.frame(
    variable    = impute_vars,
    n_missing   = vapply(impute_vars, function(v) sum(is.na(raw[[v]])), integer(1)),
    pct_missing = round(100 * vapply(impute_vars,
                        function(v) mean(is.na(raw[[v]])), numeric(1)), 2))
  for (v in impute_vars) {
    raw[[paste0(v, "_imputed")]] <- is.na(raw[[v]])
    raw[[v]] <- .impute_realm_mean(raw[[v]], raw$REALM)
  }

  ## ---- 3. Within / between decomposition + standardisation ----------------
  df <- raw %>%
    group_by(host_species) %>%
    mutate(prob_occur_within  = prob_occur_orderNorm - mean(prob_occur_orderNorm, na.rm = TRUE),
           prob_occur_between = mean(prob_occur_orderNorm, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      ## focal predictor (standardised)
      prob_occur_within  = as.numeric(scale(prob_occur_within)),
      prob_occur_between = as.numeric(scale(prob_occur_between)),
      ## nuisance covariates (standardised) for the model comparison
      NDVI_z         = as.numeric(scale(NDVI / 10000)),
      precip_ratio_z = as.numeric(scale(precip_ratio_3mo)),
      temp_anom_z    = as.numeric(scale(temp_anom_3mo)),
      travel_z       = as.numeric(scale(log1p(travel_time_health))),
      human_footprint_z = as.numeric(scale(human_footprint)),
      ## sampling-effort control
      log_tested     = log(number_tested),
      ## pathogen identity: species where known, otherwise family
      pathogen_id = ifelse(is.na(pathogen_species_cleaned),
                           paste0("fam_", pathogen_family), pathogen_species_cleaned)
    )

  ## ---- 4. Factor / index columns + month ----------------------------------
  df <- df %>%
    mutate(
      host_species_int   = factor(host_species),   # random intercept index
      host_species_slope = factor(host_species),   # random slope index (comparison only)
      pathogen_id_f      = factor(pathogen_id),
      study_id_f         = factor(study_id),
      REALM_f            = factor(REALM)
    )
  df$month <- as.integer(format(as.Date(df$midpoint_date), "%m"))

  ## Discretised climate covariates for the rw2 smooths (comparison only).
  ## inla.group lives in the INLA package; only create these if it is loaded.
  if (requireNamespace("INLA", quietly = TRUE)) {
    df$NDVI_z_disc         <- INLA::inla.group(df$NDVI_z, n = 30)
    df$precip_ratio_z_disc <- INLA::inla.group(df$precip_ratio_z, n = 30)
    df$temp_anom_z_disc    <- INLA::inla.group(df$temp_anom_z, n = 30)
  }

  attr(df, "imputation_summary") <- imp_summary
  if (verbose) {
    message(sprintf("prepare_data: %d records | %d host species | %d pathogens | %d studies | %d realms",
                    nrow(df), nlevels(df$host_species_int), nlevels(df$pathogen_id_f),
                    nlevels(df$study_id_f), nlevels(df$REALM_f)))
    message("Imputation (realm-mean) summary:")
    print(imp_summary, row.names = FALSE)
  }
  df
}
