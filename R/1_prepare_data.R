#------------------------------------------------------------------------------#
# 1_prepare_data.R
#
# Purpose: defines prepare_data(), which builds the modelling data frame used
#          by all later scripts. Sourcing this file only defines the function.
#
# Inputs:  data/df_modNDVI2.csv (serosurvey records with niche, climate and
#          site covariates)
#
# Methods: 1. filter records to complete, valid tests
#          2. impute missing nuisance covariates with the realm mean (global
#             mean if a realm has no data) and flag imputed values
#          3. split probability of occurrence into within-species (record
#             minus species mean) and between-species (species mean)
#             components, then standardise all predictors
#          4. create factor indices and month for INLA
#
# Outputs: data frame returned by prepare_data(), with the imputation summary
#          attached as attr(df, "imputation_summary")
#
# Usage:   source("R/1_prepare_data.R")
#          df <- prepare_data("data/df_modNDVI2.csv")
#------------------------------------------------------------------------------#

# 1. LOAD PACKAGES ----
suppressPackageStartupMessages({
  library(dplyr)
})

# 2. HELPER FUNCTIONS ----

## a. Realm-mean imputation ----
# Replace NAs with the group mean; use the global mean where a group is all NA
.impute_realm_mean <- function(x, grp) {
  gm <- ave(x, grp, FUN = function(v) mean(v, na.rm = TRUE))
  gl <- mean(x, na.rm = TRUE)
  filled <- ifelse(is.na(gm), gl, gm)
  ifelse(is.na(x), filled, x)
}

# 3. PREPARE DATA FUNCTION ----
prepare_data <- function(path = "data/df_modNDVI2.csv", verbose = TRUE) {

  ## a. Load and filter records ----
  raw <- read.csv(path) %>%
    filter(number_tested >= 1,
           number_positive <= number_tested,
           !is.na(host_species),
           !is.na(pathogen_family),
           !is.na(longitude), !is.na(latitude),
           !is.na(prob_occur_orderNorm),       # habitat suitability available
           !is.na(assay_group))

  ## b. Impute nuisance covariates ----
  # human_footprint (~18% missing) is used in the model comparison only
  impute_vars <- c("NDVI", "precip_ratio_3mo", "temp_anom_3mo",
                   "travel_time_health", "human_footprint")
  impute_vars <- impute_vars[impute_vars %in% names(raw)]
  imp_summary <- data.frame(
    variable    = impute_vars,
    n_missing   = vapply(impute_vars, function(v) sum(is.na(raw[[v]])), integer(1)),
    pct_missing = round(100 * vapply(impute_vars,
                        function(v) mean(is.na(raw[[v]])), numeric(1)), 2))
  for (v in impute_vars) {
    raw[[paste0(v, "_imputed")]] <- is.na(raw[[v]])   # flag before filling
    raw[[v]] <- .impute_realm_mean(raw[[v]], raw$REALM)
  }

  ## c. Within/between decomposition and standardisation ----
  df <- raw %>%
    group_by(host_species) %>%
    mutate(prob_occur_within  = prob_occur_orderNorm - mean(prob_occur_orderNorm, na.rm = TRUE),
           prob_occur_between = mean(prob_occur_orderNorm, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      prob_occur_within  = as.numeric(scale(prob_occur_within)),
      prob_occur_between = as.numeric(scale(prob_occur_between)),
      NDVI_z         = as.numeric(scale(NDVI / 10000)),
      precip_ratio_z = as.numeric(scale(precip_ratio_3mo)),
      temp_anom_z    = as.numeric(scale(temp_anom_3mo)),
      travel_z       = as.numeric(scale(log1p(travel_time_health))),
      human_footprint_z = as.numeric(scale(human_footprint)),
      log_tested     = log(number_tested),   # sampling-effort control
      # Pathogen species where known, otherwise family
      pathogen_id = ifelse(is.na(pathogen_species_cleaned),
                           paste0("fam_", pathogen_family), pathogen_species_cleaned)
    )

  ## d. Factor indices and month ----
  df <- df %>%
    mutate(
      host_species_int   = factor(host_species),   # random intercept index
      host_species_slope = factor(host_species),   # random slope index (comparison only)
      pathogen_id_f      = factor(pathogen_id),
      study_id_f         = factor(study_id),
      REALM_f            = factor(REALM)
    )
  df$month <- as.integer(format(as.Date(df$midpoint_date), "%m"))

  ## e. Discretised climate covariates ----
  # Binned versions for rw2 smooths in the model comparison; needs INLA
  if (requireNamespace("INLA", quietly = TRUE)) {
    df$NDVI_z_disc         <- INLA::inla.group(df$NDVI_z, n = 30)
    df$precip_ratio_z_disc <- INLA::inla.group(df$precip_ratio_z, n = 30)
    df$temp_anom_z_disc    <- INLA::inla.group(df$temp_anom_z, n = 30)
  }

  ## f. Return data and summary ----
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
