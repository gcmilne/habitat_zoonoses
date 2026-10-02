# A consistent signature of habitat suitability on rodent zoonotic pathogen prevalence

### Gregory C. Milne, Gonzalo Albaladejo-Robles, Artur Trebski, David Simons, Harry Gordon, Ana Martinez-Checa, David W. Redding

This project asks whether rodent populations living in better-quality habitat have different prevalence of infection with zoonotic viruses (Hantaviridae and Arenaviridae).

It combines serosurvey records with outputs from species distribution models. Each record is the number of animals of one host species tested and found positive at a site, and the distribution models give each record's probability of occurrence. Probability of occurrence is split into two parts:

- a **within-species** component: how a population's habitat compares with its species' average;
- a **between-species** component: the species' average.

Seroprevalence is modelled on both component with a Bayesian beta-binomial model, fitted with [INLA](https://www.r-inla.org/). The main quantity of interest is the within-species coefficient (`prob_occur_within`, written β<sub>within</sub> in the figures).

## Folder structure

```
habitat_zoonoses/
├── habitat_zoonoses.Rproj    RStudio project: open this first (sets the working directory)
├── R/                        analysis scripts, numbered in the order they are run
│   ├── 1_prepare_data.R
│   ├── 2_run_model_comparison.R
│   ├── 3_run_null_test.R
│   ├── 4_run_logo_cv.R
│   ├── 5_final_model.R
│   ├── 6_month_smooth_comparison.R
│   └── 7_figures.R
├── data/
│   └── df_modNDVI2.csv       input data: one row per serosurvey record
├── outputs/                  tables and diagnostic figures written by scripts 2-5
│   ├── comparison/           model comparison (script 2)
│   ├── null_test/            permutation null test (script 3)
│   ├── logo_cv/              leave-one-group-out and family checks (script 4)
│   ├── reference_coefficient.rds          full-data coefficient (script 4)
│   ├── nichepos_vs_seroprev.csv                 predicted prevalence vs within-species suitability (script 5)
│   ├── fullmodel_fixed_marginal_posteriors.csv  fixed-effect posterior densities (script 5)
│   ├── fullmodel_random_month_rw2_summary.csv   month effect, log-odds scale (script 5)
│   ├── month_effect_response_scale.csv          month effect, prevalence scale (script 5)
│   ├── samplingeffort_vs_seroprev.csv           predicted prevalence vs number tested (script 5)
│   ├── labbook_figs/         older lab-book figures (not produced by the current scripts)
│   └── logo_run.log          console log from a past run of script 4
└── plots/                    final figures (script 7)
    ├── results_multipanel.pdf / .png   main results figure
    ├── null_coefficient.png            supplementary figure
    └── posterior_predictions.png, random_effects.png, prior_width.png
                              older check plots from script 5 (their ggsave lines are
                              now commented out, so they are not updated)
```

All paths in the scripts are relative to the project root, so open `habitat_zoonoses.Rproj`, or set the working directory to this folder, before running anything.

## Scripts

Run the scripts in numbered order. Each one begins with a header giving its purpose, inputs, methods, outputs and approximate run time.

| # | Script | What it does | Writes to | Run time* |
|---|---|---|---|---|
| 1 | `1_prepare_data.R` | Defines `prepare_data()`, which every later script uses to build the modelling data. | — | seconds |
| 2 | `2_run_model_comparison.R` | Tests whether extra terms improve the model, to justify the final model. | `outputs/comparison/` | ~1–1.5 h |
| 3 | `3_run_null_test.R` | Permutation test: is the within-species effect stronger than chance? | `outputs/null_test/` | ~1.5–2 h |
| 4 | `4_run_logo_cv.R` | Checks that no single study, pathogen, realm or pathogen family drives the result. | `outputs/logo_cv/`, `outputs/reference_coefficient.rds` | ~35 min |
| 5 | `5_final_model.R` | Fits the final model, runs model checks and computes effect curves. | `outputs/*.csv` | ~1 min |
| 6 | `6_month_smooth_comparison.R` | Sensitivity check on the seasonal (month) smooth. | console only | ~40 s |
| 7 | `7_figures.R` | Builds the published figures from the saved outputs of scripts 2–5. | `plots/` | < 1 min |

\*Approximate times on a desktop PC, with INLA in single-threaded mode (see [Reproducibility](#reproducibility)).

### 1. `1_prepare_data.R`: data preparation

Sourcing this file only defines `prepare_data(path, verbose)`. Calling the function:

1. loads `data/df_modNDVI2.csv` and keeps records with valid test counts and complete host, pathogen, location, niche and assay information;
2. fills in missing values of the background (nuisance) covariates (NDVI, precipitation, temperature anomaly, travel time, human footprint) with the realm mean, and flags each filled value in a `*_imputed` column;
3. splits probability of occurrence into within-species and between-species parts, and standardises all predictors (mean 0, SD 1);
4. creates the grouping columns and month variable that INLA needs.

A summary of the imputation is attached to the result as `attr(df, "imputation_summary")`.

### 2. `2_run_model_comparison.R`: model comparison

Fits seven candidate models: the core model, with or without a realm random effect, plus versions adding climate smooths, human footprint, a host random slope or travel time. Each is scored by 10-fold cross-validation and compared with "Core + REALM" (difference in expected log predictive density, ±2 SE). DIC and WAIC are recorded but flagged as unreliable for this model, so they aren't used to choose between models.

### 3. `3_run_null_test.R`: permutation null test

Shuffles `prob_occur_within` *within each host species* 500 times and refits the final model each time, to build a null distribution for the within-species coefficient. The observed coefficient is compared with this distribution to give randomisation p-values. Progress is saved every 10 permutations; set `RESUME = TRUE` to continue an interrupted run.

### 4. `4_run_logo_cv.R`: leave-one-group-out checks

Refits the final model with each study (138), pathogen (42) and realm (5) left out in turn, and separately within each pathogen family. It records the within-species coefficient each time. It also saves the full-data reference coefficient that script 7 uses.

### 5. `5_final_model.R`: final model

Fits the final model:

```
number_positive ~ prob_occur_within + prob_occur_between + log_tested
                  + random intercepts for host species, pathogen, study and realm
                  + RW2 smooth over month
```

It then runs model checks (posterior predictive check, Moran's I for spatial autocorrelation in the residuals, random-effect estimates, and sensitivity to the prior). Finally, it draws 500 posterior samples to predict seroprevalence across the within-species gradient, by month, and by number tested. The check plots are drawn on screen but not saved.

### 6. `6_month_smooth_comparison.R`: month smooth comparison

Refits the final model with a cyclic RW1 month smooth, in which December links back to January, and compares it with the RW2 smooth by fixed effects, DIC, WAIC and CPO. This is a sensitivity check: it writes no files and isn't used by script 7.

### 7. `7_figures.R`: figures

Reads the outputs of scripts 2–5 and builds:

- `plots/results_multipanel.pdf` / `.png`: the six-panel main figure, at Nature Ecology & Evolution double-column width (183 mm):
  - fixed-effect posteriors;
  - sampling effort;
  - month;
  - within-species suitability;
  - model-comparison coefficients;
  - leave-one-group-out results;
- `plots/null_coefficient.png`: the permutation null distribution (supplementary figure).

## How the scripts depend on each other

```
                    1_prepare_data.R  (sourced by every script below)
                            │
   ┌────────────┬───────────┼────────────┬────────────┐
   2            3           4            5            6
comparison   null test    LOGO       final model   month smooth
   │            │           │            │         (stand-alone check)
   └────────────┴─────┬─────┴────────────┘
                      7_figures.R
```

Scripts 2–6 don't depend on each other, so you can run them in any order or separately. Script 7 needs the saved outputs of scripts 2–5. The files in `outputs/` are kept in the repository, so `7_figures.R` can be run straight away without refitting anything.

## Requirements

- **R** ≥ 4.4.
- **INLA** isn't on CRAN. Install it from the INLA repository:
  ```r
  install.packages("INLA", repos = c(getOption("repos"),
                   INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
  ```
- **CRAN packages:** `dplyr`, `tidyr`, `purrr`, `ggplot2`, `knitr`, `patchwork`, `ggridges`, `scales`, `RColorBrewer`, `pacman`, and optionally `spdep` (for the spatial autocorrelation check in script 5).

The analysis was last run with R 4.4.0, INLA 24.12.11 and ggplot2 4.0.3.

## Reproducibility

- **Model fits:** every INLA fit uses `num.threads = "1:1"`, `control.compute = list(internal.opt = FALSE)` and `control.inla = list(reordering = "metis")`. With INLA's default multi-threaded settings, fits can differ slightly between runs, because INLA chooses some internal settings by timing them. These settings make each fit identical on every run, at the cost of slower fitting (about 2–4× per fit).
- **Random numbers:** R's random numbers are seeded with `set.seed()` (script 3 uses a separate seed for each permutation). INLA's posterior sampler is seeded with `inla.posterior.sample(..., seed = 42L)`.
- **Saved outputs:** some files in `outputs/` and `plots/` were produced before these settings were added. Rerun scripts 2–5 and then 7 to regenerate them all consistently.

## Key data columns

`data/df_modNDVI2.csv` has one row per serosurvey record. The columns used in the models are:

| Column | Meaning |
|---|---|
| `number_tested`, `number_positive` | animals tested and seropositive (the response) |
| `host_species`, `host_family` | rodent host |
| `pathogen_species_cleaned`, `pathogen_family` | virus (species where known, otherwise family) |
| `study_id` | source study |
| `REALM` | biogeographic realm |
| `prob_occur_orderNorm` | probability of occurrence from the species distribution model; split into `prob_occur_within` / `prob_occur_between` by `prepare_data()` |
| `midpoint_date` | sampling date, used to get `month` |
| `longitude`, `latitude` | sampling location |
| `NDVI`, `precip_ratio_3mo`, `temp_anom_3mo`, `travel_time_health`, `human_footprint` | background covariates, used only in the model comparison (script 2) |
