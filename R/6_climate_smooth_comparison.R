# 1. Manually define the formulas to prevent duplication errors

# Your original formula using rw2
formula_rw2 <- number_positive ~ 1 +
  prob_occur_within  +          
  prob_occur_between +          
  log_tested  + # Using your new standardised version                  
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw2")

# Swapping out rw2 for a cyclic rw1 with a protective PC prior
formula_rw1 <- number_positive ~ 1 +
  prob_occur_within  +          
  prob_occur_between +          
  log_tested         +                  
  f(host_species_int, model = "iid", hyper = re) +
  f(pathogen_id_f,    model = "iid", hyper = re) +
  f(study_id_f,       model = "iid", hyper = re) +
  f(REALM_f,          model = "iid", hyper = re) +
  f(month,            model = "rw1", cyclic = TRUE, 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))


# 2. Fit the RW1 Model (Spelling fixed: 'optimiser')
m_rw1 <- inla(
  formula_rw1, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.inla    = list(
    strategy = "simplified.laplace", 
    # optimiser = "gsl",               # Fixed spelling
    tolerance = 1e-6
  ),
  # num.threads     = "1:1", 
  verbose         = FALSE
)

fx <- m_rw1$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
fx

# 3. Fit the RW2 Model with identical strict optimization
m_rw2 <- inla(
  formula_rw2, family = "betabinomial", Ntrials = df$number_tested, data = df,
  control.fixed   = prior_fixed,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.inla    = list(
    strategy = "simplified.laplace", 
    # optimiser = "gsl",               # Fixed spelling
    tolerance = 1e-6
  ),
  # num.threads     = "1:1", 
  verbose         = FALSE
)

# 4. Compare Fit and Penalised Information Criteria
comparison_fit <- data.frame(
  Model = c("RW1 (Cyclic)", "RW2 (Original)"),
  DIC   = c(m_rw1$dic$dic, m_rw2$dic$dic),
  WAIC  = c(m_rw1$waic$waic, m_rw2$waic$waic),
  Mean_Log_CPO = c(mean(log(m_rw1$cpo$cpo), na.rm = TRUE), 
                   mean(log(m_rw2$cpo$cpo), na.rm = TRUE))
)

print(comparison_fit)


# 1. Extract the random effect values for month
rw1_trend <- m_rw1$summary.random$month

# 2. Plot the cyclic seasonal trend
p1 <- ggplot(rw1_trend, aes(x = ID, y = mean)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), 
              alpha = 0.2, fill = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Month", 
    y = "Log-odds deviation (linear predictor scale)"
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_pub


# 1. Extract the random effect values for month
rw2_trend <- m_rw2$summary.random$month

# 2. Plot the cyclic seasonal trend
p2 <- ggplot(rw2_trend, aes(x = ID, y = mean)) +
  geom_line(color = "#1f77b4", size = 1.2) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), 
              alpha = 0.2, fill = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Month", 
    y = "Log-odds deviation (linear predictor scale)"
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_pub

wrap_plots(p1, p2, nrow=2)