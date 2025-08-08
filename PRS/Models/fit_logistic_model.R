fit_mod <- function(data, covar, outcomes) {
 
# Define formulas
    form_prs <- function(pheno) {
    as.formula(sprintf('%s ~ %s + st.prs', pheno, paste(covar, collapse = "+")))
  }
  
  form_quintile <- function(pheno) {
    as.formula(sprintf('%s ~ %s + quintile', pheno, paste(covar, collapse = "+")))
  }
  
  # Fit models
  models_prs <- map(outcomes, ~glm(form_prs(.x), data = data, family = binomial))
  models_quintile <- map(outcomes, ~glm(form_quintile(.x), data = data, family = binomial))
  
  # Build results
  bind_rows(
    map_dfr(models_prs, broom::tidy, conf.int = TRUE, .id = "model_index"), 
    map_dfr(models_quintile, broom::tidy, conf.int = TRUE, .id = "model_index")
  ) %>%
    filter(term %in% c("st.prs", "quintile")) %>%
    mutate(
      pheno = outcomes[as.numeric(model_index)],
      OR = exp(estimate),
      CI_l = exp(conf.low),  
      CI_u = exp(conf.high)
    ) %>%
    select(-model_index) %>%
    relocate(pheno, .before = term)
}
