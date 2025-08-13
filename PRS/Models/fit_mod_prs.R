# Given a set of covariables, a PRS, and a set of binary phenotype (outcomes), this function
# fits a logistic model for each of the outcomes and returns
# association statistics for the PRS:
# 1. As a continuous variable
# 2. Binary comparison between risk categories.

# Parameters:
# - data: dataframe with data
# - covar: vector of covariables 
# - outcomes: vector of binary outcomes
# - prs.variable: name of the variable with the PRS
# - compute.quintiles: TRUE by default. If TRUE, computes quintiles from the PRS variable.
# - prs.quintile.variable: NULL by default. If available in your data, you may provide the quintile variable using this argument.
# - prs.ntile.comparison.variable: Custom comparison variable (binary, 1/0).

# By default, the function computes PRS quintiles and then compares the 5th vs 4-6th (top 20% risk vs 40-60% -average- risk).
# If your data already has a quintile variable, `compute.quintiles` may be set to FALSE and the variable may be provided to `prs.quintile.variable`.
# Lastly, the function allows for a custom comparison of risk categories, and in this case, this must be provided to `prs.ntile.comparison.variable`.



fit_mod_prs <- function(data, covar, outcomes, prs.variable, compute.quintiles = TRUE, prs.quintile.variable = NULL, prs.ntile.comparison.variable = NULL) {
  
  if (is.null(prs.quintile.variable) && compute.quintiles) {
    data <- data %>%
      mutate(
        quintile = ntile(.data[[prs.variable]], 5),
        comparison = ifelse(quintile == 5, 1, ifelse(quintile == 3, 0, NA))
        #comparison = ifelse(quintile == 5, 1, 0)
        # Comparamos o top 20% co risco medio (40-60%)
      )
  } else if (!is.null(prs.quintile.variable)) {
    data <- data %>%
      mutate(
        quintile = .data[[prs.quintile.variable]],
        #comparison = ifelse(quintile == 5, 1, 0)
        comparison = ifelse(quintile == 5, 1, ifelse(quintile == 3, 0, NA))
      )
  } else if (!is.null(prs.ntile.comparison.variable)) {
   data$comparison <- data[[prs.ntile.comparison.variable]] }

 if (is.null(prs.quintile.variable) && !compute.quintiles && is.null(prs.ntile.comparison.variable)) { 
  stop("There is no variable for comparison between risk categories. Either provide one or allow computation of quintiles (compute.quintiles=T)") }
  
  data <- data %>%
    mutate(st.score = .data[[prs.variable]])
  
  form_prs <- function(pheno) {as.formula(sprintf('%s ~ %s + st.score', pheno, paste(covar, collapse = "+")))}
  
  form_quintile <- function(pheno) {as.formula(sprintf('%s ~ %s + comparison', pheno, paste(covar, collapse = "+")))}
  
  models_prs <- map(outcomes, ~ glm(form_prs(.x), data = data, family = binomial))
  models_quintile <- map(outcomes, ~ glm(form_quintile(.x), data = data, family = binomial))
  
  bind_rows(
    map_dfr(models_prs, broom::tidy, conf.int = TRUE, .id = "model_index"),
    map_dfr(models_quintile, broom::tidy, conf.int = TRUE, .id = "model_index")
  ) %>%
    filter(term %in% c("st.score", "comparison")) %>%
    mutate(
      pheno = outcomes[as.numeric(model_index)],
      OR = exp(estimate),
      CI_l = exp(conf.low),
      CI_u = exp(conf.high)
    ) %>%
    dplyr::select(-model_index) %>%
    relocate(pheno, .before = term)
}
