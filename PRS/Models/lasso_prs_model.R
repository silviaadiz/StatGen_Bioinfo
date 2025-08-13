library(glmnet)
library(pROC)

# First, we compute the AUC corrected for optimism:

boot_lasso <- function(data, pheno, covar, interactions = FALSE, interactions_covar = NULL, B = 1000) {
  set.seed(1711)
  
  dat_na <- data[complete.cases(data[, c(pheno, covar)]), ]
  
  if (interactions && !is.null(interactions_covar)) {
    f1 <- reformulate(c(covar, interactions_covar), response = pheno) } else {f1 <- reformulate(covar, response = pheno) }
  
  x <- model.matrix(f1, data = dat_na)
  n <- nrow(dat_na)
  
  coef_bootstrap_1se <- matrix(NA, nrow = B, ncol = ncol(x))
  colnames(coef_bootstrap_1se) <- colnames(x)
  
  boot_lambda <- numeric(B)
  
  preds <- matrix(NA, nrow = n, ncol = B)
  colnames(preds) <- paste0("Boot_", seq_len(B))
  
  optimism <- numeric(B)
  
  boot_cv_lasso_full <- cv.glmnet(x = x, y = dat_na[, pheno], alpha = 1, family = "binomial")
  full_model_pred <- predict(boot_cv_lasso_full, newx = x, s = "lambda.1se", type = "response")
  full_model_auc <- auc(roc(dat_na[, pheno], as.vector(full_model_pred)))
  
  for (i in seq_len(B)) {
    ind <- sample(seq_len(n), size = n, replace = TRUE)
    x_boot <- x[ind, , drop = FALSE]
    y_boot <- dat_na[ind, pheno]
    
    boot_cv_lasso <- cv.glmnet(x = x_boot, y = y_boot, alpha = 1, family = "binomial")
    
    # Now we retrieve coefficients for the best lambda
    coef_bootstrap_1se[i, ] <- as.numeric(coef(boot_cv_lasso, s = "lambda.1se"))
    boot_lambda[i] <- boot_cv_lasso$lambda.1se
    
    # Predict in bootstrap resample
    preds_boot <- predict(boot_cv_lasso, newx = x_boot, s = "lambda.1se", type = "response")
    auc_boot <- auc(roc(y_boot, as.vector(preds_boot)))
    
    # Predict in test data (full sample)
    preds_full <- predict(boot_cv_lasso, newx = x, s = "lambda.1se", type = "response")
    preds[, i] <- preds_full
    auc_test <- auc(roc(dat_na[, pheno], as.vector(preds_full)))
    
    # Compute optimism for this resample
    optimism[i] <- auc_boot - auc_test
  }
  
   # Lastly, we compute the AUC for the model corrected for optimism
  auc_opt_cor_avg <- full_model_auc - mean(optimism)

  prob_nonzero_1se <- colMeans(abs(coef_bootstrap_1se) > 1e-5)
  
  # Final results
  res <- list(
    coef = coef_bootstrap_1se,
    lambda = boot_lambda,
    prob_nonzero = prob_nonzero_1se,
    pred_cal = preds,
    optimism = optimism,
    opt_cor_auc_promediado = auc_opt_cor_avg,
    full_model_auc = full_model_auc
  )
  
  return(res)
}
     
