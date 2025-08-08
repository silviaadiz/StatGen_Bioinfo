
library(glmnet)
library(pROC)

# First, we compute the AUC corrected for optimism:

boot_lasso <- function(data, pheno, covar, interactions = FALSE, interactions_covar = NULL, B) {
  set.seed(1711)

  # Filter rows with missing values (pending imputation)
  dat_na <- data[complete.cases(data[, c(as.character(pheno), as.character(covar))]), , drop = FALSE]
 
  # Build model formula
  if (interactions && !is.null(interactions_covar)) {
    f1 <- as.formula(sprintf('%s ~ %s + %s', as.character(pheno), paste(covar, collapse = "+"), paste(interactions_covar, collapse = "+")))
  } else if (!interactions) {
    f1 <- as.formula(sprintf('%s ~ %s', as.character(pheno), paste(covar, collapse = "+")))
  }
 
  # Create model matrix 
  x <- model.matrix(f1, data = dat_na)
  n <- nrow(dat_na)
 
  # Initialize matrices for results
  coef_bootstrap_1se <- matrix(NA, nrow = B, ncol = ncol(x) + 1)
  colnames(coef_bootstrap_1se) <- c("intercept", colnames(x))

  boot.lambda <- matrix(NA, nrow = B, ncol = 1)
  colnames(boot.lambda) <- c("1se")
 
  preds <- matrix(NA, nrow = nrow(dat_na), ncol = B)
  colnames(preds) <- seq(1, B, 1)
 
  optimism <- matrix(NA, nrow = B, ncol = 1)  
 
  # Train full model
  boot_cv_lasso_full <- cv.glmnet(x = x, y = dat_na[, pheno], alpha = 1)
  full_model_auc <- auc(roc(dat_na[, pheno], predict(boot_cv_lasso_full, newx = x, s = boot_cv_lasso_full$lambda.1se)))
 
  # Bootstrapping
  for (i in 1:B) {
    ind <- sample(1:n, size = n, replace = TRUE)  
    x_boot <- x[ind, ]
    y_boot <- dat_na[ind, pheno]
   
    boot_cv_lasso <- cv.glmnet(x = x_boot, y = y_boot, alpha = 1)
   
    # Now we retrieve coefficients for the best lambda
    coef_bootstrap_1se[i, ] <- as.numeric(coef(boot_cv_lasso$glmnet.fit, s = boot_cv_lasso$lambda.1se))
    boot.lambda[i, ] <- boot_cv_lasso$lambda.1se
   
    # Predict in bootstrap resample
    preds_boot <- predict(boot_cv_lasso, newx = x_boot, type = "response", s = "lambda.1se")
   
    # Predict in test data (full sample)
    preds[, i] <- predict(boot_cv_lasso, newx = x, type = "response", s = "lambda.1se")
   
    # Compute AUC
    roc_curve_boot <- auc(roc(y_boot, as.vector(preds_boot)))
   
    roc_curve_test <- auc(roc(dat_na[, pheno], as.vector(preds[, i])))
   
    # Compute optimism for this resample
    optimism[i, ] <- roc_curve_boot - roc_curve_test
  }
 
  # Last, we compute the AUC for the model corrected for optimism
  auc_opt_cor_avg <- full_model_auc - mean(optimism)
 
  prob.nonzero.1se <- apply(abs(coef_bootstrap_1se) > 1e-5, 2, mean)
 
  # Final results
  res <- list(
    coef = coef_bootstrap_1se,  
    lambda = boot.lambda,
    prob.nonzero = prob.nonzero.1se,
    pred_cal = preds,
    optimism = optimism,  
    opt_cor_auc_promediado = auc_opt_cor_avg,  
    full_model_auc = full_model_auc  
  )
 
  return(res)
}


# Second: we would fit the full model and report the corrected AUC                  
