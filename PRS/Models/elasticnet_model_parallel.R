  
boot_elasticnet <- function(data, pheno, covar, interactions = FALSE, interactions_covar = NULL, B = 1000, cores = 4) {
  if (!requireNamespace("foreach", quietly = TRUE)) {
    install.packages("foreach")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    install.packages("doParallel")
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    install.packages("glmnet")
  }
  
  library(foreach)
  library(doParallel)
  library(glmnet)
  
  set.seed(1711)
  
  # Arrange data (remove cases with NA)
  complete_cases <- complete.cases(data[, c(as.character(pheno), as.character(covar))])
  dat_na <- data[complete_cases, , drop = FALSE]
  
  # Then we build the formula strings from the arguments, with or without interactions
  
  if (interactions && !is.null(interactions_covar)) {
    formula_str <- paste0(as.character(pheno), " ~ ",
                          paste(covar, collapse = " + "), " + ",
                          paste(interactions_covar, collapse = " + "))
  } else {
    formula_str <- paste0(as.character(pheno), " ~ ",
                          paste(covar, collapse = " + "))
  }
  
  # Building the model matrix and other parameter, such as the alpha grid for search
  f1 <- as.formula(formula_str)
  x <- model.matrix(f1, data = dat_na)
  x <- x[, !colnames(x) %in% "(Intercept)", drop = F]
  y <- dat_na[[pheno]]
  n <- nrow(dat_na)
  alphas_grid <- seq(0, 1, by = 0.1)
  
  # Set parallelization parameters 
  cl <- NULL
  tryCatch({ # error proof 
    cl <- makeCluster(cores) # create cluster with number of specified cores 
    # Note: usually optimal number of cores available - 2
    registerDoParallel(cl) # register the cluster so we can use it with foreach
    
    # Export variables to workers, it copies the specified variables (base_seed) to all workers, which
    # are independent R processes
    base_seed<-1711

    # Parallelized bootstrap
    # Foreach divides the iterations between the workers
    # %dopar% indicates that eatch iteration should be executed in parallel
    # export sends additional variables that are not inside the loop but are necessary
    # because the workers don't have access to the global environment
    # By setting a different seed for each iteration, we ensure that there are no correlations
    # between the bootrstrap resamples
    
    resultados <- foreach(
      i = 1:B,
      .combine = "rbind",
      .packages = c("glmnet"),
      .export = c("x", "y", "alphas_grid")
    ) %dopar% {
      set.seed(base_seed + i)
      ind <- sample(1:n, size = n, replace = TRUE)
      x_boot <- x[ind, ]
      y_boot <- y[ind]
      
      best_alpha <- NA
      best_lambda <- NA
      best_auc <- -Inf
      best_coef <- rep(NA, ncol(x) + 1)
      
      # Alpha grid search
      for (a in alphas_grid) {
        cv_fit <- glmnet::cv.glmnet(
          x = x_boot,
          y = y_boot,
          alpha = a,
          family = "binomial",
          nfolds = 5,
          type.measure = "auc",
          parallel = FALSE  # No parallelization within the loop
        )
        
        current_auc <- max(cv_fit$cvm, na.rm = TRUE)
        if (current_auc > best_auc) {
          best_alpha <- a
          best_lambda <- cv_fit$lambda.1se
          best_auc <- current_auc
          best_coef <- as.numeric(coef(cv_fit, s = "lambda.1se"))
        }
      }
      
      c(best_coef, best_lambda, best_alpha)
    }
  }, finally = {
    if (!is.null(cl)) {
      stopCluster(cl)
    } # we close the cluster
  })
  
  n_coefs <- ncol(x) + 1  # +1 por el intercepto
  coef_bootstrap <- resultados[,1:n_coefs]
  boot_lambda <- resultados[,n_coefs + 1]
  boot_alpha <- resultados[,n_coefs + 2]
  
  colnames(coef_bootstrap) <- c("(Intercept)", colnames(x))
  prob.nonzero.1se<-apply(abs(coef_bootstrap)>1e-5,2,FUN=mean) # Frequency of each covariable to be more than 0
  mean_coef <- apply(coef_bootstrap, 2, mean)
  sd_coef <- apply(coef_bootstrap, 2, sd)
  
  stats_table <- data.frame(
    variable = colnames(coef_bootstrap),
    frequency = prob.nonzero.1se,
    mean_coef = round(mean_coef, 4),
    sd_coefficient = round(sd_coef, 4)
  )
  
 results= list(
    coefficients = coef_bootstrap,
    lambda = boot_lambda,
    alpha = boot_alpha,
    summary_stats = stats_table,
    bootstrap_samples = B,
    n_observations = n,
    formula_used = formula_str
  )
}
