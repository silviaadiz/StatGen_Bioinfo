# Bootstrap Elastic Net para modelos de PRS con moitas covariables ou selección de SNPs (moitas variables)
#
# Para cada mostra bootstrap, busca os parámetros alpha e lambda óptimos 
# usando validación cruzada e devolve estatísticas de estabilidade para os coeficientes.
# Alpha controla o equilibrio entre penalizacións ridge (α=0) e lasso (α=1).
# Lambda controla a forza xeral da regularización (determina con que forza se shrinkean os coeficientes cara 0: lambda máis alto dará lugar a coef. máis pequenos)

# data: Dataframe con variables
# pheno: Nome da variable fenotípica
# covar: Vector cos nomes das covariables a incluír no modelo
# interactions: Lóxico, se incluír termos de interacción (por defecto: FALSE)
# interactions_covar: Vector cos termos de interacción a incluír se interactions=TRUE
# B: número de remostras bootstrap a xerar (por defecto: 1000)
# cores: número de núcleos para procesamento paralelo (por defecto: 4)
#         Recomendado: usar núcleos totais - 2 para rendemento óptimo
#
# Resultado--> Lista que contén:
#   coefficients: Matriz de coeficientes bootstrap (B x p+1)
#   lambda: valores lambda seleccionados para cada mostra bootstrap
#   alpha:  valores alpha seleccionados para cada mostra bootstrap
#   summary_stats: dataframe con frecuencia de selección de variables e estatísticas de coeficientes
#   bootstrap_samples: Número de mostras bootstrap usadas
#   n_observations: tamaño final da mostra despois de eliminar datos perdidos
#   formula_used: Fórmula usada para o axuste do modelo
#
# 

library(glmnet)
library(foreach)
library(doParallel)

boot_elasticnet <- function(data, pheno, covar, interactions = FALSE, 
                           interactions_covar = NULL, B = 1000, cores = 4) {
  
  # Verificar e instalar paquetes necesarios
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package foreach is missing!!")}
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package doParallel is missing!!")}
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package glmnet is missing!!")}
  
  library(foreach)
  library(doParallel)
  library(glmnet)
  
  set.seed(1711)
  
  complete_cases <- complete.cases(data[, c(pheno, covar)])
  dat_na <- data[complete_cases, , drop = FALSE]
  
  
  if (interactions && !is.null(interactions_covar)) {
    formula_str <- paste0(as.character(pheno), " ~ ",
                          paste(covar, collapse = " + "), " + ",
                          paste(interactions_covar, collapse = " + "))
  } else {
    formula_str <- paste0(as.character(pheno), " ~ ",
                          paste(covar, collapse = " + "))
  }
  
  # Construír matriz do modelo e extraer compoñentes
  f1 <- as.formula(formula_str)
  x <- model.matrix(f1, data = dat_na)
  x <- x[, !colnames(x) %in% "(Intercept)", drop = FALSE]  # Eliminar columna de intercepto
  y <- dat_na[[pheno]]
  n <- nrow(dat_na)
  
  # Definir cuadrícula alpha para búsqueda de parámetros elastic net
  alphas_grid <- seq(0, 1, by = 0.1)  # 0=Ridge, 1=Lasso, 0.5=Elastic Net
  
  # Configurar procesamento paralelo
  cl <- NULL
  tryCatch({
    cl <- makeCluster(cores)  # Crear cluster co número especificado de núcleos
    registerDoParallel(cl)    # Rexistrar cluster para foreach
    # Foreach divide as iteracioóns entre clústers (workers)ç
    # %dopar% indica que cada iteración debe facerse en parealelo
    # export exporta variables adicionais que NON están no loop pero se utilizan porque os workers
    # non teñen acceso ao global environment
    # establecemos unha semilla diferente por iteración para asegurarnos de que non hai correlación
    # entre remostras bootstrap
    
  
    base_seed <- 1711
    
    # Procedemento bootstrap paralelo
    # Cada iteración realiza unha mostra bootstrap con optimización de hiperparámetros
    resultados <- foreach(
      i = 1:B,
      .combine = "rbind",
      .packages = c("glmnet"),
      .export = c("x", "y", "alphas_grid", "base_seed", "n")
    ) %dopar% {
      
      set.seed(base_seed + i)
      
      # Xerar índices da mostra bootstrap (con reemplazo)
      ind <- sample(1:n, size = n, replace = TRUE)
      x_boot <- x[ind, ]
      y_boot <- y[ind]
      
      best_alpha <- NA
      best_lambda <- NA
      best_auc <- -Inf
      best_coef <- rep(NA, ncol(x))  
      
      # gride search sobre o parámetro alpha
      for (a in alphas_grid) {
        cv_fit <- glmnet::cv.glmnet(
          x = x_boot,
          y = y_boot,
          alpha = a,                    # Parámetro de mixture elastic net
          family = "binomial",          
          nfolds = 5,                   
          type.measure = "auc",         # Optimizar para AUC
          parallel = FALSE              # non se paraleliza dentro do loop
        )
        
        # Extraer mellor AUC para o alpha actual
        curr_auc <- max(cv_fit$cvm, na.rm = TRUE)
        
        # Actualizar mellor modelo se o AUC actual é mellor
        if (curr_auc > best_auc) {
          best_alpha <- a
          best_lambda <- cv_fit$lambda.1se  
          best_auc <- curr_auc
          best_coef <- as.numeric(coef(cv_fit, s = "lambda.1se"))[-1]
        }
      }
      
      c(best_coef, best_lambda, best_alpha)
    }
    
  }, finally = {
    # pechamos  o cluster para liberar recursos
    if (!is.null(cl)) {
      stopCluster(cl)
    }
  })
  
  # Procesar resultados bootstrap
  n_coefs <- ncol(x)  
  coef_bootstrap <- resultados[, 1:n_coefs]              
  boot_lambda <- resultados[, n_coefs + 1]            
  boot_alpha <- resultados[, n_coefs + 2]              
  
  colnames(coef_bootstrap) <- colnames(x)
  
  prob.nonzero.1se <- apply(abs(coef_bootstrap) > 1e-5, 2, FUN = mean)  
  mean_coef <- apply(coef_bootstrap, 2, mean)                            
  sd_coef <- apply(coef_bootstrap, 2, sd)                                
  
  stats_table <- data.frame(
    variable = colnames(coef_bootstrap),
    frequency = prob.nonzero.1se,
    mean_coef = round(mean_coef, 4),
    sd_coefficient = round(sd_coef, 4)
  )
  
  results <- list(
    coefficients = coef_bootstrap,
    lambda = boot_lambda,
    alpha = boot_alpha,
    summary_stats = stats_table,
    bootstrap_samples = B,
    n_observations = n,
    formula_used = formula_str
  )
  
  return(results)
}
