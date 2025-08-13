
# Función para calcular o incremento no pseudo-R2 asociado ao PRS con bootstrap.
# Calcula tamén o R2 do modelo nulo con covariables, do modelo co PRS só
# e do modelo completo.
# O incremento no R2 corresponde a R2_modelo.completo-R2_modelo_covariables.
# Emprega remostraxe bootstrap para a estimación. Tamén calcula os IC para todas as estimas.
#
# Data: dataframe cos datos
# pheno_var: string con nome do fenotipo
# covar: string coas covariables tipo COVAR1+COVAR2..etcc
# score_var: variable co nome do PRS
# B: número de remostras
# parallel.boot: modo de paralelización de boot ("no","multicore","snow"). Ver axuda da función boot.



library(rms)
library(boot)

calc_r2_inc_bootstrap <- function(data, pheno_var, covar, score_var = "st.score", B = 1000, parallel.boot="multicore") {

f.cov <- reformulate(covar, response = pheno_var)
f.full <- reformulate(c(covar, score_var), response = pheno_var)
f.prs <- reformulate(score_var, response = pheno_var)
  
  dif_rsq <- function(dat, idx) {
    bt <- dat[idx, ]
    mod.cov <- lrm(f.cov, data = bt)
    mod.full <- lrm(f.full, data = bt)
    mod.prs <- lrm(f.prs, data = bt)
    
    r2_cov.mod <- mod.cov$stats["R2"]
    r2_full.mod <- mod.full$stats["R2"]
    r2_prs.mod <- mod.prs$stats["R2"]
    r2_inc <- r2_full.mod - r2_cov.mod
    
    c(r2_cov.mod,r2_full.mod,r2_prs.mod,r2_inc)}
  
  set.seed(394855)
  results <- boot(data = data, statistic = dif_rsq, R = B, parallel = parallel.boot)
  
  outp <- colMeans(results$t)
  cis_list <- lapply(1:ncol(results$t), function(i) {
    ci <- tryCatch(
      boot.ci(results, type = "basic", index = i), 
      error = function(e) NULL)
    
    if (!is.null(ci) && !is.null(ci$basic)) {
      c(low = ci$basic[4], high = ci$basic[5])
    } else {
      c(low = NA, high = NA)}
  })
  
  cis.data <- do.call(rbind, cis_list)
 
  results.fin <- data.frame(
    model = c("r2.mod.cov","r2.mod.full","r2.mod.prs","r2.increase.prs"),
    r2 = outp,
    ic_low = cis.data[,"low"],
    ic_high = cis.data[,"high"]
  )

  return(results.fin)
}
