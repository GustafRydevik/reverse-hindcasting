ess_rfun<-function (sims) 
{
  if (is.vector(sims)) 
    dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  acov <- lapply(1:chains, FUN = function(i) {
    cov <- acf(sims[, i], lag.max = n_samples - 1, plot = FALSE, 
               type = c("covariance"))
    cov$acf[, , 1]
  })
  acov <- do.call(cbind, acov)
  chain_mean <- apply(sims, 2, mean)
  mean_var <- mean(acov[1, ]) * n_samples/(n_samples - 1)
  var_plus <- mean_var * (n_samples - 1)/n_samples
  if (chains > 1) 
    var_plus <- var_plus + var(chain_mean)
  rho_hat_sum <- 0
  if(var_plus>0){
  for (t in 2:nrow(acov)) {
    rho_hat <- 1 - (mean_var - mean(acov[t, ]))/var_plus
    if (rho_hat < 0) 
      break
    rho_hat_sum <- rho_hat_sum + rho_hat
  }}
  ess <- chains * n_samples
  if (rho_hat_sum > 0) 
    ess <- ess/(1 + 2 * rho_hat_sum)
  ess
}