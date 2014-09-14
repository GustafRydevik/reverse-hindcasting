monitor<-function (sims, warmup = floor(dim(sims)[1]/2), probs = c(0.025, 
    0.25, 0.5, 0.75, 0.975), digits_summary = 1, print = TRUE, 
    ...) 
{
    dim_sims <- dim(sims)
    dimnames_sims <- dimnames(sims)
    parnames <- dimnames_sims[[3]]
    if (length(dim_sims) != 3) 
        stop("'sims' is not a 3-d array")
    if (warmup > dim_sims[1]) 
        stop("warmup is larger than the total number of iterations")
    num_par <- dim_sims[3]
    if (is.null(parnames)) 
        parnames <- paste0("V", 1:num_par)
    sims_wow <- if (warmup >= 1) 
        apply(sims, c(2, 3), FUN = function(x) x[-(1:warmup)])
    else sims
    m <- apply(sims_wow, 3, mean)
    sd <- sapply(1:num_par, FUN = function(i) sd(as.vector(sims_wow[, 
        , i])))
    quan <- lapply(1:num_par, FUN = function(i) quantile(sims_wow[, 
        , i], probs = probs))
    probs_str <- names(quan[[1]])
    quan <- do.call(rbind, quan)
    rhat <- sapply(1:num_par, FUN = function(i) rstan:::split_rhat_rfun(sims_wow[, 
        , i]))
    ess <- sapply(1:num_par, FUN = function(i){if (var(c(sims_wow[,, i]))==0){
        length(sims_wow[,,i])}else{
            rstan:::ess_rfun(sims_wow[, 
        , i])}})
    sem <- sd/sqrt(ess)
    summary <- cbind(m, sem, sd, quan, ess, rhat)
    colnames(summary) <- c("mean", "se_mean", "sd", probs_str, 
        "n_eff", "Rhat")
    rownames(summary) <- parnames
    if (print) {
        cat("Inference for the input samples (")
        cat(dim_sims[2], " chains: each with iter=", dim_sims[1], 
            "; warmup=", warmup, "):\n\n", sep = "")
        summary[, "n_eff"] <- round(summary[, "n_eff"], 0)
        print(round(summary, digits_summary), ...)
        cat("\nFor each parameter, n_eff is a crude measure of effective sample size,\n", 
            "and Rhat is the potential scale reduction factor on split chains (at \n", 
            "convergence, Rhat=1).\n", sep = "")
    }
    invisible(summary)
}