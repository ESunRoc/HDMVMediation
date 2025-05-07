boot_stat <- function(data, indices) {
  # Subset the data using bootstrap indices
  U_boot <- data$U[indices, , drop = FALSE]
  X_boot <- data$X[indices]
  M_boot <- data$M[indices, , drop = FALSE]
  Y_boot <- data$Y[indices, , drop = FALSE]
  
  # Stage 1: Linear model
  XU_boot <- cbind(1, X_boot, U_boot)
  stage1_boot_mod <- .lm.fit(x = XU_boot, y = M_boot)
  stage1_boot_coef <- stage1_boot_mod$coefficients[2, ]
  
  # Stage 2: MSGLasso
  stage2_boot_mod <- MSGLasso(
    X.m = cbind(M_boot, U_boot, X_boot),
    Y.m = Y_boot,
    data$grpWTs,
    data$Pen_L,
    data$Pen_G,
    data$PQgrps,
    data$GRgrps,
    data$grp_Norm0,
    data$MSGLassolam1,
    data$MSGLassolamG.m
  )
  stage2_boot_selected <- stage2_boot_mod$Beta
  
  # Calculate effects
  boot_DE <- stage2_boot_selected[(data$p + data$l + 1):(data$p + data$l + data$q), ]
  boot_pide <- stage1_boot_coef * stage2_boot_selected[1:data$p, ]
  boot_tide <- colSums(boot_pide)
  
  # Return as a vector
  c(c(boot_pide), boot_tide, boot_DE)
}

BootSim <- function(sim_sett, nsim = 500, nB = 5e3, alpha = 0.05, ncores = 8) {
  #### initialize simulation params and results ####
  q <- sim_sett$q; p <- sim_sett$p; p0 <- sim_sett$p0; n <- sim_sett$n
  k <- sim_sett$k; l <- sim_sett$l; corr <- sim_sett$corr
  
  results <- list(
    "data_genSeeds"  = (1:nsim)^2 + (1:nsim),
    "true_coefs"     = matrix(nrow = k * p, ncol = nsim),
    "orig_coefs"     = matrix(nrow = k * (p + q + 1), ncol = nsim),
    "Boot_Draw_Seed" = ((1:nB) + 1)^2,
    # "Boot_Coefs"     = matrix(ncol = nB * nsim, nrow = k * (p + q + 1)),
    "pvals"          = matrix(ncol = nsim, nrow = k * (p + q + 1)),
    "bca_ints_lower" = matrix(ncol = nsim, nrow = k * (p + q + 1)),
    "bca_ints_upper" = matrix(ncol = nsim, nrow = k * (p + q + 1)),
    "piv_ints_lower" = matrix(ncol = nsim, nrow = k * (p + q + 1)),
    "piv_ints_upper" = matrix(ncol = nsim, nrow = k * (p + q + 1))
  )
  
  #### Stage 2 model preliminaries ####
  P <- p + l + q
  Q <- k
  
  G <- p + 2
  if (k == 5) { R <- k } else { R <- k }
  
  gmax <- 1
  cmax <- l
  
  GarrStarts <- c(0:(p - 1), p, p + l + 1)
  GarrEnds <- c(0:(p - 1), p + l, p + l + 2)
  
  if (k == 5) {
    # RarrStarts <- c(0, 2, 3, 4); RarrEnds <- c(1, 2, 3, 4)
    RarrStarts <- c(0:4); RarrEnds <- c(0:4)
  } else {
    # RarrStarts <- c(0, 3, 5, 6, 7, 8, 9); RarrEnds <- c(2, 4, 5, 6, 7, 8, 9)
    RarrStarts <- c(0:(k-1)); RarrEnds <- c(0:(k-1))
  }
  
  tmp_PQgrps <- FindingPQGrps(P = P, Q = Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQgrps <- tmp_PQgrps$PQgrps
  
  tmp_grpWts <- Cal_grpWTs(P = P, Q = Q, G, R, gmax, PQgrps)
  grpWTs <- tmp_grpWts$grpWTs
  
  tmp_GRgrps <- FindingGRGrps(P = P, Q = Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps <- tmp_GRgrps$GRgrps
  
  Pen_L <- matrix(1, P, Q)
  Pen_L[(p + 1):P, ] <- 0
  
  Pen_G <- matrix(1, G, R)
  Pen_G[(G - 1):G, ] <- 0
  
  grp_Norm0 <- matrix(1, G, R)
  
  assign("Pen_L", Pen_L, envir=.GlobalEnv); assign("Pen_G", Pen_G, envir=.GlobalEnv) # yeah. this again. still don't know why
  
  lam1.v <- seq(1e-5, 0.05, length = 15)
  lamG.v <- seq(1e-5, 0.05, length = 15)
  
  #### Run nsim datasets ####
  start_time <- proc.time()
  for (sim in 1:nsim) {
    #### Generate data ####
    orig_data <- gen_data(n = n, p = p, p0 = p0, l = l, q = q, k = k, 
                          seed = sim^2 + sim, corr = corr, quiet = T, b_scale = 2)
    X <- orig_data$data$X; M <- orig_data$data$M
    U <- orig_data$data$U; Y <- orig_data$data$Y
    
    #### Fast Stage 1 model using .lm.fit ####
    XU <- cbind(1, X, U)
    stage1_mod <- .lm.fit(x = XU, y = M)
    stage1_coefs <- stage1_mod$coefficients[2, ]
    
    #### Stage 2 model ####
    invisible(capture.output(step2_try_cv <- MSGLasso.cv(X = cbind(M, U, X), Y = Y,
                                                         grpWTs, Pen_L, Pen_G, PQgrps, 
                                                         GRgrps, lam1.v, lamG.v,
                                                         fold = 5, seed = 7^7)))
    
    best_lam <- step2_try_cv$lams.c[[which.min(step2_try_cv$rss.cv)]]
    MSGLassolam1 <- best_lam$lam1
    MSGLassolamG.m <- matrix(best_lam$lam3, G, R)
    MSGLassolamG.m[(G - 1):G, ] <- 0
    
    stage2_mod <- MSGLasso(X.m = cbind(M, U, X), Y.m = Y,
                           grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                           grp_Norm0, MSGLassolam1, MSGLassolamG.m)
    
    stage2_mod_selected <- stage2_mod$Beta
    # rownames(stage2_mod_selected) <- colnames(cbind(M, U, X))
    # colnames(stage2_mod_selected) <- colnames(Y)
    
    #### Extract effects using vectorized operations ####
    mod_DE <- stage2_mod_selected[(p+l+1):(p+l+q), ]
    stage2_mod_pide <- stage2_mod_selected[1:p, ]
    mod_pide_mat <- stage1_coefs * stage2_mod_pide  # Vectorized multiplication
    mod_tide_mat <- colSums(mod_pide_mat)
    
    mod_origFit_mat <- c(c(mod_pide_mat), mod_tide_mat, mod_DE)
    names(mod_origFit_mat) <- c(
      paste0(rownames(stage2_mod_selected)[1:p], "_ide_resp", rep(1:k, each = p)),
      paste0("TIDE_resp", 1:k),
      paste0("DE_resp", 1:k)
    )
    
    #### Calculate true PIDEs ####
    true_pides_mat <- orig_data$parameters$`True A_mat`[l + 1, ] * 
      orig_data$parameters$`True B_mat`[1:p, ]
    
    #### Bootstrap inference ####
    B_draws <- replicate(nB, sample(n, n, replace = TRUE))
    mod_bootRes <- matrix(0, nrow = length(mod_origFit_mat), ncol = nB)
    
    data_list <- list(U = U, X = X, M = M, Y = Y,
                      grpWTs = grpWTs, Pen_L = Pen_L, Pen_G = Pen_G,
                      PQgrps = PQgrps, GRgrps = GRgrps,
                      grp_Norm0 = grp_Norm0,
                      MSGLassolam1 = MSGLassolam1,
                      MSGLassolamG.m = MSGLassolamG.m,
                      p = p, l = l, q = q)
    
    
    start_time2 <- proc.time()
    boot_results <- boot(data = data_list, statistic = boot_stat, 
                         R = nB, parallel = "multicore", ncpus = ncores)
    

    mod_bootRes <- t(boot_results$t)
    
    mod_bootRes[is.nan(mod_bootRes)] <- 0
    
    #### Vectorized p-value calculation ####
    twice_orig <- 2 * mod_origFit_mat
    pos <- mod_origFit_mat >= 0
    
    count_ge <- rowSums(mod_bootRes >= twice_orig, na.rm = TRUE)
    count_lt <- rowSums(mod_bootRes < twice_orig, na.rm = TRUE)
    
    pvals <- numeric(length(mod_origFit_mat))
    pvals[pos] <- pmin(2 * count_ge[pos] / nB, 1)
    pvals[!pos] <- pmin(2 * count_lt[!pos] / nB, 1)
    
    #### Confidence intervals ####
    mod_boot_bcaCI <- t(apply(mod_bootRes, 1, coxed::bca))
    mod_boot_bcaCI[is.nan(mod_boot_bcaCI)] <- 0
    
    mod_boot_pivCI <- 2 * cbind(mod_origFit_mat, mod_origFit_mat) - 
      t(apply(mod_bootRes, 1, quantile, c(1 - alpha/2, alpha/2)))
    
    #### Store results ####
    results$orig_coefs[, sim] <- mod_origFit_mat
    results$true_coefs[, sim] <- c(true_pides_mat)
    # results$Boot_Coefs[, sim] <- rowMeans(mod_bootRes)
    results$pvals[, sim] <- pvals
    results$bca_ints_lower[, sim] <- mod_boot_bcaCI[, 1]
    results$bca_ints_upper[, sim] <- mod_boot_bcaCI[, 2]
    results$piv_ints_lower[, sim] <- mod_boot_pivCI[, 1]
    results$piv_ints_upper[, sim] <- mod_boot_pivCI[, 2]
    
    if (sim %% 5 == 0) {
      elapsed <- (proc.time() - start_time)[3]
      remaining <- (elapsed / sim) * (nsim - sim) / 60
      cat(sprintf("Sim %d/%d | Elapsed: %.1fm | Remaining: %.1fm\n",
                  sim, nsim, elapsed/60, remaining))
    }
  }
  return(results)
}