#### Simulation functions ####
# Generate a single dataset with specified sim parameters
gen_data <- function(n, p, p0, l, q, k, seed = seed, corr = c("low", "medium", "high"), quiet = T, b_scale = 1){
  #### Function information: ####
  # Requires: 
    # n, sample size
    # p, num. candidate mediators
    # p0, number of true mediators
    # l, num. of confounders
    # q, num. of treatments
    # k, num. of responses
    # corr, amount of correlation for response residuals
    # b_scale, the factor by which to adjust the signal size
  # Intermediaries:
    # W_mat, design matrix for step 1 model
    # Z_mat, design matrix for step 2 model
  # Returns:
    # test_data, a list of X, U, M, and Y needed to fit the models
    # sim_params, a list of the input parameters, drawn parameter matrices, and error matrices
  
  #### Confounders (n x l) ####
  set.seed(seed)
  U <- cbind(rbinom(n,1,0.25), rbinom(n,1,0.75), rbinom(n,1,0.50), rpois(n,15), runif(n),
             matrix(rnorm(n*(l-5)), nrow = n, ncol = l-5))
  colnames(U) <- paste0("U", 1:l)
  
  #### Treatment (n x q) ####
  set.seed(seed)
  X <- matrix(rbinom(n, 1, 0.5), ncol = q, nrow = n)
  colnames(X) <- paste0("X", 1:q)
  
  #### W_mat (n x l+q) ####
  W_mat <- cbind(U, X)
  
  #### alpha1 (l x p), confounder estimates (nuisance parameters) ####
  set.seed(1^1); alpha1 <- matrix(runif(l*p, min = -2, max = 2), nrow = l, ncol = p)
  rownames(alpha1) <- paste0("alpha1_", 1:l)
  
  #### xi1 (q x p), PIDE (trt --> mediator) ####
  set.seed(2^2)
  xi1 <- matrix(c(sample(c(-5,5), p0, replace = T), sample(c(-.05,.05), p-p0, replace = T)), nrow = q, ncol = p)
  rownames(xi1) <- paste0("xi1", 1:q)
  
  #### A_mat (l+q x p) ####
  A_mat <- rbind(alpha1, xi1)
  
  #### E1 (n x p), Mediator error matrix ####
  set.seed(seed)
  E1 <- mvrnorm(n=n, mu = rep(0,p), Sigma = diag(p))
  
  
  #### Mediators (n x p) ####
  M <- scale(W_mat%*%A_mat + E1)
  colnames(M) <- paste0("M", 1:p)
  
  # tmp <- lm(M ~ W_mat - 1); tmp_tidy <- broom::tidy(tmp); tmp_tidy[which(tmp_tidy$term=="W_matX1"),] %>% View
  
  #### Z_mat (n x q+p+l) ####
  Z_mat <- cbind(M, U, X)
  
  
  
  if(p0 == 5){
    if(quiet == T){
      suppressWarnings({
        beta <-      b_scale * rbind(rep(c(-8, -9,  4,  6,  4), k/5),        
                                     rep(c( 7, -8, -7, -4, -4), k/5),
                                     rep(c( 7, -7, -5,  5,  5), k/5),
                                     rep(c( 8, -7, -9,  4,  5), k/5),
                                     rep(c(-6, -8,  5, -4, -5), k/5),
                                     matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
      })
    } else {
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      beta <- b_scale * rbind(rep(c(-5, -3,  4,  6, 7), k/5),        
                              rep(c(-7,  6, -3, -4, 4), k/5),
                              rep(c(-4,  9, -5,  5, 3), k/5),
                              rep(c(-5,  8, -9,  4, 5), k/5),
                              rep(c(-8, -4,  3, -4, 9), k/5),
                              matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
    }
  } else{
    if(quiet == T){
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      suppressWarnings({
        beta <- b_scale * rbind(rep(c(-6, -4,  5,  7, 7), k/5),        
                                rep(c(-7,  6, -3, -4, 4), k/5),
                                rep(c(-4,  8, -5,  5, 3), k/5),
                                rep(c(-3,  8,  9,  5, 5), k/5),
                                rep(c(-3,  9, -5,  5, 3), k/5),
                                rep(c(-6,  8,  4,  7, 5), k/5),
                                rep(c( 8, -5,  3, -6, 9), k/5),
                                rep(c( 6,  9,  7, -8, 6), k/5),
                                rep(c( 5, -7,  5, -2, 5), k/5),
                                rep(c( 4, -3,  7,  5, 7), k/5),
                                matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
      })
    } else {
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      beta <- b_scale * rbind(rep(c(-6, -4,  5,  7, 7), k/5),        
                              rep(c(-7,  6, -3, -4, 4), k/5),
                              rep(c(-4,  8, -5,  5, 3), k/5),
                              rep(c(-3,  8,  9,  5, 5), k/5),
                              rep(c(-3,  9, -5,  5, 3), k/5),
                              rep(c(-6,  8,  4,  7, 5), k/5),
                              rep(c( 8, -5,  3, -6, 9), k/5),
                              rep(c( 6,  9,  7, -8, 6), k/5),
                              rep(c( 5, -7,  5, -2, 5), k/5),
                              rep(c( 4, -3,  7,  5, 7), k/5),
                              matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
    }
  }
  
  rownames(beta) <- paste0("beta_", 1:p)
  
  #### alpha2 (l x k), nuisance param(s) ####
  set.seed(4^4); alpha2 <- matrix(runif(l*k, -2, 2), nrow = l, ncol = k)
  
  #### xi2 (q x k), nuisance param(s) ####
  set.seed(5^5); xi2 <- matrix(sample(c(-1,1), q*k, replace=T), nrow = q, ncol = k)
  
  #### B_mat (p+l+q x k) ####
  B_mat <- rbind(beta, alpha2, xi2)
  
  #### E1 (n x k), Mediator error matrix ####
  set.seed(5^5)
  CorrMat_E2 <- gencor::gencor(d = k, method = corr, lim_low = 0.3, lim_medium = 0.7)$Matrix
  
  Sigma_E2 <- MBESS::cor2cov(cor.mat = CorrMat_E2, sd = rep(1,k))
  set.seed(seed)
  E2 <- mvnfast::rmvn(n=n, mu = rep(0,k), sigma = Sigma_E2)
  
  #### Y_mat (n x k), response matrix ####
  ## Generate from the model equation
  Y <- Z_mat %*% B_mat + E2 # correlation induced on residuals

  
  colnames(Y) <- paste0("Y", 1:k) 
  
  #### Define returns ####
  test_data <- list("X" = X, "U" = U, "M" = M, "Y" = Y)
  sim_params <- list("True B_mat" = B_mat, "True A_mat" = A_mat,
                     "E1" = E1, "E2" = E2, "ResponseCor" = CorrMat_E2)
  full_return <- list("data" = test_data, "parameters" = sim_params)
  return(full_return)
}

# Generate the matrices needed to initialize MSGLasso
set_MSGparams <- function(p, l, q, k){
  #### Function information: ####
  # Requires: 
    # p, num. candidate mediators
    # l, num. of confounders
    # q, num. of treatments
    # k, num. of responses
  # Intermediaries:
    # W_mat, design matrix for step 1 model
    # Z_mat, design matrix for step 2 model
  # Returns:
    # MSGparams, a list containing the needed (and intermediary) parameters for running MSGLasso
  
  P <- p + l + q
  Q <- k
  
  # Groups on X: p + 2; p mediator singleton groups + 1 confounder group + 1 treatment group
  G <- p + 2
  
  # Groups on Y: if k=2, 2; if k=5, 4; if k=10, 7
  if(k == 2) {R <- k} else if(k ==5) {R <- 4} else {R <- 7}
  
  gmax <- 1 # each variable (resp or pred) belongs to only 1 group
  cmax <- l # a group contains at most l variables (confounders form largest group)
  
  GarrStarts <- c(0:(p-1), p,  p+l+1)
  GarrEnds <-   c(0:(p-1), p+l, p+l+2)
  
  if(k == 2){
    RarrStarts <- c(0:(k-1)); RarrEnds <- c(0:(k-1))
  } else if(k==5){ # if k = 5: {(0,1), 2, 3, 4}
    # RarrStarts <- c(0,2,3,4); RarrEnds   <- c(1,2,3,4)
    RarrStarts <- c(0:(k-1)); RarrEnds <- c(0:(k-1))
  } else{ # if k > 5: {(0,1,2), (3,4), 5, ..., k}
    RarrStarts <- c(0:(k-1)); RarrEnds <- c(0:(k-1))
    # RarrStarts <- c(0,3,5:(k-1)); RarrEnds   <- c(2,4,5:(k-1)) 
  }
  
  tmp_PQgrps <- FindingPQGrps(P = P, Q = Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQgrps <- tmp_PQgrps$PQgrps
  
  tmp_grpWts <- Cal_grpWTs(P = P, Q = Q, G, R, gmax, PQgrps)
  grpWTs <- tmp_grpWts$grpWTs
  
  tmp_GRgrps <- FindingGRGrps(P = P, Q = Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps <- tmp_GRgrps$GRgrps
  
  penL <- matrix(rep(1, P*Q), P, Q, byrow=T)
  penL[(p+1):P,] <- 0 # don't penalize the confounders, (p+1):(P-1), or treatment, P
  
  penG <- matrix(rep(1,G*R),G,R, byrow=TRUE)
  penG[(G-1):G,] <- 0 # don't penalize confounder group, G-1, or treatment group, G
  
  grpNorm0 <- matrix(rep(1, G*R), nrow=G, byrow=TRUE)
  
  MSGparams <- list("P" = P, "Q" = Q, "G" = G, "R" = R, "gmax" = gmax, "cmax" = cmax, "PQgrps" = PQgrps, 
                    "grpWTs" = grpWTs, "GRgrps" = GRgrps, "Pen_L" = penL, "Pen_G" = penG, "grp.Norm0" = grpNorm0)
  return(MSGparams)
}

# Run a simulation given parameter settings and nsim generated datasets
run_sim <- function(sim_sett, nsim, b_scale = 1, drop_penn = F){
  #### Function information: ####
  # Requires: 
    # sim_sett, a list containing n, p, p0, l, k, q, corr
    # nsim, the number of data sets to simulate
    # b_scale, the factor by why to adjust signal size
    # drop_penn, a flag for whether to use smaller regularization grid
  # Intermediaries:
    # all_dataAndparams, a list containing all generated data sets and their true parameters
    # multivar_results, a list containing all step 1 models and multivariate step 2 models
    # univar_results, a list containing all step 1 models and univariate step 2 models
  # Returns:
    # Results, a list containing: all generated data sets, all step 1 model fits, all step 2 model fits
  
  #### Extract sim parameters ####
  n <- sim_sett$n; p <- sim_sett$p; p0 <- sim_sett$p0
  l <- sim_sett$l; k <- sim_sett$k; q <- sim_sett$q
  corr <- sim_sett$corr
  
  #### Initialize results storage ####
  multivar_results <- univar_results <- list()
  
  step1_univar <- step1_multivar <- matrix(nrow = 3*p, ncol = nsim)
  rownames(step1_univar) <- rownames(step1_multivar) <- c(paste0("MediatorNum", 1:p), 
                                                          paste0("Med",1:p,"TrtEst"), 
                                                          paste0("Med",1:p,"pval"))
  colnames(step1_univar) <- colnames(step1_multivar) <- paste0("sim",1:nsim)
  
  step2_multivar <- matrix(nrow = k*(p+l+q), ncol = nsim)
  if(k==2){
    rownames(step2_multivar) <- c(paste0("Beta", 1:(p+q+l), "_resp1"), paste0("Beta", 1:(p+q+l), "_resp2"))
  } else if(k==5){
    rownames(step2_multivar) <- c(paste0("Beta", 1:(p+q+l), "_resp1"), paste0("Beta", 1:(p+q+l), "_resp2"),
                                  paste0("Beta", 1:(p+q+l), "_resp3"), paste0("Beta", 1:(p+q+l), "_resp3"),
                                  paste0("Beta", 1:(p+q+l), "_resp5"))
  } else{
    rownames(step2_multivar) <- c(paste0("Beta", 1:(p+q+l), "_resp1"), paste0("Beta", 1:(p+q+l), "_resp2"),
                                  paste0("Beta", 1:(p+q+l), "_resp3"), paste0("Beta", 1:(p+q+l), "_resp4"),
                                  paste0("Beta", 1:(p+q+l), "_resp5"), paste0("Beta", 1:(p+q+l), "_resp6"), 
                                  paste0("Beta", 1:(p+q+l), "_resp7"), paste0("Beta", 1:(p+q+l), "_resp8"),
                                  paste0("Beta", 1:(p+q+l), "_resp9"), paste0("Beta", 1:(p+q+l), "_resp10"))
  }
  colnames(step2_multivar) <- paste0("sim",1:nsim)
  step2_univar <- step2_multivar
  
  pen_params_multivar <- matrix(nrow = 2, ncol = nsim); rownames(pen_params_multivar) <- c("Pen.L", "Pen.G")
  pen_params_univar <- matrix(nrow = k, ncol = nsim); rownames(pen_params_univar) <- paste0("resp",1:k)
  
  colnames(pen_params_multivar) <- colnames(pen_params_univar) <- paste0("sim", 1:nsim)
  
  all_dataAndparams <- list()
  
  #### Initialize MSGLasso parameters ####
  MSG_params <- set_MSGparams(p=p,l=l,q=q,k=k)
  P <- MSG_params$P; Q <- MSG_params$Q
  G <- MSG_params$G; R <- MSG_params$R
  
  gmax <- MSG_params$gmax ; cmax <- MSG_params$cmax
  
  PQgrps <- MSG_params$PQgrps; GRgrps <- MSG_params$GRgrps
  grpWTs <- MSG_params$grpWTs; grp_Norm0 <- MSG_params$grp.Norm0
  
  Pen_L <- MSG_params$Pen_L; Pen_G <- MSG_params$Pen_G
  
  # This is bad practice. It breaks if I try anything else; I'm not sure why.
  assign("Pen_L", Pen_L, envir=.GlobalEnv); assign("Pen_G", Pen_G, envir=.GlobalEnv)
  
  start_time <- proc.time()
  
  #### Repeat simulation nsim times ####
  for(sim in 1:nsim){
    #### Set up sim parameters ####
    seed <- sim # 1 to nsim
    simdata <- gen_data(n = n, p = p, p0 = p0, l=l, k = k, q = q, seed = seed, corr = corr, b_scale = b_scale)
    
    data <- simdata$data
    sim_params <- simdata$parameters
    
    X <- data$X; U <- data$U; M <- as.matrix(data$M); Y <- data$Y
    
    all_dataAndparams <- rlist::list.append(all_dataAndparams, simdata)
    
    #### Step 1 Model ####
    step1_model <- lm(M ~ X + U -1)
    step1_tidy <- broom::tidy(step1_model)
    step1_res <- cbind("MediatorNum" = 1:p, step1_tidy[which(step1_tidy$term == "X"),c(3,6)])
    step1_univar[,sim] <- step1_multivar[,sim] <- vec(as.matrix(step1_res))
    
    #### Step 2 multivariate ####

    if(drop_penn == T){
      lam1.v <- seq(0.0075, 0.075, length=10)
      lamG.v <- seq(0.0075, 0.075, length=10)
    } else{
      lam1.v <- seq(0.025, 0.15, length=25)
      lamG.v <- seq(0.025, 0.15, length=25)
    }
    
    # this prints (lam1.v,lamg.v) counter for each fold. This is suppressed so printing to console doesn't take too long
    invisible(capture.output(step2_mod_cv <- MSGLasso.cv(X = cbind(M,U,X), Y = Y,
                                                         grpWTs, Pen_L, Pen_G,
                                                         PQgrps, GRgrps, lam1.v, lamG.v,
                                                         fold = 5, seed = 7^7)))
    
    MSGLassolam1 <- step2_mod_cv$lams.c[which.min(as.vector(step2_mod_cv$rss.cv))][[1]]$lam1
    MSGLassolamG <- step2_mod_cv$lams.c[which.min(as.vector(step2_mod_cv$rss.cv))][[1]]$lam3
    MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
    MSGLassolamG.m[(G-1):G,] <- 0
    
    pen_params_multivar[,sim] <- c(MSGLassolam1, MSGLassolamG)
    
    step2_model <- MSGLasso(X.m = cbind(M, U, X), Y.m = Y,
                            grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                            grp_Norm0, MSGLassolam1, MSGLassolamG.m)
    
    step2_multivar[,sim] <- vec(step2_model$Beta)
    
    #### Step 2 univariate model ####
    penaltyFactor <- c(rep(1,p),rep(0,q+l)) # don't penalize coefs on U or X
    
    step2_univar_mod <- apply(Y, 2, function(val){
      step2Uni_train <- cv.glmnet(x = cbind(M,U,X), y = val, intercept = F, 
                                  penalty.factor = penaltyFactor, relax = T, 
                                  standardize = F, parallel = T)
      step2Uni <- glmnet(x = cbind(M,U,X), y = val, intercept = F, 
                         penalty.factor = penaltyFactor, 
                         relax = T, standardize = F, 
                         lambda = step2Uni_train$lambda.1se)
      return(c(step2Uni$beta,step2Uni_train$lambda.1se))
    })
    
    
    step2_univar_coef <- vector()
    for(i in 1:k){
      step2_univar_coef <- append(step2_univar_coef, as.matrix(step2_univar_mod[[i]][[1]])[c(1:p,p+q+l,(p+1):(p+q+l-1)),])
      pen_params_univar[i,sim] <- step2_univar_mod[[i]][[2]]
    }
    step2_univar_coef <- unname(step2_univar_coef)
    step2_univar[,sim] <- step2_univar_coef
    
    
    if(sim %% 5 == 0){
      cat(paste0("Done with simulation ", sim, "; ", nsim-sim, " remaining. Time elapsed: ", 
                 round(-1*(start_time[3] - proc.time()[3])/60, 3), " minutes. Apx. ",
                 round(((-1*(start_time[3] - proc.time()[3])/60)/sim)*(nsim-sim), 3), " minutes remaining\n"))
    }
  }
  (elapsed <- start_time - proc.time())
  
  #### Organize results ####
  multivar_results <- list("Step1Multivar" = step1_multivar, 
                           "Step2Multivar" = step2_multivar, 
                           "MultivarPenParams" = pen_params_multivar)
  univar_results <- list("Step1Univar" = step1_univar, 
                         "Step2Univar" = step2_univar, 
                         "UnivarPenParams" = pen_params_univar)
  results <- list("SimData" = all_dataAndparams, 
                  "MultivarResults" = multivar_results, 
                  "UnivarResults" = univar_results)
  
  return(results)
}


#### Summary functions ####
# Summarize results into object that can then be passed to table generator
summarize_results <- function(results, p, k, l = 5, q = 1, nsim = 500, n = 100, filter_step1 = F, alpha = 0.1, conf.level = 0.95){
  #### Function information ####
  # Requires:
    # results, a list output from an instance of "run_sim"
    # p, the number of candidate mediators
    # k, the number of responses
    # l, the number of confounders considered
    # q, the number of exposures considered
    # nsim, the number of datasets generated for "results"
    # n, the sample size used for "results"
    # filter_step1, a flag for whether to use a Baron-Kenney style significance check for first PIDE
    # alpha, the significance level at which to construct pseudo-CIs
    # conf.level, the confidence level at which to construct pseudo-CIs
  # Intermediaries:
    # Step2_Selection_Accuracy, data frame of how frequently MSGLasso selections matched selections of the oracle
    # Step2_Sign_Accuracy, data frame of how frequently the sign of the MSGLasso matched the sign of the oracle
    # Multivar_PIDE, data frame containing summary results for PIDEs via multivariate model
    # Univar_PIDE, data frame containing summary results for PIDEs via univariate model
    # Multivar_DE, data frame containing summary results for DEs via multivariate model
    # Univar_DE, data frame containing summary results for DEs via univariate model
  # Returns:
    # summary_results, a list containing all created summary variables
  
  simparams <- results$SimData[[1]][[2]]
  sim_TrueBmatMedis <- simparams[["True B_mat"]][1:p,]
  
  sim_TrueBmatMedis_vec <- as.matrix(vec(sim_TrueBmatMedis))
  if(k == 5){
    rownames(sim_TrueBmatMedis_vec) <- c(paste0("Beta",1:p,"_resp1"),paste0("Beta",1:p,"_resp2"),
                                         paste0("Beta",1:p,"_resp3"),paste0("Beta",1:p,"_resp4"),
                                         paste0("Beta",1:p,"_resp5"))
    ## Fix rownames for Beta()_resp4; these indices will need to be changed for different p and k. Not strictly needed.
    # rownames(results$MultivarResults$Step2Multivar)[154:204] <- paste0("Beta",1:51,"_resp4")
    # rownames(results$UnivarResults$Step2Univar)[154:204] <- paste0("Beta",1:51,"_resp4")
  } else {
    rownames(sim_TrueBmatMedis_vec) <- c(paste0("Beta",1:p,"_resp1"),paste0("Beta",1:p,"_resp2"),
                                         paste0("Beta",1:p,"_resp3"),paste0("Beta",1:p,"_resp4"),
                                         paste0("Beta",1:p,"_resp5"),paste0("Beta",1:p,"_resp6"),
                                         paste0("Beta",1:p,"_resp7"),paste0("Beta",1:p,"_resp8"),
                                         paste0("Beta",1:p,"_resp9"),paste0("Beta",1:p,"_resp10"))
  }
  
  colnames(sim_TrueBmatMedis_vec) <- "Truth"
  
  ## This only works for p=25. There will be a different sequence for p=100
  if(k == 5 & p == 25 & l == 25){
    step2_beta_indices <- c(1:25, 52:76, 103:127, 154:178, 205:229)
  } else if(k == 10 & p == 25  & l == 25){
    step2_beta_indices <- c(1:25, 52:76, 103:127, 154:178, 205:229, 256:280, 307:331, 358:382, 409:433, 460:484)
  } else if(k == 5  & p == 100 & l == 25){
    step2_beta_indices <- c(1:100, 127:226, 253:352, 379:478, 505:604)
  } else if(k == 10 & p == 100 & l == 25){
    step2_beta_indices <- c(1:100, 127:226, 253:352, 379:478, 505:604, 631:730, 757:856, 883:982, 1009:1108, 1135:1234)
  } else if(k == 5 & p == 25 & l == 5){
    step2_beta_indices <- c(1:25, 32:56, 63:87, 94:118, 125:149)
  } else if(k == 10 & p == 25 & l == 5){
    step2_beta_indices <- c(1:25, 32:56, 63:87, 94:118, 125:149, 156:180, 187:211, 218:242, 249:273, 280:304)
  } else if(k == 5 & p == 100 & l == 5){
    step2_beta_indices <- c(1:100, 107:206, 213:312, 319:418, 425:524)
  } else if(k == 10 & p == 100 & l == 5){
    step2_beta_indices <- c(1:100, 107:206, 213:312, 319:418, 425:524, 531:630, 637:736, 743:842, 849:948, 955:1054)
  }
  
  sim_multivarMedisStep2 <- results$MultivarResults$Step2Multivar[step2_beta_indices,]
  sim_univarMedisStep2 <- results$UnivarResults$Step2Univar[step2_beta_indices,]
  
  #### Step 2 selection and sign consistency ####
  Step2_sign_const <- matrix(nrow = nrow(sim_TrueBmatMedis_vec), ncol = 3); Step2_sign_const[,1] <- sign(sim_TrueBmatMedis_vec[,1])
  rownames(Step2_sign_const) <- rownames(sim_TrueBmatMedis_vec)
  colnames(Step2_sign_const) <- c("Truth", "Multivar.Prop", "Univar.Prop")  
  
  Step2_select_const <- matrix(nrow = nrow(sim_TrueBmatMedis_vec), ncol = 3); Step2_select_const[,1] <- ifelse(sim_TrueBmatMedis_vec[,1]!=0,1,0)
  rownames(Step2_select_const) <- rownames(sim_TrueBmatMedis_vec)
  colnames(Step2_select_const) <- c("Truth", "Multivar.Prop", "Univar.Prop")  
  
  for(i in 1:nrow(sim_TrueBmatMedis_vec)){
    multivar_signConst <- sum(sign(sim_multivarMedisStep2[i,])==sign(sim_TrueBmatMedis_vec[i,]))/500
    Step2_sign_const[i,2] <- multivar_signConst
    
    multivar_selected <- (sim_multivarMedisStep2[i,]!=0)
    multivar_selectConst <- sum(multivar_selected == (sim_TrueBmatMedis_vec[i,]!=0))/500
    Step2_select_const[i,2] <- multivar_selectConst
    
    univar_signConst <- sum(sign(sim_univarMedisStep2[i,])==sign(sim_TrueBmatMedis_vec[i,]))/500
    Step2_sign_const[i,3] <- univar_signConst
    
    univar_selected <- (sim_univarMedisStep2[i,]!=0)
    univar_selectConst <- sum(univar_selected == (sim_TrueBmatMedis_vec[i,]!=0))/500
    Step2_select_const[i,3] <- univar_selectConst
  }
  
  #### Full selection and sign consistency ####
  Step1MediTrueSelect <- results$SimData[[1]][[2]][["True A_mat"]][l+1,] # true trt-->medi PIDE
  
  MediTrueSelect <- rep(Step1MediTrueSelect,k)*sim_TrueBmatMedis_vec
  MediTrueSigns <- as.matrix(sign(MediTrueSelect))
  MediTrueNonzero <- ifelse(MediTrueSigns !=0, 1, 0)
  
  if(k == 5){ #  & May need: p == 25
    rownames(MediTrueNonzero) <- rownames(MediTrueSelect) <- rownames(MediTrueSigns) <- c(paste0("pide",1:p,"_resp1"),
                                                                                          paste0("pide",1:p,"_resp2"),
                                                                                          paste0("pide",1:p,"_resp3"),
                                                                                          paste0("pide",1:p,"_resp4"),
                                                                                          paste0("pide",1:p,"_resp5"))
  } else { # May need:  if(p == 25)
    rownames(MediTrueNonzero) <- rownames(MediTrueSelect) <- rownames(MediTrueSigns) <- c(paste0("pide",1:p,"_resp1"),paste0("pide",1:p,"_resp2"),
                                                                                          paste0("pide",1:p,"_resp3"),paste0("pide",1:p,"_resp4"),
                                                                                          paste0("pide",1:p,"_resp5"),paste0("pide",1:p,"_resp6"),
                                                                                          paste0("pide",1:p,"_resp7"),paste0("pide",1:p,"_resp8"),
                                                                                          paste0("pide",1:p,"_resp9"),paste0("pide",1:p,"_resp10"))
  }
  colnames(MediTrueSelect) <- colnames(MediTrueNonzero) <- colnames(MediTrueSigns) <- NULL
  
  
  ## Stack k copies of Step 1 PIDEs
  tmp_block <- matrix(rep(1,k), nrow = k, ncol = 1)
  
  ## Again, these may/will break for p=100 but there's a good amount of time before I get there so they work for now
  if(filter_step1 == F){ ## Without step 1 filtering at alpha 
    ## The first step is the same in both models, so we subset from Step1Multivar, here wlog
    step1_multivar_ide <- step1_univar_ide <- tmp_block %x% results$MultivarResults$Step1Multivar[(p+1):(2*p),] 
  } else { ## With step 1 filtering at alpha
    Step1_pide <- results$MultivarResults$Step1Multivar[(p+1):(2*p),] 
    for(j in 1:nsim){
      Step1_pide[,j] <- ifelse(results$MultivarResults$Step1Multivar[(2*p+1):(3*p),j]>alpha, 
                               0, results$MultivarResults$Step1Multivar[(2*p+1):(3*p),j])
    }
    step1_multivar_ide <- step1_univar_ide <- tmp_block %x% Step1_pide
  }
  
  
  
  univar_pide <- multivar_pide <- matrix(nrow = p*k, ncol = nsim)
  if(k == 5){
    rownames(univar_pide) <- rownames(multivar_pide) <- c(paste0("pide",1:p,"_resp1"),paste0("pide",1:p,"_resp2"),
                                                          paste0("pide",1:p,"_resp3"),paste0("pide",1:p,"_resp4"),
                                                          paste0("pide",1:p,"_resp5"))
  } else {
    rownames(univar_pide) <- rownames(multivar_pide) <- c(paste0("pide",1:p,"_resp1"),paste0("pide",1:p,"_resp2"),
                                                          paste0("pide",1:p,"_resp3"),paste0("pide",1:p,"_resp4"),
                                                          paste0("pide",1:p,"_resp5"),paste0("pide",1:p,"_resp6"),
                                                          paste0("pide",1:p,"_resp7"),paste0("pide",1:p,"_resp8"),
                                                          paste0("pide",1:p,"_resp9"),paste0("pide",1:p,"_resp10"))
  }
  colnames(univar_pide) <- colnames(multivar_pide) <- paste0("nsim",1:nsim)
  
  for(i in 1:ncol(sim_multivarMedisStep2)){
    univar_pide[,i] <- unname(sim_univarMedisStep2[,i] * step1_multivar_ide[,i])
    multivar_pide[,i] <- unname(sim_multivarMedisStep2[,i] * step1_multivar_ide[,i])
  }
  
  # the rows of univar_pide and multivar_pide are now all 500 generated pides for their respective model
  # That is, we now have empirical pseudo-distributions (each conditional on its sample) of what the pide looks like. 
  # This still isn't good for things like intervals since the data sets are all different, but conditional on the distribution(s)
  # of the data generative process, we can compare these to the truth
  
  tmp_nonzero_multi <- tmp_sign_multi <- matrix(ncol=nsim, nrow = nrow(multivar_pide))
  tmp_nonzero_uni <- tmp_sign_uni <- matrix(ncol=nsim, nrow = nrow(univar_pide))
  
  for(i in 1:ncol(multivar_pide)){
    tmp_sign_multi[,i] <- sign(multivar_pide)[,i]==MediTrueSigns
    tmp_nonzero_multi[,i] <- ifelse(multivar_pide!=0,1,0)[,i]==MediTrueNonzero
    
    tmp_sign_uni[,i] <- sign(univar_pide)[,i]==MediTrueSigns
    tmp_nonzero_uni[,i] <- ifelse(univar_pide!=0,1,0)[,i]==MediTrueNonzero
  }
  
  multivar_pide_signPercent <- rowSums(tmp_sign_multi)/nsim
  multivar_pide_selectPercent <- rowSums(tmp_nonzero_multi)/nsim
  
  multivar_pide_mean <- apply(multivar_pide, 1, mean)
  multivar_pide_stde <- apply(multivar_pide, 1, sd)
  multivar_pide_Lowconf <- multivar_pide_mean - 1.96*multivar_pide_stde
  multivar_pide_Uppconf <- multivar_pide_mean + 1.96*multivar_pide_stde
  multivar_pide_summ <- data.frame("True_PIDE" = MediTrueSelect, 
                                   "True_Sign" = MediTrueSigns,
                                   "Multivar_PIDE_Mean" = multivar_pide_mean, 
                                   "Multivar_PIDE_SE" = multivar_pide_stde,
                                   "Multivar_PIDE_Lower" = multivar_pide_Lowconf,
                                   "Multivar_PIDE_Upper" = multivar_pide_Uppconf,
                                   "Multivar_PIDE_SelectPercent" = multivar_pide_selectPercent,
                                   "Multivar_PIDE_SignPercent" = multivar_pide_signPercent)
  
  
  univar_pide_signPercent <- rowSums(tmp_sign_uni)/nsim
  univar_pide_selectPercent <- rowSums(tmp_nonzero_uni)/nsim
  
  univar_pide_mean <- apply(univar_pide, 1, mean)
  univar_pide_stde <- apply(univar_pide, 1, sd)
  univar_pide_Lowconf <- univar_pide_mean - 1.96*univar_pide_stde
  univar_pide_Uppconf <- univar_pide_mean + 1.96*univar_pide_stde
  univar_pide_summ <- data.frame("True_PIDE" = MediTrueSelect, 
                                 "True_Sign" = MediTrueSigns,
                                 "Univar_PIDE_Mean" = univar_pide_mean, 
                                 "Univar_PIDE_SE" = univar_pide_stde,
                                 "Univar_PIDE_Lower" = univar_pide_Lowconf,
                                 "Univar_PIDE_Upper" = univar_pide_Uppconf,
                                 "Univar_PIDE_SelectPercent" = univar_pide_selectPercent,
                                 "Univar_PIDE_SignPercent" = univar_pide_signPercent)
  
  ## Handle direct effects:
  MediTrueDE <- results$SimData[[1]]$parameters$`True B_mat`[p+l+1,]
  
  DE_Muiltivar_allSims <- results$MultivarResults$Step2Multivar[which(grepl(as.character(p+l+1), rownames(results$MultivarResults$Step2Multivar))),]
  rownames(DE_Muiltivar_allSims) <- paste0("tau_",1:k)
  multivar_de_mean <- apply(DE_Muiltivar_allSims, 1, mean)
  multivar_de_sd <- apply(DE_Muiltivar_allSims, 1, sd)
  multivar_de_Lowconf <- multivar_de_mean - 1.96*multivar_de_sd
  multivar_de_Uppconf <- multivar_de_mean + 1.96*multivar_de_sd
  
  multivar_de_summ <- data.frame("True_DE" = MediTrueDE,
                                 "Multivar_DE_mean" = multivar_de_mean,
                                 "Multivar_DE_SE" = multivar_de_sd,
                                 "Multivar_DE_Lower" = multivar_de_Lowconf,
                                 "Multivar_DE_Upper" = multivar_de_Uppconf)
  rownames(multivar_de_summ) <- paste0("DE_",1:k)
  
  DE_Univar_allSims <- results$UnivarResults$Step2Univar[which(grepl(as.character(p+l+1), rownames(results$UnivarResults$Step2Univar))),]
  rownames(DE_Univar_allSims) <- paste0("tau_",1:k)
  univar_de_mean <- apply(DE_Univar_allSims, 1, mean)
  univar_de_sd <- apply(DE_Univar_allSims, 1, sd)
  univar_de_Lowconf <- univar_de_mean - 1.96*univar_de_sd
  univar_de_Uppconf <- univar_de_mean + 1.96*univar_de_sd
  
  univar_de_summ <- data.frame("True_DE" = MediTrueDE,
                               "Univar_DE_mean" = univar_de_mean,
                               "Univar_DE_SE" = univar_de_sd,
                               "Univar_DE_Lower" = univar_de_Lowconf,
                               "Univar_DE_Upper" = univar_de_Uppconf)
  rownames(univar_de_summ) <- paste0("DE_",1:k)
  
  
  
  summary_results <- list("Step2_Selection_Accuracy" = Step2_select_const,
                          "Step2_Sign_Accuracy" = Step2_sign_const,
                          "Multivar_PIDE" = multivar_pide_summ,
                          "Univar_PIDE" = univar_pide_summ,
                          "Multivar_DE" = multivar_de_summ,
                          "Univar_DE" = univar_de_summ)
  return(summary_results)
}

# Helper function to calculate mean of every fixed num of entries; from https://stackoverflow.com/questions/43635846/
BinMean <- function (vec, every, na.rm = FALSE) {
  #### Function information ####
  # Requires:
    # vec, an input vector
    # every, which items to include in the mean calculation
    # na.rm, how to handle missing values
  # Intermediaries:
  # Returns:
    # x, the vector of means of "every" items
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}

# Takes output of summarize_results and formats into a matrix that resembles 
summTblFunc <- function(res_summ, p, p0, k, return_tex = T){
  # Requires:
    # res_summ, summary output from summarize_results
    # p, the number of candidate mediators
    # p0, the number of true mediators; wlog the first p0 values of each column
    # k, the number of responses
    # return_tex, a flag for whether to directly copy to clipboard the generated tex table
  # Intermediaries:
    # res_multivar, summarized p0*k matrix of selection accuracy for true signal
    # res_univar, summarized 5*k matrix of selection frequency for true noise
  # Returns:
    # res_specific, matrix form of simulation tables; passable to knitr::kable to generate tex table
  
  #### Row indices ####
  if(k == 5){
    row_idx <- c(1:p0, (p+1):(p+p0), (2*p+1):(2*p+p0), (3*p+1):(3*p+p0), (4*p+1):(4*p+p0))
    row_idx_noise <- c((p0+1):(p0+5), (p+p0+1):(p+p0+5), (2*p+p0+1):(2*p+p0+5), (3*p+p0+1):(3*p+p0+5), (4*p+p0+1):(4*p+p0+5))
  } else if(k == 10){
    row_idx <- c(1:p0, (p+1):(p+p0), (2*p+1):(2*p+p0), (3*p+1):(3*p+p0), (4*p+1):(4*p+p0), 
                 (5*p+1):(5*p+p0), (6*p+1):(6*p+p0), (7*p+1):(7*p+p0), (8*p+1):(8*p+p0), (9*p+1):(9*p+p0))
    row_idx_noise <- c((0*p+p0+1):(0*p+p0+5), (1*p+p0+1):(1*p+p0+5), (2*p+p0+1):(2*p+p0+5), (3*p+p0+1):(3*p+p0+5), (4*p+p0+1):(4*p+p0+5),
                       (5*p+p0+1):(5*p+p0+5), (6*p+p0+1):(6*p+p0+5), (7*p+p0+1):(7*p+p0+5), (8*p+p0+1):(8*p+p0+5), (9*p+p0+1):(9*p+p0+5))
  }
  
  res_multivar <- rbind(matrix(c(res_summ$Multivar_PIDE$Multivar_PIDE_SelectPercent[row_idx], 
                                 1-res_summ$Multivar_PIDE$Multivar_PIDE_SelectPercent[row_idx_noise]),
                               ncol = k, byrow = T),
                        matrix(1-BinMean(res_summ$Multivar_PIDE$Multivar_PIDE_SelectPercent[-row_idx], p-p0), nrow = 1))
  
  res_univar <- rbind(matrix(c(res_summ$Univar_PIDE$Univar_PIDE_SelectPercent[row_idx],
                               1-res_summ$Univar_PIDE$Univar_PIDE_SelectPercent[row_idx_noise]), ncol = k, byrow = T),
                      matrix(1-BinMean(res_summ$Univar_PIDE$Univar_PIDE_SelectPercent[-row_idx], p-p0), nrow = 1))
  
  res_specific <- round(matrix(rbind(res_univar, res_multivar), nrow = nrow(res_univar)), 3)
  
  colnames(res_specific) <- paste0(c("Uni_Y", "Multi_Y"), rep(1:k, each = 2))
  rownames(res_specific) <- c(paste0("p", 1:p0), paste0("p",(p0+1):(p0+5)), "AvgNoiseSelAcc")
  
  if(return_tex==T){
    res_specific %>% knitr::kable(format = "latex") %>% clipr::write_clip()
  } else {
    return(res_specific)
  }
}


#### Inference/bootstrap functions ####
# Function to compute the bootstrap statistic; for use with "boot" function in R 
boot_stat <- function(data, indices) {
  # Requires:
    # data, a simulated data set
    # indices, the samples to use for the given iteration
  # Intermediaries:
  # Returns:
    # vector of parameters to be bootstrapped
  
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

# Function to run a simulation of bootstrap inference under given parameters
BootSim <- function(sim_sett, nsim = 500, nB = 5e3, alpha = 0.05, ncores = 8) {
  # Requires:
    # sim_sett, the simulation parameters to consider
    # nsim, the number of datasets to simulate
    # nB, the number of bootstrap draws to calculate for each simulated dataset
    # alpha, the nominal significance level for bootstrap inference
    # ncores, the number of cores to use for parallel boostrapping
  # Intermediaries:
    # true_coefs, the oracle model coefficients
    # orig_coefs, the coefficients estimated by MANCOVA/MSGLasso under each simulated dataset
    # pvals, nominal p-valuesfir each effect
    # bca_ints_lower, lower CL for BCA bootstrap intervals
    # bca_ints_upper, upper CL for BCA bootstrap intervals
    # piv_ints_lower, lower CL for pivot bootstrap intervals
    # piv_ints_upper, upper CL for pivot bootstrap intervals
  # Returns:
    # results, a list containing all summary data frames/matrices
  
  #### initialize simulation params and results ####
  q <- sim_sett$q; p <- sim_sett$p; p0 <- sim_sett$p0; n <- sim_sett$n
  k <- sim_sett$k; l <- sim_sett$l; corr <- sim_sett$corr
  
  results <- list(
    # "data_genSeeds"  = (1:nsim)^2 + (1:nsim),
    "true_coefs"     = matrix(nrow = k * p, ncol = nsim),
    "orig_coefs"     = matrix(nrow = k * (p + q + 1), ncol = nsim),
    # "Boot_Draw_Seed" = ((1:nB) + 1)^2,
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

#### OMEI Application ####
# Sample code for fitting a single OMEI model with the proposed multivariate mediation model
Ca_Saliva_OHsurf <- readxl::read_xlsx(paste0(path2,"/OMEIMediationV2_datasets.xlsx"), sheet = "Ca_Saliva_OHsurf")

mod1_confounders <- Ca_Saliva_OHsurf[,2:30] %>% as.matrix()
mod1_trt <- Ca_Saliva_OHsurf$saliva_genus_trt
mod1_mediators <- Ca_Saliva_OHsurf[,32:92] %>% as.matrix()
mod1_outcomes <- Ca_Saliva_OHsurf[,93:97] %>% as.matrix()

q <- 1
p <- ncol(mod1_mediators)
l <- ncol(mod1_confounders)
k <- ncol(mod1_outcomes)
n <- length(mod1_trt)

#### Treatment --> Mediators
mod1_TrtToMedi <- lm(mod1_mediators ~ mod1_trt + mod1_confounders)
mod1_TrtToMedi_summ <- broom::tidy(mod1_TrtToMedi)

mod1_Stage1 <- mod1_TrtToMedi_summ[which(mod1_TrtToMedi_summ$term=="mod1_trt"),]

#### Mediators --> Outcomes
P <- p + l + q
Q <- k
G <- p + 2 # Groups on X: p + 2; p mediator singleton groups + 1 confounder group + 1 treatment group
R <- k # Groups on Y: k; all responses are singleton groups

gmax <- 1 # each variable (resp or pred) belongs to only 1 group
cmax <- l # a group contains at most l variables (confounders form largest group)

GarrStarts <- c(0:(p-1), p,  p+l+1)
GarrEnds <-   c(0:(p-1), p+l, p+l+2)

RarrStarts <- c(0:(k-1)); RarrEnds <- c(0:(k-1)) # singleton response groups

tmp_PQgrps <- FindingPQGrps(P = P, Q = Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
PQgrps <- tmp_PQgrps$PQgrps

tmp_grpWts <- Cal_grpWTs(P = P, Q = Q, G, R, gmax, PQgrps)
grpWTs <- tmp_grpWts$grpWTs

tmp_GRgrps <- FindingGRGrps(P = P, Q = Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
GRgrps <- tmp_GRgrps$GRgrps


Pen_L <- matrix(rep(1, P*Q), P, Q, byrow=T)
Pen_L[(p+1):P,] <- 0 # don't penalize the confounders, (p+1):(P-1), or treatment, P

Pen_G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
Pen_G[(G-1):G,] <- 0 # don't penalize confounder group, G-1, or treatment group, G

grp_Norm0 <- matrix(rep(1, G*R), nrow=G, byrow=TRUE)


lam1.v <- seq(1e-3, 0.05, length=20)
lamG.v <- seq(1e-3, 0.05, length=20)

invisible(capture.output(mod1_try.cv <- MSGLasso.cv(X = cbind(mod1_mediators,mod1_confounders,mod1_trt), 
                                                    Y = mod1_outcomes,
                                                    grpWTs, Pen_L, Pen_G,
                                                    PQgrps, GRgrps, lam1.v, lamG.v,
                                                    # grp_Norm = grp_Norm0,
                                                    fold = 5, seed = 7^7)))


MSGLassolam1 <- mod1_try.cv$lams.c[which.min(as.vector(mod1_try.cv$rss.cv))][[1]]$lam1 # 0.04742105
MSGLassolamG <- mod1_try.cv$lams.c[which.min(as.vector(mod1_try.cv$rss.cv))][[1]]$lam3 # 0.05
MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE) 
MSGLassolamG.m[(G-1):G,] <- 0


mod1_Stage2 <- MSGLasso(X.m = cbind(mod1_mediators, mod1_confounders, mod1_trt),
                        Y.m = mod1_outcomes,
                        grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                        grp_Norm0, MSGLassolam1, MSGLassolamG.m)



#### Extract PIDE and DE
mod1_Stage2_selected <- mod1_Stage2$Beta
rownames(mod1_Stage2_selected) <- colnames(cbind(mod1_mediators, mod1_confounders, mod1_trt))
colnames(mod1_Stage2_selected) <- colnames(mod1_outcomes)

mod1_DE <- mod1_Stage2_selected["mod1_trt",]

mod1_stage2_pide <- mod1_Stage2_selected[1:p,]
mod1_stage1_pide <- mod1_Stage1$estimate

mod1_pide_mat <- matrix(nrow = p, ncol = k)
for(i in 1:k) mod1_pide_mat[,i] <-   mod1_stage1_pide*mod1_stage2_pide[,i]

rownames(mod1_pide_mat) <- rownames(mod1_Stage2_selected)[1:p]; colnames(mod1_pide_mat) <- colnames(mod1_Stage2_selected)
mod1_tide_mat <- colSums(mod1_pide_mat)

mod1_origFit_mat <- matrix(nrow = p*k+2*k, ncol = 1)
mod1_origFit_mat[,1] <- c(vec(mod1_pide_mat),unname(mod1_tide_mat),unname(mod1_DE))
rownames(mod1_origFit_mat) <- c(paste0(rownames(mod1_pide_mat),"_ide_resp",rep(1:k, times = rep(p,k))),
                                paste0("TIDE_resp",1:k),
                                paste0("DE_resp", 1:k))
colnames(mod1_origFit_mat) <- "Orig_Est"


#### Bootstrap inference
nB <- 5e3

B_draws <- matrix(nrow = n, ncol = nB)
for(i in 1:nB) {set.seed((i+1)^1); B_draws[,i] <- sample(1:n, n, replace=T)}

mod1_bootRes <- matrix(nrow = p*k+2*k, ncol = nB)
rownames(mod1_bootRes) <- c(paste0(rownames(mod1_pide_mat),"_ide_resp",rep(1:k, times = rep(p,k))),
                            paste0("TIDE_resp",1:k),
                            paste0("DE_resp", 1:k))
colnames(mod1_bootRes) <- paste0("BootDraw",1:nB)

start_time <- proc.time()
for(i in 1:nB){
  ## Stage 1 fit
  mod1_confounders_boot <- mod1_confounders[B_draws[,i],]
  mod1_trt_boot <- mod1_trt[B_draws[,i]]
  mod1_mediators_boot <- mod1_mediators[B_draws[,i],]
  mod1_outcomes_boot <- mod1_outcomes[B_draws[,i],]
  
  mod1_stage1_boot_fit <- lm(mod1_mediators_boot ~ mod1_trt_boot + mod1_confounders_boot)
  mod1_stage1_boot_tidy <- broom::tidy(mod1_stage1_boot_fit)
  mod1_stage1 <- mod1_stage1_boot_tidy[which(mod1_stage1_boot_tidy$term=="mod1_trt_boot"),]
  
  
  ## Stage 2 fit
  mod1_stage2_boot_fit <- MSGLasso(X.m = cbind(mod1_mediators_boot, mod1_confounders_boot, mod1_trt_boot),
                                   Y.m = mod1_outcomes_boot,
                                   grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                                   grp_Norm0, MSGLassolam1, MSGLassolamG.m)
  mod1_stage2_boot_selected <- mod1_stage2_boot_fit$Beta
  rownames(mod1_stage2_boot_selected) <- colnames(cbind(mod1_mediators_boot, mod1_confounders_boot, mod1_trt_boot))
  colnames(mod1_stage2_boot_selected) <- colnames(mod1_outcomes_boot)
  
  ## PIDE and DE calculation
  mod1_boot_DE <- unname(mod1_stage2_boot_selected["mod1_trt_boot",])
  
  mod1_stage2_pide_boot <- mod1_stage2_boot_selected[1:p,]
  mod1_stage1_pide_boot <- mod1_Stage1$estimate
  
  mod1_pide_mat_boot <- matrix(nrow = p, ncol = k)
  for(j in 1:k) mod1_pide_mat_boot[,j] <- mod1_stage1_pide_boot*mod1_stage2_pide_boot[,j]
  
  mod1_tide_mat_boot <- colSums(mod1_pide_mat_boot)
  
  mod1_bootRes[,i] <- c(vec(mod1_pide_mat_boot),mod1_tide_mat_boot,mod1_boot_DE)
  
  if(i %% 100 == 0){
    cat(paste0("Done with boot sample ", i, "; ", nB-i, " remaining. Time elapsed: ", 
               round(-1*(start_time[3] - proc.time()[3])/60, 3), " minutes. Apx. ",
               round(((-1*(start_time[3] - proc.time()[3])/60)/i)*(nB-i), 3), " minutes remaining\n"))
  }
}

mod1_bootRes <- ifelse(is.nan(mod1_bootRes),0,mod1_bootRes)

## bootstrap p-values
mod1_boot_pvals <- matrix(nrow=p*k+2*k, ncol = 1)
for(i in 1:nrow(mod1_bootRes)){
  if(mod1_origFit_mat[i,] >= 0){
    mod1_boot_pvals[i,] <- min(2*sum(mod1_bootRes[i,]>=2*mod1_origFit_mat[i,],na.rm=T)/nB,1)
  } else{
    mod1_boot_pvals[i,] <- min(2*sum(mod1_bootRes[i,]<2*mod1_origFit_mat[i,],na.rm=T)/nB,1)
  }
}


## BCa confidence intervals
mod1_boot_bcaCI <- t(apply(mod1_bootRes, 1, coxed::bca))
mod1_boot_bcaCI[is.nan(mod1_boot_bcaCI)] <- 0

## pivotal bootstrap CIs
mod1_boot_pivCI <- 2*matrix(data = c(mod1_origFit_mat, mod1_origFit_mat), ncol = 2) - t(apply(mod1_bootRes, 1, quantile, c(1-alpha/2,alpha/2)))

## Summary table
mod1_boot_summ <- cbind(mod1_origFit_mat, apply(mod1_bootRes,1,mean), apply(mod1_bootRes,1,sd),
                        mod1_boot_pvals, mod1_boot_bcaCI, mod1_boot_pivCI)
colnames(mod1_boot_summ) <- c("OrigEst", "Mean_boot", "boot_SE", "boot_pval",
                              "bca_lowerCL", "bca_upperCL",
                              "piv_lowerCL", "piv_upperCL")
mod1_boot_summ[is.nan(mod1_boot_summ)] <- 0

save(mod1_boot_summ, file = paste0(path3, "mod1_boot_summ.Rda"))


#### Model results
load(paste0(path3,"mod1_boot_summ.Rda")); p <- 61

## PI: TIDE = 0.124229, DE = 0.2897, TE = 0.4140
mod1_resp1_bootSumm <- mod1_boot_summ[c(1:p,306,311),]           
mod1_resp1_bootSumm[which(mod1_resp1_bootSumm[,1]!=0),] %>% round(4)

## ICDAS: TIDE = 0.0099, DE = 0.2630, TE = 0.2729
mod1_resp2_bootSumm <- mod1_boot_summ[c((p+1):(2*p),307,312),]  
mod1_resp2_bootSumm[which(mod1_resp2_bootSumm[,1]!=0),] %>% round(4)

## DS: TIDE = 0.0117, DE = 0.8351, TE = 0.8467
mod1_resp3_bootSumm <- mod1_boot_summ[c((2*p+1):(3*p),308,313),] 
mod1_resp3_bootSumm[which(mod1_resp3_bootSumm[,1]!=0),] %>% round(4)

## MS: TIDE = -0.6689, DE = 0.8351, TE = 0.1662
mod1_resp4_bootSumm <- mod1_boot_summ[c((3*p+1):(4*p),309,314),] 
mod1_resp4_bootSumm[which(mod1_resp4_bootSumm[,1]!=0),] %>% round(4)

## FS: TIDE = -0.0356, DE = -0.1985, TE = -0.2341
mod1_resp5_bootSumm <- mod1_boot_summ[c((4*p+1):(5*p),310,315),] 
mod1_resp5_bootSumm[which(mod1_resp5_bootSumm[,1]!=0),] %>% round(4)

mod1_bootSumm_allRes <- rbind(mod1_resp1_bootSumm,mod1_resp2_bootSumm,
                              mod1_resp3_bootSumm,mod1_resp4_bootSumm,
                              mod1_resp5_bootSumm)

# to_latex_tab(mod1_bootSumm_allRes)

## Total Effects
mod1_TE <- c("TE_PI"    = sum(mod1_resp1_bootSumm[1:p,1])+mod1_resp1_bootSumm[p+1,1],
             "TE_ICDAS" = sum(mod1_resp2_bootSumm[1:p,1])+mod1_resp2_bootSumm[p+1,1],
             "TE_DS"    = sum(mod1_resp3_bootSumm[1:p,1])+mod1_resp3_bootSumm[p+1,1],
             "TE_MS"    = sum(mod1_resp4_bootSumm[1:p,1])+mod1_resp4_bootSumm[p+1,1],
             "TE_FS"    = sum(mod1_resp5_bootSumm[1:p,1])+mod1_resp5_bootSumm[p+1,1])