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
  
  #### beta (p x k), PIDE (mediators --> outcomes) ####
  # if(p == 25){
  #   if(quiet == T){
  #     # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
  #     suppressWarnings({
  #       # beta <- rbind(rep(c(-4, -2,  3,  5, 6), k/5),        
  #       #               rep(c(-6,  5, -2, -3, 3), k/5),
  #       #               rep(c(-3,  0, -4,  4, 2), k/5),
  #       #               rep(c(-4,  7,  0,  3, 4), k/5),
  #       #               rep(c( 0, -3,  2, -3, 0), k/5),
  #       #               matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  #       beta <- rbind(rep(c(-4, -2,  3,  5, 6), k/5),        
  #                     rep(c(-6,  5, -2, -3, 3), k/5),
  #                     rep(c(-3,  8, -4,  4, 2), k/5),
  #                     rep(c(-4,  7, -8,  3, 4), k/5),
  #                     rep(c(-7, -3,  2, -3, 8), k/5),
  #                     matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  #     })
  #   } else {
  #     # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
  #     beta <- rbind(rep(c(-4, -2,  3,  5, 6), k/5),        
  #                   rep(c(-6,  5, -2, -3, 3), k/5),
  #                   rep(c(-3,  0, -4,  4, 2), k/5),
  #                   rep(c(-4,  7,  0,  3, 4), k/5),
  #                   rep(c( 0, -3,  2, -3, 0), k/5),
  #                   matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  #   }
  # } else{
  #   if(quiet == T){
  #     # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
  #     suppressWarnings({
  #       beta <- rbind(rep(c(-4, -2,  3,  5, 6), k/5),        
  #                     rep(c(-6,  5, -2, -3, 3), k/5),
  #                     rep(c(-3,  7, -4,  4, 2), k/5),
  #                     rep(c(-4,  7,  8,  3, 4), k/5),
  #                     rep(c(-2,  0, -4,  4, 2), k/5),
  #                     rep(c(-5,  7,  0,  6, 4), k/5),
  #                     rep(c( 0, -4,  3, -5, 0), k/5),
  #                     rep(c( 0,  0,  6, -7, 5), k/5),
  #                     rep(c( 4, -6,  0, -1, 0), k/5),
  #                     rep(c( 3, -2,  0,  0, 6), k/5),
  #                     matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  #     })
  #   } else {
  #     # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
  #     beta <- rbind(rep(c(-4, -2,  3,  5, 6), k/5),        
  #                   rep(c(-6,  5, -2, -3, 3), k/5),
  #                   rep(c(-3,  7, -4,  4, 2), k/5),
  #                   rep(c(-4,  7,  8,  3, 4), k/5),
  #                   rep(c(-2,  0, -4,  4, 2), k/5),
  #                   rep(c(-5,  7,  0,  6, 4), k/5),
  #                   rep(c( 0, -4,  3, -5, 0), k/5),
  #                   rep(c( 0,  0,  6, -7, 5), k/5),
  #                   rep(c( 4, -6,  0, -1, 0), k/5),
  #                   rep(c( 3, -2,  0,  0, 6), k/5),
  #                   matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  #   }
  # }
  
  if(p0 == 5){
    if(quiet == T){
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      suppressWarnings({
        beta <-      b_scale * rbind(rep(c(-5, -3,  4,  6, 7), k/5),        
                               rep(c(-7,  6, -3, -4, 4), k/5),
                               rep(c(-4,  9, -5,  5, 3), k/5),
                               rep(c(-5,  8, -9,  4, 5), k/5),
                               rep(c(-8, -4,  3, -4, 9), k/5),
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
  set.seed(6^6)
  CorrMat_E2 <- gencor::gencor(d = k, method = corr, lim_low = 0.3, lim_medium = 0.7)$Matrix
  
  Sigma_E2 <- MBESS::cor2cov(cor.mat = CorrMat_E2, sd = rep(1,k))
  set.seed(seed)
  # E2 <- mvrnorm(n=n, mu = rep(0,k), Sigma = Sigma_E2, empirical = T)
  # E2 <- mvrnorm(n=n, mu = rep(0,k), Sigma = Sigma_E2, empirical = F)
  E2 <- mvnfast::rmvn(n=n, mu = rep(0,k), sigma = Sigma_E2)
  
  #### Y_mat (n x k), response matrix ####
  ## Sample from the matrix normal, conditional on Z_mat
  # Y <- MBSP::matrix_normal(M = Z_mat %*% B_mat, U = diag(n), V = Sigma_E2)
  
  ## Generate from the model equation
  Y <- Z_mat %*% B_mat + E2
  # Y <- M%*%beta + U%*%alpha2 + X%*%
  
  ## sample directly from the conditional distribution -- this is better but takes several minutes for large n,k
  # Sigma <- diag(n) %x% cor(E2)
  # mu <- (diag(n) %x% t(B_mat))%*%vec(t(Z_mat))
  # Y <- mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = T)
  # Y <- Y[,1:k]
  
  colnames(Y) <- paste0("Y", 1:k) 
  
  #### Define returns ####
  test_data <- list("X" = X, "U" = U, "M" = M, "Y" = Y)
  sim_params <- list("True B_mat" = B_mat, "True A_mat" = A_mat,
                     "E1" = E1, "E2" = E2, "ResponseCor" = CorrMat_E2)
  full_return <- list("data" = test_data, "parameters" = sim_params)
  return(full_return)
}

gen_data_v2 <- function(n, p, p0, l, q, k, seed = seed, corr = c("low", "medium", "high"), quiet = T, bscale = 1){
  #### Function information: ####
  # Requires: 
  # n, sample size
  # p, num. candidate mediators
  # p0, number of true mediators
  # l, num. of confounders
  # q, num. of treatments
  # k, num. of responses
  # corr, amount of correlation for response residuals
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
  rownames(xi1) <- paste0("xi1_", 1:q)
  
  #### A_mat (l+q x p) ####
  A_mat <- rbind(alpha1, xi1)
  
  #### E1 (n x p), Mediator error matrix ####
  set.seed(seed)
  E1 <- mvrnorm(n=n, mu = rep(0,p), Sigma = diag(p))
  
  
  #### Mediators (n x p) ####
  M <- scale(U%*%alpha1 + X%*%xi1 + E1) # scale(W_mat%*%A_mat + E1)
  colnames(M) <- paste0("M", 1:p)
  
  # tmp <- lm(M ~ W_mat - 1); tmp_tidy <- broom::tidy(tmp); tmp_tidy[which(tmp_tidy$term=="W_matX1"),] %>% View
  
  #### Z_mat (n x q+p+l) ####
  Z_mat <- cbind(M, U, X)
  
  #### beta (p x k), PIDE (mediators --> outcomes) ####
  if(p == 25){
    if(quiet == T){
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      suppressWarnings({
        beta <- bscale * rbind(rep(c(-5, -3,  4,  6, 7), k/5),        
                        rep(c(-7,  6, -3, -4, 4), k/5),
                        rep(c(-4,  9, -5,  5, 3), k/5),
                        rep(c(-5,  8, -9,  4, 5), k/5),
                        rep(c(-8, -4,  3, -4, 9), k/5),
                        matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
      })
    } else {
      # set.seed(3^3) # rbind(beta_1,...,beta_p0, beta_(p0+1),...,beta_p)
      beta <- bscale * rbind(rep(c(-5, -3,  4,  6, 7), k/5),        
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
        beta <- bscale * rbind(rep(c(-6, -4,  5,  7, 7), k/5),        
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
      beta <- bscale * rbind(rep(c(-6, -4,  5,  7, 7), k/5),        
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
  rownames(B_mat) <- c(paste0("beta_", 1:p),paste0("alpha2_",1:l),paste0("xi2_",1:q))
  
  #### E2 (n x k), Mediator error matrix ####
  set.seed(6^6); CorrMat_E2 <- gencor::gencor(d = k, method = corr, lim_low = 0.3, lim_medium = 0.7)$Matrix
  
  ### Ideally, empirical = T, but this was causing problems (see notes)
  Sigma_E2 <- MBESS::cor2cov(cor.mat = CorrMat_E2, sd = rep(1,k))
  set.seed(seed); E2_tmp <- mvrnorm(n=n, mu = rep(0,k), Sigma = Sigma_E2, empirical = F)
  # set.seed(seed); E2 <- mvnfast::rmvn(n=n, mu = rep(0,k), sigma = Sigma_E2)
  
  #### Y_mat (n x k), response matrix ####
  ## Sample from the matrix normal, conditional on Z_mat
  # Y <- MBSP::matrix_normal(M = Z_mat %*% B_mat, U = diag(n), V = Sigma_E2)
  
  ## Generate from the model equation
  Y1 <- Z_mat %*% B_mat + E2
  # Y2 <- M%*%beta + U%*%alpha2 + X%*%xi2 + E2
  
  ## sample directly from the conditional distribution -- this is better but takes several minutes for large n,k
  # Sigma <- diag(n) %x% cor(E2)
  # mu <- (diag(n) %x% t(B_mat))%*%vec(t(Z_mat))
  # Y <- mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = T)
  # Y <- Y[,1:k]
  
  colnames(Y) <- paste0("Y", 1:k) 
  
  #### Define returns ####
  test_data <- list("X" = X, "U" = U, "M" = M, "Y" = Y)
  sim_params <- list("True B_mat" = B_mat, "True A_mat" = A_mat,
                     "E1" = E1, "E2" = E2, "ResponseCor" = CorrMat_E2)
  full_return <- list("data" = test_data, "parameters" = sim_params)
  return(full_return)
}

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

run_sim <- function(sim_sett, nsim, dataGenerator = c("v1","v2"), b_scale = 1, drop_penn = F){
  #### Function information: ####
  # Requires: 
    # sim_sett, a list containing n, p, p0, l, k, q, corr
    # nsim, the number of data sets to simulate
  # Intermediaries:
    # all_dataAndparams, a list containing all generated data sets and their true parameters
    # multivar_results, a list containing all step 1 models and multivariate step 2 models
    # univar_results, a list containing all step 1 models and univariate step 2 models
  # Returns:
    # Results, a list containing: all generated data sets, all step 1 model fits, all step 2 model fits
  
  dataGenerator <- match.arg(dataGenerator)
  
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
  
  # it breaks if I try anything else.
  assign("Pen_L", Pen_L, envir=.GlobalEnv); assign("Pen_G", Pen_G, envir=.GlobalEnv)
  
  start_time <- proc.time()
  
  #### Repeat simulation nsim times ####
  for(sim in 1:nsim){
    #### Set up sim parameters ####
    seed <- sim # 1 to nsim
    if(dataGenerator == "v1") simdata <- gen_data(n = n, p = p, p0 = p0, l=l, k = k, q = q, seed = seed, corr = corr, b_scale = b_scale)
    if(dataGenerator == "v2") simdata <- gen_data_v2(n = n, p = p, p0 = p0, l=l, k = k, q = q, seed = seed, corr = corr, b_scale = b_scale)
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
    
    # lam1.v <- seq(1e-2, 0.15, length=25)
    # lamG.v <- seq(1e-2, 0.15, length=25)
    
    
    
    if(drop_penn == T){
      lam1.v <- seq(0.0125, 0.15, length=25)
      lamG.v <- seq(0.0125, 0.15, length=25)
    } else{
      lam1.v <- seq(0.025, 0.15, length=25)
      lamG.v <- seq(0.025, 0.15, length=25)
    }
    
    # this prints (lam1.v,lamg.v) counter for each fold. This is suppressed so printing to console doesn't take forever
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
    
    # cat(paste0("Done with simulation ", sim, "; ", nsim-sim, 
    #            " remaining. Time elapsed: ", 
    #            round(-1*(start_time[3] - proc.time()[3])/60, 3), " minutes\n"))
    
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
