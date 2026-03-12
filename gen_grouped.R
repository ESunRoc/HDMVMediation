gen_data <- function(n, p, k, l, p0, cov_groups = NULL, resp_groups = NULL, coef.seed = 823543, data.seed = 823543){
  
  ##### Generate fixed coefficients (all use the same seed)  #####
  set.seed(coef.seed)
  
  # True mediators are the first p0
  true_idx <- 1:p0
  
  # alpha: 1 x p matrix
  alpha <- matrix(0, nrow = 1, ncol = p)
  if(p0 > 0){
    alpha[1, true_idx] <- sample(c(-5, 5), p0, replace = TRUE)
  }
  if(p - p0 > 0){
    alpha[1, -true_idx] <- sample(c(-0.05, 0.05), p - p0, replace = TRUE)
  }
  
  # tau: 1 x q matrix
  tau <- matrix(sample(c(-1, 1), k, replace = TRUE), nrow = 1, ncol = k)
  
  # xi: l x p matrix
  xi <- matrix(runif(l * p, -2, 2), nrow = l, ncol = p)
  
  # eta: l x q matrix
  eta <- matrix(runif(l * k, -1, 1), nrow = l, ncol = k)
  
  # phi: p x q matrix (all zero except for true mediators)
  phi <- matrix(0, nrow = p, ncol = k)
  if(p0 > 0){
    set.seed(coef.seed + 1)   # separate seed for phi to avoid unwanted correlation
    for(i in true_idx){
      phi[i, ] <- sample(c(-5, -4, -3, -2, 2, 3, 4, 5), k, replace = TRUE)
    }
  }
  
  ##### Construct covariance matrix for mediator errors (Sigma_M) based on cov_groups #####
  # Default: independent errors (identity)
  Sigma_M <- diag(p)
  
  if(!is.null(cov_groups)){
    # Assume cov_groups is a list of index vectors; groups are disjoint
    rho <- 0.5   # fixed intra‑group correlation (can be changed if needed)
    for(grp in cov_groups){
      if(length(grp) > 1){
        # Set off‑diagonals within group to rho
        for(i in grp){
          for(j in grp){
            if(i != j) Sigma_M[i, j] <- rho
          }
        }
      }
    }
  }
  
  #### Compute required error variance for Y to satisfy correlation bound ####
  # Variance of A (Bernoulli(0.5))
  varA <- 0.25
  
  # Variance of U (confounders are independent)
  varU <- numeric(l)
  if(l >= 1) varU[1] <- 0.25 * 0.75          # Bern(0.25)
  if(l >= 2) varU[2] <- 0.5 * 0.5            # Bern(0.5)
  if(l >= 3) varU[3] <- 0.75 * 0.25          # Bern(0.75)
  if(l >= 4) varU[4] <- 15                   # Pois(15)
  if(l >= 5) varU[5] <- 1 / 12               # Unif(0,1)
  if(l > 5)  varU[6:l] <- 1                  # extra N(0,1)
  
  VarU <- diag(varU, nrow = l, ncol = l)
  
  # Compute beta_A = alpha phi + tau   (1 x q)
  beta_A <- alpha %*% phi + tau
  
  # Compute beta_U = xi phi + eta   (l x q)
  beta_U <- xi %*% phi + eta
  
  # Σ_sys = varA * (beta_A' beta_A) + beta_U' VarU beta_U + phi' Sigma_M phi
  term1 <- varA * crossprod(beta_A)                # q x q
  term2 <- t(beta_U) %*% VarU %*% beta_U           # q x q
  term3 <- t(phi) %*% Sigma_M %*% phi              # q x q (contribution from M errors)
  Sigma_sys <- term1 + term2 + term3
  
  # Function to compute maximum absolute off‑diagonal correlation for a given σ²
  max_cor <- function(s2){
    diag_vals <- diag(Sigma_sys) + s2
    cor_mat <- Sigma_sys / sqrt(outer(diag_vals, diag_vals))
    diag(cor_mat) <- 0
    max(abs(cor_mat))
  }
  
  # Find smallest common σ² such that max |cor| ≤ 0.3
  if(max_cor(0) <= 0.3){
    sigma2_Y <- 0
  } else{
    # Upper bound search
    upper <- 100 * max(diag(Sigma_sys))
    while(max_cor(upper) > 0.3){
      upper <- upper * 2
    }
    lower <- 0
    for(iter in 1:50){
      mid <- (lower + upper) / 2
      if(max_cor(mid) <= 0.3){
        upper <- mid
      } else{
        lower <- mid
      }
    }
    sigma2_Y <- upper
  }
  
  ##### Define data generator #####
  gen_samp_dat <- function(data.seed){
    set.seed(data.seed)
    
    # Generate confounders U as an n x l matrix
    U <- matrix(0, nrow = n, ncol = l)
    U[, 1] <- rbinom(n, 1, 0.25); U[, 2] <- rbinom(n, 1, 0.5)
    U[, 3] <- rbinom(n, 1, 0.75); U[, 4] <- rpois(n, 15)
    U[, 5] <- runif(n, 0, 1)
    if(l > 5){
      for(j in 6:l) U[, j] <- rnorm(n, 0, 1)
    }
    
    # Exposure A as an n x 1 matrix (Bernoulli 0.5)
    A <- matrix(rbinom(n, 1, 0.5), ncol = 1)
    
    # Generate mediator errors ε_M with covariance Sigma_M
    # Use Cholesky decomposition for efficiency
    chol_M <- chol(Sigma_M)
    epsilon_M <- matrix(rnorm(n * p), nrow = n) %*% chol_M   # n x p
    
    # Mediators M: M = A alpha + U xi + ε_M
    M <- A %*% alpha + U %*% xi + epsilon_M
    
    # Outcomes Y: Y = M phi + A tau + U eta + ε_Y, ε_Y ~ N(0, σ² I)
    Y <- M %*% phi + A %*% tau + U %*% eta +
      matrix(rnorm(n * k, 0, sqrt(sigma2_Y)), n, k)
    
    colnames(M) <- paste0("M", 1:p)
    # colnames(A) <- "trt"
    colnames(U) <- paste0("U", 1:l)
    colnames(Y) <- paste0("Y", 1:k)
    
    list(M = M, A = A, U = U, Y = Y,
         true_mediators = true_idx, alpha = alpha,
         tau = tau, xi = xi, eta = eta, phi = phi,
         sigma2_Y = sigma2_Y, Sigma_M = Sigma_M,
         cov_groups = cov_groups,
         resp_groups = resp_groups)
  }
  
  ##### Return parameters and generator #####
  list(fixed_params = list(alpha = alpha, tau = tau, xi = xi,
                           eta = eta, phi = phi, sigma2_Y = sigma2_Y,
                           Sigma_M = Sigma_M,
                           true_mediators = true_idx,
                           cov_groups = cov_groups,
                           resp_groups = resp_groups),
       generate = gen_samp_dat)
}-