# Compute column j of Theta \approx solve(t(X) %*% X)
theta_column <- function(X, j){
  p <- ncol(X); n <- nrow(X)
  
  y_j <- X[, j]
  X_minus_j <- X[, -j, drop = FALSE]
  
  cv_fit <- glmnet::cv.glmnet(X_minus_j, y_j, intercept = FALSE, 
                              standardize = FALSE, nfolds = 5)
  lambda_opt <- cv_fit$lambda.min
  
  lasso_fit <- glmnet::glmnet(X_minus_j, y_j, lambda = lambda_opt, 
                              intercept = FALSE, standardize = FALSE)
  gamma_j <- as.numeric(lasso_fit$beta)
  
  residuals <- y_j - X_minus_j %*% gamma_j
  tau_sq_j <- mean(residuals^2) + lambda_opt * sum(abs(gamma_j))
  
  theta_col <- numeric(p)
  theta_col[j] <- 1 / tau_sq_j
  
  # Fill off-diagonal elements
  k_idx <- 1
  for (k in 1:p) {
    if (k != j) {
      theta_col[k] <- -gamma_j[k_idx] / tau_sq_j
      k_idx <- k_idx + 1
    }
  }
  
  return(list(theta_col = theta_col, tau_sq = tau_sq_j))
}

# Compute, in serial, Theta \approx solve(t(X) %*% X) column-wise
theta_calc <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  X_std <- scale(X)
  X_std[is.na(X_std)] <- 0
  
  Theta <- matrix(0, p, p)
  tau_sq <- numeric(p)
  
  for (j in 1:p) {
    result <- theta_column(X_std, j, lambda, method)
    Theta[, j] <- result$theta_col
    tau_sq[j] <- result$tau_sq
  }
  
  return(Theta)
}

# Compute, in parallel, Theta \approx solve(t(X) %*% X)
theta_calc_parallel <- function(X, cores = parallel::detectCores(), folds = 5){
  n <- nrow(X)
  res_mat <- matrix(nrow = ncol(X), ncol = ncol(X)+2)
  
  clust <- makeCluster(cores)
  registerDoParallel(clust)
  
  res_mat <- foreach(j = 1:ncol(X), .combine = "rbind", .packages = "glmnet") %dopar% {
    lasso_j <- cv.glmnet(x = X[,-j], y = X[,j], nfolds = folds, intercept = F, standardize = F)
    
    # gamma_j <- lasso_j$glmnet.fit$beta[,lasso_j$index[2,]]
    
    Chat_j <- numeric(ncol(X))
    Chat_j[-j] <- lasso_j$glmnet.fit$beta[,lasso_j$index[2,]]
    Chat_j[j] <- 1
    
    lambda_1se <- lasso_j$lambda.1se
    gamma_j <- lasso_j$glmnet.fit$beta[,lasso_j$index[2,]]
    
    tau2_j <- norm(X[,j]-(X[,-j] %*% as.matrix(gamma_j)), "2")/n + lambda_1se * norm(as.matrix(gamma_j), "1")
    
    c(lambda_1se, tau2_j, Chat_j)
  }
  
  stopCluster(clust)
  colnames(res_mat) <- c("lambda.1se", "tau2", paste0("gamma", 1:(ncol(res_mat)-2)))
  rownames(res_mat) <- colnames(X)
  
  tau2_mat <- res_mat[,"tau2"]
  lambda_mat <- res_mat[,"lambda.1se"]  
  
  That2_inv <- diag(1/tau2_mat)
  Chat_mat <- res_mat[,-which(colnames(res_mat)%in%c("lambda.1se", "tau2"))]
  
  Theta <- That2_inv %*% Chat_mat
  
  return(unname(Theta))
}
