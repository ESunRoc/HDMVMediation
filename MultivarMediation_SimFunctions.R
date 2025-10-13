gen_data <- function(n, p, p0, l, q, k, seed = seed, corr = c("low", "medium", "high"), b_scale = 1){
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
  
  
  #### Z_mat (n x q+p+l) ####
  Z_mat <- cbind(M, U, X)
  
  #### beta (p x k), PIDE (mediators --> outcomes) ####
  if(p0 == 5 & k == 5 & corr == "low"){
    beta <- b_scale * rbind(c( 4,  2,  3, -1,  5),
                            c( 3, -5,  1,  4, -2),
                            c( 2,  4,  5,  4, -1),
                            c( 5,  1,  3, -3,  4),
                            c(-1, -3, -4,  2,  4),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 5 & k == 5 & corr == "medium"){
    beta <- b_scale * rbind(c( 4,  2,  3,  2,  5),
                            c( 3, -5, -1,  4, -2),
                            c( 2,  4, -5, -3,  3),
                            c(-5,  1, -2,  5,  4),
                            c( 1,  3, -4, -2, -5),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 5 & k == 5 & corr == "high"){
    beta <- b_scale * rbind(c( 4,  4,  3,  4,  5),
                            c( 3,  5,  5,  4,  3),
                            c(-2, -4, -5, -3, -1),
                            c(-5,  4, -2,  5,  4),
                            c( 1,  3, -4, -3, -5),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 5 & k == 10 & corr == "low"){
    beta <- b_scale * rbind(c( 3, -2,  4,  1,  5,  2, -3,  4,  1,  3),
                            c( 1,  5,  1, -3,  2, -4,  2, -5,  3,  5),
                            c( 2, -4, -5,  2,  1, -3,  5,  2, -4,  4),
                            c(-5,  3,  2,  4, -1,  5, -4,  1,  5, -5),
                            c( 4,  1, -3, -5,  3, -1,  1,  3, -2,  3),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 5 & k == 10 & corr == "medium"){
    beta <- b_scale * rbind(c( 4,  2,  3,  1,  5,  1, -5,  5,  3, 4),
                            c( 3, -5,  1,  4, -2,  3,  1,  4, -5, 3),
                            c( 2,  4, -5, -3,  1, -4,  2, -5,  3, 3),
                            c(-5,  1,  2,  5,  4, -2,  5, -3,  4, 2),
                            c( 1,  3, -4, -2, -5, -5,  4,  5, -2, 5),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 5 & k == 10 & corr == "high"){
    beta <- b_scale * rbind(c( 4,  2,  3,  1,  5,  1, -5,  3,  3,  4),
                            c( 3,  5,  3,  4,  2,  3,  1,  4,  5,  3),
                            c(-2, -4, -5, -3, -1, -4, -2, -5, -3, -4),
                            c(-5,  1,  5,  5,  4, -2,  5, -3,  4,  3),
                            c( 2,  3, -4, -4, -5, -5,  4, -3, -3,  5),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 10 & k == 5 & corr == "low"){
    beta <- b_scale * rbind(c( 2,  1, -3,  4, -2),
                            c(-1,  3, -2, -5,  1),
                            c( 4, -2,  1,  3, -4),
                            c(-3,  4, -5,  2,  3),
                            c( 1, -5,  4,  1,  5),
                            c( 5,  2, -1, -4, -3),
                            c(-2,  4,  3,  1,  2),
                            c( 3,  1,  2,  5,  5),
                            c(-4,  5,  5,  3, -1),
                            c( 2, -3, -4,  2,  4), 
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 10 & k == 5 & corr == "medium"){
    beta <- b_scale * rbind(c( 2,  1, -3,  4, -2),
                            c(-1, -3, -4, -5,  1),
                            c( 3, -4,  1,  3, -4),
                            c( 3,  4,  3, -1,  3),
                            c( 1,  5,  2,  1,  5),
                            c( 3,  2, -1, -4, -3),
                            c(-2,  4,  3,  4,  2),
                            c(-3,  1,  2, -3,  5),
                            c(-4,  5,  5,  3, -3),
                            c( 2, -3, -4,  2,  4), 
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 10 & k == 5 & corr == "high"){
    beta <- b_scale * rbind(c( 2,  1, -3,  4, -2),
                            c(-1, -3, -2, -5,  1),
                            c( 4, -2,  1,  3, -4),
                            c(-3,  4, -5,  2,  3),
                            c( 1, -5,  4,  1, -5),
                            c( 5,  2, -1, -4, -3),
                            c(-2,  4,  3,  1,  2),
                            c( 3,  1,  2, -5,  5),
                            c(-4,  5,  5,  3, -1),
                            c( 2, -3, -4,  2,  4), 
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 10 & k == 10 & corr == "low"){
    beta <- b_scale * rbind(c( 3, -2, -1,  3, -2,  2, -1,  1, -2, -4),
                            c(-4, -2,  4,  2,  3, -1,  2,  5,  3,  3),
                            c(-2, -3, -4,  1,  2,  5, -3,  1, -3, -5),
                            c( 1, -3, -2,  3,  4, -3, -5, -2,  1,  4),
                            c( 5,  4, -3, -5, -1,  2, -4,  3,  5, -2),
                            c(-2,  1,  5,  5, -3,  4,  1, -4, -2,  3),
                            c( 1, -3,  4,  1,  5, -2,  3,  2, -5,  2),
                            c( 5,  2, -1, -3,  4,  1, -2,  3,  1,  4),
                            c( 4, -4,  3, -3,  2, -5,  2, -1,  3,  5),
                            c( 2,  3, -2, -4,  1,  3, -5,  4, -4,  2),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else if(p0 == 10 & k == 10 & corr == "medium"){
    beta <- b_scale * rbind(c( 3, -2,  1,  4, -5,  2, -1,  1,  3,  2),
                            c(-4,  1,  5, -2,  3, -1,  4, -5,  2,  3),
                            c( 2, -3,  4,  1, -2,  5, -3,  1,  4,  5),
                            c( 1, -5, -2,  3,  4, -3, -5, -2, -1,  4),
                            c( 5,  4,  3, -5,  1,  2, -4,  3,  5, -1),
                            c(-2,  1, -5,  5, -3,  4,  1, -4, -2,  3),
                            c( 1, -3,  4, -1,  5, -2,  3,  2, -5,  1),
                            c(-5,  2, -1,  5,  4,  1, -2,  3,  1,  4),
                            c( 4, -4,  3, -3,  2, -5,  2, -1,  3,  5),
                            c( 3,  5, -2, -4, -1,  3, -5,  4, -4,  2),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } else{ # if(p0 == 10 & k == 10 & corr == "high")
    beta <- b_scale * rbind(c( 3, -2,  1,  4, -5,  2, -1,  5, -3,  2),
                            c(-4,  1,  5, -2,  3, -1,  4, -5,  2, -3),
                            c( 2,  3, -4,  1, -2,  5, -3,  1,  4, -5),
                            c(-1, -5,  2,  3,  4, -3,  5, -2, -1,  4),
                            c( 5,  4, -3, -5,  1,  2, -4,  3,  5, -1),
                            c(-2,  1, -5,  2, -3,  4,  1, -4, -2,  3),
                            c( 1, -3,  4, -1,  5, -2,  3,  2, -5,  1),
                            c(-5,  2, -1,  5, -4,  1, -2,  3,  1, -4),
                            c( 4, -4,  3, -3,  2, -5,  2, -1,  3,  5),
                            c(-3,  5, -2, -4, -1,  3, -5,  4, -4, -2),
                            matrix(rep(0,k*(p-p0)), nrow = p-p0, ncol = k)) %>% as.matrix
  } 
  rownames(beta) <- paste0("beta_", 1:p)
  
  #### alpha2 (l x k), nuisance param(s) ####
  
  # set.seed(4^4); alpha2 <- matrix(runif(l*k, -2, 2), nrow = l, ncol = k) # eta
  # alpha2 <- matrix(rep(seq(-2, 2, by = 0.2),5, byrow=T)[1:(l*k)], nrow = l, ncol = k)
  if(k==5){
    alpha2 <- matrix(c(-0.61, -0.61, -0.03, -0.27, -0.03,
                       0.35, -0.02,  0.57,  0.11, -0.56,
                       -0.16,  0.44, -0.54,  0.21, -0.33,
                       0.64, -0.39, -0.38,  0.02,  0.38,
                       -0.23,  0.58,  0.37, -0.07,  0.54), byrow = T, nrow = l, ncol = k)
  } else{
    alpha2 <- matrix(c( 0.34, -0.45,  0.56, -0.12,  0.78, -0.45,  0.12, -0.78,  0.45, -0.12,
                        -0.78,  0.12, -0.45,  0.78, -0.12,  0.56, -0.78,  0.12, -0.45,  0.78,
                        0.12, -0.78,  0.12, -0.45,  0.78, -0.12,  0.56, -0.45,  0.78, -0.12,
                        -0.45,  0.56, -0.78,  0.12, -0.45,  0.78, -0.12,  0.56, -0.12,  0.45,
                        0.56, -0.12,  0.45, -0.56,  0.12, -0.34,  0.45, -0.56,  0.12, -0.45), byrow = T, nrow = l, ncol = k)
  }
  
  
  #### xi2 (q x k), nuisance param(s) ####
  set.seed(5^5); xi2 <- matrix(sample(c(-1,1), q*k, replace=T), nrow = q, ncol = k) # tau
  # xi2 <- matrix # set manually
  
  
  # Y = M*beta + X*tau + U*eta + error
  # Y = Z*B + error
  # cor(Y) = cor(error)
  #### B_mat (p+l+q x k) ####
  B_mat <- rbind(beta, alpha2, xi2)
  
  #### E1 (n x k), Mediator error matrix ####
  set.seed(5^5); CorrMat_E2 <- gencor::gencor(d = k, method = corr, lim_low = 0.3, lim_medium = 0.7)$Matrix
  
  # Sigma_E2 <- MBESS::cor2cov(cor.mat = CorrMat_E2, sd = rep(1,k))
  # Sigma_E2 <- matrix(rep(1,k),ncol=1) %*% matrix(rep(1,k),nrow=1) * CorrMat_E2
  
  set.seed(seed); E2 <- mvnfast::rmvn(n=n, mu = rep(0,k), sigma = CorrMat_E2)
  
  #### Y_mat (n x k), response matrix ####
  ## Sample from the matrix normal, conditional on Z_mat
  # Y <- MBSP::matrix_normal(M = Z_mat %*% B_mat, U = diag(n), V = CorrMat_E2)
  
  ## Generate from the model equation
  Y <- Z_mat%*%B_mat + E2
  # Y <- matrix(nrow = n, ncol = 5)
  # for(i in 1:n){
  #   Y[i,]<-MASS::mvrnorm(n=1, mu = t((Z_mat %*% B_mat)[i,]), Sigma = CorrMat_E2, empirical = F)
  # }
  
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
      lam1.v <- seq(0.0075, 0.075, length=10)
      lamG.v <- seq(0.0075, 0.075, length=10)
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
                          seed = sim^2 + sim, corr = corr, b_scale = 1)
    # seed = sim^2 + sim, corr = corr, quiet = T, b_scale = 2)
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