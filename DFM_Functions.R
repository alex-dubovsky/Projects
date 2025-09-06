
########################################### Dynamic Factor Models ##########################################
#' @import data.table
#' @importFrom data.table :=

#' @export
Estimate_DFM <- function(X, n_f = NULL, lags_f = NULL, lags_v = NULL, max_EV = NULL, undo_scale=TRUE, factor_method="pc", VAR_method="OLS") {
  if (is.null(n_f)) {
    n_f <- IC_BaiNg(X)
  }
  if(factor_method=="pc"){
    SFM <- pc(X, undo_scale=undo_scale)
  }else if(factor_method=="sWF"){
    SFM <- sWF(X, n_f=n_f, undo_scale=undo_scale) 
  }else{
    stop("invalid factor_method")
  }
  if(undo_scale){
    v <- scale(X, center=TRUE, scale=FALSE) - SFM$F_hat[, 1:n_f] %*% t(SFM$Lambda_hat[, 1:n_f])
  }else{
    v <- scale(X, center=TRUE, scale=TRUE) - SFM$F_hat[, 1:n_f] %*% t(SFM$Lambda_hat[, 1:n_f])
  }
  if(VAR_method=="OLS"){
    VAR_factors <- VAR_est(SFM$F_hat[, 1:n_f], p = lags_f)
  }else if(VAR_method=="lasso"){
    VAR_factors <- lasso_VAR_est(SFM$F_hat[, 1:n_f], p = lags_f)
  }else if(VAR_method=="sparse_lasso"){
    VAR_factors <- sparse_lasso_VAR_est(SFM$F_hat[, 1:n_f], p = lags_f)
  }else{ 
    stop("invalid VAR_method")
  }
  
  VAR_idio <- diagVAR_est(v, lags_v)
  if (!is.null(max_EV)) {
    roots_phi <- check_roots(VAR_factors$coef)
    if (max(abs(roots_phi)) > max_EV) {
      Phi <- max_EV * VAR_factors$coef / max(abs(roots_phi))
    } else {
      Phi <- VAR_factors$coef
    }
    roots_delta <- check_roots(VAR_idio$coef)
    if (max(abs(roots_delta)) > max_EV) {
      Delta <- max_EV * VAR_idio$coef / max(abs(roots_delta))
    } else {
      Delta <- VAR_idio$coef
    }
  }
  est_factors <- list(Phi = Phi, Lambda = SFM$Lambda_hat[, 1:n_f],
                      H = focused_H(stats::cov(VAR_factors$resid), SFM$Lambda_hat[, 1:n_f],
                                    policy_var = ncol(X)), F_hat=SFM$F_hat[, 1:n_f])
  est_idio <- list(Delta = Delta, Xi = t(chol(VAR_idio$var)))
  return(list(factors = est_factors, idio = est_idio))
}

ols <- function(x, y) {
  xx_inv <- chol2inv(chol(crossprod(x)))
  b <- xx_inv %*% crossprod(x, y)
  e <- y - x %*% b
  return(list(coef = b, resid = e))
}

pc <- function(z, undo_scale = TRUE) {      #Principal Component Estimation
  z <- as.matrix(z)                         #make sure z/data is a matrix
  x <- scale(z)                             #if undo_scale = TRUE do not scale. If false, keep scale
  if (undo_scale) {
    xs <- z
  } else {
    xs <- x                                #if you do not want to undo the scale then use x, which is the scaled data.
  }
  n <- nrow(x)                             #the rows of the data matrix is the number of observation/time periods
  p <- ncol(x)                             #number of colums is the amount of variables
  evv <- eigen(x %*% t(x), symmetric = TRUE)
  F.hat <- sqrt(n) * evv$vectors
  l.hat <- crossprod(xs, F.hat) / n                                         
  scree <- evv$values
  V <- rep(NA,n)                                              #residuals of Xt = Lambda*f +V in the next lines we limit this to be a diagonal matrix (idiosyncratic)
  for (r in 1:min(n,p)) {                                      #included min(n,p) was 1:p. Mcracken and Ng use monthly dataset so have more observations n>p, whereas i use quarterly p>n
    V[r] <- sum(diag(crossprod(xs - F.hat[, 1:r, drop = FALSE] %*% t(l.hat[, 1:r, drop = FALSE])))) / (n * p);  #diagonal matrix of idiosyncratic shocks summed into RSS as seen in stock and watson 2016. also Stock and Watson Factor Model Loss Function
  }
  return(list(F_hat = F.hat, Lambda_hat = l.hat, eigenv = scree, RSS = V))
}

#' @export
pc_lean <- function(z, undo_scale = TRUE) {
  z <- as.matrix(z)
  x <- scale(z)
  if (undo_scale) {
    xs <- z
  } else {
    xs <- x
  }
  n <- nrow(x)
  p <- ncol(x)
  evv <- eigen(x %*% t(x), symmetric = TRUE)
  F.hat <- sqrt(n) * evv$vectors
  l.hat <- crossprod(xs, F.hat) / n
  return(list(F_hat = F.hat, Lambda_hat = l.hat))
}


sWF <- function(z, n_f=6, undo_scale=TRUE) {
  var_names<-colnames(z)
  z <- as.matrix(z)
  x <- scale(z)
  n <- nrow(x)
  p <- ncol(x)
  F.hat<-matrix(0,n,n)
  l.hat<-matrix(0,p,p);
  control_<-list(penA=FALSE, penD=FALSE, lam.min.factor=1e-6)
  s<-rrpack::sofar(Y=x, nrank=n_f, X=diag(n), control=control_, ic.type="GCV",)
  if(length(s$D)>1){
    Lambda<-s$V%*%diag(s$D)
  }else{
    Lambda<-s$V*s$D
  }
  l.hat[,1:ncol(Lambda)]<-Lambda
  F.hat[,1:ncol(s$U)]<-s$U
  if(undo_scale){
    l.hat<-diag(attributes(x)$'scaled:scale')%*%l.hat
  }
  rownames(l.hat)<-var_names
  return(list(F_hat = F.hat, Lambda_hat = l.hat))
}

g.BN <- function(N, T) {
  C.NT.1 <- (N * T)/(N + T)
  C.NT.2 <- min(N, T)
  g <- c(log(C.NT.1) / C.NT.1, log(C.NT.2) / C.NT.1, log(C.NT.2) / C.NT.2)
}

IC_BaiNg <- function(X, k_max = 8 * floor((min(dim(X)) / 100)^(1/4))) {
  PC.X <- pc(X)
  N <- ncol(X)
  T <- nrow(X)
  V.0 <- sum(diag(crossprod(X))) / (T * N)
  IC.m <- matrix(c(V.0, PC.X[[4]][1:k_max]), nrow = k_max + 1, ncol = 3) + (0:k_max) %o% g.BN(N, T)
  r.hat <- apply(IC.m, 2, which.min) - 1
  return(r.hat)
}

VAR_est <- function(X, p = 1) {
  X_lags <- create_lags(X, p, include.original = FALSE)
  return(ols(X_lags, X[-(1:p), ]))
}

lasso_VAR_est <- function(X, p=1){
  # Assuming your data matrix is in the object Y, please use this code then to estimate the VAR sparsely with p=1 lag:
  fit = bigtime::sparseVAR(Y = X, p = p, VARpen = 'L1', selection = 'bic', check_std = F)
  N<-ncol(X)
  Phi<-matrix(0, N*p, N)
  for(i in 1:p){
    Phi[((i-1)*N+1):(i*N),1:N]<-fit$Phihat[1:N,((i-1)*N+1):(i*N)]#Stephan stacks matrices on top of each other, Ines next to each other
  }
  X_lags <- create_lags(X, p, include.original = FALSE)
  e <- X[-(1:p), ] - X_lags %*% Phi
  return(list(coef=Phi, resid=e))
}

sparse_lasso_VAR_est <- function(X, p=1){
  # Assuming your data matrix is in the object Y, please use this code then to estimate the VAR sparsely with p=1 lag:
  fit = bigtime::sparseVAR(Y = X, p = p, VARpen = 'L1', selection = 'bic', verbose = F, check_std = F, VARgran=c(10,10))
  N<-ncol(X)
  Phi<-matrix(0, N*p, N)
  for(i in 1:p){
    Phi[((i-1)*N+1):(i*N),1:N]<-fit$Phihat[1:N,((i-1)*N+1):(i*N)]#Stephan stacks matrices on top of each other, Ines next to each other
  }
  X_lags <- create_lags(X, p, include.original = FALSE)
  e <- X[-(1:p), ] - X_lags %*% Phi
  return(list(coef=Phi, resid=e))
}

diagVAR_est <- function(X, p = 1) {
  N <- NCOL(X)
  coefs <- matrix(0, nrow = N * p, ncol = N)
  res_var <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    X_lags <- create_lags(X[, i], p, include.original = FALSE)
    AR <- ols(X_lags, X[-(1:p), i])
    coefs[1:p * N + i - N, i] <- AR$coef
    res_var[i, i] <- stats::var(AR$resid)
  }
  return(list(coef = coefs, var = res_var))
}

create_lags <- function(y, p = 1, include.original = TRUE, trim = TRUE) {     #create lags enters VAR_est and so y are our factors. This is our data
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = (p + include.original) * k)
  if (is.null(colnames(x))) {                                                #colnames are null on the factor matrix
    c.names <- rep("", k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, p + include.original)
  for (i in (1 - include.original):p) {
    cols <- k * (i - 1 + include.original) + 1:k
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, " l", i, sep = "")
  }
  return(lx[(1 + trim * p):n, , drop = FALSE])
}
################# Generate DFM / Artifical Data ########################
#' @export
generate_DFM <- function(n, factors, idio, init = 50, max_EV = 0.95,rho = c(0,0.2,0.5),sigma.omega = c(1.6,2.1,2.6)){
  n_f <- NCOL(factors$Lambda)
  p_f <- NROW(factors$Phi) / n_f
  N <- NCOL(idio$Delta)
  q_v <- NROW(idio$Delta) / N
  n1 <- n + init
 
   #check roots for stability
  roots_phi <- check_roots(factors$Phi)
  if (max(abs(roots_phi)) > max_EV) {
    Phi <- max_EV * factors$Phi / max(abs(roots_phi))
  } else {
    Phi <- factors$Phi
  }
  roots_delta <- check_roots(idio$Delta)
  if (max(abs(roots_delta)) > max_EV) {
    Delta <- max_EV * idio$Delta / max(abs(roots_delta))
  } else {
    Delta <- idio$Delta
  }
  
  #factor shock generation
  eps <- matrix(stats::rnorm(n1 * n_f), nrow = n1, ncol = n_f)
  eta <- eps %*% t(factors$H)
  f <- matrix(0, nrow = n1, ncol = n_f)
  for (t in (p_f + 1):n1) {
    f[t, ] <- c(t(f[t-1:p_f, ])) %*% Phi + eta[t, ]
  }
  
  #idiosyncratic variable generation
  e <- matrix(stats::rnorm(n1 * N), nrow = n1, ncol = N)
  u <- e %*% t(idio$Xi)
  v <- matrix(0, nrow = n1, ncol = N) 
  for (t in (q_v + 1):n1) { 
    v[t, ] <- c(t(v[t-1:q_v, ])) %*% Delta + u[t, ]
  }
  X <- f %*% t(factors$Lambda) + v
  
  #instrumental variable generation
  instrument = matrix(0,nrow=n1,ncol=1)
  R = sample(rho,1)
  Sig = sample(sigma.omega,1)
  omega = matrix(stats::rnorm(n1,mean=0,sd=Sig),nrow=n1,ncol=1)
  for (t in 2:n1){
    instrument[t,] = R*instrument[t-1,]+ omega[t,]
  }
  return(list(X = X[init + 1:n, ], obs_shock = eps[init + 1:n, ], IV = instrument[init+1:n,]))
}
################# Generate Impusle Responses via State Space modelling #######################
#' @export
impulse_response_ABCD <- function(factors, idio, S, h_max = 15, policy_var = length(S), outcome_var = 1) {       #policy_var = length(s) where S was our length of the dataset. THis assumes that our policy variable is ordered last in the dataset
  N <- NROW(factors$Lambda)
  n_w <- length(S)
  n_f <- NCOL(factors$Lambda) # number of factors
  p_f <- NROW(factors$Phi) / n_f  #lag order of our factor VAR
  
  # Reduce to selection
  q_v <- NROW(idio$Delta) / N #number of lags of the idiosyncratic components. Rows of Delta (variables) / original variables. If we have 2 lags then rows of delta 2x the number of original variables so q_v=2
  Delta_w <- idio$Delta[c(outer(S, (1:q_v - 1) * N, "+")), S] # selects the starting components of the lag blocks 
  Lambda_w <- factors$Lambda[S, ]
  Xi_w <- idio$Xi[S, ]
  
  # Find A, B, C, D
  if (p_f <= q_v + 1) {   # q_v is the number of lags in the idiosyncratic AR process (we have 1 lag). If the number of lags in our factor VAR is less than or equal to the number of lags in our idiosyncratic VAR +1
    n_s <- (1 + q_v) * n_f #n_s equal to 16 if we have one lag of idiosyncratic and one lag of factor VAR
    A <- companion(factors$Phi, q = q_v + 1)  #Produces [Phi 0 ,next row, 0 0]
    C <- cbind(diag(n_w), -t(Delta_w)) %*% (diag(1 + q_v) %x% Lambda_w) #Produces [Lambda -Delta'%*%Lambda]
  } else {
    n_s <- p_f * n_f 
    A <- companion(factors$Phi, q = p_f)
    Lambda_blocks <- matrix(0, nrow = (1 + q_v) * n_w, ncol = n_s)
    Lambda_blocks[1:((1 + q_v) * n_w), 1:((1 + q_v) * n_f)] <- diag(1 + q_v) %x% Lambda_w
    C <- cbind(diag(n_w), -t(Delta_w)) %*% Lambda_blocks
  }
  B <- rbind(cbind(factors$H, matrix(0, nrow = n_f, ncol = N)), matrix(0, nrow = n_s - n_f, ncol = N + n_f)) #companion matrix B uses contemporaneous effects stored in H
  D <- cbind(matrix(0, nrow = n_w, ncol = n_f), Xi_w)  # Produces [0 Xi] matrix
  
  out <- kalman_recursions(A, B, C, D) # Computes the steady state Kalmain Gain,Steady State error covariance matrix, innovation covariance Delta
  Sigma_uw <- out$Delta
  
  # g(L) = (I - (A - K C) L)^-1
  # A(L) = (I - C g(L) K L) (I - Delta_w(L)L) = A_w1(L) A_w2(L))]L
  
  g_L <- invert_VAR(t(A - out$K %*% C), q_max = h_max - 1)
  A_w1L <- sapply(1:h_max, function(j){C %*% g_L[, , j] %*% out$K}, simplify = "array")
  A_w1 <- array(c(diag(n_w), - A_w1L), dim = c(n_w, n_w, h_max + 1))
  A_w2 <- array(0, dim = c(n_w, n_w, h_max + 1))
  A_w2[, , 1] <- diag(n_w)
  A_w2[, , 1 + 1:q_v] <- array(-t(Delta_w), dim = c(n_w, n_w, q_v)) 
  A_w <- multiply_polynomials(A_w1, A_w2)[, , 1 + 0:h_max]
  
  C_w <- invert_VAR(A_w, q_max = h_max)
  B_w <- t(chol(Sigma_uw))
  
  CB_w <- sapply(0:h_max + 1, function(h){C_w[, , h] %*% B_w}, simplify = "array")
  
  IRF <- CB_w[outcome_var, policy_var, ] / CB_w[policy_var, policy_var, 1]
  IRF_cumsum = cumsum(IRF)
  list(IRF,IRF_cumsum)
  return(rbind(IRF,IRF_cumsum))
}

multiply_polynomials <- function(A, B) {
  p <- dim(A)[3]
  k <- dim(A)[1]
  C <- array(0, dim = c(k, k, 2*p - 1))
  for (i in 1:p) {
    for (j in 1:p) {
      C[, , i + j - 1] <- C[, , i + j - 1] + A[, , i] %*% B[, , j]
    }
  }
  return(C)
}

invert_VAR <- function(A, q_max = 20) {
  d <- length(dim(A))
  if (d == 2) { #2d array
    B <- invert_VAR_comp(A, q_max = q_max)
  } else if (d == 3) { # 3d array
    A_mat <- matrix(A[, , -1], nrow = dim(A)[1], ncol = dim(A)[2] * (dim(A)[3] - 1))
    A_mat_comp <- -t(solve(A[, , 1]) %*% A_mat) #solve produces inverse of the contemporaneous matrix A. A_mat is the matrices of future/dynamics
    B0 <- invert_VAR_comp(A_mat_comp, q_max = q_max)
    B <- sapply(0:q_max + 1, function(j){A[, , 1] %*% B0[, , j]}, simplify = "array")
  }
  return(B)
}

invert_VAR_comp <- function(A, q_max = 20) {
  k <- ncol(A)    
  p <- nrow(A) / k
  A_comp <- companion(A)
  B <- array(dim = c(k, k, 1 + q_max))
  B[, , 1] <- diag(k)
  A_powj <- diag(p * k)
  for (j in 1:q_max) {
    A_powj <- A_powj %*% A_comp
    B[, , 1 + j] <- A_powj[1:k, 1:k]
  }
  return(B)
}

companion <- function(A, q = nrow(A) / ncol(A)) {      ### A is our data input of coefficients to add into our companion form. Phi,Lambda,Xi,Delta. Default q is how many lags of the input we have. If we have Phi, q is p_v +1 which equals 2 in our case
  k <- ncol(A)   #number of variables (factors or observed variables depending on the input)
  p <- nrow(A) / k #number of lags of the input
  if (q > p) {  # If the desired number of lags q is greater than what’s in A, we need to pad with zeros to reach q.
    A_comp <- rbind(cbind(t(A), matrix(0, nrow = k , ncol = k * (q - p))), ##the matrix of zeroes needs to have k rows, not k*(q-p).   Take note of the transpose of A.    using Phi, we are returned a 16x16 matrix in companion form [Phi (8x8) 0(8x), next row, I(8x8) 0(8x8)]
                    diag(q * k)[1:((q - 1) * k), ]) # the identity part of the matrix advances lag states.  
  } else {
    if (p > 1) {
      A_comp <- rbind(t(A), diag(p * k)[1:((p - 1) * k), ])
    } else {
      A_comp <- t(A) #case for matrix already in block form
    }
  }
  return(A_comp)
}

check_roots <- function(A) {
  d <- length(dim(A))
  if (d == 2) {
    A_mat_comp <- A
  } else if (d == 3) {
    A_mat <- matrix(A[, , -1], nrow = dim(A)[1], ncol = dim(A)[2] * (dim(A)[3] - 1))
    A_mat_comp <- -t(solve(A[, , 1]) %*% A_mat)
  }
  A_comp <- companion(A_mat_comp)
  roots <- abs(eigen(A_comp)$values)
  return(roots)
}

kalman_recursions <- function(A, B, C, D, Omega0 = diag(nrow(A)), relax = 0, 
                              tol = 1e-8, max_iter = 1000) {
  Q <- B %*% t(B)
  R <- D %*% t(D)
  S <- B %*% t(D)
  
  for (i in 1:max_iter) {
    Delta <- C %*% Omega0 %*% t(C) + R
    Theta <- A %*% Omega0 %*% t(C) + S
    Delta_inv <- chol2inv(chol(Delta))
    K <- Theta %*% Delta_inv
    Omega <- A %*% Omega0 %*% t(A) + Q - K %*% t(Theta)
    if (max(abs(Omega - Omega0)) < tol) {
      break
    } else {
      Omega0 <- relax * Omega0 + (1 - relax) * Omega
    }
  }
  return(list(Omega = Omega, K = K, Delta = Delta))
}

focused_H <- function(S, Lambda, policy_var, shock_var = 1) {
  h1 <- S %*% Lambda[policy_var, ] / c(sqrt(t(Lambda[policy_var, ]) %*% S %*% Lambda[policy_var, ]))
  G <- t(chol(S))
  q1 <- solve(G, h1)
  x <- qr(q1)
  Q <- qr.Q(x, complete = TRUE)
  Q <- Q * sign(Q[1,1]) * sign(q1[1]) # fix sign
  H <- G %*% Q
  if (shock_var > 1) {
    if (shock_var < NCOL(S)) {
      H <- cbind(H[, 2:shock_var], H[, 1], H[, -(1:shock_var)])
    } else if (shock_var == NCOL(S)) {
      H <- cbind(H[, 2:shock_var], H[, 1])
    }
  }
  return(H)
}

generate_rnorm_toeplitz_factors <- function(N, n_f, p_f, phi_f, rho_f, threshold=0, diagonal=FALSE) {
  Lambda <- matrix(stats::rnorm(N * n_f), nrow = N)
  if(threshold!=0){
    Lambda<-hard_threshold(Lambda, threshold)
  }
  Phi <- matrix(0, nrow = p_f * n_f, ncol = n_f)
  if(!diagonal){
    for (j in 1:p_f) {
      Phi[n_f * (j - 1) + 1:n_f, 1:n_f] <- stats::toeplitz(phi_f[j]^(1:n_f))
    }
  }else{
    for (j in 1:p_f) {
      Phi[n_f * (j - 1) + 1:n_f, 1:n_f] <- phi_f[j]*diag(n_f) 
    }
  }
  
  ################# I took these out of the loop above
  if (check_roots(Phi)[1] > 0.999) {
    stop("Non-invertible lag polynomial")
  }
  Sigma_f <- matrix(rho_f, nrow = n_f, ncol = n_f) + diag(x = 1 - rho_f, nrow = n_f)
  H <- focused_H(Sigma_f, Lambda, N)
  return(list(Lambda = Lambda, H = H, Phi = Phi))
}

generate_multiAR <- function(N, q_v, delta_v, rho_v) { #####changed p_v to q_v just to avoid confusion with the notation
  Delta <- matrix(0, nrow = q_v * N, ncol = N)
  for (j in 1:q_v) {
    Delta[(j - 1) * N + 1:N, 1:N] <- diag(stats::runif(N, min = delta_v[1, j], max = delta_v[2, j]))
  }
  if (check_roots(Delta)[1] > 0.999) {
    stop("Non-invertible lag polynomial")
  }
  Sigma_v <- matrix(rho_v, nrow = N, ncol = N) + diag(1 - rho_v, nrow = N)
  Xi <- t(chol(Sigma_v))
  return(list(Xi = Xi, Delta = Delta))
}

hard_threshold<-function(mat, threshold){
  apply(X=mat, MARGIN=c(1,2), FUN=function(x){if(abs(x)<=threshold){0}else{x}})
}

#' @export
one_replication_lean<-function(dummy_list=list(),Ts=c(200, 400, 600), DFM, IRF){
  h_max<-length(IRF)-1
  intervals<-array(0,dim=c(h_max+1,3,length(Ts), 3, 3), dimnames=list(
    horizon=0:h_max, 
    information=c("lower","bhat","upper"), 
    T_=c(paste0("T_",Ts)), 
    DGP=c("HDLP_04", "LP", "FALP"), 
    LRV=c("intervals","intervals_EWC","intervals_NWfb")))
  
  covered<-array(FALSE,dim=c(h_max+1,length(Ts), 3, 3), dimnames=list(
    horizon=0:h_max, 
    T_=c(paste0("T_",Ts)), 
    DGP=c("HDLP_04", "LP", "FALP"), 
    LRV=c("intervals","intervals_EWC","intervals_NWfb")))
  
  width<-array(0,dim=c(h_max+1,length(Ts), 3, 3), dimnames=list(
    horizon=0:h_max, 
    T_=c(paste0("T_",Ts)), 
    DGP=c("HDLP_04", "LP", "FALP"), 
    LRV=c("intervals","intervals_EWC","intervals_NWfb")))
  
  for(ts in 1:length(Ts)){
    T_<-Ts[ts]
    S<-1:122
    n_f=6
    data <- generate_DFM(n=T_, DFM$factors, DFM$idio, init = 50, max_EV = 0.98)
    Z <- data$X[, S] #Z is the selection of data to be used
    x <- Z[, ncol(Z)] #set up the policy variable
    y <- Z[, 1] # outcome variable of our interest
    slow <- Z[, 3:ncol(Z) - 1] #every other variable is predetermined since we use only observed shocks and IV. All variables except HICP and IRT3M (non lagged)
    f<-pc(z=slow)  #not needed we are not doing FAVAR/FALP
    slow_factors<-f$F_hat[,1:n_f]
    cpi<-Z[,"HICPOV_EA"]
    gdp = Z[,"GDP_EA"]
    e <- data$obs_shock[, 1] #observed shock to the first structural shock of the Factor VAR
    IV = data$IV
    
    HDLP_04 <- desla::HDLP(x = x, y = y, r = slow, PI_constant=0.4, y_predetermined = TRUE, hmax = h_max, lags = 3, progress_bar = FALSE, parallel=FALSE)
    
    #HDLP_08 <- desla::HDLP(x = x, y = y, r = slow, PI_constant=0.8, y_predetermined = TRUE, hmax = h_max, lags = 3, progress_bar = FALSE, parallel=FALSE)
    
    LP <- desla::HDLP(x = x, y = y, r = gdp, OLS=TRUE, y_predetermined = TRUE, hmax = h_max, lags = 3, progress_bar = FALSE, parallel=FALSE)
    
    FALP <- desla::HDLP(x = x, y = y, r = slow_factors, OLS=TRUE, y_predetermined = TRUE, hmax = h_max, lags = 3, progress_bar = FALSE, parallel=FALSE)
    
    #observed <- desla::HDLP(x = e, y = y, r = slow, q = x, PI_constant=0.4, y_predetermined = TRUE, hmax = h_max, lags = 3, progress_bar = FALSE, parallel=FALSE)
    
    for(dgp_name in c("HDLP_04","LP","FALP")){
      for(lrv_name in c("intervals", "intervals_EWC","intervals_NWfb")){
        intervals[,,ts,dgp_name,lrv_name]<-get(lrv_name, get(dgp_name))
        for(h in 1:(h_max+1)){
          lower<-intervals[h,"lower",ts,dgp_name,lrv_name]
          upper<-intervals[h,"upper",ts,dgp_name,lrv_name]
          if(lower<= IRF[h] && IRF[h]<=upper){
            covered[h,ts,dgp_name,lrv_name]<-TRUE
          }else{
            covered[h,ts,dgp_name,lrv_name]<-FALSE
          }
          w<-upper-lower
          width[h,ts,dgp_name,lrv_name]<-w
        }
        
      }
    }
  }
  return(list(covered=covered, width=width))
}
#' @export
# small helper used below
`%||%` <- function(a, b) if (is.null(a)) b else a # If LHS is null return the RHS : NULL %||% 100 = 100
#' @export
grab_ci <- function(obj) {
  # prefer the per-horizon pointwise CI + theta
  if (!is.null(obj$ci_pointwise)) {
    lo <- obj$ci_pointwise[, "lower", drop = TRUE]
    up <- obj$ci_pointwise[, "upper", drop = TRUE]
    bh <- obj$theta
    return(cbind(lower = lo, bhat = bh, upper = up))
  }
  
  # otherwise use obj$intervals (can be matrix or 3-D array)
  if (!is.null(obj$intervals)) {
    A <- obj$intervals
    
    # case 1: already a (H+1) x 3 matrix
    if (is.matrix(A)) {
      cols <- colnames(A)
      if (!is.null(cols) && all(c("lower","upper") %in% cols)) {
        if (!("bhat" %in% cols) && !is.null(obj$theta)) {
          # make bhat from theta if missing
          A <- cbind(lower = A[, "lower", drop=TRUE],
                     bhat  = obj$theta,
                     upper = A[, "upper", drop=TRUE])
        } else {
          A <- A[, c("lower","bhat","upper"), drop = FALSE]
        }
      } else {
        # assume columns are in order (lower, bhat, upper)
        if (ncol(A) != 3) stop("Unknown DML intervals matrix layout")
        colnames(A) <- c("lower","bhat","upper")
      }
      return(A)
    }
    
    # case 2: a 3-D array like [h, info, something]
    if (is.array(A) && length(dim(A)) >= 2) {
      info_names <- dimnames(A)[[2]]
      li <- if (!is.null(info_names)) match("lower", info_names) else 1
      bi <- if (!is.null(info_names)) match("bhat",  info_names) else 2
      ui <- if (!is.null(info_names)) match("upper", info_names) else 3
      if (any(is.na(c(li,bi,ui)))) stop("Could not find lower/bhat/upper in intervals array")
      # take the first extra-slice if present
      slice <- if (length(dim(A)) >= 3) 1 else NULL
      lo <- A[, li, slice, drop = TRUE]
      bh <- A[, bi, slice, drop = TRUE]
      up <- A[, ui, slice, drop = TRUE]
      return(cbind(lower = lo, bhat = bh, upper = up))
    }
  }
  
  stop("No intervals in DML fit")
}

#' @export
extract_selected <- function(fit, spec_label) {
  out <- list()
  for (h in seq_along(fit$betas_g)) {
    mat <- fit$betas_g[[h]]
    if (is.null(mat)) next
    active <- rownames(mat)[rowSums(abs(mat)) > 0]
    out[[h-1]] <- data.frame(
      horizon  = h-1,
      nuisance = "g",
      spec     = spec_label,
      var      = active
    )
  }
  for (h in 1) { # m(X) is horizon-invariant
    mat <- fit$betas_m
    active <- rownames(mat)[rowSums(abs(mat)) > 1e-8]
    out[["m"]] <- data.frame(
      horizon  = NA,
      nuisance = "m",
      spec     = spec_label,
      var      = active
    )
  }
  do.call(rbind, out)
}

#' @export
hdlp_regressors <- function(r_names, lags, hmax, y_name="Y", x_name="X") {
  r_terms <- unlist(lapply(r_names, function(v) paste0(v, "_lag", 0:(lags+hmax))))
  y_terms <- paste0(y_name, "_lag", 1:(lags+hmax))
  x_terms <- paste0(x_name, "_lag", 0:(lags+hmax))
  c(r_terms, y_terms, x_terms)
}


#' @export
one_replication_lean_manual_seed_safe <- function(
    seed,
    Ts  = c(200,400,600),
    DFM, IRF,
    dml_grid = NULL
){
  tryCatch({
    one_replication_lean_manual_seed(seed, Ts = Ts, DFM = DFM, IRF = IRF, dml_grid = dml_grid)
  }, error = function(e) {
    message("⚠️ Replication dropped (seed=", seed, "): ", conditionMessage(e))
    return(NULL)   # discard weak instrument / failed replication
  })
}














#' @export
one_replication_lean_manual_seed <- function(
    seed,
    Ts  = c(200,400,600),
    DFM, IRF,
    dml_grid = NULL
){
  
  set.seed(seed)
  h_max <- length(IRF) - 1
  scale_IRF <- sqrt(mean(IRF^2))
  
  # DGP names
  base_DGP <- character(0)
  if (is.null(dml_grid) || nrow(dml_grid) == 0) { #  if dml grid is not provided or no rows are assigned, build one specification with default settings
    dml_grid <- data.frame(
      label = "DML_LP",
      nuisance_m="asw", nuisance_g="asw",
      stringsAsFactors = FALSE
    )
  }
  if (!("label" %in% names(dml_grid))) { #assign label/column for model/specification name if none
    dml_grid$label <- paste0("DML_", seq_len(nrow(dml_grid)))
  }
  dml_labels <- dml_grid$label #names of all our DML specifications
  all_DGP    <- c(base_DGP, dml_labels) # names of all our DGPs
  
  # containers
  LRV_names <- c("intervals")
  
  intervals <- array(NA_real_,
                     dim = c(h_max+1, 3, length(Ts), length(all_DGP), length(LRV_names)),
                     dimnames = list(
                       horizon     = 0:h_max,
                       information = c("lower","bhat","upper"),
                       T_          = paste0("T_", Ts),
                       DGP         = all_DGP,
                       LRV         = LRV_names
                     )
  )
  
  covered <- array(FALSE,
                   dim = c(h_max+1, length(Ts), length(all_DGP), length(LRV_names)),
                   dimnames = list(
                     horizon = 0:h_max,
                     T_      = paste0("T_", Ts),
                     DGP     = all_DGP,
                     LRV     = LRV_names
                   )
  )
  
  width <- array(NA_real_,
                 dim = c(h_max+1, length(Ts), length(all_DGP), length(LRV_names)),
                 dimnames = dimnames(covered)
  )
  
  bias_abs <- array(NA_real_,
                    dim = c(h_max+1, length(Ts), length(all_DGP), length(LRV_names)),
                    dimnames = dimnames(covered)
  )
  
  t.start = Sys.time()
  for (ts in seq_along(Ts)) {
    Tn <- Ts[ts]
    
    # -------- data draw --------
    S   <- 1:122
    sim <- generate_DFM(n = Tn, DFM$factors, DFM$idio, init = 50, max_EV = 0.98)
    Z   <- sim$X[, S, drop = FALSE]
    x <- Z[, ncol(Z)]                # policy
    y <- Z[, 1]                      # outcome
    cpi <- Z[, "CPIAUCSL", drop = FALSE]
    slow_idx <- 3:ncol(Z)-1# everything in the middle
    slow <- Z[, slow_idx, drop = FALSE]
    
    # factors for FALP
    fhat <- pc_lean(z = slow)
    n_f = 6
    slow_factors <- fhat$F_hat[, 1:n_f, drop = FALSE]
    
    # data.table for DML
    dt <- data.table::as.data.table(Z)
    y_var     <- colnames(Z)[1]
    shock_var <- colnames(Z)[ncol(Z)]
    # -------- benchmarks --------
    #HDLP_04 <- desla::HDLP(
     # x = x, y = y, r = slow,
      #PI_constant = 0.4, y_predetermined = TRUE,
      #hmax = h_max, lags = 3, progress_bar = FALSE, parallel = FALSE
    #)
    
    #LP <- desla::HDLP(
      #x = x, y = y, r = cpi,
      #OLS = TRUE, y_predetermined = TRUE,
      #hmax = h_max, lags = 3, progress_bar = FALSE, parallel = FALSE
    #)
    #FALP <- desla::HDLP(
      #x = x, y = y, r = slow_factors,
      #OLS = TRUE, y_predetermined = TRUE,
      #hmax = h_max, lags = 3, progress_bar = FALSE, parallel = FALSE
    #)
    
    # put the three LRVs in for each benchmark model
    #for (dgp_name in base_DGP) {
      #obj <- get(dgp_name)
      #for (lrv_name in LRV_names) {
        #ci <- get(lrv_name, obj)                     # matrix (H+1) x 3
        ##intervals[,, ts, dgp_name, lrv_name] <- ci   # write all 3 cols
        #for (h in 1:(h_max+1)) {
          #lo <- intervals[h,"lower",ts,dgp_name,lrv_name]
          #hi <- intervals[h,"upper",ts,dgp_name,lrv_name]
          #bh <- intervals[h,"bhat", ts,dgp_name,lrv_name]
          #covered[h,ts,dgp_name,lrv_name] <- (lo <= IRF[h] && IRF[h] <= hi)
          #width[  h,ts,dgp_name,lrv_name] <- (hi - lo)
          #bias_abs[h,ts,dgp_name,lrv_name] <- abs(bh - IRF[h]) / scale_IRF
        #}
      #}
    #}
    
    # -------- DML grid  --------
    dml_fits <- list()
    for (i in seq_len(nrow(dml_grid))) {
      cfg <- as.list(dml_grid[i, , drop = FALSE])
      dml_name <- dml_grid$label[i]
    
      fit <- dml_lp(
        dt, y_var = y_var, shock_var = shock_var,
        H = h_max, lags = 3,
        theta_kernel = cfg$theta_kernel %||% "bartlett",
        theta_bw     = cfg$theta_bw     %||% "andrews_ar1",
        
        nuisance_m       = cfg$nuisance_m   %||% "asw",
        m_asw_c_pi       = cfg$m_asw_c_pi   %||% 0.6,
        m_asw_alpha      = cfg$m_asw_alpha  %||% 0.05,
        m_asw_B          = cfg$m_asw_B      %||% 1000,
        m_asw_aggL       = cfg$m_asw_aggL   %||% "p95",
        scale_outcome_m  = cfg$scale_outcome_m %||% FALSE,
        
        nuisance_g       = cfg$nuisance_g   %||% "asw",
        g_asw_c_pi       = cfg$g_asw_c_pi   %||% 0.6,
        g_asw_alpha      = cfg$g_asw_alpha  %||% 0.05,
        g_asw_B          = cfg$g_asw_B      %||% 1000,
        g_asw_aggL       = cfg$g_asw_aggL   %||% "p95",
        scale_outcome_g  = cfg$scale_outcome_g %||% FALSE
      )
      
      fit$meta$label      <- dml_name
      ci_mat <- grab_ci(fit)  # (H+1) x 3: lower/bhat/upper
      dml_fits[[dml_name]] <- fit
      # fill ALL LRVs with the same DML CI so process_sim won't break
      for (lrv_name in LRV_names) {
        intervals[,, ts, dml_name, lrv_name] <- ci_mat
        for (h in 1:(h_max+1)) {
          lo <- intervals[h,"lower",ts,dml_name,lrv_name]
          hi <- intervals[h,"upper",ts,dml_name,lrv_name]
          bh <- intervals[h,"bhat", ts,dml_name,lrv_name]
          covered[h,ts,dml_name,lrv_name] <- (lo <= IRF[h] && IRF[h] <= hi)
          width[  h,ts,dml_name,lrv_name] <- (hi - lo)
          bias_abs[h,ts,dml_name,lrv_name] <- abs(bh - IRF[h]) / scale_IRF
        }
      }
    }
  }
  t.end = Sys.time()
  total_time = t.end - t.start
  
  list(
    covered   = covered,
    width     = width,
    bias_abs  = bias_abs,
    intervals = intervals,
    DML       = dml_fits
  )
}



#' @export
process_sim<-function(sim){
  
  DGPs<-dimnames(sim[[1]]$covered)$DGP
  LRVs<-dimnames(sim[[1]]$covered)$LRV
  M<-length(sim)
  T_strings<-dimnames(sim[[1]]$covered)$T_
  h_max<-dim(sim[[1]]$covered)[1]-1
  
  coverages<-array(0,dim=c(h_max+1,length(T_strings), length(DGPs), length(LRVs)), dimnames=list(
    horizon=0:h_max, 
    T_=T_strings, 
    DGP=DGPs, 
    LRV=LRVs))
  
  median_widths<-array(0,dim=c(h_max+1,length(T_strings), length(DGPs), length(LRVs)), dimnames=list(
    horizon=0:h_max, 
    T_=T_strings, 
    DGP=DGPs, 
    LRV=LRVs))
  
  for(dgp_name in DGPs){
    for(lrv_name in LRVs){
      for(h in 1:(h_max+1)){
        for(ts in T_strings){
          temp_widths<-rep(0,M)
          for(m in 1:M){
            if(sim[[m]]$covered[h,ts,dgp_name,lrv_name]){
              coverages[h,ts,dgp_name,lrv_name]<-coverages[h,ts,dgp_name,lrv_name]+1/M
            }
            temp_widths[m]<-sim[[m]]$width[h,ts,dgp_name,lrv_name]
          }
          median_widths[h,ts,dgp_name,lrv_name]<-stats::median(temp_widths)
        }
      }
    }
  }
  return(list(coverages=coverages, median_widths=median_widths, M=M, h_max=h_max, Ts=as.numeric(gsub("T_", "", T_strings)), DGPs=DGPs, LRVs=LRVs))
}







## Plotting Everything





#' @export
plot_processed_sim<-function(processed_sim, models=NULL, min=0){
  horizon=model=lrv=value=NULL
  DGPs=processed_sim$DGPs 
  LRVs=processed_sim$LRVs
  Ts<-processed_sim$Ts
  T_strings<-paste0("T_", Ts)
  h_max<-processed_sim$h_max
  if(is.null(models)){
    models=DGPs
  }
  for(ts in T_strings){
    temp_cov<-data.frame("value"=NULL, "horizon"=NULL, "model"=NULL, "lrv"=NULL)
    temp_mw<-data.frame("value"=NULL, "horizon"=NULL, "model"=NULL, "lrv"=NULL)
    for(dgp_name in models){
      for(lrv_name in LRVs){
        d1<-data.frame("value"=processed_sim$coverages[,ts,dgp_name,lrv_name], "horizon"=0:h_max, "model"=dgp_name, "lrv"=lrv_name)
        temp_cov<-rbind(temp_cov,d1)
        d2<-data.frame("value"=processed_sim$median_widths[,ts,dgp_name,lrv_name], "horizon"=0:h_max, "model"=dgp_name, "lrv"=lrv_name)
        temp_mw<-rbind(temp_mw,d2)
      }
    }
    #temp_cov<-cbind(temp_cov,data.frame(id=1:nrow(temp_cov)))
    #temp_mw<-cbind(temp_mw,data.frame(id=1:nrow(temp_mw)))
    assign(paste0(ts,"_cov_df"), temp_cov)
    assign(paste0(ts,"_mw_df"), temp_mw)
  }
  plot_list<-vector("list", 2*length(Ts))
  for(dgp_name in models){
    for(lrv_name in c("intervals", "intervals_EWC","intervals_NWfb")){
      for(ts in T_strings){
        p_cov<-ggplot2::ggplot(data=get(paste0(ts,"_cov_df")), mapping=ggplot2::aes(x=horizon, y=value, colour=model, linetype=lrv))+
          ggplot2::theme_bw()+
          ggplot2::geom_hline(yintercept=0.95, color="black", linetype=1, linewidth=0.3)+
          ggplot2::geom_line(linewidth=0.3)+
          ggplot2::scale_linetype_manual(values=c(1,5,3),
                                         name="LRV", 
                                         breaks=c("intervals", "intervals_EWC","intervals_NWfb"), 
                                         labels = c("NW", "EWC", "NWfb"))+
          ggplot2::labs(title=ts)+
          ggplot2::ylab("Coverage")+
          ggplot2::xlab("Horizon")
        assign(paste0("p_cov_",ts),p_cov)
        plot_list[[which(T_strings==ts)]]<-get(paste0("p_cov_",ts))
        
        p_mw<-ggplot2::ggplot(data=get(paste0(ts,"_mw_df")), mapping=ggplot2::aes(x=horizon, y=value, colour=model, linetype=lrv))+
          ggplot2::theme_bw()+
          ggplot2::geom_line(linewidth=0.3)+
          ggplot2::scale_linetype_manual(values=c(1,5,3),
                                         name="LRV", 
                                         breaks=c("intervals", "intervals_EWC","intervals_NWfb"), 
                                         labels = c("NW", "EWC", "NWfb"))+
          ggplot2::labs(title=ts)+
          ggplot2::xlab("Horizon")+
          ggplot2::ylab("Median width")
        assign(paste0("p_mw_",ts),p_mw)
        plot_list[[length(Ts)+which(T_strings==ts)]]<-get(paste0("p_mw_",ts))
      }
    }
  }
  combined<-ggpubr::ggarrange(plotlist=plot_list, common.legend=TRUE, legend="right")
  
  if(setequal(models,c("HDLP_04", "FALP", "LP"))){
    #separate only models
    for(ts in T_strings){
      temp_cov<-data.frame(value=NULL, horizon=NULL, model=NULL,lrv=NULL)
      temp_mw<-data.frame(value=NULL, horizon=NULL, model=NULL,lrv=NULL)
      for(dgp_name in models){
        d1<-data.frame(value=processed_sim$coverages[,ts,dgp_name,"intervals"], horizon=0:h_max, model=dgp_name)
        temp_cov<-rbind(temp_cov,d1)
        d2<-data.frame(value=processed_sim$median_widths[,ts,dgp_name,"intervals"], horizon=0:h_max, model=dgp_name)
        temp_mw<-rbind(temp_mw,d2)
      }
      #temp_cov<-cbind(temp_cov,data.frame(id=1:nrow(temp_cov)))
      #temp_mw<-cbind(temp_mw,data.frame(id=1:nrow(temp_mw)))
      assign(paste0(ts,"_only_models_cov_df"), temp_cov)
      assign(paste0(ts,"_only_models_mw_df"), temp_mw)
    }
    only_models<-vector("list", length(Ts))
    for(ts in T_strings){
      p_cov<-ggplot2::ggplot(data=get(paste0(ts,"_only_models_cov_df")), mapping=ggplot2::aes(x=horizon, y=value, color=model))+
        ggplot2::theme_bw()+
        ggplot2::geom_hline(yintercept=0.95, color="black", linetype=1, linewidth=0.3)+
        ggplot2::geom_line()+
        #ggplot2::labs(title=ts)+
        ggplot2::labs(title=paste0("T=", Ts[which(T_strings==ts)]))+
        ggplot2::scale_color_manual(values=c("blue","darkorange","darkgreen"),
                                    name="Model", 
                                    breaks=c("HDLP_04", "FALP", "LP"), 
                                    labels = c("HDLP", "FA-LP", "LP"))+
        #ggplot2::ylab("Coverage")+
        ggplot2::xlab("Horizon")
      assign(paste0("p_cov_only_models_",ts),p_cov)
      only_models[[which(T_strings==ts)]]<-get(paste0("p_cov_only_models_",ts))
      
      p_mw<-ggplot2::ggplot(data=get(paste0(ts,"_only_models_mw_df")), mapping=ggplot2::aes(x=horizon, y=value, color=model))+
        ggplot2::theme_bw()+
        ggplot2::geom_line()+
        #ggplot2::labs(title=ts)+
        ggplot2::labs(title=paste0("T=", Ts[which(T_strings==ts)]))+
        ggplot2::scale_color_manual(values=c("blue","darkorange","darkgreen"),
                                    name="Model", 
                                    breaks=c("HDLP_04", "FALP", "LP"), 
                                    labels = c("HDLP", "FA-LP", "LP"))+
        #ggplot2::ylab("Median width")+
        ggplot2::xlab("Horizon")
      assign(paste0("p_mw_only_models_",ts),p_mw)
      only_models[[length(Ts)+which(T_strings==ts)]]<-get(paste0("p_mw_only_models_",ts))
    }
    only_models[[1]]<-only_models[[1]]+ggplot2::ylab("Coverage")
    only_models[[length(Ts)+1]]<-only_models[[length(Ts)+1]]+ggplot2::ylab("Median width")
    combined_only_models<-ggpubr::ggarrange(plotlist=only_models, common.legend=TRUE, legend="right")
    
    #only LRVs
    for(ts in T_strings){
      temp_cov<-data.frame(value=NULL, horizon=NULL, model=NULL,lrv=NULL)
      temp_mw<-data.frame(value=NULL, horizon=NULL, model=NULL,lrv=NULL)
      for(lrv_name in c("intervals", "intervals_EWC","intervals_NWfb")){
        d1<-data.frame(value=processed_sim$coverages[,ts,"HDLP_04",lrv_name], horizon=0:h_max, lrv=lrv_name)
        temp_cov<-rbind(temp_cov,d1)
        d2<-data.frame(value=processed_sim$median_widths[,ts,"HDLP_04",lrv_name], horizon=0:h_max, lrv=lrv_name)
        temp_mw<-rbind(temp_mw,d2)
      }
      #temp_cov<-cbind(temp_cov,data.frame(id=1:nrow(temp_cov)))
      #temp_mw<-cbind(temp_mw,data.frame(id=1:nrow(temp_mw)))
      assign(paste0(ts,"_only_lrvs_cov_df"), temp_cov)
      assign(paste0(ts,"_only_lrvs_mw_df"), temp_mw)
    }
    only_lrvs<-vector("list", length(Ts))
    for(ts in T_strings){
      p_cov<-ggplot2::ggplot(data=get(paste0(ts,"_only_lrvs_cov_df")), mapping=ggplot2::aes(x=horizon, y=value, linetype=lrv))+
        ggplot2::theme_bw()+
        ggplot2::geom_hline(yintercept=0.95, color="black", linetype=1, linewidth=0.3)+
        ggplot2::geom_line()+
        #ggplot2::labs(title=ts)+
        ggplot2::labs(title=paste0("T=", Ts[which(T_strings==ts)]))+
        ggplot2::scale_linetype_manual(values=c(1,5,3),
                                       name="LRV", 
                                       breaks=c("intervals", "intervals_EWC","intervals_NWfb"), 
                                       labels = c("NW", "EWC", "NW-fb"))+
        #ggplot2::ylab("Coverage")+
        ggplot2::xlab("Horizon")
      assign(paste0("p_cov_only_lrvs_",ts),p_cov)
      only_lrvs[[which(T_strings==ts)]]<-get(paste0("p_cov_only_lrvs_",ts))
      
      p_mw<-ggplot2::ggplot(data=get(paste0(ts,"_only_lrvs_mw_df")), mapping=ggplot2::aes(x=horizon, y=value, linetype=lrv))+
        ggplot2::theme_bw()+
        ggplot2::geom_line()+
        #ggplot2::labs(title=ts)+
        ggplot2::labs(title=paste0("T=", Ts[which(T_strings==ts)]))+
        ggplot2::scale_linetype_manual(values=c(1,5,3),
                                       name="LRV", 
                                       breaks=c("intervals", "intervals_EWC","intervals_NWfb"), 
                                       labels = c("NW", "EWC", "NW-fb"))+
        #ggplot2::ylab("Median width")+
        ggplot2::xlab("Horizon")
      assign(paste0("p_mw_only_lrvs_",ts),p_mw)
      only_lrvs[[length(Ts)+which(T_strings==ts)]]<-get(paste0("p_mw_only_lrvs_",ts))
    }
    only_lrvs[[1]]<-only_lrvs[[1]]+ggplot2::ylab("Coverage")
    only_lrvs[[length(Ts)+1]]<-only_lrvs[[length(Ts)+1]]+ggplot2::ylab("Median width")
    combined_only_lrvs<-ggpubr::ggarrange(plotlist=only_lrvs, common.legend=TRUE, legend="right")
  }else{
    only_models<-NULL
    combined_only_models<-NULL
    combined_only_lrvs<-NULL
  }
  
  return(list(plot_list=plot_list, combined=combined, only_models=only_models, combined_only_models=combined_only_models, combined_only_lrvs=combined_only_lrvs))
}
