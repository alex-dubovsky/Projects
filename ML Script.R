
################################################################################
                                  # UTILITIES
################################################################################
#' @import data.table
#' @importFrom data.table :=

# -------- Data matrix: add lags for all, H leads for y --------
#' @export
build_matrix <- function(dt, y_var, shock_var, H = 20, lags = 3) {
  dt <- data.table::as.data.table(dt)
  cols <- names(dt)
  for (cl in cols) {
    for (L in seq_len(lags)) {
      nm <- paste0(cl, "_lag", L)
      dt[, (nm) := data.table::shift(get(cl), n = L, type = "lag")]
    }
  }
  y0 <- dt[[y_var]]
  for (L in seq_len(H)) {
    nm <- paste0(y_var, "_lead", L)
    dt[, (nm) := data.table::shift(y0, n = L, type = "lead")]
  }
  stats::na.omit(dt[])
}

# -------- NLO folds with buffers --------
#' @export
create_nlo_folds <- function(Tn, K = 10, buffer_folds = 1) {
  fold_id <- floor(((seq_len(Tn) - 1) * K) / Tn) + 1L
  idx_by_fold <- split(seq_len(Tn), fold_id)
  folds <- vector("list", K)
  for (k in seq_len(K)) {
    test_idx <- idx_by_fold[[k]]
    neighbours <- setdiff(seq(max(1, k - buffer_folds), min(K, k + buffer_folds)), k)
    buf_idx <- unlist(idx_by_fold[neighbours], use.names = FALSE)
    train_idx <- setdiff(seq_len(Tn), c(test_idx, buf_idx))
    folds[[k]] <- list(train = train_idx, test = test_idx)
  }
  folds
}

# -------- HAC bits --------
#' @export
nw_lag_default <- function(Tn) floor(4 * (Tn / 100)^(2/9))
#' @export
kernel_w <- function(l, L, kernel = c("bartlett","cosine","uniform")) {
  kernel <- match.arg(kernel)
  if (L <= 0) return(1)
  if (kernel == "bartlett") 1 - l / (L + 1)
  else if (kernel == "cosine") base::cos((pi/2) * l / (2 * (L + 1)))
  else 1
}


# Full HAC Ω, finite-adjust: 1/(T-l) on off-diagonals
#' @export
hac_lrv_matrix_kernel <- function(Psi, L, kernel = c("bartlett","cosine","uniform"), finite_adjust = TRUE) {
  kernel <- match.arg(kernel)
  Psi <- as.matrix(Psi); Tn <- nrow(Psi)
  L <- max(0L, min(as.integer(L), Tn - 2L))
  Omega <- crossprod(Psi) / Tn
  if (L > 0) {
    for (l in 1:L) {
      w <- kernel_w(l, L, kernel)
      denom <- if (finite_adjust) (Tn - l) else Tn
      G <- t(Psi[(l+1):Tn, , drop = FALSE]) %*% Psi[1:(Tn-l), , drop = FALSE] / denom
      Omega <- Omega + w * (G + t(G))
    }
  }
  (Omega + t(Omega)) / 2
}

# Diagonal HAC var (fast), same finite-adjust
#' @export
hac_var_diag_truncated <- function(Z, L, kernel = c("bartlett","cosine","uniform"), finite_adjust = TRUE) {
  kernel <- match.arg(kernel)
  Z <- as.matrix(Z); Tn <- nrow(Z)
  v <- colMeans(Z^2)
  if (L <= 0) return(pmax(v, 1e-12))
  for (l in 1:L) {
    w <- kernel_w(l, L, kernel)
    if (finite_adjust) {
      v <- v + 2 * w * colSums(Z[(l+1):Tn, , drop = FALSE] * Z[1:(Tn-l), , drop = FALSE]) / (Tn - l)
    } else {
      v <- v + 2 * w * colMeans(Z[(l+1):Tn, , drop = FALSE] * Z[1:(Tn-l), , drop = FALSE])
    }
  }
  pmax(v, 1e-12)
}

# Andrews (1991) AR(1) bandwidth eq.(6.4), Bartlett class
#' @export
L_eq64_bartlett_vec <- function(s) {
  s <- as.numeric(s) - mean(s)
  n <- length(s)
  if (n < 3) return(0L)
  num <- sum(s[-1] * s[-n]); den <- sum(s[-n]^2)
  rho <- num/den; rho <- pmin(pmax(rho, -0.99), 0.99)
  alpha1 <- 4 * rho^2 / (1 - rho^2)^2
  L <- floor(1.1447 * (alpha1 * n)^(1/3))
  as.integer(max(0L, min(L, n - 2L)))
}
#' @export
choose_common_L_eq64 <- function(S, agg = c("p95","max","median")) {
  agg <- match.arg(agg)
  S <- as.matrix(S)
  if (!ncol(S)) return(0L)
  Ls <- apply(S, 2, L_eq64_bartlett_vec)
  if (agg == "p95") as.integer(stats::quantile(Ls, 0.95, na.rm = TRUE))
  else if (agg == "max") as.integer(max(Ls, na.rm = TRUE))
  else as.integer(stats::median(Ls, na.rm = TRUE))
}

# Safe scaling helpers (store, then apply)
#' @export
safe_scale_fit   <- function(X) list(center = colMeans(X), scale = {s <- apply(X,2,sd); s[!is.finite(s)|s==0] <- 1; s})
#' @export
safe_scale_apply <- function(X, center, scale) sweep(sweep(X,2,center,"-"),2,scale,"/")


################################################################################
                            #NUISANCE LEARNERS
################################################################################


# ---------------- Adamek–Smeekes–Wilms: pick one Q_T for this T ----------------
#' @export
compute_QT_bandwidth <- function(X, y, mu, sdv, aggL = "p95") {
  stopifnot(is.finite(sdv), sdv > 0)
  Xs <- scale(X)
  y_std <- (as.numeric(y) - mu) / sdv
  u0 <- y_std - mean(y_std)
  S0 <- Xs * u0
  as.integer(choose_common_L_eq64(S0, agg = aggL))
}

# Fit ASW with fixed Q_T (sup-norm bootstrap for λ)
#' @export
fit_asw_with_QT <- function(X_tr, y_tr, Q_T,
                            c_pi = 0.8, alpha = 0.05, B = 1000,
                            full_Omega = FALSE,
                            y_center, y_scale) {
  stopifnot(length(Q_T) == 1L, Q_T >= 0, is.finite(y_scale), y_scale > 0)
  n <- nrow(X_tr); p <- ncol(X_tr)
  ss <- safe_scale_fit(X_tr); Xs <- safe_scale_apply(X_tr, ss$center, ss$scale)
  y_std <- (as.numeric(y_tr) - y_center) / y_scale
  
  # initial residuals at λ0 = arg from score mean (unit loadings)
  u0 <- y_std - mean(y_std)
  S0 <- Xs * u0
  grad0 <- colMeans(S0); lambda0 <- c_pi * max(abs(grad0))
  
  fit0 <- glmnet(Xs, y_std, alpha = 1, lambda = lambda0,
                 penalty.factor = rep(1, p),
                 intercept = TRUE, standardize = FALSE)
  u <- y_std - as.numeric(predict(fit0, Xs, s = lambda0))
  
  # Ω with Bartlett and fixed Q_T
  S <- Xs * u
  if (full_Omega) {
    Omega <- hac_lrv_matrix_kernel(S, L = Q_T, kernel = "bartlett", finite_adjust = TRUE)
    Omega <- (Omega + t(Omega)) / 2
    ee <- eigen(Omega, symmetric = TRUE)
    Omega <- ee$vectors %*% diag(pmax(ee$values, 1e-12), p) %*% t(ee$vectors)
    Z <- MASS::mvrnorm(B, mu = rep(0, p), Sigma = Omega)
  } else {
    sdcol <- sqrt(hac_var_diag_truncated(S, L = Q_T, kernel = "bartlett", finite_adjust = TRUE))
    Z <- matrix(stats::rnorm(B * p), nrow = B, ncol = p) * rep(sdcol, each = B)
  }
  q_sup <- unname(stats::quantile(apply(abs(Z), 1, max), 1 - alpha))
  lambda <- (c_pi * q_sup) / sqrt(n)
  
  fit <- glmnet(Xs, y_std, alpha = 1, lambda = lambda,
                penalty.factor = rep(1, p),
                intercept = TRUE, standardize = FALSE)
  
  list(
    fit = fit, lambda = lambda,
    x_center = ss$center, x_scale = ss$scale,
    y_center = y_center, y_scale = y_scale,
    diag = list(Q_T = Q_T, lambda0 = lambda0, q_sup = q_sup,
                n_active = sum(abs(as.vector(coef(fit, s = lambda))[-1]) > 0))
  )
}
#' @export
predict_asw <- function(obj, X_te) {
  Xs <- safe_scale_apply(as.matrix(X_te), obj$x_center, obj$x_scale)
  yhat_std <- as.numeric(predict(obj$fit, Xs, s = obj$lambda))
  yhat_std * obj$y_scale + obj$y_center
}

# ---------------- Ahrens HAC-LASSO: fixed q, kernel for φ ----------------
#' @export
fit_ahrens_fixed <- function(X_tr, y_tr,
                             C_PI = 0.8, q_fixed = 8,
                             iters = 3,
                             kernel_phi = c("uniform","bartlett","cosine"),
                             y_center, y_scale) {
  kernel_phi <- match.arg(kernel_phi)
  n <- nrow(X_tr); p <- ncol(X_tr)
  ss <- safe_scale_fit(X_tr); Xs <- safe_scale_apply(X_tr, ss$center, ss$scale)
  y_std <- (as.numeric(y_tr) - y_center) / y_scale
  
  # iid penalty level, then iterate φ with fixed bandwidth q_fixed
  gamma <- 0.1 / log(pmax(2, n))
  qcrit <- stats::qnorm(1 - gamma / (2 * p))
  lambda <- (2 * C_PI * qcrit) / sqrt(n)
  
  fit <- glmnet(Xs, y_std, alpha = 1, lambda = lambda,
                penalty.factor = rep(1, p),
                intercept = TRUE, standardize = FALSE)
  u <- y_std - as.numeric(predict(fit, Xs, s = lambda))
  
  for (it in seq_len(max(1, iters))) {
    phi <- sqrt(hac_var_diag_truncated(Xs * u, L = q_fixed,
                                       kernel = kernel_phi, finite_adjust = TRUE))
    fit <- glmnet(Xs, y_std, alpha = 1, lambda = lambda,
                  penalty.factor = phi,
                  intercept = TRUE, standardize = FALSE)
    u <- y_std - as.numeric(predict(fit, Xs, s = lambda))
  }
  
  list(
    fit = fit, lambda = lambda, phi = phi, q_fixed = q_fixed,
    x_center = ss$center, x_scale = ss$scale,
    y_center = y_center, y_scale = y_scale,
    diag = list(n_active = sum(abs(as.vector(coef(fit, s = lambda))[-1]) > 0),
                phi_min = min(phi), phi_med = median(phi), phi_max = max(phi))
  )
}
#' @export
predict_ahrens <- function(obj, X_te) {
  Xs <- safe_scale_apply(as.matrix(X_te), obj$x_center, obj$x_scale)
  yhat_std <- as.numeric(predict(obj$fit, Xs, s = obj$lambda))
  yhat_std * obj$y_scale + obj$y_center
}


################################################################################
                    # MODULAR CROSS FITTING WRAPPER 
################################################################################

#' @export
cf_predict <- function(X, y, folds,
                       method = c("ahrens","asw"),
                       # Ahrens
                       ah_C_PI = 0.8, ah_q_fixed = 8, ah_iters = 3, ah_kernel_phi = "uniform",
                       # ASW
                       asw_Q_T = NULL, asw_c_pi = 0.8, asw_alpha = 0.05, asw_B = 1000, asw_full_Omega = FALSE,
                       # Target scaling shared
                       scale_outcome = FALSE, mu = NULL, sdv = NULL) {
  method <- match.arg(method)
  Tn <- nrow(X); p <- ncol(X)
  oof <- rep(NA_real_, Tn)
  bet <- matrix(NA_real_, nrow = p, ncol = length(folds),
                dimnames = list(colnames(X), NULL))
  
  # one mean/sd for this target (keeps CF consistent)
  if (scale_outcome) {
    y_mu <- if (is.null(mu)) mean(y) else mu
    y_sd <- if (is.null(sdv)) sd(y)   else sdv
    if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  } else {
    y_mu <- 0; y_sd <- 1
  }
  
  diag_list <- vector("list", length(folds))
  
  for (k in seq_along(folds)) {
    tr <- folds[[k]]$train; te <- folds[[k]]$test
    Xtr <- X[tr,,drop=FALSE]; Xte <- X[te,,drop=FALSE]
    ytr <- y[tr]
    
    if (method == "ahrens") {
      f <- fit_ahrens_fixed(Xtr, ytr,
                            C_PI = ah_C_PI, q_fixed = ah_q_fixed,
                            iters = ah_iters, kernel_phi = ah_kernel_phi,
                            y_center = y_mu, y_scale = y_sd)
      yhat <- predict_ahrens(f, Xte)
    } else {
      if (is.null(asw_Q_T)) stop("ASW requires a global Q_T for this T.")
      f <- fit_asw_with_QT(Xtr, ytr, Q_T = asw_Q_T,
                           c_pi = asw_c_pi, alpha = asw_alpha, B = asw_B,
                           full_Omega = asw_full_Omega,
                           y_center = y_mu, y_scale = y_sd)
      yhat <- predict_asw(f, Xte)
    }
    
    oof[te] <- yhat
    cf <- as.vector(coef(f$fit, s = f$lambda))
    if (length(cf) > 1) {
      bet[, k] <- cf[-1]
    } else {
      bet[, k] <- rep(0, nrow(bet))  # all coefficients zero
    }
    
    diag_list[[k]] <- c(list(method = method), f$diag)
  }
  
  list(pred = oof, beta = bet,
       scale_center = y_mu, scale_scale = y_sd,
       diag = diag_list)
}



#Efficient HAC-weighted DML estimator (two-step GMM style)
#' @export
theta_hac_estimator <- function(psi_a_t, psi_b_t, L, kernel = "bartlett") {
  Psi <- cbind(psi_a_t, psi_b_t)
  Omega <- hac_lrv_matrix_kernel(Psi, L = L, kernel = kernel, finite_adjust = TRUE)
  W_hat <- tryCatch(solve(Omega), error = function(e) MASS::ginv(Omega))
  
  g_hat <- colMeans(Psi)  # (mean psi_a, mean psi_b)
  
  # θ = - (W_{12}*g_a + W_{22}*g_b) / (W_{11}*g_a + W_{12}*g_b)
  theta_hac <- - (W_hat[1,2] * g_hat[1] + W_hat[2,2] * g_hat[2]) /
    (W_hat[1,1] * g_hat[1] + W_hat[1,2] * g_hat[2])
  
  as.numeric(theta_hac)
}





################################################################################
                            # DML LOCAL PROJECTION
################################################################################

#' @export
dml_lp <- function(dt, y_var, shock_var,
                   H = 20, lags = 3, K = 10, buffer_folds = 1,
                   # θ(h) HAC settings
                   theta_kernel = c("bartlett","cosine","uniform"),
                   theta_bw     = c("fixed","andrews_ar1"),
                   L_nw = NULL,
                   # m(X) settings (choose method separately)
                   nuisance_m   = c("asw","ahrens"),
                   m_asw_c_pi   = 0.8,  m_asw_alpha = 0.05, m_asw_B = 1000, m_asw_full_Omega = FALSE, m_asw_aggL = "p95",
                   m_ah_C_PI    = 0.8,  m_ah_q_fixed = 8,   m_ah_iters = 5,   m_ah_kernel_phi = "uniform",
                   scale_outcome_m = FALSE,
                   # g(X) settings
                   nuisance_g   = c("ahrens","asw"),
                   g_asw_c_pi   = 0.4,  g_asw_alpha = 0.05, g_asw_B = 1000, g_asw_full_Omega = FALSE, g_asw_aggL = "p95",
                   g_ah_C_PI    = 0.5,  g_ah_q_fixed = 8,   g_ah_iters = 5,   g_ah_kernel_phi = "uniform",
                   scale_outcome_g = FALSE,
                   # identification
                   predetermined_y = TRUE
                   ) {
  
  theta_kernel <- match.arg(theta_kernel)
  theta_bw     <- match.arg(theta_bw)
  nuisance_m   <- match.arg(nuisance_m)
  nuisance_g   <- match.arg(nuisance_g)
  
  DATA <- build_matrix(dt, y_var, shock_var, H, lags)
  all_names <- names(DATA)
  lead_cols <- grep(paste0("^", y_var, "_lead"), all_names, value = TRUE)
  
  Y0    <- DATA[[y_var]]
  D_all <- DATA[[shock_var]]
  Tn    <- nrow(DATA)
  if (is.null(L_nw)) L_nw <- nw_lag_default(Tn)
  
  J_hat        <- numeric(H + 1)   # NEW
  omega_over_T <- numeric(H + 1)   # NEW
  theta <- numeric(H + 1)
  se_pt <- numeric(H + 1)
  ci_lo <- numeric(H + 1)
  ci_hi <- numeric(H + 1)
  
  
  folds <- create_nlo_folds(Tn, K, buffer_folds)
  
  # X designs
  exclude  <- c(lead_cols, shock_var)
  exclude0 <- c(exclude, y_var)
  X_mat_h  <- as.matrix(DATA[, setdiff(all_names, exclude),  with = FALSE])  # use for h>0 (and for g across h)
  X_mat_0  <- as.matrix(DATA[, setdiff(all_names, exclude0), with = FALSE])  # not used unless you want separate h=0
  
  # ---------- m(X): E[D|X], one fit (h-invariant design), reused for all h ----------
  if (nuisance_m == "asw") {
    d_mu <- if (scale_outcome_m) mean(D_all) else 0
    d_sd <- if (scale_outcome_m) sd(D_all)   else 1; if (!is.finite(d_sd) || d_sd==0) d_sd <- 1
    Q_T_m <- compute_QT_bandwidth(X_mat_h, D_all, mu = d_mu, sdv = d_sd, aggL = m_asw_aggL)
    M_X <- cf_predict(X = X_mat_h, y = D_all, folds = folds,
                      method = "asw",
                      asw_Q_T = Q_T_m, asw_c_pi = m_asw_c_pi,
                      asw_alpha = m_asw_alpha, asw_B = m_asw_B, asw_full_Omega = m_asw_full_Omega,
                      scale_outcome = scale_outcome_m, mu = d_mu, sdv = d_sd)
  } else {
    M_X <- cf_predict(X = X_mat_h, y = D_all, folds = folds,
                      method = "ahrens",
                      ah_C_PI = m_ah_C_PI, ah_q_fixed = m_ah_q_fixed, ah_iters = m_ah_iters,
                      ah_kernel_phi = m_ah_kernel_phi,
                      scale_outcome = scale_outcome_m)
  }
  m_hat  <- M_X$pred
  D_orth <- D_all - m_hat
  Betas_m <- M_X$beta
  
  # ---------- g(X): per horizon ----------
  theta <- J_hat <- se_pt <- ci_lo <- ci_hi <- numeric(H + 1)
  L_used <- integer(H + 1)
  Betas_g_all <- vector("list", H + 1)
  
  # If g uses ASW, compute ONE Q_T^g from (X_mat_h, Y0) and reuse across h
  if (nuisance_g == "asw") {
    g_mu0 <- if (scale_outcome_g) mean(Y0) else 0
    g_sd0 <- if (scale_outcome_g) sd(Y0)   else 1; if (!is.finite(g_sd0) || g_sd0==0) g_sd0 <- 1
    Q_T_g <- compute_QT_bandwidth(X_mat_h, Y0, mu = g_mu0, sdv = g_sd0, aggL = g_asw_aggL)
  } else {
    Q_T_g <- NULL
  }
  
  # diagnostics we’ll fill
  diag_theta <- vector("list", H + 1)
  
  for (h in 0:H) {
    if (predetermined_y && h == 0) {
      theta[1] <- 0; J_hat[1] <- 1; se_pt[1] <- 0; ci_lo[1] <- 0; ci_hi[1] <- 0
      Betas_g_all[[1]] <- matrix(0, nrow = ncol(X_mat_0), ncol = K,
                                 dimnames = list(colnames(X_mat_0), NULL))
      diag_theta[[1]] <- list(L = 0, score_var_over_T = 0, acf1 = 0, width = 0,
                              J_hat = 1, mean_Dorth2 = mean(D_orth^2))
      next
    }
    
    Yh <- if (h == 0) Y0 else DATA[[paste0(y_var, "_lead", h)]]
    
    if (nuisance_g == "ahrens") {
      G_X <- cf_predict(X = X_mat_h, y = Yh, folds = folds,
                        method = "ahrens",
                        ah_C_PI = g_ah_C_PI, ah_q_fixed = g_ah_q_fixed, ah_iters = g_ah_iters,
                        ah_kernel_phi = g_ah_kernel_phi,
                        scale_outcome = scale_outcome_g)
    } else {
      G_X <- cf_predict(X = X_mat_h, y = Yh, folds = folds,
                        method = "asw",
                        asw_Q_T = Q_T_g, asw_c_pi = g_asw_c_pi,
                        asw_alpha = g_asw_alpha, asw_B = g_asw_B, asw_full_Omega = g_asw_full_Omega,
                        scale_outcome = scale_outcome_g)
    }
    
    g_hat <- G_X$pred
    Betas_g_all[[h + 1]] <- G_X$beta
    
    # DML score pieces
    psi_a_t <- -(D_orth^2)
    psi_b_t <-  (D_orth * (Yh - g_hat))
    psi_a <- mean(psi_a_t); psi_b <- mean(psi_b_t)
    
    theta_h <- - psi_b / psi_a
    score_t <- psi_b_t + theta_h * psi_a_t
    Lh <- if (theta_bw == "fixed") L_nw else L_eq64_bartlett_vec(score_t)
    
    Om <- hac_lrv_matrix_kernel(matrix(score_t, ncol = 1), L = Lh,
                                kernel = theta_kernel, finite_adjust = TRUE)
    omega_over_T[h + 1] <- Om[1,1] / Tn
    se_h <- sqrt((Om[1,1] / Tn) / (psi_a^2))
    zc   <- qt(0.975, df = Tn - lags)
    
    theta[h + 1] <- theta_h
    J_hat[h + 1] <- -psi_a
    se_pt[h + 1] <- se_h
    ci_lo[h + 1] <- theta_h - zc * se_h
    ci_hi[h + 1] <- theta_h + zc * se_h
    L_used[h + 1] <- Lh
    
    acf1 <- if (length(score_t) > 1) stats::cor(score_t[-1], score_t[-length(score_t)]) else NA_real_
    diag_theta[[h + 1]] <- list(
      L = Lh,
      score_var_over_T = Om[1,1] / Tn,
      acf1 = acf1,
      width = (ci_hi[h + 1] - ci_lo[h + 1]),
      J_hat = J_hat[h + 1],
      mean_Dorth2 = mean(D_orth^2)
    )
  }
  
 #output
  intervals <- array(NA_real_,
                     dim = c(H + 1, 3, 1),
                     dimnames = list(horizon = as.character(0:H),
                                     information = c("lower","bhat","upper"),
                                     LRV = "intervals"))
intervals[, , 1] <- cbind(lower = ci_lo, bhat = theta, upper = ci_hi)

list(
  horizon = 0:H,
  theta   = theta,
  se      = se_pt,
  ci_pointwise = cbind(lower = ci_lo, upper = ci_hi),
  intervals = intervals,                  # <— matches your old extraction
  betas_m = Betas_m,
  betas_g = Betas_g_all,
  J_hat = J_hat,
  omega_over_T = omega_over_T,
  meta = list(
    nuisance_m = nuisance_m, nuisance_g = nuisance_g,
    scale_outcome_m = scale_outcome_m, scale_outcome_g = scale_outcome_g,
    theta_kernel = theta_kernel, theta_bw = theta_bw, L_used = L_used
  ),
  diagnostics = list(theta = diag_theta, m = M_X$diag) # g diag is per-horizon if you need it: compute similarly
)
}


#' @export
dml_applicattion <- function(dt, y_var, shock_var,
                   H = 20, lags = 3, K = 10, buffer_folds = 1,
                   B_boot = 500,
                   scale_outcome_m = FALSE,
                   scale_outcome_g = FALSE,
                   predetermined_y = TRUE) {
  
  # Build data with lags/leads
  DATA <- build_matrix(dt, y_var, shock_var, H, lags)
  all_names <- names(DATA)
  lead_cols <- grep(paste0("^", y_var, "_lead"), all_names, value = TRUE)
  
  Y0    <- DATA[[y_var]]
  D_all <- DATA[[shock_var]]
  Tn    <- nrow(DATA)
  
  folds <- create_nlo_folds(Tn, K, buffer_folds)
  
  # Exclude variables appropriately
  exclude  <- c(lead_cols, shock_var)
  exclude0 <- c(exclude, y_var)
  X_mat_h  <- as.matrix(DATA[, setdiff(all_names, exclude),  with = FALSE])
  X_mat_0  <- as.matrix(DATA[, setdiff(all_names, exclude0), with = FALSE])
  
  # ---------- m(X): E[D|X] ----------
  d_mu <- if (scale_outcome_m) mean(D_all) else 0
  d_sd <- if (scale_outcome_m) sd(D_all)   else 1; if (!is.finite(d_sd) || d_sd==0) d_sd <- 1
  Q_T_m <- compute_QT_bandwidth(X_mat_h, D_all, mu = d_mu, sdv = d_sd, aggL = "p95")
  M_X <- cf_predict(X = X_mat_h, y = D_all, folds = folds,
                    method = "asw",
                    asw_Q_T = Q_T_m, asw_c_pi = 0.8,
                    asw_alpha = 0.05, asw_B = 1000, asw_full_Omega = TRUE,
                    scale_outcome = scale_outcome_m, mu = d_mu, sdv = d_sd)
  m_hat  <- M_X$pred
  D_orth <- D_all - m_hat
  
  # ---------- g(X): setup ----------
  g_mu0 <- if (scale_outcome_g) mean(Y0) else 0
  g_sd0 <- if (scale_outcome_g) sd(Y0)   else 1; if (!is.finite(g_sd0) || g_sd0==0) g_sd0 <- 1
  Q_T_g <- compute_QT_bandwidth(X_mat_h, Y0, mu = g_mu0, sdv = g_sd0, aggL = "p95")
  
  # Storage
  theta_hat <- numeric(H + 1)
  ci_lo     <- numeric(H + 1)
  ci_hi     <- numeric(H + 1)
  
  # ---------- Loop over horizons ----------
  for (h in 0:H) {
    if (predetermined_y && h == 0) {
      theta_hat[1] <- 0
      ci_lo[1]     <- 0
      ci_hi[1]     <- 0
      next
    }
    
    Yh <- if (h == 0) Y0 else DATA[[paste0(y_var, "_lead", h)]]
    
    # g(X)
    G_X <- cf_predict(X = X_mat_h, y = Yh, folds = folds,
                      method = "asw",
                      asw_Q_T = Q_T_g, asw_c_pi = 0.4,
                      asw_alpha = 0.05, asw_B = 1000, asw_full_Omega = TRUE,
                      scale_outcome = scale_outcome_g)
    g_hat <- G_X$pred
    
    # Score components
    psi_a_t <- -(D_orth^2)
    psi_b_t <-  (D_orth * (Yh - g_hat))
    psi_a   <- mean(psi_a_t); psi_b <- mean(psi_b_t)
    
    # Point estimate
    theta_hat[h + 1] <- - psi_b / psi_a
    
    # ---------- Bootstrap ----------
    boot_thetas <- numeric(B_boot)
    for (b in 1:B_boot) {
      idx <- sample(seq_len(Tn), replace = TRUE)
      psi_a_b <- mean(psi_a_t[idx])
      psi_b_b <- mean(psi_b_t[idx])
      boot_thetas[b] <- - psi_b_b / psi_a_b
    }
    
    ci_lo[h + 1] <- quantile(boot_thetas, 0.025, na.rm = TRUE)
    ci_hi[h + 1] <- quantile(boot_thetas, 0.975, na.rm = TRUE)
  }
  
  # ---------- Return ----------
  list(
    horizon = 0:H,
    theta   = theta_hat,
    ci      = data.frame(
      lower = ci_lo,
      upper = ci_hi,
      row.names = 0:H
    )
  )
}

