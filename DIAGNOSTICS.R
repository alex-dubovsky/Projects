
.extract_intervals <- function(obj) {
  if (!is.null(obj$ci_pointwise)) {
    cbind(lower = obj$ci_pointwise[, "lower"],
          bhat  = obj$theta,
          upper = obj$ci_pointwise[, "upper"])
  } else if (!is.null(obj$intervals)) {
    arr <- obj$intervals
    if (is.array(arr) && length(dim(arr)) == 3) {
      cbind(lower = arr[, "lower", 1, drop = TRUE],
            bhat  = arr[, "bhat",  1, drop = TRUE],
            upper = arr[, "upper", 1, drop = TRUE])
    } else if (is.matrix(arr)) {
      arr[, c("lower","bhat","upper")]
    } else stop("Unknown 'intervals' format.")
  } else stop("No intervals found in object.")
}

# horizon-level rows (theta, width, se, J_hat, omega/T, L_used) — DML_LP only
.dml_diag_ci_rows <- function(fit, Tn) {
  ci <- .extract_intervals(fit)
  H1 <- nrow(ci)
  se   <- fit$se_pointwise %||% rep(NA_real_, H1)
  Jhat <- fit$J_hat        %||% (if (!is.null(fit$psi_a_bar)) -fit$psi_a_bar else rep(NA_real_, H1))
  Luse <- if (!is.null(fit$meta$L_used)) fit$meta$L_used else rep(NA_integer_, H1)
  OmT  <- if (all(is.finite(se)) && all(is.finite(Jhat))) (se * Jhat)^2 else rep(NA_real_, H1)
  
  data.frame(
    T      = Tn,
    h      = seq_len(H1) - 1L,
    theta  = if (!is.null(fit$theta)) fit$theta else ci[, "bhat"],
    se     = se,
    width  = ci[, "upper"] - ci[, "lower"],
    J_hat  = Jhat,
    omega_over_T = OmT,
    L_used = Luse,
    check.names = FALSE
  )
}

# compact summary per DML_LP fit
.dml_diag_compact <- function(fit, Tn) {
  ci <- .extract_intervals(fit)
  w  <- ci[, "upper"] - ci[, "lower"]
  base <- mean(w[2:min(9, length(w))],  na.rm = TRUE)
  mid  <- mean(w[11:min(15, length(w))], na.rm = TRUE)
  data.frame(
    T          = Tn,
    mean_width = mean(w, na.rm = TRUE),
    w_h1       = if (length(w) >= 2)  w[2]  else NA_real_,
    w_h12      = if (length(w) >= 13) w[13] else NA_real_,
    hump_mid   = is.finite(base) && is.finite(mid) && (mid > 1.2 * base)
  )
}

# simple selection tables for DML_LP betas
.sel_table <- function(beta_mat, thr = 0) {
  if (is.null(beta_mat)) return(NULL)
  vnames <- rownames(beta_mat) %||% paste0("X", seq_len(nrow(beta_mat)))
  data.frame(
    var      = vnames,
    sel_freq = rowMeans(abs(beta_mat) > thr, na.rm = TRUE),
    avg_abs  = rowMeans(abs(beta_mat), na.rm = TRUE),
    check.names = FALSE
  )
}