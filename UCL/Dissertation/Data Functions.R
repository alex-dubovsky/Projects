#### DATA PREPROCESSING ####

# Stationarity Checks

Stationarity_checks = function(data){
  ur.kpss(data,type = "tau",lags="long")
  
}







sparsity_checks <- function(
    data,
    nrows = NULL,                  
    ncols = NULL,                  
    dfm = NULL,                    # pass DFM object
    F_hat = NULL,                  # OR pass factors and loadings explicitly
    Lambda = NULL,
    make_plots = TRUE
) {
  # Factors from either DFM object or Manual imputation
  if (!is.null(dfm)) {
    F_hat  <- dfm$factors$F_hat
    Lambda <- dfm$factors$Lambda
  }
  if (is.null(F_hat) || is.null(Lambda))
    stop("Provide either `dfm` or both `F_hat` and `Lambda`.")
  
  # Optional - Use subset or full Dataset
  if (!is.null(nrows)) data <- data[seq_len(min(nrows, nrow(data))), , drop = FALSE]
  if (!is.null(ncols)) data <- data[,  seq_len(min(ncols, ncol(data))), drop = FALSE]
  
  # Safety check - keep only numeric columns
  num <- vapply(data, is.numeric, logical(1)) # applies is.numeric to each element(variable) in data. Logical(1), single True/False per variable
  if (!all(num)) data <- data[, num, drop = FALSE]
  
  # Dimensions
  obs <- nrow(data)
  p    <- ncol(data)
  k    <- ncol(F_hat)
  
  # Check observations/sample size of Factors and Dataset and align samples
  if (nrow(F_hat) != obs) {
    warn_msg <- sprintf("Aligning time dimension: data has %d rows, F_hat has %d rows. Using common min.",
                        obs, nrow(F_hat))
    warning(warn_msg)
    keep <- seq_len(min(obs, nrow(F_hat)))
    data <- data[keep, , drop = FALSE]
    F_hat <- F_hat[keep, , drop = FALSE]
  }
  
  # subset Lambda rows to variables kept (if it has more rows)
  if (nrow(Lambda) != p) {
    Lambda <- Lambda[seq_len(min(nrow(Lambda), p)), seq_len(k), drop = FALSE]
  }
  
  # --- compute R^2 for each variable against each factor ---
  r2_mat <- matrix(NA_real_, nrow = p, ncol = k,
                   dimnames = list(colnames(data),
                                   paste0("Factor", seq_len(k))))
  
  for (i in seq_len(k)) {
    for (j in seq_len(p)) {
      r2_mat[j, i] <- summary(lm(data[[j]] ~ F_hat[, i]))$r.squared
    }
    
    if (make_plots) {
      # R^2 bar plot for factor i
      plotdata <- data.frame(variable = colnames(data), R2 = r2_mat[, i])
      print(
        ggplot(plotdata, aes(x = variable, y = R2)) +
          geom_bar(stat = "identity") +
          labs(
            title = bquote(R^2 ~ " for regression of each variable on " ~ .(paste0("Factor ", i))),
            x = "Variable", y = expression(R^2)
          ) +
          theme_minimal(base_size = 14) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      )
      
      # factor loadings for factor i
      fl <- as.numeric(Lambda[seq_len(p), i])
      print(
        ggplot(data.frame(variable = colnames(data), FactorLoading = fl),
               aes(x = variable, y = FactorLoading)) +
          geom_bar(stat = "identity") +
          labs(
            title = paste0("Factor Loadings for Factor ", i),
            x = "Variable", y = "Factor Loading"
          ) +
          theme_minimal(base_size = 14) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      )
    }
  }
  
  return(r2_mat)
}
