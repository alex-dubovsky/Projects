#### DATA PREPROCESSING ####

# Reorder Data
Process <- function(data) {
  data <- data[-1] # drop time columns
  
  nm  <- names(data)
  idx <- grep("^(IR3TM|IRT3M)(_|$)", nm)   # match 3 month interest rate variables
  
  if (length(idx) >= 1) {
    # move matching column(s) to the end (keeps their order if >1)
    ord  <- c(setdiff(seq_along(nm), idx), idx)
    data <- data[, ord, drop = FALSE]
  }
  data
}






# Stationarity Checks

stationarity_table <- function(data) {
  res <- matrix(NA_character_, nrow = ncol(data), ncol = 3,
                dimnames = list(colnames(data), c("KPSS","PP","ADF")))
  
  for (j in seq_len(ncol(data))) {
    x <- data[[j]]
    
    kpss <- ur.kpss(x, type = "tau", lags = "long")                   # H0: (trend) stationarity
    pp   <- ur.pp(x,   type = "Z-tau", model = "trend", lags = "long") # H0: unit root
    adf  <- ur.df(x,   type = "drift", lags = 20, selectlags = "AIC")  # H0: unit root
    
    k_cv <- if (is.matrix(kpss@cval)) kpss@cval[1, "5pct"] else kpss@cval["5pct"]
    p_cv <- if (is.matrix(pp@cval))   pp@cval[1, "5pct"]   else pp@cval["5pct"]
    a_cv <- if (is.matrix(adf@cval))  adf@cval[1, "5pct"]  else adf@cval["5pct"]
    
    res[j, "KPSS"] <- if (as.numeric(kpss@teststat[1]) > k_cv) "Reject" else "Fail to Reject"
    res[j, "PP"]   <- if (as.numeric(pp@teststat[1])   < p_cv) "Reject" else "Fail to Reject"
    res[j, "ADF"]  <- if (as.numeric(adf@teststat[1])  < a_cv) "Reject" else "Fail to Reject"
  }
  
  as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
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
  
  # Check that each factor matches to each variable
  if (nrow(Lambda) != p) {
    Lambda <- Lambda[seq_len(min(nrow(Lambda), p)), seq_len(k), drop = FALSE]
  }
  
  # compute R^2 for each variable against each factor
  r2_mat <- matrix(NA_real_, nrow = p, ncol = k,
                   dimnames = list(colnames(data),
                                   paste0("Factor", seq_len(k))))
  
  for (i in seq_len(k)) {
    for (j in seq_len(p)) {
      r2_mat[j, i] <- summary(lm(data[[j]] ~ F_hat[, i]))$r.squared # How much variation does each factor explain in each variable
    }
    
    if (make_plots) {
      # R^2 bar plot for each factor
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
