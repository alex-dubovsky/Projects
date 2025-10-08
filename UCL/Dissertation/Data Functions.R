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

series_id = function(nm) sub("\\..*$", "", sub("_[^_]+$", "", nm))
coerce_date = function(v) {
  if (inherits(v,"Date")) return(v)
  if (inherits(v,"POSIXt")) return(as.Date(v))
  s = as.character(v)
  ym = grepl("^\\d{4}-\\d{2}$", s)
  as.Date(ifelse(ym, paste0(s,"-01"), s))
}


aggregate_monthlies_by_id = function(dt, date_col = "Time",
                                      avg_ids = monthly_stock_id,
                                      sum_ids = monthly_flow_id) {
  x = as.data.table(copy(dt))
  x[, (date_col) := coerce_date(get(date_col))] # assign time column as date object version of itself
  yy = as.integer(format(x[[date_col]], "%Y")) #extracts year numbers
  mm = as.integer(format(x[[date_col]], "%m")) # extracts month numbers
  qq = (mm - 1L) %/% 3L + 1L # sorts "months" into the quarters of the year
  x[, `:=`(qkey = paste0(yy,"Q",qq),
           qstart = as.Date(sprintf("%d-%02d-01", yy, (qq - 1L) * 3L + 1L)))] # creates two new columns, qkey "XXXXQX" and qstart "XXXX-XX-XX"
  
  series_cols = setdiff(names(x), c(date_col,"qkey","qstart")) # all columns except c(...)
  x[, (series_cols) := lapply(.SD, function(z) suppressWarnings(as.numeric(z))), .SDcols = series_cols] # replace series_cols with numeric versions
  
  long = data.table::melt(x, id.vars = c("qkey","qstart"), # change data format from wide to long (stacked series columns) using melt
                           measure.vars = series_cols,
                           variable.name = "var", value.name = "value",
                           variable.factor = FALSE)
  long[, id := series_id(var)] #add id column
  long[, rule := fifelse(id %chin% sum_ids, "sum", # identify data to keep (quarterly), or aggregate via sum or average
                         fifelse(id %chin% avg_ids, "avg", "keep"))]
  
  agg = long[, {
    v = value[!is.na(value)]
    n = length(v)
    out = if (rule[1] == "keep") { # if variable rule = keep, in accordance with quarterly data
      if (n >= 1L) v[1L] else NA_real_                        # quarterly -> keep data as is
    } else {                                                  # monthly aggregations
      if (n == 3L && data.table::uniqueN(v) == 1L) v[1L]      # guard against quarterly data misreported as monthly (same values in all months of the quarter)
      else if (n == 3L && rule[1] == "sum") sum(v)            # sum monthly data
      else if (n == 3L && rule[1] == "avg") mean(v)           # average monthly data
      else NA_real_                                           # 0/2 months -> NA
    }
    .(value = out) # return value column with values of "out"
  }, by = .(qstart, var, rule)]
  
  wide = data.table::dcast(agg, qstart ~ var, value.var = "value") # rows - time, columns - series
  data.table::setnames(wide, "qstart", "Time") # rename q start as time
  data.table::setorder(wide, Time) # sort by time
  data.table::setcolorder(wide, c("Time", intersect(series_cols, names(wide)))) # reorder columns as original
  num_cols = setdiff(names(wide), "Time") # all series columns
  wide[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols] # replace series columns with numeric versions
  wide[]
}

aggregate_list_monthlies <- function(lst, date_col = "Time",
                                     avg_ids = monthly_stock_id,
                                     sum_ids = monthly_flow_id) {
  out = lapply(lst, aggregate_monthlies_by_id,
                date_col = date_col, avg_ids = avg_ids, sum_ids = sum_ids)
  names(out) = names(lst)
  out
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

pick_slow_fast = function(df, slow_prefixes, fast_prefixes) {
  vars = setdiff(names(df), c("Time","lag_FFR"))
  
  slow = vapply(vars, function(n) any(startsWith(n, slow_prefixes)), logical(1))
  fast = vapply(vars, function(n) any(startsWith(n, fast_prefixes)), logical(1))
  
  list(
    slow       = vars[ slow & !fast ],
    fast       = vars[ fast & !slow ],
    collisions = vars[ slow &  fast ],  # matched both sets
    unknown    = vars[!slow & !fast ]   # matched neither
  )
}


###### BAI AND NG, DFM with EM Algorithm for compuation of missing values/ragged edges ##########
imputation = function(X, undo_scale = FALSE, maxiter = 1000, tol = 1e-6, verbose = FALSE) {
                         
  X = as.matrix(X)                               # ensure numeric matrix
  miss = !is.finite(X)                           # logical matrix of missing cells (NA/NaN/Inf)
  
  mu = colMeans(X, na.rm = TRUE)                # column means (ignore NAs)
  if (undo_scale) {
    sdv = rep(1, ncol(X))                        # no scaling if undo_scale=TRUE
    XS = sweep(X, 2, mu, "-")                    # center: X - mu
  } else {
    sdv = apply(X, 2, sd, na.rm = TRUE)          # column sds
    sdv[!is.finite(sdv) | sdv == 0] = 1          
    XS = sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")# standardize: (X - mu)/sd
  }
  
  # zero on working scale = mean on original scale
  XSf = XS
  XSf[miss] = 0 # put 0s only where data are missing
  rvec = IC_BaiNg(XSf)
  q  = max(1L, as.integer(rvec))
  
  # EM loop
  it = 0L
  err = Inf    
  denom = sum(XSf^2)
  
  while (it < maxiter && err > tol) {             # iterate until convergence or maxiter
    pf = pc_lean(XSf, undo_scale = TRUE)        # M-step: your PCA on the completed matrix
    chi = pf$F_hat[, 1:q, drop = FALSE] %*%      # common component: F_q * Lambda_q'
      t(pf$Lambda_hat[, 1:q, drop = FALSE])
    
    XSn = XSf                                    # start from current completed matrix
    XSn[miss] = chi[miss]                        # update ONLY missing cells to common component
    err = sum((XSn - XSf)^2) / denom             # relative change for convergence check
    XSf = XSn                                    
    it = it + 1L                                 # increment iteration
    if (verbose && it %% 25 == 0)                 
      message(sprintf("iter %d, err=%.2e", it, err))
  }
  
  # revert to original scale
  if (undo_scale) {
    sweep(XSf, 2, mu, "+")                        # add means back (center-only case)
  } else {
    sweep(sweep(XSf, 2, sdv, "*"), 2, mu, "+")    # unstandardize: *sd then +mu
  }
}

