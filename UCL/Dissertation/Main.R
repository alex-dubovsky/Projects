
install.packages(pkgs = c("dplyr", "ggpattern", "ggplot2", "ggpubr", "readxl", "reshape2","bigtime","xtable","here"))
devtools::install_github("cykbennie/fbi")
install.packages(pkgs = c("remotes","vars","readx",'ggpatterm','dplyr','ggplot2','ggpubr','reshape2','bigtime','xtable','stats','rrpack','MTS','urca','
         MCMCpack','tictoc'))
library(DoubleML)
library(devtools)
library(remotes)
library(vars)
library(readxl)
library(fbi)
library(ggpattern)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(bigtime)
library(xtable)
library(stats)
library(rrpack)
library(MTS)
library(urca)
library(Matrix)
library(MCMCpack)
library(tictoc)
library(parallel)
library(doParallel)
library(foreach)
library(tseries)
library(missMethods)
library(norm)
library(data.table)
library(glmnet)
library(here)
source(here("UCL/Dissertation", "DFM_Functions.R")) 
source(here("UCL/Dissertation", "Data Functions.R")) 
source(here("UCL/Dissertation", "PreProcessing Functions.R"))
setwd(here("UCL/Dissertation/Outputs"))

############################  European Data ####################################

path = here("UCL/Dissertation/Data")
files = list.files(path, pattern = "\\.xlsx?$", full.names = TRUE)
all_list <- lapply(files, read_excel)
stems = tools::file_path_sans_ext(basename(files)) 
names(all_list) = stems

monthly_stock_id = c( # prefixes of the variables in each dataset which are stock variables
  "IRT3M","IRT6M","LTIRT",
  # (7) Industrial production & turnover indexes
  "IPMN","IPCAG","IPDOG","IPIDCOG","IPINDCOG","IPING","IPINRG",
  "TRNNM","TRNCAG","TRNCOG","TRNDCG","TRNNDCG","TRNING","TRNNRG",
  # (8) Prices
  "PPICAG","PPICOG","PPIDCOG","PPINDCOG","PPING","PPINRG",
  "HICPOV","HICPNEF","HICPCG","HICPIN","HICPSV","HICPNG",
  # (9) Confidence indicators
  "ICONFIX","CCONFIX","ESENTX","KCONFIX","RTCONFIX","SCONFIX","BCI","CCI",
  # (10) Monetary aggregates (levels)
  "CURR","M1","M2",
  # (11) Others
  "SHIX"
)
monthly_flow_id = c(
  "CAREG"
)


# Quarterly aggregation 
Data_list =  aggregate_list_monthlies(all_list, date_col = "Time")
Mat_list <- lapply(Data_list, as.data.frame)

for (nm in names(Mat_list)) {                     # loop over each dataset in your list
  df  = Mat_list[[nm]]                          # take one data.frame
  vals = as.matrix(df[, 2:ncol(df), drop = FALSE]) # drop date col; keep numeric block
  df[, 2:ncol(df)] = suppressWarnings(imputation(               # impute missing/ragged cells in the block
    vals, undo_scale = FALSE
  )
  )
  
  Mat_list[[nm]] = df                          # write back to the list
}

na_rows = sapply(Mat_list, function(d) sum(rowSums(is.na(d[,-1,drop=FALSE])) > 0))



# Adding EA-wide 3 month interest rate to each quarterly dataset
merge_obj = Mat_list[[4]][,c("Time","IRT3M_EACC")]
for (i in seq_along(Mat_list)){
  if (i == 4) next
  else Mat_list[[i]] = merge(Mat_list[[i]],merge_obj,by = "Time",all.x = TRUE,sort=FALSE)
}

EA_LT = c(
  # (1) National Accounts (1–17)
  0,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,
  # (2) Labor Market Indicators (18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (39–64);;
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (65–72)
  2,2,2,2,2,2,2,2,
  # (5) Exchange Rates (73–74)
  2,2,
  # (6) Interest Rates (75–77)
  4,4,4,
  # (7) Industrial Production and Turnover (78–91)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (92–105)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (106–113)
  0,0,0,0,0,0,2,2,
  # (10) Monetary Aggregates (114–116)
  2,2,2,
  # (11) Others (117–118)
  2,2
)
DE_LT = c(
  # (1) National Accounts (1–17)
  0,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (39–61)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (62–68)
  2,2,2,2,2,2,2,
  # (5) Exchange Rates (69)
  2,
  # (6) Interest Rates (70)
  4,
  # (7) Industrial Production and Turnover (71–84)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (85–98)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (99–106)
  0,0,0,0,0,0,2,2,
  # (10) Others (107)
  2,4
)

FR_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–17)
  2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 39–63)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 64–71)
  2,2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 72)
  2,
  # (6) Interest Rates (Series 73)
  4,
  # (7) Industrial Production and Turnover (Series 74–87)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 88–101)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 102–109)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 110)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

IT_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–14)
  2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 15–35)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 36–52)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 61–67)
  2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 68)
  2,
  # (6) Interest Rates (Series 69)
  4,
  # (7) Industrial Production and Turnover (Series 70–82)
  2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 83–96)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 97–104)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 105)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

ES_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–14)
  2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 15–35)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 36–52)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 58–64)
  2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 65)
  2,
  # (6) Interest Rates (Series 66)
  4,
  # (7) Industrial Production and Turnover (Series 67–79)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 81–94)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 95–102)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 103)
  2,
  # Extra variable (IRT3M_EACC)
  4
)
  
PT_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–14)
  2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 15–35)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 36–56)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 61–68)
  2,2,2,2,2,2,2,2, 
  # (5) Exchange Rates (Series 69)
  2,
  # (6) Interest Rates (Series 70)
  4,
  # (7) Industrial Production and Turnover (Series 71–83)
  2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 84–90)
  2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 91–98)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 99)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

NL_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–17)
  2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 39–56)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 64–71)
  2,2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 72)
  2,
  # (6) Interest Rates (Series 73)
  4,
  # (7) Industrial Production and Turnover (Series 74–87)
  2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 88–100)
  2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 101–108)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 109)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

AT_LT <- c(
  0,  # Time
  # (1) National Accounts (1–17)
  2,2,2,2,2,2,2,2,2,2,2,2,2,  # 1–13
  0,0,0,0,                     # 14–17
  # (2) Labor Market Indicators (18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (39–62)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (63–70)
  2,2,2,2,2,2,2,2,             # 63–70
  # (5) Exchange Rates (71)
  2,                           # 71
  # (6) Interest Rates (72)
  4,                           # 72
  # (7) Industrial Production & Turnover (73–86)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2, # 73–86
  # (8) Prices (87–100)
  2,2,2,2,2,2,                 # 87–92
  2,2,2,2,2,2,                 # 93–98
  2,2,                         # 99–100
  # (9) Confidence Indicators (101–108)
  0,0,0,0,0,0,                 # 101–106
  2,2,                         # 107–108
  # (11) Others (109)
  2,                           # 109
  # Extra last code
  4                            # IRT3M_EACC
)

EL_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–12)
  2,2,2,2,2,2,2,2,2,2,0,0,
  # (2) Labor Market Indicators (Series 13–33)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 39–49)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 54–61)
  2,2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 62)
  2,
  # (6) Interest Rates (Series 63)
  4,
  # (7) Industrial Production and Turnover (Series 64–76)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 77–89)
  2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 90–97)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 98)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

IE_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–17)
  2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,
  # (2) Labor Market Indicators (Series 18–38)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 39–46)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 57–62)
  2,2,2,2,2,2,
  # (5) Exchange Rates (Series 63)
  2,
  # (6) Interest Rates (Series 64)
  4,
  # (8) Prices (Series 65–72)
  2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 73–80)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 81)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

BE_LT = c(
  # Time column
  0,
  # (1) National Accounts (Series 1–14)
  2,2,2,2,2,2,0,0,0,
  # (2) Labor Market Indicators (Series 15–35)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,2,
  # (3) Credit Aggregates (Series 35–56)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (4) Labor Costs (Series 57–64)
  2,2,2,2,2,2,2,2,
  # (5) Exchange Rates (Series 65)
  2,
  # (6) Interest Rates (Series 66)
  4,
  # (7) Industrial Production and Turnover (Series 67–79)
  2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (8) Prices (Series 80–92)
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  # (9) Confidence Indicators (Series 94–101)
  0,0,0,0,0,0,2,2,
  # (11) Others (Series 102)
  2,
  # Extra variable (IRT3M_EACC)
  4
)

slow_prefixes = c("ULC","GDP","EXPGS","GFC","HFCE","CONS","GCF","GFA","AHRDI",
                  "AHFCE","GNFC","GHIR","GHSR","TEMP","EMP","SEMP","THOURS","EMP",
                  "UNE","RPRP","WS","ESC","IP","TR","CAREG","IMP")
fast_prefixes = c("TASS.SDB","TASS.LBD","TLB","NFC","GGASS","GGLB","HHASS","HHLB","REER42","ERUS",
                  "IRT6M","LTIRT","PPI","HICP","DFGDP","HPRC","ICONFIX","CCONFIX",
                  "ESENTIX","KCONFIX","RTCONFIX","SCONFIX","BCI","CCI", "CURR","M1","M2","SHIX", "TASS.SLN","TASS.LLN")
Slow = setNames(vector("list", length(Mat_list)), names(Mat_list))
Fast = setNames(vector("list", length(Mat_list)), names(Mat_list))
for (nm in names(Mat_list)) { # IRT3M_EACC should be unclassified in all datasets
  cls = pick_slow_fast(Mat_list[[nm]], slow_prefixes, fast_prefixes)
  Slow[[nm]] = cls$slow
  Fast[[nm]] = cls$fast
  if (length(cls$unknown))
    cat("Unclassified in", nm, ":", paste(cls$unknown, collapse = ", "), "\n")
  else print("No unclassified variables")
}

Codes = list(
  EAdata = EA_LT,
  DEdata = DE_LT,
  FRdata = FR_LT,
  ITdata = IT_LT,
  ELdata = EL_LT,
  PTdata = PT_LT,
  ESdata = ES_LT,
  IEdata = IE_LT,
  NLdata = NL_LT,
  BEdata = BE_LT,
  ATdata = AT_LT
)
Cleaned = list()
for (nm in names(Mat_list)) {
  df = Mat_list[[nm]]
  codes = Codes[[nm]]
  Cleaned[[nm]] = clean_data(
    raw_data = df,
    slow_names = Slow[[nm]],
    FFR_name = "IRT3M_EACC",
    fast_names = Fast[[nm]],
    transform_codes = codes
  )
  
  cat("Finished cleaning", nm, "\n")
}





Stationary_Checks <- setNames(lapply(all_datasets, function(nm) stationarity_table(get(nm))),
                   all_datasets)
for (i in length(Stationary_Checks)){
  print(Stationary_Checks[[i]])
}

n_f <- 8 # number of factors determined by Bai and Ng information critera
p_f <-1 # vector autoregressive order of the factors, we found that by BIC this is equal to 1
q_v <- 1 # autoregressive order of the idiosyncratic errors
S<-1:119# indices of the variables included in estimating the factors. Useful for determining the position of the policy variable which has been ordered last
h_max=20 #  maximum horizon - the impulse response function is evaluated from horizon 0 to h_max
set.seed(1) # We use set.seed(9) in the calibration

########################## DENSE DYNAMIC FACTOR MODEL ##########################
#factors by pc, VAR by OLS, can take a few minutes to run
Dense_DFM<-Estimate_DFM(Data, n_f = n_f, lags_f = p_f, lags_v = q_v, max_EV = 0.98, undo_scale=TRUE, factor_method = "pc", VAR_method="OLS") 
Dense_IRF<-impulse_response_ABCD(Dense_DFM$factors, Dense_DFM$idio, S, 20, policy_var = length(S), outcome_var = 1) # Industrial Production

#evidence for sparsity
r2_mat = matrix(NA, nrow = 118, ncol = 8,
                 dimnames = list(colnames(Data), paste0("Factor", 1:8)))
for (i in 1:ncol(Dense_DFM$factors$F_hat)) {
  fac_name = colnames(r2_mat)[i]
  for (j in 1:118) {
    reg <- lm(Data[, j] ~ Dense_DFM$factors$F_hat[, i])
    r2_mat[j, i] = summary(reg)$r.squared
  }
  Plotdata = data.frame(variable = colnames(Data),R2 = as.numeric(r2_mat[,i]))
  print(
    ggplot(Plotdata, aes(x = variable, y = R2)) +
    geom_bar(stat = "identity") +
    labs(
      title = bquote(R^2 ~ " for regression of each variable on" ~.(fac_name)),
      x = "Variable",
      y = expression(R^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  print(
    ggplot(data.frame(variable = colnames(Data),FactorLoading = as.numeric((Dense_DFM$factors$Lambda[,i]))) , aes(x = variable, y = FactorLoading)) +
    geom_bar(stat = "identity") +
    labs(
      title = bquote("Factor Loadings for" ~.(fac_name)),
      x = "Variable",
      y = "Factor Loading"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
return(r2_mat)
}

#################################################### SPARSE DYNAMIC FACTOR MODEL #########################################################################

#Sparse_DFM<-Estimate_DFM(Data,n_f=n_f, lags_f = 1, lags_v = 1, max_EV = 0.98, undo_scale=TRUE, factor_method = "sWF", VAR_method="lasso")
#Sparse_IRF<-impulse_response_ABCD(Sparse_DFM$factors, Sparse_DFM$idio, S, 20, policy_var = length(S), outcome_var = 1) [1,]

Sparse_IRF = readRDS("Sparse_IRF.RData")
Sparse_DFM = readRDS("Sparse_DFM.RData")


cor_mat <- cor(Data, use = "pairwise.complete.obs")
heatmap(cor_mat)
rk(Data)
#If you see values close to 1, you have strong collinearity.


r2_sparse = matrix(NA, nrow = 119, ncol = 8,
                dimnames = list(colnames(Data), paste0("Factor", 1:8)))
for (i in 1:ncol(Dense_DFM$factors$F_hat)) {
  fac_name = colnames(r2_sparse)[i]
  for (j in 1:118) {
    reg <- lm(Data[, j] ~ Sparse_DFM$factors$F_hat[, i])
    r2_mat[j, i] = summary(reg)$r.squared
  }
  Plotdata = data.frame(variable = colnames(Data),R2 = as.numeric(r2_mat[,i]))
  print(
    ggplot(Plotdata, aes(x = variable, y = R2)) +
      geom_bar(stat = "identity") +
      labs(
        title = bquote(R^2 ~ " for regression of each variable on" ~.(fac_name)),
        x = "Variable",
        y = expression(R^2)
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  print(
    ggplot(data.frame(variable = colnames(Data),FactorLoading = as.numeric((Dense_DFM$factors$Lambda[,i]))) , aes(x = variable, y = FactorLoading)) +
      geom_bar(stat = "identity") +
      labs(
        title = bquote("Factor Loadings for" ~.(fac_name)),
        x = "Variable",
        y = "Factor Loading"
      ) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  return(r2_mat)
}

# Saving DFM Data Locally

saveRDS(Dense_DFM, file="Dense_DFM.RData")
saveRDS(Dense_IRF, file="Dense_IRF.RData")
saveRDS(Sparse_DFM, file="Sparse_DFM.RData")
saveRDS(Sparse_IRF, file="Sparse_IRF.RData")

### Agnostic Identification for Shocks ### - after 10,000 iterations, unable to find contemporaneous/rotation matrix such that Sparse CPI is negative

F_hat = Sparse_DFM$factors$F_hat
Phi_hat = Sparse_DFM$factors$Phi
Sigma_eta = cov(Sparse_DFM$VAR_resid)
neg_index = c(1:4)
shock_id = 119
n1=1   
n2=10   
lag = 1

keep_Sigma <- list()
keep_H     <- list()
keep_IRF <-list()
n_rotations = 10
seed=9
factors = sWF_lasso_DFM$factors
idio = sWF_lasso_DFM$idio
Sigma_eta = cov(sWF_lasso_DFM$VAR_resid)
agnostic_H_all <- function(n_rotations, factors, idio, Sigma_eta, S,
                           h_max = 20, policy_var, outcome_var, neg_index,seed=2){ 
  # Load parallel back-end
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required for parallel execution")
  }
  library(doParallel)
  cores <- parallel::detectCores(logical = FALSE)
  cl <-makeCluster(cores - 1)
  clusterEvalQ(cl,{source("C:/Users/AlexD/Documents/Economics/UCL/Dissertation/R Code/DFM_Functions.R")})
  clusterExport(cl,varlist=c("impulse_response_H_given","neg_index"))
  registerDoParallel(cl)
  
  m <- ncol(Sigma_eta)
  P <- t(chol(Sigma_eta)) #lower cholesky
  
  # Parallel search over rotations
  results <- foreach(s = seq_len(n_rotations), .combine = 'c') %dopar% {
    # Draw random orthogonal matrix Q with det(Q) = +1
    seeds = sample(1e6)
    set.seed(seeds)
    Z <- matrix(rnorm(m * m), m, m)
    QR_decomp <- qr(Z)
    Q_temp <- qr.Q(QR_decomp)
    R_temp <- qr.R(QR_decomp)
    Q <- Q_temp %*% diag(sign(diag(R_temp)))
    
    Hcand <- P %*% Q
    IRF0 <- impulse_response_H_given(factors, idio, S,
                                     h_max = h_max,
                                     policy_var = policy_var,
                                     outcome_var = outcome_var,
                                     H = Hcand)
    # Check sign restrictions
    if (all(IRF0[1:length(neg_index), 2:6] <= 0) &&
        all(IRF0[length(neg_index)+1, 2:(h_max + 1)] > 0)) {
      list(list(H = Hcand, IRF = IRF0))
    } else {
      NULL
    }
  }
  
  stopCluster(cl)
  
  # results is a list of lists; flatten if necessary
  valid_results <- Filter(Negate(is.null), results)
  if (length(valid_results) == 0) {
    warning(sprintf("No rotations found after %d attempts.", n_rotations))
    return(NULL)
  }
  on.exit(toc(),add=FALSE)
  # Extract H matrices and IRFs
  H_list   <- lapply(valid_results, `[[`, "H")
  IRF_list <- lapply(valid_results, `[[`, "IRF")
  
  return(list(H_all = H_list, IRF_all = IRF_list))
  
}
tic()
agnostic_H = agnostic_H_all(n_rotations = 1000,Sparse_DFM$factors,Sparse_DFM$idio,Sigma_eta = Sigma_eta,S,h_max = 20,policy_var = 119,outcome_var = 5,neg_index)
toc()

#trial and error/adapt appropriate sigma.omega (error variance of IV) until all F-statistics are approximately within F = (10,30) which are empirically valid instruments.
f = iv_calibration(rho = c(0,0.2,0.5),n_rep = 100,sigma.omega =c(1.6,2.1,2.6),init=50)
f



##### Simulations #####


# ----- DML grid -----
cpi_vals <- c(0.4,0.5,0.6,0.7)

dml_grid <- expand.grid(
  nuisance_m = "asw",
  nuisance_g = "asw",
  
  m_asw_c_pi = cpi_vals,
  g_asw_c_pi = cpi_vals,   # will filter to equal pairs below
  
  m_asw_alpha = 0.05, g_asw_alpha = 0.05,
  m_asw_B = 1000,     g_asw_B     = 1000,
  m_asw_aggL = "p95", g_asw_aggL  = "p95",
  
  scale_outcome_m = FALSE,
  scale_outcome_g = FALSE,
  
  theta_bw = "andrews_ar1",
  theta_kernel = "bartlett",
  stringsAsFactors = FALSE,
  K = 15
)

# keep only identical (m,g) c_pi pairs: 
dml_grid <- dml_grid[dml_grid$m_asw_c_pi == dml_grid$g_asw_c_pi, , drop = FALSE]

# (optional) neat labels
dml_grid$labels <- with(dml_grid, paste0(
  "ASW_c", m_asw_c_pi,
  "_sm", ifelse(scale_outcome_m, "T","F"),
  "_sg", ifelse(scale_outcome_g, "T","F"),
  "_K",K
))


# ----- simulation settings -----
M  <- 750
Ts <- c(200, 400, 600)
set.seed(2)
seeds <- data.frame(matrix(sample(x = 1:1e8, size = M * 1), ncol = 1))
colnames(seeds) <- c("Sparse")
# ----- cluster -----
threads <- parallel::detectCores() - 2
t0 <- Sys.time()
cl <- parallel::makeCluster(threads)
t1 <- Sys.time()
cat("Cluster setup:", t1 - t0, "\n")

# --- load packages on workers
parallel::clusterEvalQ(cl, {
  library(data.table)
  library(glmnet)
  library(MASS)
  library(DissertationSimPack)
})

setup <- "Sparse"
t2 <- Sys.time()
cat("Packages loaded on workers:", t2 - t1, "\n")

# export data objects
parallel::clusterExport(
  cl,
  varlist = c("dml_grid","Ts","Sparse_DFM","Sparse_IRF"),
  envir   = .GlobalEnv
)
t3 <- Sys.time()
cat("Data export:", t3 - t2, "\n")
cat("Starting", setup, "simulation...\n")
setup_start <- Sys.time()
simulation <- parallel::parLapply(
  cl,
  X   = seeds[[setup]],
  fun = one_replication_lean_manual_seed_safe,
  Ts  = Ts,
  DFM = get(paste0(setup,"_DFM")),
  IRF = get(paste0(setup,"_IRF")),
  dml_grid = dml_grid
)
setup_end <- Sys.time()
cat(setup, "simulation finished in", setup_end - setup_start, "\n")
cat("Replications returned:", length(simulation), "\n")
simulation <- Filter(Negate(is.null), simulation) 
cat("Valid replications:", length(simulation), "out of", M, "\n")
saveRDS(simulation, paste0(setup,"_sim.RData"))
parallel::stopCluster(cl)


source("C:/Users/AlexD/Documents/Economics/UCL/Dissertation/R Code/ML Script.R")

simdata <- readRDS(paste0(setup,"_sim.RData"))
length(simdata)
adamek_simdata = readRDS("C:/Users/AlexD/Documents/Economics/UCL/Dissertation/R Code/Output/sWF_lasso_sim.RData")

processed_sim = process_sim(simdata)
adamek_processed = process_sim(adamek_simdata)



your_subset <- processed_sim
your_subset$coverages     <- your_subset$coverages[, , c("DML_2","DML_3"), , drop = FALSE]
your_subset$median_widths <- your_subset$median_widths[, , c("DML_2","DML_3"), , drop = FALSE]
your_subset$DGPs          <- c("DML_2","DML_3")


needed_LRVs <- adamek_processed$LRVs
your_subset$LRVs <- needed_LRVs

pad_array <- function(arr, needed_LRVs){
  old_LRVs <- dimnames(arr)$LRV
  new_arr <- array(NA_real_,
                   dim = c(dim(arr)[1], dim(arr)[2], dim(arr)[3], length(needed_LRVs)),
                   dimnames = list(
                     horizon = dimnames(arr)$horizon,
                     T_      = dimnames(arr)$T_,
                     DGP     = dimnames(arr)$DGP,
                     LRV     = needed_LRVs
                   ))
  # Fill what exists
  for(lrv in old_LRVs){
    new_arr[,,,lrv] <- arr[,,,lrv]
  }
  new_arr
}

your_subset$coverages     <- pad_array(your_subset$coverages, needed_LRVs)
your_subset$median_widths <- pad_array(your_subset$median_widths, needed_LRVs)


combined_proc <- your_subset
combined_proc$coverages     <- abind::abind(adamek_processed$coverages,
                                            your_subset$coverages, along = 3)
combined_proc$median_widths <- abind::abind(adamek_processed$median_widths,
                                            your_subset$median_widths, along = 3)
combined_proc$DGPs <- c(adamek_processed$DGPs, your_subset$DGPs)


source("C:/Users/AlexD/Documents/Economics/UCL/Dissertation/R Code/DFM_Functions.R")

plots <- plot_processed_sim(combined_proc, models = combined_proc$DGPs)

model_colors <- c(
  "HDLP_04" = "blue",
  "FALP"    = "lightblue",
  "LP"      = "purple",
  "DML_2"   = "red",
  "DML_3"   = "orange"
)

# fixed color palette matching exactly combined_proc$DGPs
model_colors <- c(
  "HDLP_04" = "blue",
  "LP"      = "darkgreen",
  "FALP"    = "darkorange",
  "DML_2"   = "red",
  "DML_3"   = "purple"
)

plot_processed_intervals <- function(processed_sim, models = NULL) {
  horizon = model = value = NULL
  
  Ts <- processed_sim$Ts
  T_strings <- paste0("T_", Ts)
  h_max <- processed_sim$h_max
  
  if (is.null(models)) {
    models <- processed_sim$DGPs
  }
  
  plot_list <- vector("list", 2 * length(Ts))
  
  for (ts in T_strings) {
    # Coverage + median width data
    temp_cov <- data.frame()
    temp_mw  <- data.frame()
    
    for (dgp_name in models) {
      d1 <- data.frame(value = processed_sim$coverages[, ts, dgp_name, "intervals"],
                       horizon = 0:h_max, model = dgp_name)
      temp_cov <- rbind(temp_cov, d1)
      
      d2 <- data.frame(value = processed_sim$median_widths[, ts, dgp_name, "intervals"],
                       horizon = 0:h_max, model = dgp_name)
      temp_mw <- rbind(temp_mw, d2)
    }
    
    # ✅ enforce consistent factor levels
    temp_cov$model <- factor(temp_cov$model, levels = names(model_colors))
    temp_mw$model  <- factor(temp_mw$model,  levels = names(model_colors))
    
    # Coverage plot
    p_cov <- ggplot2::ggplot(temp_cov, aes(x = horizon, y = value, color = model)) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.93, ymax = 0.97),
                         fill = "grey90", alpha = 0.3, inherit.aes = FALSE) +
      ggplot2::geom_hline(yintercept = 0.95, color = "black", linetype = 2, linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 1.0, alpha = 0.9) +
      ggplot2::scale_color_manual(values = model_colors, name = "Model") +
      ggplot2::labs(title = paste0("T = ", Ts[which(T_strings == ts)]),
                    x = "Horizon", y = "Coverage")
    
    # Median width plot
    p_mw <- ggplot2::ggplot(temp_mw, aes(x = horizon, y = value, color = model)) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::geom_line(linewidth = 1.0, alpha = 0.9) +
      ggplot2::scale_color_manual(values = model_colors, name = "Model") +
      ggplot2::labs(title = paste0("T = ", Ts[which(T_strings == ts)]),
                    x = "Horizon", y = "Median Width")
    
    # Store plots
    plot_list[[which(T_strings == ts)]] <- p_cov
    plot_list[[length(Ts) + which(T_strings == ts)]] <- p_mw
  }
  
  # Arrange into grid
  combined <- ggpubr::ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "right")
  return(combined)
}




plot_processed_intervals(combined_proc,
                         models = c("HDLP_04","LP","FALP","DML_2","DML_3"))















############## VARIABLE FREQUENCY ##########################################

all_betas <- list()

for (rep in seq_along(simdata)) {
  for (spec in names(simdata[[rep]]$DML)) {
    fit <- simdata[[rep]]$DML[[spec]]
    
    # --- m(X): fold-averaged ---
    beta_m <- rowMeans(abs(fit$betas_m))
    if (is.null(names(beta_m))) names(beta_m) <- rownames(fit$betas_m)
    
    # --- g(X): fold- and horizon-averaged ---
    beta_g_list <- lapply(fit$betas_g, function(mat) {
      if (is.null(mat) || length(mat) == 0) return(NULL)   # skip empty horizons
      b <- rowMeans(abs(mat))
      if (is.null(names(b))) names(b) <- rownames(mat)
      b
    })
    beta_g_list <- beta_g_list[!sapply(beta_g_list, is.null)]
    
    if (length(beta_g_list) > 0) {
      beta_g <- Reduce(function(a, b) {
        all_vars <- union(names(a), names(b))
        a <- a[match(all_vars, names(a))]; a[is.na(a)] <- 0
        b <- b[match(all_vars, names(b))]; b[is.na(b)] <- 0
        setNames(a + b, all_vars)
      }, beta_g_list) / length(beta_g_list)
    } else {
      beta_g <- numeric(0)
    }
    
    # --- align by names ---
    all_vars <- union(names(beta_m), names(beta_g))
    beta_m <- beta_m[match(all_vars, names(beta_m))]; beta_m[is.na(beta_m)] <- 0
    beta_g <- beta_g[match(all_vars, names(beta_g))]; beta_g[is.na(beta_g)] <- 0
    
    df <- data.frame(
      replication = rep,
      spec        = spec,
      variable    = all_vars,
      beta_m      = beta_m,
      beta_g      = beta_g
    )
    
    all_betas[[length(all_betas) + 1]] <- df
  }
}

betas_df <- dplyr::bind_rows(all_betas)

## selection frequency ##
library(stringr)
library(dplyr)

# Add base_var column to collapse lags/leads
betas_df2 <- betas_df %>%
  mutate(base_var = str_replace(variable, "_lag\\d+$", ""),
         base_var = str_replace(base_var, "_lead\\d+$", ""))

# Recompute activity
activity_df <- betas_df2 %>%
  group_by(spec, base_var) %>%
  summarise(
    freq_m = mean(beta_m != 0),
    freq_g = mean(beta_g != 0),
    .groups = "drop"
  )
top_vars_m <- activity_df %>%
  group_by(spec) %>%
  slice_max(order_by = freq_m, n = 25)
library(dplyr)
library(ggplot2)
library(tidytext)  # for reorder_within()

k <- 20  # top-k variables per spec

top_vars_m <- activity_df %>%
  group_by(spec) %>%
  slice_max(order_by = freq_m, n = k, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    base_var = reorder_within(base_var, freq_m, spec),
    freq_pct = freq_m * 100
  )

p_m <- ggplot(top_vars_m, aes(x = freq_pct, y = base_var, fill = spec)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ spec, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_brewer(palette = "Set2") +   # << consistent palette
  labs(title = paste("Top", k, "Active Variables in m(X)"),
       x = "Selection Frequency (%)", y = "Variable") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 12)
  )

# ----- g(X) -----
top_vars_g <- activity_df %>%
  group_by(spec) %>%
  slice_max(order_by = freq_g, n = 40, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    base_var = reorder_within(base_var, freq_g, spec),
    freq_pct = freq_g * 100
  )

p_g <- ggplot(top_vars_g, aes(x = freq_pct, y = base_var, fill = spec)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ spec, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_brewer(palette = "Set2") +   # << same palette
  labs(title = paste("Top", k, "Active Variables in g(X)"),
       x = "Selection Frequency (%)", y = "Variable") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 12)
  )
p_m
p_g
     


##### empirical #####
source("C:/Users/AlexD/Documents/Economics/UCL/Dissertation/R Code/ML Script.R")
Data = as.data.table(Data)
y_var     <- colnames(Data)[1]
shock_var <- colnames(Data)[ncol(Data)]
set.seed(1)
other_slow<-as.matrix(CDbmedium$slow_data[,-c(which(colnames(CDbmedium$slow_data)=="INDPRO")
                                              ,which(colnames(CDbmedium$slow_data)=="CPIAUCSL"))]) 
fast <- as.matrix(CDbmedium$fast_data)

IP  <- as.matrix(CDbmedium$data_all[,"INDPRO",   drop = FALSE])
CPI <- as.matrix(CDbmedium$data_all[,"CPIAUCSL", drop = FALSE])
FFR <- as.matrix(CDbmedium$FFR, drop = FALSE)

slow_vars <- cbind(other_slow, IP, CPI)
fast_vars <- fast

dml_FFR <- dml_application_slowfast(
  data          = CDbmedium$data_all,
  slow_vars     = slow_vars,
  fast_vars     = fast_vars,
  y_var         = "FEDFUNDS",
  shock_var     = "FEDFUNDS",
  H             = 49,
  predetermined_y = FALSE,
  cumulate      = FALSE
)

dml_ahrens <- dml_application_slowfast(Data, slow_vars, fast_vars,
                                       y_var = "FEDFUNDS", shock_var = "FEDFUNDS",
                                       H = 49)

time1 = Sys.time()
DML_2_INDPRO = dml_application(dt=Data,y_var,shock_var,H=49,predetermined_y = TRUE,lags = 13, cumulate=TRUE)
DML_2_FFR = dml_application_slowfast(Data,slow_vars,fast_vars,y_var="FEDFUNDS",shock_var = "FEDFUNDS",H=49,predetermined_y = FALSE,lags = 13, cumulate=FALSE, scale_outcome_m = TRUE,scale_outcome_g = TRUE)
DML_2_CPI = dml_application(Data,y_var="CPIAUCSL",shock_var = "FEDFUNDS",H=49,predetermined_y = TRUE,lags = 13, cumulate=TRUE)
time2 = Sys.time()
time2-time1



 plot_irf <- function(irf_df, title = "Impulse Response",
                     line_color = "red", band_color = "red", alpha = 0.25) {
  ggplot2::ggplot(irf_df, ggplot2::aes(x = horizon, y = theta)) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.4, linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                         fill = band_color, alpha = alpha) +
    ggplot2::geom_line(color = line_color, linewidth = 1) +
    ggplot2::labs(title = title, x = "Horizon", y = "Response") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}
DML = dml_lp_final(Data,slow_vars,fast_vars,y_var = "FEDFUNDS", shock_var="FEDFUNDS",theta_kernel = "bartlett",theta_bw = "andrews_ar1",nuisance_m = "asw",nuisance_g="asw",scale_outcome_m = FALSE,scale_outcome_g = FALSE,predetermined_y = FALSE)
irf_df <- data.frame(
  horizon   = DML$horizon,
  theta     = DML$theta,
  ci_lower  = DML$ci_pointwise[, "lower"],
  ci_upper  = DML$ci_pointwise[, "upper"]
)

plot_irf(irf_df = DML_2_INDPRO)
plot_irf(irf_df = DML_2_CPI)
plot_irf(irf_df=DML_2_FFR)
dml_ahrens$theta
plot_irf(irf_df)








library(ggplot2); library(ggpubr); library(data.table)

as_dml_df <- function(obj) {
  # Handle either a data.frame directly or a list that contains one
  if (is.data.frame(obj) && all(c("horizon","theta","ci_lower","ci_upper") %in% names(obj))) {
    return(as.data.table(obj))
  }
  if (is.list(obj)) {
    for (el in obj) {
      if (is.data.frame(el) && all(c("horizon","theta","ci_lower","ci_upper") %in% names(el))) {
        return(as.data.table(el))
      }
    }
  }
  stop("Could not find DML output with columns horizon/theta/ci_lower/ci_upper.")
}

dml_plot <- function(df, ylab_txt) {
  ggplot(df, aes(x = horizon, y = theta)) +
    theme_bw() +
    geom_hline(yintercept = 0) +
    geom_line(color = "red") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.5, fill = "red") +
    xlab("Horizon") + ylab(ylab_txt)
}

# Expect these to be the returned data.frames from your DML calls
dml_ffr_df <- as_dml_df(DML_2_FFR)      # non-cumulative
dml_ip_df  <- as_dml_df(DML_2_INDPRO)   # cumulative (theta already cumulated in function)
dml_cpi_df <- as_dml_df(DML_2_CPI)      # cumulative

p7 <- dml_plot(dml_ffr_df, "FFR")
p8 <- dml_plot(dml_ip_df,  "IP")
p9 <- dml_plot(dml_cpi_df, "CPI")

top_row <- ggarrange(p7, p8, p9, nrow = 1)
top_row <- annotate_figure(top_row, top = text_grob("DML-2", color = "black", face = "bold", size = 14))

# --- combine all rows: DML (top), HDLP (middle), FAVAR (bottom) ------------
combined <- ggarrange(top_row, first_row, second_row, ncol = 1)

# (optional) save
ggsave(filename = "fig3_with_DML.pdf", plot = combined, device = "pdf",
       width = 18, height = 18, units = "cm", dpi = 1000)
       