### Data PreProcessing from Adamek et al - with adjustments for Barigozzi et al 2024 ###
# functions related to processing of data

individual_transform<-function(column, code){
  Tn<-length(column)
  if(code==0){   
    ret<-column 
  }else if(code==4){
    ret<-as.numeric(c(NA, base::diff(column,na.pad = T, differences=1))) 
  }else if(code==5){
    ret<-as.numeric(c(NA,NA, base::diff(column, na.pad=T, differences=2)))
  }else if(code==1){
    ret<-as.numeric(log(column)*100)
  }else if(code==2){
    ret<-as.numeric(c(NA, base::diff(log(column), na.pad=T, differences=1)*100))
  }else if(code==3){
    ret<-as.numeric(c(NA,NA, base::diff(log(column), na.pad=T, differences=2)*100))
  }else if(code==7){
    ret<-rep(NA,Tn)
    for(i in 2:Tn){
      ret[i]<-column[i]/column[i-1]-1
    }
    ret<-as.numeric(ret)
  }else{warning("not a valid code")}
  return(ret)
}

data_transform<-function(mat, transform_codes=NULL){
  if(is.null(transform_codes)){
    transform_codes<-mat[1,] #takes the first row of our data matrix as the transformation codes. 
    transform_codes[1]<-1 # assuming Time/Date variable is first, leaves this untransformed
    transform_codes<-as.numeric(transform_codes) # makes the row vector numeric
  }
  transform_codes[1]<-1 #if transform codes are not NULL skip the above, first value of transform codes = 1 (Time/Date untransformed), 
  ret<-mat[-c(1),] #removes the first row of the data matrix. Assumes the first row of the matrix is the transformation codes
  for(i in 1:ncol(mat)){ #for each variable in the data matrix
    ret[,i]<-individual_transform(column=as.numeric(unlist(mat[-1, i])), code=transform_codes[i]) #column = loop over the variables in the data
  }
  return(ret)
}

#' @export
group_which<-function(full_set, subset){
  ind<-1:length(subset)
  for(i in 1:length(subset)){
    if(length(which(full_set==subset[i]))>0){
      ind[i]<-which(full_set==subset[i])
    }else{
      ind[i]<-NA
    }
    
  }
  return(ind)
}

#' @export
clean_data<-function(raw_data, slow_names, FFR_name, fast_names, start_date, end_date, transform_codes=NULL){
  manually_transformed<-fred_transform(mat = raw_data, transform_codes = transformation_codes)
  var_names<-colnames(manually_transformed)
  dates<-raw_data[-1,1]
  variables_we_use<-c(slow_names, FFR_name, fast_names) # minus the date/Time variable it should be 118 variables
  indexes<-group_which(var_names, variables_we_use)
  FFR=manually_transformed[,FFR_name]
  data_all_dates<-data.frame(cbind(dates,manually_transformed[,indexes]), lag_FFR=dplyr::lag(FFR))
  data_all_dates$Time = as.Date(data_all_dates$Time)
  start<-which(data_all_dates$Time==start_date)
  end<-which(data_all_dates$Time==end_date)
  data_in_window<-data_all_dates[start:end,]
  
  delete<-NULL
  for(i in 1:nrow(data_in_window)){
    if(sum(anyNA(data_in_window[i,]))){
      delete<-c(delete, i)
    }else{
      break
    }
  }
  for(i in nrow(data_in_window):1){
    if(sum(anyNA(data_in_window[i,]))){
      delete<-c(delete, i)
    }else{
      break
    }
  }
  delete <- unique(delete)
  
  if(length(delete)>0){
    data_all<-data_in_window[-delete,]
    warning(paste0("Deleted ", length(delete), " rows with missing data"))
  }else{
    data_all<-data_in_window
  }
  lag_FFR<-data_all$lag_FFR
  data_sans_dates<-data_all[,-c(1, which(colnames(data_all)=="lag_FFR"))]
  missing_plot<-NULL
  missing_summary<-NULL
  if(sum(is.na(data_all))>0){
    warning("Internal missing values!")
    missing_values<-is.na(data_all)
    missing_plot<-graphics::image(missing_values, col=c("white", "red"))
    missing_summary<-apply(data_all[-delete,],FUN=function(x){sum(is.na(x))},MARGIN=2)
  }
  return(list(data_all=data_all,
              data_sans_dates=data_sans_dates,
              var_names=var_names,
              slow_data=data_sans_dates[,group_which(variables_we_use, slow_names)],
              fast_data=data_sans_dates[,group_which(variables_we_use, fast_names)],
              FFR=data_sans_dates[,group_which(variables_we_use, FFR_name)],
              missing_plot=missing_plot,
              missing_summary=missing_summary,
              start_date=data_all$dates[1],
              end_date=data_all$dates[nrow(data_all)],
              dates=data_all$dates,
              lag_FFR=lag_FFR
  ))
}
