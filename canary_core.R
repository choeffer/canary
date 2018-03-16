#########################################################################
##
## Offline mode of CANARY: A Water Quality Event Detection Tool
## https://software.sandia.gov/trac/canary
## implemented in R. At the moment its not fully implemented. Pattern
## library and baseline change are missing.
##
## Version 1.0
## Copyright (c) 2018 Christian Hoeffer
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
#########################################################################

#define function offline_mode with given default values which can be called from another R script
offline_mode <- function(loaded_data=NULL,vals=NULL,window_size=2000,treshold=1.5,used_algo=lm,bed_window=5,prob_bed=0.01,bed_tresh=0.975,debug_flag=FALSE,debug_number_iteration=2040,...){
  
  #how to pass variables to a function 
  #https://cran.r-project.org/doc/manuals/r-devel/R-intro.html#Named-arguments-and-defaults
  #https://cran.r-project.org/doc/manuals/r-devel/R-intro.html#The-three-dots-argument
  
  ### check loaded_data
  #check if data object is empty or not passed to the function
  if(is.null(loaded_data)){
    stop("No data object loaded_data is passed to the function.")
  }
  
  #check if passed data object is a data.frame
  if(!is.data.frame(loaded_data)){
    stop("Function needs data object loaded_data in the form of a data.frame to run properly")
  }

  #check if one or more NA values are in the data.frame; anyNA instead of any(is.na) is used, see
  #https://stackoverflow.com/questions/6551825/fastest-way-to-detect-if-vector-has-at-least-1-na/35713234#35713234
  if(anyNA(loaded_data)){
    warning("Given loaded_data data.frame contains one ore more NA values!")
  }
  
  #check if first column of the data object is of data type factor which indicates the timestamp column
  if(is.factor(loaded_data[0,1])){
    #remove timestamps from the data object (mostly first column), it will be added later again
    subset_loaded_data <- loaded_data[,-1]
  }
  else{
    stop("Unexpected structure in the given loaded_data data.frame. Function offline_mode of canary_core.R expects data type factor in the first column of loaded_data indicating that the first column contains the timestamps")
  }
  
  ### check vals
  #check if vals is empty or not passed to the function
  if(is.null(vals)){
    stop("No data object vals of column names is passed to the function.")
  }
  
  #check if passed data object is a vector
  if(!is.vector(vals)){
    stop("Function needs data object vals in the form of a vector to run properly")
  }
  
  #check if one or more elements in the vector vals are NOT column names of the loaded_data
  if(any(!(is.element(vals,colnames(loaded_data))))){
    stop("One or more elements in the vector vals are NOT column names of the loaded_data data.frame. In this case function offline_mode will not work properly.")
  }
  
  #check if one or more elements of the vector vals are duplicated. This will break the the creation of the data frames and also later causes trouble if assigning results to one specific column name.
  if(any(duplicated(vals))){
    stop("One or more elements of the vector vals are duplicated. In this case function offline_mode will not work properly.")
  }
  
  #check if debugging/testing mode is selected or not
  if(debug_flag==TRUE){
    #limit number_iteration to 2030 (default, set by debug_number_iteration) and just loop over the first 2030 rows of the given loaded_data data.frame 
    #for fast testing different used_algo or algoAddParam on a loaded_data or other testing of this function
    number_iteration<-debug_number_iteration
  }
  if(debug_flag==FALSE){
    #loop over the whole length of the given loaded_data data.frame which means analyze all rows
    number_iteration<-nrow(loaded_data)
  }
  
  #check if window size + bed_windows is smaller number_iteration; if not the function offline_mode will not work properly
  if((window_size+bed_window)>number_iteration){
    stop("window size + bed_window=",window_size+bed_window," is bigger than number_iteration=",number_iteration,". In this case function offline_mode will not work properly.")
  }

  ### BED (Binomial Event Discriminator)
  #define a the function for BED see: https://www.osti.gov/servlets/purl/1267022
  binevdis <- function(r,n,p){
    part1 <- factorial(n)/(factorial(r)*factorial(n-r))
    part2 <- (p^r)*(1-p)^(n-r)
    return(part1*part2)
  }

  ############### create data.frames which are needed in the function ###############

  ### create empty data.frame as a blueprint for adding rows of the loaded_data while filling the window
  subset_used <- subset_loaded_data[0,]
  
  ### create empty data.frame as a blueprint for adding the timestamps of the loaded_data while filling the window
  subset_used_timestamps <- loaded_data[0,1, drop = FALSE]
  
  #create empty data frame for residuals and standard deviation of each column, defined in vals, in the window
  temp_res_sd<-setNames(as.data.frame(matrix(ncol = length(vals), nrow = 2)),vals)
  row.names(temp_res_sd)<-c("res","sd")
  
  ### create data.frame for residuals of each selected column (defined in vals)
  #create empty variable
  used_names_res <- NULL
  #fill variable with the needed names for the residuals data.frame
  #https://stackoverflow.com/questions/6984796/how-to-paste-a-string-on-each-element-of-a-vector-of-strings-using-apply-in-r
  used_names_res <- paste0(vals, "_residual")
  #create data.frame for residuals for each selected column and assign the needed column names to the residual data.frame
  #https://stackoverflow.com/a/32712555
  calc_residuals <- setNames(as.data.frame(matrix(ncol = length(used_names_res), nrow = 0)),used_names_res)
  #fill columns XXX_residual with NA until first calculation (depending on the window size)
  calc_residuals[1:window_size,] <- as.numeric(NA)
  
  ### create data.frame for if event or not of each selected column (defined in vals)
  #create empty variable
  used_names_events <- NULL
  #fill variable with the needed names for the column specific event data.frame
  used_names_events <- paste0(vals, "_event")
  #create data.frame for outlier detection saving for each selected column and assign the needed column names to the column specific event data.frame
  events_col <- setNames(as.data.frame(matrix(ncol = length(used_names_events), nrow = 0)),used_names_events)
  #fill columns XXX_event with NA until first calculation (depending on the window size)
  events_col[1:window_size,] <- NA

  ### create data.frame for each row if event or not (like canary does it)
  #create data.frame for outlier detection saving
  events <- read.csv(text="Event")
  #fill column Event with NA until first calculation (depending on the window size)
  events[1:window_size,] <- NA
  
  ### create data.frame for BED for each selected column (defined in vals)
  #create empty variable
  used_names_bed <- NULL
  #fill variable with the needed names for the column specific probability of an event data.frame and assign the needed column 
  #names to the column specific probability of an event data.frame
  used_names_bed <- paste0(vals, "_prob_event")
  #create data.frame for BED results saving for each column used
  events_col_bed <- setNames(as.data.frame(matrix(ncol = length(used_names_bed), nrow = 0)),used_names_bed)
  #fill columns XXX_prob_event with NA until first calculation (depending on the window_size and bed_window)
  events_col_bed[1:(window_size+bed_window),] <- as.numeric(NA)
  
  ### create data.frame for each row for BED probability (like canary does it)
  #create data.frame for saving the probabilities of an event
  events_bed <- read.csv(text="Prob_event")
  #fill column Event with NA until first calculation (depending on the window_size and bed_window)
  events_bed[1:(window_size+bed_window),] <- as.numeric(NA)
  
  ### create data.frame for each row for final decision if event or not
  #create data.frame for saving the probabilities of an event
  events_final <- read.csv(text="Final_event")
  #fill column Event with NA until first calculation (depending on the window size)
  events_final[1:(window_size+bed_window),] <- NA

  ############### loop over the number of rows (defined by number_iteration) in the loaded_data data.frame ###############
  
  #fill subset_used until window is filled with enough observations
  subset_used <- subset_loaded_data[1:(window_size-1),]
  #subset_used_return will be used later for returning the original values from the loaded_data
  #because the data in subset_used will be modified during the for loop
  subset_used_return <- subset_loaded_data[1:number_iteration,]
  #fill subset_used_timestamps
  subset_used_timestamps <- loaded_data[1:number_iteration,1,drop = FALSE]

  #start loop of the canary offline mode with adding each time a new row/element to the window and perform the calculations
  for(row in window_size:number_iteration){
    
    #add new row to subset_used and subset_used_print
    subset_used <- rbind(subset_used,subset_loaded_data[row,])

    #reset all values of temp_res_sd for next loop
    temp_res_sd[]<-NA
    
    #take the last observations of subsed_used which fit the window_size, so the observations which are in the window
    defined_window <- tail(subset_used,window_size)
    
    #scale() returns NaN in some cases if window size is too small, see https://stackoverflow.com/questions/15363610/why-does-scale-return-nan-for-zero-variance-columns
    #and therefore windowssize has to be big enough to avoid this!
    
    #normalize data; scale returns a matrix object, therefore the object is transformed back to a data.frame
    scaled_defined_window <- as.data.frame(scale(defined_window))
    
    #check if one or more NaN values are in the scaled_defined_window (as.matrix is needed because this function does not work with a data.frame)
    if(any(is.nan(as.matrix(scaled_defined_window)))){
       stop("Normalized data in the window contains NaN values. This might crash the used_algo, therefore function offline_mode stops")
    }
    
    #remove last added row
    scaled_defined_window_last_row_rm <- head(scaled_defined_window, -1)

    #exclude the last added row
    scaled_defined_window_last_row <- tail(scaled_defined_window, 1)
    
    #calculations for the residuals and events for all selected columns (defined in vals) 
    for(i in 1:length(vals)){
      #assign name
      used_name <- vals[i]
      #train selected algorithm (defined by used_algo);formula has to be inserted as ONE text, therefore paste is used and in addation 
      #as.formula() is used to transform the text to a formula so that it is recognized by algorithms as type formula
      #unlist() unpacks the additional parameters which can be passed to the selected algorithm
      used_model <- used_algo(formula = as.formula(paste0(used_name, " ~ . ")), data = scaled_defined_window_last_row_rm,...)
      pred_value <- predict(used_model, scaled_defined_window_last_row)
      #calculate residual between predicted and real value. unname() is used because return value is a named num and abs() is used to get absolut value
      temp_res_sd["res",used_name] <- abs(unname(pred_value-scaled_defined_window_last_row[,used_name]))
      #assign residual to the apppropriate column in the calc_residuals data.frame
      calc_residuals[row,paste0(used_name, "_residual")] <- temp_res_sd["res",used_name]
      #calculate sd for this column in the window
      temp_res_sd["sd",used_name]<-sd(defined_window[,used_name],na.rm = TRUE)
      #print(sd(defined_window[,"B_TOC_VAL"],na.rm=TRUE))

      #check if outlier or not for each column separately (residual classification) and store it in the events_col data.frame
      #additional functionality which is added and is not in CANARY
      #check for NA values because they brake the if statement
      if(!is.na(temp_res_sd["res",used_name])){
        #outlier
        if(temp_res_sd["res",used_name]>treshold*temp_res_sd["sd",used_name]){
          #add event = TRUE to events_col data.frame
          events_col[row,paste0(used_name,"_event")]=TRUE
        }
        #backround
        if(temp_res_sd["res",used_name]<=treshold*temp_res_sd["sd",used_name]){
          #add event = FALSE to events_col data.frame
          events_col[row,paste0(used_name,"_event")]=FALSE
        }
      }
      else{
        #add event = NA to events_col data.frame
        events_col[row,paste0(used_name,"_event")]=NA
      }
      
      ### BED analysis for each column
      #check if all of the last rows in events_col, defined by the bed_window, contain no NA values -> so are valid for BED analysis
      #additional functionality which is added and is not in CANARY itself
      if(!any(is.na(tail(events_col[,paste0(used_name,"_event")],bed_window)))){
        #how often was an event detected in the bed_window? event = TRUE https://stackoverflow.com/questions/2190756/how-to-count-true-values-in-a-logical-vector
        occur_col_true<-sum(tail(events_col[,paste0(used_name,"_event")],bed_window),na.rm=TRUE)
        #Calculation of the probability of an event, therefore 1-binevdis
        events_col_bed[row,paste0(used_name, "_prob_event")]<-1-binevdis(occur_col_true,bed_window,prob_bed)
      }
    }
    #find biggest residual in temp_res_sd and compare it with treshold
    index_biggest_res<-which.max(temp_res_sd["res",])
    max_rs<-temp_res_sd["res",index_biggest_res]
    #check if event or not (Result outlier Classification)
    #check for existence of NA values because they break the if statement
    if(length(max_rs) > 0){
      #outlier
      if(max_rs>treshold*temp_res_sd["sd",index_biggest_res]){
        #add event = TRUE to events data.frame
        events[row,"Event"] <- TRUE
        #remove value from the column of subset_used which is in this row the source for the outlier classification so that it is excluded from
        #the values used to predict water quality at next time step needs to be checked if assigning there NA is the proper way to solve it or not
        subset_used[row,colnames(temp_res_sd[index_biggest_res])]<-NA
      }
      #backround
      if(max_rs<=treshold*temp_res_sd["sd",index_biggest_res]){
        #add event = FALSE to events data.frame
        events[row,"Event"] <- FALSE
        
      }
    }
    else{
      #add event = NA to events data.frame
      events[row,"Event"] <- NA
    }
    
    ### BED analysis for each row (also final decision if event or not)
    #check if all of the last rows in events, defined by the bed_window, contain no NA values -> so are valid for BED analysis
    if(!any(is.na(tail(events[,"Event"],bed_window)))){
      #how often was an event detected in the bed_window? event = TRUE https://stackoverflow.com/questions/2190756/how-to-count-true-values-in-a-logical-vector
      occur_true<-sum(tail(events[,"Event"],bed_window),na.rm=TRUE)
      #Calculation of the probability of an event, therefore 1-binevdis
      temp_prob_event <- 1-binevdis(occur_true,bed_window,prob_bed)
      events_bed[row,"Prob_event"] <- temp_prob_event
      ### check if probability of an event is bigger a chosen threshold (defined by bed_tresh)
      #if yes, mark it as TRUE -> event
      if(temp_prob_event>=bed_tresh){
        events_final[row,"Final_event"]<-TRUE
      }
      #if not, mark it as FALSE -> NO event
      if(temp_prob_event<=bed_tresh){
        events_final[row,"Final_event"]<-FALSE
      }
    }
  }
  
  #combine all data.frames to one big data.frame for the results return value
  results <- cbind(subset_used_timestamps,subset_used_return,calc_residuals,events_col,events,events_col_bed,events_bed,events_final)
  
  #return the results data frame
  return(results)
}
