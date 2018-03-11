#########################################################################
##
## Test file for canary_core.R
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

#set working directory, test_canary_core.R and canary_core.R should be in this directory
setwd("~/Dokumente/AIT_Master/3rd_semester/case_study_bartz/R_data/splitted_up")

#clear Environment (for testing/developing the code)
rm(list=ls())

#import the R script which includes the function offline_mode
source("canary_core.R")

param<-list(
  
  #load data set as a data.frame; format data.frame is mandatory
  #sample data can be found in canary-4.3-src.zip under /canarytrunk/examples/sample_data/
  #from https://software.sandia.gov/trac/canary/downloader/download
  loaded_data = read.csv("sample_data/test_station_b.csv"),
  
  #vector of column names of the loaded_data passed to the function to select which columns to be used for the analyzis
  vals = c("B_CL2_VAL","B_TURB_VAL","B_PH_VAL","B_TOC_VAL","B_COND_VAL","B_TEMP_VAL","B_PLNT_PH_VAL","B_PLNT_TURB_VAL","B_PLNT_CL2_VAL"),
  #vals = c("B_TEMP_VAL","B_PLNT_PH_VAL","B_PLNT_TURB_VAL","B_PLNT_CL2_VAL"),
  #vals = c("B_TEMP_VAL","B_PLNT_PH_VAL"),
  
  #define window size
  window_size = 2000,
  
  #define treshold for decision if outlier or not depending on residuals
  treshold = 1.5,
  
  #define which algorithm to use for prediction of last row of window
  used_algo = lm,
  
  #define optional parameters for the used algorithm defined with used_algo
  #this parameters will directly be passed to the used_algo function
  #model=TRUE,
  #singular.ok = TRUE,
  #rainer="test",
  
  ### BED (Binomial Event Discriminator) settings ###
  #define how many rows are taken for the BED function
  bed_window=10,
  #define probability for BED function
  prob_bed=0.01,
  #define treshold probability if event is detected or not
  bed_tresh=0.85,
  
  #debug/test mode; if debug_flag = FALSE, ALL rows of the data.frame are analyzed
  #if debug_flag = TRUE debug_number_iteration is used, so just the first XX rows of the data.frame are analyzed
  debug_flag=TRUE,
  debug_number_iteration=2030
)

#call of the function can be done with do.call(offline_mode,list) or offline_mode(parameters, ...)
#this are the minimal needed parameters for a proper function call
#final_results <- offline_mode(loaded_data,vals)

#call with do.cal would look like, where param is a list with all parameters, see above for an example
final_results <- do.call(offline_mode,param)

#write final result to a csv if wanted
#write.csv(final_results, "final_results.csv")

#print out the last x rows of the returned data.frame for debugging or just seeing the results
tail(final_results,35)