# ABCD 2WayMOSUM/Kuldorff Simulation R Script

# This script includes functions which will perform the simulations comparing 2WayMOSUM from Li et al. (2024) and
# the spatial scan statistic from Kuldorff et al. (2009) automatically.

# Requires the following packages to perform simulations 
library(MASS)
library(ade4)
library(L2hdchange)
library(ggplot2)
library(rsatscan)

# Generating covariance matrix for image data, iterating through each pair of pixels
vcov_mat <- matrix(NA, nrow = 100, ncol = 100)
for(i in 1:10){
  for(j in 1:10){
    index1 <- (i-1)*10 + j
    row <- i
    col <- j
    for(k in 1:10){
      for(l in 1:10){
        index2 <- (k - 1)*10 + l
        compare_row <- k
        compare_col <- l
        dist <- sqrt((row - compare_row)^2 + (col-compare_col)^2)
        vcov_mat[index1, index2] <- (0.6)^dist
      }
    }
  }
}

## Important note: #####
# We set the seed of each simulation type to get repeatable results; we started with set.seed(1) for each trial. Unfortunately 2WayMOSUM was prone to fairly 
# frequent errors, so we had to change the seed to set.seed(2) in the case of trials with t_5 and dependent normal data for a full simulation without errors from 2WayMOSUM
#####

image_sim <- function(ntrials, data_type = c("Normal", "t", "Normal_cov"), sss = TRUE, cp_time = NULL){
  # ntrials is the number of trials per setting, a total of 3*ntrials runs are performed
  # data_type controls which type of data to be analyzed: standard Gaussian, t_5 distributed or spatially dependent Gaussian
  # sss=T includes space-time scan statistic in the simulation, 
  # cp_time can be set if one wants to set a custom time of change-point for a simulation
  
  ABCD_tests_full <- list()
  scan_tests_full <- list()
  l2_tests_full<- list()
  
  n <- 200 
  if(is.null(cp_time) == F){
    cplus1 <- cp_time + 1
  } else {
    cplus1 <- 121}
  # 2WayMOSUM temporal window set to its default of 30
  MOSUM_window <- 30
  
  changed_region1 <- c(1:5, 11:15, 21:25, 31:35, 41:45)
  changed_region2 <- c(1:4, 11:14, 21:24, 31:34)
  changed_region3 <- c(1:3, 11:13, 21:23)
  
  # Ultimately we set the magnitude of mean change per component to be the same for every data type for our simulation...
  delta_t5 <- sqrt(2/9) 
  delta_norm <- sqrt(2/9)
  delta_covnorm <- sqrt(2/9) 
  
  if(data_type == "Normal"){
    set.seed(1)
    delta <- delta_norm}
  if(data_type == "Normal_cov"){
    set.seed(2)
    delta <- delta_covnorm
  }
  
  if(data_type == "t"){
    set.seed(2)
    delta <- delta_t5
  }
  print(delta)
  
  # Making a list of the 3 changed regions (5x5, 4x4, 3x3)
  changed_regions <- list(changed_region1, changed_region2, changed_region3)
  
  for(l in 1:3){
    ABCD_tests <- list()
    scan_tests<- list()
    l2_tests <- list()
    
    changed_region <- unlist(changed_regions[[l]])
    
    for(k in 1:ntrials){
      print(k)
      # Generating the data in a flattened way
      if(data_type == "Normal"){
        test_data <- matrix(rnorm(n*100), ncol = n)
      } 
      if(data_type == "t"){
        test_data <- matrix(rt(n*100, df = 5), ncol = n)
      } 
      if(data_type == "Normal_cov"){
        # using covariance matrix as defined above
        test_data <- t(MASS::mvrnorm(n, mu = rep(0, 100), Sigma=vcov_mat))
      }
      #print(var(test_data[,1]))
      
      # Randomly selecting 9 dimensions within changed region to change
      changed <- sample(changed_region, 9)
      test_data[changed,cplus1:n] <- test_data[changed, cplus1:n] + delta
      
      # Generating an image array from the flattened data
      test_array <- array(NA, dim = c(10,10, n))
      for(i in 1:10){
        for(j in 1:10){
          index <- (i-1)*10 + j
          test_array[i,j, ] <- test_data[index,]
        }
      }
      
      ### ABCD Test (ran on the image array): #######
      ABCD_test<- ABCD_img(test_array, k = floor(n/5), stat.type = "m", B = 1000, perm.p = T,
                           P_1_grid = c(1, 2, 3), P_2_grid = c(1, 2, 3))
      ABCD_tests[[k]] <- ABCD_test
      print(ABCD_test$tauhat_max)
      print(ABCD_test$p_val)
      
      ### 2WayMOSUM Test: #######
      # equivalent blocks to to P_1 = (1,2,3), P_2 = (1,2,3) for ABCD, which has a total of 1 + 4 + 9 = 14 blocks to test over
      block_groups <- list(c(1:100), #1
                           c(1:5, 11:15, 21:25, 31:35, 41:45), #2 
                           c(51:55, 61:65, 71:75, 81:85, 91:95),#3
                           c(6:10, 16:20, 26:30, 36:40, 46:50),#4
                           c(56:60, 66:70, 76:80, 86:90, 96:100),#5
                           c(1:3, 11:13, 21:23), c(4:6, 14:16, 24:26), c(7:10, 17:20, 27:30),#6-8
                           c(31:33, 41:43, 51:53), c(34:36, 44:46, 54:56), c(37:40, 47:50, 57:60),#9-11
                           c(61:63, 71:73, 81:83, 91:93), c(64:66, 74:76, 84:86, 94:96), c(67:70, 77:80, 87:90, 97:100))#12-14
      y_obj <- ts_hdchange(test_data, window_size = MOSUM_window, 
                           # These are all standard settings for the 2WayMOSUM procedure...
                           m = 8, 
                           h = 1,
                           N_rep = 999, # Generation of significance threshold 
                           alpha = 0.05, # Nominal Type 1 error rate
                           nbd_info = block_groups) 
      
      l2_test <- hdchange(y_obj)
      l2_tests[[k]] <- l2_test
      print(l2_test$nbd_and_stamps_pair)
      
      ### Space-time Scan Statistic Test #######
      # Note: One needs to have SaTScan software installed on your machine to be able use RSatscan 
      # Code adapted from example found in an RSatscan vignette on CRAN: https://cran.r-project.org/web/packages/rsatscan/vignettes/simulation.html
      if(sss == TRUE){
        mygeo = expand.grid(1:10,1:10)
        daysbase = n # number of observations 
        locid = rep(1:100, times=daysbase) # Setting 100 locations per observation
        basecas = as.numeric(test_data) # Feeding in test data
        day = rep(1:n, each = 100)
        casecount <- rep(1, 100*n) # total number of observations
        mycas = data.frame(locid, casecount, day, basecas) # formatting the specific way Rsatscan takes in data
        td = tempdir() # As a wrapper function for SaTScan, one needs to write files like this to run a trial
        # We have opted to create a temporary directory to do this to make things manageable
        write.geo(mygeo, location = td, file = "mygeo", userownames=TRUE)
        write.cas(mycas, location = td, file = "mycas")
        
        # Code to specify how SaTScan is run, using a space-time scan for Normal data
        invisible(ss.options(reset=TRUE))
        ss.options(list(CaseFile="mycas.cas", PrecisionCaseTimes=4)) #
        ss.options(list(StartDate="1", CoordinatesType=0, TimeAggregationUnits=4)) # Specifying Cartesian Coordinates
        ss.options(list(EndDate=as.character(n), CoordinatesFile="mygeo.geo", AnalysisType=3, ModelType=5)) #Space-time scan with Normal model
        ss.options(list(NonCompactnessPenalty=0, MaxTemporalSizeInterpretation=0, MaxTemporalSize=50)) # Can find at maximum a cluster half the length of the series, no penalty for non-compactness
        ss.options(list(SaveSimLLRsDBase="y"))
        
        write.ss.prm(td, "mybase")
        # Running the actual SaTScan procedure... sslocation and SaTScanBatch64 may vary based on your SaTScan installation
        mybase = satscan(td, "mybase", sslocation="C:/Program Files/SaTScan", ssbatchfilename="SaTScanBatch64")
        
        # Saving output of each test, a matrix of detected clusters with associated start/end points and p-values
        scan_tests[[k]] <- mybase$col 
        print(mybase$col)
      } else {
        scan_tests[[k]] <- NULL
      }
    }
    ABCD_tests_full[[l]] <- ABCD_tests
    scan_tests_full[[l]] <- scan_tests
    l2_tests_full[[l]] <- l2_tests
  }
  
  # Returning the full output for all tests
  return(list(ABCD_tests_full = ABCD_tests_full, scan_tests_full = scan_tests_full,
              l2_tests_full = l2_tests_full))
}

# A function to aggregate the number of significant and accurate change-points found by each method, for a given setting
sig_acc_tbl <- function(sim, ntrials, sss = T, cp_time = NULL){
  ABCD_sig_acc <- matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  l2_sig_acc <- matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  scan_sig_acc <- matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  
  ABCD_cps <-  matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  l2_cps <-  matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  scan_cps <- matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  
  ABCD_ps <-  matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  scan_ps <- matrix(rep(0, 3*ntrials), nrow = ntrials, ncol = 3)
  
  if(is.null(cp_time)==T){
    cp_time <- 120
  }
  scan_prop <- rep(NA, 3)
  for(i in 1:3){
    for(j in 1:ntrials){
      ABCD_p_val <- sim$ABCD_tests_full[[i]][[j]]$p_val
      ABCD_cp_est <- sim$ABCD_tests_full[[i]][[j]]$tauhat_max
      ABCD_cps[j,i] <- ABCD_cp_est
      ABCD_ps[j,i] <- ABCD_p_val
      ABCD_sig_acc[j,i] <- ifelse(ABCD_p_val < 0.05 & abs(ABCD_cp_est - cp_time) <= 10, 1, 0)
      if(sss == T){
        scan_test <- sim$scan_tests_full[[i]][[j]]
        scan_p_val <- as.numeric(scan_test$P_VALUE[1])
        start_est <- as.numeric(as.matrix(scan_test)[1,7])
        end_est <- as.numeric(as.matrix(scan_test)[1,8])
        ABCD_cps[j,i] <- start_est
        ABCD_ps[j,i] <- end_est
        scan_sig_acc[j,i] <- ifelse(scan_p_val < 0.05 & min(abs(start_est - (cp_time + 1)), abs(end_est - cp_time)) <= 10, 1, 0) # The algorithm is set to detect clusters at level 0.05, 
        # they will be included in this list if their critical value is significant to this level
        # NOTE: as SSS detects clusters instead of an exact change-point, the start time of the cluster found should start at time tau + 1, which is
        # the first post-change observation, hence we we assess abs(start_est - (cp_time + 1)) here ...
      }
      
      if(sim$l2_tests_full[[i]][[j]]$total_num_breaks > 0){ # The algorithm is set to detect changes at level 0.05, listing them in order of significance 
        l2_sig_acc[j,i] <- ifelse(abs(sim$l2_tests_full[[i]][[j]]$nbd_and_stamps_pair[1,2] - cp_time) <= 10, 1, 0)
        l2_cps[j,i] <- sim$l2_tests_full[[i]][[j]]$nbd_and_stamps_pair[1,2]
      } else {
        l2_sig_acc[j,i] <- 0}
    }
  }
  ABCD_prop <- colSums(ABCD_sig_acc)
  if(sss == T){
    scan_prop <- colSums(scan_sig_acc)
  }
  l2_prop <- colSums(l2_sig_acc)
  prop_tbl <- rbind(ABCD_prop, l2_prop, scan_prop)
  rownames(prop_tbl) <- c("ABCD", "2WayMOSUM","SSS")
  colnames(prop_tbl) <- c("5 X 5", "4 X 4", "3 X 3")
  # The object prob_tbl re-creates the sub-table for each data type found in section 4.2 
  return(list(prop_tbl=prop_tbl, ABCD_cps = ABCD_cps, ABCD_ps=ABCD_ps,
              l2_cps=l2_cps, scan_cps = scan_cps, scan_ps = scan_ps))
}