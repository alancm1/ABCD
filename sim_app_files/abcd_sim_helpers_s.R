# ABCD Simulation Study helper functions

# This file includes functions which help to perform all of the simulations comparing ABCD, DC, Inspect and SCAN automatically.
# The file "simulation_testing_s.Rmd" contains code to generate the exact simulations found in the paper, in sections 4.1 and 4.2

# Requires the following packages to perform simulations (factorcpt is found on the following Github https://github.com/markov10000/factorcpt)
library(MASS)
library(InspectChangepoint)
library(factorcpt)
library(ggplot2)

# Function to generate covariance matrix with local structure (i.e. Cov(h_i, h_j) = \rho^(|i-j|))
local_cov_mat <- function(rho, d){ #rho is a value between 0 and 1, d is the dimension of the time series
construction_vec <- c()
for(i in 1:(d-1)){
  construction_vec <- c(construction_vec, rho^(d - (i)))
}
construction_vec[d] <- 1
construction_vec <-c(construction_vec, rev(construction_vec[1:(d-1)]))

cov_mat <- matrix(NA, nrow = d, ncol = d)
for(i in 1:d){
  start <- ((2*d) - ((d-1) + i))
  end <- start + d - 1
  cov_mat[i, ] <- construction_vec[start:end]
}
return(cov_mat)
}

# Levels of compactness of changed region
p_cs <- c(1/4, 1/2, 1)
n_p_cs <- length(p_cs)

# Functions for manual implementation of SCAN (Enikeeva and Harchaoui (2019))
# Implemented as described in the paper

# X is a d by n matrix
# t is the time point of interest
# n is the total number of observations
Z_n_t <- function(X, t, n){
  return(sqrt(t*(n-t) / n) * (rowMeans(X[,1:t]) - rowMeans(X[,(t+1):n])) )
}

L_lin <- function(X, d, n){
  stats <- c()
  for(i in 2:(n-2)){
    Z <- Z_n_t(X = X, t= i, n = n)
    L2_Z <- sqrt(sum(Z^2))
    stats <- c(stats, (L2_Z^2 - d) /sqrt(2*d))
  }
  return(list(stat = max(stats), tauhat = which.max(stats) + 1 ))
}


T_p_n <- function(d, n, p, alpha_s){
  return(sqrt(2/p)*log(choose(d, p)*(2*n*p^2) / alpha_s) + sqrt(2*log(choose(d, p)*(2*n*p^2) / alpha_s)))
}


L_scan <- function(X, d, n, alpha_s, maxp = NULL){
  Z_s <- matrix(NA, nrow = d, ncol = n)
  
  # storing CUSUM vectors for each time point
  for(i in 2:(n-2)){
    Z_s[, i]<- Z_n_t(X = X, t= i, n = n)
  }
  
  stat_vec <- c()
  tauhats <- c()
  for(i in 1:maxp){
    p <- i
    weight <- T_p_n(d, n, p, alpha_s)
    t_stats <- c()
    for(j in 2:(n-2)){
      t_stats <- c(t_stats, (sum(sort((Z_s[,j])^2, decreasing= T)[1:p]) - p)/(sqrt(2*p)))
    }
    stat_vec <- c(stat_vec, max(t_stats) / weight)
    tauhats <- tauhats <- c(tauhats, which.max(t_stats)+ 1)
  }
  return(list(stat = max(stat_vec), tauhat = tauhats[which.max(stat_vec)], best_p = which.max(stat_vec)))
}


lin_cutoff <- function(d, n, alpha_l){
  e <- sqrt(2*log(d) / d)
  a_e <- alpha_l*(log(1 + e) / log(n))
  return( sqrt(log(d)) + (1 + e)*sqrt(2*log(1 / a_e)) + (1 + e)*sqrt(2 / d)*log(1 / a_e) )
}


# The overall SCAN test, with exact quantiles.  maxp generally is never much above 10, we set it's maximum value to be 100 for our simulations, as at most 100 dimensions change.
scan_test <- function(X, d , n, alpha_l, alpha_s, maxp){
  
  lin_stat_init <- L_lin(X, d, n)
  lin_tau <- lin_stat_init$tauhat
  lin_stat_norm <- lin_stat_init$stat / lin_cutoff(d, n, alpha_l)
  
  scan_stat_init <- L_scan(X, d, n, alpha_s, maxp)
  scan_tau <- scan_stat_init$tauhat
  scan_stat_norm <- scan_stat_init$stat
  
  final.stat <- max(scan_stat_norm, lin_stat_norm)
  tau_vec <- as.numeric(c(lin_tau, scan_tau))
  
  stat_vec <- as.numeric(c(lin_stat_norm, scan_stat_norm))
  final.tau <- tau_vec[which(stat_vec == final.stat)]
  return(list(stat = final.stat, tauhat = final.tau, best_p = scan_stat_init$best_p))}

# Function for Monte Carlo simulation, to approximate Type 1 error thresholds for varying alpha, for EACH method
# N is number of null simulations
MC_null_stats <- function(data_type = c("Normal","Normal_cov", "Lognormal", "t3"), N = 20000){
  set.seed(1234)
  if(data_type == "Normal_cov"){
    cov_mat <- local_cov_mat(0.6, 1000)
    chol_cov <- chol(cov_mat) # setting up covariance matrix and corresponding Cholesky decomp to be used 
  }
  ABCD_null_stats <- c()
  ins_null_stats <- c()
  dc_null_stats <- c()
  scan_null_stats <- c()
  for(i in 1:N){
    #print(i)
    if(data_type == "Normal"){
      y <- matrix(rnorm(200000), ncol=200)}
    if(data_type == "Normal"){
      y <- matrix(rnorm(200000), ncol=200)}
    if(data_type == "Normal_cov"){
      y <- matrix(rnorm(200000), ncol=200) 
      y <- t(t(y) %*% chol_cov)} # using Cholesky decomposition to speed up simulating correlated normal data 
    if(data_type == "Lognormal"){
      y <- matrix(rlnorm(200000), ncol=200)}
    if(data_type == "t3"){
      y <- matrix(rt(200000, df = 3), ncol=200)}
    
  ABCD_test <- ABCD(t(y), P_grid = c(1, 4, 10, 20, 40), k = 40)
  ABCD_null_stats <- c(ABCD_null_stats, ABCD_test$test_stat)
  ins_test <- locate.change(y)
  ins_null_stats <- c(ins_null_stats, ins_test$cusum)
  dc_test <- func_dc(y)
  dc_null_stats <- c(dc_null_stats, max(dc_test$res))
  scan_test <- scan_test(y, 1000 , 200, alpha_l= 0.05, alpha_s= 0.05, maxp= 100)
  scan_null_stats <- c(scan_null_stats, scan_test$stat)
  }

  list_of_null_stats <- list(ABCD_null_stats, ins_null_stats, dc_null_stats, scan_null_stats)
  return(list_of_null_stats)
}
  
# Below is a function generating plots to check the accuracy of the Type 1 
# error control for a given set of Monte Carlo null tests
# N is the number of MC tests, n is the number of new null tests to compare it to
null_plot <- function(list_of_null_stats, N, n, data_type = c("Normal", "Normal_var", "Normal_cov", "Lognormal", "t3")){
  if(data_type == "Normal_cov"){
    cov_mat <- local_cov_mat(0.6, 1000)
    chol_cov <- chol(cov_mat) # setting up covariance matrix and corresponding Cholesky decomp to be used 
  }
  mc_stats <- list_of_null_stats
  alphas <- c(0.001, seq(0.005, 0.1, by = 0.005))
n_quants <- length(alphas)
tiles <- rep(N, length(alphas)) - floor(N*alphas)

ABCD_new_null_stats <- c()
ins_new_null_stats <- c()
dc_new_null_stats <- c()
scan_new_null_stats <- c()
for(i in 1:n){
  #print(i)
  #Generating the appropriate data type
  if(data_type == "Normal"){
    y <- matrix(rnorm(200000), ncol=200)
    title_start = "Multivariate Gaussian"}
  if(data_type == "Normal_var"){
    y <- matrix(rnorm(200000), ncol=200)
    title_start = "Multivariate Gaussian"}
  if(data_type == "Normal_cov"){
    y <- matrix(rnorm(200000), ncol=200) 
    y <- t(t(y) %*% chol_cov)
    title_start = "Spatially Dependent Gaussian"} # using Cholesky decomposition to speed up simulating correlated normal data 
  if(data_type == "Lognormal"){
    y <- matrix(rlnorm(200000), ncol=200)
    title_start = "Multivariate Log-normal"}
  if(data_type == "t3"){
    y <- matrix(rt(200000, df = 3), ncol=200)
    title_start = "Multivariate t, df = 3"}
  
  # Running a test of each method on null data
  ABCD_test <- ABCD(t(y), P_grid = c(1, 4, 10, 20, 40), k = 40)
  ABCD_new_null_stats <- c(ABCD_new_null_stats, ABCD_test$test_stat)
  ins_test <- locate.change(y)
  ins_new_null_stats <- c(ins_new_null_stats, ins_test$cusum)
  dc_test <- func_dc(y)
  dc_new_null_stats <- c(dc_new_null_stats, max(dc_test$res))
  scan_test <- scan_test(y, 1000 , 200, alpha_l= 0.05, alpha_s= 0.05, maxp= 100)
  scan_new_null_stats <- c(scan_new_null_stats, scan_test$stat)
}
# putting new statistics into a list
      new_null_stats <- list(ABCD_new_null_stats, ins_new_null_stats, dc_new_null_stats, scan_new_null_stats)
      
# Generating a data frame of proportion of Type 1 error (comparing to the nominal level) for each set of new null statistics
null_dfs <- list()
method_names <- c("ABCD", "Inspect", "Double CUSUM", "SCAN")
for(i in 1:length(mc_stats)){
  method_name <-method_names[i]
  sorted_mc <- sort(as.numeric(mc_stats[[i]]))
  new_nulls <- new_null_stats[[i]]
  #print(new_nulls)
  # print(sorted_mc)  # generating indices of empirical quantiles 
  mc_cutoffs <- sorted_mc[tiles]
  props <- c()
  for(j in 1:length(mc_cutoffs)){
    prop <- length(which(new_nulls >= mc_cutoffs[j]))
    props <- c(props, prop)
  }
  #print(props)
  null_df <- cbind(props, alphas, rep(method_name, length(props)))
  #print(null_df)
  null_dfs[[i]] <- null_df
}

# Row binding each data frame to create a full data frame to plot
abcd_df <- null_dfs[[1]]
ins_df <- null_dfs[[2]]
dc_df <- null_dfs[[3]]
scan_df <- null_dfs[[4]]

full_null_df <- as.data.frame(rbind(null_dfs[[1]], null_dfs[[2]], null_dfs[[3]],null_dfs[[4]]))
colnames(full_null_df) <- c("prop", "alpha", "test_type")
full_null_df$prop <- as.numeric(full_null_df$prop)
full_null_df$alpha <- as.numeric(full_null_df$alpha)
full_null_df$prop <- full_null_df$prop/n
#head(full_null_df)

null_plot_title <-paste(title_start, "FPR under no change, 1000 trials", sep = " ")
# Final plot comparing nominal to actual Type 1 error for each method
null_t1plot <-ggplot(full_null_df, mapping = aes(x = alpha, y = prop, group = test_type)) + 
  geom_line(aes(color = test_type, linetype = test_type)) + 
  geom_point(color="black", size= 0.15, alpha = 0.5) +
  labs(title = null_plot_title, x = "Alpha", y = "False Positive Rate") + 
  geom_abline(intercept = 0, slope = 1, color="black") +theme_bw() +
  theme(panel.spacing = unit(1.3, "lines"), text=element_text(size=10.7,  family="serif"))

return(null_t1plot)
}

# Magnitudes of change for each scenario are given below
# IID Normal (Mean): 0.3= sqrt(2.25/25) (25 dims change), sqrt(2.25/50) (50 dims change), 0.15 = sqrt(2.25/100) (100 dims change), slightly altered so that 2.25 is the L2 distance for each change!!
# IID Normal (Variance): sd = 1.25 (25 dims change), sd = 1.175 (50 dims change) , sd = 1.125 (100 dims change)
# IID Lognormal: sqrt(3/25) (25 dims change), sqrt(3/50) (50 dims change) , sqrt(3/100) (100 dims change), this is a change in the mean of the normal, resulting in a distributional change for the lognormal data
# IID t3: sqrt(5/25) (25 dims change), sqrt(5/50) (50 dims change) ,sqrt(5/100) (100 dims change)!!
# Normal with local covariance: sqrt(4/25) (25 dims change), sqrt(4/50) (50 dims change) , sqrt(4/100)  (100 dims change), slightly altered so that 4 is the L2 distance for each change!!

# Function to generate data for each simulation setting; there are 3 such settings per data type/change time, we have a total 
# of 3*6 = 18 scenarios

# All changes occur at t = 120, changes detected between 110-130 are considered accurate
# Each data frame is 1000 rows, 200 columns
data_gen <- function(ch_dims = c(25, 50, 100), d = 1000, n = 200, data_type = c("Normal", "Normal_var", "Normal_cov", "Lognormal", "t3")){
  # a method to find which magnitude of change to use depending on the input
  n_dims_changed <- c(25, 50, 100)
  delta_index <- which(ch_dims == n_dims_changed)
    
  # Magnitudes of change in mean (or variance) for the different settings
  norm_mean_deltas <- c(sqrt(2.25/25), sqrt(2.25/50), sqrt(2.25/100))
  norm_var_newsds <- c(1.25, 1.175, 1.125)
  lognorm_deltas <- c(sqrt(3/25), sqrt(3/50), sqrt(3/100))
  t3_mean_deltas <- c(sqrt(5/25), sqrt(5/50), sqrt(5/100))
  norm_cov_deltas <- c(sqrt(4/25), sqrt(4/50), sqrt(4/100))
  
  p_cs <- c(1/4, 1/2, 1)
  n_p_cs <- length(p_cs)
  list_of_lists <- list()
  
# generating
  if(data_type == "Normal_cov"){
    cov_mat <- local_cov_mat(0.6, d)
  }
  
for(j in 1:n_p_cs){
  list_of_dfs <- list()
  p_c <- p_cs[j]
  print(p_c)
  
for(i in 1:100){
  # setting the size of the changed region, and sampling dimensions within this region which will change
  cr_size <- ch_dims/p_c
  changed <- sample(cr_size, size = ch_dims) #}
  
  if(data_type == "Normal" | data_type == "Normal_var"){
    data <- matrix(data = rnorm(d*n), ncol = n)
    delta <- norm_mean_deltas[delta_index]
    if(data_type == "Normal"){
    for(k in changed){
      data[k,121:200] <- rnorm(80, mean = delta)
    }
  } else {
    newsd <- norm_var_newsds[delta_index]
    for(k in changed){
      data[k,121:200] <- rnorm(80, sd = newsd)
      }
    }
  }

  if(data_type == "Normal_cov"){
    data <- t(as.matrix(mvrnorm(n = 200, mu = rep(0, 1000), Sigma = cov_mat)))
    delta <- norm_cov_deltas[delta_index]
    for(k in changed){
      data[k,121:200] <- data[k, 121:200] + delta
    }
  }
  
  if(data_type == "Lognormal"){
    data <- matrix(data = rlnorm(d*n), ncol = n)
    delta <- lognorm_deltas[delta_index]
    for(k in changed){
      data[k,121:200] <- rlnorm(80, meanlog = delta)
    }
  }
  
  if(data_type == "t3"){
    data <- matrix(data = rt(d*n, df = 3), ncol = n)
    delta <- t3_mean_deltas[delta_index]
    for(k in changed){
      data[k,121:200] <- data[k,121:200] + delta
    }
  }
  
  list_of_dfs[[i]] <- data
}
  list_of_lists[[j]] <- list_of_dfs
}
  
return(list_of_lists)
}
    
    
  
single_sim <-function(sim_dfs, method = c("ABCD", "Inspect", "DC", "SCAN"), change_type = c("mean", "variance", "both")){
  stat_mat <- matrix(NA,ncol = n_p_cs, nrow = 100)
  chpt_mat <- matrix(NA,ncol = n_p_cs, nrow = 100)
  for(j in 1:length(p_cs)){
    data_list <- sim_dfs[[j]]
    #print(j)
    for(i in 1:100){
      y <- data_list[[i]]
      
  if(method == "ABCD"){
    test <- ABCD(t(y), P_grid = c(1, 4, 10, 20, 40), k = 40)
    stat_mat[i,j] <- test$test_stat
    chpt_mat[i,j] <- test$tauhat_max
    #print(i)
  }
      
  if(method == "Inspect"){
      test <- InspectChangepoint::locate.change(y)
        stat_mat[i,j] <- test$cusum
        chpt_mat[i,j] <- test$changepoint
      }
      
  if(method == "DC"){
      test <- factorcpt::func_dc(y)
      stat_mat[i,j] <-  max(test$res)
      chpt_mat[i,j] <-  which.max(test$res)
      }
      
  if(method == "SCAN"){
        test <-  scan_test(y, d = 1000, n = 200, alpha_l = 0.05, alpha_s = 0.05, maxp = 100)
        stat_mat[i,j] <- test$stat
        chpt_mat[i,j] <- test$tauhat
      }
    
     }
  }
  return(list(stat_mat = stat_mat, chpt_mat = chpt_mat))
}

single_sim_setting <- function(sim_dfs, change_type = c("mean", "variance", "both")){
  methods <- c("ABCD", "Inspect", "DC", "SCAN")
  list_of_results <- list()
  for(i in 1:4){
    list_of_results[[i]] <- single_sim(sim_dfs = sim_dfs, method = methods[i], change_type = change_type)
  }
  
  return(list_of_results = list_of_results)
}
  
# Function to generate a facet plot of a simulation for given a data type and number of dimensions changing 
sim_plots <- function(list_of_results, mc_stats, N, plot_title){ 
  # list of results should be results taken from a call of single_sim_setting
  # The list of Monte Carlo null stats should be taken from a call from MC_null_stats.
  # N is the number of MC simulations used
  methods <- c("ABCD", "Inspect", "DC", "SCAN")
  alphas <- c(0.001, seq(0.005, 0.1, by = 0.005))
  n_quants <- length(alphas)
  #prop_dfs <- list()
  full_prop_df <- matrix(NA, ncol = 4, nrow = 4*n_p_cs*n_quants)
  
  for(i in 1:4){ # iterating through the 4 methods
    start_row_ind <- (i-1)*n_p_cs*n_quants + 1
    end_row_ind <- i*n_p_cs*n_quants
    method <- methods[i]
     stats <- as.matrix(list_of_results[[i]]$stat_mat)
     preds <- as.matrix(list_of_results[[i]]$chpt_mat)
    for(j in 1:100){
      for(k in 1:n_p_cs){
        if(preds[j,k] < 110 | preds[j,k] > 130){ # Essentially filtering out non-accurate observations
          stats[j,k] <- 0
        }
      }
    }
     #print(preds)
     #print(stats)
     
     # Creating a partial data frame of proportions of significant accurate change-points for a given method from matrices
     # of change-point estimates and corresponding statistics
     sorted_mc <- sort(as.numeric(mc_stats[[i]]))
    # print(sorted_mc)
     tiles <- rep(N, length(alphas)) - floor(N*alphas)  # generating indices of empirical quantiles 
     mc_cutoffs <- sorted_mc[tiles] # finding empirical thresholds
    # print(mc_cutoffs)
     prop_mat <- matrix(data = NA, nrow = length(mc_cutoffs), ncol = n_p_cs)
     for(i in 1:length(mc_cutoffs)){
       props <- apply(stats, MAR = 2, FUN = function(x){return(sum(ifelse(x > mc_cutoffs[i], 1, 0)))}) 
       # finds proportion of significant accurate change-points for a given empirical threshold
       prop_mat[i,] <- as.numeric(props)
     }
     #print(prop_mat)
     colnames(prop_mat) <- c("1/4",
                             "1/2",
                             "1")
     prop_df <- as.data.frame(cbind(as.numeric(c(prop_mat[,1], prop_mat[,2], prop_mat[,3] )),
                                    c(rep("1/4", n_quants),
                                      rep("1/2", n_quants),
                                      rep("1", n_quants)),
                                    as.numeric(rep(alphas, n_p_cs)),
                                    rep(method, n_p_cs*n_quants)))
     #prop_dfs[[i]] <- prop_df
     #print(dim(prop_df))
    # print(prop_df)
     #print(c(start_row_ind, end_row_ind))
     full_prop_df[start_row_ind:end_row_ind,] <- as.matrix(prop_df)
    
     #print(full_prop_df)
  }
  
  # Combining proportion matrices for each method, prepping it to be plotted
  colnames(full_prop_df) <- c("sig_acc_rate", "cr_prob", "alpha", "test_type")
  full_prop_df <- as.data.frame(full_prop_df)
  full_prop_df$cr_prob <- factor(full_prop_df$cr_prob,
                         levels =  c("1/4",
                                     "1/2",
                                     "1"))
  full_prop_df$sig_acc_rate <- as.numeric(full_prop_df$sig_acc_rate) / 100
  
  # Labels for Facet Plots
  my_labeller = as_labeller(
    c("1/4" = bquote('"p"[italic(C)]*" = 1/4"'), 
      "1/2" = bquote('"p"[italic(C)]*" = 1/2"'),
      "1" = bquote('"p"[italic(C)]*" = 1"')), 
    default = label_parsed)
  
  plot <- ggplot(full_prop_df, mapping = aes(x = as.numeric(alpha), y = as.numeric(sig_acc_rate), group = test_type)) + 
    geom_line(aes(color = test_type, linetype = test_type, show.legend = FALSE), size = 3/4) +
    geom_point(color="black", size= 0.2, alpha = 0.5) + 
    labs(title = plot_title) + 
    scale_x_continuous(limits = c(0, 0.1)) + scale_y_continuous(limits = c(0, 1)) +
    ylab("Detection Accuracy") + xlab(bquote("Nominal Type 1 Error Rate")) +
    facet_wrap(~ cr_prob,  labeller = my_labeller) + 
    labs(color = "Test", linetype = "Test") +theme_bw() + theme(panel.spacing = unit(1.3, "lines"),
                                             text=element_text(size=11.5,  family="serif"))
  return(plot)
}


# Function to return results for change-point simulation for a data type, including all simulation settings for each method, with plots
sim_full_data_type <- function(data_type = c("Normal", "Normal_var", "Normal_cov", "Lognormal", "t3"),
                     N = 10000, plot_titles = NULL, MC_nulls = NULL){
    ABCD_stype <- "m"

  if(is.null(MC_nulls) == T){
  MC_nulls <- MC_null_stats(data_type = data_type, N = N)}
    null_title <- paste(data_type, "Type 1 Error Actual vs. Nominal", sep = " ")
    null_t1plot <- null_plot(list_of_null_stats = MC_nulls, N = N, n = 1000, data_type = data_type)
    print(null_t1plot)
  ch_dims <- c(25, 50, 100) 
  
  # Setting type of change based on simulation data
  if(data_type == "Normal_var"){
    change_type <- "variance"
  } else if (data_type == "Lognormal"){
    change_type <- "both"
  } else {
    change_type <- "mean"
  }
  
    ABCD_stype <- "m"
  
  list_of_plots <- list()
  results_list_of_lists <- list()
  # iterate over the differing number of changed dimensions
  for(i in 1:3){
    sim_dfs <- data_gen(ch_dim = ch_dims[i], data_type = data_type)
    results_list <- single_sim_setting(sim_dfs, change_type = change_type)
    results_list_of_lists[[i]] <- results_list
    list_of_plots[[i]] <- sim_plots(list_of_results = results_list, mc_stats = MC_nulls, N = N, plot_title = plot_titles[i])
  }
return(list(plot25 = list_of_plots[[1]], plot50 = list_of_plots[[2]], plot100 = list_of_plots[[3]],
            results_list_of_lists = results_list_of_lists, MC_nulls = MC_nulls, null_t1plot = null_t1plot))
}


entire_sim <- function(N = 10000, use_old_mc_trials = F){
  null_title_starts <- c("IID Gaussian:", "IID Gaussian:", "Spatially Dependent Gaussian:", "IID Log-normal:", "IID t, df = 3:")
  data_types <- c("Normal", "Normal_var", "Normal_cov", "Lognormal", "t3")
  ABCD_stypes <- c("m", "m", "m", "m", "m") # Using only the max-type statistic
  plot_titles_norm_mean <-  c("Multivariate Gaussian: 25 of 1000 dimensions change in mean",
                              "Multivariate Gaussian: 50 of 1000 dimensions change in mean",
                              "Multivariate Gaussian: 100 of 1000 dimensions change in mean")
  plot_titles_norm_var <-  c("Multivariate Gaussian: 25 of 1000 dimensions change in variance",
                             "Multivariate Gaussian: 50 of 1000 dimensions change in variance",
                             "Multivariate Gaussian: 100 of 1000 dimensions change in variance")
  plot_titles_cov_norm <- c("Spatially Dependent Gaussian: 25 of 1000 dimensions change in mean",
                            "Spatially Dependent Gaussian: 50 of 1000 dimensions change in mean",
                            "Spatially Dependent Gaussian: 100 of 1000 dimensions change in mean")
  plot_titles_lognorm <- c("Multivariate Log-normal: 25 of 1000 dimensions change in distribution",
                           "Multivariate Log-normal: 50 of 1000 dimensions change in distribution",
                           "Multivariate Log-normal: 100 of 1000 dimensions change in distribution")
  plot_titles_t3 <- c("Multivariate t, df = 3: 25 of 1000 dimensions change in mean",
                      "Multivariate t, df = 3: 50 of 1000 dimensions change in mean",
                      "Multivariate t, df = 3: 100 of 1000 dimensions change in mean")
  plot_titles_list <- list(plot_titles_norm_mean, plot_titles_norm_var, plot_titles_cov_norm,
                           plot_titles_lognorm, plot_titles_t3)
  
  data_type_sim_list <- list()
  
 # methods <- c("ABCD", "Inspect", "DC", "SCAN") order of methods tested
  if(use_old_mc_trials==T){
    print("Using pre-run Monte Carlo Trials...")
    # Intended for developer use only, below code reads in pre-run Monte Carlo trials to speed up results
    ABCD_norm_m_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/ABCD_norm_MC_stats.csv")[,2]
    inspect_norm_nulls <-  read.csv("/Users/alanmoore/Documents/final_MC_stats/inspect_norm_MC_stats.csv")[,2]
    dc_norm_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/DC_norm_MC_stats.csv")[,2]
    scan_norm_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/SCAN_norm_MC_stats.csv")[,2]

    ABCD_cov_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/ABCD_covnorm_MC_stats.csv")[,2]
    inspect_cov_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/inspect_covnorm_MC_stats.csv")[,2]
    dc_cov_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/DC_covnorm_MC_stats.csv")[,2]
    scan_cov_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/SCAN_covnorm_MC_stats.csv")[,2]
    
    ABCD_log_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/ABCD_lognorm_MC_stats.csv")[,2]
    inspect_log_nulls <-  read.csv("/Users/alanmoore/Documents/final_MC_stats/inspect_lognorm_MC_stats.csv")[,2]
    dc_log_nulls <-  read.csv("/Users/alanmoore/Documents/final_MC_stats/DC_lognorm_MC_stats.csv")[,2]
    scan_log_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/SCAN_lognorm_MC_stats.csv")[,2]
    
    ABCD_t3_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/ABCD_t3_MC_stats.csv")[,2]
    inspect_t3_nulls <-  read.csv("/Users/alanmoore/Documents/final_MC_stats/inspect_t3_MC_stats.csv")[,2]
    dc_t3_nulls <-  read.csv("/Users/alanmoore/Documents/final_MC_stats/DC_t3_MC_stats.csv")[,2]
    scan_t3_nulls <- read.csv("/Users/alanmoore/Documents/final_MC_stats/SCAN_t3_MC_stats.csv")[,2]
    
    # Aggregating all Monte Carlo null trials into a list of lists
    list_of_list_of_MC_nulls <- list()
    
    iidnorm_MC_nulls <- list(ABCD_norm_m_nulls, inspect_norm_nulls, dc_norm_nulls, scan_norm_nulls)
    covnorm_MC_nulls <- list(ABCD_cov_nulls, inspect_cov_nulls, dc_cov_nulls, scan_cov_nulls)
    iidlognorm_MC_nulls <- list(ABCD_log_nulls, inspect_log_nulls, dc_log_nulls, scan_log_nulls)
    iidt3_MC_nulls <- list(ABCD_t3_nulls, inspect_t3_nulls, dc_t3_nulls, scan_t3_nulls)
    list_of_list_of_MC_nulls <- list(iidnorm_MC_nulls, iidnorm_MC_nulls, covnorm_MC_nulls, iidlognorm_MC_nulls, iidt3_MC_nulls)
    
    # Giving a list of Monte Carlo trials for each data type
    for(m in 1:5){
      data_type <- data_types[m]
      plot_titles <- plot_titles_list[[m]]
      ABCD_stype <- ABCD_stypes[m]
      data_type_sim <- sim_full_data_type(data_type = data_type, 
                                            N = N, plot_titles = plot_titles, MC_nulls = list_of_list_of_MC_nulls[[m]])
      data_type_sim_list[[m]] <- data_type_sim
    }
  } else { # If Monte Carlo trials for Type 1 error control need to be created from scratch
  # iterating through data type
  for(m in 1:5){
    data_type <- data_types[m]
    plot_titles <- plot_titles_list[[m]]
    ABCD_stype <- ABCD_stypes[m]
    if(m == 2){ # No need to re-run Monte Carlo simulations for the same data type under the null hypothesis (IID normal)
      data_type_sim <- sim_full_data_type(data_type = data_type, N = N,
                                  plot_titles = plot_titles, MC_nulls = data_type_sim_list[[1]]$MC_nulls)
    } else {
    data_type_sim <- sim_full_data_type(data_type = data_type, 
                                        N = N, plot_titles = plot_titles, MC_nulls = NULL)}
    data_type_sim_list[[m]] <- data_type_sim
    
  }
}
  return(simulation_results = data_type_sim_list)
  #plot25, plot50 and plot100 from each entry in the simulation_results list are the 3 corresponding plots for each data type
}
