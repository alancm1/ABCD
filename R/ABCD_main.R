#'  Main function implementing Adaptive Block-based Change-point Detection (ABCD)
#'
#' This function implements the change-point method ABCD from the corresponding paper
#' on a given time series, using a range of blocking structures to test for spatially
#' clustered change-points at multiple scales
#'
#' @param y An n by d time series to be analyzed
#' @param P_grid A vector of numbers of blocks associated with each blocking structure
#' @param k The value k used for constructing the k-MSTs within blocks
#' @param B The number of permutations used to generate a permutation p-value
#' @param n0 First observation which a test statistic is computed for
#' @param n1 Last observation which a test statistic is computed for
#' @param perm.p A T/F statement on whether a permutation p-value should be computed
#' @param block_list_of_lists A list of lists each indicating dimension indices
#' for utilizing more complicated blocking structures where blocks may overlap; overrides P_grid if used
#' @return A list containing:
#' \item{tauhat_max}{Final ABCD change-point estimate}
#' \item{test_stat}{Final ABCD test statistic}
#' \item{p_val}{The associated permutation p-value, if computed}
#' \item{tauhat_inds}{A length(P_grid) by 2 matrix of which block obtained the maximum statistic in each blocking structure}
#' \item{tauhat_cols}{A length(P_grid) by 2 matrix of the time of maximum scan statistic within each blocking structure}
#' \item{tauhat_maxes}{A vector of the maximum scan statistic value within each blocking structure}
#' \item{boot_stats}{A vector of B-1 bootstrap statistics}
#' @export
ABCD <-function(y, P_grid = NULL, k = NULL, B = 500,
                n0 = NULL, n1 = NULL, perm.p = F,
                block_list_of_lists = NULL){


  n <- nrow(as.matrix(y)) # "n" is the total number of nodes (observations).
  d <- ncol(as.matrix(y))# "d" is the dimension of the original series.

  if(is.null(k)==T){
    k <- floor(n/5)
  }

  if(is.null(n0) == T){n0 = floor(0.05*n)} # getting range of detection, if not pre-specified
  if(is.null(n1) == T){n1 = ceiling(0.95*n)}

  # Collecting information on statistics within each blocking structures
  block_stats_list <- list()
  tauhat_inds <- c()
  tauhat_cols <- c()
  tauhat_maxes <- c()

  # Only used if blocks are manually specified (with potential to overlap)
  if(is.null(block_list_of_lists)==F){
    P_grid <- seq(1, by = 1, length.out = length(block_list_of_lists))
  } # setting P_grid to have the same length as the list of lists given for manual blocks,
  # to keep code concise

  # Generating automatic blocking if no specifications are given for P_grid
  if(is.null(P_grid)==T){
    P_grid <- generate_P_vec(d = d)
  }

  maxes_per_col <- matrix(NA, nrow = length(P_grid), ncol = length(n0:n1))

  # Generating a list of matrices to store bootstrap statistics for each blocking structure
  if(perm.p == T){
    boot_mpc_list <- replicate((B-1),  matrix(NA, nrow = length(P_grid), ncol = length(n0:n1)), simplify = FALSE)}

  # Generating data permutations OUTSIDE of the loop, to perform the exact same test, but with permuted data
  perm_stuff <- list()
  if(perm.p == T){
    perm_stuff <- perm_func(B, n)}

  for(m in 1:length(P_grid)){
    if(is.null(block_list_of_lists) == T){
      P <- P_grid[m]
      breaks <- breaks_1d(d, P)} else {
        block_list <- block_list_of_lists[[m]]
        P <- length(block_list)
      }

    if(perm.p == T){
      Z_mats <- replicate(B, array(dim = c(P, n)), simplify = FALSE) # generating B matrices of NAs
    } else {
      Z_mats <- replicate(1, array(dim = c(P, n)), simplify = FALSE)
    }

    for(l in 1:P){ # Looping over each block
      if(is.null(block_list_of_lists) == T){
        start <- breaks[l] # beginning of block
        end <- breaks[l+1]-1 # end of block
        block <- y[,start:end]
      } else {
        block <- y[,block_list[[l]]]
      }

      # creates a k-MST, 10 by default.
     # if(is.null(graph.type) == T){
     #   E <- kmst(block, k=k) # We need to generate E for each block
    #  }  else {
     #   if(graph.type == "kmst"){
          E <- kmst(block, k=k) # We need to generate E for each block
     #   }
      # if(graph.type == "knn"){  # code for using a knn graph, ultimately unused
      #    if(is.null(knn.type) == F){
      #      E <- to_edges(get.knn(block, k = k, algorithm = knn.type)$nn.index)
      #    } else {
      #      E <- to_edges(get.knn(block, k = k)$nn.index) # We need to generate E for each block
      #    }
      #  }
      #}
      Ebynode <- Ebynode_gen(n, E)

      Ebnstars <- list()
      Ebnstars[[1]] <- Ebynode # the first entry is the actual node info

      if(perm.p == T){
        #     perm_stuff <- perm_func(B, n)
        permmatches <- perm_stuff[[1]]
        perms <- perm_stuff[[2]]
        for(b in 2:B){ # the other B-1 statistics are null generated node info, from the permutations created above
          Ebnstar =  vector("list", n)
          for(i in 1:n){
            oldlinks = Ebynode[[permmatches[[b]][i]]]
            Ebnstar[[i]] = perms[[b]][oldlinks]
          }
          Ebnstars[[b]] <- Ebnstar}}

      if(perm.p == F){
        Z <- stat_row(Ebynode = Ebynode, edge_info = Ebnstars[[1]], n = n, n0=n0, n1=n1)
        Z_mats[[1]][l,] <- Z } else {
          for(b in 1:B){
            Z <- stat_row(Ebynode = Ebynode, edge_info = Ebnstars[[b]], n = n, n0=n0, n1=n1)
            Z_mats[[b]][l,] <- Z
          }
        }
    } # end of 1:P loop

    tauhat_ind = which(Z_mats[[1]] == max(Z_mats[[1]][,n0:n1]), arr.ind = T)[1] # Index of block with largest statistic
    tauhat_col = which(Z_mats[[1]] == max(Z_mats[[1]][,n0:n1]), arr.ind = T)[2] # The estimate column (time/location)
    tauhat_max <- max(Z_mats[[1]][,n0:n1])
    block_stats <- Z_mats[[1]]

    #Special case if blocking structure consists of a single block
    if(P==1){
      max_per_col = block_stats[,n0:n1]
    } else {
      max_per_col <- apply(block_stats[,n0:n1], 2, max)}

    if(perm.p == T){
      for(b in 2:B){
        if(P==1){
          permZ_mat = Z_mats[[b]][n0:n1]
          boot_mpc_list[[(b -1)]][m,] <- permZ_mat
        } else {
          permZ_mat = Z_mats[[b]][,n0:n1]
          boot_mpc_list[[(b -1)]][m,] <- apply(permZ_mat, 2, max)}
      }
    }

    tauhat_cols <- c(tauhat_cols, tauhat_col)
    tauhat_inds <- c(tauhat_inds, tauhat_ind)
    tauhat_maxes <- c(tauhat_max, tauhat_maxes)
    block_stats_list[[m]] <- block_stats
    maxes_per_col[m, ] <- max_per_col
  }
  tauhat_max = (which.max(colMeans(maxes_per_col)) + (n0-1))
  test_stat = max(colMeans(maxes_per_col))

  # Generating p-value using bootstrap statistics
  p_val <- NULL
  boot_stats <- c()
  if(perm.p == T){
    for(b in 2:B){
      boot_stats <- c(boot_stats, max(colMeans(boot_mpc_list[[(b-1)]])))
    }

    p_val <- length(which(c(test_stat, boot_stats) >= test_stat)) / B
  }

  # Generating clean output for information from each blocking structure...
  tauhat_inds <- as.matrix(tauhat_inds)
  colnames(tauhat_inds) <- "Block w/ Max Stat"
  if(is.null(block_list_of_lists) == T){
    rownames(tauhat_inds) <- paste(rep("P = ", length(P_grid)),  as.character(P_grid), sep="")
  } else {
    rownames(tauhat_inds) <- paste(rep("List ", length(P_grid)),  as.character(P_grid), sep="")
  }

  tauhat_cols <- as.matrix(tauhat_cols)
  colnames(tauhat_cols) <- "Time of Max Stat"
  if(is.null(block_list_of_lists) == T){
    rownames(tauhat_cols) <- paste(rep("P = ", length(P_grid)),  as.character(P_grid), sep="")
  } else {
    rownames(tauhat_cols) <- paste(rep("List ", length(P_grid)),  as.character(P_grid), sep="")
  }

  tauhat_maxes <- as.matrix(tauhat_maxes)
  colnames(tauhat_maxes) <- "Max Stat"
  if(is.null(block_list_of_lists) == T){
    rownames(tauhat_maxes) <- paste(rep("P = ", length(P_grid)),  as.character(P_grid), sep="")
  } else {
    rownames(tauhat_maxes) <- paste(rep("List ", length(P_grid)),  as.character(P_grid), sep="")
  }

  return(list(tauhat_max = tauhat_max, test_stat = test_stat, p_val = p_val,
       tauhat_cols = tauhat_cols, tauhat_inds = tauhat_inds,
       tauhat_maxes = tauhat_maxes, boot_stats = boot_stats))
}
