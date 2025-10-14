# Image-focused version of ABCD
# y input is now array, with the 1st and 2nd dimensions being x and y, and the 3rd being time

#'  Main function implementing Adaptive Block-based Change-point Detection (ABCD)
#'
#' This function implements the change-point method ABCD from the corresponding paper
#' on a given time series, using a range of blocking structures to test for spatially
#' clustered change-points at multiple scales
#'
#' @param y An d_1 by d_2 by n array to be analyzed
#' @param P_1_grid A vector of numbers of row divisions associated with each blocking structure
#' @param P_2_grid A vector of numbers of column divisions associated with each blocking structure
#' @param k The value k used for constructing the k-MSTs within blocks
#' @param B The number of permutations used to generate a permutation p-value
#' @param n0 First observation which a test statistic is computed for
#' @param n1 Last observation which a test statistic is computed for
#' @param perm.p A T/F statement on whether a permutation p-value should be computed
#' for utilizing more complicated blocking structures where blocks may overlap; overrides P_grid if used
#' @return A list containing:
#' \item{tauhat_max}{Final ABCD change-point estimate}
#' \item{test_stat}{Final ABCD test statistic}
#' \item{p_val}{The associated permutation p-value, if computed}
#' \item{tauhat_inds}{A length(P_grid) by 2 matrix of which block obtained the maximum statistic in each blocking structure}
#' \item{tauhat_cols}{A length(P_grid) by 2 matrix of the time of maximum scan statistic within each blocking structure}
#' \item{boot_stats}{A vector of B-1 bootstrap statistics}
#' \item{max_stats_at_change}{A vector of the maximum two-sample test statistic at tauhat_max for each blocking structure. With max_blocks_at_change, can be used to find a locational estimate of a potential change-point detected at tauhat_max.}
#' \item{max_blocks_at_change}{A vector of the block with the maximum two-sample test statistic at tauhat_max for each blocking structure}
#' @export
#'
ABCD_img <-function(y, P_1_grid = NULL, P_2_grid = NULL, k = NULL, B = 500,
                    n0 = NULL, n1 = NULL, perm.p = F){

  if(length(P_1_grid) != length(P_2_grid)){
    print("Vertical/Horizontal Grid Vectors must be of equal length.")
    return()
  }

  flat_blocks_list <- list()
  blocks <- list()

  n <- dim(y)[3] # "n" is the total number of images (observations).
  d <- dim(y)[1]*dim(y)[2] # "d" is the dimension of the original series.
  side_length_1 <- dim(y)[1]
  side_length_2 <- dim(y)[2]

  if(is.null(k)==T){
    k <- floor(n/5)
  }

  if(is.null(n0) == T){n0 = ceiling(0.05*n)} # getting range of detection, if not pre-specified
  if(is.null(n1) == T){n1 = floor(0.95*n)}

  # Generating automatic blocking if no specifications are given for P_grid
  if(is.null(P_1_grid)==T){
    P_1_grid <- generate_P_vec(d = side_length_1, target_min_block_size = 5)
  }

  if(is.null(P_2_grid)==T){
    P_2_grid <- generate_P_vec(d = side_length_2, target_min_block_size = 5)
  }

  block_stats_list <- list()
  tauhat_inds <- matrix(NA, nrow = length(P_1_grid), ncol = 2)
  tauhat_cols <- c()
  maxes_per_col <- matrix(NA, nrow = length(P_1_grid), ncol = length(n0:n1))

  if(perm.p == T){
    boot_mpc_list <- replicate((B-1),  matrix(NA, nrow = length(P_1_grid), ncol = length(n0:n1)), simplify = FALSE)}
  if(perm.p == T){
    perm_stuff <- perm_func(B, n)}

  for(q in 1:length(P_1_grid)){
    flat_blocks <- list()
    P_1 <- P_1_grid[q]
    P_2 <- P_2_grid[q]
    break_list <- breaks_img(P_1 = P_1, P_2 = P_2, side_length_1 = side_length_1, side_length_2 = side_length_2)
    breaks_1 <- break_list[[1]]
    breaks_2 <- break_list[[2]]

    if(perm.p == T){
      Z_arrays <- replicate(B, array(dim = c(P_1,P_2,n)), simplify = FALSE) #generating an array of NAs
    } else {
      Z_arrays <- replicate(1, array(dim = c(P_1,P_2,n)), simplify = FALSE) #generating an array of NAs
    }

    for(l in 1:P_1){ # Looping over each block
      start_x <- breaks_1[l] # left x-boundary
      end_x <- breaks_1[l+1]-1 # right x-boundary
      for(m in 1:P_2){
        start_y <- breaks_2[m] # upper y-boundary
        end_y <- breaks_2[m+1]-1 # lower y-boundary
        block <- y[start_x:end_x, start_y:end_y, ]

        # Flattening the block into a matrix of n
        flat_block <- matrix(,nrow = 0, ncol = n)
        for(i in start_x:end_x){
          for(j in start_y:end_y){
            flat_block <- rbind(flat_block, y[i,j, ])
          }
        }

        list_index <- (l-1)*P_2 + m
        flat_blocks[[list_index]] <- flat_block
        blocks[[list_index]] <- block

        E <- kmst(t(flat_block), k=k) # Generating E for each block
        # Code for knn implementation unused...
       # if(is.null(graph.type) == T){
       #   E <- kmst(t(flat_block), k=k, dis.type=dis.type) # We need to generate E for each block
       # }  else {
       #   if(graph.type == "kmst"){
       #     E <- kmst(t(flat_block), k=k, dis.type=dis.type) # We need to generate E for each block
       #   }
       #   if(graph.type == "knn"){
       #     if(is.null(knn.type) == F){
       #       E <- to_edges(get.knn(t(flat_block), k = k, algorithm = knn.type)$nn.index)
       #     } else {
       #       E <- to_edges(get.knn(t(flat_block), k = k)$nn.index) # We need to generate E for each block
       #     }
       #   }
       # }

        Ebynode <- Ebynode_gen(n, E)

        Ebnstars <- list()
        Ebnstars[[1]] <- Ebynode # the first entry is the actual node info, subsequent entries are for permuted observations

        if(perm.p == T){
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
          Z_arrays[[1]][l,m,] <- Z } else {
            for(b in 1:B){
              Z <- stat_row(Ebynode = Ebynode, edge_info = Ebnstars[[b]], n = n, n0=n0, n1=n1)
              Z_arrays[[b]][l,m,] <- Z
            }
          }
      }} # end of l/m loops, looping over each image section

    block_stats <- Z_arrays[[1]]

    if(P_1 == 1 & P_2 == 1){
      max_per_col = block_stats[,,n0:n1]
    } else {
      max_per_col <- apply(block_stats[,,n0:n1], 3, max)}
    tauhat_ind = which(Z_arrays[[1]] == max(Z_arrays[[1]][,,n0:n1]), arr.ind = T)[1:2]  # The estimate index
    tauhat_col = which(Z_arrays[[1]] == max(Z_arrays[[1]][,,n0:n1]), arr.ind = T)[3]    # The estimate column (time/location)

    if(perm.p == T){
      for(b in 2:B){
        if(P_1==1 & P_2 == 1){
          permZ_array = Z_arrays[[b]][,,n0:n1]
          boot_mpc_list[[(b -1)]][m,] <- permZ_array
        } else {
          permZ_array = Z_arrays[[b]][,,n0:n1]
          boot_mpc_list[[(b - 1)]][q,] <- apply(permZ_array, 3, max)}
      }
    }

    tauhat_cols <- c(tauhat_cols, tauhat_col)
    tauhat_inds[q,] <- tauhat_ind
    block_stats_list[[q]] <- block_stats
    maxes_per_col[q, ] <- max_per_col
    flat_blocks_list[[q]] <- flat_blocks
  } # end of q loop

  tauhat_max = (which.max(colMeans(maxes_per_col)) + (n0-1))
  test_stat = max(colMeans(maxes_per_col))

  # Generation of p-value
  p_val <- NULL
  boot_stats <- c()
  if(perm.p == T){
    for(b in 2:B){
      boot_stats <- c(boot_stats, max(colMeans(boot_mpc_list[[(b-1)]])))
    }

    p_val <- length(which(c(test_stat, boot_stats) >= test_stat)) / B
  }

  # Generating clean results for information from tests in each blocking structure
  tauhat_cols <- as.matrix(tauhat_cols)
  colnames(tauhat_cols) <- "Time of Max Stat"
  p1lab <- paste(rep("P_1 = ", length(P_1_grid)),  as.character(P_1_grid), sep="")
  p2lab <- paste(rep("P_2 = ", length(P_1_grid)),  as.character(P_2_grid), sep="")
  rownames(tauhat_cols) <- paste(p1lab, p2lab)

  rownames(tauhat_inds) <- rownames(tauhat_cols)
  colnames(tauhat_inds) <- c("Row w/ Max Stat", "Column w/ Max Stat")

  # Retrieving most significant block within each blocking structure at tauhat_max
  max_block_at_change <-matrix(NA, ncol = 2, nrow = length(block_stats_list))
  max_stats_at_change <- c()

  for(i in 1:length(block_stats_list)){
    stat_array <- block_stats_list[[i]]
    if(i ==1){
      max_block_at_change[i,] <- c(1,1)
      max_stats_at_change <- c(max_stats_at_change, as.numeric(stat_array)[tauhat_max])
    } else {
      max_stat_at_change <- max(stat_array[,,tauhat_max])
      max_block_at_change[i,] <- which(stat_array[,,tauhat_max] == max_stat_at_change, arr.ind = T)[1:2]
      max_stats_at_change <- c(max_stats_at_change, max_stat_at_change)
    }
  }

  return(list(tauhat_max = tauhat_max, test_stat = test_stat, p_val = p_val,
              tauhat_cols = tauhat_cols, tauhat_inds = tauhat_inds,
              boot_stats = boot_stats,
              max_stats_at_change = max_stats_at_change, max_blocks_at_change = max_block_at_change))
}
