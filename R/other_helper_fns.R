#' Generates Ebynode object
#'
#' This function constructs a list of vectors of edges for each node
#'
#'
#' @param n Number of observations in the time series being analyzed
#' @param E Edge list of a similarity graph (i.e. can take kmst(y, k=k) as input)
#' @return A list of n vectors, each including edges for each node
#' @export
Ebynode_gen <- function(n, E){
  Ebynode <- vector("list", n)  # A list of length n, one for each node
  for(i in 1:n) {Ebynode[[i]] = rep(0,0)}
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1]) # These for loops construct Ebynode from E
  }
  return(Ebynode)
}

#' Helper function for Permuting observations
#'
#' This function helps to permute the labels of observations for running tests to
#' generate the permutation null distribution for computing a p-value
#'
#' @param B Number of permutations
#' @param n Number of observations in the time series being analyzed
#' @return A list containing:
#' \item{permmatches}{A list of B vectors of permutations of 1:n}
#' \item{perms}{A list of B sets of vectors matching pre-permuted indices to a permutation}
#' @export
perm_func <- function(B,n){
  permmatches <- list()
  perms <- list()
  for(b in 2:B){
    perm = sample(n)
    perms[[b]] <- perm
    permmatch = rep(0,n)
    for(i in 1:n) {permmatch[perm[i]] = i}
    permmatches[[b]] <- permmatch}
  return(list(permmatches, perms))
}


#' Generates a vector of graph-based scan statistics for a single block
#'
#' This function uses the similarity graph computed for a specific block
#' of dimensions to generate a vector of graph-based scan statistics
#' for all times from t = n0 to t = n1
#'
#'
#' @param Ebynode A list of length n of vectors of edges for a given node
#' @param edge_info A list of lists of links for each node in the similarity graph
#' @param n Number of observations in the time series being analyzed
#' @param n0 First observation which a test statistic is computed for
#' @param n1 Last observation which a test statistic is computed for
#' @return A vector of length n1-n0 of graph-based max-type scan statistics
#' @export
stat_row <- function(Ebynode, edge_info, n, n0, n1){

  nodedeg = rep(0,n)  # degree of each node
  for(i in 1:n){nodedeg[i] = length(Ebynode[[i]])}
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2  # Stands for number of edges

  g = rep(1,n) # a vector of 1s of length n
  R = rep(0,n) # a vector of length n for storing the between-group edge count
  R1 = rep(0,n) # a vector for storing the within-group edge count BEFORE the potential change per node
  R2 = rep(0,n) # a vector for storing the within-group edge count AFTER the potential change per node

  # Below, we generate counts of edges to find the desired test stat in each block
  for(i in 1:(n-1)){
    g[i] = 0  # update g
    links = edge_info[[i]] # the list of links for a given node

    if(i==1){ # looking at the first node specifically
      if(length(links)>0){
        R[i] = sum(rep(g[i],length(links)) != g[links]) # R[1] set to the number of links to first node
      } else {
        R[i] = 0
      }
      R1[i]=0 # the "before" group edge count
      R2[i]=nE-length(links) # the "after" group edge count, all edges minus those touching the first node
    } else { # all nodes between the 2nd and (n-1)th
      if(length(links)>0){
        add = sum(rep(g[i],length(links)) != g[links]) # counts number of edges coming after the given node
        subtract = length(links)-add # any edges between current node and any prior to it
        R[i] = R[i-1]+add-subtract # adding new edges between, and subtracting edges that are now only in the "before" group
        R1[i] = R1[i-1]+subtract # adding all edges added to the "before" group
      } else {
        R[i] = R[i-1] # this else loop is in case a node has no edges
        R1[i]=R1[i-1]
      }
    }

    R2[i] = nE-R[i]-R1[i] # by process of elimination, all edges not in between/before groups are in the after group
  }

  tt = 1:(n) # a vector from 1 to n
  temp=n0:n1 # a vector comprising the middle 80% of indices of the series
  mu.t = nE* 2*tt*(n-tt)/(n*(n-1))  # These are attained via combinatorial analysis
  p1.tt = 2*tt*(n-tt)/(n*(n-1))
  p2.tt = tt*(n-tt)*(n-2)/(n*(n-1)*(n-2))
  p3.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
  A.tt = (p1.tt-2*p2.tt+p3.tt)*nE+(p2.tt-p3.tt)*sumEisq+p3.tt*nE^2

  # Obtaining the weighted statistic and max-type statistic.
  # Ultimately ABCD has only been implemented with the max-type statistic
  stat.type <- "m"

  if(stat.type %in% c("m", "w")){
    # Code for determining the max-type test stat, taken from the original gseg function
    E_diff <- nE*((2*tt - n)/n)
    Var_diff <- ((tt*(n-tt))*(sumEisq - ((4*(nE^2))/n)))/(n*(n-1))
    Z_diff <- (R1 -R2 - E_diff)/sqrt(Var_diff)

    E_w <- nE*(((tt-1)*(n-tt-1))/((n-1)*(n-2)))
    Var_w <- ((tt*(tt-1)*(n-tt)*(n - tt - 1))/(n*(n-1)*(n-2)*(n-3)))*
      (nE - (sumEisq/(n-2)) + (2*(nE^2))/((n-1)*(n-2)))

    p_t <- (tt - 1)/ (n-2)
    q_t <- 1- p_t
    R_w <- q_t*R1 + p_t*R2

    Z_w <- (R_w - E_w) / sqrt(Var_w)

    #if(stat.type == "w"){
    #  Z <- Z_w
    #  Z[n] <- 0
    #  return(Z)
    #}

    # Statistic which is ultimately used in the overall ABCD algorithm
    if(stat.type == "m"){
      Z <- pmax(abs(Z_diff), Z_w, na.rm = T)
      Z[n] <- 0
      return(Z)
    }
    # Note that scan vector denoted as Z, as to not relabel everything, but this represents "M" in gSeg
  }
}
