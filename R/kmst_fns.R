#' Compute pairwise L2 Distance
#'
#' This function computes the L2 distance between all pairs of
#' observations in a given time series of length n
#'
#' @param y Time series of dimension n by d
#' @return An n by n matrix of pairwise L2 distances between observations
#' @export
getdis = function(y){ #Code for L2 distance
  n = dim(y)[1]
  G = y%*%t(y)
  g = diag(G)
  dis = sqrt(matrix(rep(g,n),n) + matrix(rep(g,n),n,byrow=T) - 2*G)
  dis
}

#' Construct k-MST on time series
#'
#' This function computes the L2 distance between all pairs of observations
#' in a given time series of length n, output as an n by n matrix
#'
#' @param y Time series of dimension n by d
#' @param dis n by n matrix of pairwise L2 distances between observations
#' @param k  Controls number of edges k*(n-1) of k-MST similarity graph
#' @return  An n by 2 matrix of observations pairs connected by an edge
#' @export
kmst = function(y=NULL, dis=NULL, k=1){
  if (is.null(dis) && is.null(y)){
    cat("Please input data or the distance matrix!\n")
    return(0)
  }
  if (is.null(dis)){
    dis = getdis(y)
  }

  mymst = ade4::mstree(stats::as.dist(dis),k)
  cbind(mymst[,1], mymst[,2])
}
