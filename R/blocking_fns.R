#' Generate block divisions for high-dimensional data
#'
#' This function computes divisions between blocks
#' in a single blocking structure for high-dimensional data
#'
#' @param d Dimension of time series to be analyzed
#' @param P Number of blocks to divide time series into
#' @return A vector of dimension indices used to divide time series into blocks
#' @export
breaks_1d <- function(d, P){
  block_size <- floor(d/P) # the size of the first P-r blocks
  breaks <- c(seq(from = 1, to = block_size*(P-1) + 1, by = block_size), d+1)
  return(breaks)
}


#' Generate block divisions for image data
#'
#' This function computes divisions between blocks
#' in a single blocking structure for image data
#'
#' @param side_length_1 Number of rows of image dimensions
#' @param side_length_2 Number of columns of image dimensions
#' @param P_1 Number of row divisions
#' @param P_2 Number of column divisions
#' @return A list of two vectors of row/column indices used to divide the image series into blocks
#' \item{breaks_1}{A vector of row indices used to divide the image series into blocks}
#' \item{breaks_2}{A vector of column indices used to divide the image series into blocks}
#' @export
breaks_img <- function(side_length_1, side_length_2, P_1, P_2){
  # blocking divisions for rows
  block_length_1 <- floor(side_length_1/P_1) # the number of rows of dimensions in  the first P_1-1 rows of blocks

  # blocking divisions for columns
  block_length_2 <- floor(side_length_2/P_2) # the number of rows of dimensions in  the first P_2-1 cols of blocks


  breaks_1 <- c(seq(from = 1,
                    to = side_length_1 - block_length_1 + 1,
                    by = block_length_1), side_length_1+1)

  breaks_2 <- c(seq(from = 1,
                    to = side_length_2 - block_length_2 + 1,
                    by = block_length_2), side_length_2+1)
  return(list(breaks_1 = breaks_1, breaks_2 = breaks_2))
}
