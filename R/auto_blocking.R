#' Generate automatic blocking structures
#'
#' This function implements an method for automatically generating
#' blocking structures for detecting clustered changes at multiple scales
#'
#' @param first_P Number of blocks in first blocking structure (with largest blocks)
#' @param last_P  Number of blocks in last blocking structure (with smallest blocks)
#' @param S Total number of blocking structures
#' @param d Dimension of analyzed data
#' @param extra_prop The maximum allowable proportion of extra dimensions in last block of a given blocking structure
#' @param search_range Defines how many integers are checked to be a number of blocks in a blocking structure
#' @param poly_order Order of polynomial used to generate initial vector of values to search around for numbers of blocks.
#' @param target_min_block_size Targeted size of block in P_S, used if last_P unspecified
#' @return A vector of integers defining the number of blocks in each blocking structure
#' @export
generate_P_vec <- function(first_P=1, last_P=NULL, S= 5, d = NULL, extra_prop = 0.25, search_range = 3,
                           poly_order = 1.5, target_min_block_size= 20) {

  # In case the blocking structure with smallest blocks is not given, we produce it automatically
  if(is.null(last_P)){
    last_P_target <- floor(d/target_min_block_size)

    # Search around the value v_s (for integer # of blocks with d%%P reasonably small)
    candidates <- (last_P_target - search_range):(last_P_target + search_range)
    candidates <- candidates[candidates > 1]  # Keep only values above P=1

    valid_candidates <- candidates[d %% candidates < extra_prop * floor(d/candidates)]

    # Pick the valid candidate closest to target sequence value...
    if (length(valid_candidates) > 0) {
      last_P <- valid_candidates[which.min(abs(valid_candidates - last_P_target))]
    }
  }

  # Generating sequence on low-order polynomial
  t <- seq(0, 1, length.out = S)
  values <- first_P + (last_P - first_P) * t^poly_order

  # Rounding to nearest whole number
  rounded_values <- round(values)


  if (!is.null(d) && !is.null(extra_prop)) {
    for (i in 1:length(rounded_values)) {
      target <- rounded_values[i]

      # Search around the value v_s (for integer # of blocks with d%%P reasonably small)
      candidates <- (target - search_range):(target + search_range)
      candidates <- candidates[candidates > 0]  # Keep only values above P=1

      # Filter candidates that meet the remainder constraint (don't want much larger remainder blocks)
      # Here floor(d/candidates) is the normal block size, and
      # floor(d/candidates) + d %% candidates is remainder block size...
      valid_candidates <- candidates[d %% candidates < extra_prop * floor(d/candidates)]

      # Pick the valid candidate closest to target sequence value...
      if (length(valid_candidates) > 0) {
        rounded_values[i] <- valid_candidates[which.min(abs(valid_candidates - target))]
      }
    }
  }

  return(rounded_values)
}
