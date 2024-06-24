#' Calculate the absolute probability of an index case to case interval belonging to a transmission route component
#'
#' @param dat a (numeric) vector of index case to case intervals
#' @param sigma standard deviation of underlying probability distribution; defaults to sd(dat)
#' @param mu mean of underlying probability distribution; defaults to mean(dat)
#'
#' @return matrix where each row represents a transmission route and each column represents the absolue probability of belonging to that transmission route for each index case to case interval in dat
#' @export
#' @importFrom mitey vec_int
#'
#' @examples
#' my_data<-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))
#' get_abs_prob(my_data)

get_abs_prob <- function(dat, sigma = sd(dat), mu = mean(dat)){

  result_mat <- matrix(, nrow = 7, ncol = length(dat))

  for(i in 1:7){

    tmp_zero <- vec_int(dat, l = dat, u = (dat + 1), s = sigma, m = mu, rt = i, q = "zero")
    tmp_lower <- vec_int(dat, l = (dat-1), u = dat, s = sigma, m = mu, rt = i, q = "lower")
    tmp_upper <- vec_int(dat, l = dat, u = (dat + 1), s = sigma, m = mu, rt = i, q = "upper")

    result_mat[i,] <- c(tmp_zero[which(dat == 0)], tmp_lower[which(dat != 0)] + tmp_upper[which(dat != 0)])
  }

  return(result_mat)
}
