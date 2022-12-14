#' Generalized Gompertz distribution
#'
#' This function returns the density function value given y, the parameters
#'  alpha, lambda and mu and q, considering a mixture model.
#'
#' @param y the value for which one intends to calculate the density function.
#' @param alpha scale parameter. It must be greater than zero.
#' @param lambda scale parameter. It must be greater than zero.
#' @param mu location parameter
#' @param q quantile of interest.
#' @param log if the log density is to be calculated or not. Default is
#'  FALSE.
#' @return The density value for given y, the parameters alpha, lambda and 
#'  mu and q.
#' @export

dGG_mm <- function(y, alpha, lambda, mu, q, log = FALSE){
  
  logfy  <-
    log(theta) + log(lambda) + (theta - 1) * 
    log(1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) - 
    (-alpha ^ 2 * y + lambda * exp(alpha * y) - lambda) / alpha
  
  if (!log)
    fy <- exp(logfy)
  else
    fy <- logfy
  fy
}


