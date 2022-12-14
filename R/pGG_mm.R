#' Generalized Gompertz distribution
#'
#' This function returns the cdf function value given y, the parameters
#'  alpha, lambda and mu and q, considering a mixture model.
#'
#' @param y the value for which one intends to calculate the density function.
#' @param alpha scale parameter. It must be greater than zero.
#' @param lambda scale parameter. It must be greater than zero.
#' @param mu shape parameter. It must be greater than zero.
#' @param q quantile of interest.
#' @param lower.tail if TRUE returns the lower tail of the cdf function.
#'  Otherwise, returns the probability of being greater than y.
#' @param log.p if TRUE returns the log of the cdf. Default if FALSE.
#' @return A
#' @export

pGG <- function(y, alpha, lambda, mu, q, lower.tail = TRUE, log.p = FALSE){
  
  cdf   <-
    (1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) ^ theta
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  cdf
}
