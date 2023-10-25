#' Generalized Gompertz distribution
#'
#' This function returns the cdf function value given y, the parameters
#'  alpha, lambda and mu and q.
#'
#' @param y the value for which one intends to calculate the density function.
#' @param alpha scale parameter. It must be greater than zero.
#' @param lambda scale parameter. It must be greater than zero.
#' @param mu shape parameter. It must be greater than zero.
#' @param q quantile of interest.
#' @param lower.tail if TRUE returns the lower tail of the cdf function.
#'  Otherwise, returns the probability of being greater than y.
#' @param log.p if TRUE returns the log of the cdf. Default if FALSE.
#' @return The cdf for given y, for parameters alpha, lambda and mu and q.
#' @export

pGG <- function(y, alpha, lambda, mu, q, lower.tail = TRUE, log.p = FALSE){

  check_parameters(y = y, lambda = lambda, mu = mu, q = q)

  theta <-
    -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                            (exp(
                                                              alpha * mu
                                                            ) - 1))))

  cdf <- (1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) ^ theta

  if (lower.tail) cdf <- cdf
  else cdf <- 1 - cdf
  if (!log.p) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
