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
#' @param lower.p if TRUE returns the log probability of the cdf function.
#'  Default is FALSE.
#' @return A
#' @references
#' @export

pGG <- function(y, alpha, lambda, mu, q, lower.tail = TRUE, log.p = FALSE){

  check_parameters(y = y, alpha = alpha, lambda = lambda,
                   mu = mu, q = q)

  if (any(lambda < 0))
    stop(paste("lambda must be positive", "\n", ""))

  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))

  if (any(y <= 0))
    stop(paste("y must be greater than 0 ", "\n", ""))

  if (any(q < 0) | any(q > 1))
    stop(paste("q must be in the interval (0,1) ", "\n", ""))

  theta <-
    -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                            (exp(
                                                              alpha * mu
                                                            ) - 1))))

  cdf <- (1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) ^ theta

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
