#' Generalized Gompertz distribution
#'
#' This function returns the density function value given y, the parameters
#'  alpha, lambda and mu and q.
#'
#' @param y the value for which one intends to calculate the density function.
#' @param alpha scale parameter. It must be greater than zero.
#' @param lambda scale parameter. It must be greater than zero.
#' @param mu location parameter
#' @param q quantile of interest.
#' @param log if the log density is to be calculated or not. Default is
#'  FALSE.
#' @return A
#' @export

dGG <- function(y, alpha, lambda, mu, q, log = FALSE){

  check_parameters(y = y, lambda = lambda,
                   mu = mu, q = q)

  theta <-
    -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                            (exp(
                                                              alpha * mu
                                                            ) - 1))))

  logfy <- log(theta) + log(lambda) + (theta - 1) *
    log(1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) -
    (-alpha ^ 2 * y + lambda * exp(alpha * y) - lambda) / alpha

  if (log == FALSE)
    fy = exp(logfy)
  else
    fy = logfy
  fy
}


