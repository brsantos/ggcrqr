#' Generalized Gompertz distribution
#'
#' This function returns the cure fraction value given the parameters
#'  alpha, lambda and theta.
#'
#' @param alpha scale parameter. It must be greater than zero.
#' @param lambda scale parameter. It must be greater than zero.
#' @param mu location parameter.
#' @param q quantile of interest.  
#' @return The cure fraction for given parameters alpha, lambda and mu and q.
#' @export 

p0_dGG <- function(alpha, lambda, mu, q) {

  # check_parameters(y = y, lambda = lambda, mu = mu, q = q)
  
  theta <-
    -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                            (exp(
                                                              alpha * mu
                                                            ) - 1))))

  1 - (1 - exp(lambda / alpha)) ^ theta
}
