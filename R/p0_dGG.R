#' Generalized Gompertz distribution
#'
#' This function returns the cure fraction value given the parameters
#'  alpha, lambda and theta.
#'
#' @param alpha scale parameter. It must be greater than zero and the default
#'  value is -1.
#' @param lambda scale parameter. It must be greater than zero and the default
#'  value is 2.
#' @param theta shape parameter. It must be greater than zero and the default
#'  value is 2.
#' formula
#' @return A
#' @export

p0_dGG <- function(alpha, lambda, theta) {
  check_parameters(lambda = lambda, theta = theta)

  1 - (1 - exp(lambda / alpha)) ^ theta
}
