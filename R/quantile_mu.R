#' Generalized Gompertz distribution
#'
#' This function returns the quantile function for the Generalized Gompertz
#'  distribution given the quantile q, the parameters alpha, lambda and
#'  theta.
#'
#' @param q the value for the quantile to be used in the calculation of the
#'  mu parameter. It must be between 0 and 1 and the default value is 0.5.
#' @param alpha scale parameter. It must be greater than zero and the default
#'  value is -1.
#' @param lambda scale parameter. It must be greater than zero and the default
#'  value is 2.
#' @param theta shape parameter. It must be greater than zero and the default
#'  value is 2.
#' formula
#' @return A numeric value for quantile q, given the other parameters.
#' @export

quantile_mu <- function(q = 0.5, alpha = -1, lambda = 2, theta = 2){

  check_parameters(q, lambda, theta)

  (1/alpha) *
    (log(lambda - alpha * log(1 - q^(1/theta) * (1 - exp(lambda/alpha)))) -
        log(lambda))
}
