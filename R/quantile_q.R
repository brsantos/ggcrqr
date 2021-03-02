#' Generalized Gompertz distribution
#'
#' This function returns the quantile q for the Generalized Gompertz
#'  distribution given the parameters mu_q, alpha, lambda and theta.
#'
#' @param mu parameter mu. Defatult value is 0.5
#'  mu parameter. It must be between 0 and 1 and the default value is 0.5.
#' @param alpha scale parameter. It must be greater than zero and the default
#'  value is -1.
#' @param lambda scale parameter. It must be greater than zero and the default
#'  value is 2.
#' @param theta shape parameter. It must be greater than zero and the default
#'  value is 2.
#' formula
#' @return A
#' @references
#' @export


quantile_q <- function(mu = 0.5, alpha = -1, lambda = 2, theta = 2){

  check_parameters(mu = mu, alpha = alpha, lambda = lambda, theta = theta)

  exp(theta *
        log(
          (1 - exp((lambda - exp(alpha * mu + log(lambda))) / alpha)) /
              (1 - exp(lambda/alpha))))
}
