#' Generalized Gompertz distribution
#'
#' This function returns the quantile function for the Generalized Gompertz
#'  distribution given the quantile q, the parameters alpha, lambda and
#'  theta, considering a mixture model
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

quantile_mu_mm <- function(q = 0.5, 
                           alpha = 0.5, 
                           lambda = 2, 
                           theta = 2){

  (1 / alpha) * (log(1 - (alpha / lambda) * log(1 - q ^ (1 / theta))))
}