#' Generalized Gompertz distribution
#'
#' This function returns the likelihood function value given data, the parameters
#'  alpha, lambda and mu and q.
#'
#' @param par vector with dimension equal to the
#' @param q quantile of interest.
#' @param data data to be calculated the quantile of interest. It must have a
#'  column with d value containing a censoring indicator, where 1 corresponds
#'  to a non-censored observation and 0 otherwise. It must also contain a x
#'  column with
#' @param J number of regressors in the model.
#' @param log if the log-likelihood is to be calculated or not. Default is
#'  FALSE.
#' @return the likelihood function value given the provided information.
#' @export

likGG <- function(par, q, data, J, log = FALSE) {
  alpha <- par[1]
  lambda <- par[2]
  beta <- par[3:(3 + J)]

  cens <- data$d
  t <- data$t
  xe <- data$x
  X <- stats::model.matrix( ~ 1 + xe)
  mu <- exp(tcrossprod(X, t(beta)))
  out <-
    sum(
      cens * dGG(
        t,
        alpha = alpha,
        lambda = lambda,
        mu = mu,
        q = q,
        log = TRUE
      ) +
        (1 - cens) * pGG(
          t,
          alpha = alpha,
          lambda = lambda,
          mu = mu,
          q = q,
          lower.tail = FALSE,
          log.p = TRUE
        )
    )

  if (log == TRUE)
    lik <- out
  else
    lik <- exp(out)
  return(lik)
}
