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
#' @param linear_pred the linear predictor formula to be included in the model,
#'  in the LHS one put the name of variable to considered as the response
#'  variable and in the RHS one puts the name of explanatory variables,
#'  where all variables are included in data.
#' @param J number of regressors in the model.
#' @param log if the log-likelihood is to be calculated or not. Default is
#'  FALSE.
#' @return the likelihood function value given the provided information.
#' @export

likGG <- function(par, q, data, linear_pred, J, log = FALSE) {
  alpha <- par[1]
  lambda <- par[2]
  beta <- par[3:(3 + J)]

  cens <- data$d
  # t <- data$t
  # xe <- as.matrix(data[,4:(4+J-1)])
  t <- as.numeric(
    stats::model.extract(stats::model.frame(linear_pred, data), 'response')
  )
  X <- stats::model.matrix(linear_pred, data)
  # X <- stats::model.matrix( ~ 1 + xe)
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
