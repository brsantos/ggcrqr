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
#' @param linear_pred_mu the linear predictor formula to be included in the 
#'  model for the conditional quantiles; in the LHS one put the name of 
#'  variable to considered as the response variable and in the RHS one puts the
#'  name of explanatory variables, where all variables are included in data.
#' @param linear_pred_alpha the formula for the linear predictor of alpha. 
#' There is no need to put the LHs of the formula
#' @param linear_pred_lambda the formula for the linear predictor of alpha. 
#' There is no need to put the LHs of the formula    
#' @param d name of variable containing the censoring indicator in data.
#' @param log if the log-likelihood is to be calculated or not. Default is
#'  FALSE.
#' @return the likelihood function value given the provided information.
#' @export

likGG <- function(par, q, data, 
                  linear_pred_mu,
                  linear_pred_alpha, 
                  linear_pred_lambda,
                   d, log = FALSE) {
  
  J1 <- dim(stats::model.matrix(linear_pred_mu, data))[2] 
  J2 <- dim(stats::model.matrix(linear_pred_alpha, data))[2]
  J3 <- dim(stats::model.matrix(linear_pred_lambda, data))[2] 
  
  beta <- par[1:J1]
  gama <- par[(1 + J1):(J1 + J2)]
  vi <- par[(J1 + J2 + 1):(J1 + J2 + J3)]
  
  cens <- get(d, data)
  t <- as.numeric(
    stats::model.extract(stats::model.frame(linear_pred_mu, data), 
                         'response')
  )
  X <- stats::model.matrix(linear_pred_mu, data)
  mu <- exp(tcrossprod(X, t(beta)))
  
  X2 <- stats::model.matrix(linear_pred_alpha, data)
  alpha <- -exp(tcrossprod(X2, t(gama)))
  
  X3 <- stats::model.matrix(linear_pred_lambda, data)
  lambda <- exp(tcrossprod(X3, t(vi)))
  
  
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
