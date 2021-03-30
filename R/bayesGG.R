#' Generalized Gompertz distribution
#'
#' This function returns the fit for a Bayesian regression model.
#'
#' @param linear_pred the linear predictor formula to be included in the model,
#'  in the LHS one put the name of variable to considered as the response
#'  variable and in the RHS one puts the name of explanatory variables,
#'  where all variables are included in data.
#' @param data data to be calculated the quantile of interest.
#' @param q quantile of interest.
#' @param d name of variable containing the censoring indicator in data.
#' @param burn number of iterations to be considered as the burn-in.
#' @param jump thinning number for the MCMC.
#' @param guess initial values for the MCMC algorithm.
#' @return the likelihood function value given the provided information.
#' @export
#' @import LaplacesDemon
#' @import truncnorm

bayesGG  <- function(linear_pred, data, q, d, burn, jump, guess) {
  n.size <- 1000
  mon.names <- c("LP")
  J <- dim(stats::model.matrix(linear_pred, data))[2] - 1
  parm.names <- LaplacesDemon::as.parm.names(list(
    alpha = 0,
    lambda = 0,
    beta = rep(0, (J + 1))
  ))
  pos.alpha  <- grep("alpha", parm.names)
  pos.lambda   <- grep("lambda", parm.names)
  pos.beta <- grep("beta", parm.names)

  MyData <- list(
    linear_pred = linear_pred,
    data = data,
    J = J,
    q = q,
    d = d,
    log = TRUE,
    mon.names = mon.names,
    parm.names = parm.names,
    pos.alpha = pos.alpha,
    pos.lambda = pos.lambda,
    pos.beta = pos.beta,
    N = 1
  )

  Model <- function(parm, Data) {
    ### Parameters
    alpha  <- LaplacesDemon::interval(parm[Data$pos.alpha], -Inf, -1e-100)
    parm[Data$pos.alpha]  <- alpha
    lambda  <-  LaplacesDemon::interval(parm[Data$pos.lambda], 1e-100, Inf)
    parm[Data$pos.lambda] <- lambda
    beta <-  parm[Data$pos.beta]

    ### Log(Prior Densities)
    alpha.prior <-
      log(truncnorm::dtruncnorm(
        alpha,
        a = -Inf,
        b = -1e-100,
        mean = -0.1,
        sd = 10
      ))
    lambda.prior <-
      stats::dgamma(lambda,
             shape = 0.01,
             rate = 0.01,
             log = TRUE) #mean 1 and var 100
    beta.prior <- sum(LaplacesDemon::dnormv(
      beta,
      mean = 0,
      var = 100,
      log = TRUE
    ))

    ### Log-Likelihood
    LL <-
      likGG(
        par = c(alpha, lambda, beta),
        q = Data$q,
        d = Data$d,
        data = Data$data,
        linear_pred = Data$linear_pred,
        J = Data$J,
        log = Data$log
      )

    ### Log-Posterior
    LP <- LL + alpha.prior + lambda.prior + beta.prior
    Modelout <-
      list(
        LP = LP,
        Dev = -2 * LL,
        Monitor = LP,
        yhat = 1,
        parm = parm
      )
    return(Modelout)
  }

  Fit <- LaplacesDemon::LaplacesDemon(
    Model = Model,
    Data = MyData,
    Initial.Values = guess,
    Covar = NULL,
    Iterations = burn + jump * n.size,
    Status = 20000,
    Thinning = jump,
    Algorithm = "AM",
    Specs = list(Adaptive = 500, Periodicity = 100),
    Debug=list(DB.chol = FALSE, DB.eigen = FALSE,
               DB.MCSE = FALSE, DB.Model = FALSE)
  )


  Posterior <-
    Fit$Posterior1[
      (length(seq(jump, burn, jump)) + 1):length(Fit$Posterior1[, 1]),
      ]
  burn_rec <- Fit$Rec.BurnIn.Thinned
  jump_rec <- Fit$Rec.Thinning
  AR <- Fit$Acceptance.Rate
  DIC <- Fit$DIC1[1]

  output <-
    list(
      "post" = Posterior,
      "AR" = AR,
      "DIC" = DIC,
      "rec_burnin" = burn_rec,
      "rec_jump" = jump_rec
    )
  return(output)
}
