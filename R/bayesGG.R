#' Generalized Gompertz distribution
#'
#' This function returns the fit for a Bayesian regression model.
#'
#' @param linear_pred_mu the linear predictor formula to be included in the 
#'  model for the conditional quantiles; in the LHS one put the name of 
#'  variable to considered as the response variable and in the RHS one puts the
#'  name of explanatory variables, where all variables are included in data.
#' @param linear_pred_alpha the formula for the linear predictor of alpha. 
#' There is no need to put the LHs of the formula
#' @param linear_pred_lambda the formula for the linear predictor of alpha. 
#' There is no need to put the LHs of the formula    
#' @param data data to be calculated the quantile of interest.
#' @param q quantile of interest.
#' @param d name of variable containing the censoring indicator in data.
#' @param iter number of posterior samples. Default value is 1000.
#' @param burn number of iterations to be considered as the burn-in.
#' @param jump thinning number for the MCMC.
#' @param guess initial values for the MCMC algorithm.
#' @return the likelihood function value given the provided information.
#' @export
#' @import LaplacesDemon
#' @import truncnorm

bayesGG  <- function(linear_pred_mu,
                     linear_pred_alpha,
                     linear_pred_lambda,
                     data, q, d, iter = 1000, 
                     burn, jump, guess = NULL) {
  mon.names <- c("LP")
  J1 <- dim(stats::model.matrix(linear_pred_mu, data))[2] 
  J2 <- dim(stats::model.matrix(linear_pred_alpha, data))[2] 
  J3 <- dim(stats::model.matrix(linear_pred_lambda, data))[2] 
  parm.names <- LaplacesDemon::as.parm.names(list(
    beta = rep(0, J1),
    gama = rep(0, J2), 
    vi = rep(0, J3)
  ))
  pos.beta <- grep("beta", parm.names)
  pos.gama <- grep("gama", parm.names)
  pos.vi <- grep("vi", parm.names)
  
  if (is.null(guess)) guess <- rep(0, J1 + J2 + J3)
  
  MyData <- list(
    linear_pred_mu = linear_pred_mu,
    linear_pred_alpha = linear_pred_alpha,
    linear_pred_lambda = linear_pred_lambda,
    data = data,
    d = d,
    q = q,
    log = TRUE,
    mon.names = mon.names,
    parm.names = parm.names, 
    pos.beta = pos.beta,
    pos.gama = pos.gama,
    pos.vi = pos.vi,
    N = 1  
  )

  Model <- function(parm, Data) {
    ### Parameters
    beta <-  parm[Data$pos.beta]
    gama <-  parm[Data$pos.gama]
    vi <-  parm[Data$pos.vi]

    ### Log(Prior Densities)
    beta.prior <- sum(dnormv(
      beta,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    gama.prior <- sum(dnormv(
      gama,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    vi.prior <- sum(dnormv(
      vi,
      mean = 0,
      var = 100,
      log = TRUE
    ))

    ### Log-Likelihood
    LL <-
      likGG(
        par = c(beta, gama, vi),
        q = Data$q,
        d = Data$d,
        data = Data$data,
        linear_pred_mu = Data$linear_pred_mu,
        linear_pred_alpha = Data$linear_pred_alpha,
        linear_pred_lambda = Data$linear_pred_lambda,
        log = Data$log
      )

    ### Log-Posterior
    LP <- LL + beta.prior + gama.prior + vi.prior
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
    Iterations = burn + jump * iter,
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
      "rec_jump" = jump_rec,
      "tau" = q
    )
  return(output)
}
