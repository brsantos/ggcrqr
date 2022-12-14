#' Generalized Gompertz distribution
#'
#' This function returns the fit for a Bayesian regression model, considering
#'  a mixtture model, considering a mixture model.
#'
#' @param linear_pred_mu the linear predictor formula to be included in the 
#'  model for the conditional quantiles; in the LHS one put the name of 
#'  variable to considered as the response variable and in the RHS one puts the
#'  name of explanatory variables, where all variables are included in data.
#' @param linear_pred_alpha the formula for the linear predictor of alpha. 
#' There is no need to put the LHS of the formula
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


bayesGG_mm  <- function(linear_pred_mu,
             linear_pred_alpha,
             data,
             q,
             d,
             iter = 1000,
             burn,
             jump,
             guess) {
      
    mon.names <- c("LP")
    J1 <- dim(stats::model.matrix(linear_pred_mu, data))[2] 
    J2 <- dim(stats::model.matrix(linear_pred_alpha, data))[2] 
    
    parm.names <- LaplacesDemon::as.parm.names(list(
      alpha = 0,
      lambda = 0,
      beta = rep(0, J1),
      gama = rep(0, J2)
    ))
    
    pos.alpha  <- grep("alpha", parm.names)
    pos.lambda   <- grep("lambda", parm.names)
    pos.beta <- grep("beta", parm.names)
    pos.gama <- grep("gama", parm.names)
    
    MyData <- list(
      linear_pred_mu = linear_pred_mu,
      linear_pred_alpha = linear_pred_alpha,
      data = data, 
      d = d,
      q = q,
      log = TRUE,
      mon.names = mon.names,
      parm.names = parm.names,
      pos.alpha = pos.alpha,
      pos.lambda = pos.lambda,
      pos.beta = pos.beta,
      pos.gama = pos.gama,
      N = 1
    )
    
    Model <- function(parm, Data) {
      alpha  <- parm[Data$pos.alpha]
      lambda  <- parm[Data$pos.lambda]
      beta <-  parm[Data$pos.beta]
      gama <-  parm[Data$pos.gama]
      
      ### Log(Prior Densities)
      alpha.prior <-
        stats::dnorm(alpha,
                     mean = 0,
                     sd = sqrt(100),
                     log = TRUE) 
      
      lambda.prior <-
        stats::dnorm(lambda,
                     mean = 0,
                     sd = sqrt(100),
                     log = TRUE)
      beta.prior <- sum(LaplacesDemon::dnormv(
        beta,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      
      gama.prior <- sum(LaplacesDemon::dnormv(
        gama,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      
      ### Log-Likelihood
      LL <- likGG_mm(  
        par = c(alpha, lambda, beta, gama),
        linear_pred_mu = Data$linear_pred_mu,
        linear_pred_alpha = Data$linear_pred_alpha,
        q = Data$q,
        d = Data$d,
        data = Data$data, 
        log = Data$log
      )
      
      ### Log-Posterior
      LP <- LL + alpha.prior + lambda.prior + beta.prior + gama.prior
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
    
    
    Fit <- LaplacesDemon(
      Model,
      Data = MyData,
      Initial.Values = guess,
      Covar = NULL,
      Iterations = burn + jump * iter,
      Status = burn + jump * iter,
      Thinning = jump,
      Algorithm = "AMWG",
      Specs = list(B = NULL, n = 0, Periodicity = 200)
    )
    
    Posterior <-
      Fit$Posterior1[(length(seq(jump, burn, jump)) + 
                        1):length(Fit$Posterior1[, 1]), ]
    
    Posterior[,1] <- exp(Posterior[,1])
    Posterior[,2] <- exp(Posterior[,2])
    
    burn_rec <- Fit$Rec.BurnIn.Thinned
    jump_rec <- Fit$Rec.Thinning
    AR <- Fit$Acceptance.Rate
    DIC <- Fit$DIC1[1]
    

    list(
      "post" = Posterior,
      "AR" = AR,
      "DIC" = DIC,
      "rec_burnin" = burn_rec,
      "rec_jump" = jump_rec
    )
}
