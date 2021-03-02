cpo <- function(Posterior, dados, J, q) {
  #média do log do CPO

  alpha  <- Posterior[, 1]
  lambda   <- Posterior[, 2]
  beta <- Posterior[, 3:(3 + J)]
  aux <- 0
  for (i in 1:dim(dados)[1]) {
    d <- dados$d[i]
    t <- dados$t[i]
    xe <- dados$x[i]
    X <- model.matrix( ~ 1 + xe)
    x_beta <- do.call(rbind, lapply(apply(beta, 1, list), function(b) {
      tcrossprod(X, t(b[[1]]))
    }))
    mu <- exp(x_beta)
    out <-
      sum(1 / ((
        dGG(
          t,
          alpha = alpha,
          lambda = lambda,
          mu = mu,
          q = q,
          log = FALSE
        ) ^ d
      ) *
        (
          pGG(
            t,
            alpha = alpha,
            lambda = lambda,
            mu = mu,
            q = q,
            lower.tail = FALSE,
            log.p = FALSE
          ) ^ (1 - d)
        )))

    CPO <- (out / dim(Posterior)[1]) ^ (-1)

    aux <- aux + log(CPO)

  }
  B <- aux #/dim(dados)[1]
  return(B)
}
