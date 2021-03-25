cpo <- function(posterior, data, J, q) {
  alpha  <- posterior[, 1]
  lambda   <- posterior[, 2]
  beta <- posterior[, 3:(3 + J)]
  aux <- 0
  for (i in 1:dim(data)[1]) {
    d <- data$d[i]
    t <- data$t[i]
    xe <- data$x[i]
    X <- stats::model.matrix( ~ 1 + xe)
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

    CPO <- (out / dim(posterior)[1]) ^ (-1)

    aux <- aux + log(CPO)

  }
  B <- aux
  return(B)
}
