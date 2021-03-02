theta_par <- function(mu = 0.5,
                      alpha = -1,
                      lambda = 2,
                      q = 0.5) {
  log(q) /
    (log(1 - exp((
      lambda - exp(alpha * mu + log(lambda))
    ) / alpha)) - log(1 - exp(lambda / alpha)))
}
