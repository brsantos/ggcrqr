outputs <- function(chain) {
  out <- cbind(t(summary(chain)),
               as.matrix(stats::sd(chain)),
               t(stats::quantile(chain, probs = c(0.025, 0.975))),
               t(TeachingDemos::emp.hpd(chain, conf = 0.95)))
  return(out)
}
