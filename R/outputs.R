outputs <- function(cadeia) {
  out <- cbind(t(summary(cadeia)),
               as.matrix(sd(cadeia)),
               t(quantile(cadeia, probs = c(0.025, 0.975))),
               t(emp.hpd(cadeia, conf = 0.95)))
  return(out)
}
