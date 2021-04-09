outputs <- function(chain) {
  out <- as.data.frame(cbind(t(summary(chain)),
               as.matrix(stats::sd(chain)),
               t(stats::quantile(chain, probs = c(0.025, 0.975))),
               t(TeachingDemos::emp.hpd(chain, conf = 0.95))))
  colnames(out) <-
    c(
      "min",
      "1qt",
      "med",
      "media",
      "3qt",
      "max",
      "sd",
      "emp_025",
      "emp_975",
      "HPD_025",
      "HPD_975"
    )
  out
}
