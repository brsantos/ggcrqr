matrix_outputs <- function(cadeia) {
  out <- cbind(t(apply(cadeia, 2, summary)),
               as.matrix(apply(cadeia, 2, sd)),
               t(apply(cadeia, 2, function(x)
                 quantile(x, probs = c(0.025, 0.975)))),
               t(apply(cadeia, 2, function(x)
                 emp.hpd(x, conf = 0.95))))
  out <- as.data.frame(out)
  colnames(out) <-
    c(
      "min",
      "1qt",
      "med",
      "media",
      "3qt",
      "max",
      "dp",
      "emp_2.5",
      "emp_97.5",
      "HPD_2.5",
      "HPD_97.5"
    )
  return(out)
}
