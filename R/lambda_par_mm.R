lambda_par_MM <- function(mu = 0.5,
                          alpha = 0.5,
                          theta = 2,
                          q = 0.5) {
  
  (alpha * log(1 - q ^ (1 / theta))) / (1 - exp(alpha * mu))
}