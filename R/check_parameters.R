check_parameters <- function(q = NULL, alpha = NULL, lambda = NULL,
                             theta = NULL, mu = NULL, y = NULL){

  if (any(lambda < 0))
    stop(paste("lambda must be positive", "\n", ""))

  if (any(alpha < 0))
    stop(paste("alpha must be positive", "\n", ""))

  if (any(theta < 0))
    stop(paste("mu must be positive", "\n", ""))

  if (any(q < 0) | any(q > 1))
    stop(paste("q must be in the interval (0,1) ", "\n", ""))

  if (any(y <= 0))
    stop(paste("y must be greater than 0 ", "\n", ""))
}
