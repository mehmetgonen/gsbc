library(Rcplex)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K

  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200

  result <- Rcplex(cvec = rep(-1, N), 
                   Amat = matrix(y, nrow = 1, byrow = TRUE), 
                   bvec = 0, 
                   Qmat = yyK, 
                   lb = rep(0, N),
                   ub = rep(C, N),
                   control = opts,
                   objsense = "min",
                   sense = "E")

  alpha <- result$xopt[1:N]
  alpha[alpha < +C * epsilon] <- 0
  alpha[alpha > +C * (1 - epsilon)] <- +C
  objective <- sum(alpha) - 0.5 * (t(alpha) %*% yyK) %*% (alpha)
  objective <- objective * (objective >= 0)
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & abs(alpha) < C)
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] * (1 - yyK[active_indices, support_indices] %*% alpha[support_indices]))
  } else {
    b <- 0
  }

  model <- list(alpha = alpha * y, b = b, objective = objective)
}
