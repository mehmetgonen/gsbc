library(Rmosek)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K

  problem <- list()
  problem$sense <- "min"
  problem$c <- rep(-1, N)
  problem$A <- Matrix(y, nrow = 1, byrow = TRUE, sparse = TRUE)
  problem$bc <- rbind(blc = 0, buc = 0) 
  problem$bx <- rbind(blx = rep(0, N), bux = rep(C, N))

  I <- matrix(1:N, N, N, byrow = FALSE)
  J <- matrix(1:N, N, N, byrow = TRUE)
  problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
                       j = J[lower.tri(J, diag = TRUE)],
                       v = yyK[lower.tri(yyK, diag = TRUE)])

  opts <- list()
  opts$verbose <- 0
  result <- mosek(problem, opts)

  alpha <- result$sol$itr$xx[1:N]
  alpha[alpha < +C * epsilon] <- 0
  alpha[alpha > +C * (1 - epsilon)] <- +C
  objective <- sum(alpha) - 0.5 * (t(alpha) %*% yyK) %*% (alpha)
  objective <- objective * (objective >= 0)
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & alpha < C)
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] * (1 - yyK[active_indices, support_indices] %*% alpha[support_indices]))
  } else {
    b <- 0
  }

  model <- list(alpha = alpha * y, b = b, objective = objective)
}
