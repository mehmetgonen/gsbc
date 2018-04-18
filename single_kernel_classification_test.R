single_kernel_classification_test <- function(K, state) {
  f <- K %*% state$alpha + state$b
  
  prediction <- list(f = f)
}
