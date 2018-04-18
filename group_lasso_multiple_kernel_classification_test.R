group_lasso_multiple_kernel_classification_test <- function(Km, state) {
  Keta <- calculate_Keta(Km, state$eta)
  f <- Keta %*% state$alpha + state$b

  prediction <- list(f = f)
}
