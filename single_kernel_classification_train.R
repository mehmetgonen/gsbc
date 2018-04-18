single_kernel_classification_train <- function(K, y, parameters) {
  model <- solve_classification_svm(K, y, parameters$C, parameters$epsilon)
  
  state <- list(alpha = model$alpha, b = model$b, parameters = parameters)
}
