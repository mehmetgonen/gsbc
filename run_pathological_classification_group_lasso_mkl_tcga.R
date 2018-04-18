library(AUC)
optimizer <- "mosek"

source("classification_helper.R")
if (optimizer == "cplex") {
  source("solve_classification_svm_cplex.R")
}
if (optimizer == "mosek") {
  source("solve_classification_svm_mosek.R")
}
source("group_lasso_multiple_kernel_classification_train.R")
source("group_lasso_multiple_kernel_classification_test.R")

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

data_path <- "./data"
result_path <- "./results_1_vs_234"
pathway <- "hallmark"

for (cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for (replication in 1:100) {
    if (file.exists(sprintf("%s/%s/glmkl_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
      load(sprintf("%s/%s.RData", data_path, cohort))
      
      common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$pathologic_stage) == FALSE)], rownames(TCGA$mrna))
      
      X <- log2(TCGA$mrna[common_patients,] + 1)
      y <- rep(NA, length(common_patients))
      
      y[TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
      y[TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                  "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                  "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
      
      valid_patients <- which(is.na(y) == FALSE)
      valid_features <- as.numeric(which(apply(X[valid_patients,], 2, sd) != 0))
      X <- X[valid_patients, valid_features]
      y <- y[valid_patients]
      
      negative_indices <- which(y == -1)
      positive_indices <- which(y == +1)
      
      C_set <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000)
      epsilon <- 1e-5
      fold_count <- 4
      train_ratio <- 0.8
      iteration_count <- 200
      
      pathways <- read_pathways(pathway)
      gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
      X <- X[, which(colnames(X) %in% gene_names)]
      
      auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(C_set), dimnames = list(1:fold_count, sprintf("%g", C_set)))
      
      set.seed(1606 * replication)
      train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
      train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))
      
      negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
      positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))
      for (fold in 1:fold_count) {
        train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
        test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
        
        X_train <- X[train_indices,]
        X_test <- X[test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        
        N_train <- nrow(X_train)
        N_test <- nrow(X_test)
        N_pathway <- length(pathways)
        K_train <- array(0, dim = c(N_train, N_train, N_pathway))
        K_test <- array(0, dim = c(N_test, N_train, N_pathway))
        for (m in 1:N_pathway) {
          feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
          D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
          sigma <- mean(D_train)
          K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
          K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
        }
        
        y_train <- y[train_indices]
        y_test <- y[test_indices]
        
        for (C in C_set) {
          print(sprintf("running fold = %d, C = %g", fold, C))
          parameters <- list()
          parameters$C <- C
          parameters$epsilon <- epsilon
          parameters$iteration_count <- iteration_count
          
          state <- group_lasso_multiple_kernel_classification_train(K_train, y_train, parameters)
          prediction <- group_lasso_multiple_kernel_classification_test(K_test, state)
          auroc_matrix[fold, sprintf("%g", C)] <- auc(roc(prediction$f, as.factor(1 * (y_test == +1))))
        }
      }
      
      C_star_AUROC <- C_set[max.col(t(colMeans(auroc_matrix)), ties.method = "last")]    
      
      train_indices <- c(train_negative_indices, train_positive_indices)
      test_indices <- setdiff(1:length(y), train_indices)
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      N_pathway <- length(pathways)
      K_train <- array(0, dim = c(N_train, N_train, N_pathway))
      K_test <- array(0, dim = c(N_test, N_train, N_pathway))
      for (m in 1:N_pathway) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train <- y[train_indices]
      y_test <- y[test_indices]
      
      parameters <- list()
      parameters$C <- C_star_AUROC
      parameters$epsilon <- epsilon
      parameters$iteration_count <- iteration_count
      
      state <- group_lasso_multiple_kernel_classification_train(K_train, y_train, parameters)
      prediction <- group_lasso_multiple_kernel_classification_test(K_test, state)
      result <- list()
      result$AUROC <- auc(roc(prediction$f, as.factor(1 * (y_test == +1))))
      
      save("state", file = sprintf("%s/%s/glmkl_pathway_%s_measure_AUROC_replication_%d_state.RData", result_path, cohort, pathway, replication))
      save("result", file = sprintf("%s/%s/glmkl_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication))
    } 
  }
}
