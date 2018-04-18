library(AUC)
library(randomForestSRC)

source("classification_helper.R")

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

data_path <- "./data"
result_path <- "./results_1_vs_234"
pathway <- "none" #replace with "hallmark" if you would like to restrict RF to Hallmark genes only

for (cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for (replication in 1:100){
    if (file.exists(sprintf("%s/%s/random_forest_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
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
      
      ntree_set <- c(1:5) * 500
      fold_count <- 4
      train_ratio <- 0.8
      
      if (pathway != "none") {
        pathways <- read_pathways(pathway)
        gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
        X <- X[, which(colnames(X) %in% gene_names)]
      }
      
      set.seed(1606 * replication)
      train_negative_indices <- sample(negative_indices, ceiling(train_ratio * length(negative_indices)))
      train_positive_indices <- sample(positive_indices, ceiling(train_ratio * length(positive_indices)))
      
      auroc_matrix <- matrix(NA, nrow = fold_count, ncol = length(ntree_set), dimnames = list(1:fold_count, sprintf("%g", ntree_set)))
      
      negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
      positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))
      for (fold in 1:fold_count) {
        train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
        test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
        
        X_train <- X[train_indices,]
        X_test <- X[test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        X_train[is.na(X_train)] <- 0
        X_test[is.na(X_test)] <- 0
        
        y_train <- y[train_indices]
        y_test <- y[test_indices]
        
        data <- as.data.frame(X_train)
        data$label <- as.factor(y_train)
        for (ntree in ntree_set) {
          print(sprintf("running fold = %d, ntree = %g", fold, ntree))
          state <- rfsrc(formula = label ~ ., data = data, ntree = ntree)
          prediction <- predict(state, newdata = as.data.frame(X_test))
          auroc_matrix[fold, sprintf("%g", ntree)] <- auc(roc(prediction$predicted[,2], as.factor(1 * (y_test == +1))))
        }
      }
      
      ntree_star_AUROC <- ntree_set[max.col(t(colMeans(auroc_matrix, na.rm = TRUE)), ties.method = "last")]
      
      train_indices <- c(train_negative_indices, train_positive_indices)
      test_indices <- setdiff(1:length(y), train_indices)
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      X_train[is.na(X_train)] <- 0
      X_test[is.na(X_test)] <- 0
      
      y_train <- y[train_indices]
      y_test <- y[test_indices]
      
      data <- as.data.frame(X_train)
      data$label <- as.factor(y_train)
      state <- rfsrc(formula = label ~ ., data = data, ntree = ntree_star_AUROC)
      prediction <- predict(state, newdata = as.data.frame(X_test), tree.err = TRUE)
      result <- list()
      result$AUROC <- auc(roc(prediction$predicted[,2], as.factor(1 * (y_test == +1))))

      save("result", file = sprintf("%s/%s/random_forest_pathway_%s_measure_AUROC_replication_%d_result.RData", result_path, cohort, pathway, replication))
    } 
  }
}
