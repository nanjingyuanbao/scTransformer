#####requirement
###install
# packages_needed <- c("ggplot2", "dplyr", "tidyr","Seurat","SingleCellNet","SingleR","Clustifyr","CHETAH","scmap","scLearn","scPred","SciBet","caret")
# sapply(packages_needed, function(pkg) {
#  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#    install.packages(pkg)
#  }
#})

packages <- c("ggplot2", "dplyr", "tidyr","Seurat","SingleCellNet","SingleR","Clustifyr","CHETAH","scmap","scLearn","scPred","SciBet","caret")
sapply(packages, library, character.only = TRUE)

# Define 5 k-folds
X <- as.matrix(lung1@assays$RNA@counts)
y <- lung1@meta.data$True_Class
set.seed(128) # set randoms
folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)
accuracies <- c()
f1_scores <- c()
F1_Score <- function(actual, predicted) {
  actual <- factor(actual, levels = levels(factor(predicted)))
  predicted <- factor(predicted, levels = levels(factor(predicted)))
  cm <- confusionMatrix(predicted, actual)
  f1 <- cm$byClass[4]
  return(f1)
}
calculate_metrics <- function(y_test, y_pred) {
  # Calculate the accuracy of each cell
  accuracy <- sapply(unique(y_test), function(c) {
    mean(y_test[y_pred == c] == c)
  })
  
  # Compute the f1 score
  f1 <- sapply(unique(y_test), function(c) {
    actual <- factor(y_test == c, levels = c(FALSE, TRUE))
    predicted <- factor(y_pred == c, levels = c(FALSE, TRUE))
    cm <- confusionMatrix(predicted, actual)
    f1 <- cm$byClass[4]
    return(f1)
  })
  
  # Combine accuracy and f1 score into one dataframe
  metrics <- data.frame(accuracy, f1)
  rownames(metrics) <- unique(y_test)
  
  # Returns a list of dataframes and category names with accuracy and f1 scores
  return(list(metrics, unique(y_test)))
}
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  #dim(X_train)
  #dim(X_test)
  #length(y_train)
  #length(y_test)
  train1 <- Seurat::CreateSeuratObject(counts = X_train)
  train1 <- NormalizeData(object =train1)
  train1 <- FindVariableFeatures(object = train1)
  train1<- ScaleData(object =train1)
  train1 <- RunPCA(object =train1)
  train1 <- FindNeighbors(object = train1)
  train1 <- FindClusters(object = train1)
  train1 <- RunTSNE(object =train1)
  train1<- RunUMAP(object = train1,dim=1:10 )
  train1@meta.data$True_Class1=y_train
  train1.sce= as.SingleCellExperiment(train1)
  train_seurat[[i]]=train1
  train_singlecellnet[[i]]=train1.sce
  test1 <- Seurat::CreateSeuratObject(counts = X_test)
  test1 <- NormalizeData(object =test1)
  test1 <- FindVariableFeatures(object = test1)
  test1<- ScaleData(object =test1)
  test1 <- RunPCA(object =test1)
  test1 <- FindNeighbors(object = test1)
  test1 <- FindClusters(object = test1)
  test1 <- RunTSNE(object =test1)
  test1<- RunUMAP(object = test1,dim=1:10 )
  test1@meta.data$True_Class=y_test
  test1.sce= as.SingleCellExperiment(test1)
  test_seurat[[i]]=test1
  test_singlecellnet[[i]]=test1.sce
}

###############singler########
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
all_metrics <- list()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  #X_train = as.SingleCellExperiment(X_train)
  #X_test= as.SingleCellExperiment(X_test)
  # Prediction using SingleR
  sr <- SingleR(test=test_singlecellnet[[i]],ref=train_singlecellnet[[i]],labels=y_train)
  y_pred <- sr$labels
  
  # Calculate the accuracy and f1 score
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  
  #Calculate the accuracy  and f1-score of each category
  
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}

# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_SingleR <- mean(accuracies)
mean_f1_score_Dolung_SingleR <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_SingleR * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_SingleR * 100, 2), "%\n")
###Mean accuracy: 95.56 %
###Mean F1 score: 95.28 %
# Returns a list of dataframes and category names with accuracy and f1 scores
accuracy_list <- c()
f1_list <- c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
}

# Combine accuracy and f1 scores into one dataframe
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class = c(row.names(all_metrics[[1]][[1]]),row.names(all_metrics[[2]][[1]]),row.names(all_metrics[[3]][[1]]),row.names(all_metrics[[4]][[1]]),row.names(all_metrics[[5]][[1]])))
df_all$class <- factor(df_all$class, levels = unique(y))
# Calculate average accuracy and average F1 score grouped by category name
df_avg_Dolung_SingleR <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
#Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_SingleR$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_SingleR$f1)

###########CHETAH######################
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  #Cell type of reference data
  celltypes_hn <- y_train
  #Expression matrix for reference data
  counts_hn <- assay(train_singlecellnet[[i]])
  #Expression matrix for input data
  counts_melanoma=assay(test_singlecellnet[[i]])
  #Dimensionality reduction information for input data
  tsne_melanoma <- reducedDim(test_singlecellnet[[i]],"TSNE")
  #Construct SingleCellExperiment object for reference data
  headneck_ref <- SingleCellExperiment(assays = list(counts = counts_hn),
                                       colData = DataFrame(celltypes = celltypes_hn))
  #input
  input_mel <- SingleCellExperiment(assays = list(counts = counts_melanoma),
                                    reducedDims = SimpleList(TSNE = tsne_melanoma))
  #predict
  input_mel <- CHETAHclassifier(input = input_mel,
                                ref_cells = headneck_ref)
  input_mel <- Classify(input_mel, 0)
  y_pred=input_mel$celltype_CHETAH
  # Calculate  accuracy and  F1 score 
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_CHETAH <- mean(accuracies)
mean_f1_score_Dolung_CHETAH <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_CHETAH * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_CHETAH * 100, 2), "%\n")
#Mean accuracy: 93.63 %
#Mean F1 score: 92.78 %
#Calculate the unassigned rate
mean_unassigned_Dolung_CHETAH <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_CHETAH * 100, 2), "%\n")
#Mean Unssigned score:0 %
# Save accuracy and f1 scores for all test sets to a list
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) { 
  metrics <- all_metrics[[i]]
  
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
  
  
}
# Combine accuracy and f1 scores into one dataframe
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class = row_names)
df_all$class <- factor(df_all$class, levels = unique(y))
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_CHETAH <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
#Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_CHETAH$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_CHETAH$f1)

#######################clustifyr################
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  # Get unique cell type labels
  unique_cell_types <- unique(y_train)
  # Get a list of cell types
  cell_types <-  unique(y_train)
  gene_names <- rownames(train_singlecellnet[[i]])
  gene_expression_matrix <- data.frame(Gene = gene_names)
  for (cell_type in cell_types) {
    gene_expression_matrix[cell_type] <- 0
  }
  X1=assay(train_singlecellnet[[i]])
  for (cell_type in cell_types) {
    x=(y_train == cell_type)
    subset_data <- X1[, x]
    gene_expression_values <- rowMeans(subset_data)  
    gene_expression_matrix[cell_type] <- gene_expression_values
  }
  rownames(gene_expression_matrix)=gene_expression_matrix[,1]
  gene_expression_matrix=gene_expression_matrix[,-1]
  res <- clustify(
    input = test_seurat[[i]],       # a Seurat object
    ref_mat = gene_expression_matrix,    # matrix of RNA-seq expression data for each cell type
    cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
    obj_out = TRUE,
    per_cell = FALSE,# output Seurat object with cell type inserted as "type" column
  )
  
  y_pred=res$type
  # Calculate accuracy and f1 scores
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate Save accuracy and f1 scores for each 
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_Clustifyr <- mean(accuracies)
mean_f1_score_Dolung_Clustifyr <- mean(f1_scores[!is.na(f1_scores)])
cat("Mean accuracy:", round(mean_accuracy_Dolung_Clustifyr * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_Clustifyr * 100, 2), "%\n")
#Mean accuracy: 81.21 %
#Mean F1 score: 83.19 %
#Calculate the unassigned rate
mean_unassigned_Dolung_Clustifyr <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_Clustifyr * 100, 2), "%\n")
#Mean Unssigned score: 4.79 %
# Save accuracy and f1 scores for all test sets to a list
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_Clustifyr <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_Clustifyr$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_Clustifyr$f1)

##############scibet##########################
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  train_set=t(assay(train_singlecellnet[[i]]))
  train_set=as.data.frame(train_set)
  train_set$label=as.character(y_train)
  test_set=t(assay(test_singlecellnet[[i]]))
  test_set=as.data.frame(test_set)
  y_pred=SciBet_R(train_set,test_set)
  # Calculate accuracy and f1 scores
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #计算未分配率
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_scibetr <- mean(accuracies)
mean_f1_score_Dolung_scibetr <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_scibetr * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_scibetr * 100, 2), "%\n")
#Mean accuracy: 96.17 %
#Mean F1 score: 100 %
#Calculate the unassigned rate
mean_unassigned_Dolung_scibetr <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_scibetr * 100, 2), "%\n")
#Mean Unssigned score: 0 %
# Returns a list of dataframes and category names with accuracy and f1 scores
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 score into one data frame
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Calculate average accuracy and average F1 score grouped by category name
df_avg_Dolung_Scibet <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_Scibet$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_Scibet$f1)


#######################scLearn###################
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  rawcounts=assay(train_singlecellnet[[i]])
  refe_ann<-as.character(y_train)
  names(refe_ann)<-colnames(train_singlecellnet[[i]])
  data_qc<-Cell_qc(rawcounts,refe_ann,species="Mm")
  data_type_filtered<-Cell_type_filter(data_qc$expression_profile,data_qc$sample_information_cellType,min_cell_number = 10)
  high_varGene_names <- Feature_selection_M3Drop(data_type_filtered$expression_profile)
  scLearn_model_learning_result<-scLearn_model_learning(high_varGene_names,data_type_filtered$expression_profile,data_type_filtered$sample_information_cellType,bootstrap_times=1)
  rawcounts2<-assay(test_singlecellnet[[i]])
  data_qc_query<-Cell_qc(rawcounts2,species="Mm")
  data_qc_query$expression_profile=rawcounts2
  scLearn_predict_result<-scLearn_cell_assignment(scLearn_model_learning_result,data_qc_query$expression_profile,diff=0.05,threshold_use=TRUE,vote_rate=0.6)
  y_pred=scLearn_predict_result$Predict_cell_type
  # Calculate accuracy and f1 scores
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_scLearn <- mean(accuracies)
mean_f1_score_Dolung_scLearn <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_scLearn * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_scLearn * 100, 2), "%\n")
#Mean accuracy: 76.66 %
#Mean F1 score: 71.39 %
#Calculate the unassigned rate
mean_unassigned_Dolung_scLearn <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_scLearn * 100, 2), "%\n")
#Mean Unssigned score: 15.25 %
#Returns a list of dataframes and category names with accuracy and f1 scores
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_scLearn <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_scLearn$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_scLearn$f1)

###################scmap######################
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  sce_1=test_singlecellnet[[i]]
  sce_2=train_singlecellnet[[i]]
  # Define the cell-type1 column
  colData(sce_2)$cell_type1 <- sce_2$True_Class1
  # Define feature_symbol column
  rowData(sce_2)$feature_symbol=rownames(sce_2)
  #Build the cell index
  refSCE <- scmap::selectFeatures(sce_2, suppress_plot = FALSE)
  # Single-cell-based index construction
  cell_ref <- scmap::indexCell(refSCE)
  rowData(sce_1)$feature_symbol <- rownames(sce_1)
  # Cell-type annotation based on nearest-neighbour algorithm
  nearest_neighbors <- scmap::scmapCell(
    projection = sce_1,
    index_list = list(immune1 = metadata(cell_ref)$scmap_cell_index),
    w = 10
  )
  # Get cell-type labels
  mode_label <- function(neighbors, metadata = sce_2$True_Class1){
    freq <- table(metadata[neighbors])
    label <- names(freq)[which(freq==max(freq))]
    if (length(label)>1) {
      return("ambiguous")
    }
    return(label)
  }
  cell_labs <- apply(nearest_neighbors$immune1$cells, 2, mode_label)
  y_pred=cell_labs
  #Calculate accuracy and f1 scores
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_scmap <- mean(accuracies)
mean_f1_score_Dolung_scmap <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_scmap * 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_scmap * 100, 2), "%\n")
mean_accuracy_Dolung_scmap =0.9314
mean_f1_score_Dolung_scmap=0.9418
#Mean accuracy: 93.14 %
#Mean F1 score: 94.18 %
#Calculate the unassigned rate
mean_unassigned_Dolung_scmap <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_scmap * 100, 2), "%\n")
#Mean Unssigned score: 0.99 %
#Returns a list of dataframes and category names with accuracy and f1 scores
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_scmap <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_scmap$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_scmap$f1)

#############scPred##################

all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  reference <- getFeatureSpace(train_seurat[[i]], "True_Class1")
  reference <- trainModel(reference)
  ctrl <- scPredict(test_seurat[[i]], reference)
  y_pred=ctrl@meta.data$scpred_prediction
  # Calculate average accuracy and average f1 score
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}

# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_scPred <- mean(accuracies)
mean_f1_score_Dolung_scPred <- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_scPred* 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_scPred * 100, 2), "%\n")
# Mean accuracy: 95.14 %
#Mean F1 score: 90.56 %
#Calculate the unassigned rate
mean_unassigned_Dolung_scPred <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_scPred * 100, 2), "%\n")
#Mean Unssigned score: 2.44 %
# Returns a list of dataframes and category names with accuracy and f1 scores
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_scPred <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
#Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_scPred$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_scPred$f1)

#######################seurat###############

all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  #Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  transfer.anchors <- FindTransferAnchors(reference =train_seurat[[i]], query = test_seurat[[i]], dims = 1:30)
  predictions <- TransferData(anchorset = transfer.anchors, refdata = y_train, dims = 1:30)
  y_pred=predictions$predicted.id
  # Calculate the unassigned rate
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_Seurat <- mean(accuracies)
mean_f1_score_Dolung_Seurat<- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_Seurat* 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_Seurat * 100, 2), "%\n")
#Mean accuracy: 96.5 %
#Mean F1 score: 92.5 %
#Calculate the unassigned rate
mean_unassigned_Dolung_Seurat <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_Seurat * 100, 2), "%\n")
#Mean Unssigned score: 0 %

# Returns a list of dataframes and category names with accuracy and f1
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_seurat <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_seurat$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_seurat$f1)

#####################singlecellnet#############
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
for (i in 1:length(folds)) {
  # Get training and test sets
  train_index <- unlist(folds[-i])
  test_index <- folds[[i]]
  X_train <- X[train_index, ]
  y_train <- y[train_index]
  X_test <- X[test_index, ]
  y_test <- y[test_index]
  X_train=t(X_train)
  X_test=t(X_test)
  ########Load training dataset###########
  stTM=train_seurat[[i]]@meta.data
  expTMraw=X_train
  stTM$cell=rownames(stTM)
  ########Load query dataset##########
  stQuery=test_seurat[[i]]@meta.data
  expQuery=X_test
  #If same species
  commonGenes<-intersect(rownames(expTMraw), rownames(expQuery))
  expTMraw = expTMraw[commonGenes,]
  set.seed(100)
  stList = splitCommon(sampTab=stTM,  dLevel="True_Class1")
  stTrain = stList[[1]]
  expTrain = expTMraw[,rownames(stTrain)]
  class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "True_Class1", colName_samp = "cell")
  crPBMC<-scn_predict(class_info[['cnProc']], expQuery, nrand=50)
  test = crPBMC[,colnames(crPBMC) %in% colnames(expQuery)]
  stQuery <- assign_cate(classRes = test, sampTab = stQuery, cThresh = 0.5)
  y_pred=stQuery$category
  # Calculate accuracy and f1 scores
  accuracy <- mean(y_pred == y_test)
  f1_score <- F1_Score(y_test, y_pred)
  accuracies <- c(accuracies, accuracy)
  f1_scores <- c(f1_scores, f1_score)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  cat("F1 score:", round(f1_score * 100, 2), "%\n")
  #Calculate the unassigned rate
  unassigned=0
  for (i in 1:length(y_pred)) {
    if (!(y_pred[i] %in% y_train) ) {
      unassigned=unassigned+1
    }
  }
  unassigned_rate=unassigned/length(y_test)
  unassigned_rates=c(unassigned_rates,unassigned_rate)
  cat("Uassigned_rate:", round(unassigned_rate * 100, 2), "%\n")
  #Save accuracy and f1 scores for each category
  all_metrics[[i]] <- calculate_metrics(y_test, y_pred)
}
# Calculate average accuracy and average f1 score
mean_accuracy_Dolung_SingleCellNet <- mean(accuracies)
mean_f1_score_Dolung_SingleCellNet<- mean(f1_scores)
cat("Mean accuracy:", round(mean_accuracy_Dolung_SingleCellNet* 100, 2), "%\n")
cat("Mean F1 score:", round(mean_f1_score_Dolung_SingleCellNet * 100, 2), "%\n")
#Mean accuracy: 92.89 %
#Mean F1 score: 97.78 %
#Calculate the unassigned rate
mean_unassigned_Dolung_SingleCellNet <- mean(unassigned_rates)
cat("Mean Unssigned score:", round(mean_unassigned_Dolung_SingleCellNet * 100, 2), "%\n")
#Mean Unssigned score: 1.55 %
# Save accuracy and f1 scores for all test sets to a list
accuracy_list <- c()
f1_list <- c()
row_names=c()
for (i in 1:length(all_metrics)) {
  metrics <- all_metrics[[i]]
  accuracy_list <- c(accuracy_list, metrics[[1]]$accuracy)
  f1_list <- c(f1_list, metrics[[1]]$f1)
  row_names=c(row_names,row.names(metrics[[1]]))
}
# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_singlecellnet <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_singlecellnet$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_singlecellnet$f1)


##########scTransformer##########
all_metrics <- list()
accuracies <- c()
f1_scores <- c()
unassigned_rates=c()
library(reticulate)
# Create numpy object
np <- import("numpy", convert = FALSE)
# Read .npy file
data<- np$load('/home/yuanjiaxin/ddd/Domingo-Gonzalez/log/y_predicts.npy')
#Convert data to R vector
result <- as.vector(data)
table(result)
vector1=unique(lung1@meta.data$True_Class)
vector1=as.character(vector1)
vector2 <- toupper(vector1)
vector2 =as.character(vector2 )
for (i in 1:length(result)) {
  index=which(vector2 == result[i])
  result[i]=vector1[index]
}

y_pred=result
table(y_pred)
###################
# Read .npy file
data<- np$load('/home/yuanjiaxin/ddd/Domingo-Gonzalez/log/y_tests.npy')
result <- as.vector(data)
table(result)
vector1=unique(lung1@meta.data$True_Class)
vector1=as.character(vector1)
vector2 <- toupper(vector1)
vector2 =as.character(vector2 )
for (i in 1:length(result)) {
  index=which(vector2 == result[i])
  result[i]=vector1[index]
}
y_test=result
table(y_test)
#####Calculate average accuracy and average f1 score
accuracy <- mean(y_pred == y_test)
f1_score <- F1_Score(y_test, y_pred)
print(accuracy)
print(f1_score)
#########Calculate the unassigned rate

unassigned=0
for (i in 1:length(y_pred)) {
  if (!(y_pred[i] %in% y_test) ) {
    unassigned=unassigned+1
  }
}
unassigned_rate=unassigned/length(y_test)
#Save accuracy and f1 scores for each category
all_metrics <- calculate_metrics(y_test, y_pred)
#Save accuracy and f1 scores for all test sets to a list
accuracy_list <- c()
f1_list <- c()
row_names=c()

metrics <- all_metrics[[1]]
accuracy_list <- c(accuracy_list, metrics$accuracy)
f1_list <- c(f1_list, metrics$f1)
row_names=c(row_names,row.names(metrics))

# Combine accuracy and f1 scores into one dataframe
accuracy_list[is.na(accuracy_list)]=0
f1_list[is.na(f1_list)]=0
df_all <- data.frame(accuracy = accuracy_list, f1 = f1_list, class =row_names )
# Group by category name and calculate average accuracy and average f1 score
df_avg_Dolung_scTransformer <- aggregate(df_all[, c("accuracy", "f1")], by = list(df_all$class), mean)
# Output average accuracy and average F1 score for each category
cat("Average accuracy by class:\n")
print(df_avg_Dolung_scTransformer$accuracy)
cat("Average F1 score by class:\n")
print(df_avg_Dolung_scTransformer$f1)







