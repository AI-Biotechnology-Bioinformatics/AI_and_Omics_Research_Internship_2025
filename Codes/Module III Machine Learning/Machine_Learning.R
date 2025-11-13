# =====================================================================
#                  AI and Biotechnology / Bioinformatics
# =====================================================================
#               AI and Omics Research Internship (2025)
# =====================================================================
#        Module III: Introduction to Machine Learning for Genomics 
# =====================================================================
#                      Microarray Data Analysis
# =====================================================================

# Topics Covered:
#   1. Preprocessing microarray data
#   2. Feature selection (Boruta, RFE)
#   3. Model training (RF, SVM, ANN)
#   4. Model evaluation (Confusion Matrix, AUC, ROC)

# --------------------------
# Load Required Packages
# --------------------------
# Each library below supports a specific step in our ML workflow
library(caret)         # End-to-end ML workflow (preprocessing, training, validation)
library(RANN)          # KNN-based imputation of missing values
library(randomForest)  # Random Forest algorithm
library(Boruta)        # Wrapper for feature selection
library(kernlab)       # SVM algorithms
library(pROC)          # ROC curve and AUC calculation
library(ggplot2)       # Visualization (for plots, bar charts, ROC etc.)

# ---------------------------
#### DATA PREPROCESSING ####
# ---------------------------

# In this stage, we prepare our data for modeling by cleaning, 
# transforming, and standardizing it.

data <- read.delim("Batch_Corrected_Processed_Matrix.txt")  # Input = normalized expression matrix

dim(data)  # Number of genes (columns) and samples (rows)
str(data)  # Data structure overview
View(data) # Open dataset for inspection in RStudio

# -------------------------------------
# ---- Step 1. Data Transformation ----
# --------------------------------------
# Many expression datasets are right-skewed (some genes have very large values).
# Log transformation helps reduce skewness and stabilize variance.

boxplot(data[1:30],
        outline = FALSE,
        col = "skyblue",
        main = "Before Log Transformation")

par(mfrow = c(1, 2))
hist(data[, 18], main = "Before Log Transform", col = "skyblue")
hist(log10(data + 1)[, 18], main = "After Log Transform", col = "coral")

data_log <- log10(data + 1)  # +1 avoids log(0) which is undefined

# --------------------------------------------------------------------
#  ---- Step 2. Transpose Data (Samples = Rows, Genes = Columns)  ----
# ---------------------------------------------------------------------
# ML algorithms expect data in the format: rows = samples, columns = features.
data_t <- as.data.frame(t(data_log))

# -------------------------------------
#  ---- Step 3. Feature Filtering  ----
# -------------------------------------
# Some genes show almost no variation across samples — these carry little information.
# Removing them reduces computational load and improves learning.
nzv <- preProcess(data_t, method = "nzv", uniqueCut = 15)
data_t <- predict(nzv, data_t)

# Optionally, you can filter by variance or correlation.
# to keep only the top variable genes (e.g., top 500 or 1000) that avoid redundancy.

# ----------------------------------------
#  ---- Step 4. Center and Scale Features
# ----------------------------------------
# Centering subtracts the mean; scaling divides by standard deviation.
# Ensures all genes contribute equally to the model (important for SVM/ANN).
process <- preProcess(data_t, method = c("center", "scale"))
data_t <- predict(process, data_t)

# -------------------------------------------
#  ---- Step 5. Handle Missing Values  ----
# -------------------------------------------
# Missing values are common in omics data.
# We either remove or impute them.
anyNA(data_t)
sum(is.na(data_t))

# Option 1: Remove missing values — simple but can lose data.
# Option 2: Impute missing values (recommended).
# KNN imputation replaces missing values using the mean of the k nearest samples.

knnImpute <- preProcess(data_t, method="knnImpute")
data_t <- predict(knnImpute, data_t)
sum(is.na(data_t))  # Check again to confirm all missing values are handled.

# -------------------------------------
#  ---- Step 6. Add Class Labels  ----
# -------------------------------------
# Machine Learning requires both input (X) and output (Y).
# Here, Y = class labels (e.g., cancer vs normal).
groups <- factor(brief$Characteristics.sampling.site.,
                 levels = c("cancer", "normal"),
                 labels = c(1, 0))

class(groups)
levels(groups)
table(groups)

data_t <- cbind(groups, data_t)
View(data_t)

#--------------------------------------
# ---- Step 7. DATA PARTITIONING  ----
#--------------------------------------

# Split data into Training (70%) and Testing (30%) for unbiased evaluation.
set.seed(123)
trainIndex <- createDataPartition(data_t$groups, p = 0.7, list = FALSE)
train_data <- data_t[trainIndex, ]
test_data  <- data_t[-trainIndex, ]

x_train <- train_data[, -1]  # Gene expression matrix
y_train <- train_data$groups # Class labels

x_test  <- test_data[, -1]
y_test  <- test_data$groups

#-----------------------------------------------------------------------------------------

#----------------------------#
#### FEATURE SELECTION ####
#----------------------------#

# The goal is to reduce thousands of genes to a small informative subset
# that carries the most predictive signal for the classification task.

#--------------------------------------------------
#### Method 1: Feature Selection with Boruta ####
#--------------------------------------------------

# Boruta is a wrapper around Random Forest that tests each feature’s importance.
# It keeps only those that significantly contribute to the prediction.

boruta_result <- Boruta(x = x_train, 
                        y = y_train, 
                        doTrace = 1)

final_boruta <- TentativeRoughFix(boruta_result) # fix tentative features

boruta_features <- getSelectedAttributes(final_boruta)

write.csv(boruta_features, "Results/boruta_selected_genes.csv", row.names = FALSE)

#------------------------------------------------------------------------------
#### Method 2: Feature Selection with RFE (Recursive Feature Elimination) ####
#-------------------------------------------------------------------------------
# RFE iteratively removes least important features and retrains the model.

rfe_control <- rfeControl(functions = rfFuncs, #random forest
                          method = "cv", 
                          number = 3, 
                          verbose = TRUE)

rfe_result <- rfe(x = x_train, 
                  y = y_train, 
                  sizes = c(1:5, 10, 50), 
                  rfeControl = rfe_control)

plot(rfe_result, type = c("g", "o"))

rfe_features <- predictors(rfe_result)

write.csv(rfe_features, "Results/rfe_selected_genes.csv", row.names = FALSE)

# Common features selected by both Boruta and RFE

common_features <- intersect(boruta_features, rfe_features)

write.csv(common_features, "Results/common_selected_genes.csv", row.names = FALSE)

#--------------------------------------
#### MODEL TRAINING & EVALUATION ####
#--------------------------------------

# We’ll compare three classifiers: Random Forest, SVM, and ANN.
# Each learns patterns differently.

set.seed(123)
ctrl <- trainControl(method = "cv", number = 10, verboseIter = TRUE)

                        ## ---- BORUTA Feature Set ----

train_boruta <- x_train[, boruta_features]
test_boruta  <- x_test[, boruta_features]

# ----------------------------------
# ---- Model 1: Random Forest  ----
# ----------------------------------
RF_boruta <- train(x = train_boruta, y = y_train, method = "rf", importance = TRUE, trControl = ctrl)
plot(varImp(RF_boruta), top = 10)  # Visualize top contributing genes

# -----------------------------------------
#  ---- Model 2: SVM (Radial Kernel)  ----
# -----------------------------------------
svm_boruta <- train(x = train_boruta, y = y_train, method = "svmRadial",
                    trControl = ctrl,
                    tuneGrid= data.frame(C=c(0.25,0.5,1), sigma= 0.5),
                    prob.model = TRUE)

# --------------------------
#  ---- Model 3: ANN  ----
# --------------------------
ann_boruta <- train(x = train_boruta, y = y_train, method = "nnet",
                    trControl = ctrl,
                    tuneGrid = data.frame(size = 1:2, decay = 0),
                    MaxNWts = 2000)

# -----------------------------
#  ---- Model Evaluation  ----
# -----------------------------

# Prediction
rf_pred_boruta  <- predict(RF_boruta, newdata = test_boruta)
svm_pred_boruta <- predict(svm_boruta, newdata = test_boruta)
ann_pred_boruta <- predict(ann_boruta, newdata = test_boruta)

# Confusion Matrix shows model’s classification performance
rf_conf_boruta  <- confusionMatrix(rf_pred_boruta, y_test)
svm_conf_boruta <- confusionMatrix(svm_pred_boruta, y_test)
ann_conf_boruta <- confusionMatrix(ann_pred_boruta, y_test)

# Compare accuracy
results_boruta <- data.frame(
  Model = c("Random Forest", "SVM", "ANN"),
  Accuracy = c(rf_conf_boruta$overall["Accuracy"],
               svm_conf_boruta$overall["Accuracy"],
               ann_conf_boruta$overall["Accuracy"])
)
print(results_boruta)

# training error
# 1 - accuracy(mode's accuracy)

# ROC & AUC Curves
# ROC curve = True Positive Rate vs False Positive Rate
# AUC = area under ROC; measures model discrimination

rf_prob_boruta  <- predict(RF_boruta, newdata = test_boruta, type = "prob")
svm_prob_boruta <- predict(svm_boruta, newdata = test_boruta, type = "prob")
ann_prob_boruta <- predict(ann_boruta, newdata = test_boruta, type = "prob")

rf_roc_boruta  <- roc(y_test, as.numeric(rf_prob_boruta[, 1]))
svm_roc_boruta <- roc(y_test, as.numeric(svm_prob_boruta[, 1]))
ann_roc_boruta <- roc(y_test, as.numeric(ann_prob_boruta[, 1]))

rf_auc_boruta; svm_auc_boruta; ann_auc_boruta

# Plot ROC curves for visual comparison
par(mfrow = c(1, 3))
plot(rf_roc_boruta, col = "black", lwd = 2, main = paste("RF (AUC =", round(rf_auc_boruta, 3), ")"))
plot(svm_roc_boruta, col = "green", lwd = 2, main = paste("SVM (AUC =", round(svm_auc_boruta, 3), ")"))
plot(ann_roc_boruta, col = "blue", lwd = 2, main = paste("ANN (AUC =", round(ann_auc_boruta, 3), ")"))

                    
                    ## ---- RFE Feature Set ----


train_rfe <- x_train[, rfe_features]
test_rfe  <- x_test[, rfe_features]

# Train again using RFE features to compare results
RF_rfe  <- train(x = train_rfe,  y = y_train, method = "rf", trControl = ctrl)
svm_rfe <- train(x = train_rfe, y = y_train, method = "svmRadial", trControl = ctrl)
ann_rfe <- train(x = train_rfe, y = y_train, method = "nnet", trControl = ctrl)

# Evaluate and compare
rf_pred_rfe  <- predict(RF_rfe, newdata = test_rfe)
svm_pred_rfe <- predict(svm_rfe, newdata = test_rfe)
ann_pred_rfe <- predict(ann_rfe, newdata = test_rfe)

rf_conf_rfe  <- confusionMatrix(rf_pred_rfe, y_test)
svm_conf_rfe <- confusionMatrix(svm_pred_rfe, y_test)
ann_conf_rfe <- confusionMatrix(ann_pred_rfe, y_test)

results_rfe <- data.frame(
  Model = c("Random Forest", "SVM", "ANN"),
  Accuracy = c(rf_conf_rfe$overall["Accuracy"],
               svm_conf_rfe$overall["Accuracy"],
               ann_conf_rfe$overall["Accuracy"])
)
print(results_rfe)

# AUC comparison plot
rf_prob_rfe  <- predict(RF_rfe, newdata = test_rfe, type = "prob")
svm_prob_rfe <- predict(svm_rfe, newdata = test_rfe, type = "prob")
ann_prob_rfe <- predict(ann_rfe, newdata = test_rfe, type = "prob")

rf_roc_rfe  <- roc(y_test, as.numeric(rf_prob_rfe[, 1]))
svm_roc_rfe <- roc(y_test, as.numeric(svm_prob_rfe[, 1]))
ann_roc_rfe <- roc(y_test, as.numeric(ann_prob_rfe[, 1]))

plot(rf_roc_rfe, col = "blue", lwd = 2, main = "ROC Curve - RFE Features")
plot(svm_roc_rfe, col = "red", lwd = 2, add = TRUE)
plot(ann_roc_rfe, col = "green", lwd = 2, add = TRUE)

legend("bottomright",
       legend = c(
         paste("RF (AUC =", round(auc(rf_roc_rfe), 3), ")"),
         paste("SVM (AUC =", round(auc(svm_roc_rfe), 3), ")"),
         paste("ANN (AUC =", round(auc(ann_roc_rfe), 3), ")")
       ),
       col = c("blue", "red", "green"), lwd = 2, bty = "n")
