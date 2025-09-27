# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================

# ---------------------------------------------------------------------
#              AI and Omics Research Internship (2025)
# ---------------------------------------------------------------------
#             Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------------------
#                     Microarray Data Analysis
# =====================================================================

# Topics covered in this script:
# 1. Quality Control (QC) 
# 2. RMA Normalization
# 3. Pre-processing and Filtering

#######################################################################
#### 0. Install and Load Required Packages ####
#######################################################################

# Bioconductor provides R packages for analyzing omics data (genomics, transcriptomics, proteomics etc).

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------

# Series matrix files are preprocessed text files containing 
# expression values, sample annotations, and probe information.
# Reason: Useful for quick exploratory analysis when raw CEL files are not needed.

gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)

# Extract expression data matrix (genes/probes × samples)
# Rows corresponds to probes and columns corresponds to samples
expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)


# Extract feature (probe annotation) data
# Rows corresponds to probes and columns corresponds to samples
feature_data <-  fData(gse_data$GSE79973_series_matrix.txt.gz)


# Extract phenotype (sample metadata) data
# Rows corresponds to samples and columns corresponds to probes
phenotype_data <-  pData(gse_data$GSE79973_series_matrix.txt.gz)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1)) 

# --------------------------------------
#### Download Raw Data (CEL files) ####
# --------------------------------------

# CEL files contain raw probe-level intensity values for affymetrix platform.
# Raw data required full preprocessing (e.g., RMA normalization, QC)

# CEL files are large, and downloads may fail even with a good connection. 
#It's recommended to download raw data directly from NCBI GEO

# skip this step if you already downloaded data from NCBI

# Fetch GEO supplementry files
getGEOSuppFiles("GSE79973", baseDir = "Raw_Data", makeDirectory = TRUE)

# Important Note: 
# For Affymetrix microarray data, the preprocessing pipeline is the same 
# whether raw CEL files are downloaded from NCBI GEO or ArrayExpress. 

# (See tutorial for detailed explanation of this step: https://youtu.be/DZMxkHxwWag?feature=shared) 

# Untar CEL files if compressed as .tar
untar("Raw_Data/GSE79973_RAW.tar", exdir = "Raw_Data/CEL_Files")

# Alternatively, unzip if data is compressed as .zip
unzip("Raw_Data/E-GEOD-79973.zip", exdir = "Raw_Data/E_GEOD79973")

# Read CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

raw_data   # Displays basic information about the dataset
# Note down the annotation (e.g annotation= hgu133plus2) from this output.
# You will need this in the next step to select and install the correct
# annotation package (e.g., hgu133plus2.db) for mapping probe IDs to genes.

# ---------------------------------------------------
#### Quality Control (QC) Before Pre-processing ####
# ---------------------------------------------------

# QC identifies outlier arrays, hybridization problems, or technical biases.
# arrayQualityMetrics: # This package generates automated QC reports for microarray data.
# It applies multiple complementary methods to detect technical issues:
#   - Boxplots and density plots: check distribution of intensities 
#   - MA-plots: visualize systematic biases between arrays 
#   - Heatmaps and distance matrices: identify clustering/outliers
#   - PCA: detect unusual variation/detecting outliers or batch effects
#
# The output is an interactive HTML report (index.html file) summarizing QC results.

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

# -------------------------------------------------------
#### RMA (Robust Multi-array Average) Normalization ####
# -------------------------------------------------------

# RMA is a popular method for normalizing Affymetrix microarray data by:
# 1. Background correcting, 
# 2. normalizing probe intensities using quantile normalization and 
# 3. summarizing them into gene-level expression values using a robust median polish algorithm.

# This method reduces experimental variation across multiple arrays, 
# producing more symmetrical and reliable normalized expression data 
# compared to other approaches

normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples

# ---------------------------------------------------------------------------
#### Filter Low-Variance Transcripts (“soft” intensity based filtering) ####
# ---------------------------------------------------------------------------

# Filtering removes probes with low or uninformative expression signals.
# Reason: Reduces noise and improves statistical power in differential expression & Machine Learning.


# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 

# -----------------------------------
#### Phenotype Data Preparation ####
# -----------------------------------

# Phenotype data contains sample-level metadata such as condition, 
# tissue type, or disease status.
# Required to define experimental groups for statistical analysis.

class(phenotype_data$source_name_ch1) 

# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)

#-------------
# Assignment 
#-------------

# Work through the full preprocessing workflow with your own dataset

# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

# Submission Instructions:

# Upload your R script implementing the workflow to GitHub and 
# provide the repository link in the form.

# Google Form: https://forms.gle/1e9tj2Mqf5T9FKEJ7

# Deadline: Sunday 5th October, 2025 (Midnight)


