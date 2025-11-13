# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
# ===================================================================

# --------------------
# Topics in this module
# --------------------
#   1. Operators in R
#   2. Data Structures in R
#   3. User-defined Functions in R
#   4. Automating Workflows with for-Loops
# --------------------------------------------------------------------------------------------------

#### 1. Operators in R ####

# Operators are special symbols in R.
# They instruct R to perform actions on values or variables.
# You will use them for assignment, calculations, comparisons, and logical tests.

# ----------------------
# Assignment Operators
# ----------------------
# Used to store values inside variables

# <-   (Most common, assigns value on the right to the name on the left)
height <- c(1.75, 1.76, 1.82, 1.67)

# ->   (Same but assigns left-hand values to right-hand variable name)
c(68, 78, 85, 75) -> weight

# =    (Also assigns values, often used in function arguments)
smoking_status = c("Yes", "No", "No", "Yes")


# ----------------------
# Arithmetic Operators
# ----------------------
# Perform basic math: +, -, *, /, ^
#   +  addition
#   -  subtraction
#   *  multiplication
#   /  division
#   ^  exponent (to the power of)

# Example: calculate BMI (Body Mass Index = weight / height^2)
BMI <- weight/(height^2)
BMI

# Note: R applies operations element-wise when variables are vectors.
# This is called "vectorization", every weight is divided by every height squared.


# ---------------------
# Comparison Operators
# ---------------------
# Comparison operators ask TRUE/FALSE questions about values.
# They do not calculate a number, they return a logical output.

BMI > 25                    # Is BMI greater than 25?
BMI < 18.5                  # Is BMI less than 18.5?
height >= 1.75              # IS height greater than or equal to 1.75
weight <= 65                # Is weight less than or equal to 65
smoking_status == "No"      # Is smoking status equal to No
smoking_status != "No"      # Is smoking status not equal to No


# -------------------
# Logical Operators
# -------------------
# Combine multiple conditions using:
#   &   AND
#   |   OR
#   !   NOT

# Example: is the patient overweight AND a smoker?
(BMI > 25) & (smoking_status == "Yes")

# Example: overweight OR smoker?
(BMI > 25) | (smoking_status == "Yes")

# Example: NOT a smoker?
!(smoking_status == "Yes")

# ----------------------
# Summary of Operators
# ----------------------
#   - Assignment stores values
#   - Arithmetic performs calculations
#   - Comparison returns TRUE/FALSE
#   - Logical combines conditions
# --------------------------------------------------------------------------------------------------

# ------------------------
# 2. Data Structures in R
# ------------------------
# Data structures are how R organizes information.
# We commonly use:
#   1. Vectors
#   2. Lists
#   3. Matrices
#   4. Data Frames

# ---------
# Vectors
# ---------
# Simplest structure in R. Stores a sequence of values of the SAME type.

num_vec <- c(1, 2, 3, 4)   # numeric vector
chrc_vector <- c("gene1", "gene2", "gene3")   # character vector
logical_vector <- c(TRUE, FALSE, TRUE)        # logical vector

# Important:
# Maybe R does not throw an error if you combine mixed types.
# Instead, it coerces all values into a single type.
mix_vector <- c("gene1", 1, "gene2", 2)
mix_vector  # all values become characters

# Indexing: extract elements with []
num_vec[2]       # second element
num_vec[2:4]     # elements 2 through 4

# ---------
# Lists
# ---------
# Unlike vectors, lists can store different types together.
all_vectors <- list(num_vec, chrc_vector, logical_vector)
all_vectors[[2]]   # extract the second element of the list

# ---------
# Matrices
# ---------
# A 2D table where all values must be the same type.

my_matrix <- matrix(1:9, nrow = 3, ncol = 3)  # fills column-wise by default
my_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE) # row-wise

# Access with [row, column]
my_matrix[2,3]
my_matrix[2,]

# ---------
# Data Frames
# ---------
# Most important structure for biological datasets.
# Columns can have different types (numeric, character, factor).
data <- data.frame(
  patient_id = c("P1", "P2", "P3"),
  age = c(65, 78, NA),
  diagnosis = c("cancer", "diabetes", "cancer")
)

# --------------------
# Dataset Assessment
# --------------------
# Before analyzing, inspect your dataset to understand its structure.
str(data)       # structure (column names, data types, preview)
head(data)      # first 6 rows
tail(data, 2)   # last 2 rows
dim(data)       # number of rows and columns
names(data)     # column names

# Access data using:
data$patient_id     # extract single column
data[1:2, c(1,3)]   # extract specific rows and columns

# Create new columns:
data$new_column <- c(1, 2, 3)

# ----------------
# Missing Values
# ----------------
# Real data often contains missing values (NA).
# You must check and handle them before analysis.

is.na(data)                # identify missing values
sum(is.na(data))           # total missing values
colSums(is.na(data))       # missing values per column
rowSums(is.na(data))       # missing values per row

# Ways to handle NA:

# remove rows with NA
clean_data1 <- na.omit(data)   


# remove columns with NA
clean_data_2 <- data[, colSums(is.na(data))==0]

# replace NA value with 0
clean_data_3 <- data
clean_data_3[is.na(clean_data_3)] <- 0

# replace NA value with mean
clean_data_4 <- data
clean_data_4[is.na(clean_data_4)] <- mean(data$age, na.rm = TRUE)

# ------------------------------
# Summary of Data Structures:
# ------------------------------

#   - Vectors: simple sequences of same data type
#   - Lists: mix of different data types
#   - Matrices: numeric tables
#   - Data Frames: mixed-type tables 


#--------------------
# 3. Functions in R
#--------------------
# Functions let us wrap code into reusable blocks.

# function is  reusable block of code 
# Why use functions?
#   - Avoid repetition
#   - Organize and simplify code
#   - Reuse across projects (save it for later use)
#   - Share with others

# A function in R has 4 key parts:
#   1. Name         -> the name you give to the function
#   2. Arguments    -> the inputs you provide to the function
#   3. Body         -> the set of operations the function performs
#   4. Return Value -> the output the function gives back

# Example: A function to calculate Body Mass Index (BMI)

# 1. Function Name: calculate_BMI
# 2. Arguments: weight (in kg), height (in meters)
# 3. Body: performs BMI calculation e.g   # Formula: BMI = weight / (height^2)
# 4. Return Value: the BMI value

calculate_BMI <- function(weight, height) {
  # Perform the BMI calculation
  Bmi <- weight / (height ^ 2)
  
  # Return the BMI value
  return(Bmi)
}

# Call the function by naming arguments explicitly
calculate_BMI(weight = 60, height = 1.75)

# Call the function using variables as arguments
calculate_BMI(weight = weight, height = height)

# If a function expects two arguments, you must provide both
# This will give an error because 'height' is missing
calculate_BMI(60) 

# You can assign default values to function arguments
calculate_BMI <- function(weight, height = 1.75) {
  # Perform the BMI calculation
  Bmi <- weight / (height ^ 2)
  
  # Return the BMI value
  return(Bmi)
}

# In this case, if you don’t provide height, R automatically uses 1.75 as the default.
calculate_BMI(weight = 60)

# ----------------------------
# Lazy evaluation in R
# ----------------------------
# If your function has three arguments, but the body only uses two,
# R does not force you to supply the third argument for the calculation.
# Example: 'age' is defined as an argument, but not used in the formula.

# Define the function with three arguments
calculate_BMI <- function(weight, height, age){
  # Perform the BMI calculation
  bmi <- weight / (height^2)
  
  # Return the BMI value
  return(bmi)
}

# Here we pass only 'weight' and 'height'
# Even though 'age' exists as an argument, it is ignored because it is not used
calculate_BMI(60, 1.65)

#---------
# Summary:
#---------
# Functions help us package logic once and apply it to different inputs.

# ----------------------------------
# 4. Automating Workflows with for-Loop
# ----------------------------------
# Suppose you have multiple datasets and you want to:
#   - import them,
#   - check missing values,
#   - clean columns,
#   - compute BMI,
#   - and save results.
#
# Instead of repeating steps for each file, we use loops.

# -----------------------
# Typical loop workflow:
# -----------------------

# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).
# Specify the list of files that need to be processed.
# Prepare an empty list in R to store results for later use.
# For each file in the list:
#          Import the data into R.
#          Check and handle missing values (NA).
#          Calculate BMI using the calculate_BMI function.
#          Save the processed results both as a CSV file and in the R results list.

# ----------------------------------------------------
# Calculate BMI of two dataset within loop

# Define input and output folders
input_dir <- "Raw_Data" 
output_dir <- "Results"


# create output folder if not already exist

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}


# List which files to process
files_to_process <- c("data_1.csv", "data_2.csv") 
# These must match exactly with the names in your working folder,
# otherwise R will not find them.

# For practice, you can instead use the provided datasets:
# BMI_data_1.csv
# BMI_data_2.csv
# (download from the GitHub repository).
https://github.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/tree/main

# Prepare empty list to store results in R 
result_list <- list()

# For each file with in a loop:
#       - import data
#       - handle NA values
#       - calculate BMI using calculate_BMI function
#       - save results (both as CSV and inside R list)


for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # handling missing values
  
  if("height" %in% names(data)){
    missing_count <- sum(is.na(data$height))
    
    cat("Missing values in 'height':", missing_count, "\n")
    data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
  }
  
  if("weight" %in% names(data)){
    missing_count <- sum(is.na(data$weight))
    
    cat("Missing values in 'weight':", missing_count, "\n")
    data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
  }
  # calculate BMI
  data$bmi <- calculate_BMI(data$weight, data$height)
  cat("BMI has been calculated successfully.\n")
  
  # save results in R
  result_list[[file_names]] <- data 
  
  # save results in Results folder
  output_file_path <- file.path(output_dir, paste0("BMI_results", file_names))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
}

# The loop repeats until all files are processed.

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

# --------
# Summary:
# --------
# Loops automate repetitive work 
# making your workflow faster 
# consistent, and reproducible

# --------------------------
# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv

# Each file contains three columns: 
# Gene_Id	
# padj	
# logFC


