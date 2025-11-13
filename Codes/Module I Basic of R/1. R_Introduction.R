# ===================================================================
#               AI and Biotechnology/Bioinformatics        
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Resrarch Internship (2025)
# -------------------------------------------------------------------
#               Module I: Getting Started with R (Class Ib)
# -------------------------------------------------------------------
# ==================================================================

# Topics:
# 1. Setting the working directory properly
# 2. Creating and organizing project folders
# 3. How R code works: functions, syntax, and execution
# 4. Variables and data types in R
# 6. Importing CSV files and working with categorical data
# 7. Saving scripts, outputs, and the R workspace

# Lecture Recording: https://youtu.be/YybbWfD_VjE?feature=shared
# GitHub: https://github.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025

# -------------------------------------------------------------------------------------------------------

#### Working Directory ####
# It is a main folder where R looks for files and save outputs/results

getwd()  # to check current working directory

# Set Working Directory #

# Before starting, set your working directory. 
# This helps R know where to find or save your files.

# NOTE: If you're copying a folder path from Windows,
# make sure to either:
# 1. Replace all backslashes "\" with forward slashes "/" 
# OR 
# 2. Use double backslashes "\\" 

# Example:
# If your folder path is: C:\Users\YourName\R_Projects\Biostats
# Use either of these in R:

setwd("C:/Users/YourName/Documents/Data_Analysis")
# or
setwd("C:\\Users\\YourName\\Documents\\Data_Analysi")

# --------------------------------------------------------------------------

#### Organize Project Folder ####

# First, create a main project folder:
# Go to the File menu > select 'New Project' > choose 'New Directory' > then 'New Project'
# Name your project and choose the location (your main folder)
# R will create a .Rproj file and set it as your working environment
# check the name of project on top right corner of RStudio

# Create common subfolders for project organization
dir.create("data")     # For storing raw or cleaned data files
dir.create("script")   # For saving R scripts
dir.create("results")  # For saving analysis outputs

# --------------------------------------------------------------------------

#### How R Code Works ####

# Functions in R follow the format: function_name(arguments)
print("Biotechnology")

# Create a simple numeric vector 
gene_expression <- c(2.3, 3.8, 3.9, 5.6, 9.4)

# "gene_expression" is the variable name
# "<-" is the assignment operator used to store values
# c() is the function to combine values into a vector/list


# Calculate the mean of the vector
mean(gene_expression)
mean_result <- mean(gene_expression)

# Basic visualizations
plot(gene_expression)     # Scatterplot
hist(gene_expression)     # Histogram
barplot(gene_expression)  # Bar plot

# Summary statistics
summary(gene_expression)

# --------------------------------------------------------------------------

#### R Data Types ####

# 1. Numeric (decimal numbers) 

x <- 24.5     
class(x) # to check the data type of any variable


y <- -40.4 
class(y)

# 2. Integer (whole numbers) 

z <- 24L      # Integer (note the 'L' suffix)
class(z)

# By default are treated both 'decimal' and 'integer' values as numeric unless 'L' is used

numeric <- c(2, 3, 4)
class(numeric)

integer <- c(2L, 3L, 4L)
class(integer)


# Example vector of numeric values
var_1 <- c(28.4, 23.7, 2.9)
class(var_1)

# Convert numeric to integer
var_2 <- as.integer(var_1)
class(var_2)

#R truncates the decimal part when converting numeric to integer.

# Convert integer back to numeric
var_3 <- as.numeric(var_2)
class(var_3)

# Integer uses less memory than numeric.

# 3. Character / String Data 
# Must be enclosed in single (' ') or double (" ") quotes

var_4 <- c("gene1", "gene2", "gene3")
class(var_4)

# Incorrect example (gene1 is unquoted and will throw error if not defined)
var_5 <- c(gene1, "gene2", "gene3")

# R consider unquoted strings as variable name

# Character vector with numeric values in quotes treated as text
var_6 <- c("3.8", "4.6", "5.9")
class(var_6)

mean(var_6)  # This will return an error since the values are characters

# Proper numeric vector
var_7 <- c(3.8, 4.6, 5.9)
class(var_7)

mean(var_7)  # Returns numeric mean

# 4. Factor (Categorical) Data 

# Categorical vector with string labels
disease_status <- c("cancer", "normal", "cancer", "cancer", "normal")
class(disease_status)

# Convert character vector to factor
disease_status <- as.factor(disease_status)
class(disease_status)
disease_status

# Set factor level order manually
disease_status <- factor(disease_status, 
                         levels = c("normal", "cancer"))

disease_status

# 5. Logical Data Type 

age <- c(23, 25, 28)
var_8 <- age < 25    
print(var_8   )# Returns logical vector (TRUE/FALSE)

# --------------------------------------------------------------------------
# Import Data from CSV

data <- read.csv(file.choose())

# View data in spreadsheet format
View(data)

# Check structure of the dataset
str(data)

# Convert 'height' column to factor
data$height_fac <- as.factor(data$height)
str(data)

# Reorder factor levels manually
data$height_fac <- factor(data$height_fac, 
                          levels = c("Tall", "Medium", "Short"))
levels(data$height_fac)

# Convert character/factor to numeric using factor() function
# Example # 0 = Short, 1 = Medium & 2 = Tall

data$height_num <- factor(data$height_fac,
                        levels = c("Short", "Medium", "Tall"),
                        labels = c(0, 1, 2)) 

# Convert 'gender' to factor
data$gender_fac <- as.factor(data$gender)
str(data)

# Convert factor to numeric using ifelse statement (Female = 1, Male = 0)
data$gender_num <- ifelse(data$gender_fac == "Female", 1, 0)
class(data$gender_num)

# Convert numeric gender code to factor
data$gender_num <- as.factor(data$gender_num)
class(data$gender_num)

#### Saving Outputs ####

# To save R script use keyboard shortcut 
# ctrl + s 
# Or go to the file menu & select 'Save As' option

# To open R script use keyboard shortcut 
# ctrl + Enter 
# Or go to the file menu & select 'Open File' option


# Save file  as CSV in results folder
write.csv(disease_status, file = "results/patient_data.csv")

# Save the entire R workspace
save.image(file = "full_workspace.RData")

# Save selected objects only
save(gene_expression, disease_status, file = "workspace.RData")

# Load workspace
load("workspace.RData")
load("full_workspace.RData")

# --------------------------------------------------------------------------
#### Tasks ####

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

# 2. Create Project Folder

# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

# ---------------------------------------------------------------------------
# 3. Download "patient_info.csv" dataset from GitHub repository

# load the dataset into your R environment

# Inspect the structure of the dataset using appropriate R functions

# Identify variables with incorrect or inconsistent data types.

# Convert variables to appropriate data types where needed

# Create a new variable for smoking status as a binary factor:
  
# 1 for "Yes", 0 for "No"

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
# Save your R script in your script folder with name "class_Ib"
# Upload "class_Ib" R script into your GitHub repository
