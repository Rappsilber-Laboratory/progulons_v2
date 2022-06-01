############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Create a set of random modules to test the progulonFinder
# Authors: G Kustatscher
# Date: April 2022
############################################################################################################# #

# Load the required libraries
library(data.table)

# Set random seed 
set.seed(42)

# Load ProHD
ProHD <- fread("input_files/ProteomeHD_v1_1.csv.gz") # loading ProteomeHD from NBT supplements

# Define the ratio columns
prohd_ratio_cols <- grep("^Ratio", names(ProHD), value = TRUE)

# Select relevant columns and turn into data.frame
ProHD <- data.frame( ProHD[, c("Majority_protein_IDs", prohd_ratio_cols), with = FALSE], row.names = "Majority_protein_IDs")

# Select relevant proteins with >95 features (same as in NBT paper)
count_features <- function(x){ sum( !is.na(x) ) }   # Function to calculate number of SILAC ratios per protein
feature_count <- apply(ProHD, 1, count_features)    # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 95 ,]              # Discard proteins detected in fewer than 95 experiments

# Calculate the sizes of the 72 real core modules
core_modules <- fread("output_files/core_modules.csv")
module_sizes <- core_modules[, sapply(.SD, function(x){ sum(x != "") }) ]

# Create 72 random groups of genes in sizes that match the actually tested core modules
random_modules <- data.table()
random_modules_names <- paste("RM", formatC( 1:length(module_sizes), width = 2, flag = "0"), sep = "_")
for(i in 1:length(module_sizes)){
  tmp_sample <- sample( rownames(ProHD), module_sizes[i] )   # The random sample
  tmp_sample <- tmp_sample[ 1:max(module_sizes)]             # Add NAs to extend to maximum length
  tmp_sample[ is.na(tmp_sample) ] <- ""                      # Replace NAs with ""
  random_modules[, random_modules_names[i] := tmp_sample ]   # Create column in data table with this entry
}

# Write out the result file
fwrite(random_modules, "output_files/random_modules.csv")

