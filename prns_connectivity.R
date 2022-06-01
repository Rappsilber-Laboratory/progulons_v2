############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Check if proteins in progulons are significantly inter-connected
# Authors: E Fiagbedzi and G Kustatscher
# Date: January 2022
############################################################################################################# #

# Load the required libraries
library(data.table); library(treeClust)

# Set random state
set.seed(42) 


#### Load and prep input data ####

# Progulons
prns <- fread("output_files/Progulons.csv.gz")

## Load ProteomeHD
prohd <- fread("input_files/ProteomeHD_v1_1.csv.gz")              # Read in fast using data.table
prohd <- data.frame(prohd, row.names = "Majority_protein_IDs")    # Convert to data.frame; use protein IDs as rownames
prohd <- prohd[, grep("Ratio", colnames(prohd)) ]                 # Keep only columns with SILAC ratios

## Keep only proteins that were quantified in at least 30 experiments
feature_count <- apply(prohd, 1, function(x){ sum(!is.na(x)) })
prohd_ratios_min95 <- prohd[feature_count >= 30,] 


#### Calculate protein co-regulation network using treeClust ####

## Obtain treeClust distances
tc_distances <- treeClust.dist( prohd_ratios_min95, d.num = 2, verbose = TRUE)

# convert treeClust distance object into a long table format of similarities (= edges & weights)
edges <- as.data.table(reshape2::melt(as.matrix(tc_distances)))
edges <- edges[, .(Protein_1=as.character(Var1),
                   Protein_2=as.character(Var2), 
                   edge_weight=(1-value)) ]

# keep only unique protein associations
edges <- edges[Protein_1 > Protein_2]

# identify the number of pairs that constitute 0.5%
N_edges_0.5pc <- floor(edges[,.N] / 100 * 0.5)

# get the corresponding subset of pairs
edges <- edges[order(-edge_weight)][1:N_edges_0.5pc]

# print resulting edges
edges


#### Call progulons at different cut-offs and calculate a connectivity p-value ####

# initialize data.table object
connectivity_test <- data.table()

# set range of RF score cut-offs to screen
RF_range <- seq(0.5, 1.0, 0.01)

for(i in RF_range){
  
  for (prn_id in unique(prns[, Progulon_ID])){
  
  # select proteins of current progulon
  current_prn <- prns[ Progulon_ID == prn_id & Mean_RF_score > i, Protein_IDs ]
  
  # calculate number of possible combinations
  n_possible_combinations <- choose(length(current_prn), 2)
  
  # find number of observed combinations
  n_connections <- edges[Protein_1 %in% current_prn & Protein_2 %in% current_prn, .N]
  
  # find number of random combinations
  n_random_combinations <- round( n_possible_combinations / 100 * 0.5, digits = 0 )
  
  # create a confusion matrix for fischer exact test
  conf_matrix <- matrix(c(n_connections, n_random_combinations,
                          n_possible_combinations-n_connections,
                          n_possible_combinations-n_random_combinations),
                        nrow=2,
                        dimnames=list(c("measured", "random"),
                                      c("linked", "not_linked")))
  
  # calculate pValue for each progulon connectivity
  pValue = fisher.test(conf_matrix)$p.value
  
  # create temporary data.table to store values
  tmp_dt <- data.table(
    "Progulon_ID"=prn_id,
    "No_of_connections_observed"=n_connections,
    "Connectivity_p_value"=pValue,
    "RF_cut_off" = i)
  
  # combine temporary data with initialized table
  connectivity_test <- rbind(connectivity_test, tmp_dt)
  }
}


#### Clean-up and re-name progulon data files ####

# These progulons have significant connectivity at some RF score cut-off
sign_prns <- connectivity_test[ Connectivity_p_value < 0.05, .N , Progulon_ID ][ order(Progulon_ID) , Progulon_ID ]

# Read in progulons in long format
progulons_DT <- fread("output_files/Progulons.csv.gz")

# Keep only progulons with significant connectivity (i.e. remove PRN29)
progulons_DT <- progulons_DT[ Progulon_ID %in% sign_prns ]

# Create a table with old and new desired progulon names
cross_ref <- data.table( old_names = unique(progulons_DT$Progulon_ID),
                         new_names = paste("P",
                                           formatC(x = 1:length(unique(progulons_DT$Progulon_ID)), width = 2, format = "d", flag = "0"), 
                                           sep = ""))

# Re-name the progulons by merging
progulons_DT <- merge( progulons_DT, cross_ref, by.x = "Progulon_ID", by.y = "old_names")  # Add new names
progulons_DT[, Progulon_ID := new_names ]  # Assign new names
progulons_DT[, new_names := NULL ]         # Remove new name column again


# Read in progulons in wide format
progulons_wf <- fread("output_files/Progulon_Scores.csv.gz")

# Keep only progulons with significant connectivity (i.e. remove PRN29)
progulons_wf <- progulons_wf[, c("Protein_IDs", "Feature_count", "Protein_names", "Gene_names", sign_prns), with = FALSE]

# Re-name the progulons
setnames( progulons_wf, old = cross_ref$old_names, new = cross_ref$new_names )

# Overwrite progulon files with the newly named progulons
fwrite(progulons_DT, "output_files/Progulons.csv.gz")
fwrite(progulons_wf, "output_files/Progulon_Scores.csv.gz")

# Clean up the connectivity table and write to disk
connectivity_test <- connectivity_test[ Progulon_ID %in% sign_prns ] # Keep only progulons with significant connectivity (i.e. remove PRN29)
connectivity_test <- merge( connectivity_test, cross_ref, by.x = "Progulon_ID", by.y = "old_names")  # Add new names by merging
connectivity_test[, Progulon_ID := new_names ]  # Assign new names
connectivity_test[, new_names := NULL ]         # Remove new name column again

# Write data.table file to disk
fwrite(connectivity_test, "output_files/prn_connectivity.csv.gz")



