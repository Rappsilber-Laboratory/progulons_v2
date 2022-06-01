############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Identify the core co-regulation modules in ProteomeHD by unsupervised clustering
# Authors: E Rullman & G Kustatscher
# Date: January 2022
############################################################################################################# #

# Load the required libraries
library(data.table); library(treeClust); library(dbscan); library(reshape2);library(ggplot2)

# Set random seed 
set.seed(42)


#### Calculate protein dissimilarities for ProteomeHD ####

# Load ProHD
ProHD <- fread("input_files/ProteomeHD_v1_1.csv.gz") # loading ProteomeHD from NBT supplements
ProHD2 <- ProHD[,1:4] # get geneID annotation

# Define the ratio columns
prohd_ratio_cols <- grep("^Ratio", names(ProHD), value = TRUE)

# Select relevant columns and turn into data.frame
ProHD <- data.frame( ProHD[, c("Majority_protein_IDs", prohd_ratio_cols), with = FALSE], row.names = "Majority_protein_IDs")

# Select relevant proteins with >95 features (same as in NBT paper)
count_features <- function(x){ sum( !is.na(x) ) }   # Function to calculate number of SILAC ratios per protein
feature_count <- apply(ProHD, 1, count_features)    # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 95 ,]              # Discard proteins detected in fewer than 95 experiments

# Use treeClust to learn dissimilarities
tc_dist   <- treeClust.dist(ProHD, d.num = 2, verbose = TRUE)

# Write out the top 0.5% of pairwise interactions for ClusterONE clustering (which is done in the command line using a Java tool)
edges <- as.data.table( reshape2::melt( as.matrix( tc_dist )))
edges <- edges[,.( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), edge_weight = (1-value)) ]
edges <- edges[ Protein_1 > Protein_2 ]                    # Removes duplicate pairs (incl. self-comparisons)

N_edges_0.5pc <- floor( edges[,.N] / 100 * 0.5 )              # Identify the number of pairs that constitute 0.5% 
edges   <-   edges[ order(-edge_weight) ][ 1:N_edges_0.5pc ]  # Get the corresponding subset of pairs
fwrite(edges, "output_files/edges_ProHD.txt", sep = "\t", col.names = FALSE)  # Write ProHD edges for clusterONE


#### ClusterONE clustering ####

# ClusterONE is available as a Java application from: https://paccanarolab.org/cluster-one.
# This was run in console with default setting using:
# java -jar cluster_one-1.0.jar edges_ProHD.txt -F csv > ProHD_clusterOne_output.csv

# Load the ClusterONE clusters back into R
clONE   <- fread("output_files/ProHD_clusterOne_output.csv")

# Turn into long format
clONE_long <- data.table()
for(i in 1:nrow(clONE)){ 
  tmp <- data.table( Protein = unlist( strsplit( clONE[i, Members] , split = " " )), clONE = i )
  clONE_long <- rbind(clONE_long, tmp)
}

# Prepare clONE cluster for intersection with OPTICS
clONE_long <- merge(clONE_long, ProHD2[,.(Protein = Majority_protein_IDs, Protein_names)])


#### OPTICS clustering ####

# Calculate OPTICS clustering order
OPTICS   <- optics(tc_dist)               

# Extract OPTICS clusters
xi_setting <- 0.0001 
OPTICS_xi   <- extractXi(OPTICS,   xi = xi_setting)

# Prepare output table
OPTICS_xi_clusters   <- data.table( Protein = rownames(ProHD),   optCl = OPTICS_xi$cluster )


#### Intersect all OPTICS cluster with all clONE clusters ####

Overlap_table <- data.frame(OPTICS_Cluster=numeric(),OPTICS_cl_size=numeric(),clONE_Cluster=numeric(),clONE_cl_size=numeric(),Overlap=integer(),Percentage_OPTICS = integer(), Percentage_clONE = integer(), listID = integer()) # create empty dataframe

OPTICS_cl <- sort(unique(OPTICS_xi_clusters$optCl))
clONE_cl <- sort(unique(clONE_long$clONE))

IDlist <- vector("list", length(OPTICS_cl)*length(clONE_cl))
ID <- 1

pb <- txtProgressBar(min = 0, max = length(OPTICS_cl), style = 3)
for (x in OPTICS_cl) {
  for (y in clONE_cl) {
    a <- OPTICS_xi_clusters$Protein[OPTICS_xi_clusters$optCl==x] # all IDs belonging to one cluster in subset 1
    b <- clONE_long$Protein[clONE_long$clONE==y] # all IDs belonging to one cluster in subset 2
    overlap <- intersect(a,b) # intersect all three at the same time
    perc_OPT <- length(overlap)/length(a)
    perc_clONE <- length(overlap)/length(b)
    Overlap_table <- rbind(Overlap_table,
                           data.frame(OPTICS_Cluster=x,OPTICS_cl_size=length(a),clONE_Cluster=y,clONE_cl_size=length(b),Overlap=length(overlap), Percentage_OPTICS = perc_OPT, Percentage_clONE = perc_clONE, listID = ID)) # create table with info about the clusters and their overlap
    IDlist[[ID]] <- overlap
    ID <- ID+1
  }
  setTxtProgressBar(pb, x)
}

Overlap_table <- data.table(Overlap_table) # Turn into data.table

# Extract the best-overlapping combinations
Best_overlap <- Overlap_table[0,] # create an empty table
OPTICS_cl <- sort(unique(Overlap_table$OPTICS_Cluster)) # create a vector for the loop of all optics clusters (as clONE is fuzzy)

for (i in OPTICS_cl) {
  subtable <- Overlap_table[Overlap_table$OPTICS_Cluster == i]
  best_row <- subtable[which.max(subtable$Overlap),]
  Best_overlap <- rbind(Best_overlap,best_row)
} 

# Create a longtable formatted matrix containing the overlapping IDs -> CORES
Cluster_IDs <- Best_overlap$listID
max_overlap <- max(Best_overlap$Overlap)
Cluster_prot_dt <- data.frame(Number = c(1:max_overlap))
annot <- ProHD2[,.(Majority_protein_IDs,Protein_names)]

for (y in Cluster_IDs) {
  proteins <- IDlist[[y]]
  NA_vec <- rep(NA,max_overlap-length(proteins))
  dt_col <- as.data.frame(c(proteins,NA_vec))
  dt_col <- merge(dt_col,annot, by.x = "c(proteins, NA_vec)", by.y = "Majority_protein_IDs", all.x = TRUE)
  clustname <- paste("Cluster",Best_overlap$OPTICS_Cluster[Best_overlap$listID==y],sep = "_")
  colnames(dt_col) <- c(clustname, paste(clustname,"Proteins", sep = "_"))
  Cluster_prot_dt <- cbind(Cluster_prot_dt,dt_col)
}

Cluster_prot_dt <- data.table(Cluster_prot_dt)

# Filter by minimum size of seed protein groups
min_size <- 4 # define minimum size of core
na_max <- nrow(Cluster_prot_dt) - min_size
final_cores <- Cluster_prot_dt[,colSums(is.na(Cluster_prot_dt))<=na_max, with=FALSE]


#### Write out resulting core modules ####

# Keep only protein IDs
final_cores[, Number := NULL ]
final_cores[, grep("_Proteins", names(final_cores), value = TRUE) := NULL ]

# Rename the core modules
names(final_cores) <- paste("CM", formatC( 1:ncol(final_cores), width = 2, flag = "0"), sep = "_")

# Write out the result file
fwrite(final_cores, "output_files/core_modules.csv")

