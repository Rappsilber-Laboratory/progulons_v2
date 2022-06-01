############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Polish the progulons identified by the KNIME workflow
# Authors: G Kustatscher
# Date: January 2022
############################################################################################################# #

# Load the required libraries
library(data.table)

# Load in the progulon scores from a folder than contains the individual result files obtained by running progulonFinder on the 72 core modules
filenames <- list.files("output_files/progulons/", pattern = "RF_score", full.names = FALSE)
filepaths <- list.files("output_files/progulons/", pattern = "RF_score", full.names = TRUE)
progulons <- lapply(filepaths, fread, colClasses = c("character", "numeric", "numeric", "integer", "character", "character", "character", "character", "numeric", "character"))
names(progulons) <- gsub(".csv", "", filenames, fixed = TRUE)
names(progulons) <- gsub("RF_score_", "", names(progulons), fixed = TRUE)

# Load in the stats
filenames <- list.files("output_files/progulons/", pattern = "Stats", full.names = FALSE)
filepaths <- list.files("output_files/progulons/", pattern = "Stats", full.names = TRUE)
progulon_stats <- lapply(filepaths, fread)
names(progulon_stats) <- gsub(".csv", "", filenames, fixed = TRUE)
names(progulon_stats) <- gsub("Stats_", "", names(progulon_stats), fixed = TRUE)

# Keep only progulons with AUC > 0.99
progulon_AUCs <- sapply(progulon_stats, function(x){ x[Observation == "Area Under Curve", as.numeric(Value)]})
progulons_to_keep <- names( which( progulon_AUCs > 0.99 ))
progulons <- progulons[ names(progulons) %in% progulons_to_keep ]

# Keep only progulons that have at least four (cross-validated) training protein among the top 10 hits
# (i.e. for the smallest training sets, the entire set has to be in the top 10 retrieved hits)
ts_among_top10 <- sapply( progulons, function(x){ x[ order(-`Mean RF score`) ][ 1:10 , sum( Training_label == "Positive") ]})
progulons_to_keep <- names( which( ts_among_top10 >= 4 ))
progulons <- progulons[ names(progulons) %in% progulons_to_keep ]

# Append progulon identifier
progulon_identifiers <- paste("PRN", formatC(1:length(progulons), width=2, flag="0"), sep="")
for(i in 1:length(progulons)){ progulons[[i]][, Progulon_ID := progulon_identifiers[i] ] }

# Bring all progulons into one data.table
progulons_DT <- rbindlist(progulons)

# Remove unnessary columns
progulons_DT[, Fasta.headers := NULL ]

# Rename and reorder remaining columns
progulons_DT <- progulons_DT[, .(Progulon_ID, Protein_IDs = Majority.protein.IDs, Mean_RF_score = `Mean RF score`,
                                 Score_StDev = `St Dev of score`, Training_label, Used_for_positive_training,
                                 Feature_count, Protein_names = Protein.names, Gene_names = Gene.names, Summed_intensities) ]

# Write out complete progulon table in long format
fwrite(progulons_DT, "output_files/Progulons.csv.gz")

# Create a wide-format progulon-by-progulon version of this table
progulons_wide_DT <- dcast(progulons_DT, formula = Protein_IDs + Feature_count + Protein_names + Gene_names ~ Progulon_ID,
                           value.var = "Mean_RF_score" )

# Write out the wide format
fwrite(progulons_wide_DT, "output_files/Progulon_Scores.csv.gz")

