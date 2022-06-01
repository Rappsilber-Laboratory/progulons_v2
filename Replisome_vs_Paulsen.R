############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Compare the replisome prediction with the results from Paulsen et al (Mol Cell, 2009)
# Authors: G Kustatscher
# Date: April 2022
############################################################################################################# #

# Load the required libraries
library(data.table); library(readxl); library(ggplot2)

#### Prep the data ####

# Read in the replisome results (prepared previously as table S4)
REP <- as.data.table( read_xlsx("input_files/Supplementary_Table_4.xlsx", sheet = 1))

# Keep only relevant columns
REP <- REP[ !is.na(Mean_RF_score),
            .(Majority_protein_IDs, Simplified_protein_ID, Mean_RF_score, prot_in_prn, Feature_count, Used_for_positive_training, Gene_names, Screening_type)]

# Read in ID conversion file
ids <- fread("input_files/Paulsen_RefSeq_Uniprot.tab.gz")

# Simplify ID conversion
names(ids)[1] <- "RefSeq"
ids <- unique(ids[, .(RefSeq, Entry)])

# Keep only unique IDs
ids <- ids[ !duplicated(RefSeq) ]
ids <- ids[ !duplicated(Entry) ]

# Append RefSeq IDs to REP data
REP <- merge(REP, ids, by.x = "Simplified_protein_ID", by.y = "Entry")

# Load data from Paulsen et al
PAU <- as.data.table( read_xls("input_files/1-s2.0-S1097276509004596-mmc2.xls", skip = 8))

# Keep only relevant rows and columns
PAU <- PAU[ 5:.N ]
PAU <- PAU[, .(Accession, Significance, `Average %H2AX`, `Average Cell Number`, `Average %G1`, `Average %S`, `Average %G2/M`)]

# Append UniProt IDs to PAU data
PAU <- merge(PAU, ids, by.x = "Accession", by.y = "RefSeq")

# Merge PAU and REP data
DT <- merge(REP, PAU, by.x = "Simplified_protein_ID", by.y = "Entry")


#### Analyse the data ####

# Create plotting data
plot_dt <- DT[, sum( Significance == 4) / .N * 100 , by = .(prot_in_prn == "yes") ]
names( plot_dt ) <- c("Sample", "Sign. H2AX phenotype [%]")
plot_dt[, Sample := as.character(Sample) ]
plot_dt[ is.na(Sample),  Sample := "other genes\n(RF score < 0.55)" ]
plot_dt[ Sample == TRUE, Sample := "replisome progulon\n(RF score > 0.55)" ]

# Fisher's Exact test
DT[, .N , by = .(prot_in_prn == "yes", Significance == 4) ]
test_matrix <- matrix( c(16, 206, 246, 6605), 
                       nrow = 2, 
                       dimnames = list( c("Sign", "not_Sign"),
                                        c("in_repli", "not_in_repli")))
test_matrix
fisher.test(x = test_matrix, alternative = "greater")

# Plot a bar chart
p <- ggplot( plot_dt, aes(x = Sample, y = `Sign. H2AX phenotype [%]`))+
     geom_bar( stat = "identity" )+
     ylim(0,10)+   
     theme_bw()+theme(axis.title.x = element_blank(), panel.grid = element_blank())
p

ggsave("output_files/Reviewer_H2AX_figure.pdf", p, width = 8, height = 8, units = "cm")




