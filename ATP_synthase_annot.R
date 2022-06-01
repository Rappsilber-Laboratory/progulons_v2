############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: This script matches the ATP progulon data to the manual progulon annotations
# Authors: G Kustatscher
# Date: April 2022
############################################################################################################# #

# Load the required libraries
library(ggplot2); library(data.table)

#### Annotate ATP synthase progulon ####

# Read in the manual annotation file
annot <- fread("input_files/ATPsyn_manual_annot.csv")

# Load and prep the progulon scores
prns <- fread("output_files/Progulons.csv.gz") # Load the progulon data
connectivity <- fread("output_files/prn_connectivity.csv.gz") # Load connectivity data

# Find the minimum RF cut-off that will yield a significantly connected progulon 
prns_con <- connectivity[ Connectivity_p_value < 0.05, .(target_cut_off = min(RF_cut_off)), Progulon_ID ]

# Assign proteins to progulons
for(i in prns_con$Progulon_ID){
  target_RF_cutoff <- prns_con[ Progulon_ID == i, target_cut_off ]                    # Find progulon-specific RF score cut-off
  prns[ Progulon_ID == i & Mean_RF_score <= target_RF_cutoff, prot_in_prn := "no" ]   # Proteins below cut-off are not in progulon
  prns[ Progulon_ID == i & Mean_RF_score >  target_RF_cutoff, prot_in_prn := "yes" ]  # Proteins above cut-off are in progulon
}

# Extract the ATP synthase progulon
ATP <- prns[ Progulon_ID == "P09" ]

# Keep only relevant columns
ATP <- ATP[, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Protein_names, Gene_names, prot_in_prn)]

# Merge with annotations
ATP <- merge(annot, ATP, by = "Protein_IDs", all.y = TRUE)

# Keep only proteins that were either annotated or part of the progulon
ATP <- ATP[ !is.na(Function_or_Complex) | prot_in_prn == "yes" ] 

# Rearrange columns
setcolorder(x = ATP, neworder = c(1,4:8,2,3))

# Order by score
ATP <- ATP[ order(-Mean_RF_score) ]

# Write out the resulting table
fwrite(ATP, "output_files/ATP_synthase_prn_annot.csv")


#### Create cartoon figure colour templates ####

# To assist with the creation of the cartoon figure, write out mock images of each progulon gene, colour-coded by RF score
plotting_subset <- ATP[ !is.na(Function_or_Complex) ] # Define the subset of proteins to be included in this plot
plotting_subset <- plotting_subset[order(Function_or_Complex, -Mean_RF_score)]  # Set plotting order
plotting_subset[, ymin := 1:.N]        # Add arbitray dimensions for plotting
plotting_subset[, ymax := 1:.N+1]
plotting_subset[, ymean := 1:.N+0.5]

my_colours <- colorRampPalette(c("#36a9e1", "#FFFFFF", "#e6007e"))(100)  # Define a custom colour palette

p <- ggplot(plotting_subset, aes(xmin=0, xmax=1, ymin=ymin, ymax=ymax, fill=Mean_RF_score))+
      geom_rect(colour="black")+
      geom_text(aes(x=0.5, y=ymean, label=Gene_names))+
      geom_text(aes(x=1.5, y=ymean, label=Function_or_Complex), hjust=1)+
      scale_fill_gradientn(colours=my_colours, limits=c(0,1), breaks=seq(0,1,0.2))+
      theme(panel.background = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank())

p
ggsave("output_files/colour_coding_ATP_synthase_cartoon.pdf", p, width=17, height=100, units= "cm")



