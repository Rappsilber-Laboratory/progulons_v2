############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Calculate some basic statistics for the progulons
# Authors: G Kustatscher
# Date: February 2022
############################################################################################################# #

# Load the required libraries
library(ggplot2); library(data.table)


#### Load progulon data and assign proteins ####

# Load the progulon data
prns <- fread("output_files/Progulons.csv.gz")

# Load connectivity data
connectivity <- fread("output_files/prn_connectivity.csv.gz")

# Find the minimum RF cut-off that will yield a significantly connected progulon 
prns_con <- connectivity[ Connectivity_p_value < 0.05, .(target_cut_off = min(RF_cut_off)), Progulon_ID ]

# Assign proteins to progulons
for(i in prns_con$Progulon_ID){
  target_RF_cutoff <- prns_con[ Progulon_ID == i, target_cut_off ]                    # Find progulon-specific RF score cut-off
  prns[ Progulon_ID == i & Mean_RF_score <= target_RF_cutoff, prot_in_prn := "no" ]   # Proteins below cut-off are not in progulon
  prns[ Progulon_ID == i & Mean_RF_score >  target_RF_cutoff, prot_in_prn := "yes" ]  # Proteins above cut-off are in progulon
}


#### Get some stats ####

# How many proteins are, on mean / median, in each progulon?
proteins_in_progulon <- prns[ prot_in_prn == "yes", .N, by = Progulon_ID]
mean(   proteins_in_progulon[,N] )
median( proteins_in_progulon[,N] )
min(    proteins_in_progulon[,N] )
max(    proteins_in_progulon[,N] )

# How many training proteins have been used on mean/median?
training_prots_per_prog <- prns[ Used_for_positive_training == "Used", .N, by = Progulon_ID]
mean(   training_prots_per_prog[,N] )
median( training_prots_per_prog[,N] )
min(    training_prots_per_prog[,N] )
max(    training_prots_per_prog[,N] )

# How many proteins have been assigned to how many progulons?
progs_per_protein <- prns[ prot_in_prn == "yes", .N, by = Protein_IDs]
progs_per_protein[,.N]                                      # Proteins found in any progulon
progs_per_protein[,.N] / prns[,length(unique(Protein_IDs))] # Fraction of all analysed proteins
mean(   progs_per_protein[,N] )
median( progs_per_protein[,N] )
min(    progs_per_protein[,N] )
max(    progs_per_protein[,N] )


#### Plot progulon membership ####
progs_per_protein[, category := ifelse(N <= 3, N, "> 3")]
progs_per_protein[, category := factor( category, levels = c("> 3", "3", "2", "1"))]

p <- ggplot(progs_per_protein, aes(x = 1, fill = category))+
      geom_bar()+
      ylab("# proteins")+
      scale_fill_manual(values = c("#7FDBFF", "#39CCCC", "#0074D9", "#001f3f"))+
      theme(plot.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=7), 
            axis.text.y = element_text(size=6), axis.text.x = element_blank(), axis.line.y = element_line(size=0.25),
            axis.ticks.x = element_blank())

ggsave("output_files/Progs_per_protein.pdf", p, width = 5, height = 5, units = "cm")
