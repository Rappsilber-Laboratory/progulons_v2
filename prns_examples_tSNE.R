############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Produce plots showing three example progulons
# Authors: G Kustatscher
# Date: February 2022
############################################################################################################# #

# Load required libraries
library(data.table)
library(ggplot2)
library(treeClust)
library(Rtsne)
library(GA)
library(egg)

#### Load and prep data ####

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

# Load and prep ProteomeHD
ProHD <- read.csv("input_files/ProteomeHD_v1.csv", stringsAsFactors=FALSE)
rownames(ProHD) <- ProHD$Majority.protein.IDs                   # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]                # Keep only columns with SILAC ratios

# Load the manual annotation data
annot_P09 <- fread("input_files/Manual_annotation_P09.csv")
annot_P13 <- fread("input_files/Manual_annotation_P13.csv")
annot_P28 <- fread("input_files/Manual_annotation_P28.csv")


#### treeClust, tSNE and RF score plots for ATP synthase progulon (P09) ####

# Select progulon to be plotted
my_prn <- "P09"
my_annot <- annot_P09
  
# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(1)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = FALSE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = FALSE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
        geom_vline(xintercept = prns_con[ Progulon_ID == "P09", target_cut_off], size=0.25, linetype="dashed")+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta", alpha = 0.8)+
        geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
        scale_x_continuous(limits=c(0,1), expand = c(0,0))+
        scale_y_continuous(limits=c(0,300), expand = c(0,0))+
        xlab("Random forest score")+
        ylab("# available experiments")+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              panel.grid = element_blank(),
              axis.text.x=element_text(size=5), axis.text.y=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)), axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
              axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF
ggsave("output_files/ATP_synthase_RF.png", pRF, width = 5, height = 5, units= "cm", dpi = 600)


## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <- c(`ATP synthase` = "#DC3912", `aKG depletion` = "#11069e", `Complex I` = "#FF9900",
             `Complex II` = "#109618", `Complex III` = "#990099", `Complex IV` = "#0099C6",
             `FAO` = "#DD4477", `Uncharacterised` = "#bedf00")

pSNE_09 <- ggplot(SNE, aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_09


#### treeClust, tSNE and RF score plots for prefoldin progulon (P13) ####

# Select progulon to be plotted
my_prn <- "P13"
my_annot <- annot_P13

# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(1)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = FALSE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = FALSE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                 .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
        geom_vline(xintercept = prns_con[ Progulon_ID == "P13", target_cut_off], size=0.25, linetype="dashed")+
        geom_point(size=0.05, alpha=0.4, colour="grey50")+
        geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta")+
        geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
        #geom_point(size=0.2, data= plot_DT[Used_for_positive_training =="Used"], colour="magenta")+
        scale_x_continuous(limits=c(0,1), expand = c(0,0))+
        scale_y_continuous(limits=c(0,300), expand = c(0,0))+
        xlab("Random Forest score")+
        ylab("Experiments in which protein was quantified")+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)),
              axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF

## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <-  c(`Heat shock proteins` = "#70f0ff",            
              `hnRNPs` = "peru",
              `LSm` = "coral1",
              `Sm` = "lightgoldenrod",
              `Nuclear import/export` = "blue", 
              `Prefoldin` = "#DC3912",                       
              `Prolyl isomerase` = "#f8ff3d",                
              `Proteasome, core` = "#990099",                
              `Proteasome, regulatory` = "#ef6eef",          
              `R2TP/PFDL complex` = "#ef542f",               
              `Ribosome 40S` = "#13a51c",                    
              `Ribosome 60S` = "green",                    
              `Translation factors` = "#00a09b",             
              `TRiC` = "orange",                            
              `Tubulins & tub. chaperones` = "#0087af",      
              `Uncharacterised` = "#bedf00")          


pSNE_13 <- ggplot(SNE[ V1 > -28 ], aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_13

  
#### treeClust, tSNE and RF score plots for AP3 progulon (P28) ####

# Select progulon to be plotted
my_prn <- "P28"
my_annot <- annot_P28

# Subset ProHD to proteins from the Progulon
progulon_proteins <- prns[ Progulon_ID == my_prn & prot_in_prn == "yes", Protein_IDs]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(1)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = FALSE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = FALSE)
SNE <- as.data.frame(SNE$Y)
SNE$ID <- protein_IDs

# Merge SNE data with RF score information
SNE <- merge(SNE, prns[ Progulon_ID == my_prn, .(Protein_IDs, Mean_RF_score, Used_for_positive_training, Feature_count, Protein_names)],
             by.x = "ID", by.y = "Protein_IDs")

# Merge SNE data with annotation data
SNE <- merge(SNE, my_annot[,.(Protein_IDs, Manual_annotation, Highlight_genes)], by.x = "ID", by.y = "Protein_IDs", all.x = TRUE)


## Create RF score plot ##

# Create plotting data
plot_DT <- prns[ Progulon_ID == my_prn ,
                 .(Protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, prot_in_prn)]

pRF <- ggplot(plot_DT, aes(Mean_RF_score, Feature_count))+
        geom_vline(xintercept = prns_con[ Progulon_ID == "P28", target_cut_off], size=0.25, linetype="dashed")+
        geom_point(size=0.05, alpha=0.4, colour="grey50")+
        geom_point(size=0.1, data= plot_DT[ prot_in_prn == "yes"], colour="magenta")+
        geom_point(size=0.8, data= plot_DT[Used_for_positive_training =="Used"], shape=1)+
        #geom_point(size=0.2, data= plot_DT[Used_for_positive_training =="Used"], colour="magenta")+
        scale_x_continuous(limits=c(0,1), expand = c(0,0))+
        scale_y_continuous(limits=c(0,300), expand = c(0,0))+
        xlab("Random Forest score")+
        ylab("Experiments in which protein was quantified")+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title.x=element_text(size=6, margin=margin(1.5,0,0,0)),
              axis.title.y=element_text(size=6, margin=margin(0,0,0,0)))

pRF

## Create tSNE plot ##
SNE <- as.data.table(SNE)
SNE[ Manual_annotation == "", Manual_annotation := NA ]
SNE[ Highlight_genes == "", Highlight_genes := NA ]
SNE[, Manual_annotation := as.factor(Manual_annotation) ]

my_cols <-  c(`Adaptor protein complex 3` = "#990099",
              `Coatomer` = "lightgoldenrod",
              `COG complex, lobe A` = "blue",
              `COPII coat` = "#DC3912",
              `Exocyst complex` = "#f8ff3d",
              `GIT1/2-PIX complex` = "#990099",
              `Microtubule-based vesicle motility` = "#ef6eef",
              `NRZ vesicle tethering complex` = "green",
              `Other vesicle trafficking` = "#70f0ff",
              `Regulation of actin cytoskeleton` = "orange",
              `Uncharacterised` = "#bedf00")


pSNE_28 <- ggplot(SNE, aes(x=V1, y=V2, size = Mean_RF_score))+
            geom_point( fill="black", alpha=0.5, shape=21, stroke=0)+ 
            geom_point( aes( fill = Manual_annotation ), shape=21, stroke=0)+      
            geom_point( data = SNE[ Used_for_positive_training == "Used",], colour="black", fill=NA, show.legend = FALSE, shape=21, stroke=0.2)+
            geom_text( aes(label = Highlight_genes), size = 1.5, colour = "grey40")+
            xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
            scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
            scale_fill_manual( values = my_cols )+ 
            theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
                  panel.grid = element_blank(), panel.background=element_blank(),
                  axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25),
                  axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                  legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE_28


#### Progulon lineplots: seeds + top25 co-regulated proteins in 25 experiments ####

# Select ATP synthase progulon (P09) data
P09_training_prots <- prns[ Progulon_ID == "P09" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
P09_progulon_top25 <- prns[ Progulon_ID == "P09" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# Select Prefoldin (P13) data
P13_training_prots <- prns[ Progulon_ID == "P13" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
P13_progulon_top25 <- prns[ Progulon_ID == "P13" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# Select AP3 progulon (P28) data
P28_training_prots <- prns[ Progulon_ID == "P28" & Used_for_positive_training == "Used"                                                        , Protein_IDs]    # Get seed proteins (positive training proteins)
P28_progulon_top25 <- prns[ Progulon_ID == "P28" & Used_for_positive_training != "Used" & prot_in_prn == "yes"][ order(-Mean_RF_score) ][ 1:25 , Protein_IDs]    # Get top 25 of non-training progulon proteins

# These are all the proteins that need to be in the plot
my_prots <- unique( c( P09_training_prots, P09_progulon_top25, P13_training_prots, P13_progulon_top25, P28_training_prots, P28_progulon_top25 ))      

# Fitness function for a genetic algorithm, designed to pick a useful set of experiments to display
# "Useful" is to mean that (a) good coverage and (b) they are experiments where the proteins actually show some (even modest) change
fitness_f <- function(x){ E1 <- ceiling( x[1] )         # Get column index of the currently selected experiments (ratios)
                          E2 <- ceiling( x[2] )         # Need to use ceiling to work with integers
                          E3 <- ceiling( x[3] )
                          E4 <- ceiling( x[4] )
                          E5 <- ceiling( x[5] )
                          E6 <- ceiling( x[6] )
                          E7 <- ceiling( x[7] )
                          E8 <- ceiling( x[8] )
                          E9 <- ceiling( x[9] )
                          E10 <- ceiling( x[10] )
                          E11 <- ceiling( x[11] )  
                          E12 <- ceiling( x[12] )
                          E13 <- ceiling( x[13] )
                          E14 <- ceiling( x[14] )
                          E15 <- ceiling( x[15] )
                          E16 <- ceiling( x[16] )
                          E17 <- ceiling( x[17] )
                          E18 <- ceiling( x[18] )
                          E19 <- ceiling( x[19] )
                          E20 <- ceiling( x[20] )
                          E21 <- ceiling( x[21] )  
                          E22 <- ceiling( x[22] )
                          E23 <- ceiling( x[23] )
                          E24 <- ceiling( x[24] )
                          E25 <- ceiling( x[25] )
                          ratios <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13, E14, E15, E16, E17, E18, E19, E20, E21, E22, E23, E24, E25)
                          
                          # Select the subset of ProteomeHD that covers these proteins and these experiments (ratios)
                          temp_df <- ProHD[ rownames(ProHD) %in% my_prots, ratios]
                          
                          # What's the fraction of missing values across these proteins and experiments?
                          pc_NA <- sum(is.na(temp_df)) / (ncol(temp_df)*nrow(temp_df))
                          
                          # What is the median fold-change in this experiments, and how often does it exceed a certain threshold?
                          P09_changes <- apply( temp_df[ rownames(temp_df) %in% P09_training_prots,], 2, median, na.rm = TRUE )
                          P09_changes <- sum( abs( P09_changes ) > 0.5 , na.rm = TRUE)
                          P13_changes <- apply( temp_df[ rownames(temp_df) %in% P13_training_prots,], 2, median, na.rm = TRUE )
                          P13_changes <- sum( abs( P13_changes ) > 0.5 , na.rm = TRUE)
                          P28_changes <- apply( temp_df[ rownames(temp_df) %in% P28_training_prots,], 2, median, na.rm = TRUE )
                          P28_changes <- sum( abs( P28_changes ) > 0.5 , na.rm = TRUE)
                          
                          # Calculate the output
                          if( length(unique(ratios)) < 25 ){
                            fitness_output <- 0                                                          # If some experiments were picked more than once, the solution is invalid (score zero)
                            } else if( pc_NA > 0.05 ){
                              fitness_output <- 0                                                        # If there are more than 5% missing of values missinge, the display would not be informative (score zero)
                              } else
                                fitness_output <- sum( c(P09_changes, P13_changes, P28_changes))   # In how many experiments is the result above the threshold?
                          
                          # Return the fitness output
                          return(fitness_output) 
                          }                                                                

set.seed(1)
GA <- ga(type="real-valued", fitness = fitness_f, lower = rep(1,25), upper = rep(294,25))      # Genetic algorithm searching for good combination
my_ratios <- ceiling(GA@solution)          # The 25 ratios identified by GA as good for plotting
my_ratios <- as.integer( my_ratios[1,] )   # Simplified to integer and any duplicate solutions removed     
                          
# Create data frame for plotting
plot_data <- ProHD[ rownames(ProHD) %in% my_prots , my_ratios ]        # The relevant proteins and the relevant ratios
colnames(plot_data) <- paste("exp_", 1:25, sep="")

# Median-normalise
plot_data_exp_medians <- apply(plot_data, 2, median, na.rm=TRUE)       # Median fold-change of these proteins in these experiments
plot_data <- sweep(plot_data, 2, plot_data_exp_medians, FUN="-")       # Set that to zero

# Turn into melted data.table
plot_data$Protein_IDs <- rownames(plot_data)
plot_data <- as.data.table( plot_data )
plot_data <- melt(plot_data, id.vars = "Protein_IDs")

# Make the plot
plot_data[, experiment_number := as.integer( gsub("exp_", "", variable)) ]

pBase <- ggplot(plot_data, aes(x = experiment_number, y = value, group = Protein_IDs))+
         coord_cartesian(ylim = c(-2.5,2.5))+
         xlab("Experiments")+ylab("log2 SILAC ratio")+
         scale_x_continuous( breaks = seq(1,25,1), limits = c(1, 25), expand = c(0,0))+
         scale_y_continuous( breaks = c(-2, 0, 2))+
         theme(legend.position="left", axis.text=element_text(size=5), panel.grid = element_blank(),
                panel.background=element_blank(),
                axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25), axis.ticks.x = element_line(size=0.25), axis.ticks.y = element_line(size=0.25),
                axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
                legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

p09 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% P09_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% P09_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)

p13 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% P13_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% P13_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)

p28 <- pBase + geom_line( data = plot_data[ Protein_IDs %in% P28_progulon_top25 ], colour = "#36a9e1", alpha = 0.7, size = 0.25)+
               geom_line( data = plot_data[ Protein_IDs %in% P28_training_prots ], colour = "#e6007e", alpha = 0.7, size = 0.25)


#### Output one combined plot (except RF score) ####

# Remove the legend
pSNE_09_no_leg <- pSNE_09 + theme(legend.position = "none")
pSNE_13_no_leg <- pSNE_13 + theme(legend.position = "none")
pSNE_28_no_leg <- pSNE_28 + theme(legend.position = "none")

         p <- ggarrange( pSNE_09_no_leg, pSNE_13_no_leg, pSNE_28_no_leg, p09, p13, p28, ncol = 3, heights = c(2.1,1))
p_with_leg <- ggarrange( pSNE_09,        pSNE_13,        pSNE_28,        p09, p13, p28, ncol = 3, heights = c(2.1,1))

# Save the plots with and without legend
ggsave("output_files/Progulon_example_plots.pdf",        p, width=18, height=8.3, units= "cm")
ggsave("output_files/Progulon_example_plots_legend.pdf", p_with_leg, width=18, height=8.3, units= "cm")




