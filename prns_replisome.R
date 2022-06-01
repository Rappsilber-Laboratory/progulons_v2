## This script processes the progulonFinder results for the replisome into scatter and tSNE plots,
## while highlighting the proteins that were selected for siRNA follow-up.

# Load the necessary libraries
library(data.table); library(ggplot2); library(treeClust); library(Rtsne); library(plotly); library(egg)

# Note: Set working directory to the script file location 

#### Prep the Replisome Progulon data ####

## ProteomeHD

# Load ProteomeHD (no normalisation necessary)
ProHD <- read.csv("input_files/ProteomeHD_v1_1.csv.gz", stringsAsFactors=FALSE)

# Get an "annotation frame" to which I can then map the results and further annotations
Annotation_frame <- as.data.table( ProHD[,c("Majority_protein_IDs", "Simplified_protein_ID", "Protein_names", "Gene_names")] )

# Extract SILAC ratios
rownames(ProHD) <- ProHD$Majority_protein_IDs       # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]    # Keep only columns with SILAC ratios


## DNA replication progulon

# Load the replisome progulon (result from progulonFinder)
PRN_replisome <- fread("input_files/RF_score_Replisome_all.csv")

# Polish column names
names(PRN_replisome) <- gsub(" ", "_", names(PRN_replisome), fixed = TRUE)
names(PRN_replisome) <- gsub(".", "_", names(PRN_replisome), fixed = TRUE)

# Keep only necessary columns
PRN_replisome <- PRN_replisome[, .(Majority_protein_IDs, Mean_RF_score, Feature_count, Used_for_positive_training, Training_label)]

# Turn any empty STRINGs into NAs
PRN_replisome[ Used_for_positive_training == "" , Used_for_positive_training := NA ]
PRN_replisome[ Training_label == ""             , Training_label := NA ]

# Assign proteins to the progulon: Keep all proteins that score > 0.55
PRN_replisome[ Mean_RF_score > 0.55, prot_in_prn := "yes" ]
PRN_replisome[, .N, prot_in_prn]

# Subset ProHD data to proteins from the Progulon
progulon_proteins <- PRN_replisome[ prot_in_prn == "yes", Majority_protein_IDs ]
progulon_data <- ProHD[ rownames(ProHD) %in% progulon_proteins ,]

# Use treeClust to learn a dissimilarity matrix
set.seed(42)
tc_dist <- treeClust.dist(progulon_data, d.num = 2, verbose = FALSE)
protein_IDs <- attr(tc_dist, "Labels")

# Use tSNE to reduce the dissimilarity matrix down to a 2D map
set.seed(42)
SNE <- Rtsne(tc_dist, is_distance = TRUE, theta = 0.0, verbose = TRUE)
SNE <- as.data.table(SNE$Y)
SNE$Majority_protein_IDs <- protein_IDs
names(SNE) <- c("tSNE_dim_1", "tSNE_dim_2", "Majority_protein_IDs")

# Merge replisome data (RF score and tSNE dimensions) with annotation data
PRN_replisome <- merge(PRN_replisome, Annotation_frame, by = "Majority_protein_IDs", all = TRUE)
PRN_replisome <- merge(PRN_replisome, SNE, by = "Majority_protein_IDs", all = TRUE)


#### Assign GO annotations ####

# Read in GO annotation for DNA replication (downloaded from QuickGO, reviewed human proteins,
# mapping to GO:0006260 with qualifiers "part of" or "involved in")
GO <- fread("input_files/QuickGO_DNA_replication_GO0006260.tsv")
GO[, GO_annotation := "DNA replication"]
GO <- unique( GO[,.(SimpleID = `GENE PRODUCT ID`, GO_annotation)] )

# Merge Progulon data with GO annotation data
PRN_replisome <- merge(PRN_replisome, GO, by.x = "Simplified_protein_ID", by.y = "SimpleID", all.x = TRUE)  


#### Assign data from siRNA follow-up screening ####

# Load list of proteins that we followed up
siRNA_results <- fread("input_files/Candidates_siRNA_screen.csv")

# Extract the names of the candidates
targets <- unique( siRNA_results[, .(Protein_ID_in_ProHD, Candidate, Type_of_candidate) ] )

# Assign the two protein IDs belonging to ASF1A and ASF1B, respectively
temp <- targets[ Candidate == "ASF1A+B" ] 
temp <- rbind( temp, temp )
temp[ 1, Candidate := "ASF1A" ]
temp[ 2, Candidate := "ASF1B" ]
temp[ Candidate == "ASF1A", Protein_ID_in_ProHD := PRN_replisome[ Gene_names == "ASF1A" , Majority_protein_IDs ]]
temp[ Candidate == "ASF1B", Protein_ID_in_ProHD := PRN_replisome[ Gene_names == "ASF1B" , Majority_protein_IDs ]]
targets <- targets[ Candidate != "ASF1A+B" ] 
targets <- rbind( temp, targets )

# Merge Progulon data with target annotation
targets <- targets[, .(Majority_protein_IDs = Protein_ID_in_ProHD, Screening_type = Type_of_candidate) ]
PRN_replisome <- merge(PRN_replisome, targets, by = c("Majority_protein_IDs"), all.x = TRUE)

# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "siRNA_results") ))


#### Make the plots ####

# Use NA as factor levels (exlude = NULL) and set them to come before DNA replication
PRN_replisome[, GO_annotation := factor(GO_annotation, exclude = NULL, levels = c(NA, "DNA replication")) ]

# Set fill colours
my_cols <- c(positive_control = "#A2DE00", negative_control = "#0000DE", new_replication_candidate = "#FF0051", candidate_from_NCC_exps = "#FF6800")

# A scatterplot and histogram showing enrichment of GO ontology terms
p1 <- ggplot(PRN_replisome, aes(x = Mean_RF_score, y = Feature_count, colour = Screening_type))+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_point(size=0.01, alpha=0.3, colour="grey50")+
        geom_point(data = PRN_replisome[ GO_annotation == "DNA replication" ], colour = "cyan4", size = 0.2)+
        geom_point(data = PRN_replisome[ Screening_type != "candidate_from_NCC_exps" ], size = 0.2, show.legend = FALSE)+   # Do not highlight the targets selected based solely on NCC data
        geom_point(data = PRN_replisome[ Used_for_positive_training == "Used"], colour = "black", fill=NA, size = 0.7, shape = 1, stroke=0.2)+
        xlab("Random forest score")+
        ylab("# available experiments")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(25,300), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))

p2 <- ggplot(PRN_replisome, aes(x = Mean_RF_score, fill = GO_annotation))+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+      
        geom_histogram( binwidth = 0.1, boundary = 0, position = "fill", show.legend = FALSE)+
        scale_fill_manual( values = "cyan4")+
        xlab("Random forest score")+
        ylab("GO: DNA replication [%]")+
        scale_x_continuous(limits=c(0,1), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), expand = c(0,0), breaks = c(0,0.5,1))+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text.y=element_text(size=5), axis.text.x = element_blank(),
              axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(), plot.background = element_blank(),
              axis.title.y=element_text(size=6), axis.title.x = element_blank())
        
# Combine and output plot
p3 <- ggarrange(p2, p1, heights = c(1,4))
ggsave("output_files/RF_plot.pdf", p3, width=4.7, height=5.5, units= "cm", dpi = 600)


# Plot the replisome progulon

# Need to get all fill colours into column, otherwise the size scaling is inconsistent
SNE_plot_data <- PRN_replisome[ !is.na(tSNE_dim_1) ]
SNE_plot_data[ GO_annotation == "DNA replication"           , SNE_plot_fill := "GO: DNA replication" ]
SNE_plot_data[ Screening_type == "positive_control"         , SNE_plot_fill := "screen: pos ctrl" ]
SNE_plot_data[ Screening_type == "new_replication_candidate", SNE_plot_fill := "screen: candidate" ]

pSNE <- ggplot(SNE_plot_data, aes( x = tSNE_dim_1, y = tSNE_dim_2, colour = Used_for_positive_training, fill = SNE_plot_fill, size = Mean_RF_score, label = Gene_names))+
          geom_point( shape = 21, stroke = 0.25 , alpha = 0.8)+
          scale_fill_manual( values = c("cyan4", "#FF0051", "#A2DE00"), na.value = "grey50")+
          scale_colour_manual( values = c("black"), na.value = "grey50")+
          geom_text(  data = SNE_plot_data[ !is.na(Screening_type) ], colour = "#FF0051", size = 1.7, show.legend = FALSE)+
          xlab("tSNE Dimension 1")+ylab("tSNE Dimension 2")+
          scale_size(limits=c(0.5,1), range=c(0.3,2.1), guide = guide_legend(title = "Random\nforest\nscore"))+
          theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(),
              panel.grid = element_blank(), panel.background = element_rect(fill = NA, colour = "black", size = 0.25),
              axis.title=element_text(size=6), legend.key = element_blank(), legend.text = element_text(size=6),
              legend.title = element_text(size=6), legend.key.height = unit(0.2,"cm"), legend.key.width = unit(0.2, "cm"))

pSNE
ggsave("output_files/Replisome_tSNE.pdf", pSNE, width=9.45, height=6, units= "cm")

# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "my_cols", "siRNA_results") ))


#### siRNA results bar chart ####

## Plot the results of the screen

# Get column order
my_levels <- siRNA_results[, sum(SD_score, na.rm = TRUE), by = Candidate ][ order(-V1) , Candidate ]

# Set column order
siRNA_results[, Candidate := factor( Candidate, levels = my_levels )]

pScreen <- ggplot( siRNA_results[ Type_of_candidate != "candidate_from_NCC_exps" ],
                   aes(x = Candidate, y = SD_score, fill = Process ))+
            geom_bar( stat = "identity" )+
            geom_hline( yintercept = 13/2, linetype = "dashed", size = 0.25, colour = "grey50")+
            geom_hline( yintercept = 13/3, linetype = "dashed", size = 0.25, colour = "grey50")+
            scale_y_continuous( limits = c(0,13), breaks = seq(0,13,1))+
            ylab("combined siRNA screen score")+
            scale_fill_manual( values = c("chartreuse3", "yellow2", "mediumslateblue"))+
            theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                  axis.text.y=element_text(size=5), axis.text.x=element_text(size=5, angle = 90, hjust = 1), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
                  plot.background = element_blank(), axis.title.y=element_text(size=6), axis.title.x = element_blank(),
                  legend.title = element_blank(), legend.text = element_text(size=6), legend.position = "top",
                  panel.grid.major.y = element_line(colour = "grey50", size = 0.25, linetype = "dotted"))

pScreen
ggsave("output_files/Main_screen_results.pdf", pScreen, width=9, height=5, units= "cm")

# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "my_cols", "siRNA_results") ))


#### Control 1: What about REMOVING NCC data from the prediction ####

## p4a: ProgulonFinder without the 15 NCC ratios

# Load a progulonFinder result that was run on ProteomeHD minus the 15 NCC ratios
PRN_without_NCC <- fread("input_files/RF_score_Replisome_without_NCC.csv")
PRN_without_NCC <- PRN_without_NCC[, .(Majority_protein_IDs = Majority.protein.IDs, RF_score_without_NCC = `Mean RF score`)]

# Merge with the full progulon data
PRN_replisome <- merge(PRN_replisome, PRN_without_NCC, by = "Majority_protein_IDs", all.x = TRUE)

# Calculate the R-squared
r2_full_withoutNCC <- PRN_replisome[, cor(Mean_RF_score, RF_score_without_NCC, use = "pairwise.complete.obs")^2 ]

# Make a scatterplot
p4a <- ggplot(PRN_replisome, aes(y = Mean_RF_score, x = RF_score_without_NCC, colour = Screening_type))+
        geom_hline(yintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_smooth( method = lm , size = 0.25, colour = "black")+
        geom_point(data = PRN_replisome[ Screening_type != "candidate_from_NCC_exps" ], size = 0.2, show.legend = FALSE)+ 
        annotate("text", x= 0.15, y=0.85, label = paste("R2", round(r2_full_withoutNCC, 3), sep = " "), size = 2)+
        xlab("RF score, without 15 NCC experiments")+
        ylab("RF score, full ProteomeHD")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))


## p4b: ProgulonFinder without 15 random ratios (to see if removing NCC ratios is worse than removing other ratios)

# Load a progulonFinder result that was run on ProteomeHD minus 15 ratios chose at random
PRN_without_15_ran_1 <- fread("input_files/RF_score_Replisome_without15random_1.csv")
PRN_without_15_ran_1 <- PRN_without_15_ran_1[, .(Majority_protein_IDs = Majority.protein.IDs,
                                                 RF_score_without_15_random_ratios_1 = `Mean RF score`)]

# Merge with the full progulon data
PRN_replisome <- merge(PRN_replisome, PRN_without_15_ran_1, by = "Majority_protein_IDs", all.x = TRUE)

# Calculate the R-squared
r2_full_without_15_ran_1 <- PRN_replisome[, cor(Mean_RF_score, RF_score_without_15_random_ratios_1, use = "pairwise.complete.obs")^2 ]

# Make a scatterplot
p4b <- ggplot(PRN_replisome, aes(y = Mean_RF_score, x = RF_score_without_15_random_ratios_1, colour = Screening_type))+
        geom_hline(yintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_smooth( method = lm , size = 0.25, colour = "black")+
        geom_point(data = PRN_replisome[ Screening_type != "candidate_from_NCC_exps" ], size = 0.2, show.legend = FALSE)+ 
        annotate("text", x= 0.15, y=0.85, label = paste("R2", round(r2_full_without_15_ran_1, 3), sep = " "), size = 2)+
        xlab("RF score without 15 experiments\nselected at random #1")+
        ylab("RF score, full ProteomeHD")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))


## p4c: ProgulonFinder without 15 random ratios (to see if removing NCC ratios is worse than removing other ratios)

# Load a progulonFinder result that was run on ProteomeHD minus 15 ratios chose at random
PRN_without_15_ran_2 <- fread("input_files/RF_score_Replisome_without15random_2.csv")
PRN_without_15_ran_2 <- PRN_without_15_ran_2[, .(Majority_protein_IDs = Majority.protein.IDs,
                                                 RF_score_without_15_random_ratios_2 = `Mean RF score`)]

# Merge with the full progulon data
PRN_replisome <- merge(PRN_replisome, PRN_without_15_ran_2, by = "Majority_protein_IDs", all.x = TRUE)

# Calculate the R-squared
r2_full_without_15_ran_2 <- PRN_replisome[, cor(Mean_RF_score, RF_score_without_15_random_ratios_2, use = "pairwise.complete.obs")^2 ]

# Make a scatterplot
p4c <- ggplot(PRN_replisome, aes(y = Mean_RF_score, x = RF_score_without_15_random_ratios_2, colour = Screening_type))+
        geom_hline(yintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_smooth( method = lm , size = 0.25, colour = "black")+
        geom_point(data = PRN_replisome[ Screening_type != "candidate_from_NCC_exps" ], size = 0.2, show.legend = FALSE)+ 
        annotate("text", x= 0.15, y=0.85, label = paste("R2", round(r2_full_without_15_ran_2, 3), sep = " "), size = 2)+
        xlab("RF score without 15 experiments\nselected at random #2")+
        ylab("RF score, full ProteomeHD")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))


## p4d: ProgulonFinder without 15 random ratios (to see if removing NCC ratios is worse than removing other ratios)

# Load a progulonFinder result that was run on ProteomeHD minus 15 ratios chose at random
PRN_without_15_ran_3 <- fread("input_files/RF_score_Replisome_without15random_3.csv")
PRN_without_15_ran_3 <- PRN_without_15_ran_3[, .(Majority_protein_IDs = Majority.protein.IDs,
                                                 RF_score_without_15_random_ratios_3 = `Mean RF score`)]

# Merge with the full progulon data
PRN_replisome <- merge(PRN_replisome, PRN_without_15_ran_3, by = "Majority_protein_IDs", all.x = TRUE)

# Calculate the R-squared
r2_full_without_15_ran_3 <- PRN_replisome[, cor(Mean_RF_score, RF_score_without_15_random_ratios_3, use = "pairwise.complete.obs")^2 ]

# Make a scatterplot
p4d <- ggplot(PRN_replisome, aes(y = Mean_RF_score, x = RF_score_without_15_random_ratios_3, colour = Screening_type))+
        geom_hline(yintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_smooth( method = lm , size = 0.25, colour = "black")+
        geom_point(data = PRN_replisome[ Screening_type != "candidate_from_NCC_exps" ], size = 0.2, show.legend = FALSE)+ 
        annotate("text", x= 0.15, y=0.85, label = paste("R2", round(r2_full_without_15_ran_3, 3), sep = " "), size = 2)+
        xlab("RF score without 15 experiments\nselected at random #3")+
        ylab("RF score, full ProteomeHD")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))

# Combine the output plot
p4 <- ggarrange(p4a, p4b, p4c, p4d, nrow = 1)
ggsave("output_files/Without_NCC_or_random.png", p4, width=18.3, height=5, units= "cm", dpi = 600)


# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "my_cols", "siRNA_results") ))


#### Could we have predicted the same using ONLY NCC data? ####

# Load a progulonFinder result that was run using only the 15 NCC ratios in ProHD
PRN_NCC_only <- fread("input_files/RF_score_Replisome_NCConly.csv")

# Polish column names
names(PRN_NCC_only) <- gsub(" ", "_", names(PRN_NCC_only), fixed = TRUE)
names(PRN_NCC_only) <- gsub(".", "_", names(PRN_NCC_only), fixed = TRUE)

# Keep only relevant columns and re-name them to distinguish from standard approach
PRN_NCC_only <- PRN_NCC_only[, .(Majority_protein_IDs, RF_score_NCC_only = Mean_RF_score, NCC_only_feature_count = Feature_count) ]

# Merge with the rest of the data
PRN_replisome <- merge(PRN_replisome, PRN_NCC_only, by = c("Majority_protein_IDs"), all.x = TRUE)

# Plot the difference to the full ProHD & highlight additional proteins that were followed up
p5 <- ggplot(PRN_replisome, aes(y = Mean_RF_score, x = RF_score_NCC_only, colour = Screening_type, label = Gene_names))+
        geom_hline(yintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.25)+
        geom_hline(yintercept = PRN_replisome[ !is.na(Training_label) , median(Mean_RF_score) ], colour = "green", size = 0.25)+
        geom_vline(xintercept = PRN_replisome[ !is.na(Training_label) , median(RF_score_NCC_only, na.rm = TRUE) ], colour = "green", size = 0.25)+
        geom_hline(yintercept = PRN_replisome[ GO_annotation == "DNA replication" , median(Mean_RF_score, na.rm = TRUE) ], colour = "orange", size = 0.25)+
        geom_vline(xintercept = PRN_replisome[ GO_annotation == "DNA replication" , median(RF_score_NCC_only, na.rm = TRUE) ], colour = "orange", size = 0.25)+
        geom_point(size=0.01, alpha=0.4, colour="grey50")+
        geom_point(data = PRN_replisome[ GO_annotation == "DNA replication" ], colour = "cyan4", size = 0.1)+
        geom_point(data = PRN_replisome[ Used_for_positive_training == "Used"], colour = "black", fill=NA, size = 0.7, shape = 1, stroke=0.2)+
        geom_point(data = PRN_replisome[ !is.na(Screening_type) ], size = 0.2, show.legend = FALSE)+ 
        geom_text( data = PRN_replisome[ !is.na(Screening_type) ], aes( colour = Screening_type), size = 1.7 , show.legend = FALSE)+
        xlab("RF score, using only the 15 NCC experiments")+
        ylab("RF score, full ProteomeHD")+
        scale_x_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2), expand = c(0,0))+
        scale_colour_manual( values = my_cols )+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks = element_line(size=0.25), plot.background = element_blank(),
              axis.title=element_text(size=6))


#### NCC-only prediction: siRNA results bar chart ####

pScreen2 <- ggplot( siRNA_results, aes(x = Candidate, y = SD_score, fill = Process ))+
            geom_bar( stat = "identity" )+
            geom_hline( yintercept = 13/2, linetype = "dashed", size = 0.25, colour = "grey50")+
            geom_hline( yintercept = 13/3, linetype = "dashed", size = 0.25, colour = "grey50")+
            scale_y_continuous( limits = c(0,13), breaks = seq(0,13,1))+
            ylab("combined siRNA screen score")+
            coord_flip()+
            scale_fill_manual( values = c("chartreuse3", "yellow2", "mediumslateblue"))+
            theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                  axis.text.y=element_text(size=5), axis.text.x=element_text(size=5, angle = 90, hjust = 1), axis.ticks = element_line(size=0.25),
                  plot.background = element_blank(), axis.title.y=element_text(size=6), axis.title.x = element_blank(),
                  legend.title = element_blank(), legend.text = element_text(size=6), legend.position = "top",
                  panel.grid.major.y = element_line(colour = "grey50", size = 0.25, linetype = "dotted"))

# Combine the plots and save them
pNCC_only <- ggarrange(p5, pScreen2, nrow = 1, widths = c(2.75,1))
ggsave("output_files/NCC_screen_results.pdf", pNCC_only, width=18.3, height=12, units= "cm")

# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "my_cols", "siRNA_results") ))


#### Control 2: What about using a random set of 15 ratios from the prediction ####

# We run progulonFinder three times on a random subset of 15 SILAC ratios. Load these data and merge
# them with the rest of the data
PRN_ran15_1 <- fread("input_files/RF_score_Replisome_only15random_1.csv")
PRN_ran15_2 <- fread("input_files/RF_score_Replisome_only15random_2.csv")
PRN_ran15_3 <- fread("input_files/RF_score_Replisome_only15random_3.csv")

PRN_ran15_1 <- PRN_ran15_1[, .(Majority_protein_IDs = Majority.protein.IDs, RF_score_15_random_exp_1 = `Mean RF score`)]
PRN_ran15_2 <- PRN_ran15_2[, .(Majority_protein_IDs = Majority.protein.IDs, RF_score_15_random_exp_2 = `Mean RF score`)]
PRN_ran15_3 <- PRN_ran15_3[, .(Majority_protein_IDs = Majority.protein.IDs, RF_score_15_random_exp_3 = `Mean RF score`)]

PRN_replisome <- merge(PRN_replisome, PRN_ran15_1, by = "Majority_protein_IDs", all.x = TRUE)
PRN_replisome <- merge(PRN_replisome, PRN_ran15_2, by = "Majority_protein_IDs", all.x = TRUE)
PRN_replisome <- merge(PRN_replisome, PRN_ran15_3, by = "Majority_protein_IDs", all.x = TRUE)


## Compare using boxplots, but need to quantile-normalise to get similar score distributions

# Create a new annotation column, which shows either DNA replication factors or mitochondrial proteins
PRN_replisome[ Protein_names %like% "itochondri",  Boxplot_annotation := "Mitochondrial"   ]
PRN_replisome[ GO_annotation == "DNA replication", Boxplot_annotation := "DNA replication" ]
PRN_replisome[ is.na(Boxplot_annotation),          Boxplot_annotation := "Other" ]
PRN_replisome[, Boxplot_annotation := factor( Boxplot_annotation, levels = c("Other", "DNA replication", "Mitochondrial"))]

# Extract the relevant data
bxpl_data <- PRN_replisome[, .(Majority_protein_IDs, Boxplot_annotation, Mean_RF_score, RF_score_without_NCC,
                               RF_score_NCC_only, RF_score_15_random_exp_1, RF_score_15_random_exp_2, RF_score_15_random_exp_3)]

# Remove missing values for quantile normalisation
ratio_cols <- grep("score", names(bxpl_data), value = TRUE)  # Define ratio columns
bxpl_data <- bxpl_data[ complete.cases( bxpl_data[, ratio_cols, with = FALSE ] )]


# # Define quantile normalisation function (modified from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/)
# quantile_normalisation <- function(df){
#   df_rank <- apply( df, 2, rank, ties.method = "min" )
#   df_sorted <- data.frame( apply( df, 2, sort ))
#   df_mean <- apply( df_sorted, 1, mean)
#   
#   index_to_mean <- function(my_index, my_mean){
#     return( my_mean[my_index] )
#   }
#   
#   df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
#   rownames(df_final) <- rownames(df)
#   df_final <- as.data.table( df_final )
#   return(df_final)
# }
# 
# # Quantile normalise the score to ensure even distribution
# bxpl_data[, c(ratio_cols) := quantile_normalisation(.SD), .SDcols = ratio_cols ]

# Melt the data
bxpl_data <- melt(bxpl_data, id.vars = c("Majority_protein_IDs", "Boxplot_annotation"))

# Make variable names more informative
bxpl_data[ variable == "Mean_RF_score",        variable := "Full ProteomeHD" ]
bxpl_data[ variable == "RF_score_without_NCC", variable := "ProHD without 15x NCC" ]
bxpl_data[ variable == "RF_score_NCC_only",    variable := "Only 15x NCC experiments" ]
bxpl_data[ variable == "RF_score_15_random_exp_1", variable := "15x random experiments #1" ]
bxpl_data[ variable == "RF_score_15_random_exp_2", variable := "15x random experiments #2" ]
bxpl_data[ variable == "RF_score_15_random_exp_3", variable := "15x random experiments #3" ]

# Plot the boxplot
p6 <- ggplot(bxpl_data, aes(x = variable, y = value, fill = Boxplot_annotation))+
        geom_boxplot(notch = TRUE, outlier.colour = NA, size = 0.25)+
        scale_fill_manual(values = c("grey50", "cyan4", "orange3"))+
        ylab("Random forest score")+ # Removed quantile normalisation
        scale_y_continuous( expand = c(0,0), limits = c(0,1))+
        theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
              axis.text=element_text(size=5), axis.ticks.y = element_line(size=0.25), axis.ticks.x = element_blank(),
              plot.background = element_blank(), axis.title.y=element_text(size=6), axis.title.x = element_blank(),
              legend.title = element_blank(), legend.text = element_text(size=6),
              panel.grid.major = element_line(colour = "grey50", size = 0.25, linetype = "dotted"))

p6
ggsave("output_files/NCC_only_boxplot.pdf", p6, width = 18.3, height = 5, units= "cm")

# Display protein numbers
bxpl_data[, .N, .(Boxplot_annotation, variable)]

# Clear the workspace
rm( list = setdiff( ls(), c("PRN_replisome", "my_cols", "siRNA_results") ))


#### Write out the supplementary table ####

PRN_replisome <- PRN_replisome[ order(-Mean_RF_score) ]
fwrite(PRN_replisome, "output_files/Replisome_predictions.csv")








