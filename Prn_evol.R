# Load required libraries
library(data.table); library(ggplot2); library(viridis)


#setwd("G:/Andere Computer/Mein Computer Neu/Work/Projects/ProteomeHD/Overhaul2022/")
setwd("")

#### Load and prep the data ####

# Load the mouse tissue dataset from Grabowski et al
df <- read.csv("input_files/mouse_SILAC_TPMs_log2_final_min8_features.csv", stringsAsFactors = FALSE)

# The table contains the IDs of the human one-to-one orthologs according to ENSEMBL, and I will only focus on those
df <- df[ !is.na(df$Human_Uniprot), ]

# Remove unnecessary tables
df$Ensembl <- NULL
df$Mouse_Uniprot <- NULL
df$Human_Ortholog <- NULL

# Load ProteomeHD
ProHD <- read.csv("input_files/ProteomeHD_v1_1.csv.gz", stringsAsFactors=FALSE)

# Simplify protein IDs
ProHD$SimpleID <- gsub(";.+", "", ProHD$Majority_protein_IDs)
ProHD$SimpleID <- gsub("-.+", "", ProHD$SimpleID)

# Remove duplicate IDs (isoforms)
ProHD <- ProHD[ !duplicated(ProHD$SimpleID) ,]

# Restrict both dataset to the overlapping set of proteins
   df <- df[  df$Human_Uniprot %in% ProHD$SimpleID   ,]
ProHD <- ProHD[ ProHD$SimpleID %in% df$Human_Uniprot ,]

# Keep only protein ratios of mouse tissues and median-normalise
rownames(df) <- df[, "Human_Uniprot"]                           # Set protein IDs as rownames
mouse <- df[, grep("SILAC_", colnames(df))]                     # Keep only columns with SILAC ratios
tmouse <- t(mouse)                                              # Transpose and median-normalise to remove SILAC mixing artifacts
mouse_people_medians <- apply(tmouse, 1, median, na.rm=TRUE)
tmouse_mn <- sweep(tmouse, 1, mouse_people_medians, FUN="-")

# Keep only protein ratios of ProteomeHD and median-normalise
rownames(ProHD) <- ProHD$SimpleID                               # Set protein IDs as rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]                # Keep only columns with SILAC ratios
feature_count <- apply(ProHD, 1, function(x){ sum(!is.na(x))})  # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 30 ,]                          # Only proteins with >= 30 feature counts will have RF scores
tProHD <- t(ProHD)                                              # Transpose and median-normalise to remove SILAC mixing artifacts
ProHD_people_medians <- apply(tProHD, 1, median, na.rm=TRUE)
tProHD_mn <- sweep(tProHD, 1, ProHD_people_medians, FUN="-")


#### Get correlations for all relevant protein pairs ####

# For the mouse data
mouse_cor <- cor(tmouse_mn, use = "pairwise.complete.obs", method = "spearman")  # All pairwise combinations
mouse_cor <- as.data.table( melt( mouse_cor ))                                  # Convert it to a long data table
mouse_cor <- mouse_cor[, .( Gene_1 = as.character(Var1),                        # Re-name
                            Gene_2 = as.character(Var2),
                            mouse_rho = value ) ]
mouse_cor <- mouse_cor[ Gene_1 > Gene_2 ]                                       # Remove duplicate pairs (incl self-comparisons)

# For ProteomeHD
ProHD_cor <- cor(tProHD_mn, use = "pairwise.complete.obs", method = "spearman")  # All pairwise combinations
ProHD_cor <- as.data.table( melt( ProHD_cor ))                                  # Convert it to a long data table
ProHD_cor <- ProHD_cor[, .( Gene_1 = as.character(Var1),                        # Re-name
                            Gene_2 = as.character(Var2),
                            ProHD_rho = value ) ]
ProHD_cor <- ProHD_cor[ Gene_1 > Gene_2 ]                                       # Remove duplicate pairs (incl self-comparisons)

# Combine the data
setkey(mouse_cor, Gene_1, Gene_2)
setkey(ProHD_cor, Gene_1, Gene_2)
DT <- merge(mouse_cor, ProHD_cor)

# Remove comparisons with missing values
DT <- DT[ complete.cases(DT) ]

rm( list = ls()[! ls() %in% c("DT")] )  


#### Assess distribution of ratios ####

DT[, lapply(.SD, median), .SDcols = grep("_rho", names(DT))]   # Output medians

# Plot across gene distribution
ggplot( DT[ sample(.N, 1000000) ])+
  geom_histogram( aes(mouse_rho), binwidth = 0.05, boundary = 0.025, fill = "grey80")+
  geom_histogram( aes(ProHD_rho), binwidth = 0.05, boundary = 0.025, fill = NA, colour = "red")+
  xlim(-1,1)+
  theme_bw()


#### How well is protein coexpression conserved for proteins from no, any or shared progulons? ####

# I would like to compare the coexpression in human ProHD and mouse tissues. I will break it down into three groups:
# (a) pairs where neither protein has been assigned to any progulon
# (b) pairs where the two proteins were assigned to different progulons
# (c) pairs where the both proteins were assigned to the same progulon

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

# Simplify protein IDs
prns[, SimpleID := gsub(";.+", "", Protein_IDs) ][, SimpleID := gsub("-.+", "", SimpleID)]


# Assign the coexpression pairs to the three groups
prot_in_any_prn <- prns[ prot_in_prn == "yes" , unique(SimpleID) ]   # These proteins have been assigned to any progulon
DT[  Gene_1 %in% prot_in_any_prn  &   Gene_2 %in% prot_in_any_prn,   both_in_a_prn := "yes" ]  # Annotation 1
DT[(!Gene_1 %in% prot_in_any_prn) & (!Gene_2 %in% prot_in_any_prn), neither_in_prn := "yes" ]  # Annotation 2

for(i in unique(prns$Progulon_ID)){                                                            # Annotation 3
  prots_in_current_prn <- prns[ Progulon_ID == i & prot_in_prn == "yes", unique(SimpleID) ]
  DT[ Gene_1 %in% prots_in_current_prn & Gene_2 %in% prots_in_current_prn, shared_prn := "yes" ]
}


#### Plot the results ####

# Looking at coexpressed pairs in ProHD ("coexpressed" defined as rho > 0.5),
# how many of them are also coexpressed in the mouse dataset, per category?
neither <- DT[ ProHD_rho > 0.5 & neither_in_prn == "yes"                       , sum(mouse_rho > 0.5)/.N*100 ]
notSame <- DT[ ProHD_rho > 0.5 & both_in_a_prn  == "yes" & is.na(shared_prn)   , sum(mouse_rho > 0.5)/.N*100 ]
sharedP <- DT[ ProHD_rho > 0.5 & both_in_a_prn  == "yes" & shared_prn == "yes" , sum(mouse_rho > 0.5)/.N*100 ]

dt <- rbind( data.table( type = "neither", value = neither ),
             data.table( type = "notSame", value = notSame ),
             data.table( type = "sharedP", value = sharedP ))

# Get the statistical significance
neither <- DT[ ProHD_rho > 0.5 & neither_in_prn == "yes"                       , .N, mouse_rho > 0.5 ]
notSame <- DT[ ProHD_rho > 0.5 & both_in_a_prn  == "yes" & is.na(shared_prn)   , .N, mouse_rho > 0.5 ]
sharedP <- DT[ ProHD_rho > 0.5 & both_in_a_prn  == "yes" & shared_prn == "yes" , .N, mouse_rho > 0.5 ]

neither <- neither[ order(-mouse_rho) , N ]   # Order so coexpressed (TRUE) comes before not coexpressed (FALSE)
notSame <- notSame[ order(-mouse_rho) , N ]   # Order so coexpressed (TRUE) comes before not coexpressed (FALSE)
sharedP <- sharedP[ order(-mouse_rho) , N ]   # Order so coexpressed (TRUE) comes before not coexpressed (FALSE)

neither_to_notSame <- data.frame( notSame, neither )
rownames(neither_to_notSame) <- c("coexpressed", "not_coexpressed")

notSame_to_sharedP <- data.frame( sharedP, notSame )
rownames(notSame_to_sharedP) <- c("coexpressed", "not_coexpressed")

pval_neither_to_notSame <- round(fisher.test(neither_to_notSame)$p.value,7)
pval_notSame_to_sharedP <- round(fisher.test(notSame_to_sharedP)$p.value,25)
fisher.test(notSame_to_sharedP)
pval_notSame_to_sharedP <- 2.2e-16 # updated as the $pvalue rounds to 0, 2.2e-16 is seemingly the standard lowest value given out by the fisher test

# Make the plot
p1 <- ggplot(dt, aes(x = type, y = value))+
      geom_bar(stat = "identity")+
      ylab("Human protein pairs also coexpressed in mouse [%]")+
      scale_y_continuous( limits = c(0,45), expand = c(0,0))+
      annotate(geom = "text", label = pval_neither_to_notSame, x = 1.5, y = 20, size = 2)+
      annotate(geom = "text", label = pval_notSame_to_sharedP, x = 2.5, y = 40, size = 2)+  
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 6), axis.title.x = element_blank(),
            axis.text = element_text(size = 5, colour = "black"), axis.ticks.y = element_line(size = 0.25),
            axis.ticks.x = element_blank())

p1
ggsave("output_files/Mouse_conservation.pdf", p1,
       width = 2.5, height = 3.5, units = "cm")





