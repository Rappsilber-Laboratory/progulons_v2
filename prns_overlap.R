############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Calculate the overlap between progulons
# Authors: G Kustatscher
# Date: February 2022
############################################################################################################# #

# Load the required libraries
library(ggplot2); library(data.table)

#### Load and prep input data ####

# Load the progulon data
prns <- fread("output_files/Progulon_Scores.csv.gz")

# Load connectivity data
connectivity <- fread("output_files/prn_connectivity.csv.gz")

# Find the minimum RF cut-off that will yield a significantly connected progulon 
prns_con <- connectivity[ Connectivity_p_value < 0.05, .(target_cut_off = min(RF_cut_off)), Progulon_ID ]

# Turn all values to NA that are below the progulon-specific RF cut-off
for(i in prns_con$Progulon_ID){
  target_RF_cutoff <- prns_con[ Progulon_ID == i, target_cut_off ]  # Find progulon-specific RF score cut-off
  prns[ get(i) <= target_RF_cutoff, c(i) := NA ]                    # Set to NA all values below target cut-off in that progulon
}


#### Find overlap between progulons ####

# Get handle for relevant columns
score_columns <- grep("P[[:digit:]]", names(prns), value = TRUE)

# For pairwise comparison, create all progulon combinations
prn_combinations <- as.data.table( t( combn( score_columns, 2 )))

# Calculate number of proteins overlapping between different progulons
overlap <- integer()
for(i in 1:prn_combinations[,.N]){
  prn_pair <- prns[, as.character( prn_combinations[i] ) , with = FALSE ]    # Get the scores for each pair
  overlap[i] <- prn_pair[ complete.cases(prn_pair) , .N ]                    # Count the proteins for which they both have a score (i.e. proteins are in the progulon)
}
prn_combinations[, overlap := overlap ]


#### Find the overlap in percent ####

# Calculate the number of progulon proteins
prn_N <- prns[, lapply(.SD, function(x){ sum( !is.na(x) ) }), .SDcols = score_columns]
prn_N <- melt(prn_N, variable.name = "prnID", value.name = "prn_N", measure.vars = names(prn_N) )

# Match progulon size to overlap table
prn_combinations <- merge(prn_combinations, prn_N, by.x = "V1", by.y = "prnID")
prn_combinations <- merge(prn_combinations, prn_N, by.x = "V2", by.y = "prnID", suffixes = c("_V1", "_V2"))
setcolorder(x = prn_combinations, neworder = c(2,1,3:5))

# Calculate percentage of overlap relative to the bigger of the two progulons
prn_combinations[ prn_N_V1 >= prn_N_V2 , largerPRN := prn_N_V1 ]
prn_combinations[ prn_N_V1  < prn_N_V2 , largerPRN := prn_N_V2 ]
prn_combinations[, overlap_pc := overlap / (largerPRN/100) ]


#### Get the order in which progulons are plotted in the correlation panel ####

# Load the correlation data
cor_combis <- fread("output_files/ProgulonCor.csv.gz")

# Load the progulaon annotations
prn_annot <- fread("input_files/Progulon_annotation.csv")
setnames(prn_annot, old = c("Progulon_ID", "Progulon_name"), new = c("Progulon", "Function"))

# Expand data to be able to get a complete matrix, i.e. append duplicates
cor_combis <- rbind(cor_combis[, .(PRN_A, PRN_B, RHO)],
                    cor_combis[, .(PRN_B = PRN_A, PRN_A = PRN_B, RHO)])

# Append functional annotation
cor_combis[, Function_1 := prn_annot[ match( cor_combis[, PRN_A], prn_annot[, Progulon] ), Function ]]
cor_combis[, Function_2 := prn_annot[ match( cor_combis[, PRN_B], prn_annot[, Progulon] ), Function ]]

# Cast into a correlation matrix
cor_mat <- dcast( cor_combis, Function_1 ~ Function_2, value.var = "RHO" )
my_rownames <- cor_mat[, Function_1] 
cor_mat[, Function_1 := NULL ]
cor_mat <- as.data.frame( cor_mat )
rownames(cor_mat) <- my_rownames

# Group progulons by correlation
my_dist <- as.dist( (1-cor_mat)/2 )
my_dist[ is.na(my_dist) ] <- 0        # For the few combinations that had no non-overlapping proteins, set the resulting NA to zero
my_clust <- hclust(my_dist, method = "average")
new_prn_order <- rownames(cor_mat)[ my_clust$order ]


#### Plot overlap between progulons ####

# Append functional annotation
prn_combinations <- merge(prn_combinations, prn_annot, by.x = "V1", by.y = "Progulon")
prn_combinations <- merge(prn_combinations, prn_annot, by.x = "V2", by.y = "Progulon", suffixes = c("_1", "_2"))

# Rearrange progulons (via factor levels) in clustered order
prn_combinations[, Function_1 := factor(Function_1, levels = new_prn_order )]
prn_combinations[, Function_2 := factor(Function_2, levels = new_prn_order )]

# Expand data into a complete "overlap matrix" (minus diagonal)
prn_combinations <- rbind(prn_combinations[, .(Function_1, Function_2,                           overlap, overlap_pc) ],
                          prn_combinations[, .(Function_1 = Function_2, Function_2 = Function_1, overlap, overlap_pc) ])

# Create the plot
p <- ggplot(prn_combinations, aes(Function_1, Function_2, size = overlap, colour = overlap_pc))+
     geom_point()+
     geom_vline(xintercept = seq(0.5,50.5,1), size = 0.25, colour = "grey80")+
     geom_hline(yintercept = seq(0.5,50.5,1), size = 0.25, colour = "grey80")+
     scale_size_continuous(name = "# proteins", limits = c(1,350), range = c(0,3), breaks = seq(0,500,50))+
     scale_colour_viridis_c(name = "% overlap", limits = c(1,70), option = "E", direction = -1)+
     theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black"), 
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), 
           axis.ticks = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6, angle = 90, hjust=1, vjust = 0.5),
           legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
           legend.key.size = unit(3, "mm"))

p1 <- p + theme( legend.position = "none")
p2 <- p + theme( legend.position = "right")

# Save the plot
ggsave("output_files/Progulon_overlap.pdf", p1, width = 11.5, height = 11.5, units = "cm")
ggsave("output_files/Progulon_overlap_legend.pdf", p2, width = 11.5, height = 11.5, units = "cm")

















