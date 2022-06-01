library(ggplot2); library(data.table); library(perm); library(readxl); library(grid); library(egg)

#### Get Progulon associations ####

setwd("")

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


#### Load and prep the matched mRNA - protein input data for breast cancer cell lines ####

# Load the breast cancer dataset from Lapek et al (Table S3, Nature Biotechnology, 2017)
pro <- read_xlsx("input_files/Supplementary_Table_3_Protein_mRNA_Profiles_36_Cell_Lines.xlsx", 
                 sheet = "Proteome_Profiles")
rna <- read_xlsx("input_files/Supplementary_Table_3_Protein_mRNA_Profiles_36_Cell_Lines.xlsx", 
                 sheet = "mRNA_Profiles")

# Re-name one cell line because there is a typo
colnames(rna)[ colnames(rna) == "HS578T" ] <- "Hs578T"

# Now mRNA and protein datasets have the same cell lines in the same order?
sum( colnames(pro) == colnames(rna) ) == ncol(pro)

# Rename columns
colnames(pro)[ 2:37 ] <- paste("pro_", colnames(pro)[ 2:37 ] , sep = "")
colnames(rna)[ 2:37 ] <- paste("rna_", colnames(rna)[ 2:37 ] , sep = "")

# Combine the data
df <- merge(pro, rna, by = "Uniprot Gene Name" )

# Write out Uniprot IDs and convert them to accession numbers using the Retrieve tool on Uniprot website
write.csv(df$`Uniprot Gene Name`, "output_files/temp_Uniprot_IDs.csv")
id_conversion <- fread("input_files/ID_conversion.tab")
df$UniprotIDs <- id_conversion$Entry[ match( df$`Uniprot Gene Name`, id_conversion$`Entry name` ) ]

# Remove unmatched gene / protein IDs
df <- df[ !is.na(df$UniprotIDs) ,]

# Keep only proteins which are in our progulon analysis
df <- df[ df$UniprotIDs %in% unique(prns$SimpleID) ,]

# Get the mRNA and protein ratios
rownames(df) <- df[, "UniprotIDs"]
rna <- df[, grep("rna_", colnames(df))]
pro <- df[, grep("pro_", colnames(df))]

# Transpose the matrixes and sweep out row-medians to avoid artificial correlations
trna <- t(rna)
rna_people_medians <- apply(trna, 1, median, na.rm=TRUE)
trna_mn <- sweep(trna, 1, rna_people_medians, FUN="-")

tpro <- t(pro)
pro_people_medians <- apply(tpro, 1, median, na.rm=TRUE)
tpro_mn <- sweep(tpro, 1, pro_people_medians, FUN="-")

# Restrict also the PRNS dataset to the mRNAs / proteins detected here
prns <- prns[ SimpleID %in% df$UniprotIDs ]


#### Get all possible pairwise mRNA - mRNA correlations in this dataset ####

RNA_RNA_cor <- cor(trna_mn, use = "pairwise.complete.obs", method = "spearman")   # All pairwise combinations
RNA_RNA_cor <- as.data.table( melt( RNA_RNA_cor ))                                # Convert it to a long data table
RNA_RNA_cor <- RNA_RNA_cor[, .( Gene_1 = as.character(Var1),                      # Re-name
                                Gene_2 = as.character(Var2),
                                RNA_RNA_rho = value ) ]
RNA_RNA_cor <- RNA_RNA_cor[ Gene_1 >= Gene_2 ]                                    # Remove duplicate pairs (keep self-comparisons)


#### Get all possible pairwise protein - protein correlations in this dataset ####

pro_pro_cor <- cor(tpro_mn, use = "pairwise.complete.obs", method = "spearman")  # All pairwise combinations
pro_pro_cor <- as.data.table( melt( pro_pro_cor ))                                # Convert it to a long data table
pro_pro_cor <- pro_pro_cor[, .( Gene_1 = as.character(Var1),                      # Re-name
                                Gene_2 = as.character(Var2),
                                pro_pro_rho = value ) ]
pro_pro_cor <- pro_pro_cor[ Gene_1 >= Gene_2 ]                                    # Remove duplicate pairs (keep self-comparisons)


#### Get all possible pairwise mRNA - protein correlations in this dataset ####

# Modify column names to include data type
RNA_pro_cor <- cor(trna_mn, tpro_mn, use = "pairwise.complete.obs", method = "spearman")  # All pairwise combinations
RNA_pro_cor <- as.data.table( melt( RNA_pro_cor ))                                # Convert it to a long data table
RNA_pro_cor <- RNA_pro_cor[, .( Gene_1 = as.character(Var1),                      # Re-name
                                Gene_2 = as.character(Var2),
                                RNA_pro_rho = value ) ]
RNA_pro_cor <- RNA_pro_cor[ Gene_1 >= Gene_2 ]                                    # Remove duplicate pairs (keep self-comparisons)


#### Merge all pairwise correlations into one table ####
setkey( RNA_RNA_cor , Gene_1 , Gene_2 )
setkey( pro_pro_cor , Gene_1 , Gene_2 )
setkey( RNA_pro_cor , Gene_1 , Gene_2 )

DT <- merge(RNA_RNA_cor, pro_pro_cor)
DT <- merge(DT, RNA_pro_cor)

# Clear workspace
rm( list = ls()[! ls() %in% c("DT", "prns", "trna_mn", "tpro_mn")] )  


#### Assess distribution of ratios ####

DT[ Gene_1 != Gene_2 , lapply(.SD, median), .SDcols = grep("_rho", names(DT))]   # Output across gene medians
DT[ Gene_1 == Gene_2 , lapply(.SD, median), .SDcols = grep("_rho", names(DT))]   # Output within gene medians

# Plot across gene distribution
pH1 <- ggplot( DT[ Gene_1 != Gene_2 ])+
        geom_histogram( aes(RNA_RNA_rho), binwidth = 0.05, boundary = 0.025, fill = NA, colour = "#00A0BB", size = 0.25)+
        geom_histogram( aes(pro_pro_rho), binwidth = 0.05, boundary = 0.025, fill = NA, colour = "#EC008C", size = 0.25)+
        xlim(-1,1)+
        xlab("spearman's correlation coefficient")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))

# Plot within gene distribution
pH2<- ggplot( DT[ Gene_1 == Gene_2 ] )+
       geom_histogram( aes(RNA_pro_rho), binwidth = 0.05, boundary = 0.025, fill = "grey80")+
       geom_vline( xintercept = DT[ Gene_1 == Gene_2 , median(RNA_pro_rho) ], size = 0.25, linetype = "dashed")+
       annotate("text", x = 0.6, y = 400, size = 2.5,
                label = paste("median =", round( DT[ Gene_1 == Gene_2 , median(RNA_pro_rho) ] , 2 )))+
       xlim(-1,1)+
       xlab("mRNA - protein correlation [rho]")+
       ylab("Number of gene pairs")+
       theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
             legend.position = "none", strip.background = element_blank(),
             strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))

# Combine and output plot
pH <- ggarrange(pH1, pH2, nrow = 1)
ggsave("output_files/Breast_cancer_histograms.pdf", pH, width = 14, height = 4.9, units = "cm")


#### Get pairwise correlations for each progulon ####

DT_prn <- data.table()  # Initialise the result table: This will be long-format table with pairwise rhos assigned to progulons.
                        # Within-gene comparisons removed from RNA-RNA and Pro-Pro pairs, across-gene comparisons removed from the RNA-Pro pairs

for(i in unique(prns$Progulon_ID)){      # Loop through all progulons
  
  # These are the proteins belonging to the current progulon
  progulon_proteins <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]
  
  # Get the protein (gene) pairs relevant for the current progulon
  temp_prn_DT <- DT[ Gene_1 %in% progulon_proteins & Gene_2 %in% progulon_proteins ]
  
  # These are the across-gene-pair values
  across_gene <- temp_prn_DT[ Gene_1 != Gene_2 ]
  
  # These are the same-gene pairs
  within_gene <- temp_prn_DT[ Gene_1 == Gene_2 ]
  
  # Combine into one long-format table
  temp_DT1 <- melt( across_gene[, .(RNA_RNA_rho, pro_pro_rho) ] , measure.vars = c("RNA_RNA_rho", "pro_pro_rho"))
  temp_DT2 <- melt( within_gene[, .(RNA_pro_rho) ]              , measure.vars = c("RNA_pro_rho"))
  temp_DT <- rbind( temp_DT1, temp_DT2 )
  temp_DT[, Progulon_ID := i ]
  
  # Merge  
  DT_prn <- rbind(DT_prn, temp_DT)

  # Display progess  
  print(i)
}

# Rename the columns
colnames(DT_prn) <- c("cor_type", "rho", "Progulon_ID")

# Create a summary table
DT_summary <- DT_prn[, .(median_prn_rho = median(rho)), by = .(cor_type, Progulon_ID) ]
DT_summary <- dcast(DT_summary, Progulon_ID ~ cor_type , value.var = "median_prn_rho")
DT_summary <- DT_summary[, .(Progulon_ID, med_RNA_RNA_rho = RNA_RNA_rho, med_pro_pro_rho = pro_pro_rho, med_RNA_pro_rho = RNA_pro_rho) ]
  

#### Statistical significance #1: Protein vs mRNA coexpression of each progulon ####

# For each progulon, I want to calculate if the extent of coexpression is stronger on protein than mRNA level
# I use an independent Mann Whitney U test. Remember that some gene pairs are assigend to multiple progulons.

RNAvsPROsign <- data.table()                # Initialise result table
for(i in unique(DT_prn$Progulon_ID)){       
  temp_RNA <- DT_prn[ Progulon_ID == i & cor_type == "RNA_RNA_rho" , rho ]
  temp_pro <- DT_prn[ Progulon_ID == i & cor_type == "pro_pro_rho" , rho ]
  temp_pvalue <- wilcox.test(temp_pro, temp_RNA)$p.value
  temp_dt <- data.table( Progulon_ID = i, p_PCORvsRCOR = temp_pvalue)
  RNAvsPROsign <- rbind( RNAvsPROsign, temp_dt )
  print(i)
}

# Append p-values to summary table
DT_summary <- merge(DT_summary, RNAvsPROsign, by = "Progulon_ID")


#### Statistical significance #2: mRNA to protein correlation of progulon-specific genes ####

# This is a different question: There is some correlation between the mRNA and the protein expression values of each gene
# Do genes assigned to different progulon differ in how strongly their RNA and protein expression values are correlated?
# Because of the different sizes and distributions involved I'm using a permutation test for this one

# Create the universe
all_RNA_pro <- DT[ Gene_1 == Gene_2, .(Gene_1, Gene_2, RNA_pro_rho) ]   # These are all within-gene RNA-to-protein correlations

# For each progulon, test the potential shift in mean correlation between the genes that belong to this progulon
# and the remaining genes in the "universe"

RNAtoPROsign <- data.table()             # Initialise result table
for(i in unique(prns$Progulon_ID)){      # Loop through all progulons
  progulon_proteins <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # These are the proteins belonging to the current progulon
  prn_cors <- all_RNA_pro[ Gene_1 %in% progulon_proteins        , RNA_pro_rho ]    # RNA-to-protein correlations of the current progulon genes
  remaining_cors <- all_RNA_pro[ !Gene_1 %in% progulon_proteins , RNA_pro_rho ]    # and of the remaining genes in the universe
  
  temp_pvalue <- permTS(prn_cors, remaining_cors, alternative = "two.sided",       # Calcule p-values by permutation
                        method = "exact.mc", control = permControl(nmc = 10000,    # Using 10000 Monte Carlo replications
                        setSEED = FALSE, tsmethod = "abs"))$p.value
  temp_dt <- data.table( Progulon_ID = i, p_RNAtoPRO = temp_pvalue)
  RNAtoPROsign <- rbind( RNAtoPROsign, temp_dt )
  print(i)
}
  
# Append p-values to summary table
DT_summary <- merge(DT_summary, RNAtoPROsign, by = "Progulon_ID")


#### Prepare progulon annotation and plotting order ####

# Load the prn-prn correlations
cor_combis <- fread("output_files/ProgulonCor.csv.gz")

# Load the manual prn annotations
prn_annot <- fread("input_files/Progulon_annotation.csv")
colnames(prn_annot) <- c("Progulon","Function")
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
my_dist[ is.na(my_dist) ] <- 0 # set the few NAs to 0
my_clust <- hclust(my_dist)
new_prn_order <- rownames(cor_mat)[ my_clust$order ]

# Append the progulon annotation to the data
DT_prn[, prn_function := prn_annot[ match( DT_prn$Progulon_ID , prn_annot$Progulon ), Function ]]
DT_summary[, prn_function := prn_annot[ match( DT_summary$Progulon_ID , prn_annot$Progulon ), Function ]]

# Rearrange progulons (via factor levels) in clustered order
DT_prn[, prn_function := factor(prn_function, levels = new_prn_order )]
DT_summary[, prn_function := factor(prn_function, levels = new_prn_order )]

# Clear workspace
rm( list = ls()[! ls() %in% c("DT", "prns", "DT_prn", "DT_summary", "trna_mn", "tpro_mn")] )  


#### Plot 1: Example histograms ####

# Get the relevant subset of the data 
plot_dt1 <- DT_prn[ cor_type != "RNA_pro_rho" & Progulon_ID %in% c("P25", "P29") ]
plot_dt1[, prn_function := factor(prn_function, levels = c("Ribosome", "DNA replication"))]

# Create a separate data.table for the medians
dt1_medians <- plot_dt1[, median(rho), .(cor_type, prn_function) ]

# Create a separate data.table for the number of gene pairs
dt1_N_pairs <- plot_dt1[, .N, .(cor_type, prn_function) ]

# Create the plot
p1 <- ggplot( plot_dt1 , aes(x = rho, fill = cor_type ))+
        facet_wrap(~prn_function, nrow = 2)+
        geom_vline( data = dt1_medians, aes( xintercept = V1, colour = cor_type ),
                    linetype = "dotted", size = 0.25)+
        geom_histogram(position = "identity", alpha = 0.7, binwidth = 0.05, center = 0.025)+
        geom_text( data = dt1_N_pairs, aes( label = N , x = -0.5, y = 300), size = 2)+
        scale_fill_manual( values = c("#00A0BB", "#EC008C"))+
        scale_colour_manual( values = c("#00A0BB", "#EC008C"))+
        scale_x_continuous( limits = c(-1,1))+
        xlab("Gene coexpression\n[rho]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))

p1
ggsave("output_files/Haas_Example_histograms_update.pdf", p1,
       width = 3.3, height = 6, units = "cm")


#### Plot 2: mRNA-mRNA vs protein-protein scatterplot ####
sigrhoest <- DT_summary[, cor.test(med_RNA_RNA_rho, med_pro_pro_rho, method = "spearman" )$estimate ] 
sigrhopva <- DT_summary[, cor.test(med_RNA_RNA_rho, med_pro_pro_rho, method = "spearman" )$p.value  ] 
sigRHOest <- DT_summary[, cor.test(med_RNA_RNA_rho, med_pro_pro_rho, method = "spearman")$estimate ]
sigRHOpva <- DT_summary[, cor.test(med_RNA_RNA_rho, med_pro_pro_rho, method = "spearman")$p.value  ]

sigrho <- paste("rho", round(sigrhoest, 2), ", p value", signif(sigrhopva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

p2 <- ggplot(DT_summary, aes(x = med_RNA_RNA_rho, y = med_pro_pro_rho))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_smooth( method = "lm", size = 0.25, colour = "orange", se = FALSE, fullrange = TRUE)+ 
        #annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1, linetype = "dashed", size = 0.25)+
        geom_text(data = DT_summary[ Progulon_ID %in% c("P25", "P04","P05", "P11", "P31", "P29", "P21", "P09", "P23","P15")],
                  aes(label = prn_function), size = 2, colour = "blue", hjust = -0.1)+
        scale_x_continuous( limits = c(0,1), expand = c(0,0), breaks = seq(0,1,0.2))+
        scale_y_continuous( limits = c(-0.05,1), expand = c(0,0), breaks = seq(0,1,0.2))+
        annotate(geom = "text", x = 0.02, y = 0.9, size = 2, hjust = 0, label = sigrho )+      
        annotate(geom = "text", x = 0.02, y = 0.8, size = 2, hjust = 0, label = sigRHO )+  
        xlab("median mRNA coexpression [rho]")+
        ylab("median protein coexpression [rho]")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p2 
ggsave("output_files/Haas_RNARNAvsProPro_scatter_update.pdf", p2,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 3: mRNA-to-protein boxplot for subset of progulons  ####

# Define and order a subset of progulons for plotting
plotting_subset <- DT_prn[ cor_type == "RNA_pro_rho" & Progulon_ID %in% c("P04","P05","P25","P31", "P29", "P21", "P09") ]
plotting_subset[, prn_subset_function := factor(prn_function, levels = c("mito. 39S ribosome","mito. 28S ribosome","Ribosome", "Nucleosome", "DNA replication", 
                                                                         "Exosome", "ATP Synthase"))]

p3 <- ggplot( plotting_subset, aes(x = prn_subset_function, y = rho ))+
        geom_boxplot(size = 0.25, outlier.fill = "white", outlier.shape = 21, outlier.stroke = 0.25, outlier.size = 0.5)+
        geom_hline( yintercept = DT_prn[ cor_type == "RNA_pro_rho" , median(rho), prn_function ][, mean(V1) ], linetype = "dashed", size = 0.25, colour = "grey50")+
        scale_y_continuous( limits = c(-0.4,1), expand = c(0,0), breaks = seq(-1,1,0.2))+
        ylab("mRNA to protein correlation [rho]")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 0.25, colour = "black"),
              axis.title.y = element_text(size = 6), axis.title.x = element_blank(),
              axis.text = element_text(size = 5, colour = "black"), legend.position = "none",
              axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.25))

p3
ggsave("output_files/Haas_RNAtoPro_Subset_boxplot_update.pdf", p3,
       width = 4, height = 4.8, units = "cm")


#### Plot 4: mRNA-to-protein boxplot for all progulons  ####

p4 <- ggplot( DT_prn[ cor_type == "RNA_pro_rho" ], aes(x = prn_function, y = rho ))+
        geom_boxplot(size = 0.25, outlier.fill = "white", outlier.shape = 21, outlier.stroke = 0.25, outlier.size = 0.5)+
        geom_hline( yintercept = DT_prn[ cor_type == "RNA_pro_rho" , median(rho), prn_function ][, mean(V1) ], linetype = "dashed", size = 0.25, colour = "grey50")+
        scale_y_continuous( limits = c(-0.5,1), expand = c(0,0), breaks = seq(-1,1,0.2))+
        ylab("mRNA to protein correlation [rho]")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 0.25, colour = "black"),
              axis.title.y = element_text(size = 7), axis.title.x = element_blank(),
              axis.text.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 7, colour = "black", angle = 90, hjust = 1), 
              legend.position = "none", axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 0.25))

p4
ggsave("output_files/Haas_RNAtoPro_boxplot_update.pdf", p4,
       width = 18, height = 10, units = "cm")

# Clear workspace
rm( list = ls()[! ls() %in% c("DT", "prns", "DT_prn", "DT_summary", "trna_mn", "tpro_mn", "p2")] )  


#### Is coexpression dependent on the scale of expression variation? ####  
  
# It is often assumed that large gene expression changes - typical of induced or tissue-specific genes -
# are driven by transcriptional changes, whereas smaller changes - typical of housekeeping genes - may be
# the result of post-transcriptional changes. Here I want to measure if the scale of the expression changes
# determines how well mRNA changes correlate with protein changes, or mRNA-mRNA and protein-protein

# Instead of standard variation I use a more robust measure of scale, the median absolute deviation (MAD)

# Determine mRNA and protein expression variation per gene
mRNA_mad <- apply( trna_mn,  2, mad, na.rm = TRUE )
mRNA_mad <- data.table( Gene = names(mRNA_mad), mRNA_mad = mRNA_mad) 

prot_mad <- apply( tpro_mn, 2, mad, na.rm = TRUE )
prot_mad <- data.table( Gene = names(prot_mad), prot_mad = prot_mad) 

gene_mad <- merge( mRNA_mad, prot_mad )

# Calculate the median MADs per progulon
prn_mad <- data.table()                  # Initialise result table
for(i in unique(prns$Progulon_ID)){
  temp_genes <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]
  temp_med_mRNA_mad <- gene_mad[ Gene %in% temp_genes , .(med_mRNA_mad = median(mRNA_mad)) ]
  temp_med_prot_mad <- gene_mad[ Gene %in% temp_genes , .(med_prot_mad = median(prot_mad)) ]
  prn_mad <- rbind(prn_mad, data.table( Progulon_ID = i, temp_med_mRNA_mad, temp_med_prot_mad))
}

# Append results to the summary table
DT_summary <- merge(DT_summary, prn_mad, by = "Progulon_ID")


#### Write out (prototype of) the supplementary table ####

fwrite( DT_summary , "mRNAvsProtein/Haas/Haas_Rna_Pro_PRN_update.csv" )

DT_summary <- fread("mRNAvsProtein/Haas/Haas_Rna_Pro_PRN_update.csv")

#### Plot 5: scale of mRNA expression variation vs scale of protein variation, scatterplot ####

# This shows that progulons with large mRNA expression variation will generally also have 
# larger protein expression variation

# Calculate significance of correlation and turn into plot labels
sigrhoest <- DT_summary[, cor.test(med_mRNA_mad, med_prot_mad, method = "spearman" )$estimate ] 
sigrhopva <- DT_summary[, cor.test(med_mRNA_mad, med_prot_mad, method = "spearman" )$p.value  ] 
sigRHOest <- DT_summary[, cor.test(med_mRNA_mad, med_prot_mad, method = "spearman")$estimate ]
sigRHOpva <- DT_summary[, cor.test(med_mRNA_mad, med_prot_mad, method = "spearman")$p.value  ]

sigrho <- paste("rho", round(sigrhoest, 2), ", p value", signif(sigrhopva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Make the plot
p5 <- ggplot(DT_summary, aes(x = med_mRNA_mad, y = med_prot_mad))+
        geom_smooth( method = "lm", size = 0.25, colour = "orange", se = FALSE, fullrange = TRUE)+    
        geom_point( alpha = 0.5, size = 1 )+
        geom_text(data = DT_summary[ Progulon_ID %in% c("P04","P05","P25", "P09", "P31","P02","P16","P26")], aes(label = prn_function), size = 2, colour = "blue", hjust = -0.1)+
        scale_x_continuous( limits = c(0,0.151), expand = c(0,0), breaks = seq(0,0.3,0.05))+
        scale_y_continuous( limits = c(0,1.8), expand = c(0,0), breaks = seq(0,2,0.4))+
        xlab("Scale of mRNA expression variation [median MAD]")+
        ylab("Scale of protein expression variation [median MAD]")+
        annotate(geom = "text", x = 0.02, y = 1.5, size = 2, hjust = 0, label = sigrho )+      
        annotate(geom = "text", x = 0.02, y = 1.3, size = 2, hjust = 0, label = sigRHO )+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p5 
ggsave("output_files/Haas_Scale_RNAvsPro_update.pdf", p5,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 6: scale of (protein) expression variation vs mRNA-to-protein correlation, scatterplot ####

# This shows that progulons with larger expression variation tend to have a stronger mRNA-based regulatory component

# Calculate significance of correlation and turn into plot labels
sigrhoest <- DT_summary[, cor.test(med_RNA_pro_rho, med_prot_mad, method = "spearman" )$estimate ] 
sigrhopva <- DT_summary[, cor.test(med_RNA_pro_rho, med_prot_mad, method = "spearman" )$p.value  ] 
sigRHOest <- DT_summary[, cor.test(med_RNA_pro_rho, med_prot_mad, method = "spearman")$estimate ]
sigRHOpva <- DT_summary[, cor.test(med_RNA_pro_rho, med_prot_mad, method = "spearman")$p.value  ]

sigrho <- paste("rho", round(sigrhoest, 2), ", p value", signif(sigrhopva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Make the plot
p6 <- ggplot(DT_summary, aes(x = med_RNA_pro_rho, y = med_prot_mad))+
      geom_smooth( method = "lm", size = 0.25, colour = "orange", se = FALSE, fullrange = TRUE)+
      geom_point( alpha = 0.5, size = 1 )+
      geom_text(data = DT_summary[ Progulon_ID %in% c("P04","P05","P25", "P09", "P31","P02","P16","P26")],
                aes(label = prn_function), size = 2, colour = "blue", hjust = -0.1)+
      scale_x_continuous( limits = c(0,0.85), expand = c(0,0), breaks = seq(0,1,0.2))+
      scale_y_continuous( limits = c(0,1.8), expand = c(0,0), breaks = seq(0,2,0.4))+
      xlab("mRNA-to-protein correlation [median rho]")+
      ylab("Scale of protein expression variation [median MAD]")+
      annotate(geom = "text", x = 0.1, y = 1.6, size = 2, hjust = 0, label = sigrho )+      
      annotate(geom = "text", x = 0.1, y = 1.5, size = 2, hjust = 0, label = sigRHO )+      
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
            legend.position = "none", axis.ticks = element_line(size = 0.25))

p6
ggsave("output_files/Haas_Scale_RNAtoPro_update.pdf", p6,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 7: scale of mRNA expression variation vs mRNA-to-mRNA correlation, scatterplot ####

# This shows that progulons with large mRNA expression variation do not necessarily have better mRNA-to-mRNA coordination
# For example, several progulons have larger expression changes than the ribosome but weaker mRNA-mRNA correlation. This 
# would suggest that mRNA-mRNA coordination is independent of the scale of expression changes, but perhaps dependent on 
# biological function (see protein synthesis and degradation related progulons)

# Calculate significance of correlation and turn into plot labels
sigrhoest <- DT_summary[, cor.test(med_RNA_RNA_rho, med_mRNA_mad, method = "spearman" )$estimate ] 
sigrhopva <- DT_summary[, cor.test(med_RNA_RNA_rho, med_mRNA_mad, method = "spearman" )$p.value  ] 
sigRHOest <- DT_summary[, cor.test(med_RNA_RNA_rho, med_mRNA_mad, method = "spearman")$estimate ]
sigRHOpva <- DT_summary[, cor.test(med_RNA_RNA_rho, med_mRNA_mad, method = "spearman")$p.value  ]

sigrho <- paste("rho", round(sigrhoest, 2), ", p value", signif(sigrhopva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Make the plot
p7 <- ggplot(DT_summary, aes(x = med_RNA_RNA_rho, y = med_mRNA_mad))+
      geom_point( alpha = 0.5, size = 1 )+
      geom_text(data = DT_summary[ Progulon_ID %in% c("P04","P05","P25", "P09", "P31","P02","P16","P26")],
                 aes(label = prn_function), size = 2, colour = "blue", hjust = -0.1)+
      scale_x_continuous( limits = c(0,0.5), expand = c(0,0), breaks = seq(0,0.5,0.1))+
      scale_y_continuous( limits = c(0,0.15), expand = c(0,0), breaks = seq(0,0.15,0.05))+
      xlab("median mRNA coexpression [rho]")+
      ylab("Scale of mRNA expression variation [median MAD]")+
      annotate(geom = "text", x = 0.3, y = 0.14, size = 2, hjust = 0, label = sigrho )+      
      annotate(geom = "text", x = 0.3, y = 0.13, size = 2, hjust = 0, label = sigRHO )+      
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
            legend.position = "none", axis.ticks = element_line(size = 0.25))

p7
ggsave("output_files/Haas_Scale_RNAtoRNA_update.pdf", p7,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 8: scale of protein expression variation vs protein-to-protein correlation, scatterplot ####

# This shows that progulons with larger protein expression variation do not necessarily have better protein-to-protein coordination

# Calculate significance of correlation and turn into plot labels
sigrhoest <- DT_summary[, cor.test(med_pro_pro_rho, med_prot_mad, method = "spearman" )$estimate ] 
sigrhopva <- DT_summary[, cor.test(med_pro_pro_rho, med_prot_mad, method = "spearman" )$p.value  ] 
sigRHOest <- DT_summary[, cor.test(med_pro_pro_rho, med_prot_mad, method = "spearman")$estimate ]
sigRHOpva <- DT_summary[, cor.test(med_pro_pro_rho, med_prot_mad, method = "spearman")$p.value  ]

sigrho <- paste("rho", round(sigrhoest, 2), ", p value", signif(sigrhopva, 2))
sigRHO <- paste("rho", round(sigRHOest, 2), ", p value", signif(sigRHOpva, 2))

# Make the plot
p8 <- ggplot(DT_summary, aes(x = med_pro_pro_rho, y = med_prot_mad))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_text(data = DT_summary[ Progulon_ID %in% c("P04","P05","P25", "P09", "P31","P02","P16","P26")],
                  aes(label = prn_function), size = 2, colour = "blue", hjust = -0.1)+
        scale_x_continuous( limits = c(-0.05,0.8), expand = c(0,0), breaks = seq(0,1,0.2))+
        scale_y_continuous( limits = c(0,1.2), expand = c(0,0), breaks = seq(0,1.2,0.2))+
        xlab("median protein coexpression [rho]")+
        ylab("Scale of protein expression variation [median MAD]")+
        annotate(geom = "text", x = 0.5, y = 1.1, size = 2, hjust = 0, label = sigrho )+      
        annotate(geom = "text", x = 0.5, y = 1.0, size = 2, hjust = 0, label = sigRHO )+      
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p8
ggsave("output_files/Haas_Scale_ProToPro_update.pdf", p8,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 9: Composite scale-related plots ####

# To output a properly aligned plot grid, combine the gtables first
g5 <- ggplotGrob(p5)
g6 <- ggplotGrob(p6)
g7 <- ggplotGrob(p7)
g8 <- ggplotGrob(p8)
g <- cbind( rbind(g5, g6, size = "first"), rbind(g7, g8, size = "first"), size = "first")
grid.newpage()
grid.draw(g)

ggsave("output_files/Haas_Scale_composite_update.pdf", g,
       width = 6.7, height = 6.5, units = "cm")


#### Plot 10: Volcano plot showing progulon regulation ####

p10 <- ggplot(DT_summary, aes(x = med_RNA_pro_rho, y = -log10(p_RNAtoPRO)))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_vline( xintercept = DT[ Gene_1 == Gene_2 , median(RNA_pro_rho) ], size = 0.25, linetype = "dashed")+
        geom_text(data = DT_summary[ order(med_RNA_pro_rho) ][ c(1:31) ],
                  aes(label = prn_function), size = 2, colour = "blue", hjust = 1, angle = 90)+
        scale_x_continuous( limits = c(0.1,0.8), expand = c(0,0), breaks = seq(0,0.8,0.1))+
        xlab("mRNA to protein correlation [rho]")+
        ylab("-log10 p-value")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p10 
ggsave("output_files/Haas_volcano_update.pdf", p10,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 11: mRNA-mRNA vs mRNA-to-protein correlation per progulon ####

sigRHO_11 <- paste("rho", round( DT_summary[, cor.test(med_RNA_RNA_rho, med_RNA_pro_rho, method = "spearman")$estimate ], 2),
                   ", p value", signif(DT_summary[, cor.test(med_RNA_RNA_rho, med_RNA_pro_rho, method = "spearman")$p.value  ], 2))


p11 <- ggplot(DT_summary, aes(x = med_RNA_RNA_rho, y = med_RNA_pro_rho))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_smooth( method = "lm", size = 0.25, colour = "orange", se = FALSE, fullrange = TRUE)+
        annotate(geom = "text", x = 0.1, y = 0.7, size = 2, hjust = 0, label = sigRHO_11 )+      
        geom_text(data = DT_summary[ order(med_RNA_pro_rho) ][ c(1:31) ],
                  aes(label = prn_function), size = 2, colour = "blue", hjust = 1, angle = 90)+
        scale_x_continuous( limits = c(0,0.6), expand = c(0,0), breaks = seq(0,0.6,0.2))+
        scale_y_continuous( limits = c(0,0.8), expand = c(0,0), breaks = seq(0,0.8,0.2))+
        ylab("mRNA to protein correlation [rho]")+
        xlab("mRNA coexpression [rho]")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p11
ggsave("output_files/Haas_RNARNAvsmRNAPro_update.pdf", p11,
       width = 4.9, height = 4.9, units = "cm")


#### Plot 12: pro-pro vs mRNA-to-protein correlation per progulon ####

sigRHO_12 <- paste("rho",       round( DT_summary[, cor.test(med_pro_pro_rho, med_RNA_pro_rho, method = "spearman")$estimate ], 2),
                   ", p value", signif(DT_summary[, cor.test(med_pro_pro_rho, med_RNA_pro_rho, method = "spearman")$p.value  ], 2))


p12 <- ggplot(DT_summary, aes(x = med_pro_pro_rho, y = med_RNA_pro_rho))+
        geom_point( alpha = 0.5, size = 1 )+
        geom_smooth( method = "lm", size = 0.25, colour = "orange", se = FALSE, fullrange = TRUE)+
        annotate(geom = "text", x = 0.1, y = 0.7, size = 2, hjust = 0, label = sigRHO_12 )+      
        geom_text(data = DT_summary[ order(med_RNA_pro_rho) ][ c(1:31) ],
                  aes(label = prn_function), size = 2, colour = "blue", hjust = 1, angle = 90)+
        scale_x_continuous( limits = c(-0.05,0.81), expand = c(0,0), breaks = seq(0,1,0.2))+
        scale_y_continuous( limits = c(0,0.8), expand = c(0,0), breaks = seq(0,0.8,0.2))+
        ylab("mRNA to protein correlation [rho]")+
        xlab("protein coexpression [rho]")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", axis.ticks = element_line(size = 0.25))

p12
ggsave("output_files/Haas_propro_vsmRNAPro_update.pdf", p12,
       width = 4.9, height = 4.9, units = "cm")


#### Combined plot ####

# Output combined plot to create manuscript figure
combined_p <- ggarrange(p10, p2, p11, p12, nrow = 2)

ggsave("output_files/Haas_combined_plot_update.pdf", combined_p,
       width = 8.8, height = 8.8, units = "cm")
