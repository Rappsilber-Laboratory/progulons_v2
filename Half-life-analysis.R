######################################################################################################## #
# Project: Progulons, with Groth lab     
# Purpose: Analyse and compare mRNA and protein half-lives of progulon genes
# Author: E Fiagbedzi and G Kustatscher
# Date: April 2022
######################################################################################################## #

# Load necessary libraries
library(data.table); library(readxl); library(ggplot2); library(egg)


#### Get Progulon associations ####

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

# Simplify protein IDs by removing multiple Uniprot IDs and isoform information
prns[, SimpleID := gsub(";.+", "", Protein_IDs) ][, SimpleID := gsub("-.+", "", SimpleID)]


#### Prep mRNA half-lives from Tani et al (Genome Research 2012) ####

# Read in mRNA half-lifes (hr) from Table S1 from Tani et al (Genome Research 2012), which is available here:
# https://genome.cshlp.org/content/suppl/2012/02/14/gr.130559.111.DC1/Tani_Supp_Tables_revised2.xls
# Ignore warnings - they are due to this being an Excel file 
mRNA_HL <- read_xls("input_files/Tani_Supp_Tables_revised2.xls", sheet = 1, skip = 3)

# Convert to data.table for easier handling
mRNA_HL <- as.data.table(mRNA_HL)

# Keep only mRNAs with a half-life measurement (i.e. that don't have an NA half-life)
mRNA_HL <- mRNA_HL[ !is.na(`t1/2 (h)`) ]           

# Remove commas from mRNA name (this will fuse multi-name entries, but we would discard these anyway)
mRNA_HL[, RepName := gsub(",", "", RepName) ]  

# Convert Refseq mRNA IDs to UniProt IDs
fwrite( mRNA_HL[, .(RepName)], "output_files/temp_mRNA_IDs.csv" )  # (1) Write them out as a temporary file
                                                      # (2) Upload this file to UniProts Retrieve function (https://www.uniprot.org/uploadlists/)
                                                      # (3) Select to convert from RefSeq Nucleotide to UniprotKB
                                                      # (4) Download the result as tab-separated file 
converted <- fread("input_files/uniprot-yourlist.tab.gz")         # (5) Load the file back in here
setnames(converted, old = grep("yourlist", names(converted)), new = "RefSeqID")  # (6) Convert the name of the RefSeqID column, which Uniprot calls yourlist...
setnames(converted, old = "Entry",                            new = "UniprotID") # (7) Convert the name of the UniProtID column, which Uniprot calls Entry
converted <- converted[, .(UniprotID, RefSeqID) ]     # (8) Keep only these two ID columns
converted <- converted[ !duplicated(UniprotID) ]      # (9) Remove duplicate UniProt IDs, if any
converted <- converted[ !duplicated(RefSeqID) ]       # (10) Remove duplicate RefSeq IDs, if any
mRNA_HL$UniprotID <- converted$UniprotID[ match(mRNA_HL$RepName, converted$RefSeqID) ]   # (11) Assign Uniprot IDs to the RefSeq IDs in the half-life table
mRNA_HL <- mRNA_HL[ complete.cases(mRNA_HL) ]         # (12) Remove genes which haven't been mapped

# Rename columns and remove unnecessary columns
mRNA_HL <- mRNA_HL[, .(UniprotID, mRNA_HL = `t1/2 (h)`) ]

# Restrict to genes that have been mapped to at least one progulon, in order to avoid any bias that could
# derive from a potential bias in assigning proteins to progulons
all_prn_prots <- prns[ prot_in_prn == "yes" , unique(SimpleID) ]
mRNA_HL <- mRNA_HL[ UniprotID %in% all_prn_prots ]


#### Prep protein half-lives from McShane et al (Cell 2016) ####

# Read in protein half-lifes (hr) from Table S4 from McShane et al (Cell 2016), available here:
# https://www.cell.com/cms/10.1016/j.cell.2016.09.015/attachment/ed36d0dc-573c-4f65-ba99-6ce9713d2584/mmc4.xlsx
prot_HL <- read_xlsx("input_files/mmc4.xlsx", sheet = 1)   # Read in downloaded file (ignore warning related to Excel issue)
prot_HL <- as.data.table(prot_HL)                     # Convert to data.table for easier handling

# Select and rename relevant columns
prot_HL <- prot_HL[, .(UniprotID = `Protein IDs (Uniprot)`, 
                       protein_HL = `Half-life (exponential) 1-state-model [h]`) ]

# Simplify Uniprot IDs (e.g. remove isoform information)
prot_HL[, UniprotID := gsub(";.+", "", UniprotID) ][, UniprotID := gsub("-.+", "", UniprotID)]

# Remove duplicate genes
prot_HL <- prot_HL[ !duplicated(UniprotID) ]

# Restrict to genes that have been mapped to at least one progulon, in order to avoid any bias that could
# derive from a potential bias in assigning proteins to progulons
prot_HL <- prot_HL[ UniprotID %in% all_prn_prots ]

# Half-lives of > 300 hours are not accurate, so I remove them
prot_HL[ protein_HL == "> 300",  protein_HL := NA ]

# Turn HL into numeric vector
prot_HL[, protein_HL := as.numeric(  protein_HL ) ]

# Keep only relevant columns
prot_HL <- prot_HL[, .(UniprotID, protein_HL)]

# Remove proteins with missing half-lives
prot_HL <- prot_HL[ complete.cases(prot_HL) ]

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "mRNA_HL", "prot_HL")] )  


#### Compare mRNA and protein half-lives on a per-gene basis ####

# Merge mRNA and protein HL data
HL <- merge(mRNA_HL, prot_HL, by = "UniprotID")

# Calculate Spearman's rho, including the statistical significance of the correlation between mRNA and protein HLs
mRNA_prot_RHO <- HL[, cor.test(mRNA_HL, protein_HL, method = "spearman" ) ] 

# Extract correlation coefficient and p-value and turn into plot labels
RHO_label <- paste("RHO ", 
                   round( mRNA_prot_RHO$estimate , 2),
                   "\n p-value ", 
                   signif( mRNA_prot_RHO$p.value , 2), 
                   sep = "")

# Set plotting theme
theme_set( theme( text = element_text( size = 7), line = element_line( size = 0.25),
                  panel.border = element_rect(colour = "black", size = 0.25, fill = NA), panel.background = element_blank(), panel.grid = element_blank(),
                  axis.text = element_text( colour = "black")
                  ))

# Create the plot
pHL_genes <- ggplot(HL, aes(mRNA_HL, protein_HL))+
             geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
             geom_point( alpha = 0.5, size = 0.75, shape = 16)+
             scale_x_continuous( limits = c(0,25), expand = c(0,0), breaks = seq(0,25,5))+
             scale_y_continuous( limits = c(0,300), expand = c(0,0), breaks = seq(0,300,50))+
             xlab("mRNA half-life [h] - HeLa cells")+
             ylab("protein half-life [h] - RPE1 cells")+
             annotate(geom = "text", x = 15, y = 250, size = 2, hjust = 0.5, label = RHO_label)


#### Prep mRNA half-lives per progulon ####

# This section calculates the average half-life of mRNAs in each progulon, and tests if the half-lives of genes within a progulon 
# are statistically different from genes that are not in the same progulon

HL_mRNA_prns <- data.table()                                                # Initialise result table
for(i in unique(prns$Progulon_ID) ){                                        # Loop through all progulons
   prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Capture proteins of the current progulon
      in_prn <- mRNA_HL[  UniprotID %in% prn_prots , mRNA_HL ]              # Retrieve their half-lives
  not_in_prn <- mRNA_HL[ !UniprotID %in% prn_prots , mRNA_HL ]              # Retrieve the half-lives of the non-progulon proteins
      tmp_dt <- data.table(                                                 # Aggregate the results into a temporary table
        Progulon_ID = i,                                        # ID of the current progulon
        mean_mRNA_HL = mean(in_prn),                            # Average mRNA HL of the current progulon
        N_mRNA_HL = length(in_prn),                             # Number of mRNAs HLs on which this assessment is based
        pval_mRNA_HL = wilcox.test(in_prn, not_in_prn)$p.value  # Two-sided p-value of a Mann-Whitney test - are the mRNA HLs of the progulon members different from the other proteins?
        )
  HL_mRNA_prns <- rbind(HL_mRNA_prns, tmp_dt)                               # Combine results of each progulon into the shared output table
}

# Since we test multiple progulons, also apply a p-value adjustment
HL_mRNA_prns[, BH_adj_pval_mRNA_HL := p.adjust(pval_mRNA_HL, method = "BH") ]


#### Prep protein half-lives per progulon ####

HL_prot_prns <- data.table()                                               # Initialise result table
for(i in unique(prns$Progulon_ID) ){                                       # Loop through all progulons
  prn_prots <- prns[ Progulon_ID == i & prot_in_prn == "yes" , SimpleID ]  # Capture proteins of the current progulon
     in_prn <- prot_HL[  UniprotID %in% prn_prots , protein_HL ]           # Retrieve their half-lives
 not_in_prn <- prot_HL[ !UniprotID %in% prn_prots , protein_HL ]           # Retrieve the half-lives of the non-progulon proteins
     tmp_dt <- data.table(                                                 # Aggregate the results into a temporary table
       Progulon_ID = i,                                        # ID of the current progulon
       mean_prot_HL = mean(in_prn),                            # Average prot HL of the current progulon
       N_prot_HL = length(in_prn),                             # Number of prots HLs on which this assessment is based
       pval_prot_HL = wilcox.test(in_prn, not_in_prn)$p.value  # Two-sided p-value of a Mann-Whitney test - are the prot HLs of the progulon members different from the other proteins?
       )
  HL_prot_prns <- rbind(HL_prot_prns, tmp_dt)                              # Combine results of each progulon into the shared output table
}

# Since we test multiple progulons, also apply a p-value adjustment
HL_prot_prns[, BH_adj_pval_prot_HL := p.adjust(pval_prot_HL, method = "BH") ]


#### Compare mRNA and protein half-lives on a per-progulon basis ####

# Merge mRNA and protein HL data
HL_prns <- merge(HL_mRNA_prns, HL_prot_prns, by = "Progulon_ID")

# Add progulon annotation
annot <- fread("input_files/Progulon_annotation.csv")                           # Load a table with manually created progulon annotation
annot <- annot[, .(Progulon_ID, Progulon_function = Progulon_name)]  # Keep and rename relevant columns
HL_prns <- merge(HL_prns, annot, by = "Progulon_ID")

# Tidy up workspace (this will remove everything except the listed objects)
rm( list = ls()[! ls() %in% c("prns", "mRNA_HL", "prot_HL", "HL_prns", "pHL_genes")] )  

# Calculate Spearman's rho, including the statistical significance of the correlation between mRNA and protein HLs - per progulon
mRNA_prot_RHO_prn <- HL_prns[, cor.test(mean_mRNA_HL, mean_prot_HL, method = "spearman" ) ] 

# Extract correlation coefficient and p-value and turn into plot labels
RHO_label_prn <- paste("RHO ", 
                       round( mRNA_prot_RHO_prn$estimate , 2),
                       "\n p-value ", 
                       signif( mRNA_prot_RHO_prn$p.value , 2), 
                       sep = "")

# Create the plot
pHL_prns <- ggplot(HL_prns, aes(mean_mRNA_HL, mean_prot_HL))+
            geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
            # geom_vline( xintercept = median(HL_prns$mean_mRNA_HL), size = 0.25, linetype = "dashed", colour = "grey50")+
            # geom_hline( yintercept = median(HL_prns$mean_prot_HL), size = 0.25, linetype = "dashed", colour = "grey50")+
            geom_point( alpha = 0.75, size = 1, shape = 16)+
            geom_text(aes(label = Progulon_function), size = 2, nudge_y = -2)+
            annotate(geom = "text", x = 10, y = 80, size = 2, hjust = 0.5, label = RHO_label_prn)+
            scale_x_continuous( limits = c(5.5,12.5), expand = c(0,0), breaks = seq(0,20,2))+
            scale_y_continuous( limits = c(18,88), expand = c(0,0), breaks = seq(0,100,20))+
            xlab("mean mRNA half-life [h]")+
            ylab("mean protein half-life [h]")


#### Calculate coordination of protein and mRNA half-lives of proteins in progulons ####

# Calculate the coefficient of variation of mRNA half lives within progulons
prog_corr_mRNA <- data.table()                                            # Initialise result table
for (id in unique(prns$Progulon_ID)){                                     # Loop through all progulons
  current_prn <- prns[Progulon_ID == id & prot_in_prn == "yes", SimpleID] # Capture proteins of the current progulon
  current_prn_mRNA_HL <- mRNA_HL[UniprotID %in% current_prn]              # Retrieve their Half-lives
  prn_mRNA_HL_cv <- current_prn_mRNA_HL[,sd(mRNA_HL)/mean(mRNA_HL)]       # Calculate the coefficient of variation of the half-lives
  tmp_dt <- data.table(                                                   # Aggregate the results into a data.table
    Progulon_ID = id,              # Progulon Ids
    mRNA_HL_cv = prn_mRNA_HL_cv    # Coefficient of variation of half-lives
  )
  prog_corr_mRNA <- rbind(prog_corr_mRNA, tmp_dt)                         # Combine results of each progulon into the shared data.table
}

# Calculate the coeficient of variation of protein half-lives within progulons
prog_corr_prot <- data.table()                                            #Initialise result table
for (id in unique(prns$Progulon_ID)){                                     # Loop through all progulons
  current_prn <- prns[Progulon_ID == id & prot_in_prn == "yes", SimpleID] # Capture proteins of the current progulon
  current_prn_prot_HL <-prot_HL[UniprotID %in% current_prn]               # Retrieve their Half-lives
  prn_prot_HL_cv <- current_prn_prot_HL[,sd(protein_HL)/mean(protein_HL)] # Calculate the coefficient of variation of the half-lives
  tmp_dt <- data.table(                                                   # Aggregate the results into a data.table
    Progulon_ID = id,             # Progulon Ids
    prot_HL_cv = prn_prot_HL_cv   # coefficient of variation of half-lives
  )
  prog_corr_prot <- rbind(prog_corr_prot, tmp_dt)                         # Combine results of each progulon into the shared data.table
}

# Merge the prog correlation data with the functional annotations
mRNA_prns_with_functions    <- merge( prog_corr_mRNA, HL_prns[, .(Progulon_ID, Progulon_function)], by = "Progulon_ID") # for mRNA
protein_prns_with_functions <- merge( prog_corr_prot, HL_prns[, .(Progulon_ID, Progulon_function)], by = "Progulon_ID") # for protein

# merge protein prns table and mRNA prns table
mRNA_prot_prns_with_functions <- merge( 
  mRNA_prns_with_functions,
  protein_prns_with_functions,
  by = c("Progulon_ID", "Progulon_function"))


#### Compare progulon HL, HL coordination and steady-state coordination ####

# Load table containing median RNA-RNA correlation and median protein-protein correlation per progulon
median_corr_data <- fread("input_files/RNA_protein_correlations_breast_cancer.csv")

# convert to data.table object
median_corr_data <- as.data.table(median_corr_data)

# merge median_corr_data with mRNA_prot_prns_with_functions
corr_plot_data <- merge(
  median_corr_data[, .(Progulon_ID, med_RNA_RNA_rho, med_pro_pro_rho, med_RNA_pro_rho)],
  mRNA_prot_prns_with_functions,
  by = "Progulon_ID")

# Merge with mean half-live data
corr_plot_data <- merge( HL_prns[, .(Progulon_ID, mean_mRNA_HL, mean_prot_HL)], corr_plot_data, by = "Progulon_ID" )

# Create function to format correlations
f_rho_label <- function(x){ paste( "RHO ", 
                                   round(x$estimate, 2), 
                                   "\n p-value ", 
                                   signif(x$p.value, 2),
                                   sep="") }

# Calculate the set of Spearman correlations for mRNA level
  mRNA_HL_cv_vs_mRNA_rho_corr <- corr_plot_data[, cor.test(  mRNA_HL_cv, med_RNA_RNA_rho, method = "spearman")]
mean_mRNA_HL_vs_mRNA_rho_corr <- corr_plot_data[, cor.test(mean_mRNA_HL, med_RNA_RNA_rho, method = "spearman")]

# Calculate the set of Spearman correlations for prot level
  prot_HL_cv_vs_prot_rho_corr <- corr_plot_data[, cor.test(  prot_HL_cv, med_pro_pro_rho, method = "spearman")]
mean_prot_HL_vs_prot_rho_corr <- corr_plot_data[, cor.test(mean_prot_HL, med_pro_pro_rho, method = "spearman")]

# make plot of mRNA_HL_cv vs mRNA_mRNA_correlation
mRNA_HL_cv_vs_mRNA_rho <- ggplot(corr_plot_data, aes(x=mRNA_HL_cv, y=med_RNA_RNA_rho))+
  geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
  geom_point(alpha = 0.75, size = 1, shape = 16)+
  annotate(geom = "text", x = 0.5, y = 0.4, size = 2, hjust = 0.5, label = f_rho_label(mRNA_HL_cv_vs_mRNA_rho_corr))+
  geom_text(data = corr_plot_data[ med_RNA_RNA_rho > 0.4],   # Add Labels for these outliers 
            aes(label = Progulon_function), size = 2, nudge_y = -0.05)+
  xlab("mRNA half-life coordination\n(CV values; lower is better)")+
  ylab("mRNA steady-state coordination\n(accross breast cancer cell lines)")

# make plot of mean_mRNA_HL vs mRNA_mRNA_correlation
mRNA_HL_vs_mRNA_rho <- ggplot(corr_plot_data, aes(x=mean_mRNA_HL, y=med_RNA_RNA_rho))+
  geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
  geom_point(alpha = 0.75, size = 1, shape = 16)+
  annotate(geom = "text", x = 6, y = 0.4, size = 2, hjust = 0.5, label = f_rho_label(mean_mRNA_HL_vs_mRNA_rho_corr))+
  geom_text(data = corr_plot_data[ med_RNA_RNA_rho > 0.4],   # Add Labels for these outliers 
            aes(label = Progulon_function), size = 2, nudge_y = -0.05)+
  xlab("mean mRNA half-life [h]")+
  ylab("mRNA steady-state coordination\n(accross breast cancer cell lines)")

# Plot prot_HL_cv vs protein_protein_rho
prot_HL_cv_vs_prot_rho <- ggplot(corr_plot_data, aes(x = prot_HL_cv, y = med_pro_pro_rho))+
  geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
  geom_point(alpha = 0.75, size = 1, shape = 16)+
  annotate(geom = "text", x = 1.2, y = 0.6, size = 2, hjust = 0.5, label = f_rho_label(prot_HL_cv_vs_prot_rho_corr))+
  geom_text(data = corr_plot_data[ prot_HL_cv < 0.6 | med_pro_pro_rho > 0.5],   # Add Labels for these outliers 
            aes(label = Progulon_function), size = 2, nudge_y = -0.05)+
  xlab("Protein half-life coordination\n(CV values; lower is better)")+
  ylab("Protein steady-state coordination\n(accross breast cancer cell lines)")

# Plot prot_HL_cv vs protein_protein_rho
prot_HL_vs_prot_rho <- ggplot(corr_plot_data, aes(x = mean_prot_HL, y = med_pro_pro_rho))+
  geom_smooth(method = "lm", size = 0.25, fill = "grey70")+
  geom_point(alpha = 0.75, size = 1, shape = 16)+
  annotate(geom = "text", x = 40, y = 0.6, size = 2, hjust = 0.5, label = f_rho_label(mean_prot_HL_vs_prot_rho_corr))+
  geom_text(data = corr_plot_data[ prot_HL_cv < 0.6 | med_pro_pro_rho > 0.5],   # Add Labels for these outliers 
            aes(label = Progulon_function), size = 2, nudge_y = -0.05)+
  xlab("mean protein half-life [h]")+
  ylab("Protein steady-state coordination\n(accross breast cancer cell lines)")


#### Combine and output the plots and the data ####

# Align the two plots plots
pHL <- ggarrange( pHL_genes, pHL_prns,
                  mRNA_HL_vs_mRNA_rho, mRNA_HL_cv_vs_mRNA_rho, 
                  prot_HL_vs_prot_rho, prot_HL_cv_vs_prot_rho, 
                  nrow = 3)

# Save the combined plot
ggsave("output_files/mRNA_vs_protein_HL.pdf", pHL, width = 12, height = 18, units = "cm")


# Write out the data for use as a supplementary table
fwrite(corr_plot_data[, .(Progulon_ID, mean_mRNA_HL, mean_prot_HL, mRNA_HL_cv, prot_HL_cv, Progulon_function)], 
       "output_files/mRNA_protein_HL_results.csv")








