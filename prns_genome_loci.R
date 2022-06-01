######################################################################################################## #
# Project: Progulons, with Groth lab     
# Purpose: Analyse how progulon genes are distributed across the genome
# Author: G Kustatscher
# Date: February 2022
######################################################################################################## #

# Read in the necessary libraries
library(data.table); library(ggplot2)

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

# Get a vector with all unique protein IDs that were assigned to a progulon, i.e. omitting proteins
# that are in ProteomeHD but were never assigned to any progulon
all_proteins <- prns[ prot_in_prn == "yes", unique(Protein_IDs) ]

# Simplify that vector by removing isoforms
all_proteins_simpleID <- gsub(";.+", "", all_proteins)
all_proteins_simpleID <- gsub("-.+", "", all_proteins_simpleID)
all_proteins_simpleID <- unique(all_proteins_simpleID)

# Find all pairs of proteins that are in (any) progulon together
prog_pairs <- data.table()                                                         # Initiate result table
for(i in unique(prns$Progulon_ID)){                                                # Loop over each progulon
prog_prot_IDs <- prns[ Progulon_ID == i & prot_in_prn == "yes" , Protein_IDs ]     # These are the proteins associated with the current progulon
prog_prot_IDs <- gsub(";.+", "", prog_prot_IDs)                                    # Simplify protein IDs by removing isoform information 1/2
prog_prot_IDs <- gsub("-.+", "", prog_prot_IDs)                                    # Simplify protein IDs by removing isoform information 2/2
prog_prot_IDs <- unique( prog_prot_IDs )                                           # Remove duplicates arising from isoform info removal
temp_prog_pairs <- as.data.table( t( combn( prog_prot_IDs , 2 )))                  # All possible pairwise combinations of these progulon proteins
temp_prog_pairs <- as.data.table( t( apply(temp_prog_pairs, 1, sort)))             # Sort pairs row-wise alphabetically
temp_prog_pairs <- temp_prog_pairs[, .(Protein_X = V1, Protein_Y = V2) ]           # Rename columns
temp_prog_pairs <- data.table( temp_prog_pairs, Progulon_ID = i )                  # Append ID of progulon
prog_pairs <- rbind( prog_pairs, temp_prog_pairs)                                  # Combine pairs from different progulons
print(i)                                                                           # Print progress
}

prog_pairs[, .N]                                   # Total pairs that are in any progulon together
prog_pairs[, .N, by = .(Protein_X, Protein_Y) ]    # And the unique pairs with number of progulons they are assigned to


#### Calculate gene positions and pairwise genomic distances ####

# Load ENSEMBL genome coordinates of human genes. The file "Human_gene_positions_17Feb22.txt.gz" was downloaded from www.ensembl.org using BioMart.
Genes <- fread("input_files/Human_gene_positions_17Feb22.txt.gz")

# Rename the columns
Genes <- Genes[, .(Ensembl_Gene_ID = `Gene stable ID`,          
                   Gene_start      = `Gene start (bp)`,
                   Gene_end        = `Gene end (bp)`,
                   Chromosome      = `Chromosome/scaffold name`, 
                   Strand,
                   Uniprot_ID      = `UniProtKB/Swiss-Prot ID`)]

# Keep only genes whose proteins were detected in our experiment
Genes <- Genes[ Uniprot_ID %in% all_proteins_simpleID ]

# For proteins encoded by multiple genes (e.g. histones), keep only one gene
Genes <- Genes[ !duplicated(Uniprot_ID) ]

# Calculate the transcription start sites (TSS)
Genes$Gene_TSS <- ifelse(Genes$Strand == 1, Genes$Gene_start, Genes$Gene_end)

# Calculate the transcription end sites
Genes$Gene_end <- ifelse(Genes$Strand == 1, Genes$Gene_end, Genes$Gene_start)

# Create a table with all pairwise gene combinations
gene_IDs <- unique(Genes$Ensembl_Gene_ID)
DT <- as.data.table( t( combn( gene_IDs , 2)))
DT <- DT[, .(Gene_1 = V1, Gene_2 = V2)]

# Assign gene annotation to the pairs
DT <- merge(DT, Genes, by.x="Gene_1", by.y="Ensembl_Gene_ID")                             # Add genomic annotation for Gene 1
DT <- merge(DT, Genes, by.x="Gene_2", by.y="Ensembl_Gene_ID", suffixes = c("_1", "_2"))   # Add genomic annodation for Gene 2

# For gene pairs on same chromosome, calculate distance between their TSSs
DT[ Chromosome_1 == Chromosome_2 , Distance := abs(Gene_TSS_1-Gene_TSS_2) ]   

# Assign gene pairs from bidirectional promoters
DT[ Distance < 1000 &                          # If the TSS distance between Gene 1 and Gene 2 is < 1kb ...
    Strand_1 != Strand_2 &                     # ... and they are on opposite strands
    Distance < abs(Gene_TSS_1-Gene_end_2) &    # ... and the TSS of Gene 1 is further away from the end of Gene 2 than from its TSS
    Distance < abs(Gene_TSS_2-Gene_end_1),     # ... and the TSS of Gene 2 is further away from the end of Gene 1 than from its TSS
  Bidirectional := "bidirectional" ]         # --> then call it a bidirectional gene pair

# Keep relevant columns only
DT <- DT[, .(Gene_1, Gene_2, Uniprot_ID_1, Uniprot_ID_2,
             Chromosome_1, Chromosome_2, Distance, Bidirectional)]  

# Annotate which gene pairs are actually in a (any) progulon together
DT[, Protein_1_sorted := ifelse(Uniprot_ID_1 < Uniprot_ID_2, Uniprot_ID_1, Uniprot_ID_2) ]      # Sort the IDs row-wise
DT[, Protein_2_sorted := ifelse(Uniprot_ID_1 > Uniprot_ID_2, Uniprot_ID_1, Uniprot_ID_2) ]      # Sort the IDs row-wise
DT <- DT[, .(Gene_1, Gene_2, Uniprot_ID_1 = Protein_1_sorted, Uniprot_ID_2 = Protein_2_sorted,  # Keep relevant columns
             Chromosome_1, Chromosome_2, Distance, Bidirectional)]

DT <- merge(DT, prog_pairs,                               # Merge genomic annotation with progulon annotation
            by.x = c("Uniprot_ID_1", "Uniprot_ID_2"),     # Note that pairs have been sorted
            by.y = c("Protein_X", "Protein_Y"),           # Note that pairs have been sorted
            all.x = TRUE)                                 # Keep pairs that are not in any progulon


#### Write out gene IDs belonging to the ATP synthase and Nucleosome progulons ####

P09_ATPs_genes <- DT[ Progulon_ID == "P09" , .( unique( c( Gene_1, Gene_2 ))) ]
 P31_Nuc_genes <- DT[ Progulon_ID == "P31" ,  .( unique( c( Gene_1, Gene_2 ))) ]

fwrite(P09_ATPs_genes, "output_files/P09_ATPs_gene_loci.csv")
fwrite(P31_Nuc_genes, "output_files/P31_Nuc_gene_loci.csv")

# Clear workspace
rm( list = ls()[! ls() %in% c("DT", "prns")] )   


#### Analyse genomic distribution of gene pairs: chromosomal co-localisation ####

## Question: Are gene pairs from the same progulon enriched for genes from the same chromosome?
## (analyse on a per-progulon basis because the overall calculation is biased by progulon size)

# First the get the numbers that would be expected by random chance:

N_gp_all_sChr <- DT[ Chromosome_1 == Chromosome_2 , .(Uniprot_ID_1, Uniprot_ID_2) ]  # Overall gene pairs from same chromosome
N_gp_all_sChr <- unique( N_gp_all_sChr )
N_gp_all_sChr <- N_gp_all_sChr[, .N]

N_gp_all_NsChr <- DT[ Chromosome_1 != Chromosome_2 , .(Uniprot_ID_1, Uniprot_ID_2) ] # Overall gene pairs not on the same chromosome
N_gp_all_NsChr <- unique( N_gp_all_NsChr )
N_gp_all_NsChr <- N_gp_all_NsChr[, .N]

# Enrichment in percent
pc_random_sChr <- N_gp_all_sChr / (N_gp_all_sChr + N_gp_all_NsChr)   * 100

# Then get the numbers for each progulon:
progulons_chr_enrichment <- data.table()   # Initialise result table

for( i in unique(prns$Progulon_ID)){       # Loop through all progulons

  # Number of gene pairs in the progulon that are on the same chromosome
  N_gp_sp_sChr <- DT[ Progulon_ID == i & Chromosome_1 == Chromosome_2 , .(Uniprot_ID_1, Uniprot_ID_2) ]
  N_gp_sp_sChr <- unique( N_gp_sp_sChr )
  N_gp_sp_sChr <- N_gp_sp_sChr[, .N]
  
  # Number of gene pairs in the progulon that are NOT on the same chromosome
  N_gp_sp_NsChr <- DT[ Progulon_ID == i & Chromosome_1 != Chromosome_2 , .(Uniprot_ID_1, Uniprot_ID_2) ]
  N_gp_sp_NsChr <- unique( N_gp_sp_NsChr )
  N_gp_sp_NsChr <- N_gp_sp_NsChr[, .N]
  
  # Create the contingency table
  prog_chr <- data.frame(in_same_prog = c( N_gp_sp_sChr,  N_gp_sp_NsChr),
                               random = c( N_gp_all_sChr, N_gp_all_NsChr),
                            row.names = c( "on_same_chr", "not_on_same_chr" ))
  
  # Perform Fisher's Exact test
  pvalue <- fisher.test(prog_chr)$p.value
 
  # Enrichment in percent
  pc_gp_sp_sChr <- N_gp_sp_sChr / (N_gp_sp_sChr + N_gp_sp_NsChr)   * 100
  
  # Combine output
  temp_DT <- data.table(Progulon_ID = i, pvalue, pc_gp_sp_sChr)
  progulons_chr_enrichment <- rbind(progulons_chr_enrichment, temp_DT)
  print(i)
}


##### Create the plot ####


p <- ggplot(progulons_chr_enrichment, aes(x = pc_gp_sp_sChr, y = -log10(pvalue)))+
      geom_vline(xintercept = pc_random_sChr, size = 0.25, linetype = "dashed", colour = "grey50")+
      geom_point(alpha = 0.5, size = 1)+
      geom_point(data = progulons_chr_enrichment[ Progulon_ID %in% c("P09", "P31")], size = 1, colour = "green")+
      geom_text(data = progulons_chr_enrichment[ Progulon_ID %in% c("P09", "P31")], aes(label = Progulon_ID), size = 2, hjust = -0.3)+
      scale_x_continuous(limits = c(0,25), breaks = c(0, 12.5, 25))+
      xlab("% gene pairs\non same chromosome")+
      ylab("-log10 p-value")+
      theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
            legend.position = "top", axis.ticks = element_line(size = 0.25))

p

ggsave("output_files/Gene_loci_volcano.pdf", p, width = 5, height = 5, units = "cm")


