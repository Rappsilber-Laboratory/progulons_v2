############################################################################################################# #
# Project: Progulon manuscript    
# Purpose: Plot two progulons with different RNA-to-protein behaviour as examples for Figure 2
# Authors: G Kustatscher
# Last update: March 2022
############################################################################################################# #

# Load libraries
library(ggplot2); library(data.table); library(perm); library(readxl); library(grid); library(gridExtra)

# Part 1: Schematic (mock) figure

#### Centre part of the schematic drawing (lineplots) ####

# Create 10 mock mRNA ratios (mr) and 10 similar mock protein ratios (mp)
mr <- c(0.0, 1.2, -0.6, -0.5, 0.8, -0.7, -0.9, 0.9, -0.2, 0.6)
mp <- c(0.1, 0.8, -0.6, -0.2, 1.2, 0.1, -0.9, 0.9, -0.2, 0.6)

# Create 8 mock mRNAs and proteins by modifying the ratios randomly
set.seed(123)

m_prot <- list()
for(i in 1:8){ 
  m_prot[[i]] <- mp + rnorm(10, 0, 0.25)
  }
m_prot <- as.data.table(m_prot)

m_mRNA <- list()
for(i in 1:8){ 
  m_mRNA[[i]] <- mr + rnorm(10, 0, 0.42)
}
m_mRNA <- as.data.table(m_mRNA)

# Create plotting table
m_prot$experiment <- 1:10
m_mRNA$experiment <- 1:10

mm_prot <- melt( m_prot, id.vars = "experiment")
mm_mRNA <- melt( m_mRNA, id.vars = "experiment")

mm_prot$type <- "protein"
mm_mRNA$type <- "mRNA"

plot_dt <- rbind(mm_prot, mm_mRNA)
plot_dt$type <- factor( plot_dt$type, levels = c("protein", "mRNA"))

# Create the plot
p1 <- ggplot(plot_dt, aes( x = experiment, y = value, group = variable))+
              facet_wrap( ~ type , nrow = 2 , scales = "free_y")+
              scale_x_continuous(breaks = 1:10)+
              xlab("Experiment")+
              ylab("fold-change")+
              geom_line(size = 0.25, alpha = 0.5)+
              theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    axis.title = element_text(size = 6), axis.text.y = element_blank(),
                    legend.position = "none", strip.background = element_blank(),
                    strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25),
                    axis.text.x = element_text( size = 5),
                    axis.ticks.y = element_blank())

ggsave("output_files/Mock1.pdf", p1, width = 6, height = 6, units = "cm")


#### Left part of the schematic drawing (coexpression) ####

# For proteins
p2a <- ggplot( m_prot, aes( x = V1, y = V2 ))+
        geom_point( colour = "#EC008C", alpha = 0.7)+
        xlim(-1, 1)+
        ylim(-1, 1)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "#EC008C", alpha = 0.7)+
        xlab("Gene 1")+
        ylab("Gene 2")+
        annotate( geom = "text", label = paste( "RHO", round( cor(m_prot$V1, m_prot$V2, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "#EC008C", alpha = 0.7)+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 6), axis.text = element_blank(),
          legend.position = "none", strip.background = element_blank(),
          strip.text = element_text(size = 6), axis.ticks = element_blank())

p2b <- ggplot( m_prot, aes( x = V7, y = V8 ))+
        geom_point( colour = "#EC008C", alpha = 0.7)+
        xlim(-1, 1)+
        ylim(-1, 1)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "#EC008C", alpha = 0.7)+
        xlab("Gene X")+
        ylab("Gene Y")+
        annotate( geom = "text", label = paste( "RHO", round( cor(m_prot$V7, m_prot$V8, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "#EC008C", alpha = 0.7)+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_blank())


# For mRNAs
p3a <- ggplot( m_mRNA, aes( x = V1, y = V2 ))+
  geom_point( colour = "#00A0BB", alpha = 0.7)+
  xlim(-1.4, 1.4)+
  ylim(-1.4, 1.4)+
  geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "#00A0BB", alpha = 0.7)+
  xlab("Gene 1")+
  ylab("Gene 2")+
  annotate( geom = "text", label = paste( "RHO", round( cor(m_mRNA$V1, m_mRNA$V2, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "#00A0BB", alpha = 0.7)+
  theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6), axis.text = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(size = 6), axis.ticks = element_blank())

p3b <- ggplot( m_mRNA, aes( x = V7, y = V8 ))+
  geom_point( colour = "#00A0BB", alpha = 0.7)+
  xlim(-1.4, 1.4)+
  ylim(-1.4, 1.4)+
  geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "#00A0BB", alpha = 0.7)+
  xlab("Gene X")+
  ylab("Gene Y")+
  annotate( geom = "text", label = paste( "RHO", round( cor(m_mRNA$V7, m_mRNA$V8, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "#00A0BB", alpha = 0.7)+
  theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 6), axis.text = element_blank(),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(size = 6), axis.ticks = element_blank())

p23 <- arrangeGrob(p2a, p2b, p3a, p3b, nrow = 2)

ggsave("output_files/Mock2.pdf", p23, width = 6, height = 6, units = "cm")


#### Left part of the schematic drawing (histograms) ####

# To create mock histograms, repeat the gene creation but make 800, not 8, different genes, otherwise
# the histogramms don't look nice

# Create 800 mock mRNAs and proteins by modifying the ratios randomly
m_prot800 <- list()
for(i in 1:800){ 
  m_prot800[[i]] <- mp + rnorm(10, 0, 0.35)
}
m_prot800 <- as.data.table(m_prot800)

m_mRNA800 <- list()
for(i in 1:800){ 
  m_mRNA800[[i]] <- mr + rnorm(10, 0, 0.5)
}
m_mRNA800 <- as.data.table(m_mRNA800)

# Get all pairwise protein correlations
prot_cor <- cor(m_prot800, method = "spearman")
prot_cor <- as.data.table( melt(prot_cor))
prot_cor <- prot_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
prot_cor$type <- "protein"

# Get all pairwise mRNA correlations
mRNA_cor <- cor(m_mRNA800, method = "spearman")
mRNA_cor <- as.data.table( melt(mRNA_cor))
mRNA_cor <- mRNA_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
mRNA_cor$type <- "mRNA"

# Plot the mock histograms
plot_dt_hist <- rbind(prot_cor, mRNA_cor)

pH <- ggplot( plot_dt_hist , aes( x = value , fill = type ))+
        scale_fill_manual( values = c("#00A0BB", "#EC008C"))+
        geom_histogram(position = "identity", alpha = 0.7, binwidth = 0.02, center = 0.01)+
        geom_vline( aes( xintercept = median( prot_cor$value )) , linetype = "dotted", size = 0.25, colour = "#EC008C")+
        geom_vline( aes( xintercept = median( mRNA_cor$value )) , linetype = "dotted", size = 0.25, colour = "#00A0BB")+
        scale_x_continuous( limits = c(0,1))+
        xlab("Gene coexpression\n[RHO]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text.x = element_text(size = 5, colour = "black"),
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks.x = element_line(size = 0.25))


#### Right part of the schematic drawing (mRNA contribution to protein changes) ####

plot_dt <- data.table( protein = m_prot$V1, mRNA = m_mRNA$V1 )
p5a <- ggplot( plot_dt, aes( y = protein, x = mRNA ))+
        geom_point()+
        xlim(-1.4, 1.4)+
        ylim(-1.4, 1.4)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "grey50")+
        xlab("mRNA")+
        ylab("protein")+
        annotate( geom = "text", label = paste( "RHO", round( cor(plot_dt$protein, plot_dt$mRNA, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "grey50")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_blank())

plot_dt <- data.table( protein = m_prot$V2, mRNA = m_mRNA$V2 )
p5b <- ggplot( plot_dt, aes( y = protein, x = mRNA ))+
        geom_point()+
        xlim(-1.2, 1.6)+
        ylim(-1.2, 1.6)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "grey50")+
        xlab("mRNA")+
        ylab("protein")+
        annotate( geom = "text", label = paste( "RHO", round( cor(plot_dt$protein, plot_dt$mRNA, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "grey50")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_blank())

plot_dt <- data.table( protein = m_prot$V3, mRNA = m_mRNA$V3 )
p5c <- ggplot( plot_dt, aes( y = protein, x = mRNA ))+
        geom_point()+
        xlim(-1.5, 1.6)+
        ylim(-1.5, 1.6)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "grey50")+
        xlab("mRNA")+
        ylab("protein")+
        annotate( geom = "text", label = paste( "RHO", round( cor(plot_dt$protein, plot_dt$mRNA, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "grey50")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_blank())

plot_dt <- data.table( protein = m_prot$V4, mRNA = m_mRNA$V4 )
p5d <- ggplot( plot_dt, aes( y = protein, x = mRNA ))+
        geom_point()+
        xlim(-1.5, 1.6)+
        ylim(-1.5, 1.6)+
        geom_smooth( method = "lm", size = 0.25, se = FALSE, colour = "grey50")+
        xlab("mRNA")+
        ylab("protein")+
        annotate( geom = "text", label = paste( "RHO", round( cor(plot_dt$protein, plot_dt$mRNA, method = "spearman"), 2)), x = -0.9, y = 0.9 , colour = "grey50")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_blank())

p5 <- arrangeGrob( p5a, p5b, p5c, p5d, nrow = 2)

ggsave("output_files/Mock3.pdf", p5, width = 6, height = 6, units = "cm")


#### Right part of the schematic drawing (histograms) ####

# Get the 800-gene-strong mock dataset and find mRNA-protein correlations
crosscor <- cor(m_prot800, m_mRNA800, method = "spearman")
crosscor <- as.data.table( melt(crosscor))
crosscor <- crosscor[ as.character(Var1) == as.character(Var2) ]   # Keep only within-gene (self) correlations

# Plot the mock histograms
p6 <- ggplot( crosscor , aes(x = value ))+
        geom_histogram(position = "identity", binwidth = 0.05, center = 0.025, fill = "grey50")+
        geom_vline( aes( xintercept = median( crosscor$value )) , linetype = "dotted", size = 0.25)+
        scale_x_continuous( limits = c(0,1))+
        xlab("mRNA - protein correlation\n[RHO]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text.x = element_text(size = 5, colour = "black"),
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks.x = element_line(size = 0.25))

# Combine the other histogram
p6 <- arrangeGrob(p6, pH, nrow = 1)

ggsave("output_files/Mock4.pdf", p6, width = 12, height = 6, units = "cm")


#### Part 2: Real examples ##########################################################


#### Get Progulon associations ####

# Clear workspace
rm( list = ls())  

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

# Load the breast cancer dataset from Haas et al (Table S3, Nature Biotechnology, 2017)
pro <- read_xlsx("input_files/Supplementary_Table_3_Protein_mRNA_Profiles_36_Cell_Lines.xlsx", sheet = "Proteome_Profiles")
rna <- read_xlsx("input_files/Supplementary_Table_3_Protein_mRNA_Profiles_36_Cell_Lines.xlsx", sheet = "mRNA_Profiles")

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
write.csv(df$`Uniprot Gene Name`, "~/Desktop/temp_Uniprot_IDs.csv")
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


#### Ribosome progulon (P25) - lineplot ####

# Get the progulon protein IDs
P25_IDs <- prns[ Progulon_ID == "P25" & prot_in_prn == "yes" , unique(SimpleID) ]

# Get and combine the mRNA and protein expression data of these proteins
P25_rna <- melt( trna_mn[, P25_IDs] )
P25_rna$Var1 <- gsub("rna_", "", P25_rna$Var1)
P25_rna$type <- "mRNA"

P25_pro <- melt( tpro_mn[, P25_IDs] )
P25_pro$Var1 <- gsub("pro_", "", P25_pro$Var1)
P25_pro$type <- "protein"

P25_line_dt <- rbind(P25_rna, P25_pro)
P25_line_dt$type <- factor( P25_line_dt$type, levels = c("protein", "mRNA"))

# Create the lineplot
pL1 <- ggplot(P25_line_dt, aes( x = Var1, y = value, group = Var2))+
        facet_wrap( ~ type , nrow = 2 , scales = "free_y")+
        xlab("Cell line")+
        ylab("log2 fold-change")+
        geom_line(size = 0.25, alpha = 0.1)+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25),
              axis.text = element_text( size = 5))

pL1
ggsave("output_files/P25_LP.pdf", pL1, width = 6.5, height = 6, units = "cm")


#### Ribosome progulon (P25) - mRNA & protein coexpression ####

# Get all pairwise mRNA correlations
mRNA_cor <- cor( trna_mn[, P25_IDs] , method = "spearman")
mRNA_cor <- as.data.table( melt(mRNA_cor))
mRNA_cor <- mRNA_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
mRNA_cor$type <- "mRNA"

# Get all pairwise protein correlations
prot_cor <- cor( tpro_mn[, P25_IDs] , method = "spearman")
prot_cor <- as.data.table( melt(prot_cor))
prot_cor <- prot_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
prot_cor$type <- "protein"

# Plot the histograms
plot_dt_hist_P25 <- rbind(prot_cor, mRNA_cor)

pH1 <- ggplot( plot_dt_hist_P25 , aes( x = value , fill = type ))+
        scale_fill_manual( values = c("#00A0BB", "#EC008C"))+
        geom_histogram(position = "identity", alpha = 0.7, binwidth = 0.05, center = 0.025)+
        geom_vline( aes( xintercept = median( prot_cor$value )) , linetype = "dotted", size = 0.25, colour = "#EC008C")+
        geom_vline( aes( xintercept = median( mRNA_cor$value )) , linetype = "dotted", size = 0.25, colour = "#00A0BB")+
        annotate(geom = "text", label = round( median( prot_cor$value ), 2), x = 0, y = 300, colour = "#EC008C")+
        annotate(geom = "text", label = round( median( mRNA_cor$value ), 2), x = 0, y = 200, colour = "#00A0BB")+
        scale_x_continuous( limits = c(-1,1))+
        xlab("Gene coexpression\n[RHO]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))
pH1
ggsave("output_files/P25_HG.pdf", pH1, width = 6, height = 4, units = "cm")


#### Ribosome progulon (P25) - mRNA-to-rotein correlation ####

# Get the cross-correlations
crosscor_P25 <- cor( trna_mn[, P25_IDs], tpro_mn[, P25_IDs] , method = "spearman") 
crosscor_P25 <- as.data.table( melt(crosscor_P25))
crosscor_P25 <- crosscor_P25[ as.character(Var1) == as.character(Var2) ]   # Keep only within-gene (self) correlations

# Plot the mock histograms
pHC1 <- ggplot( crosscor_P25 , aes(x = value ))+
          geom_histogram(position = "identity", binwidth = 0.05, center = 0.025, fill = "grey50")+
          geom_vline( aes( xintercept = median( crosscor_P25$value )) , linetype = "dotted", size = 0.25)+
          annotate(geom = "text", label = round( median( crosscor_P25$value ), 2), x = 0, y = 10)+
          scale_x_continuous( limits = c(-1,1))+
          xlab("mRNA - protein correlation\n[RHO]")+
          ylab("Number of gene pairs")+
          theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
                legend.position = "none", strip.background = element_blank(),
                strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))
pHC1
ggsave("output_files/P25_HC.pdf", pHC1, width = 6, height = 4, units = "cm")


#### Actin cytoskeleton progulon (P16) - lineplot ####

# Get the progulon protein IDs
P16_IDs <- prns[ Progulon_ID == "P16" & prot_in_prn == "yes" , unique(SimpleID) ]

# Get and combine the mRNA and protein expression data of these proteins
P16_rna <- melt( trna_mn[, P16_IDs] )
P16_rna$Var1 <- gsub("rna_", "", P16_rna$Var1)
P16_rna$type <- "mRNA"

P16_pro <- melt( tpro_mn[, P16_IDs] )
P16_pro$Var1 <- gsub("pro_", "", P16_pro$Var1)
P16_pro$type <- "protein"

P16_line_dt <- rbind(P16_rna, P16_pro)
P16_line_dt$type <- factor( P16_line_dt$type, levels = c("protein", "mRNA"))

# Create the lineplot
pL2 <- ggplot(P16_line_dt, aes( x = Var1, y = value, group = Var2))+
        facet_wrap( ~ type , nrow = 2 , scales = "free_y")+
        xlab("Cell line")+
        ylab("log2 fold-change")+
        geom_line(size = 0.25, alpha = 0.1)+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25),
              axis.text = element_text( size = 5))

pL2
ggsave("output_files/P16_LP.pdf", pL2, width = 6.5, height = 6, units = "cm")


#### Actin cytoskeleton progulon (P16) - mRNA & protein coexpression ####

# Get all pairwise mRNA correlations
mRNA_cor <- cor( trna_mn[, P16_IDs] , method = "spearman")
mRNA_cor <- as.data.table( melt(mRNA_cor))
mRNA_cor <- mRNA_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
mRNA_cor$type <- "mRNA"

# Get all pairwise protein correlations
prot_cor <- cor( tpro_mn[, P16_IDs] , method = "spearman")
prot_cor <- as.data.table( melt(prot_cor))
prot_cor <- prot_cor[ as.character(Var1) > as.character(Var2) ]   # Remove duplicates & self-correlations
prot_cor$type <- "protein"

# Plot the histograms
plot_dt_hist_P16 <- rbind(prot_cor, mRNA_cor)

pH2 <- ggplot( plot_dt_hist_P16 , aes( x = value , fill = type ))+
        scale_fill_manual( values = c("#00A0BB", "#EC008C"))+
        geom_histogram(position = "identity", alpha = 0.7, binwidth = 0.05, center = 0.025)+
        geom_vline( aes( xintercept = median( prot_cor$value )) , linetype = "dotted", size = 0.25, colour = "#EC008C")+
        geom_vline( aes( xintercept = median( mRNA_cor$value )) , linetype = "dotted", size = 0.25, colour = "#00A0BB")+
        annotate(geom = "text", label = round( median( prot_cor$value ), 2), x = -0.5, y = 60, colour = "#EC008C")+
        annotate(geom = "text", label = round( median( mRNA_cor$value ), 2), x = -0.5, y = 50, colour = "#00A0BB")+
        scale_x_continuous( limits = c(-1,1))+
        xlab("Gene coexpression\n[RHO]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))
pH2
ggsave("output_files/P16_HG.pdf", pH2, width = 6, height = 4, units = "cm")


#### Actin cytoskeleton progulon (P16) - mRNA-to-protein correlation ####

# Get the cross-correlations
crosscor_P16 <- cor( trna_mn[, P16_IDs], tpro_mn[, P16_IDs] , method = "spearman") 
crosscor_P16 <- as.data.table( melt(crosscor_P16))
crosscor_P16 <- crosscor_P16[ as.character(Var1) == as.character(Var2) ]   # Keep only within-gene (self) correlations

# Plot the mock histograms
pHC2 <- ggplot( crosscor_P16 , aes(x = value ))+
        geom_histogram(position = "identity", binwidth = 0.05, center = 0.025, fill = "grey50")+
        geom_vline( aes( xintercept = median( crosscor_P16$value )) , linetype = "dotted", size = 0.25)+
        annotate(geom = "text", label = round( median( crosscor_P16$value ), 2), x = 0, y = 6.5)+
        scale_x_continuous( limits = c(-1,1))+
        xlab("mRNA - protein correlation\n[RHO]")+
        ylab("Number of gene pairs")+
        theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size = 0.25),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 6), axis.text = element_text(size = 5, colour = "black"),
              legend.position = "none", strip.background = element_blank(),
              strip.text = element_text(size = 6), axis.ticks = element_line(size = 0.25))
pHC2
ggsave("output_files/P16_HC.pdf", pHC2, width = 6, height = 4, units = "cm")






