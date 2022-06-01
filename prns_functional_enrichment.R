######################################################################################################## #
# Project: Progulons, with Groth lab     
# Purpose: Test functional enrichment (GO terms, Reactome pathways and hu.MAP complexes) among progulons
# Author: E Rullmann and G Kustatscher
# Date: February 2022 
######################################################################################################## #

# Read in the necessary libraries
library(plyr); library(data.table); library(topGO); library(ggplot2); library(viridis); library(stringr); library(egg)

#### Prep the data ####
setwd("")

# Read in ProteomeHD to determine the list of proteins that need to be tested
ProHD <- fread("input_files/ProteomeHD_v1_1.csv.gz")
ProHD_ratios <- ProHD[, .SD, .SDcols = colnames(ProHD) %like% "Ratio"]      # Limit to data columns
feature_count <- apply( ProHD_ratios, 1, function(x){ sum( !is.na(x)) } )   # Calculate number of features per protein
ProHD <- ProHD[ feature_count >= 30 ,]                                      # Keep only proteins that were actually included in the Progulon analysis
ProHD_proteins <- unique( ProHD[, SimpleID_1 := gsub(";.+", "", Majority_protein_IDs)][, gsub("-.+", "", SimpleID_1) ] )  # Unique genes (protein isoforms removed)

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

prns <- prns[prot_in_prn=="yes"]

prns[, Protein_IDs := gsub(";.+", "", Protein_IDs)][, Protein_IDs := gsub("-.+", "", Protein_IDs)]           # Simplify protein IDs 
prns <- dlply( prns, "Progulon_ID", function(x){  unique( x[ x$prot_in_prn == "yes" , "Protein_IDs" ] ) })   # List of proteins in each progulon

# Download the human GO associations file (goa_human.gaf) from http://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/
# Then read them in and prep as follows:
all_GO <- fread("input_files/goa_human.gaf.gz")                             # Loading human protein GO associations
all_GO <- all_GO[, .(ID = V2, Qualifier = V4, GO_ID = V5, Aspect = V9)]     # Keep relevant columns only
all_GO <- all_GO[ !( Qualifier %like% "acts_upstream_of" | 
                     Qualifier %like% "NOT" ) ]                             # Keep only associations with reasonable qualifiers
all_GO <- all_GO[ ID %in% ProHD_proteins ]                                  # Keep only proteins that were part of our analysis

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "all_GO")] )  


#### GO enrichment analysis per progulon - Molecular function ####

# Prepare the data for topGO
GO <- all_GO[ Aspect == "F" ]                                      # Extract the relevant GO aspect
GO <- unique(GO)                                                   # Discard duplicate annotations
geneID2GO <- dlply( GO, "ID", function(x){as.character(x$GO_ID)})  # GO annotation for all proteins in the analysis, for which GO annotation is available

# Turn progulons into allGenes input for topGO
# allGenes_list is a list of progulons, with each progulon being a named factor that can be used for topGO's allGenes parameter.
# The names are all proteins in ProHD for which progulon-membership has been tested AND for which we have GO annotations (the "universe"). 
# The factor levels indicate whether proteins are in a progulon or not
geneNames <- names(geneID2GO)
allGenes_list <- lapply(prns, function(x){  geneList <- factor(as.integer(geneNames %in% x))
                                            names(geneList) <- geneNames
                                            return(geneList) })

# Create pilot topGOdata object (which will be updated with new allGenes in each iteration below)
GOdata <- new("topGOdata", ontology = "MF", allGenes = allGenes_list[[1]], annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

# Screen every progulon for GO term enrichment
pvalues <- data.frame(GO_term = usedGO(GOdata))
for(i in 1:length(allGenes_list)){
  GOdata <- updateGenes(GOdata, allGenes_list[[i]])                           # Update the pilot GOdata object with the real progulon to be tested
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")     # Calculate enrichment considering graph structure
  current_pvalue <- as.data.frame( score(result))                             # Get the pvalues for enrichment of all GO terms for this gene subset
  colnames(current_pvalue) <- names(allGenes_list)[i]                         # Attach the progulon ID
  pvalues <- merge(pvalues, current_pvalue, by.x="GO_term", by.y="row.names")
  print(i)
}

pvalues_MF <- pvalues

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "all_GO", "pvalues_MF")] )  


#### GO enrichment analysis per progulon - Cellular component ####

# Prepare the data for topGO
GO <- all_GO[ Aspect == "C" ]                                      # Extract the relevant GO aspect
GO <- unique(GO)                                                   # Discard duplicate annotations
geneID2GO <- dlply( GO, "ID", function(x){as.character(x$GO_ID)})  # GO annotation for all proteins in the analysis, for which GO annotation is available

# Turn progulons into allGenes input for topGO
# allGenes_list is a list of progulons, with each progulon being a named factor that can be used for topGO's allGenes parameter.
# The names are all proteins in ProHD for which progulon-membership has been tested AND for which we have GO annotations (the "universe"). 
# The factor levels indicate whether proteins are in a progulon or not
geneNames <- names(geneID2GO)
allGenes_list <- lapply(prns, function(x){  geneList <- factor(as.integer(geneNames %in% x))
names(geneList) <- geneNames
return(geneList) })

# Create pilot topGOdata object (which will be updated with new allGenes in each iteration below)
GOdata <- new("topGOdata", ontology = "CC", allGenes = allGenes_list[[1]], annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

# Screen every progulon for GO term enrichment
pvalues <- data.frame(GO_term = usedGO(GOdata))
for(i in 1:length(allGenes_list)){
  GOdata <- updateGenes(GOdata, allGenes_list[[i]])                           # Update the pilot GOdata object with the real progulon to be tested
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")     # Calculate enrichment considering graph structure
  current_pvalue <- as.data.frame( score(result))                             # Get the pvalues for enrichment of all GO terms for this gene subset
  colnames(current_pvalue) <- names(allGenes_list)[i]                         # Attach the progulon ID
  pvalues <- merge(pvalues, current_pvalue, by.x="GO_term", by.y="row.names")
  print(i)
}

pvalues_CC <- pvalues

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "all_GO", "pvalues_MF", "pvalues_CC")] )  


#### GO enrichment analysis per progulon - Biological process ####

# Prepare the data for topGO
GO <- all_GO[ Aspect == "P" ]                                      # Extract the relevant GO aspect
GO <- unique(GO)                                                   # Discard duplicate annotations
geneID2GO <- dlply( GO, "ID", function(x){as.character(x$GO_ID)})  # GO annotation for all proteins in the analysis, for which GO annotation is available

# Turn progulons into allGenes input for topGO
# allGenes_list is a list of progulons, with each progulon being a named factor that can be used for topGO's allGenes parameter.
# The names are all proteins in ProHD for which progulon-membership has been tested AND for which we have GO annotations (the "universe"). 
# The factor levels indicate whether proteins are in a progulon or not
geneNames <- names(geneID2GO)
allGenes_list <- lapply(prns, function(x){  geneList <- factor(as.integer(geneNames %in% x))
names(geneList) <- geneNames
return(geneList) })

# Create pilot topGOdata object (which will be updated with new allGenes in each iteration below)
GOdata <- new("topGOdata", ontology = "BP", allGenes = allGenes_list[[1]], annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

# Screen every progulon for GO term enrichment
pvalues <- data.frame(GO_term = usedGO(GOdata))
for(i in 1:length(allGenes_list)){
  GOdata <- updateGenes(GOdata, allGenes_list[[i]])                           # Update the pilot GOdata object with the real progulon to be tested
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")     # Calculate enrichment considering graph structure
  current_pvalue <- as.data.frame( score(result))                             # Get the pvalues for enrichment of all GO terms for this gene subset
  colnames(current_pvalue) <- names(allGenes_list)[i]                         # Attach the progulon ID
  pvalues <- merge(pvalues, current_pvalue, by.x="GO_term", by.y="row.names")
  print(i)
}

pvalues_BP <- pvalues

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "all_GO", "pvalues_MF", "pvalues_CC", "pvalues_BP")] )  


#### Plot GO enrichment analysis ####

# Set cut-off
p_value_cutoff <- 1e-06

# Convert to data.tables
BP <- as.data.table(pvalues_BP)
MF <- as.data.table(pvalues_MF)
CC <- as.data.table(pvalues_CC)

# Number of GO terms that are significantly enriched, per progulon
BP <- melt(BP, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
BP <- BP[, .(BP_terms = sum(pvalue < p_value_cutoff)), by = prn ]
MF <- melt(MF, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
MF <- MF[, .(MF_terms = sum(pvalue < p_value_cutoff)), by = prn ]
CC <- melt(CC, id.vars = "GO_term", variable.name = "prn", value.name = "pvalue")
CC <- CC[, .(CC_terms = sum(pvalue < p_value_cutoff)), by = prn ]

# Combine results into one table
GO <- merge(BP, MF)
GO <- merge(GO, CC)

## Re-order progulons as in the correlation figure panel

# Load the correlation data
all_prn_combinations <- fread("output_files/ProgulonCor.csv.gz")

# Expand data to be able to get a complete matrix, i.e. append duplicates
all_prn_combinations <- rbind(all_prn_combinations[, .(PRN_A, PRN_B, RHO)],
                              all_prn_combinations[, .(PRN_B = PRN_A, PRN_A = PRN_B, RHO)])

# Load the progulaon annotations
prn_annot <- fread("input_files/Progulon_annotation.csv")
setnames(prn_annot, old = c("Progulon_ID", "Progulon_name"), new = c("Progulon", "Function"))

# Append functional annotation
all_prn_combinations[, Function_1 := prn_annot[ match( all_prn_combinations[, PRN_A], prn_annot[, Progulon] ), Function ]]
all_prn_combinations[, Function_2 := prn_annot[ match( all_prn_combinations[, PRN_B], prn_annot[, Progulon] ), Function ]]

# Cast into a correlation matrix
cor_mat <- dcast( all_prn_combinations, Function_1 ~ Function_2, value.var = "RHO" )
my_rownames <- cor_mat[, Function_1] 
cor_mat[, Function_1 := NULL ]
cor_mat <- as.data.frame( cor_mat )
rownames(cor_mat) <- my_rownames

# Group progulons by correlation
my_dist <- as.dist( (1-cor_mat)/2 )
my_dist[ is.na(my_dist) ] <- 0        # For the few combinations that had no non-overlapping proteins, set the resulting NA to zero
my_clust <- hclust(my_dist, method = "average")
new_prn_order <- rownames(cor_mat)[ my_clust$order ]

# Add annotated names and reorder progulons 
GO <- merge(GO, prn_annot[, .(prn = Progulon, Function)])
GO <- GO[, .(Function, BP_terms, MF_terms, CC_terms)]
GO <- melt(GO, id.vars = "Function")
GO[, Function := factor(Function, levels = new_prn_order)]
GO[, variable := factor(variable, levels = c("MF_terms", "CC_terms", "BP_terms"))]

# Create the GO plots
pGO <- ggplot(GO, aes(x = Function, y = value, fill = variable))+
       geom_bar(stat="identity")+
       coord_flip()+
       scale_fill_viridis(discrete = TRUE, option = "D", direction = -1)+
       theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size=0.25), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(),
             axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6))

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "GO", "pGO")] )  


#### Reactome pathway enrichment analysis ####

# Read in Reactome data (available at https://reactome.org/download/current/UniProt2Reactome.txt)
R1 <- fread("input_files/UniProt2Reactome_2022.txt", header = FALSE)      # Read in lowest level Reactome pathway annotation
R1 <- unique( R1[ V6 == "Homo sapiens" , .(ID = V1, Pathway = V4) ])    # Restrict to human proteins and keep only relevant columns
R1 <- R1[ !grepl("-", R1$ID, fixed = TRUE) ]                            # Remove protein isoforms [is that really necessary?]

# Read progulon data
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

prns <- prns[prot_in_prn=="yes"]
progulons <- prns
# filter both datasets
progulons <- progulons[progulons$Protein_IDs %in% R1$ID] # reduce progulons to only IDs that are in R1
R1 <- R1[R1$ID %in% progulons$Protein_IDs]               # reduce pathways to only IDs that are in progulons

prog_uni <- unique(progulons$Protein_IDs[progulons$Protein_IDs %in% R1$ID])
path_uni <- unique(R1$ID[R1$ID %in% progulons$Protein_IDs])

selected_pathways <- R1[, .N, Pathway][ N > 20 ][, Pathway ] # minimal size of pathway
R1 <- R1[ Pathway %in% selected_pathways ]

selected_progulons <- progulons[, .N, Progulon_ID][ N > 10 ][, Progulon_ID] # minimal size of progulon
progulons <- progulons[Progulon_ID %in% selected_progulons]

prog_table <- as.data.frame(table(progulons$Progulon_ID))
hist(prog_table$Freq, breaks = c(seq(0,850,10))) # create a histogram to see the distribution of prog sizes

R1_table <- as.data.frame(table(R1$Pathway))
hist(R1_table$Freq, breaks = c(seq(0,500,10))) # create a histogram to see the distribution of pathway sizes

# Loop for contingency tables for comparing all progulons with all pathways
contig_table <- data.frame(matrix(nrow = 2,ncol = 2),row.names = c("In Progulon","Not in Progulon")) # create an empty contingency table
colnames(contig_table) <- c("In Pathway","Not in Pathway")

fishers_table <- data.frame(Pathway=character(),Progulon=character(),P_value=integer()) # create empty dataframe to append p-values
colnames(fishers_table) <- c("Pathway","Progulon","P-value")

# loop every pathway over every progulon 
for (i in selected_pathways) {
  for (x in selected_progulons) {
    contig_table[1,1] <- sum(progulons$Protein_IDs[progulons$Progulon_ID==x]%in%R1$ID[R1$Pathway==i],na.rm = TRUE) # IDs both in pathway and progulon
    contig_table[1,2] <- prog_table$Freq[prog_table$Var1==x]-contig_table[1,1] # IDs in Progulon but not in pathway
    contig_table[2,1] <- R1_table$Freq[R1_table$Var1==i]-contig_table[1,1] # IDs in Pathway but now in progulon
    contig_table[2,2] <- nrow(R1)+nrow(progulons)-contig_table[1,1]-contig_table[1,2]-contig_table[2,1]
    fish <- fisher.test(contig_table)
    fishers_table <- rbind(fishers_table,data.frame(Pathway=i,Progulon=x,P_value=fish$p.value))
  }
}

adjusted_P <- p.adjust(fishers_table$P_value, method = "bonferroni")
fishers_table$Bonferroni_P_value <- adjusted_P

bonferroni_table <- data.frame(Pathway=character(),Progulon=character(),P_value=integer(),Bonferroni_P_value=integer(),Bonferroni_P_value_PRNwise=integer())

for (y in selected_progulons) {
  PRN_fish <- fishers_table[fishers_table$Progulon==y,]
  adjusted_P <- p.adjust(PRN_fish$P_value, method = "bonferroni")
  PRN_fish$Bonferroni_P_value_PRNwise <- adjusted_P
  bonferroni_table <- rbind(bonferroni_table,PRN_fish)
}

bonferroni_table <- na.omit(bonferroni_table)

PRN_translate <- fread("input_files/Progulon_annotation.csv")
setnames(PRN_translate, old = c("Progulon_ID", "Progulon_name"), new = c("Progulon", "Function"))


p_distribution <- as.data.frame(table(bonferroni_table$Progulon[bonferroni_table$Bonferroni_P_value_PRNwise <= 1*10^-6])) # filter only significant overlaps

p_distribution <- merge(p_distribution,PRN_translate[,Progulon,Function], by.x = "Var1", by.y = "Progulon", all.y = TRUE)

p_distribution$Function <- factor(p_distribution$Function, levels = levels(GO$Function))

# Create the Reactome plot
pRE <- ggplot(p_distribution, aes(x = Function, y = Freq))+
       geom_bar(stat="identity")+ 
       coord_flip()+
       theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size=0.25), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6), legend.position = "none")

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "GO", "pGO", "pRE")] )  


#### huMAP protein complexes enrichment analysis ####

# Read in human complex data from huMAP2, available here:
# http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt
R1 <- fread("input_files/humap2_complexes_20200809.txt")                     

col_split <- str_split_fixed(R1$Uniprot_ACCs," ", n = Inf)
R1 <- cbind(R1[,HuMAP2_ID],col_split)
R1[R1==""] <- NA
R1 <- melt.data.table(as.data.table(R1),id.vars = "V1", na.rm = TRUE)
R1$variable <- NULL
colnames(R1) <- c("complex","ID")

# Read progulon data
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

prns <- prns[prot_in_prn=="yes"]
progulons <- prns

progulons$Protein_IDs <- sub("\\;.*$","",progulons$Protein_IDs)
progulons <- progulons[ !grepl("-", progulons$Protein_IDs, fixed = TRUE) ] # remove protein isoforms [is that really necessary?]

# filter both datasets
progulons <- progulons[progulons$Protein_IDs %in% R1$ID] # reduce progulons to only IDs that are in R1
R1 <- R1[R1$ID %in% progulons$Protein_IDs] # reduce pathways to only IDs that are in progulons

prog_uni <- unique(progulons$Protein_IDs[progulons$Protein_IDs %in% R1$ID])
path_uni <- unique(R1$ID[R1$ID %in% progulons$Protein_IDs])

selected_complexs <- R1[, .N, complex][ N > 1 ][, complex ] # minimal size of complex
R1 <- R1[ complex %in% selected_complexs ]

selected_progulons <- progulons[, .N, Progulon_ID][ N > 1 ][, Progulon_ID] # minimal size of progulon
progulons <- progulons[Progulon_ID %in% selected_progulons]

prog_table <- as.data.frame(table(progulons$Progulon_ID))
hist(prog_table$Freq, breaks = c(seq(0,1000,20))) # create a histogram to see the distibution of prog sizes

R1_table <- as.data.frame(table(R1$complex))
hist(R1_table$Freq, breaks = c(seq(0,100,2))) # create a histogram to see the distribution of complex sizes

# loop for contigency tables for comparing all progulons with all complexs
contig_table <- data.frame(matrix(nrow = 2,ncol = 2),row.names = c("In Progulon","Not in Progulon")) # create an empty contigency table
colnames(contig_table) <- c("In complex","Not in complex")

fishers_table <- data.frame(complex=character(),Progulon=character(),P_value=integer()) # create empty dataframe to append p-values
colnames(fishers_table) <- c("complex","Progulon","P-value")

# loop every complex over every progulon 
pb <- txtProgressBar(min = 0, max = length(selected_complexs), style = 3) 
stepi <- 0
for (i in selected_complexs) {
  for (x in selected_progulons) {
    contig_table[1,1] <- sum(progulons$Protein_IDs[progulons$Progulon_ID==x]%in%R1$ID[R1$complex==i],na.rm = TRUE) # IDs both in complex and progulon
    contig_table[1,2] <- prog_table$Freq[prog_table$Var1==x]-contig_table[1,1] # IDs in Progulon but not in complex
    contig_table[2,1] <- R1_table$Freq[R1_table$Var1==i]-contig_table[1,1] # IDs in complex but now in progulon
    contig_table[2,2] <- nrow(R1)+nrow(progulons)-contig_table[1,1]-contig_table[1,2]-contig_table[2,1]
    fish <- fisher.test(contig_table)
    fishers_table <- rbind(fishers_table,data.frame(complex=i,Progulon=x,P_value=fish$p.value))
  }
  stepi <- stepi + 1
  setTxtProgressBar(pb,stepi)
}

adjusted_P <- p.adjust(fishers_table$P_value, method = "bonferroni")
fishers_table$Bonferroni_P_value <- adjusted_P

bonferroni_table <- data.frame(complex=character(),Progulon=character(),P_value=integer(),Bonferroni_P_value=integer(),Bonferroni_P_value_PRNwise=integer())

for (y in selected_progulons) {
  PRN_fish <- fishers_table[fishers_table$Progulon==y,]
  adjusted_P <- p.adjust(PRN_fish$P_value, method = "bonferroni")
  PRN_fish$Bonferroni_P_value_PRNwise <- adjusted_P
  bonferroni_table <- rbind(bonferroni_table,PRN_fish)
}

bonferroni_table <- na.omit(bonferroni_table)

# translate IDs and create figure
PRN_translate <- fread("input_files/Progulon_annotation.csv")
setnames(PRN_translate, old = c("Progulon_ID", "Progulon_name"), new = c("Progulon", "Function"))

p_distribution <- as.data.frame(table(bonferroni_table$Progulon[bonferroni_table$Bonferroni_P_value_PRNwise <= 1*10^-6])) # filter only significant overlaps

p_distribution <- merge(p_distribution,PRN_translate[,Progulon,Function], by.x = "Var1", by.y = "Progulon", all.y = TRUE)

p_distribution$Function <- factor(p_distribution$Function, levels = levels(GO$Function))

# Create the Reactome plot
pHM <- ggplot(p_distribution, aes(x = Function, y = Freq))+
       geom_bar(stat="identity")+ 
       coord_flip()+
       theme(plot.background = element_blank(), panel.background = element_rect(fill=NA, colour="black", size=0.25), panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), axis.title = element_blank(), axis.ticks.y = element_blank(),
       axis.text.y = element_text(size=6), axis.text.x = element_text(size = 6), legend.position = "none")

# Clear up workspace
rm( list = ls()[! ls() %in% c("prns", "GO", "pGO", "pRE", "pHM")] )  


#### Arrange and save the output plots ####

# Combine plots
p <- ggarrange(pGO, pRE, pHM, nrow = 1)

# Save the plot
ggsave("output_files/GO_React_huMAP.pdf", p, width = 26, height = 10.5, units = "cm")
ggsave("output_files/GO_React_huMAP.svg", p, width = 26, height = 10.5, units = "cm")


















