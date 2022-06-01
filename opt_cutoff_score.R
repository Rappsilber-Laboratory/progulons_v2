#### load the relevant libraries ####
library(data.table); library(readxl); library(ggplot2);
library(egg); library(treeClust)
set.seed(42) # set random state

#### load the progulon data ####
prns <- fread("input_files/Progulons.csv.gz")

#### create a network of progulon proteins using STRING ####
# load proteomeHD dataset
ProHD <- fread("input_files/ProteomeHD_v1_1.csv.gz")

# select required columns
ProHD <- ProHD[, c("Majority_protein_IDs",
                   grep("^Ratio", names(ProHD), value=TRUE)), with=FALSE]

# convert to a data.frame
ProHD <- data.frame(ProHD, row.names="Majority_protein_IDs")

# keep proteins with 30 or more non-null SILAC ratios
count_features <- function(x){sum(!is.na(x))}
feature_count <- apply(ProHD, 1, count_features)
ProHD <- ProHD[feature_count>=30,]

# Use treeClust to learn dissimilarities
tc_dist <- treeClust.dist(ProHD, d.num=2, verbose=TRUE)

# convert treeClust distance object into a long table format
edges <- as.data.table(reshape2::melt(as.matrix(tc_dist)))
edges <- edges[, .(Protein_1=as.character(Var1),
                   Protein_2=as.character(Var2), edge_weight=(1-value))]

# keep only unique protein associations
edges <- edges[Protein_1 > Protein_2]


# identify the number of pairs that constitute 0.5%
N_edges_0.5pc <- floor(edges[,.N] / 100 * 0.5)

# get the corresponding ubsetof pairs
edges <- edges[order(-edge_weight)][1:N_edges_0.5pc]

edges

# create a list of progulon ids
prn_id_list <- unique(prns[, "Progulon_ID"])

# initialize data.table object
connectivity_test <- data.table()

# loop through each progulon
for (prn_id in unique(prns[, Progulon_ID])){
  
  # loop through each cut-off
  for (score in seq(0.5, 1.0, by=0.01)){
    current_prn <- prns[Progulon_ID==prn_id & Mean_RF_score>score, Protein_IDs]
    
    # calculate number of possible combinaitons
    n_possible_combinations <- choose(length(current_prn), 2)
    
    # find number of observed combinations
    n_connections <- edges[Protein_1 %in% current_prn & Protein_2 %in% current_prn, .N]
    
    # find number of random cominations
    n_random_combinations <- n_possible_combinations / 100 * 0.5
    
    # create a confusion matrix for fischer exact test
    conf_matrix <- matrix(c(n_connections, n_random_combinations,
                            n_possible_combinations-n_connections,
                            n_possible_combinations-n_random_combinations),
                          nrow=2,
                          dimnames=list(c("measured", "random"),
                                        c("linked", "not_linked")))
    
    # calculate pValue for each progulon connectivity
    pValue = fisher.test(conf_matrix)$p.value
    
    # create temporary data.table to store values
    tmp_dt <- data.table(
      "ProgulonID"=prn_id,
      "RF_score_cutoff"=score,
      "No_of_connections_observed"=n_connections,
      "Connectivity_p_value"=pValue,
      "neg_logP"= -1*log2(pValue))
    
    # combine temporary data with initialized table
    connectivity_test <- rbind(connectivity_test, tmp_dt)
  }
}


# make a plot of RF_score_cutoffs vs significance
opt_cutoff <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed"))
opt_cutoff_zoom1 <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed")) + scale_y_continuous(limits=c(-4, 250))
opt_cutoff_zoom2 <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed")) + scale_y_continuous(limits=c(-4, 100))
ggsave("output_files/opt_cutoff.pdf", opt_cutoff, width = 24, height = 12, units = "cm")
ggsave("output_files/opt_cutoff_zoom1.pdf", opt_cutoff_zoom1, width = 24, height = 12, units = "cm")
ggsave("output_files/opt_cutoff_zoom2.pdf", opt_cutoff_zoom2, width = 24, height = 12, units = "cm")


# create a list of progulon ids
prn_id_list <- unique(prns[, "Progulon_ID"])

# initialize data.table object
connectivity_test <- data.table()

# loop through each progulon
for (prn_id in unique(prns[, Progulon_ID])){
  
  # loop through each cut-off
  for (score in seq(0.0, 1.0, by=0.01)){
    current_prn <- prns[Progulon_ID==prn_id & Mean_RF_score>score, Protein_IDs]
    
    # calculate number of possible combinaitons
    n_possible_combinations <- choose(length(current_prn), 2)
    
    # find number of observed combinations
    n_connections <- edges[Protein_1 %in% current_prn & Protein_2 %in% current_prn, .N]
    
    # find number of random cominations
    n_random_combinations <- n_possible_combinations / 100 * 0.5
    
    # create a confusion matrix for fischer exact test
    conf_matrix <- matrix(c(n_connections, n_random_combinations,
                            n_possible_combinations-n_connections,
                            n_possible_combinations-n_random_combinations),
                          nrow=2,
                          dimnames=list(c("measured", "random"),
                                        c("linked", "not_linked")))
    
    # calculate pValue for each progulon connectivity
    pValue = fisher.test(conf_matrix)$p.value
    
    # create temporary data.table to store values
    tmp_dt <- data.table(
      "ProgulonID"=prn_id,
      "RF_score_cutoff"=score,
      "No_of_connections_observed"=n_connections,
      "Connectivity_p_value"=pValue,
      "neg_logP"= -1*log2(pValue))
    
    # combine temporary data with initialized table
    connectivity_test <- rbind(connectivity_test, tmp_dt)
  }
}


# make a plot of RF_score_cutoffs vs significance for full range
opt_cutoff_range <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed"))
opt_cutoff_range_zoom1 <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed")) + scale_y_continuous(limits=c(-4, 250))
opt_cutoff_range_zoom2 <- ggplot(data=connectivity_test, aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="dashed")) + scale_y_continuous(limits=c(-4, 100))
ggsave("output_files/opt_cutoff_range.pdf", opt_cutoff_range, width = 24, height = 12, units = "cm")
ggsave("output_files/opt_cutoff_range_zoom1.pdf", opt_cutoff_range_zoom1, width = 24, height = 12, units = "cm")
ggsave("output_files/opt_cutoff_range_zoom2.pdf", opt_cutoff_range_zoom2, width = 24, height = 12, units = "cm")

#### Investigate significance values between 0.5-0.6 cut-offs ####
# create a list of progulon ids
prn_id_list <- unique(prns[, "Progulon_ID"])

# initialize data.table object
connectivity_test <- data.table()

# loop through each progulon
for (prn_id in unique(prns[, Progulon_ID])){
  
  # loop through each cut-off
  for (score in seq(0.5, 0.6, by=0.01)){
    current_prn <- prns[Progulon_ID==prn_id & Mean_RF_score>score, Protein_IDs]
    
    # calculate number of possible combinaitons
    n_possible_combinations <- choose(length(current_prn), 2)
    
    # find number of observed combinations
    n_connections <- edges[Protein_1 %in% current_prn & Protein_2 %in% current_prn, .N]
    
    # find number of random cominations
    n_random_combinations <- n_possible_combinations / 100 * 0.5
    
    # create a confusion matrix for fischer exact test
    conf_matrix <- matrix(c(n_connections, n_random_combinations,
                            n_possible_combinations-n_connections,
                            n_possible_combinations-n_random_combinations),
                          nrow=2,
                          dimnames=list(c("measured", "random"),
                                        c("linked", "not_linked")))
    
    # calculate pValue for each progulon connectivity
    pValue = fisher.test(conf_matrix)$p.value
    
    # create temporary data.table to store values
    tmp_dt <- data.table(
      "ProgulonID"=prn_id,
      "RF_score_cutoff"=score,
      "No_of_connections_observed"=n_connections,
      "Connectivity_p_value"=pValue,
      "neg_logP"= -1*log2(pValue))
    
    # combine temporary data with initialized table
    connectivity_test <- rbind(connectivity_test, tmp_dt)
  }
}

# plots with constant significance between cut-offs 0.5-0.6
prn01 <- ggplot(data=connectivity_test[ProgulonID=="PRN01"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn30 <- ggplot(data=connectivity_test[ProgulonID=="PRN30"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn20 <- ggplot(data=connectivity_test[ProgulonID=="PRN20"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn29 <- ggplot(data=connectivity_test[ProgulonID=="PRN29"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))

# plots with plummeting significance between cut-offs 0.5-0.6
prn41 <- ggplot(data=connectivity_test[ProgulonID=="PRN41"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn39 <- ggplot(data=connectivity_test[ProgulonID=="PRN39"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn37 <- ggplot(data=connectivity_test[ProgulonID=="PRN37"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn36 <- ggplot(data=connectivity_test[ProgulonID=="PRN36"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))

# plots with rising significance between cut-offs 0.5-0.6
prn38 <- ggplot(data=connectivity_test[ProgulonID=="PRN38"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn23 <- ggplot(data=connectivity_test[ProgulonID=="PRN23"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn11 <- ggplot(data=connectivity_test[ProgulonID=="PRN11"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))
prn28 <- ggplot(data=connectivity_test[ProgulonID=="PRN28"], aes(x=RF_score_cutoff, y=neg_logP, color=ProgulonID)) + geom_line() + geom_hline(aes(yintercept=-1*log2(0.05), linetype="signifcance"))

# arrange plots in one figure
still <- ggarrange( prn01, prn30, prn20, prn29, nrow = 2)
down <- ggarrange( prn41, prn39, prn37, prn36, nrow = 2)
up <- ggarrange( prn38, prn23, prn11, prn28, nrow = 2)

ggsave("output_files/constant.pdf", still, width = 48, height = 24, units = "cm")
ggsave("output_files/plumetting.pdf", down, width = 48, height = 24, units = "cm")
ggsave("output_files/rising.pdf", up, width = 48, height = 24, units = "cm")