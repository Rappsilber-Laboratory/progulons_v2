# Progulons_MSB_revision
R scripts and input files for the manuscript "Higher-order modular regulation of the human proteome"

_R scripts_ are in main folder, all files that are read into R scripts are in _input_files_ folder, all files that are created by an R script are in the _output_files_ folder

### Some notes
**Half-life-analysis.R** Script that creates Extended Data Figure 3, analysing mRNA and protein half-lives of progulons

**WGCNA_proteomeHD_treeClust.R** Script that clusters the ProteomeHD + treeClust distance matrix hierachical via WGCNA package

**core_modules_clustering.R** Script that creates the seeds by overlapping OPTICS and ClusterONE (external Java application) clustering

**prn_checkup.R** Does some basic filtering, e.g. removing all progulons that have an AUC > 0.99.

**prn_connectivity** Calculates a connectivity score for each progulon at random forest score cut-offs from 0.5 to 1.0 in 0.01 increments, removing one progulon that is not a significantly connected module at any cut-off

**prns_correlation.R** Calculates and plots the Spearman's correlation between non-overlapping proteins of progulons

**prns_overlap.R** Calculates and plots the overlap in numbers and percentage between the progulons

**prns_basic_stats.R** Calculates and plots some basic statistics about each progulon (how many proteins etc)

**prns_examples_tSNE.R** Produces plots showing three example progulons

**prns_functional_enrichment.R** plots numbers of significant enrichments of GO terms (topGO), Reactome pathways and humap2 complexes in progulons
