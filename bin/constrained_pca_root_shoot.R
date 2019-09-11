#Here I will be performing a constrained PCA
library(vegan)

#the working directory ie the path to the data directory within MAE
setwd("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/MAE/data/")

#get the counts
filtered_counts <- read.delim("shoots_roots_separate/RST_mRNA_counts_filtered.txt")

#read in the metadata
metadata_design <- read.delim("shoots_roots_separate/design_RST.txt")

#function to fix groups
fix_group_for_graphing <- function(partial_db) {
  partial_db$Group <- as.character(partial_db$Group)
  for ( i in 1:nrow(partial_db) ) {
    to_paste <- ""
    #Tissue
    if (partial_db$Tissue[i] == "Root") {
      to_paste <- paste(to_paste, "Root_", sep="")
    }
    else {
      to_paste <- paste(to_paste, "Shoot_", sep="")
    }
    
    #Bacteria
    if ( partial_db$Bacteria[i] == "hk" ) {
      to_paste <- paste(to_paste, "HeatKilled_", sep="")
    }
    else if ( partial_db$Bacteria[i] == "nb" ) {
      to_paste <- paste(to_paste, "NoBacteria_", sep="")
    }
    else {
      to_paste <- paste(to_paste, "Syncom_", sep="")
    }
    
    #flg22
    if (partial_db$flg22[i] == "flg22minus") {
      to_paste <- paste(to_paste, "flg22Minus", sep="")
    }
    else {
      to_paste <- paste(to_paste, "flg22Plus", sep="")
    }
    partial_db$Group[i] <- to_paste
  }
  partial_db$Group <- as.factor(partial_db$Group)
  return(partial_db)
}

#Get the root data and time 1 data only
metadata_root <- metadata_design[metadata_design$Tissue == "Root" & metadata_design$Time == "T1",]
root_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_root$Library..))]


#This function will get the log CPM counts of the top500 most variable genes
get_cpm_data <- function(counts_mat, metadata) {
  library(edgeR)
  full_dgelist <- DGEList(counts = counts_mat[3:ncol(counts_mat)], lib.size = colSums(counts_mat[3:ncol(counts_mat)]), samples = metadata, genes=counts_mat$Geneid)
  ### filter low count reads with cpm ###
  countsPerMillion <- cpm(full_dgelist)
  #creates a T and F matrix based on CPM
  countCheck <- countsPerMillion > 1
  filter_out_crazy_genes <- countsPerMillion > quantile(countsPerMillion, probs= .99)
  #creates a vector of the rows that meet the criteria
  keep <- which(rowSums(countCheck) >= 5 & rowSums(filter_out_crazy_genes) < 1)
  full_dgelist <- full_dgelist[keep,] #filter out low count data
  
  ### get the data in CPM after adjusting for normalization factor ###
  full_dgelist <- calcNormFactors(full_dgelist)
  counts <- cpm(full_dgelist, log = T, prior.count = 1)
  row.names(counts) <- full_dgelist$genes$genes
  # select the top 500 genes with largest SD
  genes <- counts[ names(sort(apply(counts,1,sd),decreasing = TRUE)[1:500]), 1:ncol((counts))] 
  
  return(genes)
}

top_500sd_genes <- t(get_cpm_data(filtered_counts, metadata_design))

#This is taken from Isais github at https://github.com/isaisg/ohchibi/blob/master/R/oh.cap.R

# first perform a simple cap 
cap_cond_on_experiment <- capscale(top_500sd_genes ~ Tissue + Time + flg22 + Bacteria + Condition(Experiment), metadata_design, sqrt.dist = T)

cap_sum <- summary(cap_cond_on_experiment)
test_model_signif <- anova.cca(cap_cond_on_experiment, permutations = 5000)
#the model is in fact significant

test_term_signif <- anova.cca(cap_cond_on_experiment, permutations = 5000, by="terms")
#In this model, bacteria does not seem to be a significant factor

pval_model <- test_model_signif[1,4]

cop_results_meta <- cbind( metadata_design, cap_sum$sites)
#the contrabution of each variable to the variance of the data
percvar <- round(100* cap_cond_on_experiment$CCA$eig/cap_cond_on_experiment$CCA$tot.chi,2)

#Extract the proportion explained by MDS
eig <- cap_cond_on_experiment$CA$eig
mds_var <- round(100*eig/sum(eig),2) #what is the variance explained by each dimension

percvar <- c(percvar, mds_var)

#extract the variance
chi <- c(cap_cond_on_experiment$tot.chi, cap_cond_on_experiment$CCA$tot.chi, cap_cond_on_experiment$CA$tot.chi)
variability_table <- cbind(chi, chi/chi[1])
colnames(variability_table) <- c("inertia", "proportion")
rownames(variability_table) <- c("total", "constrained", "unconstrained")

vartot <- variability_table[2,2]*100

#I will then use this data for plotting

#### START ####

### THESE FUNCTIONS ARE IN OTHER R SCRIPTS ###
cap_cond_on_experiment <- perform_cap_scale(top_500sd_genes, metadata_design, "Tissue + Time + flg22 + Bacteria + Condition(Experiment)")
initial_graph <- create_figure_from_caps_data(list_ohpco = cap_cond_on_experiment, comp_a = "CAP1", comp_b = "CAP2", col_val = "Tissue", shape_val = "Time")
initial_graph

#saved as an svg at 900 width and 600 height

#Obvious effect of tissue and time
#Lets try to condition on those two things
cap_cond_on_all_but_bact_flg22 <- perform_cap_scale(top_500sd_genes, metadata_design, "flg22 + Bacteria + Condition(Tissue+Time+Experiment)")
graph <- create_figure_from_caps_data(list_ohpco = cap_cond_on_all_but_bact_flg22, comp_a = "CAP1", comp_b = "CAP2", col_val = "Tissue", shape_val = "flg22")
graph

################################
#### RUN ROOT ONLY ANALYSIS ####
################################
top500_root_genes <- t(get_cpm_data(root_counts, metadata_root))
cap_on_root <- perform_cap_scale(top500_root_genes, metadata_root, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Bacteria", shape_val = "flg22")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Bacteria", fill="flg22"),level = .7, type = "norm", linetype=1)

cap_on_root$variability_table
cap_on_root$pval_model

#saved as an svg at 1000 width and 600 height
################################
#### RUN SHOOT ONLY ANALYSIS ###
################################
metadata_shoot <- metadata_design[metadata_design$Tissue != "Root" & metadata_design$Time == "T1",]
shoot_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_shoot$Library..))]

top500_shoot_genes <- t(get_cpm_data(shoot_counts, metadata_shoot))
cap_on_root <- perform_cap_scale(top500_shoot_genes, metadata_shoot, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Bacteria", shape_val = "flg22")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Bacteria", fill="flg22"),level = .7, type = "norm", linetype=1)

###########################
#### 12 hour time point ####
###########################
#root
metadata_root <- metadata_design[metadata_design$Tissue == "Root" & metadata_design$Time == "T12",]
root_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_root$Library..))]

top500_root_genes <- t(get_cpm_data(root_counts, metadata_root))
cap_on_root <- perform_cap_scale(top500_root_genes, metadata_root, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Bacteria", shape_val = "flg22")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Bacteria", fill="flg22"),level = .7, type = "norm", linetype=1)
cap_on_root$perm_anova_terms

#shoot
metadata_shoot <- metadata_design[metadata_design$Tissue != "Root" & metadata_design$Time == "T12",]
shoot_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_shoot$Library..))]

top500_shoot_genes <- t(get_cpm_data(shoot_counts, metadata_shoot))
cap_on_root <- perform_cap_scale(top500_shoot_genes, metadata_shoot, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Bacteria", shape_val = "flg22")
root_graph +stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Bacteria", fill="flg22"),level = .7, type = "norm", linetype=1)

cap_on_root$perm_anova_terms

#####################################
#### RUN WHOLE SEEDLING ANALYSIS ####
#####################################
#get the counts
filtered_counts <- read.delim("IRM_mRNA_counts_filtered.txt")

#read in the metadata
metadata_design <- read.delim("design_filtered_full.txt")

#keep only the data in counts found in the design
filtered_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_design$Sample..))]

#run
top500_whole_genes <- t(get_cpm_data(filtered_counts, metadata_design))
cap_on_root <- perform_cap_scale(top500_whole_genes, metadata_design, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Bacteria", shape_val = "flg22")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Bacteria", fill="flg22"),level = .7, type = "norm", linetype=1)


###########################
#### All Gene Analysis ####
###########################

#I will now repeat with all genes
get_cpm_data_notop <- function(counts_mat, metadata) {
  library(edgeR)
  full_dgelist <- DGEList(counts = counts_mat[3:ncol(counts_mat)], lib.size = colSums(counts_mat[3:ncol(counts_mat)]), samples = metadata, genes=counts_mat$Geneid)
  ### filter low count reads with cpm ###
  countsPerMillion <- cpm(full_dgelist)
  #creates a T and F matrix based on CPM
  countCheck <- countsPerMillion > 1
  filter_out_crazy_genes <- countsPerMillion > quantile(countsPerMillion, probs= .99)
  #creates a vector of the rows that meet the criteria
  keep <- which(rowSums(countCheck) >= 5 & rowSums(filter_out_crazy_genes) < 1)
  full_dgelist <- full_dgelist[keep,] #filter out low count data
  
  ### get the data in CPM after adjusting for normalization factor ###
  full_dgelist <- calcNormFactors(full_dgelist)
  counts <- cpm(full_dgelist, log = T, prior.count = 1)
  row.names(counts) <- full_dgelist$genes$genes
  return(counts)
}

top500_whole_genes <- t(get_cpm_data_notop(filtered_counts, metadata_design))
cap_on_root <- perform_cap_scale(top500_whole_genes, metadata_design, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Group")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Group"),level = .7, type = "norm", linetype=1)


#roots
top500_root_genes <- t(get_cpm_data_notop(root_counts, metadata_root))
cap_on_root <- perform_cap_scale(top500_root_genes, metadata_root, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Group")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Group"),level = .7, type = "norm", linetype=1)
summary(cap_on_root$Map_cap)
#shoots
metadata_shoot <- metadata_design[metadata_design$Tissue != "Root" & metadata_design$Time == "T1",]
shoot_counts <- filtered_counts[,colnames(filtered_counts) %in% c("Geneid", "Length", as.character(metadata_shoot$Library..))]

top500_shoot_genes <- t(get_cpm_data_notop(shoot_counts, metadata_shoot))
cap_on_root <- perform_cap_scale(top500_shoot_genes, metadata_shoot, "Bacteria + flg22 + Condition(Experiment)")
root_graph <- create_figure_from_caps_data(list_ohpco = cap_on_root, comp_a = "CAP1", comp_b = "CAP2", col_val = "Group")
root_graph + stat_ellipse(inherit.aes = F,aes_string(x="CAP1", y="CAP2", color="Group"),level = .7, type = "norm", linetype=1)
