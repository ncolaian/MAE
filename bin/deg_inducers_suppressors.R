#this code will be used to perform the gene expression of the inducers vs the supressors

library(edgeR)

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)


separate_data_meta_noflg <- function( data ) {
  suppressors <- c("MF105", "MF374", "MF345", "MF314", "MF79", "MF138")
  inducers <- c("MF41", "cl18", "MF161", "MF370", "MF50", "cl69")
  samples_to_keep <- rna_metadata[rna_metadata$Bacteria %in% c(suppressors, inducers) & rna_metadata$Treatment == "flg22minus",]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
  samples_to_keep$supressor[samples_to_keep$Bacteria %in% suppressors] <- 1
  samples_to_keep$supressor[samples_to_keep$Bacteria %in% inducers] <- 0
  return(list(mm_data, samples_to_keep))
}

main_edgeR_DGE_flg <- function() {
  #first calculate the normalization facter data together so that samples can be comparable 
  full_dgelist <- DGEList(counts = rna_data[,3:ncol(rna_data)], lib.size = colSums(rna_data[,3:ncol(rna_data)]), samples = rna_metadata, genes=rna_data$Gene)
  #filter low count reads with cpm
  countsPerMillion <- cpm(full_dgelist)
  #creates a T and F matrix based on CPM
  countCheck <- countsPerMillion > 1
  #creates a vector of the rows that meet the criteria
  keep <- which(rowSums(countCheck) >= 5)
  full_dgelist <- full_dgelist[keep,] #filter out low count data
  full_dgelist <- calcNormFactors(full_dgelist)
  
  #only keep the data from the genes that pass the tests
  mm_data <- rna_data[rna_data$Gene %in% full_dgelist$genes$genes,]
  
  #pdf("/Users/nicholascolaianni/Desktop/MAE_mds_plots.pdf")
  
  mm_dataframes <- separate_data_meta_noflg(mm_data)
   #create a DGE object containing the counts and metadata
  dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
  #get the norm factors from the full dataframe
  dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
    
    #model matrix design
    # bact <- factor(dgList$samples$Bacteria, levels= c("nb", i))
  experiment <- factor(dgList$samples$Experiment, levels = c("e1", "e2", "e3"))
  treatment <- factor(dgList$samples$Treatment, levels = c("flg22minus", "flg22plus"))
  bact <- factor(dgList$samples$Bacteria)
  supp <- factor(dgList$samples$supressor, levels = c(0, 1))
  designMat <- model.matrix(~supp + experiment)
    
    #fit dispersion parameter
  dgList <- estimateGLMCommonDisp(dgList, design=designMat) #uses a common estimate across genes
  dgList <- estimateGLMTrendedDisp(dgList, design=designMat) #fits an estimate based on mean variance trends across the data
  dgList <- estimateGLMTagwiseDisp(dgList, design=designMat) #compute a genewise variance
    
    
    #perform dgexpression
  fit <- glmFit(dgList, designMat)
    
    #Perform a logratiotest with pvals and cpm
    #the coef makes sure we are testing for DGE for the addition of flg22
  lrt <- glmLRT(fit, coef="supp1") 
  return(lrt)
}

trial <- main_edgeR_DGE_flg()

#what are the differentially expressed genes?
#Are they from any particular cluster?
w <- topTags(trial, n=nrow(trial), sort.by = "none", p.value = 0.01)
w <- w[abs(w$table$logFC) > log2(1.5),]
write.csv(w, "/Users/nicholascolaianni//Desktop/comparison_degs_trial.csv")

#how many of these genes come from the clusters?
cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)

for ( i in w$table$genes ) {
  if ( i %in% cluster_info$Gene) {
    w$table$cluster[w$table$genes == i] <- cluster_info$cluster[cluster_info$Gene == i]
  }
  else {
    w$table$cluster[w$table$genes == i] <- 0
  }
}

library(ggplot2)

ggplot(w$table, aes(x=factor(cluster, levels = c(0,1,2,3,4,6))))+
  geom_bar(fill="grey12")+
  theme_minimal()+
  #ggtitle("Cluster Residence of DEGs Between Inducers and Suppressors")+
  xlab("Cluster")+
  ylab("Number of DEGs")+
  theme(plot.title = element_text(hjust = .5))

#go enrichment with this set to show its mainly defense related stuff
library(clusterProfiler)
library(AnnotationHub)

hub <- AnnotationHub()
q <- query(hub, "Arabidopsis thaliana")
id <- q$ah_id[length(q)]
arab <- hub[[id]]

ego1 <- enrichGO(gene = w$table$genes,
                 OrgDb = arab,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 keyType = 'TAIR',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05
)

dotplot(ego1, showCategory = 20)

ego2 <- enrichGO(gene = w$table$genes[w$table$cluster==0],
                 OrgDb = arab,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 keyType = 'TAIR',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05
)

dotplot(ego2, showCategory = 20)
