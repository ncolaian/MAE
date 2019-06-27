# I want to make a heatmap of the clusters Paulo defined in the pilot experiment
# I need the normalized data so that I can create heatmaps
library(edgeR)

rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)


separate_data_meta <- function( name, data ) {
  samples_to_keep <- rna_metadata[rna_metadata$Bacteria==name,]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
  return(list(mm_data, samples_to_keep))
}

get_zscore <- function(norm_rpkm) {
  return((norm_rpkm-mean(norm_rpkm))/sd(norm_rpkm))
}

main_edgeR_DGE_flg_norm <- function() {
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
  
  normalized_counts <- list()
  narmalized_medians_rpkm <- list()
  for ( i in unique(rna_metadata$Bacteria) ) {
    print(i)
    mm_dataframes <- separate_data_meta(i, mm_data) #getting each samples data frame
    #create a DGE object containing the counts and metadata
    dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
    #get the norm factors from the full dataframe
    dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
  
    count_rpkm <- rpkm(dgList, gene.length = mm_data$Length, log = T)
    
    count_rpkm <- removeBatchEffect(count_rpkm, factor(mm_dataframes[[2]]$Experiment))

    #I then want to get the median for each gene for flg22+ and flg22-
    trial_flg22_minus <- sapply(1:nrow(count_rpkm), function(x) {median(count_rpkm[x,which(colnames(count_rpkm) %in% mm_dataframes[[2]]$Sample[mm_dataframes[[2]]$Treatment == "flg22minus"])]) } )
    trial_flg22_plus <- sapply(1:nrow(count_rpkm), function(x) {median(count_rpkm[x,which(colnames(count_rpkm) %in% mm_dataframes[[2]]$Sample[mm_dataframes[[2]]$Treatment == "flg22plus"])]) } )
    medians <- matrix(c(get_zscore( trial_flg22_plus ), get_zscore(trial_flg22_minus)), ncol = 2, byrow = F)
    normalized_counts[[i]] <- count_rpkm
    narmalized_medians_rpkm[[i]] <- medians
  }
  

  return(list(mm_data$Gene,normalized_counts, narmalized_medians_rpkm))
}

normalized_counts <- main_edgeR_DGE_flg_norm()

norm_genes <- normalized_counts[[1]]
normalized_counts_all <- normalized_counts[[2]]
medians_of_norm <- normalized_counts[[3]]
#read clusters of genes
cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)
table(cluster_info$cluster)

#Read in the nb DEG for flg22+
top_nbflg22 <- read.csv("/Users/nicholascolaianni/Desktop/nb_DEG.csv", stringsAsFactors = F)
top_nbflg22 <- top_nbflg22[top_nbflg22$logFC > log2(1.5) | top_nbflg22$logFC < -(log2(1.5)),]


#filter the cluster info
cluster_info <- cluster_info[cluster_info$Gene %in% top_nbflg22$genes,]

cluster_data <- list()
for ( j in 1:6 ) {
  data_clust <- c()
  clnms <- c()
  for ( i in names(medians_of_norm) ) {
    #put the - first and then the plus
    rows_2_use <- which(norm_genes %in% cluster_info$Gene[cluster_info$cluster == j])
    #this is for all genes found DE in NB
    #rows_2_use <- which(norm_genes %in% top_nbflg22$genes)
    
    data_clust <- append(data_clust,  c(medians_of_norm[[i]][rows_2_use,2],medians_of_norm[[i]][rows_2_use,1] ))
    clnms <- append(clnms, c(paste(i, "minus", sep = "-"), paste(i, "plus", sep = "-")))
  }
  mm_mat <- matrix(data_clust, ncol = length(clnms), byrow = F)
  row.names(mm_mat) <- norm_genes[which(norm_genes %in% cluster_info$Gene[cluster_info$cluster == i])]
  colnames(mm_mat) <- clnms
  cluster_data[[j]] <- mm_mat
}



### PLOTTING
library(devtools)
library(superheat)
library(ggplot2)

for ( i in c(1,2,3,4,6)) {
  png(paste("/Users/nicholascolaianni/Desktop/NB_flg22_DEGs_all", i, ".png"), height = 800, width = 1500)
  superheat(cluster_data[[3]],
          row.dendrogram = T, 
          col.dendrogram = T,
          heat.pal = c("Dodgerblue", "Light Yellow", "Firebrick"),
          grid.hline.col = "white",
          grid.vline.col = "white",
          heat.lim = c(-2, 2),
          left.label = "none",
          bottom.label.col = "white",
          title = paste("DEGs from NB Cluster ", i, sep = ""),
          title.size = 8,
          extreme.values.na = FALSE,
          title.alignment = "center",
          bottom.label.text.alignment = "right",
          bottom.label.text.angle = 90)
  dev.off()
  
}

#interesting plot idea
trial <- log(medians_of_norm$nb[,1]/medians_of_norm$nb[,2])

plot( trial[abs(trial) > 1.5], type = "p")
plot(top_nbflg22$logFC)
awhich ( medians_of_norm$MF2[,1]/medians_of_norm$MF2[,2] > 1000 )


norm_genes[which ( medians_of_norm$MF376[,1]/medians_of_norm$MF376[,2] > 1000 )]
