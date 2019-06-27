#### This script is to analyze the clusters with a PCA analysis ####

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)

#read in the cluster data
cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)
table(cluster_info$cluster)


get_zscore <- function(count_df) {
  for ( i in 1:ncol(count_df)) {
    count_df[i,] <- (count_df[i,]-mean(count_df[i,]))/sd(count_df[i,])
  }
  return(count_df)
}

get_cpm_data_clust <- function() {
  library(edgeR)
  full_dgelist <- DGEList(counts = rna_data[,3:ncol(rna_data)], lib.size = colSums(rna_data[,3:ncol(rna_data)]), samples = rna_metadata, genes=rna_data$Gene)
  #filter low count reads with cpm
  countsPerMillion <- cpm(full_dgelist)
  #creates a T and F matrix based on CPM
  countCheck <- countsPerMillion > 1
  filter_out_crazy_genes <- countsPerMillion > quantile(countsPerMillion, probs= .99)
  #creates a vector of the rows that meet the criteria
  keep <- which(rowSums(countCheck) >= 5 & rowSums(filter_out_crazy_genes) < 1)
  full_dgelist <- full_dgelist[keep,] #filter out low count data
  
  full_dgelist <- calcNormFactors(full_dgelist)
  counts <- cpm(full_dgelist, log = T, prior.count = 1)
  
  #counts <- rpkm(full_dgelist, log = T, prior.count = )
  
  #counts <- get_zscore(counts)
  
  data <- matrix(apply(counts[,which(colnames(counts) %in% rna_metadata$Sample..[rna_metadata$Group == unique(rna_metadata$Group)[1]])],1, median), nrow = nrow(counts))
  for ( i in unique(rna_metadata$Group)[2:length(unique(rna_metadata$Group))] ) {
    medians <- apply(counts[,which(colnames(counts) %in% rna_metadata$Sample..[rna_metadata$Group == i])],1, median)
    data <- cbind(data, medians)
  }
  colnames(data) <- unique(rna_metadata$Group)
  rownames(data) <- full_dgelist$genes$genes
  data <- t(data)
  data_list <- list()
  for ( i in c(1,2,3,4,6) ) {
    mm <- data[,colnames(data) %in% cluster_info$Gene[cluster_info$cluster == i]]
    data_list[[i]] <- mm
  }
  return(data_list)
}

#I then want to perform a PCA comparing the LFC of all these genes compared to nb-
library("FactoMineR")
library("factoextra")

data_clusts <- get_cpm_data_clust()
pca_list <- list()
for ( i in c(1,2,3,4,6) ) {
  data_clusts[[i]] <- data_clusts[[i]][which(!(row.names(data_clusts[[i]]) %in% c("MF181.flg22plus","MF181.flg22minus"))),]
  rownames(data_clusts[[i]]) <- sapply(rownames(data_clusts[[i]]), function(x){strsplit(x,".flg")[[1]][1]})
  pca_res <- PCA(data_clusts[[i]], scale.unit = F, ncp = 8, graph = F)
  pca_list[[i]] <- pca_res
}

c_index <- rep("N", nrow(data_clusts[[1]]))
c_index[seq(1,nrow(data_clusts[[1]]), by=2)] <- "M"
c_index[seq(2,nrow(data_clusts[[1]]), by=2)] <- "P"

#res.pca.1 <- PCA(data_clusts[[2]], scale.unit = F, ncp = 5, graph = F)

fviz_pca_ind(pca_list[[2]], axes = c(1,2), col.ind =c_index, repel = T, invisible="quali", title="")+
  theme_minimal()

top_ax <- pca_list[[2]]$var$contrib[(pca_list[[2]]$var$contrib[,2] + pca_list[[2]]$var$contrib[,1])  >= (200/nrow(pca_list[[2]]$var$contrib)),c(1,2)]
write.csv(top_ax, file="/Users/nicholascolaianni/Desktop/clust_2_variables.csv")

suppressors <- c("MF105", "MF374", "MF345", "MF314", "MF79", "MF138")
inducers <- c("MF41", "cl18", "MF161", "MF370", "MF50", "cl69")
suppressors_flg <- c("MF105.1", "MF374.1", "MF345.1", "MF314.1", "MF79.1", "MF138.1")
inducers_flg <- c("MF41.1", "cl18.1", "MF161.1", "MF370.1", "MF50.1", "cl69.1")

pdf("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/cluster_figs/cluster_elipses_combined.pdf",
    width=6, height=3)
#produce figures
for ( i in c(1,2,3,4,6)) {
#  minimized <- as.data.frame(pca_list[[i]]$ind$coord[row.names(pca_list[[i]]$ind$coord) %in% c(suppressors, suppressors_flg,inducers, inducers_flg),])
#  minimized$Names <- row.names(minimized)
  minimized <-as.data.frame(pca_list[[i]]$ind$coord[!(row.names(pca_list[[i]]$ind$coord) %in% c("MF181", "MF181.1")),])
  minimized$Names <- row.names(minimized)
#  trial$col[trial$Names %in% suppressors] <- "No Flg Suppressors"
#  trial$col[trial$Names %in% inducers] <- "No Flg Inducers"#  trial$col[trial$Names %in% suppressors_flg] <- "Flg Suppressors"
#  trial$col[trial$Names %in% inducers_flg] <- "Flg Inducers"
  minimized$col[minimized$Names %in% c(suppressors, suppressors_flg)] <- "Suppressors"
  minimized$col[minimized$Names %in% c(inducers, inducers_flg)] <- "Inducers"
  minimized$col[is.na(minimized$col)] <- "Neutral"
  minimized$shape[minimized$Names %in% c(suppressors_flg, inducers_flg)] <- "Flg"
  minimized$shape[minimized$Names %in% c(suppressors, inducers)] <- "No Flg"
  minimized$shape[is.na(minimized$shape) & grepl(".1", minimized$Names)] <- "Flg"
  minimized$shape[is.na(minimized$shape)] <- "No Flg"
  minimized$col <- factor(minimized$col, levels = c("Suppressors", "Inducers", "Neutral"))
  minimized <- minimized[order(-as.numeric(minimized$col)),]
  #minimized$col <- as.factor(minimized$col)
  
  eig.val <- get_eigenvalue(pca_list[[i]])
  
  plot <- ggplot(minimized, aes(x=Dim.1, y=Dim.2, color=col, shape=shape), na.rm=F)+
    geom_point(na.rm = T)+
    #geom_text_repel()
    stat_ellipse(inherit.aes = F, data=minimized[minimized$col != "Neutral",],aes(x=Dim.1, y=Dim.2, color=col), level = .7, type = "t", linetype=1, na.rm=T)+
    ggtitle(paste("Cluster ", i, sep = ""))+
    labs(color="", shape="")+
    theme_minimal() +
    theme(plot.title = element_text(hjust = .5))+
    scale_color_manual(values=c(Suppressors = "firebrick", Inducers="dodgerblue", Neutral = "grey90")) +
    ylab(paste("PC2 (", round(eig.val[2,2],2), "%)", sep=""))+
    xlab(paste("PC1 (", round(eig.val[1,2],2), "%)", sep=""))
  print(plot)
}
dev.off()
# I want to combine all of the clusters into one pca and then ask which cluster
# explains the variation the best

#this is also how I found out how many dimensions that I needed to keep

eigen_mat <- matrix(nrow = nrow(pca_list[[1]]$ind$coord), ncol = 0)
for ( i in c(1,2,3,4,6) ) {
  print(i)
  eigen <- get_eigenvalue(pca_list[[i]])
  dims_keep <- sum(eigen[,1]>1)
  print(dims_keep)
  print(eigen[dims_keep,3])
  cname <- c(colnames(eigen_mat), rep(paste(i, "_", 1:dims_keep, sep = "")))
  eigen_mat <- cbind(eigen_mat, pca_list[[i]]$ind$coord[,1:dims_keep])
  colnames(eigen_mat) <- cname
}

combined.pca.res <- PCA(eigen_mat, scale.unit = F, ncp = 5, graph = F)

fviz_pca_ind(combined.pca.res, axes = c(1,2), col.ind = c_index, repel = T, invisible="quali", palette = c("Dodgerblue", "Firebrick"), title="")+
  theme_minimal()

fviz_contrib(combined.pca.res, choice = "var", axes =1)
fviz_contrib(combined.pca.res, choice = "var", axes =2)


##################
library(devtools)
library(superheat)
library(ggplot2)
superheat( data,
          row.dendrogram = T, 
          col.dendrogram = T,
          grid.hline.col = "white",
          grid.vline.col = "white",
          #heat.lim = c(min(data_clusts[[1]]), -min(data_clusts[[1]])),
          left.label = "none",
          bottom.label.col = "white",
          title = paste("DEGs from NB Cluster ", i, sep = ""),
          title.size = 8,
          extreme.values.na = FALSE,
          title.alignment = "center",
          bottom.label.text.alignment = "right",
          heat.pal = c("Dodgerblue", "Light Yellow", "Firebrick"),
          bottom.label.text.angle = 90)
library(corrplot)
corrplot(cor(t(data_clusts[[6]])))
cor(t(data_clusts[[6]]))

