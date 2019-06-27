# I want to color the DEGS of each bactria by cluster

direct_w_flg22_degs <- "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/nb_bact_data/"

#read in the cluster info 
cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)


df_bact <- list()
for (i in list.files(direct_w_flg22_degs)) {
  df <- read.csv(paste(direct_w_flg22_degs, i, sep = ""))
  df <- df[(df$logFC > log2(1.5) | df$logFC< -(log2(1.5)))& df$FDR<0.01,]
  df$cluster <- sapply(df$genes, function(x){cluster_info$cluster[cluster_info$Gene==x]})
  df$cluster[is.na(df$cluster >= 1)] <- 0
  df$cluster <- as.character(df$cluster)
  if ( nrow(df) == 0 ) {
    next
  }
  df_bact[[i]] <- df
}


library(RColorBrewer)
library(gridExtra)

colors_plot <- c("grey",brewer.pal(n = 6, name = "Set1"))
for ( i in names(df_bact) ) {
  p <- ggplot(df_bact[[i]], aes(x=genes, y=logFC, color=cluster))+
    geom_point()+
    scale_color_manual(values=colors_plot)+
    labs(x="Genes", color="Cluster")+
    ggtitle(paste(strsplit(i, "[.]")[[1]][1], " flg22 comparisons", sep = ""))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = .5))
    ggsave(paste("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/figures/bact_nb_cluster/",strsplit(i, "[.]")[[1]][1],".pdf", sep = ""), plot = p, device = "pdf" )
}

for ( i in names(df_bact) ) {
  df_bact[[i]]$bact <- strsplit(i, "[.]")[[1]][1]
  if ( i == names(df_bact)[1]) {
    next
  }
  df_bact[[names(df_bact)[1]]] <- rbind(df_bact[[names(df_bact)[1]]], df_bact[[i]])
}

ggplot(df_bact[[names(df_bact)[1]]][df_bact[[names(df_bact)[1]]]$cluster != 0,], aes(x=bact, fill=cluster))+
  geom_bar(position = "Dodge")+
  theme_classic()+
  labs(x="Bacteria", y="Number of DEGs", fill="Cluster")+
  ggtitle("Number of DEGs by Cluster in NB-Bacteria Comparison")+
  scale_color_manual(values=brewer.pal(n = 7, name = "Set1"))+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))
  #ylim(c(0,1000))
