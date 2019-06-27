#this script will look for the depletion of or enrichment of DEGS in each cluster


#diresctories with deg files
bact_comp <- "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/old_model/final_sep_data/bactm_vs_bactp/"
nb_bact <- "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/old_model/final_sep_data/nbm_vs_bactm/"

#read in the cluster info 
cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)

read_directory_files <- function(path) {
  df_bact <- list()
  for (i in list.files(path)) {
    df <- read.csv(paste(path, i, sep = ""))
    df <- df[(df$logFC > log2(1.5) | df$logFC< -(log2(1.5)))& df$FDR<0.01,]
    df$cluster <- sapply(df$genes, function(x){cluster_info$cluster[cluster_info$Gene==x]})
    df$cluster[is.na(df$cluster >= 1)] <- 0
    df$cluster <- as.character(df$cluster)
    if ( nrow(df) == 0 ) {
      next
    }
    df_bact[[i]] <- df
  }
  return(df_bact)
}


bact_list <- read_directory_files(bact_comp)
nb_list <- read_directory_files(nb_bact)

intersection_list <- list()
for (i in union( names(nb_list), names(bact_list) ) ) {
  if ( !(i %in% names(bact_list)) ) {
    intersection_list[[i]] <- as.character(nb_list[[i]]$genes)
    print(i)
  }
  else{
    intersection_list[[i]] <- union( nb_list[[i]]$genes, bact_list[[i]]$genes )
  }
  ud_v <- c()
  for ( j in intersection_list[[i]] ) {
    #handle if they both are diff expressed
    if ( j %in% nb_list[[i]]$genes & j %in% bact_list[[i]]$genes) {
      #take the higher lfc value
      if (nb_list[[i]]$logFC[nb_list[[i]]$genes == j] > bact_list[[i]]$logFC[bact_list[[i]]$genes == j]) {
        if ( nb_list[[i]]$logFC[nb_list[[i]]$genes == j] > 0 ) {
          ud_v <- c(ud_v, 1)
        }
        else {
          ud_v <- c(ud_v, 0)
        }
      }
      else {
        if ( bact_list[[i]]$logFC[bact_list[[i]]$genes == j] > 0 ) {
          ud_v <- c(ud_v, 1)
        }
        else {
          ud_v <- c(ud_v, 0)
        }
      }
    }
    #handle just no bacteria comparison
    else if ( j %in% nb_list[[i]]$genes ) {
      if ( nb_list[[i]]$logFC[nb_list[[i]]$genes == j] > 0 ) {
        ud_v <- c(ud_v, 1)
      }
      else {
        ud_v <- c(ud_v, 0)
      }
    }
    #handle just the bacteria comparison
    else {
      if ( bact_list[[i]]$logFC[bact_list[[i]]$genes == j] > 0 ) {
        ud_v <- c(ud_v, 1)
      }
      else {
        ud_v <- c(ud_v, 0)
      }
    }
  }
  intersection_list[[i]] <- matrix(c(intersection_list[[i]], ud_v), ncol = 2)
}

intersection_list[["nb.csv"]]

intersection_list[["nb_DEG.csv"]] <- NULL
#look for enrichment
#total number of arabidopsis genes = 27206
#total number used for DEG analysis = 21199

#gets the p-value for enrichment
phyper(sum(intersection_list[["cl14.csv"]] %in% cluster_info$Gene),length(cluster_info$Gene),21199-length(cluster_info$Gene), length(intersection_list[["cl14.csv"]]), lower.tail = F)


#I first want to see if the genes are enriched in the clusters
total_data <- c()
for ( i in names(intersection_list) ) {
  total_data <- c(total_data, i, "total",sum(intersection_list[[i]][,1] %in% cluster_info$Gene),
                  phyper(sum(intersection_list[[i]][,1] %in% cluster_info$Gene),length(cluster_info$Gene),21199-length(cluster_info$Gene), length(intersection_list[[i]]), lower.tail = F),
                  phyper(sum(intersection_list[[i]][,1] %in% cluster_info$Gene),length(cluster_info$Gene),21199-length(cluster_info$Gene), length(intersection_list[[i]]), lower.tail = T))
  #print the clusters that are not enriched at a threshold of 0.001
  if ( phyper(sum(intersection_list[[i]][,1] %in% cluster_info$Gene),length(cluster_info$Gene),21199-length(cluster_info$Gene), length(intersection_list[[i]]), lower.tail = F) > 0.001 ) {
    print(i)
  }
}


#Not all are inriched in the clusters
#181
#374
#136


#I next want to ask for each cluster is there an enrichment?
#going to use a p-value of 0.01

enrichment_vector <- c()
for ( i in names(intersection_list) ) {
  vect <- c()
  for( j in c(1,2,3,4,6) ) {
    #need # of genes to be greater than 10
    #if ( sum(intersection_list[[i]] %in% cluster_info$Gene[cluster_info$cluster == j]) < 5 ) {
   #   vect <- c(vect,NA)
    #}
    #else {
      #handles the decreasing clusters -> 1,2
    if ( j %in% c(1,2) ) {
      if ( sum(intersection_list[[i]][intersection_list[[i]][,2] == 0,1] %in% cluster_info$Gene[cluster_info$cluster == j]) < 1 ) {
        vect <- c(vect,i,j,0,1,phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 0,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                                     length(cluster_info$Gene[cluster_info$cluster == j]),
                                     21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                                     length(intersection_list[[i]][,1]), lower.tail = T))
      }
      else {
        vect <- c(vect, i,j, sum(intersection_list[[i]][intersection_list[[i]][,2] == 0,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                  phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 0,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                               length(cluster_info$Gene[cluster_info$cluster == j]),
                               21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                               length(intersection_list[[i]][,1]), lower.tail = F),
                  phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 0,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                         length(cluster_info$Gene[cluster_info$cluster == j]),
                         21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                         length(intersection_list[[i]][,1]), lower.tail = T)
                  )
      }
    }
    #handle the increasing clusters -> 3,4, 6
    else {
      if ( sum(intersection_list[[i]][intersection_list[[i]][,2] == 1,1] %in% cluster_info$Gene[cluster_info$cluster == j]) < 1 ) {
        vect <- c(vect,i,j,0,1, phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 1,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                                     length(cluster_info$Gene[cluster_info$cluster == j]),
                                     21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                                     length(intersection_list[[i]][,1]), lower.tail = T)
                  )
      }
      else {
        vect <- c(vect, i,j, sum(intersection_list[[i]][intersection_list[[i]][,2] == 1,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                  phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 1,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                               length(cluster_info$Gene[cluster_info$cluster == j]),
                               21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                               length(intersection_list[[i]][,1]), lower.tail = F),
                  phyper(sum(intersection_list[[i]][intersection_list[[i]][,2] == 1,1] %in% cluster_info$Gene[cluster_info$cluster == j]),
                         length(cluster_info$Gene[cluster_info$cluster == j]),
                         21199-length(cluster_info$Gene[cluster_info$cluster == j]), 
                         length(intersection_list[[i]][,1]), lower.tail = T)
                  )
      }
    }
   # }
    if ( j==3 ) {
      print(phyper(sum(intersection_list[[i]] %in% cluster_info$Gene[cluster_info$cluster == j]),length(cluster_info$Gene[cluster_info$cluster == j]),21199-length(cluster_info$Gene[cluster_info$cluster == j]), length(intersection_list[[i]]), lower.tail = F))
    }
  }
  #vect[vect < 0.01] <- 0
  #vect[vect > 0.01] <- 1
  
  enrichment_vector <- append(enrichment_vector, vect)
}



enrich <- matrix(enrichment_vector, ncol=5, byrow = T)
total_data <- matrix(total_data, ncol=5, byrow = T)
enrich <- rbind(enrich, total_data)
enrich <- as.data.frame(enrich)
colnames(enrich) <- c("Bact", "Cluster", "DEG_count", "Enrichment","Depletion")
enrich$Depletion <- as.numeric(as.character(enrich$Depletion))
enrich$Enrichment <- as.numeric(as.character(enrich$Enrichment))
enrich$DEG_count <- as.numeric(as.character(enrich$DEG_count))
enrich$Cluster <- relevel(enrich$Cluster, ref = "total")

rn <- c()
for ( i in enrich$Bact ) {
  rn <- c(rn,strsplit(i, "[.]")[[1]][1])
}

enrich$Bact <- rn
enrich$Bact <- as.factor(enrich$Bact)

enrich$border_col[enrich$Enrichment < 0.01] <- "Enrichment"
enrich$border_col[enrich$Depletion < 0.01] <- "Depletion"
enrich$border_col[enrich$Depletion < 0.01 & enrich$Enrichment < 0.01] <- "Messed Up"
enrich$border_col[enrich$Depletion > 0.01 & enrich$Enrichment > 0.01] <- "Unsignificant"

#handle 374 and 136 (cluster is insignificant)
enrich$border_col[enrich$Bact %in% c("MF374", "MF136")] <- "Unsignificant"


#superheat(enrich,
#          row.dendrogram = T,
#          heat.pal = c("black", "white"),
#          left.label.col = "white",
#          bottom.label.col = "white",
#          heat.lim = c(0,.1)
          #scale = T
#          grid.hline.col = "light grey",
#          grid.vline.col = "light grey"
#          )

library(ggplot2)

vect <- c()
for ( i in unique(enrich$Bact) ) {
  vect <- c(vect, enrich$DEG_count[enrich$Bact == i], enrich$Enrichment[enrich$Bact == i])
}

count_mat <- matrix(vect, ncol=12, byrow = T)
rownames(count_mat) <- unique(enrich$Bact)

clust <- hclust(dist(log2(count_mat[,c(1:5)]+1)))
plot(clust)
enrich$Bact <- factor(as.character(enrich$Bact), levels = clust$labels[clust$order])

cluster_clust <- hclust(cor(t(count_mat[,c(1:5,7:11)])))
plot(cluster_clust)

ggplot(enrich[enrich$Cluster != "total",], aes(x=Cluster, y=Bact)) + 
  geom_raster(aes(fill=DEG_count))+
  geom_tile(alpha=0,aes(color=border_col,width=.99, height=.85), size=.7)+
  scale_fill_gradientn(colors = c("white", "dodgerblue"))+
  scale_colour_manual(values=c( "Black", "Light Grey"))+
  theme_minimal()


tiff("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/cluster_figs/nr_hyperg_scaled_p_border.tiff",
     width = 700, height = 700)
superheat(enrich,
               row.dendrogram = T,
               heat.pal = c("Black", "White"),
               left.label.col = "white",
               bottom.label.col = "white",
               #          grid.hline.col = "light grey",
               #          grid.vline.col = "light grey"
)
ggplot(enrich[enrich$Cluster != "total",], aes(x=Cluster, y=Bact)) + 
  geom_raster(aes(fill=log2(DEG_count+1)))+
  geom_tile(alpha=0,aes(color=border_col,width=.99, height=.85), size=.7)+
  scale_fill_gradientn(colors = c("white", "dodgerblue"))+
  scale_colour_manual(values=c( "Black", "Light Grey"))+
  theme_minimal()
dev.off()

#ggpoint of the suppressors and inducers
suppressors <- c("MF105", "MF374", "MF345", "MF314", "MF79", "MF138")
inducers <- c("MF41", "cl18", "MF161", "MF370", "MF50", "cl69")

int_only <- enrich[enrich$Bact %in% c(suppressors, inducers),]
int_only$supp[int_only$Bact %in% suppressors] <- "Suppressor"
int_only$supp[int_only$Bact %in% inducers] <- "Inducer"

library(ggpubr)
ggplot(int_only[int_only$Cluster != "total",], aes(Cluster, DEG_count, color=supp)) +
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge())+
  labs(color="", y="DEG Count") +
  scale_color_manual(values = c(Inducer="dodgerblue", Suppressor="firebrick"))+
  theme_minimal()
  stat_compare_means(method = "t.test",aes(label = ..p.signif..))
 