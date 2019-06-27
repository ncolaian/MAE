#This script will be for plotting the expected number of DE genes from each cluster and how many we actually see

cluster_info <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/clusters_IDs_IRM_experiment.txt", stringsAsFactors = F)
table(cluster_info$cluster)

top_nbflg22 <- read.csv("/Users/nicholascolaianni/Desktop/nb_DEG.csv", stringsAsFactors = F)
top_nbflg22 <- top_nbflg22[top_nbflg22$logFC > log2(1.5) | top_nbflg22$logFC < -(log2(1.5)),]

#Get the original data from paulo experiment
nb_first <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/original_nbflg22degs/edgeR_glm_question2_NB_flg22_response.csv", stringsAsFactors = F)
nb_first <- nb_first[nb_first$logFC>log2(1.5) | nb_first$logFC< -(log2(1.5)),]
nb_first <- nb_first[nb_first$FDR < 0.01,]

#get expected data
nb_first$cluster <- sapply(nb_first$Gene, function(x){cluster_info$cluster[cluster_info$Gene == x]})
s <- table(nb_first$cluster) 
s <- s[c(1,2,3,4,6)]
s <- s/sum(s)

table

#get observed data
top_nbflg22$cluster <- sapply(top_nbflg22$genes, function(x){cluster_info$cluster[cluster_info$Gene == x]})
top_nbflg22$cluster[is.na(top_nbflg22$cluster >= 1)] <- 0
o <- table(as.numeric(top_nbflg22$cluster))
o <- c(17,35,124,66,46)
ex <- chisq.test(o, p=s)
ex$expected
ex$p.value
sum(o)

t <- data.frame(c(o, ex$expected),c(rep("o",5), rep("e",5)), as.character(c(1,2,3,4,6)))
colnames(t) <- c("counts", "part", "cluster")
library(ggplot2)

ggplot(t, aes(x=cluster, y=counts, fill=part))+
  geom_bar(stat = "Identity", position = "dodge")+
  labs(x="Cluster", y="Number of DEGs", fill="")+
  ggtitle("Count of DEGs From Each Cluster")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

dbinom(124, 288, s[3])

