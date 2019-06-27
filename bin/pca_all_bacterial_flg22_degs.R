#This code is the main code used to create PCAs with the monoassociation using the DEGs identified in the original experiment
# Perform a GO ennrichment on the highly contributing genes

#the working directory ie the path to the data directory within MAE
setwd("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/MAE/data/")

#Read in Paulo's data and metadata
rna_data <- read.delim("MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("design_MAE_complete.txt",sep="\t", stringsAsFactors = F)

#read in the results
all_degs <- read.csv("nb_DEG_orig.csv")
#filter the DEGs to only retain the most significant
all_degs <- all_degs[all_degs$FDR < 0.01 & (all_degs$logFC < -log2(1.5) | all_degs$logFC > log2(1.5) ),]

get_fold_change_data <- function(){
#I want cluster the log fold change seen by all these bacteria vs nb
lfc_of_all_genes <- c()
for ( i in list.files("mono_degs/nbm_vs_bactm")) {
  mm_file <- read.csv(paste("mono_degs/nbm_vs_bactm/", i, sep = ""))
  mm_file <- mm_file[mm_file$gene %in% all_degs$genes,]
  #makes sure the genes are all in order
  mm_file <- mm_file[all_degs$genes,]
  #get the gene values
  lfc_of_all_genes <- append(lfc_of_all_genes, c(strsplit(i, ".cs")[[1]][1], mm_file$logFC) )
}

#create a matrix out the resulting data where each row is a bacterial strain
lfc_mat <- matrix(lfc_of_all_genes, ncol = length(all_degs$genes)+1, byrow = T)
colnames(lfc_mat) <- c("Bacteria", all_degs$genes)
rownames(lfc_mat) <- lfc_mat[,1]
lfc_mat <- lfc_mat[2:ncol(lfc_mat)]
return(lfc_mat)
}
get_zscore <- function(count_df) {
  for ( i in 1:ncol(count_df)) {
    count_df[i,] <- (count_df[i,]-mean(count_df[i,]))/sd(count_df[i,])
  }
  return(count_df)
}

get_cpm_data <- function() {
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
  data <- data[,colnames(data) %in% all_degs$genes]
  
  return(data)
}

#I then want to perform a PCA comparing the LFC of all these genes compared to nb-
library("FactoMineR")
library("factoextra")
#install.packages("")

#for fold change analysis
#get_fold_change_data()

#For cpm analysis
lfc_mat <- get_cpm_data()

lfc_mat <- lfc_mat[which(!(row.names(lfc_mat) %in% c("MF181.flg22minus", "MF181.flg22plus"))),]
rownames(lfc_mat) <- sapply(rownames(lfc_mat), function(x){strsplit(x,".flg")[[1]][1]})
res.pca <- PCA(lfc_mat, scale.unit = F, ncp = 13, graph = F)
print(res.pca)

# We want to keep the PCs that have a eigenvalue > 1 or you can find the PC's that account
# for a certain percentage of variance
eig.val <- get_eigenvalue(res.pca)

eig.val

#extract the variable results which you can make a couple plots with
#var$coord: coordinates of variables to create a scatter plot
#var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
#var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
var <- get_pca_var(res.pca)
var$contrib

#
fviz_screeplot(res.pca)

#Look at the contributions of each variable to each PC
#Change the axes to get a different PC
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes =1, top = 50)
#var$contrib can be used to extract the contributions of each gene to the PC
var$contrib

fviz_contrib(res.pca, choice = "var", axes =2, top = 32)
#var$contrib can be used to extract the contributions of each gene to the PC
top_ax2 <- var$contrib[var$contrib[,2] >= (100/480),2]

write.csv(top_ax1, file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/dim_1_top.csv")
write.csv(top_ax2, file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/dim_2_top.csv")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 50)

fviz_cos2(res.pca, choice = "var", axes = 1:2)

#This gets the most significant variables for each of the dimensions in axes
res.desc <- dimdesc(res.pca, axes = c(1:3), proba = 0.01)
res.desc$Dim.1

library(kernlab)
sc <- specc(res.pca, centers=2)

c_index <- rep("N", nrow(lfc_mat))
row.names(lfc_mat)
c_index[seq(1,nrow(lfc_mat), by=2)] <- "M"
c_index[seq(2,nrow(lfc_mat), by=2)] <- "P"
row.names(res.pca$ind)
fviz_pca_ind(res.pca, axes = c(1,2), col.ind = c_index, repel = T, invisible="quali", palette = c("Firebrick", "Dodgerblue"), title="", addEllipses = T, ellipse.type = "confidence", elipse.level=.5, mean.point=F)+
  theme_minimal()
fviz_pca_ind(res.pca, axes = c(1,3), col.ind = c_index, repel = T, invisible="quali", palette = c("Firebrick", "Dodgerblue"), title="", addEllipses = T, ellipse.type = "confidence", elipse.level=.5, mean.point=F)+
  theme_minimal()

fviz_pca_ind(res.pca, axes = c(1,2), col.ind = as.factor(sc), repel = T, invisible="quali", title="", mean.point=F)+
  theme_minimal()

## Try using mclust
library(mclust)

l <- Mclust(res.pca$ind$coord[,1:2], G=3)


#plot a beautiful PCA
library(ggfortify)
library(cluster)
library(ggplot2)
library(ggrepel)
#plotting_data
rownames(lfc_mat) <- sapply(rownames(lfc_mat), function(x){strsplit(x,".flg")[[1]][1]})
plot_metad <- data.frame(rownames(lfc_mat))
plot_metad$flg22 <- c_index
colnames(plot_metad) <- c("Names", "Flg22")
row.names(plot_metad) <- plot_metad$Names

ggplot(res.pca$ind$coord, aes(x=Dim.1, y=Dim.3, label=row.names(res.pca$ind$coord), color=c_index))+
  geom_point()+
  geom_text_repel()+
  stat_ellipse(inherit.aes = F, data=res.pca$ind$coord,aes(x=Dim.1, y=Dim.2, color=c_index), level = .7, type = "t", linetype=1)


autoplot(prcomp(lfc_mat), data = plot_metad, shape=F, label = FALSE, label.size = 3)+
  geom_text_repel(aes(color=plot_metad$Flg22, label=plot_metad$Names), position=position_jitter())+
  theme_bw()

write.csv(res.desc$Dim.1, "/Users/nicholascolaianni/Desktop/pca_dimension1_cpm.csv", row.names = T)
write.csv(res.desc$Dim.2, "/Users/nicholascolaianni/Desktop/pca_dimension2_cpm.csv", row.names = T)
write.csv(res.desc$Dim.3, "/Users/nicholascolaianni/Desktop/pca_dimension3_cpm.csv", row.names = T)


### PLOTTING ###


### Warning -> load is actually the contribution 
### I did not change it due to having to change so much code. real_load is 
### actually the load
#BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "tair_locus", "tair_symbol"), mart = mart)


#get the top contributors of the 1st axis
top_bot_1 <- sort(res.pca$var$contrib[,1])
top_bot_1 <- as.data.frame(top_bot_1)
colnames(top_bot_1) <- "load"
#top_bot_1$load <- top_bot_1$load/sum((top_bot_1$load))
top_bot_1$names <- row.names(top_bot_1)

#get load values
top_bot_1$real_load <- sapply(top_bot_1$names, function(x){res.pca$var$coord[which(row.names(res.pca$var$coord) == x),1]})
top_bot_1$real_load[top_bot_1$real_load>0] <- "Pos"
top_bot_1$real_load[top_bot_1$real_load<0] <- "Neg"

#calculate_pvals with z score and standard error
top_bot_1$zscore <- (mean(top_bot_1$load) - top_bot_1$load)/(sd(top_bot_1$load)/sqrt(nrow(top_bot_1)))
#two-tailed test
top_bot_1$pval <- 2*pt(-abs(top_bot_1$zscore), df=nrow(top_bot_1)-1)
top_bot_1$gname <- sapply(top_bot_1$names, function(x){gene_info$tair_symbol[which(gene_info$ensembl_gene_id == x)]})

top_bot_1$gname[top_bot_1$gname == ""] = top_bot_1$names[which(top_bot_1$gname=="")]

#getting my significant vals based on the data
qnorm(.95, mean=(mean(top_bot_1$load)),sd = sd(top_bot_1$load) )
qnorm(.05, mean=(mean(top_bot_1$load)),sd = sd(top_bot_1$load) )
library(ggrepel)
library(RColorBrewer)

#I am usinng a log normal distribution to calculate significance at a 95% connfidence interval
ln_mean <- log(mean(top_bot_1$load)^2 / sqrt(sd(top_bot_1$load)^2 + mean(top_bot_1$load)^2))
ln_sd <- sqrt(log(1 + (sd(top_bot_1$load)^2 / mean(top_bot_1$load)^2)))


ggplot(top_bot_1, aes(x=names,y=load))+
  geom_point()+
  geom_hline(yintercept=qlnorm(.95, mean=ln_mean,sd = ln_sd ), linetype="dashed", color = "red")+
  geom_text_repel(data = top_bot_1[top_bot_1$load>qlnorm(.95, mean=ln_mean,sd = ln_sd),], aes(x=names, label=gname, color=real_load))+
  geom_hline(yintercept = (sum(abs(top_bot_1$load))/nrow(top_bot_1)), linetype="dashed", color="blue")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Contribution Percentage")

#second axis plot
#this performs the same as above but on the second axis
top_bot_2 <- sort(res.pca$var$contrib[,2])
top_bot_2 <- as.data.frame(top_bot_2)
colnames(top_bot_2) <- "load"
#top_bot_1$load <- top_bot_1$load/sum((top_bot_1$load))
top_bot_2$names <- row.names(top_bot_2)

top_bot_2$real_load <- sapply(top_bot_2$names, function(x){res.pca$var$coord[which(row.names(res.pca$var$coord) == x),2]})
top_bot_2$real_load[top_bot_2$real_load>0] <- "Pos"
top_bot_2$real_load[top_bot_2$real_load<0] <- "Neg"

top_bot_2$gname <- sapply(top_bot_2$names, function(x){gene_info$tair_symbol[which(gene_info$ensembl_gene_id == x)]})

top_bot_2$gname[top_bot_2$gname == ""] = top_bot_2$names[which(top_bot_2$gname=="")]

ln_mean_2 <- log(mean(top_bot_2$load)^2 / sqrt(sd(top_bot_2$load)^2 + mean(top_bot_2$load)^2))
ln_sd_2 <- sqrt(log(1 + (sd(top_bot_2$load)^2 / mean(top_bot_2$load)^2)))

ggplot(top_bot_2, aes(x=names,y=load))+
  geom_point()+
  geom_hline(yintercept=qlnorm(.95, mean=ln_mean_2,sd = ln_sd_2 ), linetype="dashed", color = "red")+
  geom_text_repel(data = top_bot_2[top_bot_2$load>qlnorm(.95, mean=ln_mean_2,sd = ln_sd_2 ),], aes(x=names, label=gname, color=real_load))+
  geom_hline(yintercept = (sum(abs(top_bot_2$load))/nrow(top_bot_2)), linetype="dashed", color="blue")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Contribution Percentage")

#get the total % contribution of the genes greater than the blue lines and red lines
#1st dimension top normal
sum(res.pca$var$contrib[rownames(res.pca$var$contrib) %in% top_bot_1$names[top_bot_1$load>qlnorm(.95, mean=ln_mean,sd = ln_sd)],1])
#34.76502% of dim 1 is explained by the top variables represented by the red line
sum(res.pca$var$contrib[rownames(res.pca$var$contrib) %in% top_bot_1$names[top_bot_1$load > (sum(abs(top_bot_1$load))/nrow(top_bot_1))],1])
#71.958% of dim1 is explained by the variables greater/less that the uniform (blue dashed lines)


#2nd dimension top normal
sum(res.pca$var$contrib[rownames(res.pca$var$contrib) %in% top_bot_2$names[top_bot_2$load>qlnorm(.95, mean=ln_mean_2,sd = ln_sd_2 )],2])
#53.168% of dim2 is explained by the top variables represented by the red lines
sum(res.pca$var$contrib[rownames(res.pca$var$contrib) %in% top_bot_2$names[top_bot_2$load > (sum(abs(top_bot_2$load))/nrow(top_bot_2))],2])
#84.10899% of dim2 is explained by the variables greater/less that the uniform (blue dashed lines)

#Get the ids of the genes for a go analysis of the genes greater than the blue lines
write.csv(top_bot_1$names[top_bot_1$load > (sum(abs(top_bot_1$load))/nrow(top_bot_1))], file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/1st_pca_uniform_ids.csv")
write.csv(top_bot_2$names[top_bot_2$load > (sum(abs(top_bot_2$load))/nrow(top_bot_2))], file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/2nd_pca_uniform_ids.csv")

write.csv(top_bot_1$names[top_bot_1$load>qlnorm(.95, mean=ln_mean,sd = ln_sd)], file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/1st_pca_lnorm_ids.csv")
write.csv(top_bot_2$names[top_bot_2$load>qlnorm(.95, mean=ln_mean_2,sd = ln_sd_2 )], file = "/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/Paulo_MAE_figs/pca_data_figs/2nd_pca_lnorm_ids.csv")

###################################################################
### Is it worth to perform the same analysis on the third axis? ###
###################################################################

#Re-make the PCA with the top componants from the 1st PCA
new_pca <- lfc_mat[,colnames(lfc_mat) %in% top_bot_1$names[top_bot_1$load > (sum(abs(top_bot_1$load))/nrow(top_bot_1))] ]

new.res.pca <- PCA(new_pca, scale.unit = F, ncp = 5, graph = F)

fviz_pca_ind(new.res.pca, axes = c(1,2), col.ind = c_index, repel = T, invisible="quali", palette = c("Dodgerblue", "Firebrick"), title="")+
  theme_minimal()


#Go enrichment analysis
library(clusterProfiler)
library(AnnotationHub)

# org.At.tair.db arabidopsis db

#first we build the arabidopsis dataset
hub <- AnnotationHub()
q <- query(hub, "Arabidopsis thaliana")
id <- q$ah_id[length(q)]
arab <- hub[[id]]

ego1 <- enrichGO(gene = top_bot_1$names[top_bot_1$load > (sum(abs(top_bot_1$load))/nrow(top_bot_1))],
         OrgDb = arab,
         ont = "BP",
         pAdjustMethod = "BH",
         keyType = 'TAIR',
         pvalueCutoff = 0.01,
         qvalueCutoff = 0.05
         )
ego2 <- enrichGO(gene = top_bot_2$names[top_bot_2$load > (sum(abs(top_bot_2$load))/nrow(top_bot_2))],
                 OrgDb = arab,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 keyType = 'TAIR',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01
)

dotplot(ego1, showCategory = 20)
dotplot(ego2, showCategory = 20)

ego1_lnorm <- enrichGO(gene = top_bot_1$names[top_bot_1$load>qlnorm(.95, mean=ln_mean,sd = ln_sd)],
                 OrgDb = arab,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 keyType = 'TAIR',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05
)
ego2_lnorm <- enrichGO(gene = top_bot_2$names[top_bot_2$load>qlnorm(.95, mean=ln_mean_2,sd = ln_sd_2 )],
                 OrgDb = arab,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 keyType = 'TAIR',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05
)

dotplot(ego1_lnorm)
dotplot(ego2_lnorm)

  