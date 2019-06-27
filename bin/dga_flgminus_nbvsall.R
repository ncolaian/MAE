#this edgeR script will create a model with all the data and make the comparisons of interest

#This particular script only uses the flg22- data


library(edgeR)

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)


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

#trying to only analyze the bacteria without flg22 data
rna_metadata <- rna_metadata[rna_metadata$Bacteria!="Hksyncom" & rna_metadata$Treatment!="flg22plus",]
cols_to_keep <- which(colnames(mm_data) %in% rna_metadata[,1])
mm_data <- mm_data[,c(1,2,cols_to_keep)]

full_dgelist <- DGEList(counts = mm_data[,3:ncol(mm_data)], lib.size = colSums(mm_data[,3:ncol(mm_data)]), samples = rna_metadata, genes=mm_data$Gene)


#I am going to make the model matrix portions
bacteria <- factor(rna_metadata$Bacteria)
bacteria <- relevel(bacteria, "nb")
experiment <- factor(rna_metadata$Experiment, levels = c("e1", "e2", "e3"))

designMat <- model.matrix(~bacteria+experiment)

full_dgelist <- estimateDisp(full_dgelist, design = designMat)
#dgList <- estimateGLMCommonDisp(full_dgelist, design=designMat) #uses a common estimate across genes
#dgList <- estimateGLMTrendedDisp(full_dgelist, design=designMat) #fits an estimate based on mean variance trends across the data
#dgList <- estimateGLMTagwiseDisp(full_dgelist, design=designMat) #compute a genewise variance

#perform dgexpression
fit <- glmFit(full_dgelist, designMat)
colnames(designMat)

lrt_list <- list()
for ( i in unique(bacteria) ) {
  if ( i == "nb") {
    next
  }
  print(i)
  lrt_list[[i]] <- glmLRT(fit, coef = paste("bacteria", i, sep = ""))
}
#lrt <- glmLRT(fit, coef = "bacteriaMF181")

library(ggplot2)
num_degs <- sapply(names(lrt_list), function(x){length(rownames(lrt_list[[x]])[as.logical(decideTestsDGE(lrt_list[[x]], p=0.01, lfc = log2(1.5)))])})

degenes_btwn_samples <- as.data.frame(matrix(data=c(names(num_degs), num_degs),ncol=2 ))
degenes_btwn_samples$V2 <- as.numeric(as.character(degenes_btwn_samples$V2))
ggplot(degenes_btwn_samples, aes(x=V1, y=V2)) +
  geom_bar(stat = "identity")+
  ylab("Number of DEGs")+
  xlab("Bacteria") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))
colnames(degenes_btwn_samples) <- c("Bacteria", "DEGs")
write.csv(degenes_btwn_samples, file="/Users/nicholascolaianni/Desktop/nb_vs_bacteria_split_together.csv", quote = F, row.names = F)

w <- topTags(lrt, n=nrow(lrt), sort.by = "none", p.value = 0.01)
hist(w$table$FDR)

