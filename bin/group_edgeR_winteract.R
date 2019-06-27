#this edgeR script will create a model with all the data and make the comparisons of interest
library(edgeR)

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)

#trying the long grouping
rna_metadata <- rna_metadata[rna_metadata$Bacteria!="Hksyncom",]
cols_to_keep <- which(colnames(rna_data) %in% rna_metadata[,1])
rna_data <- rna_data[,c(1,2,cols_to_keep)]


full_dgelist <- DGEList(counts = rna_data[,3:ncol(rna_data)], lib.size = colSums(rna_data[,3:ncol(rna_data)]), samples = rna_metadata, genes=rna_data$Gene)
#filter low count reads with cpm
countsPerMillion <- cpm(full_dgelist)
#creates a T and F matrix based on CPM
countCheck <- countsPerMillion > 1
#creates a vector of the rows that meet the criteria
keep <- which(rowSums(countCheck) >= 8)
full_dgelist <- full_dgelist[keep,] #filter out low count data
full_dgelist <- calcNormFactors(full_dgelist)


#I am going to make the model matrix portions
flg22 <- factor(rna_metadata$Treatment, levels = c("flg22minus", "flg22plus"))
bacteria <- factor(rna_metadata$Bacteria)
bacteria <- relevel(bacteria, ref = "nb")
experiment <- factor(rna_metadata$Experiment, levels = c("e1", "e2", "e3"))
group <- factor(rna_metadata$Group)

designMat <- model.matrix(~bacteria*flg22+experiment:bacteria)

full_dgelist <- estimateDisp(full_dgelist, design = designMat)
#dgList <- estimateGLMCommonDisp(full_dgelist, design=designMat) #uses a common estimate across genes
#dgList <- estimateGLMTrendedDisp(full_dgelist, design=designMat) #fits an estimate based on mean variance trends across the data
#dgList <- estimateGLMTagwiseDisp(full_dgelist, design=designMat) #compute a genewise variance

#perform dgexpression
fit <- glmQLFit(full_dgelist, designMat)

my_cont_vect <- rep(0, length(colnames(designMat)))

nb_vect <- my_cont_vect
colnames(designMat)
nb_vect[39:74] <- 1

lrt <- glmQLFTest(fit, contrast = nb_vect)

#load up the r object fit
fit <- readRDS("/Users/nicholascolaianni/Desktop/glmqlfit_fullmatrix")

degenes_btwn_samples <- c()
for ( i in unique(bacteria) ) {
  my_cont_vect <- rep(0, length(colnames(designMat)))
  my_cont_vect[which(colnames(designMat) ==  paste("bacteria", i, ":flg22flg22plus", sep = ""))] <- 1
  if ( i == "nb" ) {
    next
  }
  lrt <- glmQLFTest(fit, contrast = my_cont_vect)
  w <- topTags(lrt, n=nrow(lrt), sort.by = "none")
  degenes_btwn_samples <- append( degenes_btwn_samples, sum(as.logical(decideTestsDGE(lrt, p.value = .01, lfc = log2(1.5)))) )
  write.csv(w$table, file = paste("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/flg22plus_bactvsflg22/", i, ".csv", sep = ""), quote = F, row.names = T)
}
names(degenes_btwn_samples) <- unique(bacteria)[unique(bacteria) != "nb"]

write.csv(matrix(c(degenes_btwn_samples, names(degenes_btwn_samples)), ncol = 2, byrow = F),"/Users/nicholascolaianni/Desktop/flg22plus_bactvsflg22.csv" )


w <- topTags(lrt, n=nrow(lrt), sort.by = "none", p.value = 0.01)
write.csv(w$table, file = "/Users/nicholascolaianni/Desktop/anova_flg22_fullflg22plus.csv", quote = F, row.names = T)

