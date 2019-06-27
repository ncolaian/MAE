# This script is to model nb+ vs bact+

library(edgeR)

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)

separate_data_plus_compare <- function(name, data) {
  samples_to_keep <- rna_metadata[(rna_metadata$Bacteria == name & rna_metadata$Treatment == "flg22plus") | (rna_metadata$Bacteria == "nb" & rna_metadata$Treatment == "flg22plus"),]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
  return(list(mm_data, samples_to_keep))
}

main_edgeR_DGE_plus_compare <- function() {
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
  
  lrt_list <- list()
  for ( i in unique(rna_metadata$Bacteria) ) {
    if ( i %in% c( "nb" )) {
      next
    }
    print(i)
    mm_dataframes <- separate_data_plus_compare(i, mm_data) #getting each samples data frame
    #create a DGE object containing the counts and metadata
    dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
    #get the norm factors from the full dataframe
    dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
    
    #model matrix design
    bact <- factor(dgList$samples$Bacteria, levels= c("nb", i))
    experiment <- factor(dgList$samples$Experiment, levels = c("e1", "e2", "e3"))
    designMat <- model.matrix(~bact + experiment)
    
    #must do this for HKsyncom because we do not have 3rd experiment
    #we can also assum the randomness will be the same to no bacteria
    if ( i == "Hksyncom" ) {
      designMat <- model.matrix(~bact + experiment)
    }
    
    #fit dispersion parameter
    dgList <- estimateGLMCommonDisp(dgList, design=designMat) #uses a common estimate across genes
    dgList <- estimateGLMTrendedDisp(dgList, design=designMat) #fits an estimate based on mean variance trends across the data
    dgList <- estimateGLMTagwiseDisp(dgList, design=designMat) #compute a genewise variance
    
    #perform dgexpression
    fit <- glmFit(dgList, designMat)
    fit$coefficients
    #Perform a logratiotest with pvals and cpm
    #the coef makes sure we are testing for DGE for the addition of bact in flg+ conditions
    bact_comp <- paste("bact", i, sep = "")
    lrt <- glmLRT(fit, coef=bact_comp) 
    
    
    #same the resulting data in a list
    lrt_list[[i]] <- lrt
    #lrt_list[[i]] <- lrt
    
  }
  
  #dev.off()
  return(lrt_list)
}


#### MAIN ####
output_lrt <- main_edgeR_DGE_plus_compare()
#### PLOTTING EdgeR RESULTS ####
library(ggplot2)

num_degs <- sapply(names(output_lrt), function(x){length(rownames(output_lrt[[x]])[as.logical(decideTestsDGE(output_lrt[[x]], p=0.01, lfc = log2(1.5)))])})

degenes_btwn_samples <- as.data.frame(matrix(data=c(names(num_degs), num_degs),ncol=2 ))
degenes_btwn_samples$V2 <- as.numeric(as.character(degenes_btwn_samples$V2))
ggplot(degenes_btwn_samples, aes(x=V1, y=V2)) +
  geom_bar(stat = "identity")+
  ylab("Number of DEGs")+
  xlab("Bacteria") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))

for ( i in names(output_lrt) ) {
  w <- topTags(output_lrt[[i]], nrow(output_lrt[[i]]), sort.by = "none")
  write.csv(w$table, file = paste("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/old_model/final_sep_data/nbp_vs_bactp/", i, ".csv", sep = ""), quote = F, row.names = T)
}

write.csv(degenes_btwn_samples, file="/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/old_model/final_sep_data/nbp_vs_bactp/numDEGS_plus.csv", quote = F)

