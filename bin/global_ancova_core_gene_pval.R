#this will be an attempt to perform bootstrapping on the core gene data. I will be attempting to compare
#the medians of groups that may be biased due to correllation

library(GlobalAncova)
library(edgeR)

#Read in Paulo's data and metadata
rna_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/MAE_mRNA_counts_filtered.txt", sep="\t", stringsAsFactors = F)
rna_metadata <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/design_MAE_complete.txt",sep="\t", stringsAsFactors = F)
core_glist <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/cor_flg22/core_flg22.txt", stringsAsFactors = F)

#this function will use rna_data and rna_metadata put in the path above
separate_data_meta <- function( name, data ) {
  samples_to_keep <- rna_metadata[rna_metadata$Bacteria==name,]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
  return(list(mm_data, samples_to_keep))
}
separate_data_no_bact <- function(name, data) {
  samples_to_keep <- rna_metadata[(rna_metadata$Bacteria == name & rna_metadata$Treatment == "flg22minus") | (rna_metadata$Bacteria == "nb" & rna_metadata$Treatment == "flg22minus"),]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
  return(list(mm_data, samples_to_keep))
}
separate_data_nbp_bactm <- function(name, data) {
  samples_to_keep <- rna_metadata[(rna_metadata$Bacteria == name & rna_metadata$Treatment == "flg22minus") | (rna_metadata$Bacteria == "nb" & rna_metadata$Treatment == "flg22plus"),]
  cols_to_keep <- which(colnames(data) %in% samples_to_keep[,1])
  mm_data <- data[,c(1,2,cols_to_keep)]
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
  
  lrt_list <- list()
  for ( i in unique(rna_metadata$Bacteria) ) {
    print(i)
    mm_dataframes <- separate_data_meta(i, mm_data) #getting each samples data frame
    #create a DGE object containing the counts and metadata
    dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
    #get the norm factors from the full dataframe
    dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
   
    gene_table <- cpm(dgList)
    row.names(gene_table) <- dgList$genes$genes
    ancova_res <- GlobalAncova(xx=gene_table, formula.full = ~Treatment+Experiment, formula.red = ~Experiment, model.dat = dgList$samples, test.genes = core_glist$Core, method = "both", perm = 10000) 
    
    #same the resulting data in a list
    lrt_list[[unique(dgList$samples$Bacteria)[1]]] <- ancova_res
    
  }
  
  #dev.off()
  return(lrt_list)
}

main_edgeR_DGE_bact <- function() {
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
    print(i)
    if ( i == "nb" ) {
      next
    }
    mm_dataframes <- separate_data_no_bact(i, mm_data) #getting each samples data frame
    #create a DGE object containing the counts and metadata
    dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
    #get the norm factors from the full dataframe
    dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
    
    gene_table <- cpm(dgList)
    row.names(gene_table) <- dgList$genes$genes
    if ( i == "Hksyncom" ) {
      ancova_res <- GlobalAncova(xx=gene_table, formula.full = ~Bacteria+Experiment, formula.red = ~Experiment, model.dat = dgList$samples, test.genes = core_glist$Core, method = "both", perm = 10000) 
    }
    else{
      ancova_res <- GlobalAncova(xx=gene_table, formula.full = ~Bacteria+Experiment, formula.red = ~Experiment, model.dat = dgList$samples, test.genes = core_glist$Core, method = "both", perm = 10000) 
    }
    #same the resulting data in a list
    lrt_list[[unique(dgList$samples$Bacteria)[1]]] <- ancova_res
    
  }
  
  #dev.off()
  return(lrt_list)
}

main_edgeR_DGE_nbplus_bactminus <- function() {
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
    print(i)
    if ( i == "nb" ) {
      next
    }
    mm_dataframes <- separate_data_nbp_bactm(i, mm_data) #getting each samples data frame
    #create a DGE object containing the counts and metadata
    dgList <- DGEList(counts = mm_dataframes[[1]][,3:ncol(mm_dataframes[[1]])], samples = mm_dataframes[[2]], genes = mm_dataframes[[1]][,1], group = mm_dataframes[[2]]$Treatment )
    
    #get the norm factors from the full dataframe
    dgList$samples$norm.factors <- sapply(dgList$samples$Sample.., function(x){full_dgelist$samples$norm.factors[full_dgelist$samples$Sample.. == x]})
    
    gene_table <- cpm(dgList)
    row.names(gene_table) <- dgList$genes$genes
    ancova_res <- GlobalAncova(xx=gene_table, formula.full = ~Treatment+Experiment, formula.red = ~Experiment, model.dat = dgList$samples, test.genes = core_glist$Core, method = "both", perm = 10000) 
    #same the resulting data in a list
    lrt_list[[unique(dgList$samples$Bacteria)[1]]] <- ancova_res
    
  }
  
  #dev.off()
  return(lrt_list)
}

ancova_out <- main_edgeR_DGE_flg()
ancova_out <- main_edgeR_DGE_bact()
ancova_out <- main_edgeR_DGE_nbplus_bactminus()

#get the pvalues - 2 is the permutation and 3 is the 
for ( i in names(ancova_out) ) {
  print(i)
  print(ancova_out[[i]]$test.result[2])
  print(ancova_out[[i]]$test.result[3])
}
