# I want to test if flg22 response is correlated with root length data

root_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_mono/root_length_all_reps_filtered.txt", stringsAsFactors = F)
meta_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_mono/paulo_35_experiment.csv", stringsAsFactors = F)
meta_data$taxon_oid <- as.character(meta_data$taxon_oid)

root_data$Bacteria <- factor(root_data$Bacteria)
root_data$Bacteria <- relevel(root_data$Bacteria, "NB")
root_data$Flg <- factor(root_data$Flg)
root_data$Flg <- relevel(root_data$Flg, "minus")
root_data$Replicate <- factor(root_data$Replicate)
root_data$Replicate <- relevel(root_data$Replicate, "R1")


#build a linear model
root_lm <- lm(Length_plant~Bacteria+Replicate, data = root_data[root_data$Flg == "minus",])
anova(root_lm)
summary(root_lm)

#data from inital lm
root_lm$coefficients

coefs <- c()
for ( i in as.character(num_degs[,3]) ) {
  if ( i == "NB") {
    next
  }
  if ( strsplit(i, "l")[[1]][1] == "c") {
    i <- toupper(i)
  }
  else if (strsplit(i, "c")[[1]][1] == "e." ) {
    i <- "E.coli"
  }
  else if (strsplit(i, "y")[[1]][1] == "s" ) {
    i <- "Syncom"
  }
  else if (strsplit(i, "k")[[1]][1] == "H" ) {
    i <- "HK"
  }
  coef_label <- paste("Bacteria", i, sep = "")
  coefs <- append(coefs, root_lm$coefficients[coef_label])
}

coefs
num_degs <- read.csv("/Users/nicholascolaianni/Desktop/glmqlfit_nb_bact.csv")
bum_degs2 <- read.csv("/Users/nicholascolaianni/Desktop/bact_vs_nb_numDEGS_plus.csv")
cor(num_degs[,3], bum_degs2[,3], method = "spearman")
cor(num_degs[,2], abs(coefs), method = "spearman")
plot(log(num_degs[,3]), abs(coefs))
#This tells me that not all large root differences are due to large rnaseq changes in the plant
