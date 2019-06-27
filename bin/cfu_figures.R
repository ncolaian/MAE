#This script will look at differences in CFU counts

library(phylolm)
library(ape)

cfu_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/cfu_data/growth_data_mono_rep1.txt", sep = "\t")

tree_syn <- read.tree("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/35_mem_tree.txt")

cfu_data <- cfu_data[cfu_data$class != "control",]
cfu_data$class <- as.character(cfu_data$class)
cfu_data$class[cfu_data$class == "ns"] <- "Inducer"
cfu_data$class[cfu_data$class == "s"] <- "Supressor"
cfu_data$class <- as.factor(cfu_data$class)
#first run a simple LM to see if there is any significance between class and cfu
simple_lm <- lm(Log.CFUs.mg. ~ class + taxa, data=cfu_data)
summary(simple_lm)
anova(simple_lm)
#There does not seem to be any significance


#run a phylogenetic model
#prune tree
mean_data <- c()
for ( i in tree_syn$tip.label ) {
  if (!( i %in% cfu_data$ID ) ) {
    tree_syn <- drop.tip(tree_syn, i)
  }
  else {
    mean_data <- append(mean_data, c(i, mean(cfu_data$Log.CFUs.mg.[cfu_data$ID == i]), as.character(cfu_data$class[cfu_data$ID == i][1])))
  }
}
t <- matrix(mean_data, ncol = 3, byrow = T)
t <- as.data.frame(t)
t$V2 <- as.numeric(as.character(t$V2))
rownames(t) <- t$V1
phylo_lm <- phylolm(V2 ~ V3, data = t, tree_syn)
summary(phylo_lm)

#No significance


#plot the data
library(ggplot2)
library(ggpubr)
cfu_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/cfu_data/growth_data_mono_rep1.txt", sep = "\t")
cfu_data <- rbind(cfu_data, read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/cfu_data/growth_data_mono_rep2.csv"))
cfu_data$batch <- "One"
cfu_data$batch[40:78] <- "Two"

cfu_data <- cfu_data[cfu_data$class != "control",]
cfu_data$class <- as.character(cfu_data$class)
cfu_data$class[cfu_data$class == "ns"] <- "Non-Supressor"
cfu_data$class[cfu_data$class == "s"] <- "Supressor"
cfu_data$class <- as.factor(cfu_data$class)
cfu_data$class <- factor(cfu_data$class)
cfu_data$taxa <- as.character(cfu_data$taxa)
cfu_data$taxa[cfu_data$taxa == "actino"] <- "Actinobacteria"
cfu_data$taxa[cfu_data$taxa == "alpha"] <- "Alphaproteobacteria"
cfu_data$taxa[cfu_data$taxa == "firmicute"] <- "Firmicute"
cfu_data$taxa[cfu_data$taxa == "gamma"] <- "Gammaproteobacteria"
cfu_data$taxa <- as.factor(cfu_data$taxa)

#Naive box plot
ggplot(cfu_data, aes(class, Log.CFUs.mg.)) +
  geom_boxplot()
cchange <- colnames(cfu_data)
cchange[4] <- "CFU"
colnames(cfu_data) <- cchange
ggplot(data = cfu_data, mapping=aes(x=class, y=CFU))+
  geom_boxplot(outlier.colour = NA)+
  facet_wrap(facets = "taxa", scales="free_x", ncol = 2)+
  #stat_compare_means( comparisons = list(c("ns", "s")), method = "kruskal", label.y = 6)+
  geom_jitter(aes(color=Bacteria, shape=batch), size=2.5 )+
  stat_compare_means(aes(group=class),label.y = 6, method = "kruskal", label.x=1.2)+
  ylab("Log10(CFU)")+
  xlab("Class")
##trial w/o firmucutes
ggplot(data = cfu_data[cfu_data$taxa!="Firmicute",], mapping=aes(x=class, y=CFU))+
  geom_boxplot(outlier.colour = NA)+
  #facet_wrap(facets = "taxa", scales="free_x", ncol = 2)+
  #stat_compare_means( comparisons = list(c("ns", "s")), method = "kruskal", label.y = 6)+
  geom_jitter(aes(color=Bacteria, shape = batch), size=2.5 )+
  stat_compare_means(aes(group=class),label.y = 6, method = "kruskal", label.x=1.2)+
  ylab("Log10(CFU)")+
  xlab("Class")
  
kruskal.test(CFU ~ class, cfu_data[cfu_data$taxa!="Firmicute",])

  
ggboxplot(cfu_data, x="class", y="CFU", facet.by = "taxa")+
  stat_compare_means(label = "p.signif")

#I will now be bringing in the relative change between roots and agar
phosphate_rel <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/old_colo_data/rel_abundance.csv")
meta_35 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_mono/paulo_35_experiment.csv")
phosphate_rel <- na.omit(phosphate_rel)
phosphate_rel$lfc_agar_root <- as.character(as.numeric(phosphate_rel$lfc_agar_root))

data <- c()
for (i in unique(phosphate_rel$Taxa) ) {
  rel_d <- phosphate_rel$lfc_agar_root[phosphate_rel$Taxa==i]
  if ( !grepl("CL", i) ) {
    i <- paste("MF", i, sep="")
  }
  if ( i %in% cfu_data$Bacteria ) {
    data <- append(data, c(i, rel_d, mean(cfu_data$CFU[cfu_data$Bacteria == i]),
                           as.character(unique(cfu_data$class[cfu_data$Bacteria==i])[1]),
                           as.character(unique(cfu_data$taxa[cfu_data$Bacteria==i])[1]),
                           as.character(meta_35$Census.Family[meta_35$Isolate..==i]),
                           as.character(meta_35$Colonization[meta_35$Isolate..==i])))
  }
}

data <- matrix(data, ncol = 7, byrow = T)
data <- as.data.frame(data)
data$V2 <- as.numeric(as.character(data$V2))
data$V3 <- as.numeric(as.character(data$V3))

ggplot(data, aes(V4, V3, color=V6, shape=V7))+
  geom_point(size=2.5)+
  facet_wrap("V5")+
  ylab("log10(CFU)") + 
  xlab("")+
  labs(color="", shape="")

ggplot(data, aes(V3, V2, color=V1, shape=V4))+
  geom_point(size=2)+
  facet_wrap("V5")+
  ggtitle("Agar Vs Root")+
  labs(x="Log10(CFU)", y="logFC", color="", shape="")+
  theme(plot.title = element_text(hjust = .5))

#### ATTEMPT TO PUT THE DATA WITH THE COLONIZATION DATA

#this is for rep2 only
#cfu_data <- na.omit(cfu_data[cfu_data$batch == "Two",])

data <- c()
for (i in unique(phosphate_rel$Taxa) ) {
  rel_d <- phosphate_rel$lfc_agar_root[phosphate_rel$Taxa==i]
  if ( !grepl("CL", i) ) {
    i <- paste("MF", i, sep="")
  }
  if ( i %in% unique(cfu_data$Bacteria) ) {
    mm_data <- cfu_data$CFU[cfu_data$Bacteria == i] 
    for ( j in 1:length(mm_data) ) {
      data <- append(data, c(i, rel_d, mm_data[j],
                             as.character(unique(cfu_data$class[cfu_data$Bacteria==i])[1]),
                             as.character(unique(cfu_data$taxa[cfu_data$Bacteria==i])[1]),
                             as.character(meta_35$Census.Family[meta_35$Isolate..==i]),
                             as.character(meta_35$Colonization[meta_35$Isolate..==i])))
    }
  }
}

data <- matrix(data, ncol = 7, byrow = T)
data <- as.data.frame(data)
data$V2 <- as.numeric(as.character(data$V2))
data$V3 <- as.numeric(as.character(data$V3))
data$combine_col <- sapply(1:nrow(data), function(x){paste(data[x,6], " ",data[x,7], sep="")})
data$combine_col <- factor(data$combine_col)

ggplot(data = data, mapping=aes(x=V1, y=V3))+
  geom_boxplot(outlier.colour = NA)+
  facet_wrap(facets = "V5", scales="free_x", ncol = 2)+
  #stat_compare_means( comparisons = list(c("ns", "s")), method = "kruskal", label.y = 6)+
  geom_jitter(aes(color=combine_col, shape = V4), size=3 )+
  #stat_compare_means(aes(group=V4),label.y = 6, method = "kruskal", label.x=1.2)+
  ylab("Log10(CFU)")+
  xlab("Class")+
  labs(color="Castrillo et al. (2017) Definitions")+
  theme_bw()

#### I will be working with Isai and Omri's Dataset now ####
isai_data <- readRDS("/Users/nicholascolaianni/Documents/dangl_lab/Isai_data/stress_experiment_data/dat_amplicon_useq97_4stresses.RDS")
nrow(isai_data$RelativeAbundance$Tab)
colnames(isai_data$RelativeAbundance$Tab)
isai_data$Rarefied$Map$Tissue
#1000Pi,31C and ph 5.5

#what do we want to do?
# We want to display the rarefied counts for the control conditions

#pull in the sequence to id data 
map_seq_name <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Omri_Data/no_dup_clusters.csv")

#get the sequences in the syncom
links <- c()
for ( i in unique(cfu_data$Bacteria) ) {
  links <- append(links, c(i, as.character(map_seq_name$Useq[grepl(i,map_seq_name$Freezer_Id)])))
}
links <- matrix(links, ncol = 2, byrow = T)
links <- links[-3,]
links[2,2] <- "Sequence_96"
links <- as.data.frame(links)

#I now want to pull out the rows of data that are in our small dataset
rare_data <- isai_data$Rarefied$Tab[row.names(isai_data$Rarefied$Tab) %in% links[,2],]
#now just get the columns of interest
which(isai_data$Rarefied$Map$condition %in% c("1000Pi", "5.5pH", "31C") & isai_data$Rarefied$Map$Tissue == "Root")
rare_data <- rare_data[,which(isai_data$Rarefied$Map$condition %in% c("1000Pi", "5.5pH", "31C") & isai_data$Rarefied$Map$Tissue == "Root")]
rare_data <- as.data.frame(rare_data)
for ( i in row.names(rare_data) ) {
    rare_data$name[row.names(rare_data) == i] <- as.character(links$V1[links$V2 == i])
    rare_data$class[row.names(rare_data) == i] <- as.character(unique(cfu_data$class[cfu_data$Bacteria == rare_data$name[row.names(rare_data) == i]]))
    rare_data$taxa[row.names(rare_data) == i] <- as.character(unique(cfu_data$taxa[cfu_data$Bacteria == rare_data$name[row.names(rare_data) == i]]))
    
}

library(tidyr)
detach("package:ncf", unload = T)

rare_data <- rare_data %>% gather(key="exp", value="rare_counts", -c(class,name,taxa) )

ggplot(data = rare_data, aes(class, rare_counts))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=name))+
  facet_wrap(facets = "taxa", scales = "free_y")

#Duo-association
duo_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/data/cfu_data/duo_association_rep2.csv")

#turn the colonization data to log10
duo_data$CFU.mg_bacteria_1 <- log10(duo_data$CFU.mg_bacteria_1)
duo_data$CFU.ml_bacteria_2 <- log10(duo_data$CFU.ml_bacteria_2)

#Need to keep only the cl18 and MF50 data for plotting
plot_duo <- duo_data[duo_data$Bacteria_1_name %in% c("CL18", "MF50"),]

ggplot(data = plot_duo, aes(Bacteria_2_name, CFU.mg_bacteria_1))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  facet_wrap(facets = "Bacteria_1_name")+
  ylab("Log10(CFU/mg)")+
  theme_bw()

ggplot(data = plot_duo, aes(Bacteria_2_name, root_weight))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  facet_wrap(facets = "Bacteria_1_name")+
  ylab("Root Weight (mg)")+
  theme_bw()
