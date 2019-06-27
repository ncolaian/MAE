#This script analyzes the syncom data

#This code compares the full syncom cfu experiments

#If you clone the repository, make sure to change this to the path of the data directory
setwd("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/MAE/data/")
cfu_data <- read.csv("cfu_data/full_syncom_rep_2.csv")

library(ggplot2)
library(agricolae)

cfu_data$Syncom <- factor(cfu_data$Syncom, levels = c("Suppressor", "Non-Suppressor", "Full", "NB"))

ggplot(cfu_data, aes(Syncom, CFU))+
  geom_point()

trial <- aov(CFU ~ Syncom, cfu_data)
summary(trial)
posthoc_Tuk <- HSD.test(trial, "Syncom")
posthoc_Tuk


