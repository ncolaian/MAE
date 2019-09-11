#This is a capscale function taken from Isai

#It will be used on RNA seq data

#data points need to have rownames as the sample names and the values for each gene

#it returns a list with everthing you will ned for graphing
perform_cap_scale <- function(data_points, mapping_data, formula, distfun="vegdist", perms=5000, sqrt=T) {
  library(vegan)
  formula <- as.formula(paste("data_points ~ ", formula) )
  
  cap <- capscale(formula, data = mapping_data, dfun=distfun, sqrt.dist = sqrt, distance = "euclidean")
  
  cap_sum <- summary(cap)
  perm_anova_model <- anova.cca(cap, permutations = 5000)
  #the model is in fact significant
  
  perm_anova_terms <- anova.cca(cap, permutations = 5000, by="terms")
  #In this model, bacteria does not seem to be a significant factor
  
  pval_model <- perm_anova_model[1,4]
  cap_sum$constraints
  mapping_data <- cbind( mapping_data, cap_sum$sites)
  #the contrabution of each variable to the variance of the data
  percvar <- round(100* cap$CCA$eig/cap$CCA$tot.chi,2)
  
  #Extract the proportion explained by MDS
  eig <- cap$CA$eig
  mds_var <- round(100*eig/sum(eig),2) #what is the variance explained by each dimension
  
  percvar <- c(percvar, mds_var)
  
  #extract the variance
  chi <- c(cap$tot.chi, cap$CCA$tot.chi, cap$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  
  vartot <- variability_table[2,2]*100
  
  toret=list(Map_cap = mapping_data, cap=cap,cap_sum=cap_sum,perm_anova_model=perm_anova_model,perm_anova_terms=perm_anova_terms,
             pval_model=pval_model,variance_explained_axis=percvar,variability_table=variability_table,
             total_var=vartot)
  return(toret)
}
