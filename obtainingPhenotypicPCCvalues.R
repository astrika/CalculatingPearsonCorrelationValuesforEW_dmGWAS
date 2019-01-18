# Astrid M Manuel 01/18/19
# Obtaining Pearson correlation coefficients for genes corresponding to each phenotype
# Run after handlingDuplicateProbes.R

# Creating phenotype gene expression matrices
Control_GEmatrix <- nodupGeneRecords[,1:10]
RIM_ChronicActive_GEmatrix <- nodupGeneRecords[,11:17]
PL_ChronicActive_GEmatrix <- nodupGeneRecords[,18:24]
RIM_Inactive_GEmatrix <- nodupGeneRecords[,25:32]
PL_Inactive_GEmatrix <- nodupGeneRecords[,33:40]

# Aquiring co-expression matrices:
Control_coexpression <- cor(t(Control_GEmatrix))
RIM_ChronicActive_coexpression <- cor(t(RIM_ChronicActive_GEmatrix))
PL_ChronicActive_coexpression <- cor(t(PL_ChronicActive_GEmatrix))
RIM_Inactive_coexpression <- cor(t(RIM_Inactive_GEmatrix))
PL_Inactive_coexpression <- cor(t(PL_Inactive_GEmatrix))

# Re-formatting for java version of EW_dmGWAS
Control_pccValues <- data.frame()
colnames(Control_pccValues) <- c("Gene1","Gene2","PCC")
for(i in seq(1:nrow(Control_coexpression))){
  Control_pccValues[i,1] <- rownames(Control_coexpression)[i]
  for(j in seq(2:ncol(Control_coexpression)+1)){
    Control_pccValues[i,2] <- rownames(Control_coexpression)[j]
    Control_pccValues[i,3] <- Control_coexpression[i,j]
  }
}