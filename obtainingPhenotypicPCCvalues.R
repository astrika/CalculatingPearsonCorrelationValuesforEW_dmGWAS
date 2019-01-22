# Astrid M Manuel 01/18/19
# Obtaining Pearson correlation coefficients for genes corresponding to each phenotype
# Run after handlingDuplicateProbes.R
# *** Because of R workspace size limitations, one phenotype is done at a time. Code below obtains phenotypic PCC values for Control samples.

# Creating phenotype gene expression matrices
Control_GEmatrix <- nodupGeneRecords[,1:10]
#RIM_ChronicActive_GEmatrix <- nodupGeneRecords[,11:17]
#PL_ChronicActive_GEmatrix <- nodupGeneRecords[,18:24]
#RIM_Inactive_GEmatrix <- nodupGeneRecords[,25:32]
#PL_Inactive_GEmatrix <- nodupGeneRecords[,33:40]

# Aquiring co-expression matrices:
Control_coexpression <- cor(t(Control_GEmatrix))
#RIM_ChronicActive_coexpression <- cor(t(RIM_ChronicActive_GEmatrix))
#PL_ChronicActive_coexpression <- cor(t(PL_ChronicActive_GEmatrix))
#RIM_Inactive_coexpression <- cor(t(RIM_Inactive_GEmatrix))
#PL_Inactive_coexpression <- cor(t(PL_Inactive_GEmatrix))

# Re-formatting for java version of EW_dmGWAS, taking into consideration PPI network:
network <- read.table("C:\\Users\\amanuel1\\Dev\\ZhaoLab\\MS\\HPRD_Release9_062910\\network.txt")
Control_PPI_pccValues <- data.frame()

for(i in seq(1:nrow(network))){
  a <- toString(network[i,1])
  b <- toString(network[i,2])
  Control_PPI_pccValues[i,1] <- a
  Control_PPI_pccValues[i,2] <- b
  Control_PPI_pccValues[i,3] <- tryCatch(print(Control_coexpression[a,b]),
                                         error = function(e) print(NA))
}

Control_PPI_pccValues_omitNA <- na.omit(Control_PPI_pccValues)
colnames(Control_PPI_pccValues_omitNA) <- c("Gene1","Gene2","PCC")
write.table(Control_PPI_pccValues_omitNA, "Control_PPI_pccValues_omitNA.txt", quote = F, row.names = F, col.names = F)