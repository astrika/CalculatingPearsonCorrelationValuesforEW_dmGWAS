# Calculating Phenotypic Pearson Correlation Coefficiencies for EW_dmGWAS
Dealing with duplicate gene probes from DNA microarray analysis and obtaining gene co-expression values by Pearson correlation (preparation for MS GWAS and GEO data analysis by Java version of EW_dmGWAS).

# GSE108000_GeneExprsionMatrix_Normalized_Matched.txt
Loaded raw data from GEO108000 MS brain tissue data set, normalized the values, and matched gene symbols to Agilent DNA microarray probe IDs (from specified platform: GPL13497)
  
# handlingDuplicateProbes.R
Handled any duplicate gene probes by calculating average of the gene expression measures of duplicates

# obtainingPhenotypicPCCvalues.R
i. Obtained the co-expression gene matrices for case/control phenotypes

ii. Reformatted the co-expression gene matrices for indicated EW_dmGWAS parameter - only providing gene interactions from HPRD PPI network
