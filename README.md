# CalculatingPearsonCorrelationValuesforEW_dmGWAS
Dealing with duplicate gene probes from DNA microarray analysis and obtaining gene co-expression values by Pearson correlation.


Preparing for MS GWAS and GEO data analysis by EW_dmGWAS:

  i.	Loaded raw data from GEO108000 MS brain tissue data set, normalized the values, and matched gene symbols to Agilent DNA microarray probe IDs (from specified platform: GPL13497)
  
  ii. Handled any duplicate gene probes by calculating average of the gene expression measures of duplicates

  iii.	Obtained the co-expression gene matrices for case/control phenotypes

  iv.	Reformatted the co-expression gene matrices for indicated EW_dmGWAS parameter - only providing gene interactions from HPRD PPI network
