if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenVisR")

library("GenVisR")

mutationData <- read.csv2("~/Documents/Master_Thesis/matrix_tstv.csv")
#no editing to the file needed, it has already been provided as sample-reference-variant

mutationData$sample <- factor(mutationData$sample, levels=unique(mutationData$sample))
                                    
# run TvTi
TvTi(mutationData, fileType="MGI", palette = c("#E69F00", "#56B4E9", "#009E73", 
                                               "#F0E442", "#0072B2", "#CC79A7")) 

