library("ggplot2")
library("reshape2") 
library("tidyverse")


###ASEX
matrix_fst_asexMAL <- read.csv("~/Documents/Master_Thesis/asex_matrix_FST_MAL", sep=";")
###The data needs to be in a specific order, otherwise the plot looks weeeeeird
ggplot(data = matrix_fst_asexMAL, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "blue", name="FST")  +labs( x = "MAL", y = "MAL") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + ggtitle("Asexual MAL - PS1159 cov corrected")+  theme(plot.title = element_text(hjust = 0.5))


###SEX
matrix_fst_sexMAL <- read.csv("~/Documents/Master_Thesis/sex_matrix_FST_MAL.csv", sep=";")
###The data needs to be in a specific order, otherwise the plot looks weeeeeird
ggplot(data = matrix_fst_sexMAL, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "orange", name="FST")  +labs( x = "MAL", y = "MAL") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + ggtitle("Sexual MAL - JU765 cov corrected")+  theme(plot.title = element_text(hjust = 0.5))

##elegans
matrix_fst_sexMAL <- read.csv("~/Documents/Master_Thesis/elegans_matrix_FST_MAL.csv", sep=";")
###The data needs to be in a specific order, otherwise the plot looks weeeeeird
ggplot(data = matrix_fst_sexMAL, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "deeppink4", name="FST")  +labs( x = "MAL", y = "MAL") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + ggtitle("Sexual MAL - C. elegans cov corrected")+  theme(plot.title = element_text(hjust = 0.5))



###SEXPOPS
matrix_fst_pops <- read.csv("~/Documents/Master_Thesis/Pop_analysis/sexpop_fst_matrix.csv", sep=";")
###The data needs to be in a specific order, otherwise the plot looks weeeeeird
ggplot(data = matrix_fst_pops, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "paleturquoise4", name="FST")  +labs( x = "Population", y = "Population") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + ggtitle("Sexual panagrolaimus populations")+  theme(plot.title = element_text(hjust = 0.5))

###ASEXPOPS



matrix_fst_asexpops <- read.csv("~/Documents/Master_Thesis/Pop_analysis/asexpop_fst_matrix.csv", sep=";")
###The data needs to be in a specific order, otherwise the plot looks weeeeeird
ggplot(data = matrix_fst_asexpops, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient(low = "white", high = "sienna2", name="FST")  +labs( x = "Population", y = "Population") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),axis.text.y = element_text(size = 12)) + coord_fixed() + ggtitle("Asexual panagrolaimus populations")+  theme(plot.title = element_text(hjust = 0.5))


