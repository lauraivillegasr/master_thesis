###box plot of pi diversity#####
JU1646_pi <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_sex_cov_corrected/JU1646.corrected.pi", header=FALSE)
P_bornheim_pi <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_sex_cov_corrected/P_bornheim.corrected.pi", header=FALSE)
P_bromber_pi <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_sex_cov_corrected/P_brombeer.corrected.pi", header=FALSE)

#convert it into a data frame just with our column of interest

JU1646<- data.frame(JU1646_pi$V5)
pborn <- data.frame(P_bornheim_pi$V5)
pbrom <- data.frame(P_bromber_pi$V5)

library(tidyverse)

#specifiy that the values of WT are numbers and not characters
JU1646$pi <- as.numeric(JU1646$JU1646_pi.V5)
pborn$pi <- (pborn$P_bornheim_pi.V5)
pbrom$pi <- (pbrom$P_bromber_pi.V5)

#re-arange the files so they have a common structure and can be merged to make 
#the plot

JU1646 <- JU1646 %>%
  mutate(Strain="JU1646") %>%
  select(Strain, pi)

pborn <- pborn  %>%
  mutate(Strain="P. bornheim") %>%
  select(Strain, pi)

pbrom <- pbrom  %>%
  mutate(Strain="P. brombeer") %>%
  select(Strain, pi)

#merge all files by the matching names - strain and WT
sexpop<- rbind(JU1646, pbrom, pborn)

#export file into a CSV so fields with "na" can easily be discarded 
#--> NEED TO FIND NICER SOLUTION!!!

write.csv(sexpop,"piW1kb_allsex.csv", row.names = FALSE)

pi_sex_informative <-  read.csv("~/Documents/Master_Thesis/Pop_analysis/piW1kb_allsex.csv", sep=";")
#make violin plot including the mean value of pi for each strain.
library("ggplot2")
library("scales")
library("ggthemes")


p1 <- ggplot(pi_sex_informative, aes(y=pi, x=Strain, fill=Strain)) + 
      geom_violin(trim=FALSE)
p1 + coord_flip()
p1 + stat_summary(fun.y=mean, geom="point", size=2, color="gray")
p1 + theme_bw() + stat_summary(fun.y=mean, geom="point", size=2, color="plum4")
##boxplot(pi.all$tajima,ylab="diversity") "general script"

PS1162_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_PS1162.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
PS1579_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_PS1579.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
PS1159_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_PS1159_refpool_merged.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
PS1806_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_PS1806pear.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
DL137G2_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_DL137G2.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
JB051_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_JB051.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)
P_davidi_pi4 <- read.delim("~/Documents/Master_Thesis/Pop_analysis/pi_asex_cov_corrected/pi_p_davidi_reduced.sort.rmd.q30.bam.corrected.pileup.pi", header=FALSE)

Dl137G2 <- data.frame(DL137G2_pi4$V5)
JB051 <- data.frame(JB051_pi4$V5)
PS1159 <- data.frame(PS1159_pi4$V5)
PS1579 <- data.frame(PS1579_pi4$V5)
PS1162 <- data.frame(PS1162_pi4$V5)
PS1806 <- data.frame(PS1806_pi4$V5)
P_davidi <- data.frame(P_davidi_pi4$V5)

Dl137G2$pi <- as.numeric(Dl137G2$DL137G2_pi4.V5)
JB051$pi <- as.numeric(JB051$JB051_pi4.V5)
PS1159$pi <- as.numeric(PS1159$PS1159_pi4.V5)
PS1579$pi <- as.numeric(PS1579$PS1579_pi4.V5)
PS1162$pi <- as.numeric(PS1162$PS1162_pi4.V5)
PS1806$pi <- as.numeric(PS1806$PS1806_pi4.V5)
P_davidi$pi <- as.numeric(P_davidi$P_davidi_pi4.V5)


library(tidyverse)


#re-arange the files so they have a common structure and can be merged to make 
#the plot

Dl137G2 <- Dl137G2 %>%
  mutate(Strain="Dl137G2") %>%
  select(Strain, pi)

JB051 <- JB051  %>%
  mutate(Strain="JB051") %>%
  select(Strain, pi)

PS1159 <- PS1159  %>%
  mutate(Strain="PS1159") %>%
  select(Strain, pi)

PS1579 <- PS1579  %>%
  mutate(Strain="PS1579") %>%
  select(Strain, pi)

PS1162 <- PS1162  %>%
  mutate(Strain="PS1162") %>%
  select(Strain, pi)

PS1806 <- PS1806  %>%
  mutate(Strain="PS1806") %>%
  select(Strain, pi)
 
P_davidi <- P_davidi %>%
  mutate(Strain="P. davidi") %>%
  select(Strain, pi)

piasexpop<- rbind(Dl137G2, JB051, PS1159, PS1579, PS1162, PS1806, P_davidi)
write.csv(piasexpop,"pi_refpool_Asex.csv", row.names = FALSE)
pi_asex_informative <- read.csv("~/Documents/Master_Thesis/Pop_analysis/pi_refpool_Asex.csv", sep=";")


p1 <- ggplot(pi_asex_informative, aes(y=pi, x=Strain, fill=Strain)) + 
  geom_violin(trim=FALSE)
p1 + coord_flip()
p1 + stat_summary(fun.y=mean, geom="point", size=2, color="gray")
p1 + theme_bw() + stat_summary(fun.y=mean, geom="point", size=2, color="plum4")
  
  