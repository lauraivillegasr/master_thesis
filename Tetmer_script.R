library(Tetmer)
P_davidi = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/P_davidi", "P_davidi_k27")
tetmer(P_davidi)

PS1162 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/PS1162", "PS1162_k27")
tetmer(PS1162)
JB051 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/JB051", "JB051_k27")
tetmer(JB051)
DL137G2 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/DL137G2", "DL137G2_k27")
tetmer(DL137G2)

PS1159 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/PS1159", "PS1159_k27")
tetmer(PS1159)

PS1579 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/PS1579", "PS1579_k27")
tetmer(PS1579)

PS1806 = read.spectrum("~/Documents/Master_Thesis/kmerspectra_asex/PS1806", "PS1806_k27")
tetmer(PS1806)

###sex_pops

P_bornheim = read.spectrum("~/Documents/Master_Thesis/kmerspectra_sex/P_bornheim", "P_bornheim_k27")
tetmer(P_bornheim)

P_brombeer = read.spectrum("~/Documents/Master_Thesis/kmerspectra_sex/P_brombeer", "P_brombeer_k27")
tetmer(P_brombeer)




library("findGSE")

findGSE(histo="ES5_kat31", sizek=27, outdir="ES5_27mer")


