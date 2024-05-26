library(Seurat)
setwd("data/test_csv/")

X0 = read.csv2("scRNAseq.csv", sep=",", skip=5)
X0 <-  X0[-c(1854,1855),] # 1-mar und 2-mar erscheint doppelt, was Probleme bereiten wird. Entferne Dublets ohen Werte.
X = (X0[1:21633,-1])
gnames = X0[1:21633,1]

row.names(X)<- gnames

X<- data.frame(lapply(X, as.numeric))

SeuratObj <- CreateSeuratObject(counts = X, project = "Seurat_test_1", min.cells = 0, min.features = 0)

# nFeatures ist einfach die Anzahl an Rows ohne 0-Werte