
library(data.table)
library(tidyverse)
library(stringr)
library(Seurat)

Bcelldat <- readRDS("./result/Bcelldat.rds")

Idents(Bcelldat) <- "group"
output <- as.data.frame(FindMarkers(Bcelldat, ident.1="COVID19", ident.2="Healthy", test.use="MAST"))

save.image("./result/Bcell_MildvsHC.RData")
write.table(output, "./result/Bcell_MildvsHC.txt", sep="\t", quote=FALSE)
