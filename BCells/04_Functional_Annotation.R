
library(data.table)
library(tidyverse)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(grid)
library(gridExtra)

CD4_mild <- na.omit(fread("./result/CD4_MildvsHC.txt", sep="\t"))
CD8_mild <- na.omit(fread("./result/CD8_MildvsHC.txt", sep="\t"))
Bcell <- fread("./result/Bcell_MildvsHC.txt", sep="\t")
colnames(Bcell)[1] <- "gene"

a <- inner_join(CD4_mild, CD8_mild, by="gene") %>% inner_join(Bcell, by="gene")

ego_bp <- enrichGO(gene = a$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.05)

ego_bp <- clusterProfiler::simplify(
    ego_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

ego_mf <- enrichGO(gene = a$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.05)

ego_mf <- clusterProfiler::simplify(
    ego_mf,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

ego_cc <- enrichGO(gene = a$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "CC", 
                   pAdjustMethod = "fdr", 
                   qvalueCutoff = 0.05)
ego_cc <- clusterProfiler::simplify(
    ego_cc,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp <- dotplot(ego_bp, showCategory=30, font.size=15, title="GO Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

mf <- dotplot(ego_mf, showCategory=30, font.size=15, title="GO Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

cc <- dotplot(ego_cc, showCategory=30, font.size=15, title="GO Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))



geneid <- bitr(a$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kk <- enrichKEGG(gene         = geneid$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod= "fdr",
                 qvalueCutoff = 0.05)

kegg <- dotplot(kk, showCategory=30, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))


grid.arrange(bp, cc, kegg, layout_matrix=rbind(c(1, 2), c(1, 3)))