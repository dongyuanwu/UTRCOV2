
library(data.table)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(grid)
library(gridExtra)

mat0 <- readRDS("./CD4network/mat0.rds")

hmat0 <- mat0[[1]]
mmat0 <- mat0[[2]]
smat0 <- mat0[[3]]


min.module.size <- 5
epsilon <- 0.999
healthy_module <- dna::network.modules(hmat0, m=min.module.size, epsilon=epsilon)
mild_module <- dna::network.modules(mmat0, m=min.module.size, epsilon=epsilon)
severe_module <- dna::network.modules(smat0, m=min.module.size, epsilon=epsilon)
F1 <- dna::get.modules(healthy_module)
F2 <- dna::get.modules(mild_module)
F3 <- dna::get.modules(severe_module)

ego_bp <- kklist <- vector(mode="list", length=10)

# Healthy

hmod1 <- names(F1)[F1 == 1]
hmod2 <- names(F1)[F1 == 2]
hmod3 <- names(F1)[F1 == 3]
hmod4 <- names(F1)[F1 == 4]
hmod5 <- names(F1)[F1 == 5]
hmod6 <- names(F1)[F1 == 6]

# healthy module 1
ego_bp[[1]] <- enrichGO(gene = hmod1,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[1]] <- clusterProfiler::simplify(
    ego_bp[[1]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

bp1 <- dotplot(ego_bp[[1]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(hmod1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[1]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg1 <- dotplot(kklist[[1]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# healthy module 2
ego_bp[[2]] <- enrichGO(gene = hmod2,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[2]] <- clusterProfiler::simplify(
    ego_bp[[2]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

bp2 <- dotplot(ego_bp[[2]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(hmod2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[2]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg2 <- dotplot(kklist[[2]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))


# healthy module 3
ego_bp[[3]] <- enrichGO(gene = hmod3,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[3]] <- clusterProfiler::simplify(
    ego_bp[[3]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp3 <- dotplot(ego_bp[[3]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(hmod3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[3]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg3 <- dotplot(kklist[[3]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# healthy module 4
ego_bp[[4]] <- enrichGO(gene = hmod4,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[4]] <- clusterProfiler::simplify(
    ego_bp[[4]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp4 <- dotplot(ego_bp[[4]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(hmod4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[4]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg4 <- dotplot(kklist[[4]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# healthy module 5: no information
# healthy module 6
ego_bp[[5]] <- enrichGO(gene = hmod6,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[5]] <- clusterProfiler::simplify(
    ego_bp[[5]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp5 <- dotplot(ego_bp[[5]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(hmod6, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[5]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg5 <- dotplot(kklist[[5]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))


ego_bp[[1]] <- na.omit(as.data.frame(ego_bp[[1]])[, -7])
ego_bp[[2]] <- na.omit(as.data.frame(ego_bp[[2]])[, -7])
ego_bp[[3]] <- na.omit(as.data.frame(ego_bp[[3]])[, -7])
ego_bp[[4]] <- na.omit(as.data.frame(ego_bp[[4]])[, -7])
ego_bp[[5]] <- na.omit(as.data.frame(ego_bp[[5]])[, -7])
kklist[[1]] <- na.omit(as.data.frame(kklist[[1]])[, -7])
kklist[[2]] <- na.omit(as.data.frame(kklist[[2]])[, -7])
kklist[[3]] <- na.omit(as.data.frame(kklist[[3]])[, -7])
kklist[[4]] <- na.omit(as.data.frame(kklist[[4]])[, -7])
kklist[[5]] <- na.omit(as.data.frame(kklist[[5]])[, -7])

ego_bp[[1]]$module <- kklist[[1]]$module <- "Healthy Module 1"
ego_bp[[2]]$module <- kklist[[2]]$module <- "Healthy Module 2"
ego_bp[[3]]$module <- kklist[[3]]$module <- "Healthy Module 3"
ego_bp[[4]]$module <- kklist[[4]]$module <- "Healthy Module 4"
ego_bp[[5]]$module <- kklist[[5]]$module <- "Healthy Module 6"
ego_bp[[1]]$GeneRatio <- ego_bp[[1]]$Count / length(hmod1)
ego_bp[[2]]$GeneRatio <- ego_bp[[2]]$Count / length(hmod2)
ego_bp[[3]]$GeneRatio <- ego_bp[[3]]$Count / length(hmod3)
ego_bp[[4]]$GeneRatio <- ego_bp[[4]]$Count / length(hmod4)
ego_bp[[5]]$GeneRatio <- ego_bp[[5]]$Count / length(hmod6)
kklist[[1]]$GeneRatio <- kklist[[1]]$Count / length(hmod1)
kklist[[2]]$GeneRatio <- kklist[[2]]$Count / length(hmod2)
kklist[[3]]$GeneRatio <- kklist[[3]]$Count / length(hmod3)
kklist[[4]]$GeneRatio <- kklist[[4]]$Count / length(hmod4)
kklist[[5]]$GeneRatio <- kklist[[5]]$Count / length(hmod6)


# COVID

mmod1 <- names(F2)[F2 == 1]
mmod2 <- names(F2)[F2 == 2]
mmod3 <- names(F2)[F2 == 3]
mmod4 <- names(F2)[F2 == 4]
smod1 <- names(F3)[F3 == 1]

# mild module 1
ego_bp[[6]] <- enrichGO(gene = mmod1,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[6]] <- clusterProfiler::simplify(
    ego_bp[[6]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

bp6 <- dotplot(ego_bp[[6]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(mmod1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[6]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg6 <- dotplot(kklist[[6]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# mild module 2
ego_bp[[7]] <- enrichGO(gene = mmod2,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[7]] <- clusterProfiler::simplify(
    ego_bp[[7]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

bp7 <- dotplot(ego_bp[[7]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(mmod2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[7]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg7 <- dotplot(kklist[[7]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))


# mild module 3
ego_bp[[8]] <- enrichGO(gene = mmod3,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[8]] <- clusterProfiler::simplify(
    ego_bp[[8]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp8 <- dotplot(ego_bp[[8]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(mmod3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[8]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg8 <- dotplot(kklist[[8]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# mild module 4
ego_bp[[9]] <- enrichGO(gene = mmod4,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[9]] <- clusterProfiler::simplify(
    ego_bp[[9]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp9 <- dotplot(ego_bp[[9]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(mmod4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[9]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg9 <- dotplot(kklist[[9]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

# severe module 1
ego_bp[[10]] <- enrichGO(gene = smod1,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05)
ego_bp[[10]] <- clusterProfiler::simplify(
    ego_bp[[10]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)


bp10 <- dotplot(ego_bp[[10]], showCategory=10, font.size=15, title="Biological Process") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))

geneid <- bitr(smod1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kklist[[10]] <- enrichKEGG(gene         = geneid$ENTREZID,
                  organism     = 'hsa',
                  pAdjustMethod= "fdr",
                  qvalueCutoff = 0.05)
kegg5 <- dotplot(kklist[[10]], showCategory=10, font.size=15, title="KEGG Pathway") +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=15))



ego_bp[[6]] <- na.omit(as.data.frame(ego_bp[[6]])[, -7])
ego_bp[[7]] <- na.omit(as.data.frame(ego_bp[[7]])[, -7])
ego_bp[[8]] <- na.omit(as.data.frame(ego_bp[[8]])[, -7])
ego_bp[[9]] <- na.omit(as.data.frame(ego_bp[[9]])[, -7])
ego_bp[[10]] <- na.omit(as.data.frame(ego_bp[[10]])[, -7])
kklist[[6]] <- na.omit(as.data.frame(kklist[[6]])[, -7])
kklist[[7]] <- na.omit(as.data.frame(kklist[[7]])[, -7])
kklist[[8]] <- na.omit(as.data.frame(kklist[[8]])[, -7])
kklist[[9]] <- na.omit(as.data.frame(kklist[[9]])[, -7])
kklist[[10]] <- na.omit(as.data.frame(kklist[[10]])[, -7])


ego_bp[[6]]$module <- kklist[[6]]$module <- "Mild Module 1"
ego_bp[[7]]$module <- kklist[[7]]$module <- "Mild Module 2"
ego_bp[[8]]$module <- kklist[[8]]$module <- "Mild Module 3"
ego_bp[[9]]$module <- kklist[[9]]$module <- "Mild Module 4"
ego_bp[[10]]$module <- kklist[[10]]$module <- "Severe Module 1"
ego_bp[[6]]$GeneRatio <- ego_bp[[6]]$Count / length(mmod1)
ego_bp[[7]]$GeneRatio <- ego_bp[[7]]$Count / length(mmod2)
ego_bp[[8]]$GeneRatio <- ego_bp[[8]]$Count / length(mmod3)
ego_bp[[9]]$GeneRatio <- ego_bp[[9]]$Count / length(mmod4)
ego_bp[[10]]$GeneRatio <- ego_bp[[10]]$Count / length(smod1)
kklist[[6]]$GeneRatio <- kklist[[6]]$Count / length(mmod1)
kklist[[7]]$GeneRatio <- kklist[[7]]$Count / length(mmod2)
kklist[[8]]$GeneRatio <- kklist[[8]]$Count / length(mmod3)
kklist[[9]]$GeneRatio <- kklist[[9]]$Count / length(mmod4)
kklist[[10]]$GeneRatio <- kklist[[10]]$Count / length(smod1)

bp <- rbind(ego_bp[[1]][1:5, ], ego_bp[[2]][1:5, ], ego_bp[[3]][1:5, ], ego_bp[[4]][1:5, ], ego_bp[[5]][1:5, ], 
            ego_bp[[6]][1:5, ], ego_bp[[7]][1:5, ], ego_bp[[8]][1:5, ], ego_bp[[9]][1:5, ], ego_bp[[10]][1:5, ])
kk <- rbind(kklist[[1]][1:5, ], kklist[[2]][1:5, ], kklist[[3]][1:5, ], kklist[[4]][1:5, ], kklist[[5]][1:5, ], 
            kklist[[6]][1:5, ], kklist[[7]][1:5, ], kklist[[8]][1:5, ], kklist[[9]][1:5, ], kklist[[10]][1:5, ])

for(i in 1:10) {
    bp <- rbind(bp, ego_bp[[i]][ego_bp[[i]]$ID %in% bp$ID, ])
    kk <- rbind(kk, kklist[[i]][kklist[[i]]$ID %in% kk$ID, ])
}
bp <- na.omit(bp[!duplicated(bp), ])
kk <- na.omit(kk[!duplicated(kk), ])

bp <- bp[order(bp$GeneRatio, bp$p.adjust, decreasing=c(FALSE, TRUE)), ]
kk <- kk[order(kk$GeneRatio, kk$p.adjust, decreasing=c(FALSE, TRUE)), ]
bp$Description <- factor(bp$Description, level=rev(unique(rev(bp$Description))))
kk$Description <- factor(kk$Description, level=rev(unique(rev(kk$Description))))


# GO BP
ggplot(bp, aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw() +
    scale_colour_gradient(limits=c(0, 0.05), low="red") +
    ylab(NULL) +
    ggtitle("GO Biological Process for CD4+ T Cells") + facet_wrap(.~module, ncol=5) +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))

# KEGG
ggplot(kk, aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw() +
    scale_colour_gradient(limits=c(0, 0.05), low="red") +
    ylab(NULL) +
    ggtitle("KEGG Pathway for CD4+ T Cells") + facet_wrap(.~module, ncol=5) +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=20))