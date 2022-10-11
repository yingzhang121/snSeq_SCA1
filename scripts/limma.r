library(limma)
library(SingleCellExperiment)
library(scran)
library(scuttle)

limma.proc <- function(sce) {
    clusters <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters=clusters)
    sce <- logNormCounts(sce)
    expr <- logcounts(sce)
    keep <- rowSums(expr) > 0
    expr <- expr[keep,]
    mm <- model.matrix(~ 0 + factor(sce$genotype) + factor(sce$batch))
    colnames(mm) <- c("KI", "WT", "b2", "b3")
    fit <- lmFit(expr, design=mm)
    contr <- makeContrasts(KI-WT, levels = colnames(coef(fit)))
    fit2 <- contrasts.fit(fit, contrasts = contr)
    fit2 <- eBayes(fit2, trend=TRUE)
    res <- topTable(fit2, coef=1, adjust="BH", n=nrow(fit2))
    return(res)
}

sce <- readRDS("allsample.purkinje.Rds")
res <- limma.proc(sce)
write.table(res, file="limma.purkinje.sc.result.txt", sep="\t", quote=F)
