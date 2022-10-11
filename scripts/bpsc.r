library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(doParallel)
library(BPSC)
registerDoParallel(cores=8)

bpsc_proc <- function(sce) {
    clusters <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters=clusters)
    keep <- rowSums(counts(sce)) > 0
    sce <- sce[keep,]
    a <- calculateCPM(sce)
    message("NORMALIZATION (SCRAN) DONE.")

    groups <- factor(sce$genotype)
    batch <- factor(sce$batch)
    controlIds = which(groups == "WT")
    design <- model.matrix(~groups + batch)
    coef = 2

    message("BEGIN BPSC process.")
    res <- BPglm(data = a, controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
    message("DONE BPSC process.")

    ss <- summary(res)
    sst <- as.data.frame(ss$topTable)
    return(sst)
}

sce <- readRDS("allsample.oligo.Rds")
message("READ IN DATA DONE.")
sst <- bpsc_proc(sce)
write.table(sst, file="bpsc.oligo.sc.result.txt", sep="\t", quote=F)
