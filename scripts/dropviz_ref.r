library(scater)
library(SingleR)
library(SingleCellExperiment)

install.packages('DropSeq.util_2.0.tar.gz', repos=NULL)
library(DropSeq.util)

counts <- readRDS("metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
annot <- readRDS("annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
rownames(annot) <- paste(annot$tissue, annot$subcluster, sep="_")
# colnames(annot) = c("tissue", "class", "class_marker", "type_marker", "full_name", "common_name", "subcluster", "tissue_subcluster")

drop.meta.ref <- SummarizedExperiment(assays=SimpleList(counts = counts), colData=as.data.frame(annot))
drop.meta.ref <- logNormCounts(drop.meta.ref)
saveRDS(drop.meta.ref, "dropviz.metacell.singleR.ref.RDS")
