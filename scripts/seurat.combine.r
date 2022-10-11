library(Seurat)

load("sample7643.RData")
load("sample7644.RData")
load("sample7830.RData")
load("sample7831.RData")
load("sample7832.RData")
load("sample7833.RData")

rm(sn7830, sn7831, sn7832, sn7833, sn7643, sn7644)

levels(sce7643$orig.ident) <- "p7643"
levels(sce7644$orig.ident) <- "p7644"
levels(sce7830$orig.ident) <- "p7830"
levels(sce7831$orig.ident) <- "p7831"
levels(sce7832$orig.ident) <- "p7832"
levels(sce7833$orig.ident) <- "p7833"

sn.list <- list(sce7643, sce7644, sce7830, sce7831, sce7832, sce7833)
features <- SelectIntegrationFeatures(object.list = sn.list, nfeatures = 3000)

sn.list <- PrepSCTIntegration(object.list = sn.list, anchor.features = features)

sn.anchors <- FindIntegrationAnchors(object.list = sn.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
# the following code can't run on a pc, using too many memory, < 128G
sn.combined.sct <- IntegrateData(anchorset = sn.anchors, normalization.method = "SCT", dims = 1:30)

sn.combined.sct <- RunPCA(sn.combined.sct, verbose = FALSE)
sn.combined.sct <- RunUMAP(sn.combined.sct, reduction = "pca", dims = 1:30)
save(sn.combined.sct, file="sn.combined.sct.RData")
