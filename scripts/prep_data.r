library(Seurat)
library(SingleCellExperiment)

load("sample7643.RData")
load("sample7644.RData")
load("sample7830.RData")
load("sample7831.RData")
load("sample7832.RData")
load("sample7833.RData")

get_cells <- function(sce, cluster, type){
	idx <- sce$seurat_clusters %in% cluster & sce$singler.celltype %in% type
	return(sce[,idx])
}

bg.7643 <- get_cells(sce7643, c(11), c("Astrocyte.Gja.Slc1a3"))
bg.7644 <- get_cells(sce7644, c(11), c("Astrocyte.Gja.Slc1a3"))
bg.7830 <- get_cells(sce7830, c(11), c("Astrocyte.Gja.Slc1a3"))
bg.7831 <- get_cells(sce7831, c(11), c("Astrocyte.Gja.Slc1a3"))
bg.7832 <- get_cells(sce7832, c(11), c("Astrocyte.Gja.Slc1a3"))
bg.7833 <- get_cells(sce7833, c(11), c("Astrocyte.Gja.Slc1a3"))

cnt <- cbind(bg.7643@assays$RNA@counts, bg.7644@assays$RNA@counts, bg.7830@assays$RNA@counts,
	     bg.7831@assays$RNA@counts, bg.7832@assays$RNA@counts, bg.7833@assays$RNA@counts)
sce.bg <- SingleCellExperiment(assays=list(counts=cnt))
sce.bg$sampleID <- c(rep("sample7643", ncol(bg.7643)), rep("sample7644", ncol(bg.7644)),
			rep("sample7830", ncol(bg.7830)), rep("sample7831", ncol(bg.7831)),
			rep("sample7832", ncol(bg.7832)), rep("sample7833", ncol(bg.7833)))
sce.bg$genotype <- c(rep("WT", ncol(bg.7643)), rep("KI", ncol(bg.7644)),
                        rep("KI", ncol(bg.7830)), rep("KI", ncol(bg.7831)),
                        rep("WT", ncol(bg.7832)), rep("WT", ncol(bg.7833)))
sce.bg$batch <- c(rep(1, ncol(bg.7643)), rep(1, ncol(bg.7644)),
                        rep(2, ncol(bg.7830)), rep(3, ncol(bg.7831)),
                        rep(2, ncol(bg.7832)), rep(3, ncol(bg.7833)))
saveRDS(sce.bg, file="allsample.bg.Rds")

va.7643 <- get_cells(sce7643, c(12,14), c("Gja1.Htra1"))
va.7644 <- get_cells(sce7644, c(12,15), c("Gja1.Htra1"))
va.7830 <- get_cells(sce7830, c(8), c("Gja1.Htra1"))
va.7831 <- get_cells(sce7831, c(12), c("Gja1.Htra1"))
va.7832 <- get_cells(sce7832, c(10), c("Gja1.Htra1"))
va.7833 <- get_cells(sce7833, c(10), c("Gja1.Htra1"))

cnt <- cbind(va.7643@assays$RNA@counts, va.7644@assays$RNA@counts, va.7830@assays$RNA@counts,
	     va.7831@assays$RNA@counts, va.7832@assays$RNA@counts, va.7833@assays$RNA@counts)
sce.va <- SingleCellExperiment(assays=list(counts=cnt))
sce.va$sampleID <- c(rep("sample7643", ncol(va.7643)), rep("sample7644", ncol(va.7644)),
			rep("sample7830", ncol(va.7830)), rep("sample7831", ncol(va.7831)),
			rep("sample7832", ncol(va.7832)), rep("sample7833", ncol(va.7833)))
sce.va$genotype <- c(rep("WT", ncol(va.7643)), rep("KI", ncol(va.7644)),
                        rep("KI", ncol(va.7830)), rep("KI", ncol(va.7831)),
                        rep("WT", ncol(va.7832)), rep("WT", ncol(va.7833)))
sce.va$batch <- c(rep(1, ncol(va.7643)), rep(1, ncol(va.7644)),
                        rep(2, ncol(va.7830)), rep(3, ncol(va.7831)),
                        rep(2, ncol(va.7832)), rep(3, ncol(va.7833)))
saveRDS(sce.va, file="allsample.va.Rds")

oligo.7643 <- get_cells(sce7643, c(9), c("Oligodendroyte.Trf.Il33"))
oligo.7644 <- get_cells(sce7644, c(8), c("Oligodendroyte.Trf.Il33"))
oligo.7830 <- get_cells(sce7830, c(10), c("Oligodendroyte.Trf.Il33"))
oligo.7831 <- get_cells(sce7831, c(10), c("Oligodendroyte.Trf.Il33"))
oligo.7832 <- get_cells(sce7832, c(9), c("Oligodendroyte.Trf.Il33", "Oligodendrocyte.Trf.Klk6"))
oligo.7833 <- get_cells(sce7833, c(9), c("Oligodendroyte.Trf.Il33"))

cnt <- cbind(oligo.7643@assays$RNA@counts, oligo.7644@assays$RNA@counts, oligo.7830@assays$RNA@counts,
	     oligo.7831@assays$RNA@counts, oligo.7832@assays$RNA@counts, oligo.7833@assays$RNA@counts)
sce.oligo <- SingleCellExperiment(assays=list(counts=cnt))
sce.oligo$sampleID <- c(rep("sample7643", ncol(oligo.7643)), rep("sample7644", ncol(oligo.7644)),
			rep("sample7830", ncol(oligo.7830)), rep("sample7831", ncol(oligo.7831)),
			rep("sample7832", ncol(oligo.7832)), rep("sample7833", ncol(oligo.7833)))
sce.oligo$genotype <- c(rep("WT", ncol(oligo.7643)), rep("KI", ncol(oligo.7644)),
                        rep("KI", ncol(oligo.7830)), rep("KI", ncol(oligo.7831)),
                        rep("WT", ncol(oligo.7832)), rep("WT", ncol(oligo.7833)))
sce.oligo$batch <- c(rep(1, ncol(oligo.7643)), rep(1, ncol(oligo.7644)),
                        rep(2, ncol(oligo.7830)), rep(3, ncol(oligo.7831)),
                        rep(2, ncol(oligo.7832)), rep(3, ncol(oligo.7833)))
saveRDS(sce.oligo, file="allsample.oligo.Rds")

granule.7643 <- get_cells(sce7643, c(0:7), c("Neuron.Slc17a7.Gabra6"))
granule.7644 <- get_cells(sce7644, c(0:7), c("Neuron.Slc17a7.Gabra6"))
granule.7830 <- get_cells(sce7830, c(0:6,9), c("Neuron.Slc17a7.Gabra6"))
granule.7831 <- get_cells(sce7831, c(0:8), c("Neuron.Slc17a7.Gabra6"))
granule.7832 <- get_cells(sce7832, c(0:6,8), c("Neuron.Slc17a7.Gabra6"))
granule.7833 <- get_cells(sce7833, c(0:6), c("Neuron.Slc17a7.Gabra6"))

cnt <- cbind(granule.7643@assays$RNA@counts, granule.7644@assays$RNA@counts, granule.7830@assays$RNA@counts,
             granule.7831@assays$RNA@counts, granule.7832@assays$RNA@counts, granule.7833@assays$RNA@counts)
sce.granule <- SingleCellExperiment(assays=list(counts=cnt))
sce.granule$sampleID <- c(rep("sample7643", ncol(granule.7643)), rep("sample7644", ncol(granule.7644)),
                        rep("sample7830", ncol(granule.7830)), rep("sample7831", ncol(granule.7831)),
                        rep("sample7832", ncol(granule.7832)), rep("sample7833", ncol(granule.7833)))
sce.granule$genotype <- c(rep("WT", ncol(granule.7643)), rep("KI", ncol(granule.7644)),
                        rep("KI", ncol(granule.7830)), rep("KI", ncol(granule.7831)),
                        rep("WT", ncol(granule.7832)), rep("WT", ncol(granule.7833)))
sce.granule$batch <- c(rep(1, ncol(granule.7643)), rep(1, ncol(granule.7644)),
                        rep(2, ncol(granule.7830)), rep(3, ncol(granule.7831)),
                        rep(2, ncol(granule.7832)), rep(3, ncol(granule.7833)))
saveRDS(sce.granule, file="allsample.granule.Rds")

purkinje.7643 <- get_cells(sce7643, c(15), c("Neuron.Gad1Gad2.Pcp2"))
purkinje.7644 <- get_cells(sce7644, c(16), c("Neuron.Gad1Gad2.Pcp2"))
purkinje.7830 <- get_cells(sce7830, c(15), c("Neuron.Gad1Gad2.Pcp2"))
purkinje.7831 <- get_cells(sce7831, c(19), c("Neuron.Gad1Gad2.Pcp2"))
purkinje.7832 <- get_cells(sce7832, c(15), c("Neuron.Gad1Gad2.Pcp2"))
purkinje.7833 <- get_cells(sce7833, c(7), c("Neuron.Gad1Gad2.Pcp2"))

cnt <- cbind(purkinje.7643@assays$RNA@counts, purkinje.7644@assays$RNA@counts, purkinje.7830@assays$RNA@counts,
             purkinje.7831@assays$RNA@counts, purkinje.7832@assays$RNA@counts, purkinje.7833@assays$RNA@counts)
sce.purkinje <- SingleCellExperiment(assays=list(counts=cnt))
sce.purkinje$sampleID <- c(rep("sample7643", ncol(purkinje.7643)), rep("sample7644", ncol(purkinje.7644)),
                        rep("sample7830", ncol(purkinje.7830)), rep("sample7831", ncol(purkinje.7831)),
                        rep("sample7832", ncol(purkinje.7832)), rep("sample7833", ncol(purkinje.7833)))
sce.purkinje$genotype <- c(rep("WT", ncol(purkinje.7643)), rep("KI", ncol(purkinje.7644)),
                        rep("KI", ncol(purkinje.7830)), rep("KI", ncol(purkinje.7831)),
                        rep("WT", ncol(purkinje.7832)), rep("WT", ncol(purkinje.7833)))
sce.purkinje$batch <- c(rep(1, ncol(purkinje.7643)), rep(1, ncol(purkinje.7644)),
                        rep(2, ncol(purkinje.7830)), rep(3, ncol(purkinje.7831)),
                        rep(2, ncol(purkinje.7832)), rep(3, ncol(purkinje.7833)))
saveRDS(sce.purkinje, file="allsample.purkinje.Rds")
