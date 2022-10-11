# https://github.com/marcalva/diem/blob/master/vignettes/diem.Rmd
# srun -N 1 -n 1 --x11 -t 3:00:00 --mem=32G --pty bash

library(diem)
library(ggplot2)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)

diem_to_score <- function(counts, min.count) {
  	message("\nStart DIEM Cleaning at: ", Sys.time())
	sce <- create_SCE(counts)
	mt_genes <- grep(pattern = "^mt-", x = rownames(sce@gene_data),
			 ignore.case = TRUE, value = TRUE)
	sce <- get_gene_pct(x = sce, genes = mt_genes, name = "pct.mt")
	genes <- grep(pattern = "^malat1$", x = rownames(sce@gene_data),
		      ignore.case = TRUE, value = TRUE)
	sce <- get_gene_pct(x = sce, genes = genes, name = "MALAT1")
	sce <- set_debris_test_set(sce, min_counts = min.count)
	sce <- filter_genes(sce, cpm_thresh = 0)
	genes <- gene_data(sce)
	sce <- get_pcs(sce, n_var_genes = 50, n_pcs = 10, threads = 8)
	sce <- init(sce, k_init = 20, nstart_init = 30, min_size_init = 2,
			 threads = 8)
	sce <- run_em(sce, threads = 8) # longest step, can use multi-threads
	sce <- assign_clusters(sce)
	sce <- estimate_dbr_score(sce, thresh_genes = 50, thresh_p = 0.5)
	return(sce)
  	message("\nEnd DIEM Cleaning at: ", Sys.time())
}

standard_seurat <- function(counts.clean) {
  message("\nStart Seurat Process at: ", Sys.time())
  sce <- CreateSeuratObject(counts = counts.clean)
  sce <- SCTransform(sce, method = "glmGamPoi")
  sce <- RunPCA(sce)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:30, return.model = TRUE)
  sce <- FindNeighbors(sce, dims = 1:30, verbose = FALSE)
  sce <- FindClusters(sce, verbose = FALSE)
  message("Done Seurat Process at: ", Sys.time())
  return(sce)
}

singler_proc <- function(sce.seurat, ref){
  message("\nStart SingleR Annotation: ", Sys.time())
  se <- as.SingleCellExperiment(sce.seurat)
  common <- intersect(rownames(se), rownames(ref))
  se <- se[common,]
  ref <- ref[common,]
  pred <- SingleR(test = se, ref = ref, labels = ref$full_name,
                  assay.type.ref = "logcounts",
                  BPPARAM = BiocParallel::registered()$MulticoreParam)
  sce.seurat$singler.celltype <- pred$labels
  message("Done SingleR Annotation: ", Sys.time())
  return(sce.seurat)
}

drop.meta.ref <- readRDS("dropviz.metacell.singleR.ref.RDS")
cb.meta.ref <- drop.meta.ref[, colData(drop.meta.ref)$tissue == "CB"]
