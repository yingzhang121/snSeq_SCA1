source("pipeline.r")
counts <- read_10x("P7644/outs/filtered_feature_bc_matrix")
sn7644 <- diem_to_score(counts, min.count=3000)

# Manual check the result to decide threshold for data cleaning
plot_dbr_score(sn7644, feature="n_genes") + geom_vline(xintercept=0.) + geom_hline(yintercept=00)
sn7644 <- call_targets(sn7644, thresh_score = 0.8, min_genes = 600)
counts <- raw_counts(sn7644)
clean_drops <- get_clean_ids(sn7644)
counts.clean <- counts[,clean_drops]

sce7644 <- standard_seurat(counts.clean)

sce7644 <- singler_proc(sce7644, cb.meta.ref)

save(sce7644, sn7644, file="sample7644.RData")

# Plotting
pdf("plots_sample7644.pdf", width=8, height=6)
plot_dbr_score(sn7644, feature="n_genes") +
        geom_vline(xintercept=0.5) +
        geom_hline(yintercept=2000) +
        ggtitle("Debries Scores: sample 7644, x=0.5, y=2000")
DimPlot(sce7644) + ggtitle("Seurat Clustering: sample 7644")
DimPlot(sce7644, group.by="singler.celltype") + guides(col=guide_legend(ncol=1)) +
        ggtitle("Seurat Clustering, colored by SingleR annotation")
dev.off()

write.table(table(sce7644$singler.celltype, sce7644$seurat_clusters), file="P7644.celltyping.txt", sep="\t", quote=F)
