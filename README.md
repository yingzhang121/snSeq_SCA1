A repository for manuscript: https://www.biorxiv.org/content/10.1101/2021.10.28.466301v2 (Now Accepted by Frontiers in Cellular Neuroscience)

Project Description:


Outline of Data Analysis:

1. Quality Control:  
   We used DIEM ([Enhancing droplet-based single-nucleus RNA-seq resolution using the semi-supervised machine learning classifier DIEM](https://www.nature.com/articles/s41598-020-67513-5)) to clean cells from "Raw_feature_bc_matrix" output from Cellranger.  
   This steps involves manual check of the debries score to determine a meaningful threshold that could separate valid cells from debris.  
   Please note, we used an earlier version of the tool, so we manually set the initial cutoff of debris as cells with less than N reads. The default N=100 never worked in our analysis, and we used N=1000 for smaller samples, and N=3000 for larger samples.  
   
2. Cell Typing:  
   The cleaned count matrix was imported into Seurat for dimension deduction and clustering analysis. At the same time, an independent cell typing was done by SingleR algorithm using DropViz Cerebellum meta cell reference.  
   So the final cell groups include 1) cells grouped into the same Seurat cluster and 2) cells with the same SingleR annotation tag.
   

3. Differential Gene Expression (DGE):  
   It is still arguable for the best way to run differential gene expression using single-cell/necleus datasets. We decided to use the scran (normalization) - limma_tread (testing for DE) given the following publications:  
   1. [A systematic evaluation of single cell RNA-seq analysis pipelines](https://www.nature.com/articles/s41467-019-12266-7) 
      - [Statistics or biology: the zero-inflation controversy about scRNA-seq data](https://pubmed.ncbi.nlm.nih.gov/35063006/)
   2. [Reproducibility of Methods to Detect Differentially Expressed Genes from Single-Cell RNA Sequencing](https://www.frontiersin.org/articles/10.3389/fgene.2019.01331/full#h4)
   
