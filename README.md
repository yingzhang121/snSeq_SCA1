A repository for manuscript: https://www.biorxiv.org/content/10.1101/2021.10.28.466301v2 (Now Accepted by Frontiers in Cellular Neuroscience)

Project Description:


Outline of Data Analysis:

1. Quality Control:  
   We used DIEM ([Enhancing droplet-based single-nucleus RNA-seq resolution using the semi-supervised machine learning classifier DIEM](https://www.nature.com/articles/s41598-020-67513-5)) to clean cells from "Raw_feature_bc_matrix" output from Cellranger.  
   This steps involves manual check of the debries score to determine a meaningful threshold that could separate valid cells from debris.  
   Please note, we used an earlier version of the tool, so we manually set the initial cutoff of debris as cells with less than N reads. The default N=100 never worked in our analysis, and we used N=1000 for smaller samples, and N=3000 for larger samples.  
   
2. Cell Typing:  
   The cleaned count matrix was imported into Seurat for dimension deduction and clustering analysis. At the same time, an independent cell typing was done by SingleR algorithm using DropViz CB meta cell reference.  
   

3. Differential Gene Expression:  
   
