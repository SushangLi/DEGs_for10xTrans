# DEGs_for10xTrans
Codes to analyze Differential Expressed Genes (DEGs) between old and young mice brain which RNA molecules are sequenced by 10X Genomics spatial transcriptome.

<i><b>DEG_spatial_DESeq2_main</b></i>: compute the DEGs and return the result tables. Data used are generated by <i>preprocess</i>.

<i><b>DEG_spatial_DESeq2_branch1</b></i>: the branch of <i>DEG_spatial_DESeq2_main</i>, compute the normalized counts of all samples and all sample layers. 

<i><b>preprocess</b></i>: generate count matrix for <i>DEG_spatial_DESeq2</i>; <i><b>ppfun</b></i>: provides functions for <i>preprocess</i>.

