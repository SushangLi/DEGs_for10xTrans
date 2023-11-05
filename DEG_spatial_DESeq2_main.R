###
# Combine gene count matrixs and row-related gene name file to analyze Differential Expressed Genes (DEGs) with DESeq2.
# Visualize the results by 
# 1. Gene expression heatmap (scaled counts by sample size or z-score in each sample),  
# 2. 
# Return files: 
# DEG.csv, all the DEGs with gene name, scaled mean group counts ($baseMean), log2FoldChange and adjusted p value, based on the treatment group. 
# DEG_up.csv, the subset of DEG.csv which the genes are up-regulated in the treatment group. 
# DEG_dw.csv, the subset of DEG.csv which the genes are down-regulated in the treatment group. 
# DEG_for_enrichment.csv, only contains the gene name and log2FoldChange column of the DEG.csv, is designed for down-stream enrichment analysis codes. 
###
## Set the working directory
setwd("/work/path/total_Slice-MEC/results")

## Import dependencies
library(DESeq2)
library(pheatmap)
library(dplyr)

# Export analysis reports or not to the disk #
write_file <- T


## Prepare the datasets
# gene.csv: contains $gene_name column which refers to the rownames of the count matrix(s) exp1(2,3,...) specificly
# exp1(2,3...): count matrix(s) begain from the 2nd column
gene <- read.csv('../data/MEC_layer2_mean.csv')
exp1 <- read.csv('../data/MEC_layer2_mean.csv')

# generate the integrated dataset
gene_name <- gene$gene_name
exp <- exp1[,2:dim(exp1)[2]]  # remove unused columns
rownames(exp) <- gene_name

# And rearrange the dataset by sample groups
# g1 = exp[,1:4]     # old
# g2 = exp[,5:8]     # young
# g3 = exp[,9:12]    # old
# g4 = exp[,13:16]   # young
# exp.raw = round(cbind.data.frame(g2,g4,g1,g3))
exp_int = round(exp) + 0    # add 1 to avoid 0 exists in every line


## DEG workflow with DESeq2
# DEG parameters for DESeq2 #
# thr.p: p value of the multiple t-test; thr.padj: adjusted p value of the multiple t-test
# thr.exp: the minimum requirement of the adjusted gene expression level ($baseMean) in one sample group, always be 1 * (sample repeats)
# thr.log2fc: the minimum requirement of the absolute log2(Fold Change) of $baseMean in the comparing sets
thr.p <- 0.01
thr.padj <- 0.05
thr.exp <- 8
thr.log2fc <- 0.5

# initialize DESeq2 object with dataset conditions
cond <- factor(c(rep("young",8), rep("old",12)))
dds <- DESeqDataSetFromMatrix(exp_int, data.frame(cond), ~cond)
dds <- DESeq(dds)

# do DEG analysis 
res <- results(dds)
resultsNames(dds)
resLFC <- results(dds, contrast=c("cond", "old", "young"))      # Treatment vs Control

# generate results dataframe
resLFCd <- as.data.frame(resLFC)
resLFCd$gene <- gene_name
norm_counts <- counts(dds, normalized = TRUE)

# filter results by pre-set thresholds
DEG <- filter(resLFCd, padj < thr.padj & abs(log2FoldChange) > thr.log2fc & baseMean > thr.exp)
DEG_up <- filter(DEG, log2FoldChange > 0)
DEG_dw <- filter(DEG, log2FoldChange < 0)

# arrange the filtered results by something like log2FoldChange or padj      ## Arrange the DEGs ##
DEG_up <- arrange(DEG_up, desc(log2FoldChange))
DEG_dw <- arrange(DEG_dw, log2FoldChange)
DEG_top <- function(DEG_up_num, DEG_dw_num) {
    DEG_top <- rbind(DEG_up[1:DEG_up_num,],DEG_dw[1:DEG_dw_num,])
    DEG_top <- na.omit(DEG_top)     # remove blank lines if the selected display number > existed DEG number
    return(DEG_top)
}

# export all the results tables to the working dir
if (write_file == T) {
    write.csv(norm_counts, "normalized_counts.csv")
    write.csv(DEG, "DEG.csv")
    write.csv(DEG_up, "DEG_up.csv")
    write.csv(DEG_dw, "DEG_dw.csv")
    write.csv(data.frame(DEG$gene, DEG$log2FoldChange), "DEG_for_enrichment.csv")
    write.csv(resLFCd, "All_genes_statistics.csv")
}

## Visualizations
# Specify a list of genes of interest (GOIs) for admitted legend(s)
gene_mod <- filter(resLFCd, grepl("some_gene", gene))         ## Add your interested genes here ##

# heatmap settings
heatmap_data <- norm_counts
color = colorRampPalette(c('#436eee', 'white', '#ee0000'))(100)
anno_col = data.frame(Group = cond)
row.names(anno_col) = colnames(heatmap_data)
# heatmap function
pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = anno_col,
         show_colnames = FALSE, fontsize = 8, fontsize_row = 6, color = color,
         border_color = NA)



