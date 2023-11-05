###
# The branch 1 of DEG_spatial_DESeq2_main.R, designed to compute the normalized counts of all samples and all sample layers. 
# Raw counts from the 4 slices which contain 2 brain slices respectively are taken as the arthimetic mean of the 2 brain slices. 
# See normalized_counts.csv in MEC_total folder or run the DEG_spatial_DESeq2_main.R 
# with MEC_total_mean.csv to achieve the normalized counts across all the 20 = 3*4 + (1+1)*4 samples.  
# normalized_counts_integrated_MEC_total.csv: normalized counts across the 16 slices.
# normalized_counts_integrated_MEC_layers.csv: normalized counts across the 16 slices of the 3 layers (layer1, layer23, deep layer).
###
## Set the working directory
setwd("F:/LJY/SANDBOX/codes/DEA/DEG/loupelike/total_Slice-MEC/results")

## Import dependencies
library(DESeq2)
library(dplyr)

# Export analysis reports or not to the disk #
write_file <- T

### for the MEC total
## Prepare the datasets
# gene.csv: contains $gene_name column which refers to the rownames of the count matrix(s) exp1(2,3,...) specificly
# exp1(2,3...): count matrix(s) begain from the 2nd column
gene <- read.csv('../data/MEC_total_mean.csv')
exp1 <- read.csv('../data/MEC_total_mean.csv')

# generate the integrated dataset
gene_name <- gene$gene_name

# remove unused columns
exp1 <- exp1[,2:dim(exp1)[2]]         
  
# calculate the arthimetic mean of the g4- (13:16) and g5- (17:20) columns
arth_mean1 <- (exp1[,13:16] + exp1[,17:20]) / 2

# replace the initial columns and drop unused columns
exp1[,13:16] <- arth_mean1
exp1 <- exp1[,1:16]

# generate the final dataset
exp <- exp1
rownames(exp) <- gene_name
colnames(exp) <- c("Young_WTY2-A", "Young_WTY2-B", "Young_WTY2-C", "Young_WTY2-D",
                  "Young_P14-A", "Young_P14-B", "Young_P14-C", "Young_P14-D", 
                  "Old_LLS20201222A", "Old_LLS20201222B", "Old_LLS20201222C", "Old_LLS20201222D", 
                  "Old_LLSA", "Old_LLSB", "Old_LLSC", "Old_LLSD")

# add something like 1 to avoid 0 exists in every line
exp_int = round(exp) + 0    

## DEG workflow with DESeq2
# initialize DESeq2 object with dataset conditions
cond <- factor(c(rep("young",8), rep("old",8)))
dds <- DESeqDataSetFromMatrix(exp_int, data.frame(cond), ~cond)
dds <- DESeq(dds)

# generate the result dataframe
norm_counts <- counts(dds, normalized = TRUE)

# export all the results tables to the working dir
if (write_file == T) {
    write.csv(norm_counts, "normalized_counts_integrated_MEC_total.csv")
}
### End MEC total

### for all the MEC layers
## Prepare the datasets
# gene.csv: contains $gene_name column which refers to the rownames of the count matrix(s) exp1(2,3,...) specificly
# exp1(2,3...): count matrix(s) begain from the 2nd column
gene <- read.csv('../data/MEC_layer1_mean.csv')

exp1 <- read.csv('../data/MEC_layer1_mean.csv')
exp2 <- read.csv('../data/MEC_layer2_mean.csv')
exp3 <- read.csv('../data/MEC_layer3_mean.csv')
exp4 <- read.csv('../data/MEC_deep_layer_mean.csv')

# generate the integrated dataset
gene_name <- gene$gene_name

# remove unused columns
exp1 <- exp1[,2:dim(exp1)[2]]      
exp2 <- exp2[,2:dim(exp2)[2]]     
exp3 <- exp3[,2:dim(exp3)[2]]     
exp4 <- exp4[,2:dim(exp4)[2]]     
  
# calculate the arthimetic mean of the g4- (13:16) and g5- (17:20) columns
arth_mean1 <- (exp1[,13:16] + exp1[,17:20]) / 2
arth_mean2 <- (exp2[,13:16] + exp2[,17:20]) / 2
arth_mean3 <- (exp3[,13:16] + exp3[,17:20]) / 2
arth_mean4 <- (exp4[,13:16] + exp4[,17:20]) / 2

# replace the initial columns and drop unused columns
exp1[,13:16] <- arth_mean1
exp2[,13:16] <- arth_mean2
exp3[,13:16] <- arth_mean3
exp4[,13:16] <- arth_mean4
exp1 <- exp1[,1:16]
exp2 <- exp2[,1:16]
exp3 <- exp3[,1:16]
exp4 <- exp4[,1:16]

# take the layer2 and layer3 as a whole
exp23 <- exp2 + exp3

# generate the final dataset
exp <- cbind(exp1, exp23, exp4) 
rownames(exp) <- gene_name
#colnames(exp) <- c("Young_WTY2-A", "Young_WTY2-B", "Young_WTY2-C", "Young_WTY2-D",
#                  "Young_P14-A", "Young_P14-B", "Young_P14-C", "Young_P14-D", 
#                  "Old_LLS20201222A", "Old_LLS20201222B", "Old_LLS20201222C", "Old_LLS20201222D", 
#                  "Old_LLSA", "Old_LLSB", "Old_LLSC", "Old_LLSD")

# add something like 1 to avoid 0 exists in every line
exp_int = round(exp) + 1    

## DEG workflow with DESeq2
# initialize DESeq2 object with dataset conditions
cond <- factor(rep(c(rep("young",8), rep("old",8)),3))
dds <- DESeqDataSetFromMatrix(exp_int, data.frame(cond), ~cond)
dds <- DESeq(dds)

# generate the result dataframe
norm_counts <- counts(dds, normalized = TRUE)

# export all the results tables to the working dir
if (write_file == T) {
    write.csv(norm_counts, "normalized_counts_integrated_MEC_layers.csv")
}
### End MEC layers

