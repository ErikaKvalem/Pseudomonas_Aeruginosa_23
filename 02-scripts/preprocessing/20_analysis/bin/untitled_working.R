#!/usr/bin/env Rscript


'
Usage:
  untitled_working.R --count_mat=<count_mat> --de_res=<de_res> --prefix=<prefix> [options]

Mandatory arguments:
  --count_mat=<count_mat>                   count mat 
  --de_res=<de_res>                         de_res all genes 
  --prefix=<prefix>                   Prefix for output filenames

' -> doc

# Libraries
library(docopt)
library(dplyr)
library(tidyr)
library(conflicted)
library(circlize)
library(ggplot2)
library(docopt)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(tibble)
library(readr)
library(stringr)

# Input from command line
arguments <- docopt(doc, version = "0.1")
count_mat <- read_csv(arguments$count_mat)
de_res <- read_tsv(arguments$de_res)
prefix <- arguments$prefix
prefix = gsub("_count_mat","",prefix)


#count_mat = "/data/projects/2023/Pseudomonas_aeruginosa/30_downstream_analysis/003_batch_effect/countmat/response_to_bacterium_count_mat.csv"
#prefix="response_to_bacterium"
#de_res = "/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/04_deseq2_batch_effect/hAO_inf_hAO_ctrl_IHWallGenes.tsv"
#count_mat <- read_csv(count_mat)
#de_res <- read_tsv(de_res)
#prefix = gsub("_count_mat","",prefix)

# Create & reshape data
data <- as.matrix(count_mat)
rownames(data) <- data[,2]
data <- data[,-1:-2]
data<- data[, c( "hAO_r1_inf","hAO_r2_inf", "hAO_r3_inf",  "hAO_r4_inf",  "hAO_r1_ctrl",  "hAO_r2_ctrl",  "hAO_r3_ctrl",  "hAO_r4_ctrl" )] #Reorder columns 
data %>% data.frame() %>% mutate(across(where(is.character), as.numeric)) %>% as.matrix() -> data_num #Convert all characters to number
count_mat <- select(count_mat, -1) # remove the first column of count mat for later step 

#Metadata & Condition (ctrl , inf) 
split_names <- strsplit(colnames(data), "_")
condition<- sapply(split_names, function(x) x[3])
metadata <-  data.frame(condition)
rownames(metadata) <- colnames(data_num)

# Normalize data 
# Z score function calculation
data_norm <-  t(scale(t(data_num))) # Apply across all rows (1) the Z-score function (Substract mean and divide by sd)
data_norm_gene_name <- data_norm
data_norm_gene_name <- data.frame(gene_name = row.names(data_norm_gene_name), data_norm_gene_name)


count_mat_merged_de_res <- merge(data_norm_gene_name,de_res,by="gene_name") # Get the genes from de_res all genes 
significant_count_mat_merged_de_res <- subset(count_mat_merged_de_res, padj <0.01)
significant_genes = c(significant_count_mat_merged_de_res$gene_name)
sorted_count_mat_merged_de_res <- count_mat_merged_de_res[order(count_mat_merged_de_res$padj,decreasing=FALSE),]
top_sorted_count_mat_merged_de_res = head(sorted_count_mat_merged_de_res,40) 
data_norm <- data_norm[rownames(data_norm) %in% top_sorted_count_mat_merged_de_res$gene_name, ] 
heat_mat <- as.matrix(data_norm) # Matrix from normalized matrix

# Plotting
marker_genes <- significant_count_mat_merged_de_res$gene_name
myRows <- marker_genes
row_idx <- which(rownames(heat_mat) %in% myRows)
fontsizes <- rep(8, nrow(heat_mat))
fontsizes[row_idx] <- 9
fontcolors <- rep('black', nrow(heat_mat))
fontcolors[row_idx] <- 'black'
fontfaces <- rep('plain',nrow(heat_mat))
pvalue_col = top_sorted_count_mat_merged_de_res$padj
pvalue_col_fun = colorRamp2(c(0, 1, 2),
                            c("green", "white", "red"))
rowAnno <- rowAnnotation(
  pvalue = anno_simple(-log10(top_sorted_count_mat_merged_de_res$padj), col = pvalue_col_fun, pch = "*"),
  gap = unit(2, "mm")
)
lgd_pvalue = Legend(title = "p-adj", col_fun = pvalue_col_fun, at = c(0, 1, 2), labels = c("0.1", "0.01", "0.001"))
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.001")

color_Palette <- c("inf"="#c7eae5","ctrl" = "#f6e8c3")
ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color_Palette),labels = unique(rev(condition))))
dend <- reorder(as.dendrogram(hclust(dist(heat_mat))), -rowMeans(heat_mat), agglo.FUN = mean)
custom_colors <- colorRampPalette(c("#6b85bd",'#7fb4cd', "#eff5d0","#fcd598",'#e7664b'))(100)

color_mapping <- colorRamp2(seq(-2, 2, length.out = length(custom_colors)), custom_colors)
prefix_spaces <- gsub("_", " ", prefix)
prefix_spaces <- paste(toupper(substr(prefix_spaces, 1, 1)), substr(prefix_spaces, 2, nchar(prefix_spaces)), sep = "")

hm <- Heatmap(heat_mat, name = "ht", column_split = metadata$condition, 
        cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8), 
        top_annotation = ha, col = color_mapping, column_title_gp = gpar(fontsize = 10),
        column_title = prefix_spaces,   heatmap_legend_param = list(title = "z-score"), 
        show_column_dend = FALSE, show_row_dend=FALSE,
        right_annotation = rowAnno)
draw(hm,annotation_legend_list = list(lgd_pvalue, lgd_sig))
