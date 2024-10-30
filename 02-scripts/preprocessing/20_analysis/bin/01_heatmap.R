#!/usr/bin/env Rscript


'
Usage:
  01_heatmap.R --count_mat=<count_mat> --de_res=<de_res> --prefix=<prefix> [options]

Mandatory arguments:
  --count_mat=<count_mat>                   count mat 
  --de_res=<de_res>                         de_res all genes 
  --prefix=<prefix>                   Prefix for output filenames

' -> doc

#BE AWARE THE COUNT MAT is coming from log1p_cpm_




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
conflicts_prefer(tidygraph::select)

# Input from command line
arguments <- docopt(doc, version = "0.1")
count_mat <- read_csv(arguments$count_mat)
de_res <- read_tsv(arguments$de_res)
prefix <- arguments$prefix
prefix = gsub("_count_mat","",prefix)

# DEBUG

#"cilium organization","cilium movement"
#hAO_inf_hAO_ctrl_IHWsigFCgenes_2_fold.tsv


#count_mat = "/data/projects/2023/Pseudomonas_aeruginosa/30_downstream_analysis/003_batch_effect/countmat/cilium_assembly_count_mat.csv"
#prefix="cilium_assembly"
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

#Metadata & Condition (inf , ctrl) 
split_names <- strsplit(colnames(data), "_")
condition<- sapply(split_names, function(x) x[3])
metadata <-  data.frame(condition)
rownames(metadata) <- colnames(data_num)

# Normalize data 
# Z score function calculation
data_norm <-  t(scale(t(data_num))) # Apply across all rows (1) the Z-score function (Substract mean and divide by sd)
data_norm_gene_name <- data_norm
data_norm_gene_name <- data.frame(gene_name = row.names(data_norm_gene_name), data_norm_gene_name)

#Function for plotting heatmap
function_plotheatmap <- function(heat_mat, metadata,prefix,case,count_mat_merged_de_res) {
  
  if (nrow(heat_mat)<=1){
    print(paste0("For this pathway there is only 1 DEG with padj<0.1. gene_name: ", count_mat_merged_de_res$gene_name))
    return()
  }
  
  #Stablish plot size 
  num_rows <- nrow(heat_mat)
  num_cols <- ncol(heat_mat)
  cell_width <- 40
  cell_height <- 40
  total_heatmap_width <- cell_width * num_cols
  total_heatmap_height <- cell_height * num_rows
  heatmap_width <- total_heatmap_width + 150  # Add extra space for row names and annotations
  heatmap_height <- total_heatmap_height + 150 # Add extra space for column names and annotations
    

  custom_colors <- colorRampPalette(c("#6b85bd",'#7fb4cd', "#eff5d0","#fcd598",'#e7664b'))(100)
  
  color_mapping <- colorRamp2(seq(-2, 2, length.out = length(custom_colors)), custom_colors)
  
  # Rows to highlight
  #To be input from config (pending)
  if (case==2){

  #myRows <- marker_genes
  heat_mat.df <- as.data.frame(heat_mat)
  myRows <-  heat_mat.df$gene_name

  # Set stylings for row names and make our selected rows unique
  row_idx <- which(rownames(heat_mat) %in% myRows)
  fontsizes <- rep(8, nrow(heat_mat))
  fontsizes[row_idx] <- 9
  fontcolors <- rep('black', nrow(heat_mat))
  fontcolors[row_idx] <- 'black'
  fontfaces <- rep('plain',nrow(heat_mat))
  pvalue_col = heat_mat.df$padj
  #row_idx_goi <- which(rownames(heat_mat) %in% myrows_goi)
  fontfaces[row_idx] <- 'bold'
  
  # Create text annotation object for displaying row names
  pvalue_col_fun = colorRamp2(c(0, 1, 2),c("green", "white", "red"))


  count_mat_merged_de_res <-count_mat_merged_de_res %>% mutate_if(is.character,as.numeric)
  rowAnno <- rowAnnotation(
   pvalue = anno_simple(-log10(count_mat_merged_de_res[,"padj"]), col = pvalue_col_fun, pch = "*"),
    gap = unit(2, "mm")
  )
  lgd_pvalue = Legend(title = "p-adj", col_fun = pvalue_col_fun, at = c(0, 1, 2), labels = c("1", "0", "0.1"))
  #lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.1")
  
  # Create our own row dendrogram (ComplexHeatmap orders rows by mean by default)
  dend <- reorder(as.dendrogram(hclust(dist(heat_mat))), -rowMeans(heat_mat), agglo.FUN = mean)
  
  # Find rows in dendrogram
  dend_idx <- which(order.dendrogram(dend) %in% which(rownames(heat_mat) %in% myRows))
  
  
  # Find bottom and top of each row on heatmap (x and y axes go from 0 to 1)
  btm <- 1 - (dend_idx / nrow(heat_mat))
  top <- btm + (1/nrow(heat_mat))
  
  #Stablish plot size 
  num_rows <- nrow(heat_mat)
  num_cols <- ncol(heat_mat)
  cell_width <- 20
  cell_height <- 20
  total_heatmap_width <- cell_width * num_cols
  total_heatmap_height <- cell_height * num_rows
  heatmap_width <- total_heatmap_width + 100  # Add extra space for row names and annotations
  heatmap_height <- total_heatmap_height + 100 # Add extra space for column names and annotations
  
  
  #Color palette for top annotation
  color_Palette <- c("ctrl"="#c7eae5","inf" = "#f6e8c3")
  
  ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color_Palette),labels = unique(rev(condition))))
  
  # Draw the heatmap using our own row clustering and text decorations
  # Capitalize the first word

  prefix_spaces <- gsub("_", " ", prefix)
  prefix_spaces <- paste(toupper(substr(prefix_spaces, 1, 1)), substr(prefix_spaces, 2, nchar(prefix_spaces)), sep = "")
  prefix_spaces = paste(prefix_spaces," (padj < 0.1) ",case)
  #ht <- Heatmap(heat_mat, name = "ht", column_split = metadata$condition, cluster_rows = dend, cluster_columns = FALSE, top_annotation = ha, col = color_mapping,column_title_gp = gpar(fontsize = 10),column_title = prefix_spaces,    heatmap_legend_param = list(title = "z-score"),right_annotation = rowAnno, show_column_dend = FALSE, show_row_dend=FALSE,show_row_names = FALSE, right_annotation = rowAnno)
  hm <- Heatmap(heat_mat, name = "ht", column_split = metadata$condition, 
                cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8), 
                top_annotation = ha, col = color_mapping, column_title_gp = gpar(fontsize = 10),
                column_title = prefix_spaces,   heatmap_legend_param = list(title = "z-score"), 
                show_column_dend = FALSE, show_row_dend=FALSE,
                right_annotation = rowAnno
                ) 
   }else{
    #myRows <- marker_genes
     heat_mat.df <- as.data.frame(heat_mat)
     myRows <-  heat_mat.df$gene_name

     
     # Set stylings for row names and make our selected rows unique
    row_idx <- which(rownames(heat_mat) %in% myRows)
    fontsizes <- rep(8, nrow(heat_mat))
    fontsizes[row_idx] <- 9
    fontcolors <- rep('black', nrow(heat_mat))
    fontcolors[row_idx] <- 'black'
    fontfaces <- rep('plain',nrow(heat_mat))
    #row_idx_goi <- which(rownames(heat_mat) %in% myrows_goi)
    fontfaces[row_idx] <- 'bold'
    color_Palette <- c("inf"="#c7eae5","ctrl" = "#f6e8c3")
    ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color_Palette),labels = unique(rev(condition))))
    dend <- reorder(as.dendrogram(hclust(dist(heat_mat))), -rowMeans(heat_mat), agglo.FUN = mean)
    custom_colors <- colorRampPalette(c("#6b85bd",'#7fb4cd', "#eff5d0","#fcd598",'#e7664b'))(100)
    
    color_mapping <- colorRamp2(seq(-2, 2, length.out = length(custom_colors)), custom_colors)

    prefix_spaces <- gsub("_", " ", prefix)
    prefix_spaces <- paste(toupper(substr(prefix_spaces, 1, 1)), substr(prefix_spaces, 2, nchar(prefix_spaces)), sep = "")
    prefix_spaces = paste(prefix_spaces," (padj < 0.1) ",case)
    # Create text annotation object for displaying row names
    #pvalue_col = top_sorted_count_mat_merged_de_res$padj
    pvalue_col = heat_mat.df$padj
    pvalue_col_fun = colorRamp2(c(0, 1, 2), c("green", "white", "red"))
    count_mat_merged_de_res <-count_mat_merged_de_res %>% mutate_if(is.character,as.numeric)
    rowAnno <- rowAnnotation(
      pvalue = anno_simple(-log10(count_mat_merged_de_res[,"padj"]), col = pvalue_col_fun, pch = "*"),
      gap = unit(2, "mm")
    )
    lgd_pvalue = Legend(title = "p-adj", col_fun = pvalue_col_fun, at = c(0, 1, 2), labels = c("1", "0", "0.1"))
    #lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.1")
    #ht <- Heatmap(heat_mat, name = "ht", column_split = metadata$condition, cluster_rows = dend, cluster_columns = FALSE,row_names_gp = gpar(fontsize = 8), top_annotation = ha, col = color_mapping,column_title_gp = gpar(fontsize = 10),column_title = prefix_spaces,   heatmap_legend_param = list(title = "z-score"), show_column_dend = FALSE, show_row_dend=FALSE, right_annotation = rowAnno)
    hm <- Heatmap(heat_mat, name = "ht", column_split = metadata$condition, 
                  cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8), 
                  top_annotation = ha, col = color_mapping, column_title_gp = gpar(fontsize = 10),
                  column_title = prefix_spaces,   heatmap_legend_param = list(title = "z-score"), 
                  show_column_dend = FALSE, show_row_dend=FALSE,
                  right_annotation = rowAnno
                  )
  }
  
  calc_ht_size = function(hm, unit = "inch") {
    pdf(NULL)
    hm = draw(hm)
    w = ComplexHeatmap:::width(hm)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(hm)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()
    c(w, h)
  }

  save_plot = function(hm, filename, format = "all", res = 300) {
  
  hm_size = calc_ht_size(hm)
  
  if (format == "all" | format == "pdf") {
    pdf(file = paste0(filename, ".pdf"), width = hm_size[1], height = hm_size[2], bg = "white")
    #draw(hm,annotation_legend_list = list(lgd_pvalue, lgd_sig))
    draw(hm,annotation_legend_list = list(lgd_pvalue))
    
    
    
    dev.off()
  }
  
  if (format == "all" | format == "png") {
    png(filename = paste0(filename, ".png"), width = hm_size[1], height = hm_size[2], units = "in", res = res, bg = "white")
    #draw(hm,annotation_legend_list = list(lgd_pvalue, lgd_sig))
    draw(hm,annotation_legend_list = list(lgd_pvalue))
    
    
    dev.off()
  }
  
  if (format == "all" | format == "svg") {
    svg(filename = paste0(filename, ".svg"), width = hm_size[1], height = hm_size[2], bg = "white")
    #draw(hm,annotation_legend_list = list(lgd_pvalue, lgd_sig))
    draw(hm,annotation_legend_list = list(lgd_pvalue))
    
    dev.off()
  }
  }
  


  filename = paste0(prefix, "_heatmap")
  save_plot(hm, filename, format = "svg", res = 300)

}

max_heatmap_rows = 40 # Number of genes per heatmap 
#To be input from config (pending)
goi = c("LCN2", "NGAL","FPN1","TFR1","FTH1","FTL","NOS1","NOS2","NOSIP","NOSTRIN","NOS3","NOS2P2","NOX1","NOX5","HAMP","SLC2A6","SLC25A37","SLCO4A1","SLC16A9","SLC7A5","SLC16A1","SLC35F3","SLC5A8","SLC25A25-AS1","SLC6A14","SLCO3A1","SLC45A3","SLC43A2","SLC25A22","SLC11A2","SLC7A11","SLC8B1","SLC30A7","SLC43A3","SLC35E1","SLC38A2","DUOXA2","DUOX2","NDUFAF6","NDUFV2","IL17C","IL19","IL4I1","IL32","IRAK3","IL1RN","IL4R","IL36G","IRAK2","IL20RA","IRAK1BP1","CXCL8","CXCL6","CXCL3","CXCL5")

# Subset the dataframe count mat with goi
count_mat_subset_goi <- subset(count_mat, gene_name %in% goi)



if (nrow(count_mat_subset_goi)<1){
  print("GOI in Gene set 0 ") # No GOI in gene set, all genes for the heatmap come from de_res sorted by padj 
  count_mat_merged_de_res <- merge(data_norm_gene_name,de_res,by="gene_name") # Get the genes from de_res all genes 
  #sorted_count_mat_merged_de_res <- count_mat_merged_de_res[order(count_mat_merged_de_res$padj,decreasing=FALSE),]
  #top_sorted_count_mat_merged_de_res = head(sorted_count_mat_merged_de_res,max_heatmap_rows) 
  #data_norm <- data_norm[rownames(data_norm) %in% top_sorted_count_mat_merged_de_res$gene_name, ] 
  data_norm <- data_norm[rownames(data_norm) %in% count_mat_merged_de_res$gene_name, ] 
  
  heat_mat <- as.matrix(data_norm) # Matrix from normalized matrix

  case = 1 # indicate that there are no marker genes
  #function_plotheatmap(heat_mat, metadata,prefix,case,top_sorted_count_mat_merged_de_res = top_sorted_count_mat_merged_de_res)
  function_plotheatmap(heat_mat, metadata,prefix,case,count_mat_merged_de_res)
}else {
  num_rows = nrow(count_mat_subset_goi)
  num_genes_de_res = max_heatmap_rows - num_rows 
  #goi_merged_de_res <- merge(count_mat_subset_goi,de_res,by="gene_name")
  
  # create and allocate your new variable to Eval
  
  count_mat_merged_de_res <- merge(data_norm_gene_name,de_res,by="gene_name") # Get the genes from de_res all genes 
  #significant_count_mat_merged_de_res <- subset(count_mat_merged_de_res, padj <0.01)
  #significant_genes = c(significant_count_mat_merged_de_res$gene_name)
  #sorted_sig_count_mat_merged_de_res <- significant_count_mat_merged_de_res[order(significant_count_mat_merged_de_res$padj,decreasing=FALSE),]
  count_mat_merged_de_res <- count_mat_merged_de_res[order(count_mat_merged_de_res$padj,decreasing=FALSE),]
  top_sorted_count_mat_merged_de_res = head(count_mat_merged_de_res,max_heatmap_rows) 
  
  #sorted_count_mat_merged_de_res <- count_mat_merged_de_res[order(count_mat_merged_de_res$padj,decreasing=FALSE),] #ordering by padj 
  #top_sorted_count_mat_merged_de_res = head(sorted_count_mat_merged_de_res,num_genes_de_res) #getting the remaining genes from de_res 
  
  #goi_de_res_concat = rbind(top_sorted_count_mat_merged_de_res,goi_merged_de_res) # the goi are on top and then the rest are ordered by padj 

  #data_norm <- data_norm[rownames(data_norm) %in% significant_genes, ] 
  #data_norm <- data_norm[rownames(data_norm) %in% goi_de_res_concat$gene_name, ] 
  #rownames(data_norm) = paste(rownames(data_norm), "*",round(significant_count_mat_merged_de_res$"padj",digits=5))
  data_norm <- data_norm[rownames(data_norm) %in% count_mat_merged_de_res$gene_name, ] 
  heat_mat <- as.matrix(data_norm) # Matrix from normalized matrix
  case = 2 # case where there are marker genes 
  function_plotheatmap(heat_mat, metadata,prefix,case,count_mat_merged_de_res)
}
  
  




