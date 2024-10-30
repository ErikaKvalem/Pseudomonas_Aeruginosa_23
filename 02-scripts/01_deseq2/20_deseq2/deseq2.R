#!/usr/bin/env Rscript

#Libraries
library("DESeq2")
library("BiocParallel")
library("org.Hs.eg.db")
library("dplyr")
library("IHW")
library("tibble")
library("readr")
library(data.table)
library("argparser", quietly = TRUE)
library(conflicted)
library(tidyverse)
conflict_prefer("select", "dplyr")
conflicts_prefer(dplyr::filter)
conflicts_prefer(data.table::transpose)
conflicts_prefer(dplyr::rename)
library(EnhancedVolcano)

# Load data
readCountFile = "/data/projects/2023/Pseudomonas_aeruginosa/deseq2icbi/tmp_rnaseq_pipeline/rnaseq_out/salmon.merged.gene_counts.tsv"
sampleAnnotationCSV = "/data/projects/2023/Pseudomonas_aeruginosa/deseq2icbi/tmp_rnaseq_pipeline/tables/samplesheet.csv"
sample_col = "sample"
cond_col = "condition"

fdr_cutoff = 0.1
fc_cutoff = 0.5

# Reading the Annotation sample csv file
sampleAnno <- read_csv(sampleAnnotationCSV)
rownames(sampleAnno) <- sampleAnno$sample

# Reading the Count matrix tsv file
count_mat <- read_tsv(readCountFile)
count_mat <- count_mat |>
  select(c(gene_id, sampleAnno[[sample_col]])) |>
  column_to_rownames("gene_id") |>
  round()

# Check if names are the same
if (!all(rownames(sampleAnno) %in% colnames(count_mat))) {
  print('Row names in sample annotation and column names in count matrix are not the same')
  break
}

# Start processing
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = sampleAnno,
                              design = ~ group + condition + group:condition)

## Keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep, ]

# Save normalized filtered count file
dds <- estimateSizeFactors(dds)

# Run DESeq
dds <- DESeq(dds, parallel = TRUE, minRep=Inf)

# get normalized counts
nc <- counts(dds, normalized=T)
nc <- as.data.frame(nc)
logcounts <- log2(nc + 1)

file_names = c("IHWallGenes_fc"="resIHWsig_fc","IHWallGenes_sig"="resIHWsig","IHWallGenes"="resIHW")


for (interaction_term in resultsNames(dds)[-1]) {
  print(interaction_term)
  
  resIHW <- results(dds, filterFun=ihw, contrast = list(c(interaction_term)),  lfcThreshold=0.5)
  resIHW <- as.data.frame(resIHW) |>
    rownames_to_column(var = "symbol") |>
    as_tibble() |>
    arrange(padj)
  resSHRINK  <- lfcShrink(dds, contrast= list(c(interaction_term)), type ="ashr", lfcThreshold=1) #specifying "normal" because "apeglm" need coef instead of contrast. 
  resSHRINK <- as.data.frame(resSHRINK) |>
    rownames_to_column(var = "symbol") |>
    as_tibble() |>
    arrange(padj) |>
    rename(lfcSE_shrink = lfcSE) |>
    rename(log2FoldChange_shrink = log2FoldChange)
  
  resIHW <- left_join(resIHW, select(resSHRINK, c(symbol,log2FoldChange_shrink,lfcSE_shrink)), by="symbol")
  resIHW$symbol <- gsub("\\..*","",resIHW$symbol)
  resIHW$gene_name <- mapIds(org.Hs.eg.db, keys = resIHW$symbol, keytype = "ENSEMBL", column = "SYMBOL")
  resIHW <- arrange(resIHW, padj)
  
  # Filter for adjusted p-value < fdr_cutoff
  resIHWsig <- resIHW %>% filter(padj < fdr_cutoff)
  # significant genes as DE gene FDR < fdr_cutoff & abs(logfoldchange) > fc_cutoff , all genes as background
  resIHWsig_fc <- resIHWsig %>% filter(abs(log2FoldChange) > fc_cutoff)
  
  for (i in 1:5) {
    
    rows <- c(resIHWsig_fc[i,][1])
    last_name <- c(resIHWsig_fc[i,][11])
    countmat_rowvector <- c(rownames(count_mat))
    gene_name <- grep(rows, countmat_rowvector, value = TRUE)
    data_mod <- logcounts[rownames(logcounts) %in% gene_name, ] 
    t_data_mod  <- transpose(data_mod)
    colnames(t_data_mod) <- rownames(data_mod)
    rownames(t_data_mod) <- colnames(data_mod)
    t_data_mod$mix <- rownames(t_data_mod)
    t_data_mod  <- separate(data = t_data_mod, col = mix, into = c("condition", "number","group"), sep = "\\_")
    t_data_mod$mix <- paste(t_data_mod$group, " - ", t_data_mod$condition)
    
    
    # Check distributions of samples using boxplots
    # Plot the chart.
    p <- ggplot(t_data_mod, aes(x = condition, y = get(gene_name))) + 
      geom_boxplot(aes(fill = group), alpha = .2) +
      geom_line(aes(group = group)) + 
      geom_point(size = 2) + 
      facet_wrap(~ group) 
    title = paste("Gene Expression:", paste0(gene_name,"-",last_name))
    
    p <- p +  labs(y = "Norm log2(CPM+1)", x = "Condition", fill = "Infected/Ctrl" ) +  ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    save_plot <- function(filename, p, width=NULL, height=NULL) {
      if (!is.null(width) && !is.null(height)) {
        ggsave(file.path(paste0(filename, ".png")), plot = p, width = width, height = height)
        ggsave(file.path(paste0(filename, ".svg")), plot = p, width = width, height = height)
      } else {
        ggsave(file.path(paste0(filename, ".png")), plot = p)
        ggsave(file.path(paste0(filename, ".svg")), plot = p)
      }
    }
    
    
    path= path = gsub(" ", "", paste("/data/projects/2023/Pseudomonas_aeruginosa/deseq2icbi/results/",interaction_term,"/figures"))
    #save_plot(file.path(path, paste0(gene_name,"_pairedplot")), p, width = 15, height = 12)
  }
  
  # Stop here if we do not have any DE genes
  if(nrow(resIHWsig_fc) < 1) {
    stop("NO significant DE genes found: check fc_cutoff and fdr_cutoff!")
  }
  
  #write 3 files 
  for (name in names(file_names)){
    path= gsub(" ", "", paste("/data/projects/2023/Pseudomonas_aeruginosa/deseq2icbi/results/",interaction_term,"/",name,".tsv"))
    print(path)
    
    for (res_name in file_names[name]){
     
      #write.csv(get(res_name),path)
      
        
        
      
      }


}
}






