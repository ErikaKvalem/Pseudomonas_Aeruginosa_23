#!/usr/bin/env Rscript


library("DESeq2")
#library("docopt")
#arguments = docopt

# Input and output
#sampleAnnotationCSV <- arguments$sample_sheet
#readCountFile <- arguments$count_table
#cond_col = "group"
#sample_col = "sample"
#covariate_formula = ""
#contrast = c(cond_col, arguments$c1, arguments$c2)

# Testdata
sampleAnnotationCSV = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp_out/deseq2_out/samplesheet.csv"
readCountFile = "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp_out/deseq2_out/counts.tsv"
#result_dir ="/data/projects/2023/atlas_protocol/results/differential_expression"
#cond_col = "group"
#sample_col = "sample"
#contrast = c("group", "LUAD", "LUSC")
#covariate_formula = ""


count_mat <- as.matrix(read.csv(readCountFile,sep="\t", row.names = 1))
#count_mat <- as.data.frame(count_mat)
#count_mat$gene_name <- NULL
#countDataMatrix <- as.matrix(all_counts[ , -1])
allSampleAnno <- read.csv(sampleAnnotationCSV, row.names = 1)

row.names(allSampleAnno) <- colnames(count_mat)

#colnames(count_mat) <- gsub("\\..*", "", colnames(count_mat) ) 
#rownames(count_mat) <- gsub("\\..*", "", rownames(count_mat) ) 
#library(tidyverse)
#allSampleAnno %>% remove_rownames %>% column_to_rownames(var="sample")



all(rownames(allSampleAnno) %in% colnames(count_mat))


################# Start processing
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = allSampleAnno,
                              design = ~ disease+ tumor_stage + group + cell_type + group:cell_type)

# run DESeq

design(dds) <- formula(~ disease+ tumor_stage + group + cell_type + group:cell_type)
dds <- DESeq(dds,useT = TRUE, minmu = 1e-06, minReplicatesForReplace = Inf, parallel = TRUE)
#contrast=c("group", "male", "female")
res <- results(dds)


#Heatmap 
ntd <- normTransform(dds)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("group","cell_type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)

head(coef(dds))

#Heatmap of the sample-to-sample distances
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, vsd$cell_type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Principal component plot of the samples

plotPCA(vsd, intgroup=c("group", "cell_type"))

#volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'pvalue')

plotPCA(res)


library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res), keytype = "ENSEMBL", column = "SYMBOL")

#bcell
celltype1_res <- results(dds, contrast = list(c("group_male_vs_female","groupmale.cell_typeB.cell")))
celltype1_res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(celltype1_res), keytype = "ENSEMBL", column = "SYMBOL")
celltype1_res_df <- as.data.frame(celltype1_res)
celltype1_res_ordered<- celltype1_res_df%>% arrange(desc(padj))
head(celltype1_res_ordered, n =100)

#cd4
celltype2_res <- results(dds, contrast = list(c("group_male_vs_female","groupmale.cell_typeCD4.positive..alpha.beta.T.cell")))
celltype2_res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(celltype1_res), keytype = "ENSEMBL", column = "SYMBOL")
celltype2_res_df <- as.data.frame(celltype2_res)
celltype2_res_ordered<- celltype2_res_df%>% arrange(desc(padj))
head(celltype2_res_ordered, n =100)

#cd8
celltype3_res <- results(dds, contrast = list(c("group_male_vs_female","groupmale.cell_typeCD8.positive..alpha.beta.T.cell")))
celltype3_res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(celltype1_res), keytype = "ENSEMBL", column = "SYMBOL")
celltype3_res_df <- as.data.frame(celltype3_res)
celltype3_res_ordered<- celltype3_res_df%>% arrange(desc(padj))
head(celltype3_res_ordered, n =100)

#Treg
celltype4_res <- results(dds, contrast = list(c("group_male_vs_female","groupmale.cell_typeregulatory.T.cell")))
celltype4_res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(celltype1_res), keytype = "ENSEMBL", column = "SYMBOL")
celltype4_res_df <- as.data.frame(celltype3_res)
celltype4_res_ordered<- celltype4_res_df%>% arrange(desc(padj))
head(celltype4_res_ordered, n =100)

#Malignant
celltype5_res <- results(dds, contrast = list(c("group_male_vs_female","groupmale.cell_typemalignant.cell")))
celltype5_res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(celltype1_res), keytype = "ENSEMBL", column = "SYMBOL")
celltype5_res_df <- as.data.frame(celltype3_res)
celltype5_res_ordered<- celltype5_res_df%>% arrange(desc(padj))
head(celltype4_res_ordered, n =100)

### IHW

# use of IHW for p value adjustment of DESeq2 results
library("IHW")
library("dplyr")
library(tidyverse)
resIHW <- results(dds, filterFun=ihw)

resIHW <- as.data.frame(resIHW ) |>
  rownames_to_column(var = "symbol") |>
  as_tibble() |>
  arrange(padj)



#### write results to TSV and XLSX files
write.csv(res , "/data/projects/2023/LCBiome/nsclc_gender_atlas_tmp_out/deseq2_out/IHWallGenes_00.tsv")
