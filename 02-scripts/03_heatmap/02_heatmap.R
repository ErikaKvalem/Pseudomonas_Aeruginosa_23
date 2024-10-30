#!/usr/bin/env Rscript

#Libraries
library(dplyr)
library(tidyr)
library(pheatmap)
# load package
library(dendextend)
library(readr)
library(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(ComplexHeatmap)
ht_opt$TITLE_PADDING = unit(c(40, 40), "points")
# Load data
count_mat_heatmap_path=  "/data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/count_mat_heatmap.csv"
count_mat_heatmap= read_csv(count_mat_heatmap_path)

data <- as.matrix(count_mat_heatmap)
rownames(data) <- data[,2]
data <- data[,-1:-2]
data<- data[, c( "hAO_r1_inf","hAO_r2_inf", "hAO_r3_inf",  "hAO_r4_inf",  "hAO_r1_ctrl",  "hAO_r2_ctrl",  "hAO_r3_ctrl",  "hAO_r4_ctrl" )]

data %>% data.frame() %>% mutate(across(where(is.character), as.numeric)) %>% as.matrix() -> data_num

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
colours <- list(
  condition = c('inf' = '#f1a340', 'ctrl' = '#998ec3' ))


data_norm <- t(apply(data_num, 1, cal_z_score))

heat <- t(scale(t(data_num)))

split_names <- strsplit(colnames(data), "_")
replicates <- sapply(split_names, function(x) x[2])
condition<- sapply(split_names, function(x) x[3])
metadata <-  data.frame(replicates,condition)
rownames(metadata) <- colnames(data_num)

ann <- data.frame(
  condition = metadata$condition,
  stringsAsFactors = FALSE)


colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  #col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(0.5, 'mm'),
  annotation_legend_param = list(
    condition = list(
      nrow = 2, # number of rows across which the legend will be arranged
      title = 'condition',
      title_position = 'topcenter',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize = 10, fontface = 'plain'),
      labels_gp = gpar(fontsize = 10, fontface = 'plain'))))

heat_df <- as.data.frame(heat)
heat_mat <- as.matrix(heat_df)


ha = HeatmapAnnotation(foo = runif(4), bar = sample(c("ctrl", "inf"), 4, replace = TRUE),
                       annotation_legend_param = list(
                         bar = list(
                           title = "Condition",
                           at = c("ctrl", "inf"),
                           labels = c("ctrl", "inf"),labels_gp = gpar(fontsize = 10, fontface = 'plain')
                         )
                       ))

hmap <- Heatmap(heat_mat,
                heatmap_legend_param = list(
                  title = "z-score",
                  at = c(-2, 0, 2),
                  direction="vertical",
                  title_gp = gpar(fontsize = 10, fontface = 'plain')
                ),
                column_title_gp = gpar(fill = c("#f1a340", "#998ec3"),fontsize = 10, fontface = "plain"),
                rect_gp = gpar(col = "white",lwd = 1),
                row_names_gp = gpar(fontsize = 7),
                column_split =   condition,
                cluster_rows = FALSE,
                show_column_dend = FALSE,
                width = unit(7, "cm"), 
                height = unit(20, "cm"),
         
                #top_annotation = colAnn,
                )

draw(hmap)
