#Libraries
library(conflicted)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflicts_prefer(edgeR::cpm)
conflicts_prefer(ComplexHeatmap::pheatmap)
library(org.Hs.eg.db)
library(edgeR)
library(decoupleR)
library(OmnipathR)
library(ComplexHeatmap)
library(ggrepel)
library(ggplot2)

#Paths 
count_mat = "/data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/salmon.merged.gene_counts.tsv"
de_res = "/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/04_deseq2_batch_effect/hAO_inf_hAO_ctrl_IHWallGenes.tsv"

# read counts
countdata <- read_tsv(count_mat)


#count_mat <- read_csv(count_mat)
filterValues <- countdata$gene_name
de_res <- read_tsv(de_res)



# remove entries with non-unique gene names:
duplicates <- countdata[duplicated(countdata$gene_name), ]
count.mat <- countdata |>
  filter(! gene_name %in% duplicates$gene_name) |>
  column_to_rownames(var = "gene_name")  |>
  select(- gene_id)

count.mat = subset(count.mat, select = -c(1) )

# get normalized, log-transformed counts:
count.mat.norm <- count.mat |>
  round() |>
  cpm(normalized.lib.sizes = FALSE, log = TRUE)

count.mat.norm.df <- as.data.frame(count.mat.norm)

count.mat.norm.filtered <- count.mat.norm.df %>% dplyr:: select(starts_with("hAO"))



# get TF regulons from dorothea
net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))


# Run wmean

activity <- run_wmean(mat=count.mat.norm.filtered, net=net, .source='source', .target='target', .mor='mor', times=100, minsize=5)

## Visualization
n_tfs <- 30


# Transform to wide matrix
activity_mat <- activity |>
  filter(statistic == 'norm_wmean') |>
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') |>
  column_to_rownames('condition') |>
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- activity %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

activity_mat <- activity_mat[,tfs]


# Scale per sample
activity_mat <- scale(activity_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

activity_mat<- activity_mat[c( "hAO_r1_inf","hAO_r2_inf", "hAO_r3_inf",  "hAO_r4_inf",  "hAO_r1_ctrl",  "hAO_r2_ctrl",  "hAO_r3_ctrl",  "hAO_r4_ctrl" ),] #Reorder rows
df <- data.frame(activity_mat)
target <- c("hAO_r1_inf","hAO_r2_inf", "hAO_r3_inf",  "hAO_r4_inf",  "hAO_r1_ctrl",  "hAO_r2_ctrl",  "hAO_r3_ctrl",  "hAO_r4_ctrl")
df<- df[match(target,rownames(df)),]
df_matrix = data.matrix(df)

rownames(df_matrix) <- as.factor(relevel(c("hAO_r1_inf","hAO_r2_inf", "hAO_r3_inf",  "hAO_r4_inf",  "hAO_r1_ctrl",  "hAO_r2_ctrl",  "hAO_r3_ctrl",  "hAO_r4_ctrl")))

# Plot
pheatmap(df_matrix, border_color = NA, color=my_color, breaks = my_breaks,  show_column_dend = FALSE, show_row_dend=FALSE,
         cluster_rows = F) 

######################################
# Plot
duplicates_de_res <- de_res[duplicated(de_res$gene_name), ]
count.mat.de_res <- de_res |>
  filter(! gene_name %in% duplicates_de_res$gene_name) |>
  column_to_rownames(var = "gene_name")  |>
  select(- gene_id)


deg <-count.mat.de_res %>%
  filter(!is.na(t)) %>% 
  as.matrix()

# removing the first row an first column of the matrix
deg <- deg[, -c(10)]

deg <- as.data.frame(deg)
deg <-deg %>% mutate_if(is.character,as.numeric)


contrast_acts <- run_wmean(mat=deg[, 'padj', drop=FALSE], net=net, .source='source', .target='target',
                           .mor='mor', times = 100, minsize = 5)

# Filter norm_wmean
f_contrast_acts <- contrast_acts %>%
  filter(statistic == 'wmean') %>%
  mutate(rnk = NA)

# Filter top TFs in both signs
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)


ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")


