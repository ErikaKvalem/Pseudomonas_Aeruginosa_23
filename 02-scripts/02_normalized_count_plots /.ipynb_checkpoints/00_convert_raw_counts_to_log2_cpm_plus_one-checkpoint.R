library(edgeR)
raw_counts = read.csv( "/data/projects/2023/Pseudomonas_aeruginosa/10_rnaseq_pipeline/star_salmon/salmon.merged.gene_counts.tsv",sep = "\t")


only_raw_counts = raw_counts[,-1:-2]

cpm_counts = cpm(only_raw_counts)

log2_cpm_counts_plus_one = log2(cpm_counts+1)

log2_cpm_counts_plus_one_df =as.data.frame(log2_cpm_counts_plus_one)

gene_id = raw_counts$gene_id
gene_name =  raw_counts$gene_name
log2_cpm_counts_plus_one_df <- cbind(gene_name,log2_cpm_counts_plus_one_df)
log2_cpm_counts_plus_one_df <- cbind(gene_id,log2_cpm_counts_plus_one_df)

write.csv(log2_cpm_counts_plus_one_df, "/data/projects/2023/Pseudomonas_aeruginosa/deseq2icbi/rnaseq_pipeline/rnaseq_out/log2_cpm_counts_plus_one.tsv")