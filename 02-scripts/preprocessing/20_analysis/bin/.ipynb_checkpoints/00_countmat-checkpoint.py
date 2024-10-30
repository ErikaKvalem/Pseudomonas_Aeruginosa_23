#!/usr/bin/env python3


"""
Usage:
  00_countmat.py --cpm_data=<cpm_data> --gsea_bp=<gsea_bp> --gene_set=<gene_set> [options]

Mandatory arguments:
--cpm_data=<cpm_data>                 cpm adata
--gsea_bp=<gsea_bp>                   gsea_bp
--gene_set=<gene_set>                 gene_set
"""

from docopt import docopt
import pandas as pd

args = docopt(__doc__)

log1p_cpm_data = args["--log1p_cpm_data"]
gsea_bp = args["--gsea_bp"]
gene_set = args["--gene_set"]
gene_set = str(gene_set)
gene_set_name = gene_set.replace(" ","_")


""" /data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/analysis/02_normalized_count_plots/00_convert_raw_counts_to_log2_cpm_plus_one.R 
This is the code used to create the cpm object that is inputed --cpm_data"

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

write.csv(log2_cpm_counts_plus_one_df, "/data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/log2_cpm_counts_plus_one.tsv")
"""

count_mat = pd.read_csv(cpm_data,index_col=0 ) #count per million data
bioprod_gsea = pd.read_csv(gsea_bp , sep = "\t") 

count_mat[['gene_ensemlb','transcript']] = count_mat['gene_id'].str.split('.', expand=True)
count_mat = count_mat[count_mat.columns.drop(list(count_mat.filter(regex="COPD")))] #removing columns from gene expression matrix for the COPD samples. for specific use case

#For the gene ensemlb that are the same
#They are grouped by mean
col_dict  = {}
for i in count_mat.columns[2:-2]:
    col_dict[i] = "mean"


count_mat_mean = count_mat.groupby(count_mat['gene_ensemlb']).aggregate(col_dict)

gene_set = bioprod_gsea[bioprod_gsea["Description"] == gene_set]["core_enrichment"]
gene_set_updated_list = [item.replace('/', ',') for item in gene_set.to_list()]



if len(gene_set_updated_list)>0:
    gene_set_updated_list = gene_set_updated_list[0].split(',')

else:
    print("Empty gene set: "+ str(gene_set))


#print("Length of the gene set" + str(len(gene_set_updated_list)))
#genes_of_interest = ["LCN2", "NGAL","FPN1","TFR1","FTH1","FTL","NOS1","NOS2","NOSIP","NOSTRIN","NOS3","NOS2P2","NOX1","NOX5","HAMP","SLC2A6","SLC25A37","SLCO4A1","SLC16A9","SLC7A5","SLC16A1","SLC35F3","SLC5A8","SLC25A25-AS1","SLC6A14","SLCO3A1","SLC45A3","SLC43A2","SLC25A22","SLC11A2","SLC7A11","SLC8B1","SLC30A7","SLC43A3","SLC35E1","SLC38A2","FTH1","NOS2","DUOXA2","DUOX2","NDUFAF6","NDUFV2","IL17C","IL19","IL4I1","IL32","IRAK3","IL1RN","IL4R","IL36G","IRAK2","IL20RA","IRAK1BP1","CXCL8","CXCL6","CXCL3","CXCL5"]
#common_elements = [element for element in genes_of_interest if element in gene_set_updated_list]
#print("Length of goi in gene set " + str(len(common_elements)))
#print("These are the goi in the gene set"+ str(common_elements))


count_mat_subset_dict={}

gene_set_updated_list = [item.replace('/', ',') for item in gene_set.to_list()]
if len(gene_set_updated_list)>0:
    
    gene_set_updated_list = gene_set_updated_list[0].split(',')
    
    
    count_mat_subset = pd.DataFrame()
    for gene in gene_set_updated_list:
        
        if gene in list(count_mat["gene_name"]):            
            count_mat_subset = pd.concat([count_mat_subset,count_mat[count_mat["gene_name"]== gene]])
            
            del count_mat_subset["gene_ensemlb"]
            del count_mat_subset['transcript']
            del count_mat_subset['gene_id']
   
    count_mat_subset.to_csv(f"{gene_set_name}_count_mat.csv")

else:
    print("Empty gene set: "+ str(gene_set))







