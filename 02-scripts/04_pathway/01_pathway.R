
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
library("pathview")
library(org.Hs.eg.db)

ora =  read_tsv("/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/01_deseq2/hAO_inf_vs_hAO_ctrl/hAO_inf_hAO_ctrl_ORA_KEGG.tsv")
gsea = read_tsv("/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/01_deseq2/hAO_inf_vs_hAO_ctrl/hAO_inf_hAO_ctrl_GSEA_KEGG.tsv")
deseq_IHWallgenes = read_tsv("/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/01_deseq2/hAO_inf_vs_hAO_ctrl/hAO_inf_hAO_ctrl_IHWallGenes.tsv")


##Get the Entrez gene IDs associated with those symbols
EG_IDs = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
##Then get the KEGG IDs associated with those entrez genes.
KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA) #not used in this case 

#Subset the meaningfull columns 
deseq_subset <- deseq_IHWallgenes[,c("gene_name","gene_id","log2FoldChange")]


##Get the Entrez gene IDs associated with those symbols
sym = deseq_subset$gene_name #gene symbols 
EG_IDs = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
#Transform list to dataframe 
EG_IDs_df <- do.call(rbind, lapply(EG_IDs, function(x) data.frame(Value = as.numeric(x))))
# Rename the row names to 'GeneName'
EG_IDs_df$gene_name = rownames(EG_IDs_df)


# Merge the numerical_df with the existing_df based on gene_name
merged_df <- merge(deseq_subset, EG_IDs_df, by = "gene_name")

#Merge only with the mapped genes that have ENTREZ_id therefore remove NA values
merged_nonna= merged_df[!is.na(merged_df$Value), ]

#Select only the required columns 
merged_nonna = merged_nonna %>% select('Value','log2FoldChange')

merged_nonna <- merged_nonna %>% 
  rename("gene_id" = "Value")

#Tuberculosis
hsa05152_Tuberculosis <- pathview(gene.data  = merged_nonna,
                     pathway.id = "hsa05152", 
                     species    = "hsa",
                     low = list(gene = "cyan", cpd = "blue"), 
                     mid = list(gene = "gray", cpd = "gray"), 
                     high = list(gene = "red", cpd ="yellow")
                     #limit      = list(gene=max(abs(merged_nonna)),cpd=1)
                     )

#Leishmaniasis
hsa05140_Leishmaniasis <- pathview(gene.data  = merged_nonna,
                     pathway.id = "hsa05140", 
                     species    = "hsa",
                     #limit      = list(gene=max(abs(merged_nonna)),cpd=1)
)

#IL17-signaling pathway
hsa04657_IL17 <- pathview(gene.data  = merged_nonna,
                     pathway.id = "hsa04657", 
                     species    = "hsa",
                     #limit      = list(gene=max(abs(merged_nonna)),cpd=1)
)

#NOD like receptor signaling pathway
hsa04621_NODLR <- pathview(gene.data  = merged_nonna,
                     pathway.id = "hsa04621", 
                     species    = "hsa",
                     #limit      = list(gene=max(abs(merged_nonna)),cpd=1)
)



