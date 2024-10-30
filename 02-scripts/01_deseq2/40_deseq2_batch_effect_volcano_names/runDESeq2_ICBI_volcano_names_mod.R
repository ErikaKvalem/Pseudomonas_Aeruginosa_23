#!/usr/bin/env Rscript
'runDESeq2_ICBI.R

Usage:
  runDESeq2_ICBI.R <sample_sheet> <count_table> --result_dir=<res_dir> --c1=<c1> --c2=<c2> [options]
  runDESeq2_ICBI.R --help

Arguments:
  <sample_sheet>                CSV file with the sample annotations.
  <count_table>                 TSV file with the read counts

Mandatory options:
  --result_dir=<res_dir>        Output directory
  --c1=<c1>                     Contrast level 1 (perturbation). Needs to be contained in condition_col.
  --c2=<c2>                     Contrast level 2 (baseline). Needs to be contained in condition_col.

Optional options:
  --nfcore                      Indicate that the input samplesheet is from the nf-core RNA-seq ppipeline.
                                Will merge entries from the same sample and infer the sample_id from `group` and `replicate`.
                                If set, this option overrides `sample_col`.
  --condition_col=<cond_col>    Column in sample annotation that contains the condition [default: group]
  --sample_col=<sample_col>     Column in sample annotation that contains the sample names
                                (needs to match the colnames of the count table). [default: sample]
  --paired_grp=<paired_grp>     Column that conatins the name of the paired samples, when dealing with
                                paired data.
  --remove_batch_effect         Indicate that batch effect correction should be applied [default: FALSE]
                                If batch effect correction should be performed, a batch column is needed in the
                                samplesheet (see also --batch_col below)
  --batch_col=<batch_col>       Optional: column in sample annotation that contains the batch [default: batch]
  --covariate_formula=<formula> Formula to model additional covariates (need to be columns in the samplesheet)
                                that will be appended to the formula built from `condition_col`.
                                E.g. `+ age + sex`. Per default, no covariates are modelled.
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --fc_cutoff=<log2 fc cutoff>  Fold change (log2) cutoff for volcano plots [default: 1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
  --gene_id_type=<id_type>      Type of the identifier in the `gene_id` column compatible with AnnotationDbi [default: ENSEMBL]
  --n_cpus=<n_cpus>             Number of cores to use for DESeq2 [default: 1]
  --skip_gsea                   Skip Gene-Set-Enrichment-Analysis step
  --genes_of_interest=<genes>   File (tsv) containing a list of genes to highlight in the volcano plot (column must be named
                                "gene_name").
                                If an optional column named "group" is present, each gene will be associated with the corresponding
                                gene group and a separate volcanon plot for each gene group will be generated (e.g. cytokines).
  --organism=<human|mouse>      Ensebml annotation db [default: human]
  --save_workspace              Save R workspace for this analysis [default: FALSE]
' -> doc

library("conflicted")
library("docopt")

arguments = docopt(doc, version = "0.1")

print(arguments)

library("BiocParallel")
library("DESeq2")
library("IHW")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("clusterProfiler")
library("ReactomePA")
library("writexl")
library("readr")
library("dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
library("EnhancedVolcano")
library("ggpubr")
library("tibble")
library("stringr")
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)
options(stringsAsFactors = FALSE)

#### Get parameters from docopt



# Input and output
sampleAnnotationCSV <- arguments$sample_sheet
readCountFile <- arguments$count_table
results_dir = arguments$result_dir
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
paired_grp <- arguments$paired_grp

# prefix and plot title
prefix <- arguments$prefix
plot_title <- arguments$plot_title

# Sample information and contrasts
nfcore = arguments$nfcore
cond_col = arguments$condition_col
sample_col = arguments$sample_col
contrast = c(cond_col, arguments$c1, arguments$c2)
gene_id_type = arguments$gene_id_type
covariate_formula = arguments$covariate_formula
remove_batch_effect = arguments$remove_batch_effect
batch_col = arguments$batch_col

# Cutoff
fdr_cutoff = as.numeric(arguments$fdr_cutoff)
fc_cutoff = as.numeric(arguments$fc_cutoff)

# GTF for Biotype QC
gtf_file = arguments$gtf_file

# Other
n_cpus = as.numeric(arguments$n_cpus)
skip_gsea = arguments$skip_gsea

genes_of_interest = arguments$genes_of_interest

# set organism (human or mouse)
organism = arguments$organism

# save R workspace
save_ws = arguments$save_workspace

if (organism == "human") {
    anno_db = "org.Hs.eg.db"
    org_kegg = "hsa"
    org_reactome = "human"
    org_wp = "Homo sapiens"
} else if (organism == "mouse") {
    anno_db = "org.Mm.eg.db"
    org_kegg = "mmu"
    org_reactome = "mouse"
    org_wp = "Mus musculus"
} else {
    msg <- paste0("Organism not implemented: ", organism)
    stop(msg)
}
library(anno_db, character.only = TRUE)

# Testdata
# Example1
sampleAnnotationCSV = "/data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/samplesheet_deseq_batch_effect.csv"
readCountFile = "/data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/salmon.merged.gene_counts.tsv"
genes_of_interest = "/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/04_deseq2_batch_effect/volcano_names/genes_of_interest_volcano.csv"
c1="hAO_inf"
c2="hAO_ctrl"
sample_col="sample"
cond_col="group"
batch_col="batch"
remove_batch_effect=TRUE
fdr_cutoff=0.1
fc_cutoff=0.5
contrast = c(cond_col, c1, c2)
gene_id_type = "ENSEMBL"
organism = "human"
paired_grp = NULL
covariate_formula = ""
n_cpus = 1
prefix = "example1"
results_dir = "/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/"
plot_title = NULL
skip_gsea = FALSE

# results_dir = "/data/projects/2021/MicrobialMetabolites/bacterial-supernatant/20_deseq2icbi/debug"
# paired_grp = NULL
# prefix = "example1"
# plot_title = NULL
# nfcore=FALSE
# nfcore=TRUE
# cond_col = "group"
# sample_col = "sample"
# contrast = c("group", "11mix", "10mix")
# gene_id_type = "ENSEMBL"
# covariate_formula = ""
# fdr_cutoff = 0.1
# fc_cutoff = 1
# fc_cutoff = 0.585
# gtf_file = "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf"
# n_cpus = 1
# n_cpus = 8
# skip_gsea = FALSE
# remove_batch_effect=FALSE#
# organism="mouse"

# sampleAnnotationCSV = "testdata/example1/sampleTableN.csv"
# readCountFile = "testdata/example1/merged_gene_counts.txt"
# results_dir = "out"
# paired_grp = NULL
# prefix = "example1"
# plot_title = NULL
# nfcore=FALSE
# nfcore=TRUE
# cond_col = "treatment"
# sample_col = "sample"
# contrast = c("treatment", "PFK158", "DMSO")
# gene_id_type = "ENSEMBL"
# covariate_formula = ""
# fdr_cutoff = 0.1
# fc_cutoff = 1
# fc_cutoff = 0.585
# gtf_file = "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf"
# n_cpus = 1
# n_cpus = 8
# skip_gsea = FALSE

# example_nfcore
# sampleAnnotationCSV = "testdata/example_nfcore/rnaseq_samplesheet.csv"
# readCountFile = "testdata/example_nfcore/salmon.merged.gene_counts.subset.tsv"
# results_dir = "/home/sturm/Downloads/tmp_out"
# results_dir = "./test"
# nfcore = TRUE
# paired_grp = "donor"
# prefix = NULL
# plot_title = NULL
# cond_col = "group"
# sample_col = NULL
# contrast = c("group", "grpA", "grpB")
# fdr_cutoff = 0.1
# fc_cutoff = 1

############### Sanitize parameters and read input data
register(MulticoreParam(workers = n_cpus))

if (is.null(plot_title)) {
  plot_title = paste0(contrast[[2]], " vs. ", contrast[[3]])
}
if (is.null(prefix)) {
  prefix = paste0(contrast[[2]], "_", contrast[[3]])
}


allSampleAnno <- read_csv(sampleAnnotationCSV)

sampleAnno <- allSampleAnno %>%
  filter(get(cond_col) %in% contrast[2:3])

sampleAnno <- data.frame(sampleAnno,stringsAsFactors = FALSE)

# Let's see if we use all samples from samplesheet in DESeq 
sampleSubset = FALSE
if (length(base::setdiff(allSampleAnno[[cond_col]], contrast[2:3])) > 0) {
  sampleSubset = TRUE
}

# Add sample col based on condition and replicate if sample col is not explicitly specified
# and make samplesheet distinct (in case the 'merge replicates' functionality was used).
if(nfcore) {
  sample_col = "sample"
  sampleAnno = sampleAnno %>%
    select(-fastq_1, -fastq_2) %>%
    distinct()
  allSampleAnno = allSampleAnno %>%
    select(-fastq_1, -fastq_2) %>%
    distinct()
}

if (is.null(covariate_formula)) {
  covariate_formula = ""
}
if (remove_batch_effect) {
  if (batch_col %in% names(sampleAnno)) {
    print("Correcting possible batch effects")

    if (! grepl(paste0("+", batch_col), covariate_formula)) {
      covariate_formula = paste0("+", batch_col, covariate_formula)
    }
  } else {
    stop("No batch_col found in sampleSheet, please check")
  }
}
if(is.null(paired_grp)) {
  design_formula <- as.formula(paste0("~", cond_col, covariate_formula))
} else {
  design_formula <- as.formula(paste0("~", paired_grp , " +", cond_col, covariate_formula))
}

count_mat <- read_tsv(readCountFile)
if (gene_id_type == "ENSEMBL") {
  count_mat = count_mat %>% mutate(gene_id= remove_ensg_version(gene_id))
}


ensg_to_genesymbol = count_mat %>% select(gene_id, gene_name)
ensg_to_desc = AnnotationDbi::select(get(anno_db), count_mat$gene_id %>% unique(), keytype = gene_id_type, columns = c("GENENAME")) %>%
  distinct(across(!!gene_id_type), .keep_all = TRUE)

# if we do DESeq on sampleSubset we save also the full count mat for generating a full PCA plot
if (sampleSubset) {
  count_mat_full = count_mat %>%
    select(c(gene_id, allSampleAnno[[sample_col]])) %>%
    column_to_rownames("gene_id") %>%
    round() # salmon does not necessarily contain integers
}

count_mat = count_mat %>%
  select(c(gene_id, sampleAnno[[sample_col]])) %>%
  column_to_rownames("gene_id") %>%
  round() # salmon does not necessarily contain integers

count_mat <- data.frame(count_mat,stringsAsFactors = FALSE)

save_plot <- function(filename, p, width=NULL, height=NULL) {
  if (!is.null(width) && !is.null(height)) {
    ggsave(file.path(paste0(filename, ".png")), plot = p, width = width, height = height)
    ggsave(file.path(paste0(filename, ".svg")), plot = p, width = width, height = height)
  } else {
    ggsave(file.path(paste0(filename, ".png")), plot = p)
    ggsave(file.path(paste0(filename, ".svg")), plot = p)
  }
}

################# Start processing
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = sampleAnno,
                              design = design_formula)

# if we use only a subset of samples for DEseq, make also full dds for a generating a full PCA plot
if (sampleSubset) {
  dds_full <- DESeqDataSetFromMatrix(countData = count_mat_full,
                                colData = allSampleAnno,
                                design = as.formula(paste0("~", cond_col)))

  dds_full <- DESeq(dds_full, parallel = (n_cpus > 1))
}

# count number of detected genes
gene_count <- sapply(
  sampleAnno[[sample_col]], function(s) {
    c <- length(count_mat[[s]][(count_mat[[s]] >10)])
    } 
  ) |>
  enframe() |>
  mutate(group=sampleAnno[[cond_col]][sampleAnno[[sample_col]] == name]) |>
  dplyr::rename(sample=name, genes=value)

p <- ggplot(gene_count, aes(sample, genes, fill=group)) + 
  geom_bar(stat = "identity", color="black") +
  scale_color_brewer(type="qual", palette="Set1") +
  ggtitle("Detected genes") +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")
        )

save_plot(file.path(results_dir, paste0(prefix, "_number_of_detected_genes")), p, width=10, height=7)
write_tsv(gene_count, file.path(results_dir, paste0(prefix, "_number_of_detected_genes.tsv")))

## keep only genes where we have >= 10 reads in total
# keep <- rowSums(counts(dds)) >= 10

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep,]

# save filtered count file
write_tsv(counts(dds) %>% 
            as_tibble(rownames = "gene_id") %>%
            left_join(ensg_to_genesymbol) %>%
            left_join(ensg_to_desc, by = c("gene_id" = gene_id_type) ) %>%
            rename(genes_description = GENENAME),
          file.path(results_dir, paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv")))

# save normalized filtered count file
dds <- estimateSizeFactors(dds)
write_tsv(counts(dds, normalized=TRUE) %>%
            as_tibble(rownames = "gene_id") %>%
            left_join(ensg_to_genesymbol) %>%
            left_join(ensg_to_desc, by = c("gene_id" = gene_id_type) ) %>%
            rename(genes_description = GENENAME),
          file.path(results_dir, paste0(prefix, "_detectedGenesNormalizedCounts_min_10_reads_in_one_condition.tsv")))

# Set the reference to the contrast level 2 (baseline) given by the --c2 option
dds[[cond_col]] = relevel( dds[[cond_col]], contrast[[3]])

# run DESeq
dds <- DESeq(dds, parallel = (n_cpus > 1))

# get normalized counts
nc <- counts(dds, normalized=T)

### IHW

# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun=ihw, contrast=contrast)

resIHW <- as.data.frame(resIHW ) |>
  rownames_to_column(var = "gene_id") |>
  as_tibble() |>
  arrange(padj)

resSHRINK  <- lfcShrink(dds, contrast= contrast, type ="normal") #specifying "normal" because "apeglm" need coef instead of contrast. 
resSHRINK <- as.data.frame(resSHRINK) |>
  rownames_to_column(var = "gene_id") |>
  as_tibble() |>
  arrange(padj) |>
  rename(lfcSE_shrink = lfcSE) |>
  rename(log2FoldChange_shrink = log2FoldChange)
 
resIHW <- left_join(resIHW, select(resSHRINK, c(gene_id,log2FoldChange_shrink,lfcSE_shrink)), by="gene_id")
resIHW  <- resIHW |>
  left_join(ensg_to_genesymbol) |>
  left_join(ensg_to_desc, by = c("gene_id" = gene_id_type) ) |>
  rename(genes_description = GENENAME) |>
  arrange(pvalue)
  
summary(resIHW)
sum(resIHW$padj < fdr_cutoff, na.rm=TRUE)
##################### PCA plots

ddsl <- list(contrast = dds)

# run PCA also for full sample set
if (sampleSubset) {
  ddsl <- append(ddsl, list(full = dds_full))
}

bplapply(names(ddsl), function(dds_name) {

  pca_prefix = ifelse(dds_name == "contrast", prefix, "all_samples")

  vsd <- vst(ddsl[[dds_name]], blind=FALSE)

  if (batch_col %in% names(sampleAnno)) {
    intgroup <- c(cond_col, batch_col)
    shape <- batch_col
    my_aes <- aes(PC1, PC2, color=get(paste0(cond_col,".1")), shape=get(shape))
  } else {
    intgroup <- c(cond_col)
    my_aes <- aes(PC1, PC2, color=get(cond_col))
  }

  pcaData <- plotPCA(vsd, intgroup=intgroup, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  p <- ggplot(pcaData, my_aes) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #geom_text(aes(label=name),vjust=2) +
    scale_color_brewer(type="qual", palette="PiYG") +
    labs(colour= "Group") +
    labs(shape= "Batch") +
    theme_bw() + 
    theme(
      plot.title = element_text(color="black", size=13, face="bold", hjust = 0.5), 
      axis.text=element_text(size=12), legend.text = element_text(size=12))

  save_plot(file.path(results_dir, paste0(pca_prefix, "_PCA")), p, width=7, height=3)

  # PCA plot after removing batch effects
  if (remove_batch_effect) {
    assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd[[batch_col]])

    pcaData <- plotPCA(vsd, intgroup=intgroup, returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    p <- ggplot(pcaData, my_aes) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      #geom_text(vjust = 0,hjust = 0.2, nudge_x = -1, nudge_y = 0.5, aes(label = name)) +
      ggtitle(paste0("PCA: ", plot_title, " after batch effect correction")) +
      scale_color_brewer(type="qual", palette="PRGn") +
      labs(colour="Group" ) +
      labs(shape= "Batch") +
      theme_bw() + 
      theme(
        plot.title = element_text(color="black", size=13, face="bold",hjust = 0.5), 
        axis.text=element_text(size=12), legend.text = element_text(size=12))
  


    save_plot(file.path(results_dir, paste0(pca_prefix, "_PCA_after_batch_effect_correction")), p, width=7, height=3)

  }
})

# Filter for adjusted p-value < fdr_cutoff
resIHWsig <- resIHW %>% filter(padj < fdr_cutoff)

# significant genes as DE gene FDR < fdr_cutoff & abs(logfoldchange) > fc_cutoff , all genes as background
resIHWsig_fc <- resIHWsig %>% filter(abs(log2FoldChange) > fc_cutoff)

# Stop here if we do not have any DE genes
if(nrow(resIHWsig_fc) < 1) {
  stop("NO significant DE genes found: check fc_cutoff and fdr_cutoff!")
}

###### Perform Biotype QC
if(!is.null(gtf_file)) {
  gtf = rtracklayer::import(gtf_file, feature.type="gene") %>%
    as_tibble() %>%
    mutate(gene_id = remove_ensg_version(gene_id))

  count_before = nrow(resIHW)
  resIHW = resIHW %>% left_join(select(gtf, gene_id, gene_type), by=c("gene_id"="gene_id"))
  stopifnot("Number of genes should be the same after adding biotypes"= count_before == nrow(resIHW))
  resIHWsig = resIHWsig %>% left_join(select(gtf, gene_id, gene_type), by=c("gene_id"="gene_id"))

  biotype_counts = resIHWsig %>% group_by(gene_type) %>% count()

  p = biotype_counts %>%
    ggplot(aes(x=gene_type, y=n)) + geom_bar(stat='identity') + theme_bw() + coord_flip()

  save_plot(file.path(results_dir, paste0(prefix, "_biotype_counts")), p)
  write_tsv(biotype_counts, file.path(results_dir, paste0(prefix, "_biotype_counts.tsv")))

}

#### result list
de_res_list <- list(IHWallGenes = resIHW, IHWsigGenes = resIHWsig, IHWsigFCgenes = resIHWsig_fc)

#### write results to TSV and XLSX files
lapply(names(de_res_list), function(res) {
  fc_suffix <- ifelse(res == "IHWsigFCgenes", paste0("_", 2^fc_cutoff, "_fold"), "")
  write_tsv(de_res_list[[res]], file.path(results_dir, paste0(prefix, "_", res, fc_suffix, ".tsv")))
  write_xlsx(de_res_list[[res]], file.path(results_dir, paste0(prefix, "_" , res, fc_suffix, ".xlsx")))
})

###### Run TOPGO analysis
de_symbols <- resIHWsig_fc$gene_id
bg_symbols <- rownames(dds)[rowSums(counts(dds)) > 0]

lapply(c("BP", "MF"), function(ontology) {
  topgoDE <- topGOtable(de_symbols, bg_symbols,
                        ontology = ontology,
                        mapping = anno_db,
                        geneID = gene_id_type)
  write_tsv(topgoDE, file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".tsv")))
  write_xlsx(topgoDE %>% select(-genes), file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".xlsx")))
})


##### Pathway enrichment analysis
hgnc_to_entrez = AnnotationDbi::select(get(anno_db), resIHW %>% pull("gene_name") %>% unique(), keytype="SYMBOL", columns=c("ENTREZID"))

# full list with ENTREZIDs added
resIHW_entrez = resIHW %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
universe = resIHW_entrez %>% pull("ENTREZID") %>% unique()

# list of significant genes with ENTREZIDs added
resIHWsig_fc_entrez <- resIHWsig_fc %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
de_foldchanges <- resIHWsig_fc_entrez$log2FoldChange
names(de_foldchanges) <- resIHWsig_fc_entrez$ENTREZID

## ORA
ora_tests = list(
  "KEGG" = function(genes, universe) {
    enrichKEGG(
      gene         = genes,
      universe     = universe,
      organism     = org_kegg,
      pvalueCutoff = 0.05
    )
  },
  "Reactome" = function(genes, universe) {
    enrichPathway(
      gene = genes,
      organism = org_reactome,
      universe = universe,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  },
  "WikiPathway" = function(genes, universe) {
    enrichWP(
      gene = genes,
      universe     = universe,
      organism     = org_wp,
      pvalueCutoff = 0.05
    )
  },
  "GO_BP" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_MF" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "MF",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  }
)

# Warmup GO database - work around https://github.com/YuLab-SMU/clusterProfiler/issues/207
._ = enrichGO(universe[1], OrgDb = get(anno_db), keyType = "ENTREZID", ont = "BP", universe = universe)

get_heatplot_dims <- function(p) {
  nr_gene <- length(unique(p$data$Gene))
  nr_cat <- length(unique(p$data$categoryID))

  hp_width = min(nr_gene * 0.25, 40)
  hp_height = min(nr_cat * 0.25, 40)

  return(c(hp_width, hp_height))
}

bplapply(names(ora_tests), function(ora_name) {
  message(paste0("Performing ", ora_name, " ORA-test..."))

  test_fun = ora_tests[[ora_name]]
  ora_res = test_fun(resIHWsig_fc_entrez$ENTREZID, universe)

  if (!is.null(ora_res)) {
    ora_res = setReadable(ora_res, OrgDb = anno_db, keyType="ENTREZID")
    res_tab = as_tibble(ora_res@result)
    write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_ORA_", ora_name, ".tsv")))

    if (min(res_tab$p.adjust) < 0.05 & length(unique(res_tab$geneID)) > 1) {
      #p = dotplot(ora_res, showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=25)) +
      #ggtitle("ORA GO BP") + 
      #theme(plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5)) 
      
      p = barplot(ora_res, showCategory=10, font.size = 12) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=25)) + 
        ggtitle("ORA GO BP") + 
        theme(plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5)) +
        scale_colour_gradient2()
      
      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_barplot")), p, width = 6, height = 11)
      
      #save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_dotplot")), p, width = 6, height = 11)

      #p <- cnetplot(ora_res,
      #             categorySize="pvalue",
      #              showCategory = min(5, length(ora_res@result$ID)),
      #              foldChange=de_foldchanges,
      #              vertex.label.font=6)

      #save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_cnetplot")), p, width = 15, height = 12)

      #p <- heatplot(ora_res, foldChange=de_foldchanges, showCategory=40) +
      #  scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
      #hp_dims <- get_heatplot_dims(p)

      #save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_heatplot")), p, width = hp_dims[1], height = hp_dims[2])

    } else {
      message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
    }
  } else {
    message(paste0("Warning: No gene can be mapped in ", ora_name, " ORA analysis. "))
  }
})


## GSEA
if(!skip_gsea) {
  # for GSEA use genes ranked by test statistic
  res_ihw_ranked = resIHW_entrez %>%
    arrange(-stat) %>%
    select(ENTREZID, stat) %>%
    na.omit() %>%
    distinct(ENTREZID, .keep_all=TRUE)
  ranked_gene_list = res_ihw_ranked$stat
  names(ranked_gene_list) = res_ihw_ranked$ENTREZID

  gsea_tests = list(
    "KEGG"=function(ranked_gene_list) {
      gseKEGG(geneList = ranked_gene_list, organism = org_kegg, pvalueCutoff = 1)
    },
    "Reactome"=function(ranked_gene_list) {
      gsePathway(geneList = ranked_gene_list, organism = org_reactome, pvalueCutoff = 1)
    },
    "WikiPathway"=function(ranked_gene_list) {
      gseWP(geneList = ranked_gene_list, organism = org_wp, pvalueCutoff = 1)
    },
    "GO_BP"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
               keyType = "ENTREZID",
               OrgDb = anno_db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               minGSSize = 10)
    },
    "GO_MF"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
               keyType = "ENTREZID",
               OrgDb = anno_db,
               ont = "MF",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               minGSSize = 10)
    }
  )

  bplapply(names(gsea_tests), function(gsea_name) {
    message(paste0("Performing ", gsea_name, " GSEA-test..."))

    test_fun = gsea_tests[[gsea_name]]
    gsea_res = test_fun(ranked_gene_list)

    if (!is.null(gsea_res)) {
      gsea_res = setReadable(gsea_res, OrgDb = get(anno_db), keyType="ENTREZID")
      res_tab = gsea_res@result %>% as_tibble()
      res_mat = as.matrix(gsea_res@result)
      gene_set_oi =c("motile cilium assembly", 
                     "intraciliary transport",
                     "cilium organization", 
                     "cilium assembly",
                     "cilium movement" ,      
                     "non-motile cilium assembly" ,     
                     "epithelial cilium movement involved in extracellular fluid movement",       
                     "cilium movement involved in cell motility",
                     "regulation of cilium movement",
                     "epithelial cilium movement involved in determination of left/right asymmetry",
                     "cilium or flagellum-dependent cell motility",
                     "cilium-dependent cell motility")

      write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, ".tsv")))
      if (min(res_tab$p.adjust) < 0.05) {
        
      
  
        
        
        p = dotplot(gsea_res, showCategory=gene_set_oi) +
          scale_y_discrete(labels=function(x) str_wrap(x, width=25)) + 
          ggtitle("GSEA GO BP - Supressed Gene sets") + 
          theme(plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5)) 
        
        
  
        save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_dotplot")), p, width = 6, height = 11)


        #p <- cnetplot(gsea_res,
        #              categorySize="pvalue",
        #              showCategory = 5,
        #              foldChange=de_foldchanges,
        #              vertex.label.font=6)

        #save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_cnetplot")), p, width = 15, height = 12)

        # GSEA generates to long gene lists so that the heatplot gets to overloaded
        # p <- heatplot(gsea_res, foldChange=de_foldchanges, showCategory=40) +
        #   scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
        #
        # hp_dims <- get_heatplot_dims(p)
        #
        # ggsave(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_heatplot.png")), plot = p, width = hp_dims[1], height = hp_dims[2])
      } else {
        message(paste0("Warning: No significant enrichment in ", gsea_name, " GSEA analysis. "))
      }
    } else {
      message(paste0("Warning: No gene can be mapped in ", gsea_name, " GSEA analysis. "))
    }
  })
}



########## Volcano plot
p <- EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 1e-6,
                FCcutoff = fc_cutoff,
                subtitle = "",
                legendPosition = "right",
                caption = paste0("fold change cutoff: ", round(2**fc_cutoff, 1), ", p-value cutoff: ", 1e-6),
                title = plot_title)

save_plot(file.path(results_dir, paste0(prefix, "_volcano")), p, width = 9, height = 7)


p <- EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = fdr_cutoff,
                FCcutoff = fc_cutoff,
                subtitle = "",
                legendPosition = "top",
                caption = paste0("fold change cutoff: ", round(2**fc_cutoff, 1), ", adj.p-value cutoff: ", fdr_cutoff),
                title = plot_title)

save_plot(file.path(results_dir, paste0(prefix, "_volcano_padj")), p, width = 9, height = 7)


if(!is.null(genes_of_interest)) {
  goi_tab = read_tsv(genes_of_interest, comment = "#")

  goi_groups = list("genes_of_interest")
  group = "gene_name"

  if("group" %in% colnames(goi_tab)) {
    goi_groups = as.list(unique(goi_tab$group))
  } else {
    goi_tab$group <- "genes_of_interest"
  }

  lapply(goi_groups, function(grp) {
    goi = goi_tab |> filter(group == grp)
    goi = goi %>% filter(gene_name %in% resIHW$gene_name)

    goi_data = resIHW %>% filter(gene_name %in% goi$gene_name)

    p <- EnhancedVolcano(resIHW,
                    lab = resIHW$gene_name,
                    selectLab = goi$gene_name,
                    labSize = 5,
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = fdr_cutoff,
                    FCcutoff = fc_cutoff,
                    subtitle = "",
                    legendPosition = "top",
                    drawConnectors = TRUE,
                    widthConnectors = 0.25,
                    legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p-value ~ and~ log[2] ~ FC)),
                    caption = paste0("fold change cutoff: ", round(2**fc_cutoff, 1), ", adj.p-value cutoff: ", fdr_cutoff),
                    maxoverlapsConnectors = Inf,
                    #title = plot_title
                    ) + 
                    ggtitle("hAO_inf_hAO_ctrl") + 
                    theme(plot.title = element_text(color="black", size=15, face="bold",hjust = 0.5)) 
    
                    
                    
               

    grp_fname = gsub("[^[:alnum:]]+", "_", grp)
    save_plot(file.path(results_dir, paste0(prefix, "_", grp_fname, "_volcano_padj")), p, width = 8, height = 8)

    write_tsv(goi_data, file.path(results_dir, paste0(prefix, "_", grp_fname, ".tsv")))
    write_xlsx(goi_data, file.path(results_dir, paste0(prefix, "_" , grp_fname, ".xlsx")))

  })

}

# Save R ws
if (save_ws) {
  save.image(file = file.path(results_dir, paste0(prefix, "_ws.RData")))
}
