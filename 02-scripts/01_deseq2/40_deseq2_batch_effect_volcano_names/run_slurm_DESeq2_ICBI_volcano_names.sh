#!/bin/bash
#SBATCH --job-name=DESeq2_volcano
#SBATCH --output=DESeq2_volcano.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=00:02:00
module load R
Rscript runDESeq2_ICBI_volcano_names.R /data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/samplesheet_deseq_batch_effect.csv  /data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/salmon.merged.gene_counts.tsv --result_dir=/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/08_cilium_volcano_names --genes_of_interest=/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/04_deseq2_batch_effect/volcano_names/genes_of_interest_volcano.csv --c1=hAO_inf --c2=hAO_ctrl --sample_col=sample --condition_col=group --remove_batch_effect --batch_col=batch --fdr_cutoff=0.1 --fc_cutoff=0.5
