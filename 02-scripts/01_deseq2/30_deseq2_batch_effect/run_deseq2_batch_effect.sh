Rscript runDESeq2_ICBI.R /data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/samplesheet_deseq_batch_effect.csv  /data/scratch/kvalem/projects/2023/pseudomonas_aeruginosa/tables/salmon.merged.gene_counts.tsv --result_dir=/data/projects/2023/Pseudomonas_aeruginosa/20_deseq2icbi/04_deseq2_batch_effect --c1=hAO_inf --c2=hAO_ctrl --sample_col=sample --condition_col=group --remove_batch_effect --batch_col=batch --fdr_cutoff=0.1 --fc_cutoff=0.5