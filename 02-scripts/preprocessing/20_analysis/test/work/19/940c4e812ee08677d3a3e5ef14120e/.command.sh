#!/bin/bash -ue
00_countmat.py \
--log1p_cpm_data=log2_cpm_counts_plus_one.tsv \
--gsea_bp=hAO_inf_hAO_ctrl_GSEA_GO_BP.tsv \
--gene_set='reactive oxygen species metabolic process'
