#!/bin/bash -ue
01_heatmap_copy.R \
--count_mat=antimicrobial_humoral_response_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigGenes.tsv \
--prefix=antimicrobial_humoral_response_count_mat
