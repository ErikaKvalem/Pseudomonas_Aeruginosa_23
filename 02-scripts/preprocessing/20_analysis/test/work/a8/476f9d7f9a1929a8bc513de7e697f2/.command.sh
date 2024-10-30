#!/bin/bash -ue
01_heatmap_copy.R \
--count_mat=response_to_bacterium_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigGenes.tsv \
--prefix=response_to_bacterium_count_mat
