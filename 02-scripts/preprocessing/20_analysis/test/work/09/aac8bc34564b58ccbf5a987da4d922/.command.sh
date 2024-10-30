#!/bin/bash -ue
01_heatmap.R \
--count_mat=antimicrobial_humoral_response_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes.tsv \
--prefix=antimicrobial_humoral_response_count_mat
