#!/bin/bash -ue
01_heatmap.R \
--count_mat=nitric_oxide_metabolic_process_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes.tsv \
--prefix=nitric_oxide_metabolic_process_count_mat
