#!/bin/bash -ue
01_heatmap.R \
--count_mat=cilium_movement_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes_2_fold.tsv \
--prefix=cilium_movement_count_mat
