#!/bin/bash -ue
01_heatmap.R \
--count_mat=cilium_organization_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes.tsv \
--prefix=cilium_organization_count_mat
