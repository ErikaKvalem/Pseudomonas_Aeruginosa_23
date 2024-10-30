#!/bin/bash -ue
01_heatmap.R \
--count_mat=response_to_lipopolysaccharide_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes.tsv \
--prefix=response_to_lipopolysaccharide_count_mat
