#!/bin/bash -ue
01_heatmap.R \
--count_mat=reactive_oxygen_species_metabolic_process_count_mat.csv \
--de_res=hAO_inf_hAO_ctrl_IHWsigFCgenes.tsv \
--prefix=reactive_oxygen_species_metabolic_process_count_mat
