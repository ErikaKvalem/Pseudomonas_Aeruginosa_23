#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COUNTMAT} from "./modules/countmat"
include { HEATMAP} from "./modules/heatmap"

workflow {
    log1p_cpm_data= Channel.fromPath(params.log1p_cpm_input_path)
    gsea_bp = Channel.fromPath(params.gsea_input_path)
    gene_set = Channel.value(params.geneset_of_interest).flatten()
    //genes = Channel.value(params.genes_of_interest)
    //de_res_sig = Channel.fromPath(params.de_res_sig_input_path)
    de_res = Channel.fromPath(params.de_res_input_path)

    ch_combined = gene_set.combine(log1p_cpm_data)

    ch_input_countmat = ch_combined.combine(gsea_bp)

    COUNTMAT(ch_input_countmat)

    ch_out_count_mat = COUNTMAT.out.count_mat.flatten().map { 
         it -> [it.baseName.replace("_count_mat.csv", ""), it]
    }

    ch_input_heatmap = ch_out_count_mat.combine(de_res)

    HEATMAP(ch_input_heatmap)

    
}
