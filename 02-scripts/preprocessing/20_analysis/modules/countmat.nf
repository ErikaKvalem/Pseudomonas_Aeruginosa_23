nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process COUNTMAT{
    publishDir "${out_dir}", mode: "$mode"

    input:

    tuple  val (gene_set), path(log1p_cpm_data), path(gsea_bp)
    


    output:
    path("*count_mat.csv"), emit: count_mat optional true

    script:   
    """
    00_countmat.py \\
    --log1p_cpm_data=${log1p_cpm_data} \\
    --gsea_bp=${gsea_bp} \\
    --gene_set='${gene_set}'
    """

  


}