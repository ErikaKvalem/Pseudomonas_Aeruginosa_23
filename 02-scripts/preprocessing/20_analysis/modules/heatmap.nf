nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process HEATMAP{
    publishDir "${out_dir}", mode: "$mode"

    input:
    tuple val(id), path(count_mat), path(de_res) //, path(de_res_all)


    output:
    path("*heatmap.svg"), emit: heatmap_pdf, optional: true

    script:    
    """
    01_heatmap.R \\
    --count_mat=${count_mat} \\
    --de_res=${de_res} \\
    --prefix=${id} 
    """


}