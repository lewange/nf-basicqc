/*
========================================================================================
    MULTIQC Module
========================================================================================
    MultiQC - Aggregate results from bioinformatics analyses
*/

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path multiqc_files
    path generated_config

    output:
    path "*_multiqc_report.html", emit: report
    path "*_data"               , emit: data
    path "*_plots"              , emit: plots, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def title = params.multiqc_title ? "--title \"${params.multiqc_title}\"" : ''
    // Use generated config, but allow user override
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : "--config ${generated_config}"
    def prefix = params.project_name ? params.project_name.replaceAll('\\s+', '_') : (params.multiqc_title ? params.multiqc_title.replaceAll('\\s+', '_') : 'basicqc')
    """
    multiqc \\
        --force \\
        $title \\
        $config \\
        --filename ${prefix}_multiqc_report \\
        $args \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    def prefix = params.multiqc_title ? params.multiqc_title.replaceAll('\\s+', '_') : 'basicqc'
    """
    mkdir ${prefix}_data
    touch ${prefix}_multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.25.1
    END_VERSIONS
    """
}
