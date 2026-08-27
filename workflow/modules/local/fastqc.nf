process FASTQC {
    tag "$sample_id"

    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip
    
    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}