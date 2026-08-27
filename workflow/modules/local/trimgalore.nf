process TRIMGALORE {
    tag "$sample_id"
    
    conda 'bioconda::trim-galore=0.6.10'
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    
    publishDir "${params.outdir}/trimming", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*_trimming_report.txt"),                 emit: reports
    
    script:
    """
    trim_galore \\
        --paired \\
        --quality 20 \\
        --length 36 \\
        --cores ${task.cpus} \\
        ${reads}
    """
}