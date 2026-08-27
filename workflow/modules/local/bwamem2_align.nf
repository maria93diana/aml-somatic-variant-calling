process BWAMEM2_ALIGN {
    tag "$sample_id"
    
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21'
    container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:1a23bc9ce572f1073c21c84a58e10a04b2498d1f-0'
    
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    path fasta
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai"), emit: bai
    
    script:
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${fasta} \\
        ${r1} ${r2} | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    
    samtools index ${sample_id}.bam
    """
}