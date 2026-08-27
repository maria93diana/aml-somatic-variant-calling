process SAMTOOLS_FAIDX {
    tag "hg38"

    conda 'bioconda::samtools=1.21'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}",      emit: fasta
    path "${fasta}.fai",  emit: fai
    path "${fasta}.gzi",  emit: gzi

    script:
    """
    samtools faidx ${fasta}
    """
}