process VEP {
    tag "$sample_id"

    // conda package segfaults on Apple Silicon (osx-arm64) — use -profile docker on Linux
    conda 'bioconda::ensembl-vep=113.0'
    container 'quay.io/biocontainers/ensembl-vep:113.0--pl5321h2a08d89_0'

    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz.tbi"), emit: tbi
    tuple val(sample_id), path("${sample_id}.vep_summary.html"), emit: summary

    script:
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --everything \\
        --offline \\
        --cache \\
        --cache_version 113 \\
        --assembly GRCh38 \\
        --fork ${task.cpus} \\
        --stats_file ${sample_id}.vep_summary.html

    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
}
