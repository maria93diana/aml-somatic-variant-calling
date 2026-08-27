process MUTECT2 {
    tag "$sample_id"

    conda 'bioconda::gatk4=4.6.1.0'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py36hdfd78af_0'

    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    tuple val(sample_id), path(tumor_bam), path(tumor_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
    tuple path(fasta), path(fai), path(gzi), path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"),      emit: vcf
    tuple val(sample_id), path("${sample_id}.vcf.gz.tbi"),  emit: tbi
    tuple val(sample_id), path("${sample_id}.vcf.gz.stats"), emit: stats
    tuple val(sample_id), path("${sample_id}.f1r2.tar.gz"), emit: f1r2

    script:
    """
    # GATK expects the dict named <ref>.dict (no .fa extension).
    # PICARD_DICT produces <ref>.fa.dict, so we symlink to the expected name.
    ln -sf ${dict} ${fasta.simpleName}.dict

    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal ${normal_id} \\
        --f1r2-tar-gz ${sample_id}.f1r2.tar.gz \\
        -O ${sample_id}.vcf.gz
    """
}
