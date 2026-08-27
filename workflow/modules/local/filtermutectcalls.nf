process FILTERMUTECTCALLS {
    tag "$sample_id"

    conda 'bioconda::gatk4=4.6.1.0'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py36hdfd78af_0'

    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi), path(stats)
    tuple val(sample_id2), path(f1r2)
    tuple path(fasta), path(fai), path(gzi), path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz.tbi"), emit: tbi
    tuple val(sample_id), path("${sample_id}.filtering.stats"),      emit: stats

    script:
    """
    ln -sf ${dict} ${fasta.simpleName}.dict

    gatk LearnReadOrientationModel \\
        -I ${f1r2} \\
        -O artifact_prior.tar.gz

    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${vcf} \\
        --ob-priors artifact_prior.tar.gz \\
        --filtering-stats ${sample_id}.filtering.stats \\
        -O ${sample_id}.filtered.vcf.gz
    """
}
