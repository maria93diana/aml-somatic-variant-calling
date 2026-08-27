process PICARD_DICT {
    tag "hg38"

    conda 'bioconda::picard=3.4.0'
    container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "*.dict", emit: dict

    script:
    """
    picard CreateSequenceDictionary \
        -R ${fasta} \
        -O ${fasta.baseName}.dict
    """
}