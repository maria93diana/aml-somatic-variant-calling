process BWAMEM2_INDEX {
    tag "hg38"
    
    conda 'bioconda::bwa-mem2=2.2.1'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'
    
    publishDir "${params.outdir}/reference", mode: 'copy'
    
    input:
    path fasta
    
    output:
    path "*.{0123,amb,ann,bwt.2bit.64,pac}", emit: index
    
    script:
    """
    bwa-mem2 index ${fasta}
    """
}