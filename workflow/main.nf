#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { TRIMGALORE             } from './modules/local/trimgalore'
include { BWAMEM2_INDEX          } from './modules/local/bwamem2_index'
include { BWAMEM2_ALIGN          } from './modules/local/bwamem2_align'
include { SAMTOOLS_FAIDX         } from './modules/local/samtools_faidx'
include { PICARD_DICT            } from './modules/local/picard_dict'
include { MUTECT2                } from './modules/local/mutect2'
include { FILTERMUTECTCALLS      } from './modules/local/filtermutectcalls'
include { VEP                    } from './modules/local/vep'

workflow {

    // ── Reference preparation ─────────────────────────────────────────────────
    fasta_ch = Channel.fromPath( params.fasta, checkIfExists: true ).first()

    BWAMEM2_INDEX( fasta_ch )
    SAMTOOLS_FAIDX( fasta_ch )
    PICARD_DICT( fasta_ch )

    // Bundle all reference files into one tuple — GATK tools need fasta + fai + gzi + dict
    reference_ch = SAMTOOLS_FAIDX.out.fasta
        .combine( SAMTOOLS_FAIDX.out.fai  )
        .combine( SAMTOOLS_FAIDX.out.gzi  )
        .combine( PICARD_DICT.out.dict    )
        .first()

    // ── Input reads ───────────────────────────────────────────────────────────
    // Tumor and normal are defined explicitly (not discovered by glob)
    // because we need to label which sample is which for Mutect2
    tumor_ch  = Channel.of( tuple(params.tumor_name,  file(params.tumor_r1),  file(params.tumor_r2))  )
    normal_ch = Channel.of( tuple(params.normal_name, file(params.normal_r1), file(params.normal_r2)) )

    // Mix both into one channel so QC and trimming run on both samples together
    all_reads_ch = tumor_ch.mix( normal_ch )

    // ── QC and trimming ───────────────────────────────────────────────────────
    // FastQC expects: tuple val(sample_id), path(reads)  — reads as a list
    FASTQC_RAW(
        all_reads_ch.map { id, r1, r2 -> tuple(id, [r1, r2]) }
    )

    TRIMGALORE(
        all_reads_ch.map { id, r1, r2 -> tuple(id, [r1, r2]) }
    )

    FASTQC_TRIM(
        TRIMGALORE.out.trimmed_reads.map { id, r1, r2 -> tuple(id, [r1, r2]) }
    )

    // ── Alignment ─────────────────────────────────────────────────────────────
    // Both tumor and normal go through the same process — they are two items
    // in the trimmed_reads channel, processed in parallel
    BWAMEM2_ALIGN(
        TRIMGALORE.out.trimmed_reads,
        fasta_ch,
        BWAMEM2_INDEX.out.index.collect()
    )

    // Join BAM and BAI by sample_id, then split tumor from normal by name
    bam_ch = BWAMEM2_ALIGN.out.bam.join( BWAMEM2_ALIGN.out.bai )

    tumor_bam_ch  = bam_ch.filter { it[0] == params.tumor_name  }
    normal_bam_ch = bam_ch.filter { it[0] == params.normal_name }

    // ── Somatic variant calling ───────────────────────────────────────────────
    MUTECT2(
        tumor_bam_ch,
        normal_bam_ch,
        reference_ch
    )

    FILTERMUTECTCALLS(
        MUTECT2.out.vcf.join( MUTECT2.out.tbi ).join( MUTECT2.out.stats ),
        MUTECT2.out.f1r2,
        reference_ch
    )

    // ── Variant annotation ────────────────────────────────────────────────────
    // VEP conda package does not run on Apple Silicon (osx-arm64).
    // Skip locally with --run_vep false (default); enable with -profile docker on Linux.
    if ( params.run_vep ) {
        VEP(
            FILTERMUTECTCALLS.out.vcf.join( FILTERMUTECTCALLS.out.tbi )
        )
    }
}
