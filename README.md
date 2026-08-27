# Project 02 — Somatic Variant Calling Pipeline (Mutect2)

A Nextflow DSL2 pipeline for somatic SNV/indel detection from tumor-normal paired WES data,
using GATK Mutect2 with orientation artifact correction.

**Dataset:** HCC1143 breast cancer cell line (tumor) + HCC1143BL matched normal (Illumina WES, chr22 subset)  
**Reference:** GRCh38/hg38 chromosome 22  

---

## Pipeline overview

```
FASTQ (tumor + normal)
    │
    ├── FastQC (raw QC)
    ├── Trim Galore (adapter trimming)
    ├── FastQC (post-trim QC)
    │
    ├── BWA-MEM2 index (reference indexing)
    ├── BWA-MEM2 align (tumor → BAM, normal → BAM)
    │
    ├── GATK Mutect2 (tumor-normal somatic calling)
    │       outputs: raw VCF + f1r2 orientation counts + stats
    │
    ├── GATK LearnReadOrientationModel (artifact prior from f1r2)
    ├── GATK FilterMutectCalls (apply filters → PASS/FAIL)
    │
    └── VEP (variant annotation — Linux/Docker only; skipped on Apple Silicon)
```

**Results on chr22:** 193 raw candidates → 54 PASS somatic variants

---

## Why tumor-normal paired calling?

Sequencing a matched normal from the same patient lets Mutect2 subtract germline variants
(variants the patient was born with) from the callset. Without a normal, germline heterozygous
SNPs — of which there are ~3 million per person — would all appear as candidates.
The normal sample is NOT used as a control for quality; it is used as a **germline reference**.

---

## Key biological concepts

### Somatic variant allele frequency (VAF)
A somatic mutation is not present in every cell. Tumors are heterogeneous:
- A clonal mutation (present in all tumor cells) appears at ~50% VAF in a pure tumor sample
- A subclonal mutation (in only a fraction of cells) appears at lower VAF (5–30%)
- Normal contamination dilutes VAF further

In this dataset (chr22), PASS variants have VAFs ranging from ~10% to ~40%, consistent
with a heterogeneous tumor sample with some normal cell admixture.

### Orientation artifacts (OxoG / FFPE)
DNA damage during library preparation causes systematic base changes:
- G→T transversions on one strand (oxidative damage, "OxoG")
- C→T transitions (FFPE formalin fixation)

These look like real variants but appear almost exclusively on reads in one orientation
(F1R2 or F2R1). Mutect2 collects F1R2/F2R1 counts, then LearnReadOrientationModel
builds a prior probability of artifact vs. true variant for each trinucleotide context.
FilterMutectCalls uses these priors to flag orientation artifacts.

---

## Understanding the output files

### `HCC1143_tumor.vcf.gz` — raw Mutect2 calls

VCF format: 8 fixed columns + FORMAT definition + one column per sample (NORMAL first, TUMOR second).

**Key INFO fields:**

| Field | Good value | Meaning |
|-------|-----------|---------|
| `TLOD` | > 6.3 (default threshold) | Tumor log-odds: how much evidence FOR the somatic mutation. Higher = more confident. This is Mutect2's primary signal. |
| `NLOD` | > 2.2 | Normal log-odds: how "clean" the normal is. Low NLOD = normal also has alt reads = likely germline or artifact. |
| `NALOD` | > 0 (positive) | Evidence this is NOT a germline variant. Negative values mean the normal carries some alt allele — not an automatic fail, but a warning. |
| `POPAF` | > 5 (rare) | Population allele frequency in gnomAD (log scale). 6.0 = not in databases = novel somatic candidate. Low values (< 2) = common in population = likely germline. |
| `ECNT` | any | Number of variants in the same local cluster. High ECNT can indicate a complex event or a misalignment artifact. |
| `MBQ` | > 20 | Median base quality of alt reads. Low values suggest poor sequencing quality supporting the alt. |
| `MMQ` | 60 | Median mapping quality of alt reads. 60 = uniquely mapped. Low values (< 30) = alt reads map ambiguously = unreliable. |

**Key FORMAT fields (per sample):**

| Field | Example | Meaning |
|-------|---------|---------|
| `GT` | `0/1` | Genotype: 0=ref, 1=alt. Tumor is usually 0/1 (heterozygous somatic). Normal should be 0/0. |
| `AD` | `22,7` | Allele depth: ref reads, alt reads. Tumor 22,7 means 7 out of 29 reads support the mutation. |
| `AF` | `0.265` | Allele frequency = alt/(ref+alt) = 7/29 = 26.5%. This is the VAF. |
| `DP` | `29` | Total depth at this site. Low depth (< 10) = unreliable call. |
| `F1R2` | `11,4` | Reads with alt on forward strand (F1) or reverse strand (R2). Used for orientation artifact detection. Ideally balanced between F1R2 and F2R1. |
| `SB` | `3,19,2,5` | Strand bias counts: ref-fwd, ref-rev, alt-fwd, alt-rev. Severe imbalance (all alt reads on one strand) = likely artifact. |

### `HCC1143_tumor.vcf.gz.stats`

```
callable   903201.0
```

Number of genomic positions where Mutect2 had enough coverage to make a call (~903 kb).
Used to compute **tumor mutational burden (TMB)** = PASS variants / callable Mb.
Here: 54 / 0.9 Mb ≈ 60 mut/Mb (chr22 is capture-enriched; do not extrapolate to whole genome).

### `HCC1143_tumor.f1r2.tar.gz`

Binary archive of per-base orientation counts (F1R2 vs F2R1 reads at each alt allele).
Input to LearnReadOrientationModel. Not human-readable, but essential for artifact correction.

### `HCC1143_tumor.filtering.stats` — what was filtered and why

**Metadata section** (lines starting with `#`):
- `threshold=0.444` — variants with somatic posterior probability < 44.4% are filtered
- `sensitivity=0.939` — the model estimates it kept 93.9% of true somatic variants
- `fdr=0.077` — estimated 7.7% false discovery rate among PASS variants

**Filter table:**

| Filter | What it catches | FP removed | FN cost |
|--------|----------------|-----------|---------|
| `weak_evidence` | Low TLOD — not enough supporting reads. The largest filter in low-depth data. | ~3.1 | ~2.0 |
| `haplotype` | Variants inconsistent with the surrounding haplotype — often clustered sequencing errors. | ~0.9 | ~0.7 |
| `strand_bias` | Alt allele only seen on reads in one orientation (forward or reverse). | ~0.2 | ~0.0 |
| `germline` | Posterior probability of being germline is too high. | ~0.2 | ~0.0 |
| `orientation` | Caught by LearnReadOrientationModel — oxidative or FFPE damage pattern. | ~0.02 | ~0.0 |
| `normal_artifact` | Artifact also present in the normal sample. | ~0.0 | ~0.7 |

**FP** = estimated false positives correctly removed  
**FN** = estimated real somatic variants accidentally removed (false negatives)  
The model optimizes this tradeoff; you can shift it with `--max-false-positive-rate`.

### `HCC1143_tumor.filtered.vcf.gz` — final callset

Same format as the raw VCF, with two differences:
- FILTER column is now `PASS` (kept) or a filter name like `weak_evidence` (removed)
- Two new INFO fields added by FilterMutectCalls:
  - `GERMQ` — germline quality score. Higher = more confident it is somatic, not germline. Values > 30 are reliable.
  - `ROQ` — read orientation quality. Higher = less likely to be an orientation artifact.

---

## How to run

### Prerequisites

- Nextflow ≥ 24.0
- conda (for `-profile local`)
- Docker (for `-profile docker`, required on Linux for VEP)

### Download input data

```bash
# Reference: chr22 from hg38 (BGZF format required by GATK)
mkdir -p data/reference
curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
    | gzcat | bgzip > data/reference/chr22.fa.gz

# HCC1143 WES reads (tumor and normal, chr22 region)
# Download scripts in scripts/part1_variant_calling/
mkdir -p data/hcc1143
```

### Run the pipeline

```bash
# Local (conda, Apple Silicon / Intel Mac)
nextflow run workflow/main.nf -profile local

# Docker (Linux, CI, cloud) — enables VEP annotation
nextflow run workflow/main.nf -profile docker

# Resume after interruption
nextflow run workflow/main.nf -profile local -resume
```

### Enable VEP annotation (Linux/Docker only)

VEP conda package does not run on Apple Silicon (osx-arm64 segfault).
On Linux with Docker:

```bash
nextflow run workflow/main.nf -profile docker --run_vep true
```

---

## Portability

| Profile | Execution | Notes |
|---------|-----------|-------|
| `local` | conda | Works on macOS (Intel + Apple Silicon), Linux |
| `docker` | containers | Recommended for production; enables VEP |
| `singularity` | containers | HPC clusters without Docker |

---

## Project structure

```
aml-somatic-variant-calling/
├── workflow/
│   ├── main.nf                    # Main workflow
│   ├── nextflow.config            # Parameters, resources, profiles
│   └── modules/local/
│       ├── fastqc.nf
│       ├── trimgalore.nf
│       ├── bwamem2_index.nf
│       ├── bwamem2_align.nf
│       ├── samtools_faidx.nf
│       ├── picard_dict.nf
│       ├── mutect2.nf
│       ├── filtermutectcalls.nf
│       └── vep.nf
├── data/
│   ├── reference/chr22.fa.gz      # hg38 chr22 (BGZF)
│   └── hcc1143/                   # Tumor + normal FASTQ
└── results/
    ├── fastqc/                    # Raw and post-trim QC reports
    ├── trimgalore/                # Trimmed reads + trim reports
    ├── bwamem2/                   # Aligned BAMs
    └── mutect2/                   # VCFs, filtering stats, f1r2
```

---

## Known limitations / Apple Silicon notes

- `libgkl_compression.dylib` warnings from GATK are harmless on Apple Silicon — GATK falls back to Java implementations automatically and produces correct results.
- `IntelInflater`/`IntelDeflater` warnings are similarly harmless.
- VEP (`ensembl-vep`) conda package segfaults on osx-arm64. Use `-profile docker` on Linux for annotation.
