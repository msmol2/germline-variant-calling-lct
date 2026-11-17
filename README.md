# Germline variant calling: the lactase persistence SNP

A short-read germline variant-calling pipeline (fastp → BWA-MEM → samtools → GATK HaplotypeCaller) run on chromosome 2 of a low-coverage **1000 Genomes** genome. It ends by genotyping **rs4988235**, the regulatory SNP upstream of *LCT* that is the main determinant of lactase persistence in Europeans.

Written for the Genomics & Transcriptomics course of my MSc (IMC Krems, 2026). The pipeline ran on the university HPC server. This repository holds the pipeline as one script, plus notes on each step and on reading the output in IGV.

## Background

Lactase (*LCT*) is normally switched off after weaning. In many Europeans it stays on, because of a C>T change 13,910 bp upstream of *LCT*, inside an intron of the neighbouring gene *MCM6*. The T allele creates a stronger enhancer, and it spread with dairy farming over the last ~10,000 years.

| | |
|---|---|
| SNP | rs4988235 (LCT −13910 C>T) |
| Position (GRCh37) | 2:136,608,646 |
| Alleles on the + strand | G (ref) > A (alt). This is the C>T change read on the *LCT* strand. |
| CC (GG) | lactase non-persistent |
| CT / TT (GA / AA) | lactase persistent (dominant) |

Note that the VCF reports the plus-strand alleles **G/A**, not the C/T used in the literature.

## Data

| | |
|---|---|
| Reads | ENA run [ERR229776](https://www.ebi.ac.uk/ena/browser/view/ERR229776): Illumina HiSeq 2000, paired-end WGS, 111 M reads, low coverage (~4–8×) |
| Individual | 1000 Genomes **HG00101** (GBR, British in England and Scotland) |
| Reference | 1000 Genomes GRCh37 `human_g1k_v37.fasta`, chromosome 2 only |

Only chromosome 2 was used, to keep alignment and calling to about an hour each. The class read-group string set `SM:HG00096`, but the ENA record for ERR229776 belongs to HG00101. The script takes the sample name as an argument, so the label can be set correctly.

## Pipeline

```
reads/*.fastq.gz ─ fastp ─▶ trimmed reads ─ bwa mem (+ read group) ─▶ samtools sort/index ─▶ BAM
ref/human_g1k_v37.fasta ─ samtools faidx (chr 2) ─▶ bwa index, .fai, .dict ─────────────────┘
                                                                                             │
                                               gatk HaplotypeCaller ◀────────────────────────┘
                                                        │
                                                        ▼
                                              VCF ─ bcftools / IGV ─▶ genotype at rs4988235
```

| Step | Tool | Notes |
|---|---|---|
| Read QC and trimming | fastp | default thresholds; adapters auto-detected; HTML/JSON report |
| Reference prep | samtools faidx, bwa index, gatk CreateSequenceDictionary | chromosome 2 extracted (contig `2`, no `chr` prefix) |
| Alignment | bwa mem | read group `ID/SM/PL/LB/PU` added with `-R`, since HaplotypeCaller refuses BAMs without one |
| Sort, index, stats | samtools sort / index / flagstat | output piped straight into `sort`, so no intermediate SAM is written |
| Variant calling | gatk HaplotypeCaller | single sample, default settings; JVM temp dir set to avoid `/tmp` quotas |
| Look-up | bcftools view, samtools mpileup, IGV | call and raw read evidence at 2:136,608,646 |

## Run

```bash
# reads/ERR229776_{1,2}.fastq.gz and ref/human_g1k_v37.fasta in place
bash scripts/variant_calling.sh ERR229776 HG00101 8
```

Long steps were started on the server with `nohup … > log 2>&1 &` so they kept running after the SSH session closed.

## Checking the call in IGV

Load `ref/human_g1k_v37_chr2.fasta`, the sorted BAM (with its `.bai`) and the VCF, then go to `2:136608646`.
- About half the reads showing A: heterozygous (persistent).
- Almost all A: homozygous alternative (persistent).
- Almost all G: homozygous reference (non-persistent).

At 4–8× coverage a heterozygous site may be supported by only two or three reads, so check the depth (`DP`) and genotype quality (`GQ`) in the VCF before trusting the call.

[`docs/notes.md`](docs/notes.md) covers file formats, flagstat, the VCF columns and reading pile-ups.

## Limitations

- Low-coverage data: genotypes at single sites are uncertain. The 1000 Genomes project used population-level imputation for that reason.
- No duplicate marking or base-quality recalibration (BQSR). For a clinical-grade call, add `MarkDuplicates` and `BaseRecalibrator` before HaplotypeCaller.
- Chromosome 2 only. Reads from other chromosomes can mis-map onto chromosome 2.

## References

- Enattah N.S. *et al.* Identification of a variant associated with adult-type hypolactasia. *Nat Genet* 30:233–237 (2002).
- The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 526:68–74 (2015).
- Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997 (2013).
- Poplin R. *et al.* Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv 201178 (2018).
