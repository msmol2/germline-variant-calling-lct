#!/usr/bin/env bash
# Germline short-variant calling on chromosome 2 (GRCh37) for one 1000 Genomes run,
# then a look-up of the lactase persistence SNP rs4988235.
#
# Usage:  bash scripts/variant_calling.sh <run_accession> <sample_name> [threads]
# e.g.    bash scripts/variant_calling.sh ERR229776 HG00101 8
#
# Expects, relative to the working directory:
#   reads/<run>_1.fastq.gz, reads/<run>_2.fastq.gz
#   ref/human_g1k_v37.fasta      (1000 Genomes GRCh37 reference, contigs named 1..22,X,Y,MT)
# Tools: fastp, bwa, samtools, gatk (4.x), bcftools

set -euo pipefail

RUN=${1:?run accession}
SM=${2:?sample name}
T=${3:-8}
REF_FULL=ref/human_g1k_v37.fasta
REF=ref/human_g1k_v37_chr2.fasta
TMP=${TMPDIR:-$HOME/scratch/java_tmp}
mkdir -p "$TMP" fastp bam vcf logs

# 1. Trim and QC reads
fastp -i reads/${RUN}_1.fastq.gz -I reads/${RUN}_2.fastq.gz \
      -o fastp/${RUN}_1.trimmed.fastq.gz -O fastp/${RUN}_2.trimmed.fastq.gz \
      --html fastp/${RUN}_fastp.html --json fastp/${RUN}_fastp.json -w "$T" \
      > logs/${RUN}_fastp.log 2>&1

# 2. Reference: extract chromosome 2 and build the indexes (once)
if [[ ! -s ${REF}.bwt ]]; then
    samtools faidx "$REF_FULL"
    samtools faidx "$REF_FULL" 2 > "$REF"
    samtools faidx "$REF"
    bwa index "$REF" > logs/bwa_index.log 2>&1
    gatk CreateSequenceDictionary --java-options "-Djava.io.tmpdir=$TMP" -R "$REF" \
         > logs/gatk_dict.log 2>&1
fi

# 3. Align with a read group (required by GATK) and sort
bwa mem -t "$T" \
    -R "@RG\tID:${RUN}\tSM:${SM}\tPL:ILLUMINA\tLB:${RUN}\tPU:${RUN}" \
    "$REF" fastp/${RUN}_1.trimmed.fastq.gz fastp/${RUN}_2.trimmed.fastq.gz \
    2> logs/${RUN}_bwa.log \
  | samtools sort -@ "$T" -o bam/${RUN}_chr2.sorted.bam -
samtools index bam/${RUN}_chr2.sorted.bam
samtools flagstat bam/${RUN}_chr2.sorted.bam > logs/${RUN}_flagstat.txt

# 4. Call variants
gatk HaplotypeCaller --java-options "-Djava.io.tmpdir=$TMP" \
     -R "$REF" -I bam/${RUN}_chr2.sorted.bam -O vcf/${RUN}_chr2.vcf.gz \
     > logs/${RUN}_gatk.log 2>&1

# 5. Lactase persistence SNP rs4988235 (LCT -13910 C>T; GRCh37 2:136608646, G>A on the + strand)
echo "== rs4988235 call"
bcftools view -H -r 2:136608646 vcf/${RUN}_chr2.vcf.gz || true
echo "== read evidence at 2:136608646"
samtools mpileup -r 2:136608646-136608646 -f "$REF" bam/${RUN}_chr2.sorted.bam 2>/dev/null
