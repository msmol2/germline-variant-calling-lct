# Notes: formats, QC and reading the output

Condensed from my course notes.

## File Format Quick Reference

| Format | Extension | Human readable? | Contents | Created by |
|--------|-----------|-----------------|----------|------------|
| FASTA | `.fasta`, `.fa` | ✅ Yes | Reference genome sequences | Public databases |
| FASTA index | `.fai` | ✅ Yes | Contig lengths and byte offsets | `samtools faidx` |
| FASTQ | `.fastq`, `.fastq.gz` | ✅ (when unzipped) | Sequencing reads + Phred quality scores | Sequencer / fastp (trimmed) |
| SAM | `.sam` | ✅ Yes | Aligned reads, tab-separated, with header | `bwa mem` |
| BAM | `.bam` | ❌ No (binary) | Compressed binary SAM | `samtools sort` |
| BAM index | `.bam.bai` | ❌ No (binary) | Random-access index for BAM | `samtools index` |
| Sequence dict | `.dict` | ✅ Yes | SAM header with contig names/lengths | `gatk CreateSequenceDictionary` |
| VCF | `.vcf` | ✅ Yes | Called variants (position, ref, alt, quality) | `gatk HaplotypeCaller` |

---


## Checking the BAM File Structure with SAMtools


```bash
# View the BAM header (metadata section at the top of the file)
samtools view -H ERR229776_chr2.sorted.bam

# View the first 20 alignment records (reads) — skips the header
samtools view ERR229776_chr2.sorted.bam | head -20

# View header + first 20 alignment lines combined
samtools view -H ERR229776_chr2.sorted.bam \
  && samtools view ERR229776_chr2.sorted.bam | head -20

# Get alignment statistics (mapped %, duplicate %, etc.)
samtools flagstat ERR229776_chr2.sorted.bam
```

**SAM/BAM column layout** (each alignment line, tab-separated):

| Column | Field | Example | Description |
|--------|-------|---------|-------------|
| 1 | QNAME | `ERR229776.1` | Read name |
| 2 | FLAG | `99` | Bitwise flag (paired, mate reversed, etc.) |
| 3 | RNAME | `2` | Reference chromosome name |
| 4 | POS | `136600000` | 1-based leftmost mapping position |
| 5 | MAPQ | `60` | Mapping quality (−10 log₁₀ probability wrong position) |
| 6 | CIGAR | `100M` | CIGAR string (match/ins/del operations) |
| 7 | RNEXT | `=` | Mate's chromosome (`=` if same) |
| 8 | PNEXT | `136600150` | Mate's position |
| 9 | TLEN | `250` | Observed template/insert length |
| 10 | SEQ | `ACGT...` | Read sequence |
| 11 | QUAL | `IIII...` | ASCII-encoded Phred quality scores |
| 12+ | Optional tags | `RG:Z:ERR229776` | Read group, alignment score, etc. |

**FLAG values** (common examples):

| FLAG | Meaning |
|------|---------|
| 99 | Read is paired, proper pair, mate reverse strand, read 1 |
| 147 | Read is paired, proper pair, read reverse strand, read 2 |
| 4 | Read is unmapped |

> The SAM flag is a bitfield; see `https://broadinstitute.github.io/picard/explain-flags.html` to decode any flag value.

---


## SAMtools flagstat — Alignment Summary Statistics


`samtools flagstat` parses the FLAG field of every read and produces a human-readable summary.

```bash
samtools flagstat ERR229776_chr2.sorted.bam
```

**Example output (illustrative values):**
```
4000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
3980000 + 0 mapped (99.50% : N/A)
4000000 + 0 paired in sequencing
2000000 + 0 read1
2000000 + 0 read2
3960000 + 0 properly paired (99.00% : N/A)
3980000 + 0 with itself and mate mapped
20000 + 0 singletons (0.50% : N/A)
0 + 0 with mate mapped to a different chr
```

**Key metrics:**

| Line | Meaning |
|------|---------|
| `in total` | Total reads in the BAM |
| `mapped` | Reads that aligned to the reference (% in brackets) |
| `properly paired` | Both mates mapped in expected orientation and distance |
| `singletons` | One mate mapped, the other did not |
| `with mate mapped to a different chr` | Discordant pairs — may indicate translocations or mapping errors |

---


## Understanding the VCF Format


The VCF (Variant Call Format) stores the results of variant calling.

**Structure:**

**1. Header lines** (start with `##`):
```
##fileformat=VCFv4.2
##GATK version
##reference=human_g1k_v37_chr2.fasta
##contig=<ID=2,length=243199373>
##INFO=<ID=DP,...>
##FORMAT=<ID=GT,...>
```
The header defines the reference genome, tool used, and describes all INFO and FORMAT fields used in the data lines.

**2. Column header line** (starts with `#CHROM`):
```
#CHROM  POS     ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  ERR229776
```

**3. Variant records** (data lines, one per called variant position):
```
2   136608646   .   G   A   856.77   PASS   DP=50;AF=0.5   GT:GQ:DP   0/1:99:50
```

**Column descriptions:**

| Column | Name | Example | Description |
|--------|------|---------|-------------|
| 1 | CHROM | `2` | Chromosome |
| 2 | POS | `136608646` | 1-based position of the variant |
| 3 | ID | `.` | dbSNP identifier (`.` if not in dbSNP) |
| 4 | REF | `G` | Reference allele at this position |
| 5 | ALT | `A` | Alternative allele(s) called in the sample |
| 6 | QUAL | `856.77` | Phred-scaled variant quality score |
| 7 | FILTER | `PASS` | `PASS` = passed all filters; otherwise filter name(s) |
| 8 | INFO | `DP=50;...` | Variant-level statistics (depth, allele frequency, etc.) |
| 9 | FORMAT | `GT:GQ:DP` | Keys describing the sample column(s) |
| 10+ | Sample(s) | `0/1:99:50` | Per-sample genotype data, keyed by FORMAT column |

**Genotype (GT) encoding:**

| GT | Meaning |
|----|---------|
| `0/0` | Homozygous reference (same as reference at both alleles) |
| `0/1` | Heterozygous (one allele = reference, one = alternative) |
| `1/1` | Homozygous alternative (both alleles carry the variant) |
| `./` . | Missing data (insufficient coverage to call) |

**Key VCF INFO fields (GATK HaplotypeCaller):**

| Tag | Meaning |
|-----|---------|
| `DP` | Total read depth at this position |
| `AF` | Allele frequency of the alternative allele |
| `MQ` | Mean mapping quality |
| `QD` | Quality by depth (QUAL/DP — normalises quality for depth) |
| `FS` | Fisher's exact test for strand bias |

> **Every row in the VCF is a position where at least one non-reference allele was called.** Positions where all reads match the reference are not recorded. The reference state is always given in the REF column so results can be interpreted relative to the known reference.

---


## Interpreting the IGV Read Display


When zoomed in to the read level, each horizontal bar represents one sequencing read. Key visual features:

**Coverage track (top):**  
A bar chart showing how many reads overlap each position. Height = depth of coverage. The Y-axis scale is shown on the left. Regions with no coverage appear as gaps.

**Read arrows:**  
Individual reads are shown as horizontal arrows pointing left (←) or right (→):
- Arrow direction indicates the read **strand**: right = forward strand (5′→3′), left = reverse strand
- **Paired-end reads** from the same fragment appear in pairs, often connected by a thin line spanning the insert

**Base colouring:**  
By default, reads matching the reference are shown in grey. Any base that **differs from the reference** is shown in the colour of that base (A=green, T=red, G=orange, C=blue). This immediately highlights polymorphic positions.

**Coverage bar colours:**  
In the coverage bar, positions with mismatches are shown with a colour band at the top of the bar, proportional to the frequency of the non-reference allele.

**Interpreting a candidate SNP:**  
If you see ~half the reads showing a non-reference base at a position, and the other half matching the reference, this is characteristic of a **heterozygous germline variant** (one allele from each parent). If all reads show the alternative base, it suggests homozygous alternative.

---

## Viewing Mate Pairs in IGV


By default, IGV shows forward and reverse reads independently. You can display **read pairs** together to visualise the insert size:

1. Right-click on the reads track
2. Select **View as pairs**

**What changes:**  
Each forward read and its corresponding reverse read are visually connected, spanning the full insert. The grey connecting line shows the inferred insert region not covered by either read.

**What you can infer:**
- **Insert size:** The total span of the pair = fragment length. Typical values: 150–400 bp for standard Illumina library prep.
- **Orientation:** In a normal pair-end library, reads point **inward** toward each other (→ ←). Reads pointing the same direction (→ →) or outward (← →) suggest structural variants (inversions, translocations).

---

## Identifying Variants Visually in IGV


**To investigate a known SNP position (e.g. the LCT −13910 C/T SNP, rs4988235):**

IGV and the VCF show the **plus strand**, where this SNP is **G (ref) / A (alt)**. The C/T in the literature is the same change read on the *LCT* strand.

Navigate to the exact coordinate:
```
2:136,608,646
```

**What to look for:**

| Observation | Interpretation |
|-------------|---------------|
| ~50% of reads show G, ~50% show A | Heterozygous G/A (C/T): one non-persistent and one persistent allele |
| ~100% of reads show A | Homozygous A/A (T/T): lactase persistent |
| ~100% of reads show G | Homozygous G/G (C/C): lactase non-persistent (reference) |
| 1–2 reads show a different base | Likely sequencing error (not a true variant) |

**Distinguishing true variants from errors:**  
A real germline variant typically appears in ~50% (het) or ~100% (hom alt) of reads consistently. Sequencing errors are random and appear in very few reads at a position — the variant caller uses Phred quality scores and read depth to assign confidence.

**Clicking on a read:**  
Click any read to open a detail pop-up showing the full read sequence, mapping quality, read name, flags, and all aligned bases.

---


## Structural Variants — What IGV Can Reveal


Beyond SNPs, IGV can reveal **structural variants** through the pattern of paired-end read orientations.

**Normal paired-end orientation:**  
In a standard insert-size library, paired reads point **inward** toward each other:
```
→ ...(insert)... ←
```

**Inversion signal:**  
If a region has been inverted on the sample's chromosome, reads spanning the inversion boundary will align in **same-direction orientation** (both pointing right or both pointing left):
```
→ ........→     (inversion signal)
```
This is because one read crosses the inversion breakpoint and its direction is flipped relative to the reference.

**Deletion signal:**  
Paired reads spanning a deletion in the sample will appear to span a larger-than-expected insert size (pairs are farther apart than normal).

**Duplication signal:**  
Paired reads spanning a tandem duplication will appear to span a smaller-than-expected insert, or reads may cluster with abnormally high coverage.

> IGV does not call structural variants automatically but visualises the raw read evidence. Dedicated tools (e.g. Delly, Lumpy, Manta) use this paired-end signal statistically to call SVs.

---

