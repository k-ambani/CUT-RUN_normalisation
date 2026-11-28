# CUT&RUN Data Processing Workflow

Workflow used for processing CUT&RUN sequencing data using a "tribrid" reference
genome combining human (hg38), yeast (sacCer3), and Drosophila (dm6) sequences.
Both sacCer3 and dm6 serve as spike-in controls for normalisation.

## Overview

This workflow processes paired-end CUT&RUN sequencing data through the following steps:

1. Quality control (FastQC/MultiQC)
2. Adapter trimming (BBDuk)
3. Alignment to tribrid genome (Bowtie2)
4. BAM processing and filtering
5. Signal track generation (BigWig)
6. Peak calling (MACS2)

## Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.11.5 | Quality control |
| MultiQC | 1.8+ | Report aggregation |
| BBMap/BBDuk | — | Adapter trimming |
| Bowtie2 | 2.3.4.1 | Read alignment |
| SAMtools | 1.9 | BAM manipulation |
| Sambamba | 0.6.7 | BAM sorting/indexing |
| Picard | 3.0.0 | Duplicate marking |
| BEDTools | 2.27.1 | BED file operations |
| deepTools | 3.5.2 | BigWig generation |
| MACS2 | 2.2.7.1 | Peak calling |

## Reference Files Required

- **Tribrid genome index**: Combined hg38 + sacCer3 + dm6 Bowtie2 index
- **Blacklist regions**: ENCODE hg38 blacklist (`hg38-blacklist.v2.bed`) and dm6 blacklist (`dm6-blacklist.v2.bed`)
- **Adapter sequences**: TruSeq adapter reference (`truseq.fa.gz`)
- **FastQ Screen config**: Pre-configured contamination reference genomes
- **IgG control BAM**: Rabbit IgG control for MACS2 peak calling

---

## 1. Directory Setup

```bash
project_dir="/path/to/project"
scratch_dir="/path/to/scratch"

mkdir -p $project_dir/{raw_data,trim_data,results,scripts}
mkdir -p $project_dir/results/{bam_files,bw_files,macs2,fastq_screen}
mkdir -p $project_dir/results/bam_files/final_bams/{hg38,sacCer3,dm6}
mkdir -p $project_dir/results/bam_files/final_bams/duplicates/{hg38,sacCer3,dm6}
mkdir -p $project_dir/results/bam_files/final_bams/bowtie2_out
mkdir -p $scratch_dir/bam_files/intermediate_bams/{hg38,sacCer3,dm6}
```

---

## 2. Quality Control (FastQC)

Assess raw sequencing data quality before processing.

**`fastqc.sh`**

```bash
#!/bin/bash
# =============================================================================
# Quality Control: FastQC and MultiQC
# Usage: sbatch fastqc.sh
# =============================================================================

#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --cpus-per-task=10
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=fastqc-%j.out
#SBATCH --error=fastqc-%j.err

# --- Modules -----------------------------------------------------------------

module load fastqc/0.11.5
module load multiqc/1.8

# --- Paths -------------------------------------------------------------------

raw_data_dir="/path/to/raw/fastq/files"
fastqc_output_dir="$raw_data_dir/fastqc_out"

# --- Run FastQC --------------------------------------------------------------

set -xe

mkdir -p $fastqc_output_dir

cd $raw_data_dir
fastqc -t 9 -o $fastqc_output_dir *.fastq.gz

# --- Aggregate reports with MultiQC ------------------------------------------

cd $fastqc_output_dir
multiqc --filename multiqc_report.html *fastqc.zip
```

---

## 3. Adapter Trimming (BBDuk)

Remove adapter sequences and low-quality bases from reads.

**`trim.sh`**

```bash
#!/bin/bash
# =============================================================================
# Adapter Trimming: BBDuk
# Usage: sbatch trim.sh
# =============================================================================

#SBATCH --job-name=trim
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=trim-%j.out
#SBATCH --error=trim-%j.err

# --- Modules -----------------------------------------------------------------

module load fastqc/0.11.5
module load multiqc/1.8
module load java/1.8.0_312-jdk

# --- Paths -------------------------------------------------------------------

rawdata="/path/to/raw/fastq/files"
trim_output="/path/to/project/trim_data"
bbduk="/path/to/bbmap/bbduk.sh"
adapter_ref="/path/to/bbmap/resources/truseq.fa.gz"

# --- Trim reads --------------------------------------------------------------

set -xe

mkdir -p $trim_output
cd $rawdata

# Parameters:
#   ktrim=r: trim adapters from right end
#   k=23/mink=11/hdist=1: kmer matching settings
#   tpe/tbo: trim paired-end reads consistently
#   qtrim=r trimq=5: quality trim from right, threshold q=5
#   minlen=20: discard reads shorter than 20 bp after trimming

for trimfastq in *_R1_001_combined.fastq.gz; do
  base=$(basename ${trimfastq} _R1_001_combined.fastq.gz)
  echo "Sample: ${base}"

  $bbduk \
    in1=${trimfastq} \
    in2=${base}_R2_001_combined.fastq.gz \
    out1=$trim_output/${base}_R1.trim.fastq.gz \
    out2=$trim_output/${base}_R2.trim.fastq.gz \
    ref=$adapter_ref \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=5 minlen=20 \
    stats=$trim_output/${base}_bbduk.trimstats.txt \
    refstats=$trim_output/${base}_bbduk.refstats.txt \
    &> $trim_output/${base}_bbduk.stdout_stats.txt
done

# --- QC on trimmed reads -----------------------------------------------------

mkdir -p $trim_output/fastqc_out_trim
cd $trim_output

fastqc -o fastqc_out_trim *fastq.gz

cd $trim_output/fastqc_out_trim
multiqc --filename multiqc_trimmed.html *fastqc.zip
```

---

## 4. Alignment to Tribrid Genome (Bowtie2)

Align reads to the combined hg38/sacCer3/dm6 reference. Reads are then split by
genome and filtered to produce final BAMs for the human (hg38), yeast (sacCer3),
and Drosophila (dm6) fractions.

Chromosome names in the hybrid index carry species prefixes (`hg38_`, `sacCer3_`,
`dm6_`); these are stripped after splitting to restore standard chromosome names.

**`align_input.sh`** — single-sample script called by the batch wrapper.

```bash
#!/bin/bash
# =============================================================================
# Alignment: Bowtie2 to tribrid hg38/sacCer3/dm6 genome (single sample)
# Usage: sbatch align_input.sh <R1.trim.fastq.gz> <R2.trim.fastq.gz>
# =============================================================================

#SBATCH --job-name=align
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=15
#SBATCH --time=7:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=align-%j.out
#SBATCH --error=align-%j.err

# --- Setup -------------------------------------------------------------------

set -xe
source /etc/profile.d/modules.sh

# --- Modules -----------------------------------------------------------------

module load bowtie2/2.3.4.1
module load samtools/1.9
module load deeptools/3.5.0
module load sambamba/0.6.7
module load picard/3.0.0
module load bedtools/2.27.1

# --- Paths -------------------------------------------------------------------

fq_R1=$1
fq_R2=$2

data_dir="/path/to/project/trim_data"
output_path="/path/to/scratch/bam_files/intermediate_bams"
output_finalbam_path="/path/to/project/results/bam_files/final_bams"
genome="/path/to/tribrid_genome/hg38_sacCer3_dm6_index"
bl_hg38="/path/to/hg38-blacklist.v2.bed"
bl_dm6="/path/to/dm6-blacklist.v2.bed"

# --- Sample name -------------------------------------------------------------

cd $data_dir

base=$(basename ${fq_R1} _R1.trim.fastq.gz)
echo "Sample: ${base}"
echo "Input files: ${fq_R1} ${fq_R2}"

# --- Output directories ------------------------------------------------------

mkdir -p $output_path/{hg38,sacCer3,dm6}
mkdir -p $output_finalbam_path/{hg38,sacCer3,dm6}
mkdir -p $output_finalbam_path/bowtie2_out
mkdir -p $output_finalbam_path/duplicates/{hg38,sacCer3,dm6}

# --- File name variables -----------------------------------------------------

align_out=$output_path/${base}_allgenomes_unsorted.sam
align_bam_sorted=$output_path/${base}_allgenomes_sorted.bam

# hg38
align_bam_hg38=$output_path/hg38/${base}_onlyhg38.bam
align_bam_chrmrm_hg38=$output_path/hg38/${base}_onlyhg38_chrmrm.bam
align_bam_blrm_hg38=$output_path/hg38/${base}_onlyhg38_chrmrm_blrm.bam
align_bam_dupmarked_hg38=$output_path/hg38/${base}_onlyhg38_chrmrm_blrm_dupmarked.bam
align_bam_multimaprm_hg38=$output_finalbam_path/hg38/${base}_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam

# dm6
align_bam_dm6=$output_path/dm6/${base}_onlydm6.bam
align_bam_chrmrm_dm6=$output_path/dm6/${base}_onlydm6_chrmrm.bam
align_bam_blrm_dm6=$output_path/dm6/${base}_onlydm6_chrmrm_blrm.bam
align_bam_dupmarked_dm6=$output_path/dm6/${base}_onlydm6_chrmrm_blrm_dupmarked.bam
align_bam_multimaprm_dm6=$output_finalbam_path/dm6/${base}_onlydm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam

# sacCer3
align_bam_sacCer3=$output_path/sacCer3/${base}_onlysacCer3.bam
align_bam_dupmarked_sacCer3=$output_path/sacCer3/${base}_onlysacCer3_chrmrm_blrm_dupmarked.bam
align_bam_multimaprm_sacCer3=$output_finalbam_path/sacCer3/${base}_onlysacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam

# --- Alignment ---------------------------------------------------------------
# Parameters:
#   --dovetail/--local/--very-sensitive: recommended CUT&RUN settings
#   --no-mixed/--no-discordant: suppress unpaired and discordant alignments
#   -I 10 -X 700: fragment size range (10-700 bp)

bowtie2 --dovetail --local --very-sensitive \
  --no-mixed --no-discordant --phred33 \
  -I 10 -X 700 -p 14 \
  -x $genome \
  -1 $fq_R1 -2 $fq_R2 \
  -S $align_out &> $output_finalbam_path/bowtie2_out/${base}_bowtie2.txt

samtools view -h -b -@ 14 $align_out | samtools sort -@ 14 -o $align_bam_sorted -
samtools index -@ 14 $align_bam_sorted

# --- Split by genome: hg38 ---------------------------------------------------

# Extract hg38 reads and strip "hg38_" chromosome prefix
samtools idxstats $align_bam_sorted | cut -f 1 | grep "hg38_" | \
  xargs samtools view -@ 14 -hb $align_bam_sorted | \
  samtools sort -@ 14 -o $output_path/hg38/${base}_hg38_temp.bam -
samtools view -h -@ 14 $output_path/hg38/${base}_hg38_temp.bam | \
  sed 's/hg38_//' | samtools view -b -@ 14 -o $align_bam_hg38 -
samtools index -@ 14 $align_bam_hg38
rm $output_path/hg38/${base}_hg38_temp.bam

# Remove mitochondrial reads
samtools idxstats $align_bam_hg38 | cut -f 1 | grep -v chrM | \
  xargs samtools view -@ 14 -hb $align_bam_hg38 > $align_bam_chrmrm_hg38

# Remove blacklisted regions
bedtools intersect -v -abam $align_bam_chrmrm_hg38 -b $bl_hg38 \
  > $output_path/hg38/${base}_hg38_temp.bam
sambamba sort -t 14 -o $align_bam_blrm_hg38 $output_path/hg38/${base}_hg38_temp.bam
rm $output_path/hg38/${base}_hg38_temp.bam

# Mark PCR duplicates
java -Xmx8g -jar /path/to/picard.jar MarkDuplicates \
  INPUT=$align_bam_blrm_hg38 \
  OUTPUT=$align_bam_dupmarked_hg38 \
  METRICS_FILE=$output_finalbam_path/duplicates/hg38/${base}_dup_metrics_hg38.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads (-q 2: MAPQ ≥ 2; -F 0x04: exclude unmapped)
samtools view -q 2 -F 0x04 -b -@ 14 $align_bam_dupmarked_hg38 \
  > $output_path/hg38/${base}_hg38_temp.bam
sambamba sort -t 14 -o $align_bam_multimaprm_hg38 $output_path/hg38/${base}_hg38_temp.bam
rm $output_path/hg38/${base}_hg38_temp.bam

# --- Split by genome: dm6 ----------------------------------------------------

# Extract dm6 reads and strip "dm6_" chromosome prefix
samtools idxstats $align_bam_sorted | cut -f 1 | grep "dm6_" | \
  xargs samtools view -@ 14 -hb $align_bam_sorted | \
  samtools sort -@ 14 -o $output_path/dm6/${base}_dm6_temp.bam -
samtools view -h -@ 14 $output_path/dm6/${base}_dm6_temp.bam | \
  sed 's/dm6_//' | samtools view -b -@ 14 -o $align_bam_dm6 -
samtools index -@ 14 $align_bam_dm6
rm $output_path/dm6/${base}_dm6_temp.bam

# Remove mitochondrial reads
samtools idxstats $align_bam_dm6 | cut -f 1 | grep -v chrM | \
  xargs samtools view -@ 14 -hb $align_bam_dm6 > $align_bam_chrmrm_dm6

# Remove blacklisted regions
bedtools intersect -v -abam $align_bam_chrmrm_dm6 -b $bl_dm6 \
  > $output_path/dm6/${base}_dm6_temp.bam
sambamba sort -t 14 -o $align_bam_blrm_dm6 $output_path/dm6/${base}_dm6_temp.bam
rm $output_path/dm6/${base}_dm6_temp.bam

# Mark PCR duplicates
java -Xmx8g -jar /path/to/picard.jar MarkDuplicates \
  INPUT=$align_bam_blrm_dm6 \
  OUTPUT=$align_bam_dupmarked_dm6 \
  METRICS_FILE=$output_finalbam_path/duplicates/dm6/${base}_dup_metrics_dm6.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads
samtools view -q 2 -F 0x04 -b -@ 14 $align_bam_dupmarked_dm6 \
  > $output_path/dm6/${base}_dm6_temp.bam
sambamba sort -t 14 -o $align_bam_multimaprm_dm6 $output_path/dm6/${base}_dm6_temp.bam
rm $output_path/dm6/${base}_dm6_temp.bam

# --- Split by genome: sacCer3 (spike-in) -------------------------------------

# Extract sacCer3 reads and strip "sacCer3_" chromosome prefix
samtools idxstats $align_bam_sorted | cut -f 1 | grep "sacCer3_" | \
  xargs samtools view -@ 14 -hb $align_bam_sorted | \
  samtools sort -@ 14 -o $output_path/sacCer3/${base}_sacCer3_temp.bam -
samtools view -h -@ 14 $output_path/sacCer3/${base}_sacCer3_temp.bam | \
  sed 's/sacCer3_//' | samtools view -b -@ 14 -o $align_bam_sacCer3 -
samtools index -@ 14 $align_bam_sacCer3
rm $output_path/sacCer3/${base}_sacCer3_temp.bam

# Mark PCR duplicates
java -Xmx8g -jar /path/to/picard.jar MarkDuplicates \
  INPUT=$align_bam_sacCer3 \
  OUTPUT=$align_bam_dupmarked_sacCer3 \
  METRICS_FILE=$output_finalbam_path/duplicates/sacCer3/${base}_dup_metrics_sacCer3.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads
samtools view -q 2 -F 0x04 -b -@ 14 $align_bam_dupmarked_sacCer3 \
  > $output_path/sacCer3/${base}_sacCer3_temp.bam
sambamba sort -t 14 -o $align_bam_multimaprm_sacCer3 \
  $output_path/sacCer3/${base}_sacCer3_temp.bam
rm $output_path/sacCer3/${base}_sacCer3_temp.bam

# --- Cleanup -----------------------------------------------------------------

rm $align_out
```

**`align_batch.sh`** — submits one job per sample pair.

```bash
#!/bin/bash
# =============================================================================
# Alignment: Bowtie2 batch submission
# Usage: bash align_batch.sh
# =============================================================================

# --- Paths -------------------------------------------------------------------

data_dir="/path/to/project/trim_data"
script_dir="/path/to/project/scripts"

# --- Submit per-sample jobs --------------------------------------------------

cd $data_dir

for infile in *_R1.trim.fastq.gz; do
  base=$(basename ${infile} _R1.trim.fastq.gz)
  sbatch $script_dir/align_input.sh ${infile} ${base}_R2.trim.fastq.gz
  sleep 1
done
```

---

## 5. BigWig Generation

Generate RPKM-normalized signal tracks from the final hg38 BAMs for visualization.

**`bigwig_RPKM_input.sh`** — single-sample script called by the batch wrapper.

```bash
#!/bin/bash
# =============================================================================
# BigWig Generation: RPKM normalization (single sample)
# Usage: sbatch bigwig_RPKM_input.sh <sample.bam>
# =============================================================================

#SBATCH --job-name=bigwig_RPKM
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5
#SBATCH --time=4:00:00
#SBATCH --output=bigwig_RPKM-%j.out
#SBATCH --error=bigwig_RPKM-%j.err

# --- Modules -----------------------------------------------------------------

module load deepTools/3.5.2

# --- Paths -------------------------------------------------------------------

bam_file=$1

input_dir="/path/to/project/results/bam_files/final_bams/hg38"
output_dir="/path/to/project/results/bw_files/RPKM"

# --- Generate RPKM-normalized BigWig -----------------------------------------

cd $input_dir

base=$(basename ${bam_file} _chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
echo "Sample: ${base}"

bamCoverage \
  -b $bam_file \
  -o $output_dir/${base}.bw \
  --binSize 1 \
  --normalizeUsing RPKM \
  --extendReads \
  --numberOfProcessors 5
```

**`bigwig_RPKM_batch.sh`** — submits one job per BAM file.

```bash
#!/bin/bash
# =============================================================================
# BigWig Generation: RPKM normalization batch submission
# Usage: bash bigwig_RPKM_batch.sh
# =============================================================================

# --- Paths -------------------------------------------------------------------

input_dir="/path/to/project/results/bam_files/final_bams/hg38"
script_dir="/path/to/project/scripts"
output_dir="/path/to/project/results/bw_files/RPKM/error_files"

mkdir -p $output_dir

# --- Submit per-sample jobs --------------------------------------------------

cd $input_dir

for infile in *_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam; do
  sbatch $script_dir/bigwig_RPKM_input.sh ${infile}
  sleep 1
done
```

---

## 6. Peak Calling (MACS2)

Peaks are called with an IgG negative control using an FDR threshold of q < 0.05.

**`macs2_withinput.sh`**

```bash
#!/bin/bash
# =============================================================================
# Peak Calling: MACS2 with IgG input control (q < 0.05)
# Usage: sbatch macs2_withinput.sh
# =============================================================================

#SBATCH --job-name=macs2
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=macs2-%j.out
#SBATCH --error=macs2-%j.err
#SBATCH --container=/config/spack/containers/centos7/container.sif

# --- Setup -------------------------------------------------------------------

set -xe
source /etc/profile.d/modules.sh

# --- Modules -----------------------------------------------------------------

module load macs/2.2.7.1
module load bedtools/2.27.1

# --- Paths -------------------------------------------------------------------

input_dir="/path/to/project/results/bam_files/final_bams/hg38"
igg_control="/path/to/project/results/bam_files/final_bams/hg38/IgG_control_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam"
output_dir="/path/to/project/results/macs2/withinput/macs2_q0.05"

mkdir -p $output_dir

# --- Call peaks --------------------------------------------------------------

cd $input_dir

for bam in *_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam; do
  # Skip IgG control
  if [[ $bam == *"IgG"* ]]; then
    continue
  fi

  base=$(basename ${bam} _chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  echo "Sample: $base"

  macs2 callpeak \
    -t $bam \
    -c $igg_control \
    -f BAMPE \
    -g hs \
    -q 0.05 \
    -n $base \
    --outdir $output_dir
done
```

---

## Output Files

| Stage | Output | Description |
|-------|--------|-------------|
| QC | `*_fastqc.html` | Per-sample quality reports |
| QC | `multiqc_*.html` | Aggregated quality reports |
| Trimming | `*_R1.trim.fastq.gz` | Trimmed forward reads |
| Trimming | `*_R2.trim.fastq.gz` | Trimmed reverse reads |
| Trimming | `*_bbduk.trimstats.txt` | BBDuk trimming statistics |
| Alignment | `*_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam` | Processed human alignments |
| Alignment | `*_onlydm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam` | Drosophila spike-in alignments |
| Alignment | `*_onlysacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam` | Yeast spike-in alignments |
| Alignment | `*_dup_metrics_*.txt` | Picard duplicate metrics |
| Alignment | `*_bowtie2.txt` | Bowtie2 alignment statistics |
| BigWig | `*.bw` | RPKM-normalized signal tracks |
| Peaks | `*_peaks.narrowPeak` | Called peaks |
| Peaks | `*_summits.bed` | Peak summit positions |

---

## Notes

- **Tribrid genome construction**: The hg38/sacCer3/dm6 hybrid index was created
  by concatenating all three reference genomes before indexing with Bowtie2.
  Chromosome names are prefixed with `hg38_`, `sacCer3_`, and `dm6_` to allow
  genome-of-origin splitting post-alignment; these prefixes are stripped after
  splitting.

- **Spike-in normalization**: Both sacCer3 (yeast) and dm6 (Drosophila) alignments
  can be used to compute spike-in scaling factors for cross-sample normalization.
  See downstream analysis scripts for scale factor computation.

- **BAM filtering steps**: The suffix `_chrmrm_blrm_duprm_unmappedrm_multimaprm`
  encodes the sequential filtering applied: mitochondrial reads removed (chrmrm),
  blacklist regions removed (blrm), duplicates marked (duprm), unmapped reads
  removed (unmappedrm), multi-mappers removed (multimaprm). Note that sacCer3
  BAMs skip the blacklist removal step.

- **IgG control**: The MACS2 peak calling step uses a matched IgG negative control.
  Update the `igg_control` path to point to the correct IgG BAM for your experiment.

- **SLURM configuration**: The `#SBATCH` directives shown are examples; adjust
  partition names and resource limits according to your HPC environment.
