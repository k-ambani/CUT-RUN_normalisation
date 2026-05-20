# CUT&RUN Data Processing Workflow

Workflow used for processing CUT&RUN sequencing data using a human (hg38) and
Drosophila (dm6) hybrid reference genome for spike-in normalization.

## Overview

This workflow processes paired-end CUT&RUN sequencing data through the following steps:

1. Quality control (FastQC/MultiQC)
2. Contamination screening (FastQ Screen)
3. Alignment to hybrid genome (Bowtie2)
4. BAM processing and filtering
5. Signal track generation (BigWig)
6. Peak calling (MACS2)

## Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.11.5 | Quality control |
| MultiQC | 1.8+ | Report aggregation |
| FastQ Screen | 0.15.3 | Contamination detection |
| Bowtie2 | 2.3.4.1 | Read alignment |
| SAMtools | 1.9 | BAM manipulation |
| Sambamba | 0.6.7 | BAM sorting/indexing |
| Picard | 3.0.0 | Duplicate marking |
| BEDTools | 2.27.1 | BED file operations |
| deepTools | 3.5.0 | BigWig generation |
| MACS2 | 2.2.7.1 | Peak calling |

## Reference Files Required

- **Hybrid genome index**: Combined hg38 + dm6 Bowtie2 index
- **Blacklist regions**: ENCODE hg38 blacklist (`hg38-blacklist.v2.bed`) and dm6 blacklist (`dm6-blacklist.v2.bed`)
- **FastQ Screen config**: Pre-configured contamination reference genomes

---

## 1. Directory Setup

```bash
project_dir="/path/to/project"
scratch_dir="/path/to/scratch"

mkdir -p $project_dir/{raw_data,results,scripts}
mkdir -p $project_dir/results/{bam_files,bw_files,macs2,multiqc,fastq_screen}
mkdir -p $project_dir/results/bam_files/final_bams/{hg38,dm6}
mkdir -p $scratch_dir/bam_files/intermediate_bams/{hg38,dm6}
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
#SBATCH --mail-user=user@institute.edu
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
multiqc *fastqc.zip
```

---

## 3. Contamination Screening (FastQ Screen)

Screen for potential contamination from other organisms.

**`fastq_screen_input.sh`** — single-sample script called by the batch wrapper.

```bash
#!/bin/bash
# =============================================================================
# Contamination Screening: FastQ Screen (single sample)
# Usage: sbatch fastq_screen_input.sh <sample.fastq.gz>
# =============================================================================

#SBATCH --job-name=fastq_screen
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=6
#SBATCH --time=10:00:00
#SBATCH --mail-user=user@institute.edu
#SBATCH --mail-type=ALL
#SBATCH --output=fastq_screen-%j.out
#SBATCH --error=fastq_screen-%j.err

# --- Paths -------------------------------------------------------------------

fq=$1
fastq_screen_dir="/path/to/FastQ-Screen"

# --- Run FastQ Screen --------------------------------------------------------

cd $fastq_screen_dir
./fastq_screen $fq
```

**`fastq_screen_batch.sh`** — submits one job per FASTQ file.

```bash
#!/bin/bash
# =============================================================================
# Contamination Screening: FastQ Screen batch submission
# Usage: bash fastq_screen_batch.sh
# =============================================================================

# --- Paths -------------------------------------------------------------------

project_dir="/path/to/project"
data_dir="/path/to/raw/fastq/files"
script_dir="$project_dir/scripts"

mkdir -p "$data_dir/fastq_screen/error_files"

# --- Submit per-sample jobs --------------------------------------------------

cd $data_dir

for infile in *.fastq.gz; do
  sbatch $script_dir/fastq_screen_input.sh $data_dir/$infile
  sleep 1
done
```

---

## 4. Alignment to Hybrid Genome (Bowtie2)

Align reads to the combined hg38/dm6 reference. Reads are then split by genome
and filtered to produce final BAMs for both the human (hg38) and spike-in (dm6)
fractions.

**`align_input.sh`** — single-sample script called by the batch wrapper.

```bash
#!/bin/bash
# =============================================================================
# Alignment: Bowtie2 to hybrid hg38/dm6 genome (single sample)
# Usage: sbatch align_input.sh <R1.fastq.gz> <R2.fastq.gz>
# =============================================================================

#SBATCH --job-name=align
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=15
#SBATCH --time=4:00:00
#SBATCH --mail-user=user@institute.edu
#SBATCH --mail-type=ALL
#SBATCH --output=align-%j.out
#SBATCH --error=align-%j.err

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

data_dir="/path/to/raw/fastq/files"
output_path="/path/to/scratch/bam_files/intermediate_bams"
output_finalbam_path="/path/to/project/results/bam_files/final_bams"
genome="/path/to/hybrid_genome/hg38_dm6_index"
bl_hg38="/path/to/hg38-blacklist.v2.bed"
bl_dm6="/path/to/dm6-blacklist.v2.bed"

# --- Sample name -------------------------------------------------------------

cd $data_dir

base=$(basename ${fq_R1} _1.fastq.gz)
base=${base:34}
echo "Sample: ${base}"
echo "Input files: ${fq_R1} ${fq_R2}"

# --- Output directories ------------------------------------------------------

mkdir -p $output_path/{hg38,dm6}
mkdir -p $output_finalbam_path/{hg38,dm6}

# --- File name variables -----------------------------------------------------

align_out=$output_path/${base}_hg38_dm6_unsorted.sam
align_bam_sorted=$output_path/${base}_hg38_dm6_sorted.bam

# hg38
align_bam_hg38=$output_path/hg38/${base}_hg38_nodm6.bam
align_bam_chrmrm_hg38=$output_path/hg38/${base}_hg38_nodm6_chrmrm.bam
align_bam_blrm_hg38=$output_path/hg38/${base}_hg38_nodm6_chrmrm_blrm.bam
align_bam_dupmarked_hg38=$output_path/hg38/${base}_hg38_nodm6_chrmrm_blrm_dupmarked.bam
align_bam_multimaprm_hg38=$output_finalbam_path/hg38/${base}_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam

# dm6
align_bam_dm6=$output_path/dm6/${base}_nohg38_dm6.bam
align_bam_chrmrm_dm6=$output_path/dm6/${base}_nohg38_dm6_chrmrm.bam
align_bam_blrm_dm6=$output_path/dm6/${base}_nohg38_dm6_chrmrm_blrm.bam
align_bam_dupmarked_dm6=$output_path/dm6/${base}_nohg38_dm6_chrmrm_blrm_dupmarked.bam
align_bam_multimaprm_dm6=$output_finalbam_path/dm6/${base}_nohg38_dm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam

set -ex

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
  -S $align_out &> $output_path/${base}_bowtie2.txt

samtools view -h -b -@ 14 $align_out | samtools sort -@ 14 -o $align_bam_sorted -
samtools index -@ 14 $align_bam_sorted

# --- Split by genome: hg38 ---------------------------------------------------

# Extract hg38 reads and strip "hg38_" chromosome prefix
samtools idxstats $align_bam_sorted | cut -f 1 | grep -v "dm6_" | \
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
  METRICS_FILE=$output_path/hg38/${base}_dup_metrics_hg38.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads (-q 2: MAPQ ≥ 2; -F 0x04: exclude unmapped)
samtools view -q 2 -F 0x04 -b -@ 14 $align_bam_dupmarked_hg38 \
  > $output_path/hg38/${base}_hg38_temp.bam
sambamba sort -t 14 -o $align_bam_multimaprm_hg38 $output_path/hg38/${base}_hg38_temp.bam
rm $output_path/hg38/${base}_hg38_temp.bam

# --- Split by genome: dm6 (spike-in) -----------------------------------------

# Extract dm6 reads and strip "dm6_" chromosome prefix
samtools idxstats $align_bam_sorted | cut -f 1 | grep -v "hg38_" | \
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
  METRICS_FILE=$output_path/dm6/${base}_dup_metrics_dm6.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads
samtools view -q 2 -F 0x04 -b -@ 14 $align_bam_dupmarked_dm6 \
  > $output_path/dm6/${base}_dm6_temp.bam
sambamba sort -t 14 -o $align_bam_multimaprm_dm6 $output_path/dm6/${base}_dm6_temp.bam
rm $output_path/dm6/${base}_dm6_temp.bam

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

data_dir="/path/to/raw/fastq/files"
script_dir="/path/to/project/scripts"

# --- Submit per-sample jobs --------------------------------------------------

cd $data_dir

for infile in *_1.fastq.gz; do
  base=$(basename ${infile} _1.fastq.gz)
  sbatch $script_dir/align_input.sh ${infile} ${base}_2.fastq.gz
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

module load deeptools/3.5.0

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

Peaks are called without a negative control using an FDR threshold of q < 0.01.

**`macs2_q0.01.sh`**

```bash
#!/bin/bash
# =============================================================================
# Peak Calling: MACS2 without input control (q < 0.01)
# Usage: sbatch macs2_q0.01.sh
# =============================================================================

#SBATCH --job-name=macs2
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --mail-user=user@institute.edu
#SBATCH --mail-type=ALL
#SBATCH --output=macs2-%j.out
#SBATCH --error=macs2-%j.err

# --- Modules -----------------------------------------------------------------

module load macs/2.2.7.1
module load bedtools/2.27.1

# --- Paths -------------------------------------------------------------------

input_dir="/path/to/project/results/bam_files/final_bams/hg38"
output_dir="/path/to/project/results/macs2/macs2_q0.01/with_no_neg_control"

mkdir -p $output_dir

# --- Call peaks --------------------------------------------------------------

set -xe

cd $input_dir

bam_suffix="_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam"

for bam in *.bam; do
  base=$(basename ${bam} ${bam_suffix})
  echo "Sample: $base"

  macs2 callpeak \
    -t $bam \
    -f BAMPE \
    -g hs \
    -q 0.01 \
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
| Alignment | `*_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam` | Processed human alignments |
| Alignment | `*_nohg38_dm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam` | Spike-in (dm6) alignments |
| Alignment | `*_dup_metrics_hg38.txt` | Picard duplicate metrics (hg38) |
| Alignment | `*_dup_metrics_dm6.txt` | Picard duplicate metrics (dm6) |
| BigWig | `*.bw` | RPKM-normalized signal tracks |
| Peaks | `*_peaks.narrowPeak` | Called peaks |
| Peaks | `*_summits.bed` | Peak summit positions |

---

## Notes

- **Spike-in normalization**: The dm6 alignments are used to calculate spike-in
  scaling factors for cross-sample normalization. See downstream analysis scripts
  for scale factor computation.

- **Hybrid genome construction**: The hg38/dm6 hybrid index was created by
  concatenating the human (hg38) and Drosophila (dm6) reference genomes before
  indexing with Bowtie2. Chromosome names are prefixed with `hg38_` and `dm6_`
  respectively to allow genome-of-origin splitting post-alignment; these prefixes
  are stripped after splitting.

- **BAM filtering steps**: The suffix
  `_chrmrm_blrm_duprm_unmappedrm_multimaprm` encodes the sequential filtering
  applied: mitochondrial reads removed (chrmrm), blacklist regions removed
  (blrm), duplicates marked (duprm), unmapped reads removed (unmappedrm),
  multi-mappers removed (multimaprm).

- **Resource allocation**: Memory and CPU requirements may need adjustment based
  on sequencing depth and available compute resources.

- **SLURM configuration**: The `#SBATCH` directives shown are examples; adjust
  partition names and resource limits according to your HPC environment.
