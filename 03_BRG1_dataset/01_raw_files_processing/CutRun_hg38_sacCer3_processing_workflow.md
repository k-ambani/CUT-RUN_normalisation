# CUT&RUN Data Processing Workflow

Workflow used for processing CUT&RUN sequencing data using a human (hg38) and yeast (sacCer3) hybrid reference genome.

## Overview

This workflow processes paired-end CUT&RUN sequencing data through the following steps:

1. Directory Setup
2. Adapter trimming (BBDuk)
3. Alignment to hybrid genome (Bowtie2), BAM processing and filtering
4. Signal track generation (BigWig)
5. Peak calling (MACS2)

## Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.11.5 | Quality control |
| BBMap/BBDuk | 39.13 | Adapter trimming |
| Bowtie2 | 2.3.4.1 | Read alignment |
| SAMtools | 1.9 | BAM manipulation |
| Sambamba | 0.6.7 | BAM sorting/indexing |
| Picard | 3.0.0 | Duplicate marking |
| BEDTools | 2.27.1 | BED file operations |
| deepTools | 3.5.0 | BigWig generation |
| MACS2 | 2.2.7.1 | Peak calling |

## Reference Files Required

- **Hybrid genome index**: Combined hg38 + sacCer3 Bowtie2 index
- **Blacklist regions**: ENCODE hg38 blacklist (hg38-blacklist.v2.bed)
- **Adapter sequences**: TruSeq adapter reference (truseq.fa.gz)
---

## 1. Directory Setup

```bash
# Create project directory structure
PROJECT_DIR="/path/to/project"
mkdir -p ${PROJECT_DIR}/{raw_data,trim_data,results,scripts,workflow}
mkdir -p ${PROJECT_DIR}/results/{bam_files,bigwig,macs2,multiqc,fastq_screen}
mkdir -p ${PROJECT_DIR}/results/bam_files/{intermediate,final_bams}/{hg38,sacCer3}
```

---


## 2. Adapter Trimming (BBDuk)

Remove adapter sequences and low-quality bases from reads.

```bash
#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00

# Define paths
RAW_DATA_DIR="/path/to/raw/fastq/files"
TRIM_OUT="${PROJECT_DIR}/trim_data"
ADAPTER_REF="/path/to/bbmap/resources/truseq.fa.gz"

mkdir -p ${TRIM_OUT}

# Process each sample
for sample_dir in ${RAW_DATA_DIR}/*; do
    cd ${sample_dir}
    
    # Get R1 file and derive sample name
    R1_FILE=$(ls *R1_001.fastq.gz)
    SAMPLE=$(basename ${R1_FILE} _R1_001.fastq.gz)
    
    # Run BBDuk trimming
    # Parameters:
    #   ktrim=r: Trim adapters from right end
    #   k=23: Kmer length for matching
    #   mink=11: Minimum kmer length at read ends
    #   hdist=1: Hamming distance for adapter matching
    #   tpe/tbo: Trim paired-end reads to same length
    #   qtrim=r: Quality trim from right
    #   trimq=5: Quality threshold
    #   minlen=20: Minimum read length after trimming
    
    bbduk.sh \
        in1=${R1_FILE} \
        in2=${SAMPLE}_R2_001.fastq.gz \
        out1=${TRIM_OUT}/${SAMPLE}_R1.trim.fastq.gz \
        out2=${TRIM_OUT}/${SAMPLE}_R2.trim.fastq.gz \
        ref=${ADAPTER_REF} \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=5 minlen=20 \
        stats=${TRIM_OUT}/${SAMPLE}_bbduk.trimstats.txt \
        refstats=${TRIM_OUT}/${SAMPLE}_bbduk.refstats.txt
done

# QC on trimmed reads
mkdir -p ${TRIM_OUT}/fastqc_trimmed
fastqc -o ${TRIM_OUT}/fastqc_trimmed ${TRIM_OUT}/*.fastq.gz
cd ${TRIM_OUT}/fastqc_trimmed
multiqc --filename multiqc_trimmed.html *_fastqc.zip
```

---

## 3. Alignment to hybrid genome (Bowtie2), BAM processing and filtering

Align reads to combined hg38/sacCer3 reference for spike-in normalization.

```bash
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --mem=25G
#SBATCH --cpus-per-task=15
#SBATCH --time=10:00:00

set -xe

# Input files (passed as arguments)
FQ_R1=$1
FQ_R2=$2

# Define paths
TRIM_DIR="${PROJECT_DIR}/trim_data"
INTERMEDIATE_DIR="${PROJECT_DIR}/results/bam_files/intermediate"
FINAL_BAM_DIR="${PROJECT_DIR}/results/bam_files/final_bams"
GENOME_INDEX="/path/to/hybrid_genome/hg38_sacCer3_index"
BLACKLIST="/path/to/hg38-blacklist.v2.bed"

cd ${TRIM_DIR}

# Extract sample name
SAMPLE=$(basename ${FQ_R1} _R1.trim.fastq.gz)

# Create output directories
mkdir -p ${INTERMEDIATE_DIR}/{hg38,sacCer3}
mkdir -p ${FINAL_BAM_DIR}/{hg38,sacCer3}

#------------------------------------------------------------------------------
# ALIGNMENT
#------------------------------------------------------------------------------
# Bowtie2 parameters for CUT&RUN:
#   --dovetail: Allow dovetailing read pairs
#   --local: Local alignment (soft-clip ends)
#   --very-sensitive: Thorough search for alignments
#   --no-mixed: Suppress unpaired alignments for paired reads
#   --no-discordant: Suppress discordant alignments
#   -I 10 -X 700: Fragment size range (10-700 bp)

bowtie2 --dovetail --local --very-sensitive \
    --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 14 \
    -x ${GENOME_INDEX} \
    -1 ${FQ_R1} -2 ${FQ_R2} \
    -S ${INTERMEDIATE_DIR}/${SAMPLE}_unsorted.sam \
    2> ${INTERMEDIATE_DIR}/${SAMPLE}_bowtie2.txt

# Convert to sorted BAM
samtools view -hb -@ 14 ${INTERMEDIATE_DIR}/${SAMPLE}_unsorted.sam | \
    samtools sort -@ 14 -o ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam -
samtools index -@ 14 ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam

#------------------------------------------------------------------------------
# SPLIT BY GENOME (hg38)
#------------------------------------------------------------------------------
# Extract human-aligned reads (exclude sacCer3 chromosomes)
samtools idxstats ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam | \
    cut -f 1 | grep -v "sacCer3" | \
    xargs samtools view -@ 14 -hb ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam | \
    samtools sort -@ 14 -o ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38.bam -
samtools index -@ 14 ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38.bam

# Remove mitochondrial reads
samtools idxstats ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38.bam | \
    cut -f 1 | grep -v chrM | \
    xargs samtools view -@ 14 -hb ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38.bam \
    > ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_chrMrm.bam

# Remove blacklisted regions
bedtools intersect -v \
    -abam ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_chrMrm.bam \
    -b ${BLACKLIST} | \
    sambamba sort -t 14 -o ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_filtered.bam /dev/stdin

# Mark PCR duplicates
picard MarkDuplicates \
    INPUT=${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_filtered.bam \
    OUTPUT=${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_dupmarked.bam \
    METRICS_FILE=${FINAL_BAM_DIR}/hg38/${SAMPLE}_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads
# -q 2: Minimum mapping quality
# -F 0x04: Exclude unmapped reads
samtools view -q 2 -F 0x04 -b -@ 14 ${INTERMEDIATE_DIR}/hg38/${SAMPLE}_hg38_dupmarked.bam | \
    sambamba sort -t 14 -o ${FINAL_BAM_DIR}/hg38/${SAMPLE}_hg38_final.bam /dev/stdin

#------------------------------------------------------------------------------
# SPLIT BY GENOME (sacCer3 - Spike-in)
#------------------------------------------------------------------------------
# Extract yeast-aligned reads for spike-in normalisation
samtools idxstats ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam | \
    cut -f 1 | grep "sacCer3" | \
    xargs samtools view -@ 14 -hb ${INTERMEDIATE_DIR}/${SAMPLE}_sorted.bam | \
    samtools sort -@ 14 -o ${INTERMEDIATE_DIR}/sacCer3/${SAMPLE}_sacCer3.bam -
samtools index -@ 14 ${INTERMEDIATE_DIR}/sacCer3/${SAMPLE}_sacCer3.bam

# Mark and filter duplicates for spike-in
picard MarkDuplicates \
    INPUT=${INTERMEDIATE_DIR}/sacCer3/${SAMPLE}_sacCer3.bam \
    OUTPUT=${INTERMEDIATE_DIR}/sacCer3/${SAMPLE}_sacCer3_dupmarked.bam \
    METRICS_FILE=${FINAL_BAM_DIR}/sacCer3/${SAMPLE}_dup_metrics_sacCer3.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT

# Remove multi-mappers and unmapped reads
samtools view -q 2 -F 0x04 -b -@ 14 ${INTERMEDIATE_DIR}/sacCer3/${SAMPLE}_sacCer3_dupmarked.bam | \
    sambamba sort -t 14 -o ${FINAL_BAM_DIR}/sacCer3/${SAMPLE}_sacCer3_final.bam /dev/stdin

# Clean up intermediate SAM file
rm ${INTERMEDIATE_DIR}/${SAMPLE}_unsorted.sam
```

**Batch submission wrapper:**

```bash
#!/bin/bash
# Submit alignment jobs for all samples

TRIM_DIR="${PROJECT_DIR}/trim_data"
SCRIPT_DIR="${PROJECT_DIR}/scripts"

cd ${TRIM_DIR}

for R1_file in *_R1.trim.fastq.gz; do
    SAMPLE=$(basename ${R1_file} _R1.trim.fastq.gz)
    sbatch ${SCRIPT_DIR}/align.sh ${R1_file} ${SAMPLE}_R2.trim.fastq.gz
    sleep 1
done
```

---

## 4. BigWig Generation

Generate normalised signal tracks for visualisation.

```bash
#!/bin/bash
#SBATCH --job-name=bigwig
#SBATCH --mem=10G
#SBATCH --cpus-per-task=5
#SBATCH --time=4:00:00

# Input BAM file (passed as argument)
BAM_FILE=$1

# Define paths
INPUT_DIR="${PROJECT_DIR}/results/bam_files/final_bams/hg38"
OUTPUT_DIR="${PROJECT_DIR}/results/bigwig"

mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

SAMPLE=$(basename ${BAM_FILE} _hg38_final.bam)

# Generate RPKM-normalised BigWig
# Parameters:
#   --binSize 1: 1 bp resolution
#   --normalizeUsing RPKM: Reads per kilobase per million normalisation

bamCoverage \
    -b ${BAM_FILE} \
    -o ${OUTPUT_DIR}/${SAMPLE}.bw \
    --binSize 1 \
    --normalizeUsing RPKM \
    --numberOfProcessors 5
```

---

## 5. Peak Calling (MACS2)

### Without Input Control 

```bash
#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00

set -xe

# Define paths
INPUT_DIR="${PROJECT_DIR}/results/bam_files/final_bams/hg38"
OUTPUT_DIR="${PROJECT_DIR}/results/macs2/no_input"

mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

# Call peaks for all samples
# Parameters:
#   -f BAMPE: Paired-end BAM format
#   -g hs: Human genome size
#   -q 0.05: FDR threshold

for bam_file in *_hg38_final.bam; do
    SAMPLE=$(basename ${bam_file} _hg38_final.bam)
    
    macs2 callpeak \
        -t ${bam_file} \
        -f BAMPE \
        -g hs \
        -q 0.05 \
        -n ${SAMPLE} \
        --outdir ${OUTPUT_DIR}
done
```

### With IgG Input Control

```bash
#!/bin/bash
#SBATCH --job-name=macs2_input
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00

set -xe

# Define paths
INPUT_DIR="${PROJECT_DIR}/results/bam_files/final_bams/hg38"
OUTPUT_DIR="${PROJECT_DIR}/results/macs2/with_input"
IgG_CONTROL="${INPUT_DIR}/IgG_control_hg38_final.bam"

mkdir -p ${OUTPUT_DIR}
cd ${INPUT_DIR}

# Call peaks with IgG control
for bam_file in *_hg38_final.bam; do
    # Skip the IgG control file itself
    if [[ ${bam_file} == *"IgG"* ]]; then
        continue
    fi
    
    SAMPLE=$(basename ${bam_file} _hg38_final.bam)
    
    macs2 callpeak \
        -t ${bam_file} \
        -c ${IgG_CONTROL} \
        -f BAMPE \
        -g hs \
        -q 0.05 \
        -n ${SAMPLE} \
        --outdir ${OUTPUT_DIR}
done
```

---

## Output Files

| Stage | Output | Description |
|-------|--------|-------------|
| Trimming | `*_R1.trim.fastq.gz` | Trimmed forward reads |
| Trimming | `*_R2.trim.fastq.gz` | Trimmed reverse reads |
| Alignment | `*_hg38_final.bam` | Processed human alignments |
| Alignment | `*_sacCer3_final.bam` | Spike-in alignments |
| BigWig | `*.bw` | Signal tracks |
| Peaks | `*_peaks.narrowPeak` | Called peaks |
| Peaks | `*_summits.bed` | Peak summit positions |

---

## Notes

- **Spike-in normalisation**: The sacCer3 (yeast) alignments can be used to calculate spike-in normalisation factors for cross-sample comparisons.
  
- **Hybrid genome construction**: The hg38_sacCer3 hybrid index was created by concatenating the human (hg38) and yeast (sacCer3) reference genomes before indexing with Bowtie2.

- **SLURM configuration**: The `#SBATCH` directives shown are examples; adjust partition names and resource limits according to your HPC environment.
