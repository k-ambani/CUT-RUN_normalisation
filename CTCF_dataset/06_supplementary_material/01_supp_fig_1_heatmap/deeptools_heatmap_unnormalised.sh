#!/bin/bash
#SBATCH --job-name=unnormalised_bigwigs_deeptools
#SBATCH -t 2:00:00
#SBATCH --output=DT_bigwigs_heatmap-%j.out
#SBATCH --error=DT_bigwigs_heatmap-%j.err


# script to make unnormalised bigwigs and deeptools heatmap
# saved as deeptools_heatmap_unnormalised.sh

module load deepTools/3.5.2

set -xe 
#!/bin/bash
# ==============================================================================
# generate_raw_heatmap.sh
#
# Generates UNNORMALISED bigwig tracks from BAM files, then produces a
# deepTools heatmap showing raw signal (baseline before any normalisation).
# ==============================================================================

# ------------------------------------------------------------------------
# CONFIG - edit these paths/settings for your data
# ------------------------------------------------------------------------

# Directory containing indexed BAM files (must have matching .bam.bai)
BAM_DIR="/path/to/project/results/bam_files/final_bams/hg38"

# Output directory for bigwigs and heatmap
OUT_DIR="/path/to/project/results/bw_files/raw_signal"

# BED file of peak regions to plot (e.g. consensus peak set)
PEAKS_BED="/path/to/project/results/bedfile/consensus_peakset.bed"

HEATMAP_DIR="/path/to/project/results/Deeptools/output"
# Bin size for bigwig generation (bp)
BIN_SIZE=1

# Number of threads
THREADS=8

# Window around reference point (bp)
UPSTREAM=1000
DOWNSTREAM=1000


# ------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------


BIGWIGS=()
LABELS=()

# ------------------------------------------------------------------------
# Step 1: Generate unnormalised bigwigs
# ------------------------------------------------------------------------
echo "=== Generating unnormalised bigwigs ==="
cd $BAM_DIR
# 
for bam in *bam; do
   echo $bam
   sample_name=$(basename "${bam}" _Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
   sample_name=$(echo "$sample_name" | sed -E 's/20160.*HsDm_CTCF_//; s/CTCF_//')
   bw_out="${OUT_DIR}/${sample_name}_raw.bw"
 
   echo "Processing ${sample_name}..."
 
   bamCoverage \
     --bam "${BAM_DIR}/${bam}" \
     --outFileName "${bw_out}" \
     --normalizeUsing None \
     --binSize "${BIN_SIZE}" \
     --numberOfProcessors "${THREADS}" \
     --extendReads

   BIGWIGS+=("${bw_out}")
   LABELS+=("${sample_name}")
   echo $BIGWIGS
   echo $LABELS
 done
 
echo "Done. Bigwigs written to ${OUT_DIR}/bigwigs/"
 
# ------------------------------------------------------------------------
# Step 2: Compute matrix over peak regions
# ------------------------------------------------------------------------

echo "=== Computing matrix ==="

MATRIX_OUT="${OUT_DIR}/matrix_raw.gz"
# 
computeMatrix reference-point \
   --scoreFileName "${BIGWIGS[@]}" \
   --regionsFileName "${PEAKS_BED}" \
   --referencePoint center \
   --upstream "${UPSTREAM}" \
   --downstream "${DOWNSTREAM}" \
   --binSize "${BIN_SIZE}" \
   --numberOfProcessors "${THREADS}" \
   --samplesLabel "${LABELS[@]}" \
   --outFileName "${MATRIX_OUT}" \
   --missingDataAsZero

echo "Matrix written to ${MATRIX_OUT}"
 
# ------------------------------------------------------------------------
# Step 3: Plot heatmap
# ------------------------------------------------------------------------
 
echo "=== Plotting heatmap ==="

HEATMAP_OUT="${HEATMAP_DIR}/heatmap_raw_signal.pdf"

plotHeatmap \
  --matrixFile "${MATRIX_OUT}" \
  --outFileName "${HEATMAP_OUT}" \
  --colorMap Greens \
  --heatmapHeight 12 \
  --heatmapWidth 4 \
  --xAxisLabel "distance from peak center (bp)" \
  --refPointLabel "center" \
  --plotTitle "Raw (unnormalised) signal"

echo "Heatmap written to ${HEATMAP_OUT}"
echo "=== Done ==="