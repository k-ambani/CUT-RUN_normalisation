#!/bin/bash
#SBATCH --job-name=unnormalised_bigwigs_deeptools
#SBATCH -t 10:00:00
#SBATCH --output=DT_bigwigs_heatmap-%j.out
#SBATCH --error=DT_bigwigs_heatmap-%j.err

module load deepTools/3.5.2
set -xe

# ==============================================================================
# generate_raw_heatmap_per_comparison.sh
# ==============================================================================

OUT_DIR="/path/to/project/results/bw_files/raw_signal"
PEAKS_BED_DIR="/path/to/project/R_analysis_reviewerresponse/4.Unique_diff_peaks_analysis/peaksets"
HEATMAP_DIR="/path/to/project/R_analysis_reviewerresponse/4.Unique_diff_peaks_analysis/deeptools_of_unique_peaksets/output"

BIN_SIZE=1
THREADS=8
UPSTREAM=1000
DOWNSTREAM=1000

mkdir -p "${HEATMAP_DIR}"

# ------------------------------------------------------------------------
# Group -> bigwig(s) and label(s). Keys must match the tokens used in the
# BED filenames (E2, noE2, sgControl, sgER).
# ------------------------------------------------------------------------

declare -A GROUP_BW
GROUP_BW[E2]="${OUT_DIR}/CSS-E2-R1-ER_raw.bw ${OUT_DIR}/CSS-E2-R2-ER_raw.bw"
GROUP_BW[noE2]="${OUT_DIR}/CSS-R1-ER_raw.bw ${OUT_DIR}/CSS-R2-ER_raw.bw"
GROUP_BW[sgControl]="${OUT_DIR}/Complete-R1-ER_raw.bw ${OUT_DIR}/Complete-R2-ER_raw.bw"
GROUP_BW[sgER]="${OUT_DIR}/sgER-Complete-R1-ER_raw.bw ${OUT_DIR}/sgER-Complete-R2-ER_raw.bw"

declare -A GROUP_LABELS
GROUP_LABELS[E2]="CSS-E2-Rep1 CSS-E2-Rep2"
GROUP_LABELS[noE2]="CSS-noE2-Rep1 CSS-noE2-Rep2"
GROUP_LABELS[sgControl]="sgControl-Rep1 sgControl-Rep2"
GROUP_LABELS[sgER]="sgER-Rep1 sgER-Rep2"

# ------------------------------------------------------------------------
# Loop over each peak BED file
# ------------------------------------------------------------------------

for BED_FILE in "${PEAKS_BED_DIR}"/*_bg_rpkm_not_dm6_sacCer3.bed; do

  [ -e "${BED_FILE}" ] || { echo "No BED files found in ${PEAKS_BED_DIR}"; exit 1; }

  BASENAME=$(basename "${BED_FILE}")
  COMPARISON_NAME="${BASENAME%_bg_rpkm_not_dm6_sacCer3.bed}"

  IFS='_' read -r GROUP1 GROUP2 <<< "${COMPARISON_NAME}"

  if [ -z "${GROUP_BW[${GROUP1}]:-}" ] || [ -z "${GROUP_BW[${GROUP2}]:-}" ]; then
    echo "WARNING: could not map groups for '${COMPARISON_NAME}' (parsed as '${GROUP1}' / '${GROUP2}') — skipping ${BED_FILE}"
    continue
  fi

  BW_FILES="${GROUP_BW[${GROUP1}]} ${GROUP_BW[${GROUP2}]}"
  LABELS="${GROUP_LABELS[${GROUP1}]} ${GROUP_LABELS[${GROUP2}]}"

  echo "=== Comparison: ${COMPARISON_NAME} (${GROUP1} vs ${GROUP2}) ==="
  echo "BED:    ${BED_FILE}"
  echo "BWs:    ${BW_FILES}"
  echo "Labels: ${LABELS}"

  MATRIX_OUT="${OUT_DIR}/matrix_raw_${COMPARISON_NAME}.gz"
  HEATMAP_OUT="${HEATMAP_DIR}/heatmap_raw_signal_${COMPARISON_NAME}.pdf"

  computeMatrix reference-point \
    --scoreFileName ${BW_FILES} \
    --regionsFileName "${BED_FILE}" \
    --referencePoint center \
    --upstream "${UPSTREAM}" \
    --downstream "${DOWNSTREAM}" \
    --binSize "${BIN_SIZE}" \
    --numberOfProcessors "${THREADS}" \
    --samplesLabel ${LABELS} \
    --outFileName "${MATRIX_OUT}" \
    --missingDataAsZero

  echo "Matrix written to ${MATRIX_OUT}"

  plotHeatmap \
    --matrixFile "${MATRIX_OUT}" \
    --outFileName "${HEATMAP_OUT}" \
    --colorMap Reds \
    --heatmapHeight 12 \
    --heatmapWidth 4 \
    --xAxisLabel "distance from peak center (bp)" \
    --refPointLabel "center" \
    --plotTitle "Raw (unnormalised) signal: ${COMPARISON_NAME}"

  echo "Heatmap written to ${HEATMAP_OUT}"
  echo "=== Done with ${COMPARISON_NAME} ==="

done

echo "=== All comparisons complete ==="