#!/bin/bash
# =============================================================================
# Greenlist Normalisation: Read Quantification and Size Factor Estimation
# Usage: sbatch greenlist_norm.sh
# =============================================================================

#SBATCH --job-name=greenlist_norm
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=greenlist_norm-%j.out
#SBATCH --error=greenlist_norm-%j.err

# --- Setup -------------------------------------------------------------------

source /etc/profile.d/modules.sh

# --- Modules -----------------------------------------------------------------

module load deeptools/3.5.0
module load R/4.2.0

# --- Paths -------------------------------------------------------------------

bam_dir="/path/to/project/results/bam_files/final_bams/hg38"
bed_dir="/path/to/CUT-RUN_greenlist-main"
out_dir="/path/to/project/R_analysis/1.1Greenlist_norm/greenlist_norm"

# --- Step 1: Quantify reads over greenlist regions ---------------------------

cd $bam_dir

bam_suffix="onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam"

multiBamSummary BED-file \
  --BED $bed_dir/hg38_CUTnRUN_greenlist.v1.bed \
  --smartLabels \
  --extendReads \
  --centerReads \
  --outFileName glist_quant.npz \
  --outRawCounts $out_dir/output \
  --bamfiles \
    MCF7-Complete-R1-ER_S14_${bam_suffix} \
    MCF7-Complete-R2-ER_S15_${bam_suffix} \
    MCF7-CSS-R1-ER_S18_${bam_suffix} \
    MCF7-CSS-R2-ER_S19_${bam_suffix} \
    MCF7-CSS-E2-R1-ER_S20_${bam_suffix} \
    MCF7-CSS-E2-R2-ER_S21_${bam_suffix} \
    sgER-Complete-R1-ER_S16_${bam_suffix} \
    sgER-Complete-R2-ER_S17_${bam_suffix}

# --- Step 2: Reformat raw counts table ---------------------------------------

cat $out_dir/output | tr -d "'#" | sed $'s/\t/_/1' | cut -f 1,3- \
  > $out_dir/glist_quant.tsv

echo "Read quantification complete."

# --- Step 3: Compute size factors --------------------------------------------

Rscript --vanilla $bed_dir/original_scripts/get_sizeFactors.R \
  $out_dir/glist_quant.tsv \
  $out_dir/glist_sizeFactors.tsv

echo "Size factors written to: $out_dir/glist_sizeFactors.tsv"
