#!/bin/bash
# =============================================================================
# Greenlist Normalization: Read Quantification and Size Factor Estimation
# Usage: sbatch greenlist_norm.sh
# =============================================================================

#SBATCH --job-name=greenlist_norm
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=20
#SBATCH --time=10:00:00
#SBATCH --mail-user=krutika.ambani@unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH --output=greenlist_norm-%j.out
#SBATCH --error=greenlist_norm-%j.err

# --- Modules -----------------------------------------------------------------

module load deepTools/3.5.2
module load R/4.5.0

# --- Paths -------------------------------------------------------------------

project_dir="/path/to/project"
bam_dir="$project_dir/results/bam_files/final_bams/hg38"
bed_dir="/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/greenlist_norm_scripts/CUT-RUN_greenlist-main"
out_dir="$project_dir/R_analysis/5.Diffbind_comparisons_greenlist_normalisation_incl/greenlist_normfacs"

# --- Quantify reads over greenlist regions -----------------------------------

cd $bam_dir

multiBamSummary BED-file \
  --BED $bed_dir/hg38_CUTnRUN_greenlist.v1.bed \
  --smartLabels -e --centerReads \
  -o glist_quant.npz \
  -b 001-MCF7-Par-DMSO-R1-BRG1-native_S11_hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam \
     002-MCF7-Par-DMSO-R2-BRG1-native_S22_hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam \
     003-MCF7-Par-Abema-R1-BRG1-native_S33_hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam \
     004-MCF7-Par-Abema-R2-BRG1-native_S44_hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam \
  --outRawCounts $out_dir/output

# --- Reformat counts table ---------------------------------------------------

cat $out_dir/output | tr -d "'#" | sed $'s/\t/_/1' | cut -f 1,3- > $out_dir/glist_quant.tsv

echo "Counts were calculated!"

# --- Compute size factors ----------------------------------------------------

Rscript --vanilla $bed_dir/original_scripts/get_sizeFactors.R \
  $out_dir/glist_quant.tsv \
  $out_dir/glist_sizeFactors.tsv
