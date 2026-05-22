#!/bin/bash
#SBATCH --job-name=greenlist-norm
#SBATCH -N 1
#SBATCH --partition=your_partition        # adjust for your cluster
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=20
#SBATCH -t 10:00:00
#SBATCH --mail-user=your@email.com        # adjust or remove
#SBATCH --mail-type=ALL
#SBATCH --output=logs/greenlist-norm-%j.out
#SBATCH --error=logs/greenlist-norm-%j.err
#SBATCH --container=/path/to/container.sif  # adjust for your cluster

# ── Environment ───────────────────────────────────────────────────────────────
source /etc/profile.d/modules.sh
module load deeptools/3.5.0
module load R/4.2.0

# ── Directories ───────────────────────────────────────────────────────────────
bam_dir="path/to/bam_files/hg38"
bed_dir="path/to/CUT-RUN_greenlist-main"
out_dir="path/to/output/greenlist_norm"
script_dir="path/to/CUT-RUN_greenlist-main/original_scripts"

# ── BAM files ─────────────────────────────────────────────────────────────────
# CUT&RUN CTCF time-course BAMs (Skene & Henikoff 2017, hg38)
# Naming convention: CTCF_{condition}_{date}_{...}_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
bams=(
    CTCF_0_6million_20160315_PS_HsDm_CTCF_0_6_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_2_5million_20160315_PS_HsDm_CTCF_2_5_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_10million_20160315_PS_HsDm_CTCF_10_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_5s_20160628_PS_HsDm_CTCF_0602_A_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_15s_20160628_PS_HsDm_CTCF_0602_B_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_45s_20160628_PS_HsDm_CTCF_0602_C_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_3min_20160628_PS_HsDm_CTCF_0602_D_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_1min_20160414_PS_HsDm_CTCF_1m_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_2min_20160414_PS_HsDm_CTCF_2m_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_4min_20160414_PS_HsDm_CTCF_4m_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_8min_20160414_PS_HsDm_CTCF_8m_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_9min_20160628_PS_HsDm_CTCF_0602_E_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_7_5min_20160401_PS_HsDm_CTCF_7_5_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_15min_20160401_PS_HsDm_CTCF_15_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_27min_20160628_PS_HsDm_CTCF_0602_F_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_30min_20160401_PS_HsDm_CTCF_30_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_45min_20160401_PS_HsDm_CTCF_45_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
    CTCF_disrupted_20160315_PS_HsDm_CTCF_expl_UB_Homo_sapiens_OTHER_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam
)

# ── Step 1: Quantify reads over greenlist regions ─────────────────────────────
cd "$bam_dir"

multiBamSummary BED-file \
    --BED "$bed_dir/hg38_CUTnRUN_greenlist.v1.bed" \
    --smartLabels \
    -e \
    --centerReads \
    -o glist_quant.npz \
    -b "${bams[@]}" \
    --outRawCounts "$out_dir/output"

# ── Step 2: Tidy output table ─────────────────────────────────────────────────
# Remove comment characters, merge chr:start columns, drop start column
cat "$out_dir/output" \
    | tr -d "'#" \
    | sed $'s/\t/_/1' \
    | cut -f 1,3- \
    > "$out_dir/glist_quant.tsv"

echo "Counts calculated successfully."

# ── Step 3: Calculate greenlist size factors ──────────────────────────────────
Rscript --vanilla \
    "$script_dir/get_sizeFactors.R" \
    "$out_dir/glist_quant.tsv" \
    "$out_dir/glist_sizeFactors.tsv"
