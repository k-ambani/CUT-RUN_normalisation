#!/bin/bash
# =============================================================================
# BAM and Peak Count Summary
# Usage: sbatch bam_bw_peaks_summary.sh
# =============================================================================

#SBATCH --job-name=bam_bw_peaks_summary
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=user@institute.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=rhel_short
#SBATCH --output=bam_bw_peaks_summary-%j.out
#SBATCH --error=bam_bw_peaks_summary-%j.err

# --- Setup -------------------------------------------------------------------

module load samtools/1.19.2-gcc-13.2.0
set -xe

# --- Paths -------------------------------------------------------------------

project_dir="/path/to/project"
scratch_dir="/path/to/scratch/bam_files/intermediate_bams"

summary_out_dir="$project_dir/results"
bam_hg38_final="$project_dir/results/bam_files/final_bams/hg38"
bam_dm6_final="$project_dir/results/bam_files/final_bams/dm6"
bam_sacCer3_final="$project_dir/results/bam_files/final_bams/sacCer3"
bam_hg38_intermediate="$scratch_dir/hg38"
bam_dm6_intermediate="$scratch_dir/dm6"
bam_sacCer3_intermediate="$scratch_dir/sacCer3"
peaks_withinput_dir="$project_dir/results/macs2/withinput/macs2_q0.05/narrowPeak"

# --- Summarize BAM counts, scale factors, and peak counts --------------------

cd $bam_hg38_final

echo "basename hg38_dm6_sacCer3_initial_bam_counts hg38_final_bam_counts hg38_intermediate_bam_counts dm6_final_bam_counts dm6_intermediate_bam_counts sacCer3_final_bam_counts sacCer3_intermediate_bam_counts dm6_scale_factor sacCer3_scale_factor peak_counts_withinput" \
  >> $summary_out_dir/summary_bam_counts.txt

for bam in *.bam; do

  base=$(basename ${bam} _onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)

  # BAM read counts
  initial_bam_counts_hg38_dm6_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${scratch_dir}/../${base}_allgenomes_sorted.bam)
  final_bam_counts_hg38=$(samtools view -c -F 260 \
    ${bam_hg38_final}/${base}_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  intermediate_bam_counts_hg38=$(samtools view -c -q 2 -F 260 \
    ${bam_hg38_intermediate}/${base}_onlyhg38.bam)
  final_bam_counts_dm6=$(samtools view -c -F 260 \
    ${bam_dm6_final}/${base}_onlydm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  intermediate_bam_counts_dm6=$(samtools view -c -q 2 -F 260 \
    ${bam_dm6_intermediate}/${base}_onlydm6.bam)
  final_bam_counts_sacCer3=$(samtools view -c -F 260 \
    ${bam_sacCer3_final}/${base}_onlysacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  intermediate_bam_counts_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${bam_sacCer3_intermediate}/${base}_onlysacCer3.bam)

  # Spike-in scale factors (1,000,000 / spike-in read count)
  total_number_of_reads_dm6=$(samtools view -c -q 2 -F 260 \
    ${bam_dm6_intermediate}/${base}_onlydm6.bam)
  dm6_scalefactor=$(echo "scale=3; 1000000 / ${total_number_of_reads_dm6}" | bc)

  total_number_of_reads_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${bam_sacCer3_intermediate}/${base}_onlysacCer3.bam)
  sacCer3_scalefactor=$(echo "scale=3; 1000000 / ${total_number_of_reads_sacCer3}" | bc)

  # Peak counts (with input)
  peak_withinput_file=${peaks_withinput_dir}/${base}_onlyhg38_peaks.narrowPeak
  if [ ! -f ${peak_withinput_file} ]; then
    peaks_counts_withinput="no_peaks"
  else
    peaks_counts_withinput=$(wc -l ${peak_withinput_file} | cut -d ' ' -f 1)
  fi

  echo "${base} ${initial_bam_counts_hg38_dm6_sacCer3} ${final_bam_counts_hg38} ${intermediate_bam_counts_hg38} ${final_bam_counts_dm6} ${intermediate_bam_counts_dm6} ${final_bam_counts_sacCer3} ${intermediate_bam_counts_sacCer3} ${dm6_scalefactor} ${sacCer3_scalefactor} ${peaks_counts_withinput}" \
    >> $summary_out_dir/summary_bam_counts.txt

done
