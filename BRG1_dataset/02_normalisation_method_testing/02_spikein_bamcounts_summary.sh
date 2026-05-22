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
#SBATCH --output=bam_bw_peaks_summary-%j.out
#SBATCH --error=bam_bw_peaks_summary-%j.err

# --- Setup -------------------------------------------------------------------

module load samtools/1.9

# --- Paths -------------------------------------------------------------------

working_dir="/path/to/project"
scratch_dir="/path/to/scratch"

summary_out_dir="$working_dir/results"
bam_hg38_final="$working_dir/results/bam_files/final_bams/hg38"
bam_sacCer3_final="$working_dir/results/bam_files/final_bams/sacCer3"
bam_hg38_intermediate="$scratch_dir/hg38"
bam_sacCer3_intermediate="$scratch_dir/sacCer3"
bam_hg38_sacCer3_initial="$scratch_dir"
peaks_withinput_dir="$working_dir/results/macs2/with_input/macs2_q0.05/narrowPeak"
peaks_withoutinput_dir="$working_dir/results/macs2/with_no_input/macs2_q0.05/narrowPeak"

# --- Summarize BAM counts, scale factors, and peak counts --------------------

cd $bam_hg38_final

echo "basename hg38_sacCer3_initial_bam_counts hg38_final_bam_counts hg38_intermediate_bam_counts sacCer3_final_bam_counts sacCer3_intermediate_bam_counts scale_factor peak_counts_withinput peak_counts_withoutinput" \
  >> $summary_out_dir/summary_hg38_sacCer3_bam_counts.txt

for bam in *.bam; do

  base=$(basename ${bam} _hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)

  # BAM read counts
  final_bam_counts_hg38=$(samtools view -c -F 260 \
    ${bam_hg38_final}/${base}_hg38_nosacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  intermediate_bam_counts_hg38=$(samtools view -c -q 2 -F 260 \
    ${bam_hg38_intermediate}/${base}_hg38_nosacCer3.bam)
  final_bam_counts_sacCer3=$(samtools view -c -F 260 \
    ${bam_sacCer3_final}/${base}_nohg38_sacCer3_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam)
  intermediate_bam_counts_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${bam_sacCer3_intermediate}/${base}_nohg38_sacCer3.bam)
  initial_bam_counts_hg38_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${bam_hg38_sacCer3_initial}/${base}_hg38_sacCer3_sorted.bam)

  # Spike-in scale factor (1,000,000 / sacCer3 read count)
  total_number_of_reads_sacCer3=$(samtools view -c -q 2 -F 260 \
    ${bam_sacCer3_intermediate}/${base}_nohg38_sacCer3.bam)
  scalefactor=$(echo "scale=3; 1000000 / ${total_number_of_reads_sacCer3}" | bc)

  # Peak counts (with input)
  peak_withinput_file=${peaks_withinput_dir}/${base}_hg38_nosacCer3_peaks.narrowPeak
  if [ ! -f ${peak_withinput_file} ]; then
    peaks_counts_withinput="no_peaks"
  else
    peaks_counts_withinput=$(wc -l ${peak_withinput_file} | cut -d ' ' -f 1)
  fi

  # Peak counts (without input)
  peak_withoutinput_file=${peaks_withoutinput_dir}/${base}_hg38_nosacCer3_peaks.narrowPeak
  peak_counts_withoutinput=$(wc -l ${peak_withoutinput_file} | cut -d ' ' -f 1)

  echo "${base} ${initial_bam_counts_hg38_sacCer3} ${final_bam_counts_hg38} ${intermediate_bam_counts_hg38} ${final_bam_counts_sacCer3} ${intermediate_bam_counts_sacCer3} ${scalefactor} ${peaks_counts_withinput} ${peak_counts_withoutinput}" \
    >> $summary_out_dir/summary_hg38_sacCer3_bam_counts.txt

done
