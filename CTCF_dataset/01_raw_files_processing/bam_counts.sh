#!/bin/bash
#SBATCH --job-name=bam-summary
#SBATCH -N 1
#SBATCH --cpus-per-task=15
#SBATCH -t 4:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --partition=your_partition        # adjust for your cluster
#SBATCH --output=logs/bam-summary-%j.out
#SBATCH --error=logs/bam-summary-%j.err

# ── Environment ───────────────────────────────────────────────────────────────
module load samtools/1.9
set -euo pipefail

# ── Directories ───────────────────────────────────────────────────────────────
working_dir="path/to/2017_Skene_Henikoff"
scratch_dir="path/to/scratch/2017_Skene_Henikoff"

summary_out_dir="$working_dir/results"

# Final BAMs (processed, filtered)
bam_hg38_final="$working_dir/results/bam_files/final_bams/hg38"
bam_dm6_final="$working_dir/results/bam_files/final_bams/dm6"

# Intermediate BAMs (on scratch)
bam_hg38_dm6_initial="$scratch_dir/results/bam_files/intermediate_bams"
bam_hg38_intermediate="$scratch_dir/results/bam_files/intermediate_bams/hg38"
bam_dm6_intermediate="$scratch_dir/results/bam_files/intermediate_bams/dm6"

# Peak calls (MACS2, q < 0.01)
peaks_withinput_dir="$working_dir/results/macs2/macs2_q0.01/with_neg_control/narrowPeak"
peaks_withoutinput_dir="$working_dir/results/macs2/macs2_q0.01/with_no_neg_control/narrowPeak"

# Output summary file
summary_file="$summary_out_dir/summary_hg38_dm6_bam_counts.txt"

# ── Suffix definitions ────────────────────────────────────────────────────────
# Defined once here so they're easy to update if filenames change
hg38_final_suffix="_hg38_nodm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam"
dm6_final_suffix="_nohg38_dm6_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam"
hg38_inter_suffix="_hg38_nodm6.bam"
dm6_inter_suffix="_nohg38_dm6.bam"
initial_suffix="_hg38_dm6_sorted.bam"

# ── Header ────────────────────────────────────────────────────────────────────
echo "basename hg38_dm6_initial_bam_counts hg38_final_bam_counts hg38_intermediate_bam_counts \
dm6_final_bam_counts dm6_intermediate_bam_counts scale_factor \
peak_counts_withinput peak_counts_withoutinput" \
    > "$summary_file"   # note: > not >> so the file is fresh each run

# ── Per-sample loop ───────────────────────────────────────────────────────────
cd "$bam_hg38_final"

for bam in *.bam; do

    base=$(basename "$bam" "$hg38_final_suffix")

    # BAM read counts (primary mapped reads: -F 260 excludes unmapped + secondary)
    final_bam_counts_hg38=$(samtools view -c -F 260 \
        "$bam_hg38_final/${base}${hg38_final_suffix}")

    intermediate_bam_counts_hg38=$(samtools view -c -q 2 -F 260 \
        "$bam_hg38_intermediate/${base}${hg38_inter_suffix}")

    final_bam_counts_dm6=$(samtools view -c -F 260 \
        "$bam_dm6_final/${base}${dm6_final_suffix}")

    intermediate_bam_counts_dm6=$(samtools view -c -q 2 -F 260 \
        "$bam_dm6_intermediate/${base}${dm6_inter_suffix}")

    initial_bam_counts_hg38_dm6=$(samtools view -c -q 2 -F 260 \
        "$bam_hg38_dm6_initial/${base}${initial_suffix}")

    # Scale factor: 1,000,000 / dm6 mapped reads (spike-in normalisation)
    dm6_reads=$(samtools view -c -q 2 -F 260 \
        "$bam_dm6_intermediate/${base}${dm6_inter_suffix}")
    scale_factor=$(echo "scale=3; 1000000 / ${dm6_reads}" | bc)

    # Peak counts — with input (negative control)
    peak_withinput_file="${peaks_withinput_dir}/${base}_hg38_nodm6_peaks.narrowPeak"
    if [[ -f "$peak_withinput_file" ]]; then
        peaks_counts_withinput=$(wc -l < "$peak_withinput_file")
    else
        peaks_counts_withinput="no_peaks_file"
    fi

    # Peak counts — without input
    peak_withoutinput_file="${peaks_withoutinput_dir}/${base}_hg38_nodm6_peaks.narrowPeak"
    if [[ -f "$peak_withoutinput_file" ]]; then
        peak_counts_withoutinput=$(wc -l < "$peak_withoutinput_file")
    else
        peak_counts_withoutinput="no_peaks_file"
    fi

    echo "$base $initial_bam_counts_hg38_dm6 $final_bam_counts_hg38 \
$intermediate_bam_counts_hg38 $final_bam_counts_dm6 $intermediate_bam_counts_dm6 \
$scale_factor $peaks_counts_withinput $peak_counts_withoutinput" \
        >> "$summary_file"

done

echo "Summary written to: $summary_file"
