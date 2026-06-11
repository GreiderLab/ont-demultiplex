#!/usr/bin/env python3
"""
evaluate_error_rates.py — Error-rate sensitivity vs. specificity for ONT demultiplexing

Tests a range of edit-distance thresholds on a FASTQ file and reports, for each
threshold, how many reads would be unambiguously assigned vs. ambiguous.
Use this to choose the --max_errors value before running demultiplex_ont.py.

Reads are cached in memory so each error rate requires only one alignment pass
(no repeated I/O). For very large files (>10 M reads) RAM may be the bottleneck.

Usage:
    python evaluate_error_rates.py \
        -i reads.fastq.gz \
        -b barcodes.csv \
        -d 200 \
        -r 0-10
"""

import argparse
import os
import sys
import gzip
import glob
from Bio import SeqIO
from Bio.Seq import Seq

try:
    import edlib
except ImportError:
    print("Error: 'edlib' is required. Install with: pip install edlib")
    sys.exit(1)


def reverse_complement(seq_string):
    return str(Seq(seq_string).reverse_complement())


def find_barcode_matches(sequence_segment, barcodes_dict, max_errors):
    valid_matches = []
    for label, barcode_seq in barcodes_dict.items():
        alignment = edlib.align(barcode_seq, sequence_segment, mode="HW", task="distance")
        distance = alignment['editDistance']
        if distance <= max_errors:
            valid_matches.append((label, distance))
    return valid_matches


def analyze_fastq_error_rates(fastq_file, barcodes_to_check, error_range, args):
    print(f"\n--- Analyzing Error Rates for: {os.path.basename(fastq_file)} ---")
    print("Caching reads in memory for faster processing...")

    try:
        with gzip.open(fastq_file, "rt") as fq_in:
            records = list(SeqIO.parse(fq_in, "fastq"))
    except Exception as e:
        print(f"Error reading FASTQ file: {e}")
        return

    total_reads = len(records)
    print(f"Found {total_reads} total reads.")

    results = {}
    for e in error_range:
        print(f"  Testing error rate: {e}...")
        unambiguous_count = 0
        ambiguous_count = 0

        for record in records:
            seq = str(record.seq)
            head_matches = find_barcode_matches(seq[:args.search_dist], barcodes_to_check, e)
            tail_matches = find_barcode_matches(seq[-args.search_dist:], barcodes_to_check, e)

            all_matches = head_matches + tail_matches
            if not all_matches:
                continue

            best_distance = min(m[1] for m in all_matches)
            num_best = sum(1 for m in all_matches if m[1] == best_distance)

            if num_best > 1:
                ambiguous_count += 1
            else:
                unambiguous_count += 1

        results[e] = {'unambiguous': unambiguous_count, 'ambiguous': ambiguous_count}

    print("\n--- Error Rate Sensitivity vs. Specificity ---")
    print(f"File: {os.path.basename(fastq_file)}")
    print(f"Total Reads: {total_reads}")
    print("=" * 65)
    print(f"{'Error Rate':<12} | {'Unambiguous Matches':<22} | {'Ambiguous Reads':<20}")
    print(f"{'------------':<12}-+-{'----------------------':<22}-+-{'--------------------':<20}")

    for e, counts in results.items():
        ua_pct = counts['unambiguous'] / total_reads * 100 if total_reads > 0 else 0
        am_pct = counts['ambiguous'] / total_reads * 100 if total_reads > 0 else 0
        print(f"{e:<12} | {counts['unambiguous']} ({ua_pct:.2f}%){'':<{10 - len(str(counts['unambiguous']))}} | {counts['ambiguous']} ({am_pct:.2f}%)")

    print("=" * 65)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate demultiplexing sensitivity vs. specificity across a range of\n"
            "edit-distance thresholds. Helps you choose --max_errors for demultiplex_ont.py.\n\n"
            "A cliff in ambiguous reads (e.g., from 0.02% to 0.64%) marks the threshold\n"
            "above which the error rate is too permissive."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input_pattern", required=True,
                        help="Input FASTQ file(s); glob patterns accepted.")
    parser.add_argument("-b", "--barcodes", required=True,
                        help="Barcode file (CSV: label,sequence, one per line).")
    parser.add_argument("-d", "--search_dist", type=int, default=200,
                        help="Bases from each end to search. Default: 200.")
    parser.add_argument("-r", "--error_range", type=str, default="0-10",
                        help="Range of error rates to test, e.g. '0-10'. Default: 0-10.")
    args = parser.parse_args()

    try:
        start, end = map(int, args.error_range.split('-'))
        error_range = range(start, end + 1)
    except ValueError:
        print("Error: --error_range must be in 'start-end' format, e.g. '0-10'.")
        sys.exit(1)

    fastq_files = glob.glob(args.input_pattern)
    if not fastq_files:
        print(f"Error: No files found matching '{args.input_pattern}'")
        sys.exit(1)

    if not os.path.exists(args.barcodes):
        print(f"Error: Barcode file not found at '{args.barcodes}'")
        sys.exit(1)

    barcodes_to_check = {}
    try:
        with open(args.barcodes, 'r') as bc_file:
            for line in bc_file:
                if line.strip():
                    label, seq = line.strip().split(',')
                    barcodes_to_check[f"{label}_fwd"] = seq.upper()
                    barcodes_to_check[f"{label}_rev"] = reverse_complement(seq.upper())
        print(f"Successfully loaded {len(barcodes_to_check)//2} barcodes.")
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        sys.exit(1)

    for fq_file in fastq_files:
        analyze_fastq_error_rates(fq_file, barcodes_to_check, error_range, args)

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
