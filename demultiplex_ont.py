#!/usr/bin/env python3
"""
demultiplex_ont.py — ONT long-read FASTQ demultiplexer

Searches for barcodes at the 5' and 3' ends of each read (forward and
reverse-complement orientations) using edit distance via edlib.

Reads are assigned to the barcode with the single lowest edit distance.
If two or more barcodes match with the same best distance, the read is
written to a separate *_ambiguous.fastq.gz file rather than being
misassigned.

Usage:
    python demultiplex_ont.py \
        -i "path/to/reads/*.fastq.gz" \
        -b barcodes.csv \
        -o output_dir \
        -d 200 \
        -e 4

Note: quote glob patterns so the shell passes them to the script unexpanded.
"""

import argparse
import os
import sys
import gzip
import glob
from collections import defaultdict
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
    """Return all (label, distance) pairs within max_errors in the segment."""
    valid_matches = []
    for label, barcode_seq in barcodes_dict.items():
        alignment = edlib.align(barcode_seq, sequence_segment, mode="HW", task="distance")
        distance = alignment['editDistance']
        if distance <= max_errors:
            valid_matches.append((label, distance))
    return valid_matches


def process_single_fastq(fastq_file, barcodes_to_check, barcode_labels, args):
    print(f"\n--- Processing file: {os.path.basename(fastq_file)} ---")

    file_basename = os.path.basename(fastq_file).replace(".fastq.gz", "").replace(".fq.gz", "")
    output_files = {}

    for label in barcode_labels:
        filename = f"{file_basename}_{label}_d{args.search_dist}_e{args.max_errors}.fastq.gz"
        output_files[label] = gzip.open(os.path.join(args.output_dir, filename), "wt")

    unmatched_filename = f"{file_basename}_unmatched_d{args.search_dist}_e{args.max_errors}.fastq.gz"
    output_files['unmatched'] = gzip.open(os.path.join(args.output_dir, unmatched_filename), "wt")

    ambiguous_filename = f"{file_basename}_ambiguous_d{args.search_dist}_e{args.max_errors}.fastq.gz"
    output_files['ambiguous'] = gzip.open(os.path.join(args.output_dir, ambiguous_filename), "wt")

    stats = {
        'total_reads': 0,
        'unambiguous_matches': 0,
        'unmatched': 0,
        'ambiguous': 0,
        'barcodes': {
            label: {'total': 0, 'head_fwd': 0, 'head_rev': 0, 'tail_fwd': 0, 'tail_rev': 0}
            for label in barcode_labels
        }
    }

    with gzip.open(fastq_file, "rt") as fq_in:
        for record in SeqIO.parse(fq_in, "fastq"):
            stats['total_reads'] += 1
            if stats['total_reads'] % 50000 == 0:
                print(f"  ...processed {stats['total_reads']} reads")

            seq = str(record.seq)
            head_segment = seq[:args.search_dist]
            tail_segment = seq[-args.search_dist:]

            head_matches = find_barcode_matches(head_segment, barcodes_to_check, args.max_errors)
            tail_matches = find_barcode_matches(tail_segment, barcodes_to_check, args.max_errors)

            all_possible_matches = head_matches + tail_matches

            if not all_possible_matches:
                stats['unmatched'] += 1
                SeqIO.write(record, output_files['unmatched'], "fastq")
                continue

            # Best distance per base barcode (ignoring orientation and end)
            per_bc = {}
            for label, d in all_possible_matches:
                base = label.rsplit('_', 1)[0]
                if base not in per_bc or d < per_bc[base]:
                    per_bc[base] = d
            sorted_bc = sorted(per_bc.items(), key=lambda x: x[1])

            # Ambiguous if the top two barcodes are within 1 edit distance of each other.
            if len(sorted_bc) > 1 and sorted_bc[1][1] - sorted_bc[0][1] <= 1:
                stats['ambiguous'] += 1
                SeqIO.write(record, output_files['ambiguous'], "fastq")
            else:
                stats['unambiguous_matches'] += 1
                base_label = sorted_bc[0][0]
                best_distance = sorted_bc[0][1]
                # Find the match for location/orientation stats
                winning = next(m for m in all_possible_matches
                               if m[0].rsplit('_', 1)[0] == base_label and m[1] == best_distance)
                _, orientation = winning[0].rsplit('_', 1)
                location = 'head' if winning in head_matches else 'tail'

                stats['barcodes'][base_label][f"{location}_{orientation}"] += 1
                stats['barcodes'][base_label]['total'] += 1
                SeqIO.write(record, output_files[base_label], "fastq")

    for f in output_files.values():
        f.close()

    total = stats['total_reads']
    print("\n--- Demultiplexing Summary ---")
    print(f"File: {os.path.basename(fastq_file)}")
    print(f"Total reads processed: {total}")
    print(f"Unambiguous matches: {stats['unambiguous_matches']} ({stats['unambiguous_matches']/total*100:.2f}%)")
    print(f"Ambiguous reads:     {stats['ambiguous']} ({stats['ambiguous']/total*100:.2f}%)")
    print(f"Unmatched reads:     {stats['unmatched']} ({stats['unmatched']/total*100:.2f}%)")
    print("-" * 35)

    for label in sorted(barcode_labels):
        c = stats['barcodes'][label]
        print(f"Barcode: {label} (Total: {c['total']})")
        print(f"  ├─ Found at Head: {c['head_fwd']} (fwd) + {c['head_rev']} (rev) = {c['head_fwd'] + c['head_rev']}")
        print(f"  └─ Found at Tail: {c['tail_fwd']} (fwd) + {c['tail_rev']} (rev) = {c['tail_fwd'] + c['tail_rev']}")

    print("------------------------------------")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Demultiplex long-read FASTQ files using barcode edit-distance matching.\n"
            "Searches the 5' and 3' ends of each read in both orientations.\n"
            "Ambiguous reads (equal best distance to multiple barcodes) are written\n"
            "to a separate file rather than misassigned.\n\n"
            "Quote glob patterns: -i \"data/*.fastq.gz\""
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input_pattern", required=True,
                        help="Input FASTQ file(s); glob patterns accepted (quote them).")
    parser.add_argument("-b", "--barcodes", required=True,
                        help="Barcode file (CSV: label,sequence, one per line).")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory for demultiplexed output files.")
    parser.add_argument("-d", "--search_dist", type=int, default=200,
                        help="Bases from each end to search for a barcode. Default: 200.")
    parser.add_argument("-e", "--max_errors", type=int, default=4,
                        help="Maximum edit distance (mismatches + indels) for a match. Default: 4.")
    args = parser.parse_args()

    fastq_files = glob.glob(args.input_pattern)
    if not fastq_files:
        print(f"Error: No files found matching '{args.input_pattern}'")
        sys.exit(1)

    if not os.path.exists(args.barcodes):
        print(f"Error: Barcode file not found at '{args.barcodes}'")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Found {len(fastq_files)} file(s) to process.")
    print(f"Output will be saved to: {args.output_dir}")

    barcodes_to_check = {}
    barcode_labels = set()
    try:
        with open(args.barcodes, 'r') as bc_file:
            for line in bc_file:
                if line.strip():
                    label, seq = line.strip().split(',')
                    barcode_labels.add(label)
                    barcodes_to_check[f"{label}_fwd"] = seq.upper()
                    barcodes_to_check[f"{label}_rev"] = reverse_complement(seq.upper())
        print(f"Successfully loaded {len(barcode_labels)} barcodes.")
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        sys.exit(1)

    for fq_file in fastq_files:
        process_single_fastq(fq_file, barcodes_to_check, barcode_labels, args)

    print("\nAll files processed. Demultiplexing complete.")


if __name__ == "__main__":
    main()
