#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict
import os

# --- Input check ---
if len(sys.argv) < 2 or not sys.argv[1].endswith(('.fastq.gz', '.fq.gz')):
    sys.exit("Error: Please provide at least one input FASTQ file (R1) ending with .fastq.gz or .fq.gz\nUsage: python3 script.py read1.fastq.gz")

# --- Configurable paths ---
INDEX1_BARCODE_FILE = "index1.coding.txt"  # Tab-delimited file: label \t sequence
READ1_FASTQ = sys.argv[1]  # Input R1 fastq.gz file

filename_prefix = os.path.basename(READ1_FASTQ).split(".")[0]

debug = False

# Constant sequence anchors
primer_indexes = ['AGATCG']
primer_constant_5prime_alpha = "AGCAG"
primer_constant_5prime_alpha_distance = 31
primer_constant_5prime_beta = "GCTAGTGCTAG"
primer_constant_5prime_beta_distance = 11
barcode_3prime_constant = "AAGTA"

# Load Index1 barcodes from file
def load_index1_barcodes(filename):
    barcode_dict = dict()
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            label, seq = parts
            barcode_dict[label] = seq
    return barcode_dict

# Find all occurrences of a subsequence
def find_all_positions(line=None, index=None):
    positions = []
    if line is None or index is None:
        return positions
    search_start = 0
    while True:
        position = line.find(index, search_start)
        if position == -1:
            break
        positions.append(position)
        search_start = position + 1
    return positions

# Extract barcode, middle fragment, and i7 index from a read
def extract_index(line, barcode_dict):
    index_set = dict()
    found_primer_indexes = dict()
    found_barcode_indexes = dict()

    found_alpha = find_all_positions(line, primer_constant_5prime_alpha)
    found_beta = find_all_positions(line, primer_constant_5prime_beta)

    for primer_index in primer_indexes:
        positions = find_all_positions(line, primer_index)
        if positions:
            found_primer_indexes[primer_index] = positions

    for barcode_label, barcode_seq in barcode_dict.items():
        pattern = barcode_seq + barcode_3prime_constant
        positions = find_all_positions(line, pattern)
        if positions:
            found_barcode_indexes[barcode_label] = (barcode_seq, positions)

    for barcode_label in found_barcode_indexes:
        barcode_seq, positions = found_barcode_indexes[barcode_label]
        for barcode_pos in positions:
            best_primer_pos = None
            best_padding = len(line)

            for primer, primer_positions in found_primer_indexes.items():
                for primer_pos in primer_positions:
                    alpha_found = any(abs(primer_pos - alpha - primer_constant_5prime_alpha_distance) < 3 for alpha in found_alpha)
                    beta_found = any(abs(primer_pos - beta - primer_constant_5prime_beta_distance) < 3 for beta in found_beta)
                    if not (alpha_found and beta_found):
                        continue
                    padding = 0  # Can be refined if needed
                    if padding < best_padding:
                        best_padding = padding
                        best_primer_pos = primer_pos

            if best_primer_pos is not None:
                middle_start = barcode_pos + len(barcode_seq)
                middle_seq = line[middle_start:middle_start + 33] if middle_start + 33 <= len(line) else ""
                i7_start = best_primer_pos + 34
                i7_seq = line[i7_start:i7_start + 10] if i7_start + 10 <= len(line) else ""
                return ",".join([barcode_seq, middle_seq, i7_seq])

    return None

# --- Main ---
barcode_dict = load_index1_barcodes(INDEX1_BARCODE_FILE)
if debug:
    print("Loaded barcodes:", barcode_dict)

index_counts = defaultdict(int)
missing_reads = []
total_reads = 0

with gzip.open(READ1_FASTQ, 'rt') as f1:
    count = 1
    for line1 in f1:
        if count == 2:
            total_reads += 1
            read1_seq = line1.strip()
            result = extract_index(read1_seq, barcode_dict)
            if result:
                index_counts[result] += 1
            else:
                missing_reads.append(read1_seq)
        count += 1
        if count == 5:
            count = 1

# Write output files
with open(f"{filename_prefix}.output.txt", "w") as out_f:
    for key in sorted(index_counts, key=index_counts.get, reverse=True):
        out_f.write(f"{key}\t{index_counts[key]}\n")

with open(f"{filename_prefix}.missing.txt", "w") as missing_f:
    for line in missing_reads:
        missing_f.write(f"{line}\n")

# Report summary
missing_pct = (len(missing_reads) / total_reads) * 100 if total_reads > 0 else 0
print(f"{filename_prefix}: total_reads={total_reads}, matched={total_reads - len(missing_reads)}, missing={len(missing_reads)} ({missing_pct:.2f}%)")

