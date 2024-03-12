from collections import defaultdict
import csv
import os
import argparse

def read_fasta(file_path):
    headers = []
    sequences = []
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ""
                headers.append(line[1:])  # Remove the '>' character
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return headers, sequences

def sanitize_filename(header):
    return header.replace('|', '_')

# Command-line argument parsing
parser = argparse.ArgumentParser(description="Process a FASTA file and output kmers csv files.")
parser.add_argument('--input', required=True, help="Path to the input FASTA file.")
parser.add_argument('--output', required=True, help="Path to the output directory where CSV files will be saved.")
args = parser.parse_args()

# Read FASTA file and prepare output directory using the arguments
transcript_headers, transcripts = read_fasta(args.input)
output_directory = args.output
os.makedirs(output_directory, exist_ok=True)

kmer_length = 50
global_kmer_counts = defaultdict(int)
kmer_transcript_sets = defaultdict(set)

# Kmer counting and tracking transcripts that contain them
for isoform_index, sequence in enumerate(transcripts):
    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i + kmer_length]
        global_kmer_counts[kmer] += 1
        kmer_transcript_sets[kmer].add(isoform_index)

# Creating CSV files for each isoform with their kmers
for isoform_index, header in enumerate(transcript_headers):
    output_csv_path = os.path.join(output_directory, sanitize_filename(header) + '_kmers.csv')

    # First pass to determine the minimum global frequency for this specific isoform
    min_global_frequency = float('inf')
    for i in range(len(transcripts[isoform_index]) - kmer_length + 1):
        kmer = transcripts[isoform_index][i:i + kmer_length]
        global_freq = global_kmer_counts[kmer]
        min_global_frequency = min(min_global_frequency, global_freq)

    # Prepare for the second pass: collect all potential rows to write
    rows_to_write = []
    local_kmer_counts = defaultdict(int)

    for i in range(len(transcripts[isoform_index]) - kmer_length + 1):
        kmer = transcripts[isoform_index][i:i + kmer_length]
        local_kmer_counts[kmer] += 1
        global_freq = global_kmer_counts[kmer]
        # Only consider kmers with the minimum global frequency for this isoform
        if global_freq == min_global_frequency:
            # Format 'Present_in_Transcripts' field based on the global frequency
            if global_freq > 1:
                transcripts_containing_kmer = '-'.join(sanitize_filename(transcript_headers[transcript_idx])
                                                       for transcript_idx in sorted(kmer_transcript_sets[kmer]))
            else:
                transcripts_containing_kmer = sanitize_filename(transcript_headers[list(kmer_transcript_sets[kmer])[0]])

            rows_to_write.append((kmer, local_kmer_counts[kmer], global_freq, transcripts_containing_kmer))

    # Identify the highest local frequency content within "Present_in_Transcripts"
    max_transcript_content_frequency = 0
    transcript_content_freq = defaultdict(int)
    for _, _, _, transcripts_content in rows_to_write:
        transcript_content_freq[transcripts_content] += 1
        max_transcript_content_frequency = max(max_transcript_content_frequency,
                                               transcript_content_freq[transcripts_content])

    # Write to CSV only for rows with highest "Present_in_Transcripts" local frequency
    with open(output_csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])

        for row in rows_to_write:
            if transcript_content_freq[row[3]] == max_transcript_content_frequency:
                csv_writer.writerow(row)

print(f"Kmers for each transcript have been saved as CSV files in the directory: {output_directory}")
