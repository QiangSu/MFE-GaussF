import pandas as pd
import os
import glob
import gzip
from Bio import SeqIO
import argparse

# Normalization function
def normalize_frequency(kmer_count, total_csv_kmers, total_read_count, read_length, k):
    multiplier = 2 * 1000 * 1000000 / (total_csv_kmers * total_read_count * (read_length - k + 1))
    return kmer_count * multiplier

# Function to get the total count of reads in a FASTQ file
def get_total_read_count(fastq_filepath):
    with gzip.open(fastq_filepath, 'rt') as f:
        return sum(1 for _ in SeqIO.parse(f, 'fastq'))

# Function to merge k-mer count files to original k-mer CSV and add normalized count
def merge_kmer_counts_to_kmers(kmers_filepath, kmer_counts_filepath, total_csv_kmers, total_read_count, read_length, k, output_filepath):
    # Load the original kmer data and kmer counts
    kmers_df = pd.read_csv(kmers_filepath)
    kmer_counts_df = pd.read_csv(kmer_counts_filepath)

    # Merge the dataframes on the 'kmer' column
    merged_df = kmers_df.merge(kmer_counts_df, left_on='kmer', right_on='K-mer')

    # Drop the redundant 'K-mer' column from the merged dataframe
    merged_df.drop('K-mer', axis=1, inplace=True)

    # Add normalized k-mer counts
    merged_df['Normalized_K-mer_Count'] = merged_df['Count'].apply(normalize_frequency, args=(total_csv_kmers, total_read_count, read_length, k))

    # Add the 'Mini_Shared_Length' column, which contains the total number of kmers for each row
    merged_df['Mini_Shared_Length'] = total_csv_kmers

    # Write the merged data to a new CSV file
    merged_df.to_csv(output_filepath, index=False)

# Main function
def main(directory, output_directory, fastq_filepath, read_length, k):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Calculate the total number of reads in the FASTQ file
    total_read_count = get_total_read_count(fastq_filepath)
    print(f"Total number of reads in FASTQ file: {total_read_count}")

    # Process all the files
    for kmer_filepath in glob.glob(os.path.join(directory, "*_kmers.csv")):
        # Calculate the number of k-mers in the current CSV file
        total_csv_kmers = pd.read_csv(kmer_filepath).shape[0]
        print(f"Number of k-mers in {os.path.basename(kmer_filepath)}: {total_csv_kmers}")

        # Compute the filename of the corresponding kmer counts file
        basename = os.path.splitext(os.path.basename(kmer_filepath))[0]
        kmer_counts_filename = f"{basename}_counts.csv"
        kmer_counts_filepath = os.path.join(directory, kmer_counts_filename)

        # Define the output file path in the output directory
        output_filename = f"{basename}_merged_normalized.csv"
        output_filepath = os.path.join(output_directory, output_filename)

        # Merge the kmer counts into the kmer file and add normalized counts
        merge_kmer_counts_to_kmers(kmer_filepath, kmer_counts_filepath, total_csv_kmers, total_read_count, read_length, k, output_filepath)
        print(f"Merged and normalized file created at: {output_filepath}")

# Run script with command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge k-mer counts into original k-mer CSV files and add normalized counts.')
    parser.add_argument('--directory', type=str, required=True, help='Directory containing the k-mer and k-mer counts CSV files.')
    parser.add_argument('--output_directory', type=str, required=True, help='Output directory for storing merged CSV files.')
    parser.add_argument('--fastq', type=str, required=True, help='Path to the specific FASTQ.GZ file.')
    parser.add_argument('--read_length', type=int, default=150, help='Read length of FASTQ sequences.')
    parser.add_argument('--k', type=int, default=50, help='K-mer length.')
    args = parser.parse_args()
    main(args.directory, args.output_directory, args.fastq, args.read_length, args.k)
