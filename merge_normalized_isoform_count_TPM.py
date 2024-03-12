import os
import pandas as pd
import glob


def initial_normalize_frequency_RPK(kmer_count, csv_kmers, read_length, k):
    """
    Normalizes the k-mer count based on the number of k-mer occurrences in a single CSV file.
    """
    multiplier = 2 * 1000 / ((csv_kmers + 49) * (read_length - k + 1))
    return kmer_count * multiplier


def normalize_frequency(kmer_count, individual_csv_kmers, total_normalized_sum, read_length, k):
    """
    Normalizes the k-mer count based on the TPM formula, with csv_kmers being recalculated for each file.
    """
    multiplier = 2 * 1000 * 1000000 / ((individual_csv_kmers + 49) * total_normalized_sum * (read_length - k + 1))
    return kmer_count * multiplier


directory_path = './gene_folder/example_GaussF/kmer_distribution_frequency_result_1042353-N-1/testing_data/'
read_length = 150
k = 50

total_normalized_sum = 0

# First loop: Calculate total_normalized_sum
for filename in glob.glob(os.path.join(directory_path, '*counts.csv')):
    df = pd.read_csv(filename)
    csv_kmers = len(df)  # The total number of unique k-mers for this file

    df['Normalized'] = df.apply(
        lambda row: initial_normalize_frequency_RPK(row['Count'], csv_kmers, read_length, k), axis=1)

    # Here, Mini_Shared_Length is equal to the total csv_kmers for the file
    df['Mini_Shared_Length'] = csv_kmers

    file_normalized_sum = df['Normalized'].sum()
    total_normalized_sum += file_normalized_sum

# Print the total normalized sum after the first loop
print(f"Total normalized sum: {total_normalized_sum}")

# Second loop: Apply final TPM normalization and add Mini_Shared_Length specific to each CSV's k-mer count, then save files
for filename in glob.glob(os.path.join(directory_path, '*counts.csv')):
    df = pd.read_csv(filename)
    individual_csv_kmers = len(df)  # Recalculate for the specific file, as it might be different for each file

    df['TPM'] = df.apply(
        lambda row: normalize_frequency(row['Count'], individual_csv_kmers, total_normalized_sum, read_length, k),
        axis=1)

    # Adding Mini_Shared_Length based on individual_csv_kmers for consistency
    df['Mini_Shared_Length'] = individual_csv_kmers

    new_filename = filename.replace('.csv', '_with_TPM.csv')
    df.to_csv(new_filename, index=False)

print(f"TPM values and Mini_Shared_Length have been calculated and saved for all files.")

import os
import glob
import pandas as pd
import argparse


def initial_normalize_frequency_RPK(kmer_count, csv_kmers, read_length, k):
    multiplier = 2 * 1000 / ((csv_kmers + 49) * (read_length - k + 1))
    return kmer_count * multiplier


def normalize_frequency(kmer_count, individual_csv_kmers, total_normalized_sum, read_length, k):
    multiplier = (2 * 1000 * 1000000) / ((individual_csv_kmers + 49) * total_normalized_sum * (read_length - k + 1))
    return kmer_count * multiplier


def merge_kmer_counts_to_kmers(total_csv_kmers, total_read_count, read_length, k, kmers_filepath, kmer_counts_filepath,
                               output_filepath):
    kmers_df = pd.read_csv(kmers_filepath)
    kmer_counts_df = pd.read_csv(kmer_counts_filepath)

    merged_df = kmers_df.merge(kmer_counts_df, left_on='kmer', right_on='K-mer')
    merged_df.drop('K-mer', axis=1, inplace=True)

    merged_df['Normalized_K-mer_Count'] = merged_df['Count'].apply(
        lambda count: normalize_frequency(count, total_csv_kmers, total_read_count, read_length, k))
    merged_df['Mini_Shared_Length'] = total_csv_kmers

    merged_df.to_csv(output_filepath, index=False)


def main(directory, output_directory, read_length, k):
    total_normalized_sum = 0

    for filename in glob.glob(os.path.join(directory, '*counts.csv')):
        df = pd.read_csv(filename)
        df['Normalized'] = df.apply(lambda row: initial_normalize_frequency_RPK(row['Count'], len(df), read_length, k),
                                    axis=1)
        total_normalized_sum += df['Normalized'].sum()

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for kmer_filepath in glob.glob(os.path.join(directory, "*_kmers.csv")):
        basename = os.path.splitext(os.path.basename(kmer_filepath))[0]
        kmer_counts_filename = f"{basename}_counts_with_TPM.csv"
        kmer_counts_filepath = os.path.join(directory, kmer_counts_filename)
        total_csv_kmers = pd.read_csv(kmer_counts_filepath).shape[0]

        output_filename = f"{basename}_merged_normalized.csv"
        output_filepath = os.path.join(output_directory, output_filename)

        if os.path.exists(kmer_counts_filepath):
            merge_kmer_counts_to_kmers(total_csv_kmers, total_normalized_sum, read_length, k, kmer_filepath,
                                       kmer_counts_filepath, output_filepath)
            print(f"Merged and normalized file created at: {output_filepath}")
        else:
            print(f"No matching k-mer count file found for {basename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Merge k-mer counts into original k-mer CSV files and add normalized counts.')
    parser.add_argument('--directory', type=str, required=True,
                        help='Directory containing the k-mer and k-mer counts CSV files.')
    parser.add_argument('--output_directory', type=str, required=True,
                        help='Output directory for storing merged CSV files.')

    parser.add_argument('--read_length', type=int, default=150, help='Read length of sequences.')
    parser.add_argument('--k', type=int, default=50, help='K-mer length.')

    args = parser.parse_args()
    main(args.directory, args.output_directory, args.read_length, args.k)
