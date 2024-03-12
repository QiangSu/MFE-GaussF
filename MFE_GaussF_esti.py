import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import argparse
import re
import subprocess


# Define a function to judge the suitability of the mean before fitting
def suitable_criteria_for_mfe(xc_fitted_local, w_fitted_local):
    return abs(xc_fitted_local) > 0.2 * w_fitted_local


# Parse command line arguments
parser = argparse.ArgumentParser(description="Analyze k-mer MFE and fit Gaussian CDF.")
parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the CSV files.")
parser.add_argument('--output', type=str, required=True, help="Path and name of the output file to save the results.")
parser.add_argument('--threshold', type=int, default=10,
                    help="Minimum number of k-mers required for fitting. Default is 10.")
args = parser.parse_args()


# Function to calculate MFE using RNAfold
def calculate_mfe(sequence):
    sequence = sequence.replace('T', 'U')  # Convert T to U for RNA
    try:
        process = subprocess.run(['RNAfold', '--noPS'], input=sequence, encoding='utf-8', capture_output=True,
                                 check=True)
        output = process.stdout.strip()
        mfe_line = output.splitlines()[-1]
        mfe_str = mfe_line.split()[-1]
        mfe_str = mfe_str.strip('()')
        mfe = float(mfe_str)
        return mfe
    except subprocess.CalledProcessError as exc:
        raise ValueError(f"RNAfold failed, error message: {exc.stderr}")
    except ValueError as e:
        raise ValueError(f"Error processing RNAfold output for sequence '{sequence}': {e}")


# Gaussian CDF definition
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)


# Extract gene name and transcript ID from filename
def extract_gene_transcript_id(filename):
    match = re.search(r'(\w+)_(\w+)_kmers', filename)
    if match:
        return match.group(1), match.group(2)
    else:
        return "Unknown", "Unknown"


# List to store the results
results = []

# Loop through each file in the directory
for filename in os.listdir(args.input):
    if filename.endswith("merged_normalized.csv"):
        filepath = os.path.join(args.input, filename)
        gene_name, transcript_id = extract_gene_transcript_id(filename)
        df = pd.read_csv(filepath)

        if 'kmer' not in df.columns:
            print(f"'kmer' column is missing in file {filename}. Skipping this file.")
            continue

        df['MFE'] = df['kmer'].apply(calculate_mfe)
        mini_shared_length = df.at[0, 'Mini_Shared_Length']
        mfe_data = df.groupby('MFE').agg({
            'Local_Frequency': 'sum',
            'Normalized_K-mer_Count': 'sum',
            'Count': 'sum'  # Including the Count column for aggregation
        }).reset_index()

        global_frequency = df.at[0, 'Global_Frequency']
        present_in_transcripts = df.at[0, 'Present_in_Transcripts']

        if len(mfe_data['MFE'].unique()) < args.threshold:
            sum_normalized_kmer_count = mfe_data['Normalized_K-mer_Count'].sum()
            sum_kmer_count = mfe_data['Count'].sum()  # Sum of the Count column
            results.append({
                'File': filename,
                'Gene_Name': gene_name,
                'Transcript_ID': transcript_id,
                'Global_Frequency': global_frequency,
                'Present_in_Transcripts': present_in_transcripts,
                'Mini_Shared_Length': mini_shared_length,
                'Sum or Fitted A (Abundance) for Normalized Count': f'{sum_normalized_kmer_count:.2f}',
                'Sum or Fitted A (Abundance) for Count': f'{sum_kmer_count:.2f}',
                'Fixed Mean (xc)': 'N/A',
                'Fixed Standard Deviation (w)': 'N/A',
                'Report': 'Insufficient Data'
            })
            continue

        mfe_data_sorted = mfe_data.sort_values(by='MFE')
        mfe_data_sorted['Cumulative_Local_Frequency'] = mfe_data_sorted['Local_Frequency'].cumsum()
        mfe_data_sorted['Cumulative_Normalized_Count'] = mfe_data_sorted['Normalized_K-mer_Count'].cumsum()
        mfe_data_sorted['Cumulative_Count'] = mfe_data_sorted[
            'Count'].cumsum()  # Calculating the cumulative sum for the Count column

        x_data = mfe_data_sorted['MFE']
        y_data_local = mfe_data_sorted['Cumulative_Local_Frequency']
        y_data_normalized = mfe_data_sorted['Cumulative_Normalized_Count']
        y_data_count = mfe_data_sorted['Cumulative_Count']  # y-data for the Count column
        initial_guesses_local = [min(y_data_local), max(y_data_local) - min(y_data_local), x_data.mean(), x_data.std()]

        try:
            popt_local, pcov_local = curve_fit(gaussian_cdf, x_data, y_data_local, p0=initial_guesses_local)
            A0_fitted_local, A_fitted_local, xc_fitted_local, w_fitted_local = popt_local

            if suitable_criteria_for_mfe(xc_fitted_local, w_fitted_local):
                initial_guesses_normalized = [min(y_data_normalized), max(y_data_normalized) - min(y_data_normalized)]
                popt_normalized, pcov_normalized = curve_fit(
                    lambda x, A0, A: gaussian_cdf(x, A0, A, xc_fitted_local, w_fitted_local),
                    x_data, y_data_normalized, p0=initial_guesses_normalized
                )
                A0_fitted_normalized, A_fitted_normalized = popt_normalized

                # Fit for the 'Count' data using the fixed Gaussian parameters
                initial_guesses_count = [min(y_data_count), max(y_data_count) - min(y_data_count)]
                popt_count, pcov_count = curve_fit(
                    lambda x, A0, A: gaussian_cdf(x, A0, A, xc_fitted_local, w_fitted_local),
                    x_data, y_data_count, p0=initial_guesses_count
                )
                A0_fitted_count, A_fitted_count = popt_count

                # Append successful fitting result
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Mini_Shared_Length': mini_shared_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': f'{A_fitted_normalized:.2f}',
                    'Sum or Fitted A (Abundance) for Count': f'{A_fitted_count:.2f}',
                    'Fixed Mean (xc)': f'{xc_fitted_local:.2f}',
                    'Fixed Standard Deviation (w)': f'{w_fitted_local:.2f}',
                    'Report': 'OK'
                })
            else:
                # Handle case where fitting criteria are not met
                sum_normalized_kmer_count = mfe_data['Normalized_K-mer_Count'].sum()
                sum_kmer_count = mfe_data['Count'].sum()
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Mini_Shared_Length': mini_shared_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': f'{sum_normalized_kmer_count:.2f}',
                    'Sum or Fitted A (Abundance) for Count': f'{sum_kmer_count:.2f}',
                    'Fixed Mean (xc)': 'N/A',
                    'Fixed Standard Deviation (w)': 'N/A',
                    'Report': 'Unsuitable Fit - Absolute mean value is smaller than one fifth standard deviation, did not fit'
                })
        except Exception as e:
            # Handle any exceptions during the fitting process
            # Calculate and report the sum of the individual values instead of N/A
            sum_normalized_kmer_count = mfe_data['Normalized_K-mer_Count'].sum()
            sum_kmer_count = mfe_data['Count'].sum()
            results.append({
                'File': filename,
                'Gene_Name': gene_name,
                'Transcript_ID': transcript_id,
                'Global_Frequency': global_frequency,
                'Present_in_Transcripts': present_in_transcripts,
                'Mini_Shared_Length': mini_shared_length,
                'Sum or Fitted A (Abundance) for Normalized Count': f'{sum_normalized_kmer_count:.2f}',
                'Sum or Fitted A (Abundance) for Count': f'{sum_kmer_count:.2f}',
                'Fixed Mean (xc)': 'N/A',
                'Fixed Standard Deviation (w)': 'N/A',
                'Report': f'Unsuitable Fit - {str(e)}'
            })

# Specify the column order with the additional fields for the non-normalized count
column_order = [
    'File',
    'Gene_Name',
    'Transcript_ID',
    'Global_Frequency',
    'Present_in_Transcripts',
    'Mini_Shared_Length',
    'Sum or Fitted A (Abundance) for Normalized Count',
    'Sum or Fitted A (Abundance) for Count',
    'Fixed Mean (xc)',
    'Fixed Standard Deviation (w)',
    'Report'
]

# Create results DataFrame with the specified column order
results_df = pd.DataFrame(results, columns=column_order)

# Save to CSV file specified by the command-line argument
results_df.to_csv(args.output, index=False)
