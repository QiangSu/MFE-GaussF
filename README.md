# MFE-GaussF
MFE-GaussF Transcript Quantification Pipeline

## Overview:

The MFE-GaussF pipeline represents a cutting-edge approach to accurately quantify transcript abundance at the isoform level from RNA-seq data. This algorithm capitalizes on distinctive k-mer signatures to address prevalent biases in RNA-seq analyses. By implementing a parametric Gaussian model, it performs sophisticated bias correction, ensuring accurate k-mer counting, as well as TPM (Transcripts Per Million) or RPKM (Reads Per Kilobase of transcript, per Million mapped reads) estimates. This enhancement allows for more detailed and reliable transcriptomic analyses.

```bash
pip install MFE-GaussF
```

***
### **Step 1: `minimal-shared region filtering`**

kmer_frequency_distribution_mini_shared.py (isoform_unique_sequence-loop.py)

This tool processes a FASTA file containing transcript sequences and outputs a set of CSV files that summarize the k-mer content for each transcript. Each CSV file contains a list of k-mers of specified length that are present in the transcript, along with their local and global frequency, and information on which transcripts each k-mer is present in.

Features
K-mer Counting: For a given transcriptome FASTA file, count all k-mers of a specified length (default is set to 50).<br>
Minimal Shared K-mers Output: For each transcript, output the k-mers that have the minimum global frequency—the smallest number of transcripts in which the k-mer appears.<br>
CSV Output Content: Generate a CSV file for each isoform with the following columns:<br>

'kmer': The k-mer sequence.<br>
'Local_Frequency': The number of times the k-mer appears in the specific isoform.<br>
'Global_Frequency': The number of transcripts that contain the k-mer across the entire transcriptome.<br>
'Present_in_Transcripts': A list of transcript identifiers that share the k-mer if its global frequency is more than 1. For unique k-mers, the identifier of the single transcript is given.

Installation

To install the specific version (0.1.1) of the package `minimal-shared-kmers` using pip, run the following command in your terminal:

```bash
pip install MFE-GaussF
```

vim analyze_kmers.py
```
import argparse
from MFE_GaussF import kmer_frequency_distribution_mini_shared

# Set up argument parsing
parser = argparse.ArgumentParser(description="Process some integers.")
parser.add_argument('--input', required=True, help='Input file path')
parser.add_argument('--output', required=True, help='Output file path')

args = parser.parse_args()

# Assuming 'my_function' takes two arguments: input path and output path
result = kmer_frequency_distribution_mini_shared.my_function(input=args.input, output=args.output)
```
running minimal-shared script
```python
python analyze_kmers.py --input "./kmer_isoform_gene/mart_export_fasta_iosform.txt" --output "./kmer_isoform_gene/isoform_minimal_shared_kmers"
```

Or other Usage
To use this tool, you need to have Python installed on your system. The script requires a FASTA file with the transcript sequences as input and a directory path where the CSV files will be saved as output.

Execute the script with the necessary arguments from the command line. For example:<br>
```python
python kmer_frequency_distribution_mini_shared.py --input path/to/your/ACTB_reference/mart_export_ACTB.txt --output path/to/output/directory/
```
Command-Line Arguments<br>
--input: Path to the input FASTA file containing transcript sequences (https://useast.ensembl.org/biomart/martview/aeb3390f02325ab7951be9a7d6daaa42).<br> 
--output: Path to the output directory where CSV files for each transcript will be saved.

Output File Details
For each transcript in the input FASTA file, the script will create a corresponding CSV file in the output directory with a name derived from the transcript header, sanitized to be filesystem-friendly.

In the output CSV files for each transcript, only k-mers that have the smallest global frequency for that transcript are included. If multiple k-mers share the same smallest global frequency, then all such k-mers are included in the CSV file. The 'Present_in_Transcripts' field in the CSV may include multiple transcript names, indicating that those transcripts share the k-mer.

If the global frequency of a k-mer is 1, indicating that it is unique to a single transcript, then the 'Present_in_Transcripts' field will only contain the identifier of that specific transcript.
***
### **Step 2.1: `k-mer Counting and Normalization`**

> [!important]  
Kmer Counter(recommended)

Overview
Kmer Counter is a high-performance C++ application designed to rapidly count k-mer frequencies in large genomic datasets. Implementing optimized algorithms and leveraging the power of multi-threading, this tool offers significant speed enhancements over traditional k-mer counting techniques typically found in higher-level languages.

The tool reads k-mers from CSV files and processes gzipped FASTQ files to provide comprehensive k-mer counts. It has been carefully engineered to handle large sequencing data efficiently, making the best use of modern multi-core processors.

Features<br>
High-throughput k-mer counting from gzipped FASTQ files.<br>
Multithreading support to take advantage of multiple CPU cores.<br>
Atomic operations for thread-safe counting without significant locking overhead.<br>
Streamlined I/O with buffered reads and writes for optimal performance.<br>
User-friendly command-line interface for easy integration into genomics workflows.<br>

Dependencies
zlib: For reading gzipped FASTQ files.<br>
C++17 filesystem library: For convenient file and directory manipulation.<br>
POSIX threads (pthreads): For multi-threading support.<br>
Usage
To use Kmer Counter, compile the source code and run the resulting binary with the required arguments specifying the k-mer CSV directory, the input FASTQ gzipped file, the output directory for counts, and the number of threads.
```cpp
./kmer_counter --kmer_dir <kmer_csv_dir> --fastq_file <fastq_gz_file> --output_dir <output_csv_dir> --threads <number_of_threads>
```
// Compile this with:<br>
// g++ -std=c++17 -o kmer_counter kmer_counter.cpp -lz -lpthread

Kmer counting by python<br>

To install the package `MFE-GaussF` using pip, run the following command in your terminal:
```bash
pip install MFE-GaussF
```
vim kmer_counting.py
```
from MFE_GaussF import kmer_counting_loop
import argparse
import subprocess

def main():
    # Initialize Argument Parser
    parser = argparse.ArgumentParser(description="Wrapper for running k-mer counting on fastq files.")
    parser.add_argument('--k', type=int, required=True, help="K-mer size.")
    parser.add_argument('--threads', type=int, required=True, help="Number of threads to use.")
    parser.add_argument('--chunk_size', type=int, required=True, help="Chunk size for processing.")
    parser.add_argument('--fastq', required=True, help="Path to the fastq.gz file.")
    parser.add_argument('--kmer_dir', required=True, help="Directory for k-mer output.")
    parser.add_argument('--output', required=True, help="Output directory for final results.")
    
    args = parser.parse_args()
    
    # Construct command
    command = [
        'python', 'kmer_counting_loop.py', 
        '--k', str(args.k),
        '--threads', str(args.threads),
        '--chunk_size', str(args.chunk_size),
        '--fastq', args.fastq,
        '--kmer_dir', args.kmer_dir,
        '--output', args.output
    ]
    
    # Execute command
    subprocess.run(command)

if __name__ == '__main__':
    main()
```
running counting script<br>
```python
python kmer_counting.py --k 50 --threads 30 --chunk_size 10000000 --fastq /path/to/your.fastq.gz --kmer_dir /path/to/kmer_output --output /path/to/final_output
```

kmer_counting_loop.py

Introduction
This tool is designed for bioinformatics analysis to count k-mer frequencies in sequencing data stored in FASTQ format. It is particularly useful when dealing with large datasets as it leverages Python's multiprocessing capabilities for parallel processing, thus enhancing performance and reducing computation time.
Features
Count specified k-mer sizes in FASTQ files (compressed with gzip).<br>
Use input CSV files containing k-mer sequences to filter and count only relevant k-mers.<br>
Handle large datasets efficiently with chunk-based parallel processing.<br>
Utilize multiple CPU cores for faster computation.<br>
Generate output CSV files containing the count of each k-mer.<br>

Or other usage
To use this tool, you need to provide several command-line arguments. Here is the syntax for running the script:
```python
python kmer_counting_loop.py --k <kmer_size> --chunk_size <chunk_size> --fastq <fastq_file_path> --kmer_dir <kmer_directory> --output <output_directory> [--threads <number_of_threads>]
```
Command-Line Arguments<br>
--k: Size of the k-mer you wish to count (required).<br>
--chunk_size: Number of records from the FASTQ file to be processed in each parallel chunk (required).<br>
--fastq: Path to the compressed FASTQ file that contains the sequencing data (required).<br>
--kmer_dir: Directory path containing input CSV files with k-mer sequences; each CSV must have k-mers listed in the first column (required).<br>
--output: Path to the output directory where the k-mer count CSV files will be saved (required). The output directory should be same as the input directory (kmer_dir)!!!<br>
--threads: Number of threads to use for processing; defaults to the number of CPU cores available on the system.<br>
Example usage:
python kmer_counting_loop.py --k 50 --threads 30 --chunk_size 10000000 --fastq /path/to/data.fastq.gz --kmer_dir /path/to/directory --output /path/to/directory

Output
The script will output CSV files in the specified output directory, with each file named according to the original k-mer CSV file but appended with _kmer_counts. For example, if there is an input file named sample_kmers.csv, the output file will be sample_kmers_kmer_counts.csv. Each output file will contain two columns: K-mer and Count, where Count is the frequency of that k-mer in the FASTQ file.
Performance
This tool is designed to handle large FASTQ files efficiently. By using parallel processing, the script splits the FASTQ file into chunks and processes each chunk in a separate CPU core, speeding up the counting operation significantly. The time taken will be printed at the end of the execution for each k-mer count task.


### **Step 2.2: `K-mer Counts Merging and Normalization`**

To install the package `MFE-GaussF` using pip, run the following command in your terminal:
```bash
pip install MFE-GaussF
```
vim kmer_count_normalizing.py
```
from MFE_GaussF import merge_normalize_isoform_count_v1
from MFE_GaussF import merge_normalized_isoform_count_TPM

import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Wrapper script for count normalization tasks.")
    parser.add_argument('--directory', required=True, help="Directory for processing.")
    parser.add_argument('--output_directory', required=True, help="Output directory for processed files.")
    parser.add_argument('--read_length', type=int, required=True, help="Read length.")
    parser.add_argument('--k', type=int, required=True, help="K-mer size.")
    parser.add_argument('--fastq', help="Optional: Fastq file path for isoform count v1 script.")
    args = parser.parse_args()

    if args.fastq:
        # If a fastq file is provided, execute merge_normalize_isoform_count_v1
        subprocess.run([
            'merge_normalize_isoform_count_v1',
            '--directory', args.directory,
            '--fastq', args.fastq,
            '--output_directory', args.output_directory,
            '--read_length', str(args.read_length),
            '--k', str(args.k)
        ])
    else:
        # Otherwise, execute merge_normalized_isoform_count_TPM
        subprocess.run([
            'merge_normalized_isoform_count_TPM',
            '--directory', args.directory,
            '--output_directory', args.output_directory,
            '--read_length', str(args.read_length),
            '--k', str(args.k)
        ])

if __name__ == '__main__':
main()

```
running count normalizing script<br>
for TPM<br>
```python
python count_normalizing.py --directory ./ --output_directory ./merge_data_all1 --read_length 150 --k 50 --fastq /path/to/fastq.gz
```
for RPKM<br>
```python
python count_normalizing.py --directory ./gene_folder/example --output_directory ./gene_folder/example/12 --read_length 150 --k 50
```


merge_mormalize_isoform_count_v1.py

This script is designed to further process the output of a previous k-mer counting script. Its purpose is to merge the k-mer count data into the original k-mer CSV files and to normalize these counts to account for differences in the total number of k-mers and read counts. This is a necessary step in many bioinformatics workflows, particularly those involving comparative genomics or quantitative assessment of sequence representation.

Features<br>
Merges k-mer count data with the original k-mer list CSV files.<br>
Normalizes k-mer frequencies using the total k-mer counts and read lengths.<br>
Supports input from gzipped FASTQ files for read count determination.<br>
Efficiently calculates normalization factors and processes large datasets.<br>

Example usage:
This script accepts command-line arguments to specify the input and output directories, the FASTQ file path, the read length, and the k-mer size. Here's how to run the script:
For TPM
```python
python ./scripts/merge_normalized_isoform_count_TPM.py --directory ./data/input --output_directory ./data/output --read_length 150 --k 50
```
For RPKM
```python
python merge_normalize_isoform_count_v1.py --directory <input_directory as the output directory of last kmer_counter.py script> --output_directory <output_directory new directory> --fastq <path_to_FASTQ.GZ> --read_length 150 --k 50
```
Command-Line Arguments
--directory: The directory containing the *_kmers.csv and corresponding *_kmer_counts.csv files (required). This directory is same as the output directory from the last script (kmer_counting_loop.py).
--output_directory: The directory where the merged and normalized CSV files will be saved (required). The output directory should be to a new directory for further MFE-GaussF workflow.
--fastq: The path to the gzipped FASTQ file for which k-mer counts were computed (required).
--read_length: The length of the reads in the FASTQ sequences, necessary for normalization (default is 150).
--k: The length of the k-mers used during the counting process (default is 50).
Output

> [!note] 
> For each *_kmers.csv file in the input directory, the script will save a corresponding *_merged_normalized.csv file in the output directory. This file will contain the original k-mer data, the raw count, and > an additional column with normalized k-mer counts.
***
### **Step 3: `Gaussian CDF Fitting`**

pipeline_abundance_GaussF_esti_loop.py

Introduction
This Python script is designed to analyze MFE distribution in sequence data and estimate the sequence abundance by fitting a cumulative distribution function (CDF) of a Gaussian to the MFE profile. It serves as a post-processing tool following k-mer counting, allowing researchers to derive meaningful biological insights based on the MFE composition and k-mer abundance patterns.

Features<br>
Analyzes the MFE of sequences represented by k-mers.<br>
Performs fitting of a Gaussian CDF to the sum of normalized k-mer counts grouped by MFE percentage.<br>
Extracts gene and transcript information from the input CSV filenames.<br>
Produces structured output for quick assessment of fit success and estimated parameters.<br>
Offers flexibility through user-defined minimum thresholds for k-mer counts appropriate for fitting.

Installation

To install the package `MFE-GaussF` using pip, run the following command in your terminal:
```bash
pip install MFE-GaussF
```

vim abundance_estimation_TPM.py
```
import argparse
from MFE_GaussF import MFE_GaussF_esti

def estimate_abundance(input_file, output_file, threshold):
    # Load your data
    # Assuming the input file is a CSV, adjust according to your needs
    data = pd.read_csv(input_file)
    
    # You would replace the following lines with the actual estimation logic provided by abundance-GaussF-esti
    # For example, if the package requires you to initialize an estimator and fit the data
    estimator = SomeEstimatorClass(threshold=threshold)
    results = estimator.estimate(data)  # Pretend 'estimate' is a method to perform the estimation - adjust according to actual package capabilities
    
    # Save the results to the specified output file, again assuming CSV as output
    results.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run abundance estimation.')
    parser.add_argument('--threshold', type=int, help='Threshold for estimation', required=True)
    parser.add_argument('--input', type=str, help='Input file path', required=True)
    parser.add_argument('--output', type=str, help='Output file path', required=True)
    
    args = parser.parse_args()
    
    # Call the estimation function with parsed CLI arguments
estimate_abundance(args.input, args.output, args.threshold)
```
running minimal-shared script
```python
python abundance_estimation_TPM.py --input "./normalized_count/merged_data/" --output "./target/directory/"
```

vim abundance_estimation_RPKM.py
```
import argparse
from MFE_GaussF import pipeline_abundance_GaussF_esti_loop

def estimate_abundance(input_file, output_file, threshold):
    # Load your data
    # Assuming the input file is a CSV, adjust according to your needs
    data = pd.read_csv(input_file)
    
    # You would replace the following lines with the actual estimation logic provided by abundance-GaussF-esti
    # For example, if the package requires you to initialize an estimator and fit the data
    estimator = SomeEstimatorClass(threshold=threshold)
    results = estimator.estimate(data)  # Pretend 'estimate' is a method to perform the estimation - adjust according to actual package capabilities
    
    # Save the results to the specified output file, again assuming CSV as output
    results.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run abundance estimation.')
    parser.add_argument('--threshold', type=int, help='Threshold for estimation', required=True)
    parser.add_argument('--input', type=str, help='Input file path', required=True)
    parser.add_argument('--output', type=str, help='Output file path', required=True)
    
    args = parser.parse_args()
    
    # Call the estimation function with parsed CLI arguments
estimate_abundance(args.input, args.output, args.threshold)
```
running minimal-shared script
```python
python abundance_estimation_RPKM.py --input "./normalized_count/merged_data/" --output "./target/directory/"
```


Or other usage:
```
python pipeline_abundance_GaussF_esti_loop.py --threshold 5 --input /path/to/merge_data --output / path/to/merge_data/results_file.csv
```
Command-Line Arguments<br>
--input: The path to the input folder containing the k-mer CSV files where each file should have a filename format including gene and transcript IDs (e.g., GENE_ENST00001234567_kmers.csv) (required).<br>
--output: The full path and name of the output CSV file where the results will be saved (required).<br>
--threshold: The minimum number of k-mers required for performing fitting; the default value is 10 if not specified.<br>
Output
The script will output a CSV file containing the following columns:

> [!note] 
File: The name of the input CSV file that was processed to generate the corresponding row in the output file. It is used to trace back the results to the input data.<br>

Gene_Name: This is the name of the gene that the k-mers are associated with, typically extracted from the filename of the input CSV file according to a predetermined naming convention.<br>

Transcript_ID: The identifier for the specific transcript from which the k-mers were derived. Like the gene name, this is also extracted from the filename of the input CSV file.<br>

Global_Frequency: The frequency of the k-mer's occurrence across all transcripts in the dataset. This gives an idea of how common a particular k-mer sequence is overall.<br>

Present_in_Transcripts: An identifier indicating which transcripts include the k-mer. This can be a single transcript ID or a list of IDs, depending on k-mer representation in the data.<br>

Mini_Shared_Length: The minimum shared length between the input k-mer sequence and any of the transcripts. This value provides insight into the minimum overlap a k-mer has with known transcripts.<br>

Sum or Fitted A (Abundance) for Normalized Count: For each k-mer fitting, this field either contains the sum of the normalized k-mer counts (if curve fitting fails or is not applicable) or the value 'A' from the successfully fitted Gaussian cumulative distribution function, which represents the abundance of the k-mer after normalization for transcript length.<br>

Sum or Fitted A (Abundance) for Count: Similar to the above field, but for the raw k-mer count data. It contains the sum total of the raw counts (if curve fitting fails or is not applicable) or the value 'A' from the fitted Gaussian cumulative distribution function, indicating the overall abundance of the k-mer before normalization.<br>

Fixed Mean (xc): The mean (or center) of the k-mer distribution, denoted by 'xc', as estimated from the Gaussian CDF fitting process. It is fixed based on an initial fitting of the local frequency data and used for subsequent fittings. If fitting was not performed, this field will be 'N/A'.<br>

Fixed Standard Deviation (w): The standard deviation of the k-mer distribution, denoted by 'w', as estimated from the Gaussian CDF fitting process. It describes the spread or dispersion of the distribution. Similar to the fixed mean, this value is determined from an initial fit and used consistently for subsequent data. If fitting was not performed or failed, this field will be 'N/A'.<br>

Report: A text field containing messages about the status of the data processing and any curve fitting processes. It can include messages such as 'OK' to indicate successful processing, 'Insufficient Data' if there isn't enough data to perform the fitting, or a detailed error message if fitting failed.<br>


The interpretation of other files:

The 'human_pc_gene.tsv' file contains an annotated list of human protein-coding genes, each associated with a unique sequence region at the isoform level. This dataset comprises a total of 20,818 distinct genes derived from Homo sapiens, systematically cataloged to facilitate research on genetic variation and isoform expression.

The 'hm_mrna_isoform_summary.txt' details a comprehensive inventory of isoforms, including <br>
Transcript ID: Extracted from the filename by removing the '_kmers.csv' part, representing the gene name and transcript ID.<br>
Number of Kmers: The total count of unique kmers found in each corresponding CSV file.<br>
Global_Frequency: The frequency of the kmers across all sequences as found in the first data row of each corresponding CSV file.<br>
Present_in_Transcripts: A list of transcripts in which the kmers are present as found in the first data row of each corresponding CSV file.

The summary encompasses 175,572 isoforms subjected to unique region screening. 

We welcome contributions to this project and encourage users to submit issues or pull requests on our GitHub repository. Additionally, if you have any questions or comments about this script, please feel free to contact us via email.

