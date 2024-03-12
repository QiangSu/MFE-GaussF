MFE-GaussF Transcript Quantification Pipeline for GAPDH 11 isoforms

### **Step 1: `minimal-shared region filtering`**
Installation

To install the specific version (0.1.1) of the package `MFE-GaussF` using pip, run the following command in your terminal:

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
python analyze_kmers.py --input "./mart_export_GAPDH.txt" --output "./kmer_isoform_gene/GAPDH_isoform"
```

### **Step 2.1: `k-mer Counting`**


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
python kmer_counting.py --k 50 --threads 30 --chunk_size 10000000 --fastq /path/to/your.fastq.gz --kmer_dir ./kmer_isoform_gene/GAPDH_isoform --output ./kmer_isoform_gene/GAPDH_isoform
```

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
python count_normalizing.py --directory ./kmer_isoform_gene/GAPDH_isoform --output_directory ./GAPDH_merge_data_all --read_length 150 --k 50 --fastq /path/to/fastq.gz
```
for RPKM<br>
```python
python count_normalizing.py --directory ./kmer_isoform_gene/GAPDH_isoform --output_directory ./GAPDH_merge_data_all --read_length 150 --k 50
```

### **Step 3: `Gaussian CDF Fitting`**

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
python abundance_estimation_TPM.py --input "./GAPDH_merge_data_all" --output "./target/directory/"
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
python abundance_estimation_RPKM.py --input "./GAPDH_merge_data_all/" --output "./target/directory/"
```

