#!/bin/bash

#1. Set up required paths
SAMPLELIST="/path/to/samplelist.txt"
OUTPUT_DIR="/path/to/output"
TRIMMOMATIC_PATH="/path/to/trimmomatic"
HISAT2_INDEX="/path/to/hisat2/index"
STAR_INDEX="/path/to/STAR/index"
ADAPTERS="/path/to/adapters.fa"
REF_GENOME="/path/to/reference/genome"
CUFFLINKS_OUTPUT_DIR="/path/to/cufflinks_output"
STRINGTIE_OUTPUT_DIR="/path/to/stringtie_output"
RSEM_OUTPUT_DIR="/path/to/rsem_output"

# Always ensure that the paths are correctly assigned

# 2.Adapter Trimming with Cutadapt for different samples, the adaptor sequence is different.
while read sample
do
    # Replace raw_1.fq.gz and raw_2.fq.gz with the actual paired filenames from your dataset.
    cutadapt -j 5 \
        -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        -o ${OUTPUT_DIR}/${sample}_cutadapt_R1.fq \
        -p ${OUTPUT_DIR}/${sample}_cutadapt_R2.fq \
        ${sample}_R1.fq.gz ${sample}_R2.fq.gz
done < "$SAMPLELIST"
echo "Adapter trimming completed for all samples."

# 3.Quality Control with Trimmomatic
while read sample
do
    java -jar ${TRIMMOMATIC_PATH}/trimmomatic.jar PE -threads 5 \
        ${OUTPUT_DIR}/${sample}_cutadapt_R1.fq \
        ${OUTPUT_DIR}/${sample}_cutadapt_R2.fq \
        -baseout ${OUTPUT_DIR}/${sample}_cutadapt_trim.fq \
        LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20
done < "$SAMPLELIST"
echo "Trimmomatic quality trimming completed for all samples."

#The upper script is for cleaning raw data for GaussF model or other alignment free methods

#4. Alignment with HISAT2 or STAR
while read i
do
    name=$(echo ${i%R1_001.fastq})
    hisat2 -p 60 --dta -x ${HISAT2_INDEX} \
        -1 ${OUTPUT_DIR}/${sample}_cutadapt_trim_1P.fq \
        -2 ${OUTPUT_DIR}/${sample}_cutadapt_trim_2P.fq \
        -S ${OUTPUT_DIR}/${sample}_cutadapt_trim.sam
done < "$SAMPLELIST"
echo "HISAT2 alignment done for all samples."

#5. Conversion of SAM to BAM with Samtools
while read sample
do
    samtools sort -@ 10 -o ${OUTPUT_DIR}/${sample}_cutadapt_trim.bam \
        ${OUTPUT_DIR}/${sample}_cutadapt_trim.sam
done < "$SAMPLELIST"
echo "Samtools conversion to BAM and sorting completed for all samples."

# 6.Indexing BAM files with Samtools
# Convert SAM to BAM, sort, and index using Samtools
while read sample
do
    # Convert SAM to Sorted BAM
    samtools sort -@ 10 -o "${OUTPUT_DIR}/${sample}_sorted.bam" "${OUTPUT_DIR}/${sample}.sam"
    # Index the sorted BAM
    samtools index -@ 10 "${OUTPUT_DIR}/${sample}_sorted.bam"
done < "$SAMPLELIST"
echo "Samtools SAM to Sorted BAM conversion and indexing done for all samples."

#7. Cufflinks for transcript assembly
while read sample
do
    cufflinks -p 10 -o "${CUFFLINKS_OUTPUT_DIR}/${sample}" "${OUTPUT_DIR}/${sample}_sorted.bam"
done < "$SAMPLELIST"
echo "Cufflinks transcript assembly done for all samples."

# 8.StringTie for transcript assembly and quantification
while read sample
do
    stringtie -p 10 -o "${STRINGTIE_OUTPUT_DIR}/${sample}.gtf" "${OUTPUT_DIR}/${sample}_sorted.bam"
done < "$SAMPLELIST"
echo "StringTie transcript assembly and quantification done for all samples."

#9. Use `rsem-prepare-reference` if you haven't done so already.
RSEM_REF="/path/to/rsem/reference"

# Convert SAM to BAM, sort, and index using Samtools
# Calculate expression levels with RSEM
while read sample
do
    # Assuming HISAT2 alignment output is in SAM format, convert to BAM
    samtools view -@ 10 -Sb "${OUTPUT_DIR}/${sample}.sam" \
        | samtools sort -@ 10 -o "${OUTPUT_DIR}/${sample}_sorted.bam"
        
    # Index the sorted BAM
    samtools index -@ 10 "${OUTPUT_DIR}/${sample}_sorted.bam"
    
    # Run RSEM to calculate expression levels
    # Assuming paired-end data; for single-end, remove `--paired-end` and provide only one BAM file
    rsem-calculate-expression --bam --paired-end --no-bam-output -p 10 \
        "${OUTPUT_DIR}/${sample}_sorted.bam" \
        "${RSEM_REF}" \
        "${RSEM_OUTPUT_DIR}/${sample}"
done < "$SAMPLELIST"
echo "RSEM quantification done for all samples."


#10. Directory for Salmon and Kallisto output
SALMON_OUTPUT_DIR="/path/to/salmon_output"
KALLISTO_OUTPUT_DIR="/path/to/kallisto_output"

# Ensure output directories exist
mkdir -p "$SALMON_OUTPUT_DIR"
mkdir -p "$KALLISTO_OUTPUT_DIR"

# Path to transcriptome index for Salmon and Kallisto
# You have to prepare the index before running this pipeline using `salmon index` and `kallisto index`.
SALMON_INDEX="/path/to/salmon/index"
KALLISTO_INDEX="/path/to/kallisto/index"

# Quantify expression levels with Salmon
while read sample
do
    # Assuming 'trimmomatic' produced trimmed FASTQ files
    salmon quant -i "$SALMON_INDEX" -l A \
          -1 "${OUTPUT_DIR}/${sample}_cutadapt_trim_1P.fq" \
          -2 "${OUTPUT_DIR}/${sample}_cutadapt_trim_2P.fq" \
          -p 8 --validateMappings \
          -o "${SALMON_OUTPUT_DIR}/${sample}"
done < "$SAMPLELIST"
echo "Salmon quantification done for all samples."

# Quantify expression levels with Kallisto
while read sample
do
    # Assuming 'trimmomatic' produced trimmed FASTQ files
    kallisto quant -i "$KALLISTO_INDEX" -o "${KALLISTO_OUTPUT_DIR}/${sample}" \
          -b 100 \
          --threads=8 \
          "${OUTPUT_DIR}/${sample}_cutadapt_trim_1P.fq" \
          "${OUTPUT_DIR}/${sample}_cutadapt_trim_2P.fq"
done < "$SAMPLELIST"
echo "Kallisto quantification done for all samples."

# Clean temporary files
# Uncomment the following line if you wish to remove intermediate files to save space
# rm ${OUTPUT_DIR}/*.fq ${OUTPUT_DIR}/*.sam

