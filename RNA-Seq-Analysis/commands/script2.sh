#!/bin/bash

# set input and output directory
# INPUT_DIR=$(dirname "${BASH_SOURCE[0]}")
# OUTPUT_DIR="${INPUT_DIR}"

INPUT_DIR="processed_fastq_data"
OUTPUT_DIR="${INPUT_DIR}"

# loop through all fastq files in the input directory
for FASTQ_FILE in ${INPUT_DIR}/*.fastq
do
    # get the filename without the extension
    FILENAME=$(basename ${FASTQ_FILE} .fastq)
    
    # align the reads and sort the resulting BAM file
    hisat2 -q --rna-strandness R -x grch38/genome -U ${FASTQ_FILE} -p 32 | samtools sort -o ${OUTPUT_DIR}/${FILENAME}.bam
    
done