#!/bin/bash

# set input and output directory
INPUT_DIR=$(dirname "${BASH_SOURCE[0]}")
OUTPUT_DIR="${INPUT_DIR}"

# loop through all fastq files in the input directory
for BAM_FILE in ${INPUT_DIR}/*.bam
do
    # get the filename without the extension
    FILENAME=$(basename ${BAM_FILE} .bam)
    
    # Get Feature Counts Matrix
    featureCounts -T 32 -S 2 -a Homo_sapiens.GRCh38.106.gtf -o ${OUTPUT_DIR}/${FILENAME}.txt ${BAM_FILE}
    
done

