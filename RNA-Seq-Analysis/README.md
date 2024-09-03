

# Required Software

Necessary Programs In Linux Environment:

```bash
sudo apt-get -y install fastqc		
sudo apt-get -y	install	hisat2		
sudo apt-get -y install	subread		
sudo apt-get -y install	samtools	
```

fastQC will be used for quality check of the raw reads.

HISAT2 will be used for the alignment of the reads with reference genome.

Subread will be used for building of Count Matrix for DEG analysis.

[Timmomatic](http://www.usadellab.org/cms/?page=trimmomatic) tool will be used to perform trimming of low-quality reads and adapters.

Data: FASTQ Files (Containing Reads)

## Use Linux command lines

1. Quality check of the raw reads using FastQC
2. Alignment of the reads to the reference genome using HISAT2
3. Extraction of the reads mapping to the reference genome using SAMtools

## Analyze data using R:

Feature Count Matrix (Counting Reads) 

# Workflow under Linux Commands

1. Quality Check of the Reads with FASTQC . Run the command and two files are generated: `test_fastqc.html` and `test_fastqc.zip`

```bash
fastqc fasq_data/test.fastq
```

2. Use of Timmomatic Tool to Remove Poor Quality Reads

```bash
java -jar trimmomatic-0.39.jar SE -threads 4 fasq_data/test.fastq fasq_data/test_trimmed.fastq TRAILING:10 -phred33
```

3. Use of HISAT2 for Alignment of Reads with Reference Genome (the output sam file from hisat2 will be input to samtools, resulting in bam file.) [Download | HISAT2 (daehwankimlab.github.io)](https://daehwankimlab.github.io/hisat2/download/)

```
hisat2 -q --rna-strandness R -x grch38/genome -U fasq_data/test_trimmed.fastq | samtools sort -o fasq_data/test_trimmed.bam
```

`R: The reads are expected to be reversed stranded.`

`F: The reads are expected to be forward stranded.`

`U: The reads are UN stranded.`

`x: Reference genome.`

Output:

```
1249800 reads; of these: # 1.2 million reads which are unpaired present in a sequencing data.
  1249800 (100.00%) were unpaired; of these:
    86563 (6.93%) aligned 0 times # They are not aligned.
    1082380 (86.60%) aligned exactly 1 time
    80857 (6.47%) aligned >1 times
93.07% overall alignment rate
```

4. Use [GTF File](https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz) to Build the Feature Count Matrix with Subread Tool

```
featureCounts -S 1 -a Homo_sapiens.GRCh38.106.gtf -o test_featurecounts.txt test_trimmed.bam
```

Output:

- 1st column: **Gene id**
- 2nd: chromosome number.
- 3rd: the start position of this gene on the chromosome
- 4th:  end position
- 5th: strand
- 6th: the length of the transcript
- 7th: **the counts of the reads** of a particular gene present in the RNA SEQ data.

# Multipipe FastQ Files Processing Using Bash Scripts

There are 4 scripts used to process FastQ files.

script1.sh (Remove poor quality reads)

```bash
#!/bin/bash

# Set the path to the input directory as the current directory
input_dir="original_fastq_data"
output_dir="processed_fastq_data"
# Loop over all fq.gz files in the input directory
for file in "$input_dir"/*.fastq; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fastq)

    # Run the trimmomatic command on the file
    java -jar trimmomatic-0.39.jar SE -threads 32 "$file" "$output_dir/$filename"_trimmed.fastq TRAILING:10 -phred33 
done
```

script2.sh (Alignment)

```bash
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

```

script3.sh (Create count matrix)

```bash
#!/bin/bash

# set input and output directory
# INPUT_DIR=$(dirname "${BASH_SOURCE[0]}")
# OUTPUT_DIR="${INPUT_DIR}"

INPUT_DIR="processed_fastq_data"
OUTPUT_DIR="count_matrix"

# loop through all fastq files in the input directory
for BAM_FILE in ${INPUT_DIR}/*.bam
do
    # get the filename without the extension
    FILENAME=$(basename ${BAM_FILE} .bam)
    
    # Get Feature Counts Matrix
    featureCounts -T 32 -S 2 -a Homo_sapiens.GRCh38.106.gtf -o ${OUTPUT_DIR}/${FILENAME}.txt ${BAM_FILE}
    
done

```

script4.sh

```bash
#!/bin/bash

# Set the directory
DIR="count_matrix"
# Step 1: Extract the first and last columns from the first file and save to output.csv
FIRST_FILE=$(ls ${DIR}/*.txt | head -n 1)
awk '{print $1 "," $NF}' $FIRST_FILE > ${DIR}/output.csv

# Step 2: Loop through the remaining txt files and extract the last column, adding it to output.csv
for file in ${DIR}/*.txt; do
    if [ "$file" != "$FIRST_FILE" ]; then
        awk '{print $NF}' "$file" | paste -d, ${DIR}/output.csv - > ${DIR}/output2.csv
        mv ${DIR}/output2.csv ${DIR}/output_tmp.csv
    fi
done

# Step 3: Remove folder names from the columns in output.csv
awk -F',' '{
    for(i=1; i<=NF; i++) {
        sub(".*/", "", $i)
    }
    print $0
}' ${DIR}/output_tmp.csv > ${DIR}/output.csv

# Step 4: Remove the temporary file
rm ${DIR}/output_tmp.csv
```

