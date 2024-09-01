

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



# Workflow

Quality Check of the Reads with FASTQC 

Run the command and two files are generated: `test_fastqc.html` and `test_fastqc.zip`

```bash
fastqc fasq_data/test.fastq
```



Use of Timmomatic Tool to Remove Poor Quality Reads

```bash
java -jar trimmomatic-0.39.jar SE -threads 4 fasq_data/test.fastq fasq_data/test_trimmed.fastq TRAILING:10 -phred33
```

