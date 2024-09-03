

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

Use of HISAT2 for Alignment of Reads with Reference Genome (the output sam file from hisat2 will be input to samtools, resulting in bam file.) [Download | HISAT2 (daehwankimlab.github.io)](https://daehwankimlab.github.io/hisat2/download/)

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

Use [GTF File](https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz) to Build the Feature Count Matrix

Building of Feature Count Matrix With Subread Tool

```
featureCounts -S 1 -a Homo_sapiens.GRCh38.106.gtf -o test_featurecounts.txt test_trimmed.bam
```

