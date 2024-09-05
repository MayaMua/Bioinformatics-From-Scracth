# ChIP Next Generation Sequencing Analysis

## 0. Workflow

#### 

#### Required software:

```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install macs3 -y
conda install fastqc -y
conda install bowtie2 -y
conda install samtools -y
conda install samtools -y
```



## 1. Quality assessment with FastQC

Generate a HTML report of ChIPseq library QC:

1. Run the command `original_fastq_data/fastqc fastinput.fastq.gz` . You could also use a wildcard command `fastqc *.gz`.

2. FastQC report will be stored under `input_fastqc.html`.



