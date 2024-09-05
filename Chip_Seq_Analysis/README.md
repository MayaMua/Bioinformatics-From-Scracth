# ChIP Next Generation Sequencing Analysis

## 0. Workflow

#### 

#### Required software:

```bash
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

1. Run the command `fastqc original_fastq_data/fastinput.fastq.gz`. You could also use a wildcard command `fastqc *.gz`.

2. FastQC report will be stored under `input_fastqc.html`.

## 2. Alignment

### Alignment to the genome: creating a Bowtie2 index

```bash
# If you have hg.19.fa file, create indexed genome.
bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>
wget -c https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
mkdir ./hg19
bowtie2-build --threads 5 hg19.fa ./hg19/hg19
```

Or use other versions of reference human or other species genomes (all Bowtie and Bowtie2 index compatible) and go to the [index zone](https://benlangmead.github.io/aws-indexes/bowtie), skipping building your indexes by relying on these resources.  

```bash
unzip hg19.zip -d hg19
bowtie2 -p 5 -q --local -x hg19/hg19 -U Input.fastq.gz -S Input_aln_unsorted.sam
```

The basic options for aligning reads to the genome using `Bowtie2` are:

- `-p`: number of processors / cores

Use `lscpu` to check the number of CPU cores. Leave 1 or 2 cores for the OS. This means that if you have 10 CPUs, use `-p 8`.

- `-q`: indicates reads are in FASTQ format
- `--local`: local alignment feature to perform [soft-clipping](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/#:~:text=Aligners reaching such high mapping,for the alignment as such.)
- `-x`: /path/to/genome_indices_directory *and* genome index base name

A genome base name is also required when we specify the path `-x hg19/hg19` rather than `-x hg19` which won't work!

- `-U`: /path/to/FASTQ_file (this can be compressed)
- `-S`: /path/to/output/SAM_file (inc. of file *name*)



## 3. Filtering

### SAM/BAM Explanation

[Sequence Alignment/Map format and SAMtools | Bioinformatics | Oxford Academic (oup.com)](https://academic.oup.com/bioinformatics/article/25/16/2078/204688?login=false)

*To include a multiple mapped read in the analysis, or not to include - that is the question?*

One of the key questions that a ChIPseq data analyst needs to consider is whether to include multiple mapped reads, i.e. reads mapped to multiple locations on the reference genome. According to [this](https://pubmed.ncbi.nlm.nih.gov/21779159/) paper, allowing for non-uniquely mapped reads can boost the sensitivity of peak detection but with an added cost of enhanced false positive rate. Therefore practitioners would need to contemplate whether:

- It is more important to detect as many true peaks as possible (high [sensitivity](https://en.wikipedia.org/wiki/Sensitivity_and_specificity)), or  
- It is more important to avoid reporting false peaks as much as possible (high [specificity](https://en.wikipedia.org/wiki/Sensitivity_and_specificity))

### Changing file format from SAM to BAM

SAM format is human readable. But it is very big in size (in our case for a library size of 2.1GB the alignment file is 12GB), and not useable for downstream analyses. Convert alignment file from SAM format to a binary BAM format using `samtools` utility.

```bash
samtools view -h -S -b -o Input_aln_unsorted.bam Input_aln_unsorted.sam
```

In terms of the options used in the command, here is their explanation:

- `-h`: include header in output
- `-S`: input is in SAM format
- `-b`: output BAM format
- `-o`: /path/to/output/file *AND* file name

[samtools(1) manual page (htslib.org)](https://www.htslib.org/doc/1.2/samtools.html)
