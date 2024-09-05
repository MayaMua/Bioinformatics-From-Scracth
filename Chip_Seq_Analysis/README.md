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

1. Run the command `original_fastq_data/fastqc fastinput.fastq.gz` . You could also use a wildcard command `fastqc *.gz`.

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

```
unzip hg19.zip -d hg19
bowtie2 -p 3 -q --local -x hg19/hg19 -U Input.fastq.gz -S Input_aln_unsorted.sam
```

The basic options for aligning reads to the genome using `Bowtie2` are:

- `-p`: number of processors / cores

This you may have decided to run with number different than 3 as in the command above. The more cores you allocate the quicker should you be able to finish the alignment procedure. I recommend to leave our 1 or 2 cores for the OS. This means that if you have 10 CPUs, use `-p 8`

- `-q`: indicates reads are in FASTQ format
- `--local`: local alignment feature to perform [soft-clipping](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/#:~:text=Aligners reaching such high mapping,for the alignment as such.)
- `-x`: /path/to/genome_indices_directory *and* genome index base name

Be aware here. It isn;t enough to specify path to the genome repository. You also have to give it a genome base name. That is why in the command above we have `-x hg19/hg19` rather than `-x hg19` which won't work!

- `-U`: /path/to/FASTQ_file (this can me compressed)
- `-S`: /path/to/output/SAM_file (inc. of file *name*)



