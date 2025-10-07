

### 1.1. hg38

```bash
# Download the reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Unzip the file
gunzip hg38.fa.gz

# Index the reference genome using samtools
module load SAMtools/1.10-GCC-8.3.0
samtools faidx hg38.fa

```

### 1.2. hg19

```bash
# Download the reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Unzip the file
gunzip hg19.fa.gz

# Index the reference genome using samtools
module load SAMtools/1.10-GCC-8.3.0
samtools faidx hg19.fa
