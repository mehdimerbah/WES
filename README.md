# Functional Genomics
## Gene Expression Profiling
Study for GEP Project: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0045200

## Whole Exome Sequencing

### Quality Control and Assessment

```
fastqc 392_1.fastq.gz
fastqc 392_2.fastq.gz
```



### Trimming

```
java -jar trimmomatic.jar PE -threads 4 392_1.fastq.gz 392_2.fastq.gz forward_paried.fastq.gz forward_unparied.fastq.gz reverse_paired.fastq.gz reverse_upaired.fastq.gz ILLUMINACLIP:/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:10:30

```