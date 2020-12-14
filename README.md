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
java -jar trimmomatic.jar PE -threads 4 392_1.fastq.gz 392_2.fastq.gz forward_paried.fastq.gz \
forward_unparied.fastq.gz reverse_paired.fastq.gz reverse_upaired.fastq.gz \
ILLUMINACLIP:/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:10:30

```

### Read Mapping to Reference Chromosome
```
 bwa mem -t 8 -R "@RG\tID:rg1\tSM:foo" hg38_chr7 forward_paired.fastq.gz reverse_paired.fastq.gz > 392_aln.sam

```
> The subsequent data processing steps will be done using Picard and GATK according to the best practices outlined in the documentation.

### Picard Workflow
**Sam file validation step**
```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar ValidateSamFile -INPUT 392_aln.bam -MODE SUMMARY 

```
**Sorting the bam file**

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar SortSam -INPUT 392_aln.bam -OUTPUT 392_aligned_sorted.bam -SORT_ORDER coordinate

```

**Removing Duplicate Reads**

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar MarkDuplicates -INPUT 392_aligned_sorted.bam -OUTPUT 392_dup_marked.bam -METRICS_FILE 392_metrics.metrics

```

