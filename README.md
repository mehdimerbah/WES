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

**Editing the Read Group Name to add Platform**

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar AddOrReplaceReadGroups  -I 392_aligned_sorted.bam  -O 392_aligned_sorted_RGcorr.bam  RGID=8  RGLB=lib1  RGPL=ILLUMINA  RGPU=unit1  RGSM=20 

```

**Marking Duplicates on Edited File**

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar MarkDuplicates -INPUT 392_aligned_sorted_RGcorr.bam -OUTPUT 392_duplicate_marked.bam -METRICS_FILE markDuplicate_metrics.metrics

```
**Base Quality Score Recalibration**
First we build the model:
```
gatk BaseRecalibrator -I 392_duplicate_marked.bam -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa --known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf -O recalibrated_data.table 
 
```
**Applying BQSR**
```
gatk ApplyBQSR -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -I 392_duplicate_marked.bam --bqsr-recal-file recalibrated_data.table -O 392_recalibrated.bam 
```
**Second Pass Recalibration**
```
gatk BaseRecalibrator -I 392_recalibrated.bam -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa --known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf -O secondPass.table
```
**Covariate analysis: Before and After BQSR**
```
 gatk AnalyzeCovariates -before recalibrated_data.table -after secondPass.table -plots covAnalysis.pdf

```

**Calling Variants**
```
gatk --java-options "-Xmx4g" HaplotypeCaller -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -I 392_recalibrated.bam -O 392_GvarCall.g.vcf.gz -ERC GVCF

```

```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_GvarCall.g.vcf.gz -O 392_varCall.vcf.gz

```
**Filtering the Variants**

```
gatk --java-options "-Xms5g -Xmx15g" VariantFiltration -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_varCall.vcf.gz -O 392_varCall_filtered.vcf.gz --filter-name "lowGQ"     --filter-expression "GQ < 20.0"     --filter-name "lowMQ"     --filter-expression "MQ < 40.0"     --filter-name 'lowQD'  --filter-expression "QD < 2.0"     --filter-name "lowMQRankSum"     --filter-expression "MQRankSum < -12.5" 

```
**Selecting the Variants**

```
gatk --java-options "-Xms5g -Xmx15g" SelectVariants -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_varCall_filtered.vcf.gz -O 392_varCall_filteredExcluded_nonVar.vcf.gz --exclude-filtered true --exclude-non-variants true 

```



