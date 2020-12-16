# Whole Exome Sequencing Project

This project aims to construct a full pipeline for the identification and annotation of genetics variants given Paired End sequencing data: `392_1.fastq.gz` and `392_2.fastq.gz` as the forward and reverse reads respectively.

### Basic File Exploration
From the file we could first notice that the file uses Phred33 character encoding for the quality scores and so we can conclude that the Illumina version used must be v 1.8 or later.

`@A00721:81:HNLHYDSXX:1:1101:24361:1877 1:N:0:GCCGGACA+TGTAAGAG`

`+`

`CCTCCCACCTCAGGTCTGTGTTCATGCTCTTTTATCCCCATAACTAAAACTCTCTTCAAACT`

`::FF,FF,,,,FFF,FFF:FF:FFFF,FF:FFF,F,,:F,:,F:,,,::,FFF:F,,,F,,,`

From an overall view of both fastq files, it seems that the quality of the reads is very good with an average score well above 36. The subsequebt quality assessment steps should confirm this.
### Quality Control and Assessment
We will use FastQC to assess the quality of the data and generate our quality reports.


```
fastqc 392_1.fastq.gz
fastqc 392_2.fastq.gz
```
**Basic Statistics**  

![Image of Stats](/WES_Workflow/images/BasicStats.png "Basic Stats")

**Per Tile Quality**  
The per tile sequence quality is really good and so there are no specific tile reads to remove.  
![Tiles](/WES_Workflow/images/Tiles.png "Tiles")

### Trimming
To make sure we are only working with the pure exome data, we need to trim the adapters that were used in the sequencing procedure. 
We don't know anything about the experimental protocol used to generate the data and so we have no information about the type type of adapters that hybridized to our DNA fragments. 

I first assumed that Nextera-PE adapters were used but another pass of FastQC on the trimmed data revealed that the adapters were still present. Adapters that remained were Illumina Universal Adapters.
```
java -jar trimmomatic.jar PE -threads 4 392_1.fastq.gz 392_2.fastq.gz forward_paried.fastq.gz \
forward_unparied.fastq.gz reverse_paired.fastq.gz reverse_upaired.fastq.gz \
ILLUMINACLIP:/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:36

```   
![Image](/WES_Workflow/images/AdaptersTrim1.png "AT1")


And so I did a second trimming pass with the TruSeqPE adapters as a referce adapter file and managed to eliminate a big portion of the adapters. This did leave a small percentage of adapter content but FastQC report validated this as passable so I carried on with the analysis. 


```
java -jar trimmomatic.jar PE -threads 4 392_1.fastq.gz 392_2.fastq.gz forward_paried.fastq.gz \
forward_unparied.fastq.gz reverse_paired.fastq.gz reverse_upaired.fastq.gz \
ILLUMINACLIP:/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 MINLEN:36

```
![Image](/WES_Workflow/images/AT2.png "AT2")


**Trimmomatic Options and Arguments**  
In running `trimmomatic` I used mainly 2 trimming options:  
`ILLUMINACLIP`: This is to trim the adapter sequences given the TruSeq3 adapter fasta file. This in itself specifies the following column sperated arguments:
- `fastaAdapters`: The fasta file with adapters
- `seedMismatches`: Looking for matches to the adapter sequence and allowing only 2 mismatches in the seed.
- `palindromeClipThreshold`: threshold score for alignment of palindromic matches to adapter sequences. I picked 30 as recommended by Trimmomatic docs.
- `simpleClipThreshold`: threshold for simple matching. usually between 7 and 15 are recommended so I took 10 as an average.
      
`MINELN`: a minimum read length of 36 for all the reads after trimming. I chose to stick with 36 as a baseline minimum length so as not to eliminate possibly good reads.

### Read Mapping to Reference Chromosome
We now can safely map our reads back into a reference genome, in our case we will b using Chromosome 7 of hg38. The chromosome fasta files had been concatinated.
To use it the reference we have to index it using `bwa index` that uses the *Burrows Wheeler Transform* and *Smith Waterman* specified in the `-a` algorithm argument.    
`bwa index -p hg38_chr7 -a bwtsw hg38_chrom7.fa`   
Now we are ready to align our reads to the reference chromosome. For the RG argument, I basically used a dummy read group as I do not know about the experiment run cycles, but it would have been better to atleast specify the platform.  

```
 bwa mem -t 8 -R "@RG\tID:rg1\tSM:foo" hg38_chr7 forward_paired.fastq.gz reverse_paired.fastq.gz > 392_aln.sam

```    
We generate a SAM file, but it is more useful and less space consuming to work with binary formats and so we convert it. and so we use the fixmate function from `samtools` to convert and adjust the reads.    
`samtools fixmate -O bam 392_aln.sam 392_aln.bam`


> The subsequent data processing steps will be done using Picard and GATK according to the best practices outlined in the documentation.

### Picard Workflow
**Sam File Validation**   
To make sure we have no issues in our alignment map files, we pass them through a validation step using `picard`.
```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar ValidateSamFile -INPUT 392_aln.bam -MODE SUMMARY 

```
We have a single warning regarding the ReadGroup information as we had not specified the sequencing platform (Illumina) information in the alignment step. This could be easily corrected later on.

**Sorting the BAM File**  
It is important to sort the bam file. First the read alignments are sorted by the reference sequence name (RNAME) field using the reference sequence dictionary. These are then sorted according to mapping position

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar SortSam -INPUT 392_aln.bam -OUTPUT 392_aligned_sorted.bam -SORT_ORDER coordinate

```

**Removing Duplicate Reads**  
Removing PCR duplicates and optical contaminants.
```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar MarkDuplicates -INPUT 392_aligned_sorted.bam -OUTPUT 392_dup_marked.bam -METRICS_FILE 392_metrics.metrics

```
The input here is simply our sorted bam file and we get a deduplicated bam file.

**Editing the Read Group Name to add Platform**  

Since we had an error previously regarding the read group, it is important to edit it so as to avoid errors later on with data recalibration. I had not noticed that I did not do this before Marking the duplicates and this had caused me some trouble when doing the BQSR. 

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar AddOrReplaceReadGroups  -I 392_aligned_sorted.bam  -O 392_aligned_sorted_RGcorr.bam  RGID=8  RGLB=lib1  RGPL=ILLUMINA  RGPU=unit1  RGSM=20 

```

**Marking Duplicates on Edited File**  

```
java -jar /home/mohamed.mehdi/picard/build/libs/picard.jar MarkDuplicates -INPUT 392_aligned_sorted_RGcorr.bam -OUTPUT 392_duplicate_marked.bam -METRICS_FILE markDuplicate_metrics.metrics

```
### Data Recalibration

**Base Quality Score Recalibration**  
First we build the BQSR model based on our data. This would give us a table that we could use to do the acual recalibration with `ApplyBQSR`.
```
gatk BaseRecalibrator -I 392_duplicate_marked.bam -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa --known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf -O recalibrated_data.table 
 
```
**Applying BQSR**  
Now that we have built the model, we simply do the BQSR.
```
gatk ApplyBQSR -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -I 392_duplicate_marked.bam --bqsr-recal-file recalibrated_data.table -O 392_recalibrated.bam 
```
**Second Pass Recalibration**  
To compare the data before and after doing the recalibration let's generate another table with a second bass of BQSR on our recalibrated data.
```
gatk BaseRecalibrator -I 392_recalibrated.bam -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa --known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf -O secondPass.table
```
**Covariate Analysis: Before and After BQSR**  
Using the two tables we have we can do a co-variate analysis to see how our base substitution data varies.
```
 gatk AnalyzeCovariates -before recalibrated_data.table -after secondPass.table -plots covAnalysis.pdf

```   
<img src="/WES_Workflow/images/CoVariates.png"  width="1024" height="520">

### Variant Calling  
Our processing steps are done and we have analysis ready files. We can do variant calling using GATK through two steps:
1. `HaplotypeCaller`
```
gatk --java-options "-Xmx4g" HaplotypeCaller -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -I 392_recalibrated.bam -O 392_GvarCall.g.vcf.gz -ERC GVCF

```
2. `GenotypeGVCFs`

```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_GvarCall.g.vcf.gz -O 392_varCall.vcf.gz

```
 
 <img src="/WES_Workflow/images/Types_of_Variants.pie.png"  width="520" height="520">
 

**Filtering the Variants**  
Filtering the variants based on certain criteria. The options for the command are the ones suggested in the standard pipeline.
```
gatk --java-options "-Xms5g -Xmx15g" VariantFiltration -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_varCall.vcf.gz -O 392_varCall_filtered.vcf.gz --filter-name "lowGQ"     --filter-expression "GQ < 20.0"     --filter-name "lowMQ"     --filter-expression "MQ < 40.0"     --filter-name 'lowQD'  --filter-expression "QD < 2.0"     --filter-name "lowMQRankSum"     --filter-expression "MQRankSum < -12.5" 

```
**Selecting the Variants**  
We filter out INDELs and keep the SNPs only.
```
gatk --java-options "-Xms5g -Xmx15g" SelectVariants -R /home/mohamed.mehdi/WholeExomeProject/chrom7/hg38_chr7.fa -V 392_varCall_filtered.vcf.gz -O 392_varCall_filteredExcluded_nonVar.vcf.gz --exclude-filtered true --exclude-non-variants true 

```

<img src="/WES_Workflow/images/Types_of_Variants_Excluded.pie.png"  width="520" height="520">

**Variant Annotation**  
Variant annotation was not finished due to a problem with loading the hg38 snp database.


### General Stats  
| ID | Stats |
| ------ | ------ |
| Number of mapped reads | 8371923 (13.4%)  |  
| Avg Mapping Quality| 25.5262 |
| Number of Supplementary Reads | 670477 |
| Number of Secondary Reads | 0 |
| Number of Reads Without Pair Complement | 58186023 |
| Number of Duplicates | 1296614 |
| Number of Reads with no Indels | 61300855 |
| Number of INDELs | 7203 |
| Number of SNPs | 122768 |


#### Tools Used:
`fastqc v 0.11.9` 
`bwa v 0.7.17`
`samtools v 1.4`
`bcftools v 1.4`
`GATK v 4.1.9.0`
`snpEff v 4.3`



