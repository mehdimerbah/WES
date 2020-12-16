#!/usr/bin/env bash


####################### FASTQ ############################

#This Script can be used to filter out bad tiles manually, below is an example only, edit the range + ID 
grep -P -A 4 '@HWI-D00119:50:H7AP8ADXX:1:(?!1[1-2][0-9][0-9])' sample.fastq

#####################    SAM/BAM    ######################

samtools flagstat 392_aln.bam 

#################### Mean MAPQ #######################
#Average MAPQ for mapped reads: Filter out mapped reads and add up their quality then divide by the count.

samtools view -F 0x04 392_aln.bam | awk '{s+=$5} END {print "Avg Quality: ", s/NR}'
#Avg Quality:  25.5262


################## No Indel Reads ####################

#Reads with no insertions or deletions:

samtools view 392_aligned_sorted_RGcorr.bam | cut -f 6 | awk  '{if($0 !~ /[ID]/) sum+=1} END {print "No Ins or Del:", sum }' 

#No Ins or Del: 61300855


######################### VCF/BCF ####################





