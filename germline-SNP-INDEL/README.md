
## Introduction:  
This WDL pipeline implements data pre-processing and initial variant calling (GVCF generation) according to the 
GATK Best Practices (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) for 
germline SNP and Indel discovery in human whole-genome sequencing (WGS) data. 

Runtime parameters are optimized for Rivanna.

Version: 0.1
Date: 05/14/2018
Author: Hardik I. Parikh

## Inputs:
	- Human Reference Genome - hg38
	-	Resouce files - dbSNP, hapmap, 1000G

	- WholeGenome/WholeExome paired-end sequencing data in FASTQ format

			---------
			IMPORTANT:  
			---------
			This workflow assumes that reads were sequenced from ONE DNA library,
			sequenced on ONE lane of Illumina

			MarkDuplicates and BQSR, important pre-processing steps, are heavily 
			dependent on this information. 

			If your experimental design was different, i.e. multiple libraries were
			were prepared for the sample, and sequenced across multiple lanes of 
			Illumina HiSeq platform, please consult with us before running the ananlysis!


## Outputs:
	- Single Sample variant calls in GVCF format. 


## Rivanna Modules:
	- module load gatk/4.0.0.0
	- module load bwa/0.7.17
	- module load samtools/1.7
	- module load wdltool/0.14
	- module load cromwell/30.1
