##
##    ___               _ _             ___ _  _ ___   _  ___ _  _ ___  ___ _    
##   / __|___ _ _ _ __ | (_)_ _  ___   / __| \| | _ \_| ||_ _| \| |   \| __| |   
##  | (_ / -_) '_| '  \| | | ' \/ -_)  \__ \ .` |  _/_   _| || .` | |) | _|| |__ 
##   \___\___|_| |_|_|_|_|_|_||_\___|  |___/_|\_|_|   |_||___|_|\_|___/|___|____|
## 
##            ___ _   ___ _____ ___          _                _____   _____ ___ 
##   _ __ ___| __/_\ / __|_   _/ _ \        | |_ ___         / __\ \ / / __| __|
##  | '_ / -_| _/ _ \\__ \ | || (_) |       |  _/ _ \       | (_ |\ V | (__| _| 
##  | .__\___|_/_/ \_|___/ |_| \__\_\  ___   \__\___/  ___   \___| \_/ \___|_|  
##  |_|                               |___|           |___|                     
## 
## 
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF generation) according to the 
## GATK Best Practices (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) for 
## germline SNP and Indel discovery in human whole-genome sequencing (WGS) data. 
## 
## Runtime parameters are optimized for Rivanna.
##
## Version: 0.1
## Date: 05/14/2018
## Author: Hardik I. Parikh
##
##  
## Inputs:
##	- Human Reference Genome - hg38
##	-	Resouce files - dbSNP, hapmap, 1000G
##
##	- WholeGenome/WholeExome paired-end sequencing data in FASTQ format
##
##			---------
##			IMPORTANT:  
##			---------
##			This workflow assumes that reads were sequenced from ONE DNA library,
##			sequenced on ONE lane of Illumina
##
##			MarkDuplicates and BQSR, important pre-processing steps, are heavily 
##			dependent on this information. 
##
##			If your experimental design was different, i.e. multiple libraries were
##			were prepared for the sample, and sequenced across multiple lanes of 
##			Illumina HiSeq platform, please consult with us before running the ananlysis!
##
##
## Outputs:
##	- Single Sample variant calls in GVCF format. 
##
##
## Rivanna Modules:
##	- module load gatk/4.0.0.0
##	- module load bwa/0.7.17
##	- module load samtools/1.7
##	- module load wdltool/0.14
##	- module load cromwell/30.1
##
##	
## LICENSING: 
##
## Copyright (C) 2018 Hardik I Parikh
## 
## This script is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This script is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this script. If not, see <http://www.gnu.org/licenses/>.
##
## For any bugs or problems found, please contact us at:
## Hardik I. Parikh <hiparikh@virginia.edu>

#######################
# WORKFLOW DEFINITION #
#######################
workflow peFastq_to_GVCF {

	String sample_name
	File r1fastq
	File r2fastq
	Int? read_length
	
	File ref_fasta
	File ref_fasta_2bit
	File ref_fasta_index
	File ref_dict
	File ref_alt
	File ref_bwt
	File ref_sa
	File ref_amb
	File ref_ann
	File ref_pac	
	
	File dbSNP_vcf
	File dbSNP_vcf_index
	Array[File] known_indels_sites_VCFs
	Array[File] known_indels_sites_indices
	
	# map to reference
	## create the READ-GROUP TAG
	String rgtag="'@RG\\tID:" + sample_name + "_" + sample_name + "\\tLB:" + sample_name + "\\tSM:" + sample_name + "\\tPL:ILLUMINA'" 
	call bwamem {
		input:
			sample_name = sample_name,
			r1fastq = r1fastq,
			r2fastq = r2fastq,
			ref_fasta = ref_fasta,
			ref_bwt = ref_bwt,
			ref_sa = ref_sa,
			ref_amb = ref_amb,
			ref_ann = ref_ann,
			ref_pac = ref_pac,
			rgtag = rgtag
	}

	call sortbam {
		input:
			sample_name = sample_name, 
			inbam = bwamem.out
	}

	# mark duplicates
	call markdups {
		input:
		sample_name = sample_name,
		inbam = sortbam.outbam,
		inbamidx = sortbam.outbamidx
	}

	# bqsr
	call bqsr {
		input:
			sample_name = sample_name,
			inbam = markdups.outbam,
			ref_fasta_2bit = ref_fasta_2bit,
			dbSNP_vcf = dbSNP_vcf, 
			dbSNP_vcf_index = dbSNP_vcf_index,
			known_indels_sites_VCFs = known_indels_sites_VCFs,
			known_indels_sites_indices = known_indels_sites_indices	
	}
		
	call applybqsr {
		input:
			sample_name = sample_name,
			inbam = markdups.outbam,
			ref_fasta_2bit = ref_fasta_2bit,
			recaltable = bqsr.recaltable
	}

	# haplotype caller
	call haplotypecaller {
		input:
			sample_name = sample_name,
			inbam = applybqsr.outbam,
			ref_fasta_2bit = ref_fasta_2bit
	}

	# outputs
	output {
		File sortedbam = sortbam.outbam
		File sortedbamidx = sortbam.outbamidx

		File mkdupbam = markdups.outbam
		File mkdupmetrics = markdups.outtxt
		File mkduplog = markdups.outlog

		File bqsrtab = bqsr.recaltable
		File bqsrlog = bqsr.outlog
		File applybqsrbam = applybqsr.outbam
		File applybqsrlog = applybqsr.outlog

		File hcgvcf = haplotypecaller.outgvcf
		File hclog = haplotypecaller.outlog
	}

}


###################
# TASK DEFINITION #
###################

# This task will align paired-end reads to ref genome using bwa-mem
# and will write the output in BAM format
# It will also add the RG tag, required for Picard to MarkDuplicates
task bwamem {
	String sample_name
	File r1fastq
	File r2fastq
	File ref_fasta
	File ref_bwt
	File ref_sa
	File ref_amb
	File ref_ann
	File ref_pac
	String rgtag

	command <<<
		bwa mem -M -R ${rgtag} -t 16 ${ref_fasta} ${r1fastq} ${r2fastq} | \
			samtools view -bS - > ${sample_name}.hg38-bwamem.bam
	>>>
	
	runtime {
		cpu: 16
		requested_memory_mb: 16000
	}

	output {
		File out = "${sample_name}.hg38-bwamem.bam"
	}
	
}

# This task will sort the input BAM file by coordinates
task sortbam {
	File inbam
	String sample_name
	command <<<
		gatk --java-options "-Xms20g" SortSam \
			--INPUT ${inbam} \
			--OUTPUT ${sample_name}.hg38-bwamem.sorted.bam \
			--SORT_ORDER coordinate \
			--CREATE_INDEX
	>>>
	runtime{
		requested_memory_mb: 24000
	}
	output {
		File outbam = "${sample_name}.hg38-bwamem.sorted.bam"
		File outbamidx = "${sample_name}.hg38-bwamem.sorted.bai"
	}
}

# This task will mark duplicates
task markdups {
	File inbam 
	File inbamidx 
	String sample_name
	command <<<
		gatk --java-options "-Xms80g" MarkDuplicatesSpark \
			--spark-master local[16] \
			--input ${inbam} \
			--METRICS_FILE ${sample_name}.mkdup-spark.metrics.txt \
			--output ${sample_name}.hg38-bwamem.sorted.mkdup-spark.bam \
			--read-validation-stringency SILENT &> \
			gatk-MarkDuplicatesSpark.log
	>>>
	runtime {
		cpu: 16
		requested_memory_mb: 80000
	}
	output {
		File outbam = "${sample_name}.hg38-bwamem.sorted.mkdup-spark.bam"
		File outtxt = "${sample_name}.mkdup-spark.metrics.txt"
		File outlog = "gatk-MarkDuplicatesSpark.log"
	}
}


#BaseRecalibrator
task bqsr {
	String sample_name
	File inbam
	File ref_fasta_2bit
	File dbSNP_vcf
	File dbSNP_vcf_index
	Array[File] known_indels_sites_VCFs
	Array[File] known_indels_sites_indices
	
	command <<<
		gatk --java-options "-Xms80g" BaseRecalibratorSpark \
			--spark-master local[16] \
			--input ${inbam} \
			--reference ${ref_fasta_2bit} \
			--use-original-qualities \
			--known-sites ${dbSNP_vcf} \
			--known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
			--output ${sample_name}.hg38-bwamem.sorted.mkdup-spark.recal-spark.table &> \
			gatk-BaseRecalibratorSpark.log
	>>>
	runtime {
		cpu: 16
		requested_memory_mb: 80000
	}
	output {
		File recaltable = "${sample_name}.hg38-bwamem.sorted.mkdup-spark.recal-spark.table"
		File outlog = "gatk-BaseRecalibratorSpark.log"
	}
}


#Apply BQSR
task applybqsr {
	String sample_name
	File inbam
	File ref_fasta_2bit
	File recaltable

	command <<<
		gatk --java-options "-Xms80g" ApplyBQSRSpark \
			--spark-master local[16] \
			--input ${inbam} \
			--reference ${ref_fasta_2bit} \
			--bqsr-recal-file ${recaltable} \
			--use-original-qualities \
			--static-quantized-quals 10 \
			--static-quantized-quals 20 \
			--static-quantized-quals 30 \
			--output ${sample_name}.hg38-bwamem.sorted.mkdup-spark.bqsr-spark.bam &> \
			gatk-ApplyBQSRSpark.log
	>>>
	runtime {
		cpu: 16
		requested_memory_mb: 80000
	}
	output {
		File outbam = "${sample_name}.hg38-bwamem.sorted.mkdup-spark.bqsr-spark.bam"
		File outlog = "gatk-ApplyBQSRSpark.log"
	}	
}


#HaplotypeCaller
task haplotypecaller {
	String sample_name
	File inbam
	File ref_fasta_2bit

	command <<<
		gatk --java-options "-Xmx80g" HaplotypeCallerSpark \
			--spark-master local[16] \
			--input ${inbam} \
			--reference ${ref_fasta_2bit} \
			-ERC GVCF \
			-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
			--max-alternate-alleles 3 \
			--read-filter OverclippedReadFilter \
			--output ${sample_name}.hg38-bwamem.sorted.mkdup-spark.bqsr-spark.gvcf &> \
			gatk-HaplotypeCallerSpark.log

	>>>
	runtime {
		cpu: 16
		requested_memory_mb: 80000
	}
	output {
		File outgvcf = "${sample_name}.hg38-bwamem.sorted.mkdup-spark.bqsr-spark.gvcf"
		File outlog = "gatk-HaplotypeCallerSpark.log"
	}
}
