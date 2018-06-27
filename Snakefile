configfile: "Brain_v1.8.config.json"

import os

# Tools
# Variables should point to full paths to tools. Currently assumes all are available in one's PATH (e.g., installed with Bioconda)
fastqc_path = "fastqc"
multiqc_path = "multiqc"
trimmomatic_path = "trimmomatic"
hisat2_build_path = "hisat2-build"
hisat2_path = "hisat2"
samtools_path = "samtools"
FEATURECOUNTS = "featureCounts"
xyalign_env_path = "xyalign_env"
xyalign_path = "xyalign"
stringtie_path = "stringtie"

# Directory variables
fastq_directory = "/mnt/storage/public/dbgap-8834/brain_amygdala/"
fastqc_directory = "fastqc_results/"

# # Reference genome files: XX with Y chromosome masked, XY with Y PAR regions masked
# XX_HISAT2_INDEX_WITH_VIRAL_REF = config["XX_GRCh38_ref_with_viral_genomes_HISAT2_index"]
# XY_HISAT2_INDEX_WITH_VIRAL_REF = config["XY_withoutYpar_GRCh38_ref_with_viral_genomes_HISAT2_index"]
# GTF = config["GRCh38_gtf_path"]
# GTF_PROTEINCODING = config["GRCh38_proteincoding_gtf_path"]
# GTF_AUTOSOMAL = config["GRCh38_gtf_withoutXY_path"]

# # Directories
# FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Brain_Amygdala/RNAseq/trimmed_fastqs/" # path to directory with fastq files
# SAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Brain_Amygdala/RNAseq/GRCh38_SAM/" # path to directory for SAM alignment files
# BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Brain_Amygdala/RNAseq/GRCh38_BAM/" # path to directory for BAM alignment files
# SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Brain_Amygdala/RNAseq/GRCh38_sorted_BAM/" # path to directory for sorted BAM alignment files
# FEATURECOUNTS_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Brain_Amygdala/RNAseq/GRCh38_featurecounts/" # path to directory for FeatureCounts gene count files


# Samples
XX_SAMPLES = config["Brain_Amygdala_Female_RNA"]
XY_SAMPLES_wo = config["Brain_Amygdala_Male_RNA_wo_ZAB4"]
XY_SAMPLES = config["Brain_Amygdala_Male_RNA"]
SAMPLES = XX_SAMPLES + XY_SAMPLES_wo


rule all:
	input:
		expand(
			"fastqc_results/{sample}_fixed_{num}_fastqc.html",
			sample=SAMPLES, num=[1,2]),
		"multiqc_results/multiqc_report.html",
		expand(
			"logfiles/{sample}_trimmomatic.log",
			sample=SAMPLES),
		expand(
			"fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_{num}_fastqc.html", 
			sample=SAMPLES, num=[1,2]),
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"xyalign/reference/{assembly}_{chr}.fa",
			assembly=["hg38"], chr=["xx", "xy"]),
		expand(
			"hisat2_index/{assembly}_{chr}.6.ht2",
			assembly=["hg38"], chr=["xx", "xy"]),
		expand(
			"bams/{sample}_{assembly}.sorted.bam",
			assembly=["hg38"], sample=SAMPLES),
		expand(
			"stringtie_results_postmerge/{sample}_{assembly}.assembled_transcripts.gtf",
			assembly=["hg38"], sample=SAMPLES)

rule fastqc_analysis:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}/{sample}_fixed_1.fastq"),
		fq2 = os.path.join(fastq_directory, "{sample}/{sample}_fixed_2.fastq")
	output:
		ofq1= "fastqc_results/{sample}_fixed_1_fastqc.html",
		ofq2= "fastqc_results/{sample}_fixed_2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input.fq1} {input.fq2}"
		
rule multiqc:
	input:
		expand("fastqc_results/{sample}_fixed_{num}_fastqc.html", sample=SAMPLES, num=[1,2])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"{params.multiqc} --interactive fastqc_results -o multiqc_results"	

rule trimmomatic:
 	input:
 		fq1 = os.path.join(fastq_directory, "{sample}/{sample}_fixed_1.fastq"),
		fq2 = os.path.join(fastq_directory, "{sample}/{sample}_fixed_2.fastq"),
		ADAPTER_FASTA = "misc/adapter_sequences.fa"
	output:
		paired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq",
		paired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq",
		unpaired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_unpaired_1.fastq",
		unpaired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_unpaired_2.fastq",
		logfile = "logfiles/{sample}_trimmomatic.log"
	params:
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 3,
		trailing = 3,
		winsize = 4,
		winqual = 30,
		minlen = 50
	shell:
		"trimmomatic PE -threads {params.threads} -trimlog {output.logfile} {input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} {output.paired_2} {output.unpaired_2} ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}"

rule fastqc_analysis_trimmed:
	input:
		expand("trimmed_fastqs/{{sample}}_trimmomatic_trimmed_paired_{num}.fastq",num=[1,2])
	output:
		ofq1= "fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_1_fastqc.html",
		ofq2= "fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input}"

rule multiqc_trimmed_paired:
	input:
		expand("fastqc_trimmed_results/{sample}_trimmomatic_trimmed_paired_{num}_fastqc.html", sample=SAMPLES, num=[1,2])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"{params.multiqc} --interactive fastqc_trimmed_results -o multiqc_trimmed_results"

rule xyalign_prepare_ref:
	input:
		lambda wildcards: config["ref_genome"][wildcards.assembly]
	output:
		xx= "xyalign/reference/{assembly}_xx.fa",
		xy= "xyalign/reference/{assembly}_xy.fa"
	params:
		xyalign_env = xyalign_env_path, 
		xyalign = xyalign_path 
	shell:
		"source activate {params.xyalign_env} && "
		"xyalign --PREPARE_REFERENCE --ref {input} --output_dir xyalign "
		"--reference_mask hg38_PAR_mask.bed --x_chromosome chrX --y_chromosome chrY "
		"--xx_ref_out {wildcards.assembly}_xx.fa --xy_ref_out {wildcards.assembly}_xy.fa"

rule hisat2_index:
	input:
		xx= "xyalign/reference/{assembly}_xx.fa",
		xy= "xyalign/reference/{assembly}_xy.fa"
	output:
		xx= "hisat2_index/{assembly}_xx.6.ht2",
		xy= "hisat2_index/{assembly}_xy.6.ht2"
	params:
		hisat2=hisat2_build_path,
		base_name1="hisat2_index/{assembly}_xx",
		base_name2="hisat2_index/{assembly}_xy"
	run:
		shell("{params.hisat2} {input.xx} {params.base_name1}")
		shell("{params.hisat2} {input.xy} {params.base_name2}")

rule hisat2_align_reads:
	input:
		paired_1 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_1.fastq",
		paired_2 = "trimmed_fastqs/{sample}_trimmomatic_trimmed_paired_2.fastq",
		xx= "hisat2_index/{assembly}_xx.6.ht2",
		xy= "hisat2_index/{assembly}_xy.6.ht2"
	output:
		"bams/{sample}_{assembly}.sorted.bam"
	params:
		hisat2=hisat2_path,
		threads= 4,
		base_name_xx="hisat2_index/{assembly}_xx",
		base_name_xy="hisat2_index/{assembly}_xy",
		samtools= samtools_path
	threads: 4
	run:
		if wildcards.sample in XY_SAMPLES_wo:
			shell(
				"{params.hisat2} -p {params.threads} -x {params.base_name_xy} -1 {input.paired_1} -2 {input.paired_2} | "
				"{params.samtools} view -b - | {params.samtools} sort -O bam -o {output} -")
		else:
			shell(
				"{params.hisat2} -p {params.threads} -x {params.base_name_xx} -1 {input.paired_1} -2 {input.paired_2} | "
				"{params.samtools} view -b - | {params.samtools} sort -O bam -o {output} -")

rule stringtie:
	input:
		aligned_reads="bams/{sample}_{assembly}.sorted.bam"
	output:
		assembled_transcripts="stringtie_results/{sample}_{assembly}.assembled_transcripts.gtf",
		gene_abundances="stringtie_results/{sample}_{assembly}.gene_abundances.txt"
		# fully_covered_transcripts="stringtie_results/{sample}_{assembly}.fully_covered_transcripts.gtf"
	params:
		stringtie= stringtie_path
	shell:
		"{params.stringtie} {input.aligned_reads} -o {output.assembled_transcripts} "
		"-A {output.gene_abundances}"
		# -C {output.fully_covered_transcripts} 

rule stringtie_merge:
	input:
		assembled_transcripts=lambda wildcards: expand(
			"stringtie_results/{sample}_{assembly}.assembled_transcripts.gtf", 
			sample=SAMPLES, assembly=[wildcards.assembly])
	output:
		"stringtie_results/{assembly}.stringtie_merge.gtf"
	params:
		stringtie= stringtie_path
	shell:
		"{params.stringtie} --merge {input.assembled_transcripts} -o {output} "

rule stringtie_2:
	input:
		G="stringtie_results/{assembly}.stringtie_merge.gtf",
		aligned_reads="bams/{sample}_{assembly}.sorted.bam"
	output:
		assembled_transcripts="stringtie_results_postmerge/{sample}_{assembly}.assembled_transcripts.gtf",
		gene_abundances="stringtie_results_postmerge/{sample}_{assembly}.gene_abundances.txt",
		fully_covered_transcripts="stringtie_results_postmerge/{sample}_{assembly}.fully_covered_transcripts.gtf"
	params:
		stringtie= stringtie_path
	shell:
		"{params.stringtie} {input.aligned_reads} -G {input.G} "
		"-B -e "
		"-o {output.assembled_transcripts} "
		"-A {output.gene_abundances} "
		"-C {output.fully_covered_transcripts} "
