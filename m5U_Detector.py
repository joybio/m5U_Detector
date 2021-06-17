configfile:  "m5U_Detector.yaml"
#__date__ = "2021-5-27"
#__author__ = "Junbo Yang"
#__email__ = "yang_junbo_hi@126.com"
#__license__ = "PKU.jia.group"


#Snakefile
#Snakefile中的每一个rule其实都可以看作是一个简单的shell脚本，通过Snakefile将多个rule组织在一起并按照
#我们定义的顺序来执行。

#from snakemake.utils import makedirs
#from snakemake.utils import listfiles

#import numpy as np
#import os
sample = config["samples"]
index = config["hisat2_idx"]
strandness = config["strandness"]
rule all:
	## LOCAL ##
#	'''
#	Defines the target rule by specifying all final output files.
#	Additionally, the cluster_logs dir is deleted if
#	snakemake was run locally.
#	'''
	input: 
		#expand(config["results_dir"] + "/{sample}.R1.test.uniq.sorted.3pend",sample=sample)
		config["results_dir"] + "/passthrough.xls"		

#	output:
#		expand(config["results_dir"] + "/{sample}.test.uniq.3pend.bed",sample=sample)

#rule cluster:
#    input:
#	script = 'python/dbscan.py',
#	path   = 'data_files/umap/{sample}_umap.csv'
#    output:
#	path = 'output/{sample}'
#    shell:
#	"python {input.script} -data {input.path} -eps '0.3' -min_samples '10' "


#-------------------------------------------------------------------------------
# gunzip rule1: Dependency packages - None
#-------------------------------------------------------------------------------
rule gunzip:
	input: 
		config["input_dir"] + "/{sample}.combined_R1.fastq.gz"
	output:
		config["results_dir"] + "/{sample}.combined_R1.fastq"
	message:
		"starting unzip ..."
	shell:
		'''
		gunzip -c {input} > {output}
		'''
#-------------------------------------------------------------------------------
# UMI_rmdup rule2: Dependency packages - seqkit
# remove PCR duplication by UMI
#-------------------------------------------------------------------------------
rule UMI_rmdup:
	input:
		config["results_dir"] + "/{sample}.combined_R1.fastq"
	output:
		config["results_dir"] + "/{sample}.combined_R1.rmdup.fq",
		config["results_dir"] + "/{sample}.combined_R1.duplicated.fq.gz",
		config["results_dir"] + "/{sample}.combined_R1.duplicated.txt"
	message:
		"starting remove PCR duplication ..."
	shell:
		'''
		seqkit rmdup -s {input[0]} -o {output[0]} -d {output[1]} -D {output[2]} -j 10
		'''
#-------------------------------------------------------------------------------
# cutadapt rule3: Dependency packages - cutadapt
#-------------------------------------------------------------------------------
rule cutadapt:
	input:  
		config["results_dir"] + "/{sample}.combined_R1.rmdup.fq"
	output: 
		config["results_dir"] + "/{sample}.trimmed.R1.fq.gz"
	message: "starting cutadaptor ..."
	params:
		error_rate=config['error_rate'],
		minLen=config['min_length'],
		overlap=config['overlap']
	shell:
		'''
		cutadapt -j 10 -a NNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			-e {params.error_rate} -O {params.overlap} -m {params.minLen} \
			-o {output[0]} {input[0]}
		'''
#-------------------------------------------------------------------------------
# map rule4: Dependency packages - hisat2
#-------------------------------------------------------------------------------

rule hisat2_mapping:
	input:
		config["results_dir"] + "/{sample}.trimmed.R1.fq.gz"
##因为这里的变量名是局部变量，仅在每个rule定义的子部分起作用，因此不同子部分可以用相同变量名，但内容取决于'='后面的部分
	output:
		config["results_dir"] + "/{sample}.R1.sam"
 #temp()用于标记临时文件，最后会自动删除，以免sam文件过大挤兑储存空间
	log:
		config["log_dir"] + "/{sample}.hisat2.log"
	message:
		"start hisat2 mapping..."
	shell:
		'''
		hisat2 -p 10 -k 1 --rna-strandness {strandness} --pen-noncansplice 1000000 \
		-x {index} \
		-U {input} \
		-S {output} \
		 > {log} 2>&1
		'''
#-------------------------------------------------------------------------------
# sam2bam rule5: Dependency packages - samtools
#-------------------------------------------------------------------------------

rule sam2bam:
	input:
		config["results_dir"] + "/{sample}.R1.sam"
	output:
		config["results_dir"] + "/{sample}.R1.bam"
	shell:
		"""
		samtools view -@ 10 -bS -o {output} {input}
		"""
#-------------------------------------------------------------------------------
# resort rule6: Dependency packages - None
#-------------------------------------------------------------------------------
rule bam_resort:
	input:
		config["results_dir"] + "/{sample}.R1.bam"
	output:
		config["results_dir"] + "/{sample}.R1.resort.bam"
	shell:
		"""
		samtools sort -@ 10 -O BAM -o {output} {input}
                """
#-------------------------------------------------------------------------------
# depth rule7: Dependency packages - None
#-------------------------------------------------------------------------------
rule bam2depth:
	input:
		bam = expand(config["results_dir"] + "/{sample}.R1.resort.bam",sample = sample)
	output:
		config["results_dir"] + "/merge.depth"
	shell:
		"""
		ls {input.bam} > depth.bam.txt | samtools depth -aa -b hg38.tRNA.strand.bed -f depth.bam.txt > {output}
		"""
#-------------------------------------------------------------------------------
# depth2 rule8: Dependency packages - None
#-------------------------------------------------------------------------------

rule depth2bed:
	input:
		config["results_dir"] + "/merge.depth"
	output:
		config["results_dir"] + "/merge.depth.bed"
	shell:
		"""
		python ./scripts/depth2bed.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# extract_tRNA_info rule9: Dependency packages - intersectBed
#-------------------------------------------------------------------------------
rule extract_tRNA_info:
	input:
		config["tRNA_database"] + "/homo-mature-tRNA.format.fa"
	output:
		config["results_dir"] + "/hg38.tRNA.strand.bed"
	shell:
		"""
		python ./scripts/extract.bed.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# depth_tRNA rule10: Dependency packages - intersectBed
#-------------------------------------------------------------------------------
rule depth_tRNA:
	input:
		config["results_dir"] + "/merge.depth.bed",
		config["results_dir"] + "/hg38.tRNA.strand.bed"
	output:
		config["results_dir"] + "/tRNA.merge.depth.bed"
	shell:
		"""
		intersectBed -a {input[0]} -b {input[1]} -wo > {output[0]}
		"""
#-------------------------------------------------------------------------------
# prepare_tRNA_sequence rule11: Dependency packages - None
#-------------------------------------------------------------------------------

rule prepare_tRNA_seq:
	input:
		bed = config["results_dir"] + "/hg38.tRNA.strand.bed",
		fa = config["fa_file"]
	output:
		config["results_dir"] + "/hg38.tRNA.seq.fa"
	shell:
		"""
		bedtools getfasta -fi {input.fa} -bed {input.bed} -s -fo {output}
		"""
#-------------------------------------------------------------------------------
# prepare_single_nucleotide_bed rule12: Dependency packages - None
#-------------------------------------------------------------------------------

rule extract_single_nucleotide:
	input:
		config["results_dir"] + "/tRNA.merge.depth.bed"
	output:
		config["results_dir"] + "/tRNA.slice.merge.depth.bed"
	shell:
		"""
		python ./extract_single_nucleotide.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# getfasta rule13: Dependency packages - None
#-------------------------------------------------------------------------------
rule getfasta:
	input:
		fa = config["fa_file"],
		bed = config["results_dir"] + "/tRNA.slice.merge.depth.bed"
	output:
		config["results_dir"] + "/hg38.tRNA.slice.strand.fa"
	shell:
		"""
		bedtools getfasta -fi {input.fa} -bed {input.bed} -s -fo {output}
		"""
#-------------------------------------------------------------------------------
# merge_info rule14: Dependency packages - None
#-------------------------------------------------------------------------------

rule merge_info:
	input:
		tRNA_seq = config["results_dir"] + "/hg38.tRNA.seq.fa",
		tRNA_slice = config["results_dir"] + "/hg38.tRNA.slice.strand.fa",
		depth = config["results_dir"] + "/tRNA.merge.depth.bed"
	output:
		config["results_dir"] + "/tRNA.merge.seq.slice.depth.bed"
	shell:
		"""
		python ./scripts/merge.py -d {input.depth} -s {input.slice} -f {input.tRNA_seq} -o {output}
		"""
#-------------------------------------------------------------------------------
# orientation_adjust rule15: Dependency packages - None
#-------------------------------------------------------------------------------
rule orientation_adjust:
	input:
		config["results_dir"] + "/tRNA.merge.seq.slice.depth.bed"
	output:
		config["results_dir"] + "/new_tRNA.merge.seq.slice.depth.bed"

	shell:
		"""
		python ./scripts/strand_adjust.py -i {input} -o {output}
		"""


#-------------------------------------------------------------------------------
# passthrough rule16: Dependency packages - None
#-------------------------------------------------------------------------------

rule passthrough:
	input:
		config["results_dir"] + "/new_tRNA.merge.seq.slice.depth.bed"
	output:
		config["results_dir"] + "/passthrough.xls"
	shell:
		"""
		python ./scripts/pass_through.py -i input -o output
		"""




