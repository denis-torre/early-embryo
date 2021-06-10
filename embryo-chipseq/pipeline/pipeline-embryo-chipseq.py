#################################################################
#################################################################
############### Embryo ChIP-Seq ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Guccione Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import ruffus.cmdline as cmdline
import sys
import os
import json
import glob
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/embryo-chipseq.R'
py_source = 'pipeline/scripts/EmbryoChipseq.py'
P = 'acc_GuccioneLab'
q = 'express'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import EmbryoChipseq as P

# 3.2 R
# r.source(r_source)

#############################################
########## 2. General Setup
#############################################
##### 1. Pipeline running #####
# Pipeline args
options = cmdline.get_argparse().parse_args()

##### 2. Genome indices #####
# Open JSON
with open('/sc/arion/projects/GuccioneLab/genome-indices/genome-indices.json') as openfile:
	genome_indices = json.load(openfile)

##### 3. Variables #####
# All FASTQ
all_illumina_fastq = glob.glob('arion/chipseq/s01-fastq.dir/*/*/*/*.f*q.gz')

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():
	fastq_files = glob.glob('arion/datasets/xia/rawdata/*/*.fastq.gz')
	for fastq_file in fastq_files:
		sample_name = fastq_file.split('/')[-2]
		fastq_name = os.path.basename(fastq_file)
		read_mate = '_'+fastq_name.split('_')[-1][0] if '_' in fastq_name else ''
		infile = os.path.join(os.getcwd(), fastq_file)
		outfile = 'arion/chipseq/s01-fastq.dir/human/raw/{sample_name}/{sample_name}{read_mate}.fastq.gz'.format(**locals())
		yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/chipseq/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

# @follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	if len(infiles) == 1:
		cmd_str = '''trim_galore --nextera --cores 6 --output_dir {outdir} {infiles[0]}'''.format(**locals())
	elif len(infiles) == 2:
		cmd_str = '''trim_galore --nextera --paired --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	print(cmd_str)
	# run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='02:00', GB=6, n=6, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

# @follows(linkFASTQ, trimIlluminaAdapters)

@transform(all_illumina_fastq,
		   regex(r'(.*)/s01-fastq.dir/(.*).f.*q.gz'),
		   r'\1/s02-fastqc.dir/\2_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=True)

# ls arion/chipseq/s02-fastqc.dir/human/raw | wc -l

#############################################
########## 2. MultiQC
#############################################

# @follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/chipseq/s02-fastqc.dir/human/raw', 'arion/chipseq/multiqc/human_fastqc/multiqc_report.html'],
		# ['arion/chipseq/s02-fastqc.dir/human/trimmed', 'arion/chipseq/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		# ['arion/chipseq/s03-alignment.dir/human/bowtie2', 'arion/chipseq/multiqc/human_alignment/multiqc_report.html'],
		# ['arion/illumina/s02-fastqc.dir/human/trimmed', 'arion/illumina/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		# ['arion/illumina/s02-fastqc.dir/mouse/trimmed', 'arion/illumina/multiqc/mouse_fastqc_trimmed/multiqc_report.html'],
		# ['arion/illumina/s04-alignment.dir/human', 'arion/illumina/multiqc/human_alignment_trimmed/multiqc_report.html'],
		# ['arion/illumina/s04-alignment.dir/mouse', 'arion/illumina/multiqc/mouse_alignment_trimmed/multiqc_report.html'],
	]
	for files in filelist:
		yield [files[0], files[1]]

@files(qcJobs)

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="01:00", GB=10, n=1, print_outfile=False, run_locally=False, stdout=outfile.replace('.html', '.log'))

#######################################################
#######################################################
########## S3. Alignment
#######################################################
#######################################################

#############################################
########## 1. Bowtie
#############################################

# @follows(trimIlluminaAdapters, buildBowtieIndex)

# @collate('arion/chipseq/s01-fastq.dir/human/trimmed/*/*_trimmed.fq.gz',
@collate('arion/chipseq/s01-fastq.dir/human/raw/*/*.fastq.gz',
 		 regex(r'(.*)/s01-fastq.dir/(.*)/raw/(.*)/.*.fastq.gz'),
		 add_inputs(r'arion/atacseq/s03-alignment.dir/\2/bowtie2/index/*primary_assembly.1.bt2'),
		 r'\1/s03-alignment.dir/\2/bowtie2/results/\3/\3.bam')

def runBowtie(infiles, outfile):

	# Infiles
	fastq = [x[0] for x in infiles]
	bowtie_index = infiles[0][1].replace('.1.bt2', '')
	outname = outfile.replace('.bam', '')

	# FASTQ string
	if len(fastq) == 1:
		fastq_str = '-U {fastq[0]}'.format(**locals())
	elif len(fastq) == 2:
		fastq_str = '-1 {fastq[0]} -2 {fastq[1]}'.format(**locals())

	# Command
	cmd_str = '''bowtie2 -x {bowtie_index} {fastq_str} -q -p 40 | samtools view -bS --threads 40 | samtools sort --threads 40 -o {outfile} && \
		samtools index -b {outfile} && samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f5 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > {outname}.mapq'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=10, GB=5, modules=['bowtie2/2.4.1', 'samtools/1.9'], stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'), print_cmd=False, ow=False)

# find arion/chipseq/s03-alignment.dir/human/bowtie2/results -name "*.log" | js

#############################################
########## 2. Filter
#############################################
# Removes duplicates, and any reads not mapping to chromosomes 1-19, X, Y.
# "[XS] == null" - multimappers
# "not unmapped"
# "not duplicate"
# "ref_id <= 21 and ref_id != 19" - chromosomes 1-19,X,Y (no chrM)
# "mapping_quality >= 30"
# fragment lengths is empty for single end data

# @transform(runBowtie,
@transform('arion/chipseq/s03-alignment.dir/human/bowtie2/results/*/*0.bam',
		   suffix('.bam'),
		   '_filtered.bam')

def filterBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "ref_id <= 21 and ref_id != 19 and not unmapped and not duplicate and mapping_quality >= 30" {infile} > {outfile} && \
		samtools index -b {outfile} && samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | awk '$9>0' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > {outname}.fragment_lengths.txt'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'], print_outfile=True)

#############################################
########## 3. Create BigWig
#############################################

@transform(filterBam,
		   suffix('.bam'),
		   add_inputs('/sc/arion/projects/GuccioneLab/genome-indices/hg38/blacklists/hg38-blacklist.v2.bed'),
		   '.bw')

def createBigWig(infiles, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=10 --skipNonCoveredRegions --numberOfProcessors=48 --normalizeUsing RPKM --blackListFileName {infiles[1]} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=8, GB=4, print_outfile=False)

#######################################################
#######################################################
########## S4. Peaks
#######################################################
#######################################################

#############################################
########## 1. MACS2
#############################################

# @collate(filterBam,
# 		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/(.*)_Rep.*/.*_filtered.bam'),
# 		 r'\1/s04-peaks.dir/\2/macs2/\3/\3_peaks.xls')

# def runMacs2(infiles, outfile):

# 	# Infiles
# 	infiles_str = ' '.join(infiles)
# 	basename = outfile.replace('_peaks.xls', '')
# 	organism = outfile.split('/')[-4].replace('human', 'hs').replace('mouse', 'mm')

# 	# Command
# 	cmd_str = ''' macs2 callpeak \
# 		-t {infiles_str}\
# 		-f BAM \
# 		--nomodel \
# 		--nolambda \
# 		--extsize 250 \
# 		-q 0.01 \
# 		-n {basename} \
# 		-g {organism}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['macs/2.1.0'], W='06:00', n=1, GB=30, print_cmd=False, stdout=outfile.replace('.xls', '.log'), stderr=outfile.replace('.xls', '.err'))

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
if __name__ == '__main__':
	cmdline.run(options)
print('Done!')