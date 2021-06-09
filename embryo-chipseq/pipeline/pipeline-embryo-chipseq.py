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
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

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