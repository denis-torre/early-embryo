#################################################################
#################################################################
############### Embryo methylation ################
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
r_source = 'pipeline/scripts/embryo-methylation.R'
py_source = 'pipeline/scripts/EmbryoMethylation.py'
P = 'acc_apollo'
q = 'sla'
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
#import EmbryoMethylation as P

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
sample_file = 'arion/datasets/guo/guo-samples.csv'
name_file = 'arion/datasets/guo/guo-sample_names.csv'
all_illumina_fastq = glob.glob('arion/methylation/s01-fastq.dir/human/*/*/*.f*q.gz')

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():
	sample_dataframe = pd.read_csv(sample_file, comment='#').rename(columns={'GEO_Accession (exp)': 'geo_sample_id'})[['Run', 'geo_sample_id']]
	name_dataframe = pd.read_table(name_file)
	name_dataframe['sample_name'] = [x.replace('RRBS', 'human').replace('-cell', 'C_').replace('ICM', 'ICM_').replace('TE', 'TE_').replace('Zygote', '1C_').replace('MII_Oocyte', 'oocyte_').replace('Morula', 'morula_').replace('Postimplantation_embryo', 'postimplantation_').replace('rep', 'Rep') for x in name_dataframe['sample_name']]
	sample_dict = sample_dataframe.merge(name_dataframe, on='geo_sample_id').set_index('Run')['sample_name'].to_dict()
	fastq_files = glob.glob('arion/datasets/guo/rawdata/*.fastq.gz')
	for fastq_file in fastq_files:
		fastq_name = os.path.basename(fastq_file)
		run = fastq_name[:-len('_1.fastq.gz')]
		sample_name = sample_dict[run]
		infile = os.path.join(os.getcwd(), fastq_file)
		outfile = 'arion/methylation/s01-fastq.dir/human/raw/{sample_name}/{fastq_name}'.format(**locals())
		yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	if not os.path.exists(outfile):
		print(outfile)
		# os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/methylation/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

# @follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	cmd_str = '''trim_galore --paired --rrbs --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='05:00', GB=6, n=6, print_outfile=False, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

# @follows(linkFASTQ, trimIlluminaAdapters, trimReads)

@transform(all_illumina_fastq,
		   regex(r'(.*)/s01-fastq.dir/(.*).f.*q.gz'),
		   r'\1/s02-fastqc.dir/\2_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

# find arion/chipseq/s02-fastqc.dir/human -name "*fastqc.html" | wc -l

#############################################
########## 2. MultiQC
#############################################

@follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/methylation/s02-fastqc.dir/human/raw', 'arion/methylation/multiqc/human_fastqc/multiqc_report.html'],
		['arion/methylation/s02-fastqc.dir/human/trimmed', 'arion/methylation/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/methylation/s02-fastqc.dir/human', 'arion/methylation/multiqc/human_fastqc_all/multiqc_report.html'],
		# ['arion/chipseq/s02-fastqc.dir/human/trimmed_50bp', 'arion/chipseq/multiqc/human_fastqc_trimmed_50bp/multiqc_report.html'],
		# ['arion/chipseq/s03-alignment.dir/human', 'arion/chipseq/multiqc/human_alignment_trimmed_50bp/multiqc_report.html']
		# ['arion/chipseq/s03-alignment.dir/human', 'arion/chipseq/multiqc/human_alignment_trimmed/multiqc_report.html']
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
########## S3. Bismark
#######################################################
#######################################################

#############################################
########## 1. Prepare genome
#############################################

@transform('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
		   regex(r'.*/reference_genomes/(.*)/.*.fa'),
		   r'arion/methylation/s03-bismark.dir/\1/genome/Bisulfite_Genome')

def prepareBismarkGenome(infile, outfile):

	# Directory
	genome_dir = os.path.dirname(outfile)

	# Command
	cmd_str = ''' cp {infile} {genome_dir} && bismark_genome_preparation {genome_dir} --verbose '''.format(**locals()) #--parallel 4

	# Run
	run_job(cmd_str, outfile, print_outfile=True, modules=['bismark/0.22.3', 'bowtie2/2.4.1'], W='06:00', GB=50, n=1, stdout=outfile.replace('.fa', '.log'), stderr=outfile.replace('.fa', '.err'))#, print_outfile=False, print_cmd=False)

#############################################
########## 2. Run Bismark
#############################################

@collate('arion/methylation/s01-fastq.dir/human/trimmed/*/*.fq.gz',
		 regex(r'(.*)/s01-fastq.dir/(.*?)/.*/(.*)/.*.fq.gz'),
		 add_inputs(r'\1/s03-bismark.dir/\2/genome/'),
		 r'\1/s03-bismark.dir/\2/results/\3/\3_pe.bam')

def runBismark(infiles, outfile):

	# Directory
	output_dir = os.path.dirname(outfile)
	temp_dir = os.path.join(output_dir, 'tmp')
	basename = os.path.basename(outfile)[:-len('_pe.bam')]

	# Command
	cmd_str = '''  bismark --genome_folder {infiles[0][1]} -1 {infiles[0][0]} -2 {infiles[1][0]} --output_dir {output_dir} --temp_dir {temp_dir} --basename {basename} '''.format(**locals()) # --parallel 4 

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bismark/0.22.3', 'bowtie2/2.4.1', 'samtools/1.13'], W='10:00', GB=50, n=1, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))#, print_outfile=False, print_cmd=False)

# find arion/methylation/s03-bismark.dir/human/results -name "*.log" | js

#############################################
########## 3. Extract methylation
#############################################

@transform('arion/methylation/s03-bismark.dir/human/results/human_8C_1_Rep1/human_8C_1_Rep1_pe.bam',
		   regex(r'(.*)/(.*)'),
		   r'\1/methylation/\2_splitting_report.txt')

def extractMetylation(infile, outfile):
	
	# Directory
	output_dir = os.path.dirname(outfile)

	# Command
	cmd_str = '''  bismark_methylation_extractor --gzip --bedGraph --output {output_dir} {infile} '''.format(**locals()) # --parallel 4 

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bismark/0.22.3', 'bowtie2/2.4.1', 'samtools/1.13'], W='06:00', GB=30, n=1, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))#, print_outfile=False, print_cmd=False)

#############################################
########## 4. Report
#############################################

# @follows(runBismark)

@transform('arion/methylation/s03-bismark.dir/human/results/*/*_report.txt',
		   regex(r'(.*)/results/(.*)/.*.txt'),
		   r'\1/summary/reports/\2_report.txt')

def createBismarkReport(infile, outfile):

	# Output dir
	output_dir = os.path.dirname(outfile)

	# Command
	cmd_str = ''' bismark2report --alignment_report {infile} --dir {output_dir} '''.format(**locals()) # --parallel 4 

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bismark/0.22.3', 'bowtie2/2.4.1', 'samtools/1.13'], W='02:00', GB=30, n=1, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))#, print_outfile=False, print_cmd=False)

# #############################################
# ########## 5. Summary
# #############################################

# @transform('arion/methylation/s03-bismark.dir/human/summary/reports',
# 		   regex(r'(.*).suffix'),
# 		   r'\1.new_suffix')

# def createBismarkSummary(infile, outfile):
# 	print(infile, outfile)

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