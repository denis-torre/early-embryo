#################################################################
#################################################################
############### Embryo Data Download ################
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
r_source = 'pipeline/scripts/data-download.R'
py_source = 'pipeline/scripts/DataDownload.py'
P = 'acc_GuccioneLab'
q = 'express'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import DataDownload as P

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
# Qiao
qiao_pacbio = '../arion/datasets/qiao/qiao-pacbio.csv'
qiao_illumina = '../arion/datasets/qiao/qiao-illumina.csv'

# Human
human_illumina = '../arion/datasets/human/human_embryo_ilmn/Groups.txt'
human_pacbio_flnc = '../arion/datasets/human/human_embryo_pacb/*.flnc.bam'

# Reference genomes
ensembl_json = '../arion/datasets/reference_genomes/ensembl-links.json'

#######################################################
#######################################################
########## S2. Mouse Embryo IsoSeq (Qiao et al., 2020)
#######################################################
#######################################################

#############################################
########## 1. PacBio
#############################################

@subdivide(qiao_pacbio,
		   regex(r'(.*)/.*-(.*).csv'),
		   r'\1/rawdata/\2/*/*.subreads.bam',
		   r'\1/rawdata/\2/{Run}')

def qiaoPacbio(infile, outfiles, outfileRoot):

	# Read metadata file
	metadata_dataframe = pd.read_csv(infile, comment='#')

	# Loop
	for index, rowData in metadata_dataframe.iterrows():

		# Get outdir
		outdir = outfileRoot.format(**rowData)

		# Download
		os.system('wget -P {outdir} {URL}'.format(**locals(), **rowData))
		# Manually renamed two bams from .bam.1 to .bam

#############################################
########## 2. Illumina
#############################################

@subdivide(qiao_illumina,
		   regex(r'(.*)/.*-(.*).csv'),
		   r'\1/rawdata/\2/*.fastq.gz',
		   r'\1/rawdata/\2/{run}*.fastq.gz')

def qiaoIllumina(infile, outfiles, outfileRoot):

	# Read metadata file
	metadata_dataframe = pd.read_csv(infile, comment='#').query('Instrument=="Illumina NovaSeq 6000"')

	# Loop
	for index, rowData in metadata_dataframe.iterrows():

		# Get platform
		run = rowData['Run']

		# Get outdir
		outfile_pattern = outfileRoot.format(**locals())
		outdir = os.path.dirname(outfile_pattern)

		# Run
		# if len(glob.glob(outfile_pattern)) == 0:

			# Create directory
			# if not os.path.exists(outdir):
			# 	os.makedirs(outdir)

		# 	# Download
		# 	os.system('ml sratoolkit/2.10.5 && cd {outdir} && fasterq-dump -e 6 -p --split-files {run}'.format(**locals())) # && gzip {run}*.fastq

		# 	# Get files
		# 	uncompressed_files = glob.glob(outfile_pattern.replace('.gz', ''))

		# 	# Loop
		# 	for uncompressed_file in uncompressed_files:

		# 		# Compress
		# 		cmd_str = 'gzip {uncompressed_file}'.format(**locals())

		# 		# Run
		# 		outfile = uncompressed_file+'.gz'
		# 		run_job(cmd_str, outfile, W='03:00', GB=6, n=1)

#######################################################
#######################################################
########## S. Public data download
#######################################################
#######################################################

#############################################
########## 1. SRA
#############################################

@subdivide(('../arion/datasets/blakeley/blakeley-samples.csv', '../arion/datasets/xue/xue-samples.csv', '../arion/datasets/yan/yan-samples.csv', '../arion/datasets/deng/deng-samples.csv'),
		   regex(r'(.*)/.*.csv'),
		   r'\1/rawdata/*.fastq.gzA',
		   r'\1/rawdata/{run}*.fastq.gz')

def downloadFASTQ(infile, outfiles, outfileRoot):

	# Read metadata file
	metadata_dataframe = pd.read_csv(infile, comment='#')

	# Loop
	for index, rowData in metadata_dataframe.iterrows():

		# Get platform
		run = rowData['Run']

		# Get outdir
		outfile_pattern = outfileRoot.format(**locals())
		outdir = os.path.dirname(outfile_pattern)

		# Run
		if len(glob.glob(outfile_pattern)) == 0:
			
			# Create directory
			if not os.path.exists(outdir):
				os.makedirs(outdir)

			# Download
			print('Doing {outfile_pattern}...'.format(**locals()))
			os.system('ml sratoolkit/2.10.5 && cd {outdir} && fasterq-dump -e 6 -p --split-files {run}'.format(**locals())) # && gzip {run}*.fastq

			# Get files
			uncompressed_files = glob.glob(outfile_pattern.replace('.gz', ''))

			# Loop
			for uncompressed_file in uncompressed_files:

				# Compress
				cmd_str = 'gzip {uncompressed_file}'.format(**locals())

				# Run
				outfile = uncompressed_file+'.gz'
				run_job(cmd_str, outfile, W='03:00', GB=6, n=1)

#######################################################
#######################################################
########## S. Xia
#######################################################

#############################################
########## 1. Get SRP IDs
#############################################
# SRP tables looked inconsistent and lacked sample names
# Manually selected samples in arion/datasets/xia/xia-samples_v2.csv from GEO 
# at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124718 on 2021/05/27

@transform('arion/datasets/xia/xia-samples_v2.csv',
		   suffix('.csv'),
		   '_metadata.tsv')

def convertSRX(infile, outfile):

	# Read
	sample_dataframe = pd.read_csv(infile)

	# Get IDs
	gsm_ids = ' '.join(sample_dataframe['geo_sample_id'])

	# Command
	os.system('pysradb gsm-to-srx --detailed --saveto {outfile} {gsm_ids}'.format(**locals()))

#############################################
########## 2. SRA
#############################################

# scr
# ml sratoolkit/2.10.5 
# fasterq-dump -e 6 -p --split-files SRR9131738

@subdivide(convertSRX,
		   regex(r'(.*)/.*.tsv'),
		   add_inputs(),
		   r'\1/rawdata/*.fastq.gz',
		   r'\1/rawdata/{experiment_accession}*.fastq.gz')

def downloadSRX(infile, outfiles, outfileRoot):

	# Read metadata file
	metadata_dataframe = pd.read_table(infile)

	# Loop
	for index, rowData in metadata_dataframe.iterrows():

		# Get platform
		experiment_accession = rowData['experiment_accession']

		# Get outdir
		outfile_pattern = outfileRoot.format(**locals())
		outdir = os.path.dirname(outfile_pattern)

		print(outfile_pattern)

		# # Run
		# if len(glob.glob(outfile_pattern)) == 0:
			
		# 	# Create directory
		# 	if not os.path.exists(outdir):
		# 		os.makedirs(outdir)

		# 	# Download
		# 	print('Doing {outfile_pattern}...'.format(**locals()))
		# 	os.system('ml sratoolkit/2.10.5 && cd {outdir} && fasterq-dump -e 6 -p --split-files {run}'.format(**locals())) # && gzip {run}*.fastq

		# 	# Get files
		# 	uncompressed_files = glob.glob(outfile_pattern.replace('.gz', ''))

		# 	# Loop
		# 	for uncompressed_file in uncompressed_files:

		# 		# Compress
		# 		cmd_str = 'gzip {uncompressed_file}'.format(**locals())

		# 		# Run
		# 		outfile = uncompressed_file+'.gz'
		# 		run_job(cmd_str, outfile, W='03:00', GB=6, n=1)

#######################################################
#######################################################
########## S. Evolutionary conservation
#######################################################
#######################################################

#############################################
########## 1. PhyloP
#############################################

def phylopJobs():
	bw_files = ['http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/mm10.60way.phyloP60way.bw', 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw']
	for url_path in bw_files:
		filename = os.path.basename(url_path)
		genome = filename.split('.')[0]
		organism = genome.replace('mm10', 'mouse').replace('hg38', 'human')
		outfile = 'arion/datasets/phylop/{organism}/{filename}'.format(**locals())
		yield [url_path, outfile]

@files(phylopJobs)

def downloadPhyloP(url_path, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Download
	os.system('cd {outdir} && wget {url_path}'.format(**locals()))

#######################################################
#######################################################
########## S. Xue
#######################################################
#######################################################

#############################################
########## 1. 
#############################################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183


#######################################################
#######################################################
########## S. Yan
#######################################################
#######################################################

#############################################
########## 1. 
#############################################

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552

#######################################################
#######################################################
########## S2. Reference Genomes
#######################################################
#######################################################

#############################################
########## 1. Download
#############################################

@subdivide(ensembl_json,
		   regex(r'(.*)/.*.json'),
		   r'\1/*/*',
		   r'\1/{organism}/{filename}')

def downloadGenomes(infile, outfiles, outfileRoot):

	# Read
	with open(infile) as openfile:
		ensembl_dict = json.load(openfile)

	# Loop
	for organism, file_dict in ensembl_dict.items():
		for file_type, url in file_dict.items():
			
			# Get outfile
			filename = os.path.basename(url)
			outfile = outfileRoot.format(**locals())

			# Get directory
			outdir = os.path.dirname(outfile)
			if not os.path.exists(outdir):
				os.makedirs(outdir)

			# # Check
			# if len(outfiles) == 0:
			
			# 	# Command
			# 	os.system('wget -P {outdir} {url}'.format(**locals()))

			# 	# Command
			# 	cmd_str = 'gunzip {outfile}'.format(**locals())
			
			# 	# Run
			# 	run_job(cmd_str, outfile.replace('.gz', ''), W='00:30', GB=15)

#############################################
########## 2. Create JSON
#############################################	

@merge(downloadGenomes,
	   '../arion/datasets/reference_genomes/reference-genomes.json')

def createReferenceJson(infiles, outfile):

	# Initialize
	result_dict = {x:{} for x in ['human', 'mouse']}

	# Loop
	for infile in infiles:

		# Get organism
		organism = infile.split('/')[-2]

		# Get file type
		if 'gtf' in infile:
			file_type = 'gtf'
		elif 'cdna' in infile:
			file_type = 'cdna'
		elif 'primary_assembly' in infile:
			file_type = 'primary_assembly'

		# Add
		result_dict[organism][file_type] = os.path.join(os.path.dirname(os.getcwd()), infile.replace('../', ''))

	# Write
	with open(outfile, 'w') as openfile:
		json.dump(result_dict, openfile, indent=4)

#######################################################
#######################################################
########## S3. Sample renaming
#######################################################
#######################################################

#############################################
########## 1. Qiao
#############################################

@transform(qiao_illumina,
		   suffix('-illumina.csv'),
		   add_inputs(qiao_pacbio),
		   '-sample_names.csv')

def renameQiao(infiles, outfile):
	
	# Read data
	sample_dataframe = pd.read_csv(infiles[0], comment='#').sort_values(['Platform', 'Cell_type', 'Run'])
	sample_dataframe = sample_dataframe[sample_dataframe['Center Name']=='GEO'][['Run', 'Platform', 'Cell_type', 'Strain']]

	# Rename
	sample_dataframe['count'] = sample_dataframe.groupby(['Platform', 'Cell_type']).cumcount()+1
	sample_dataframe['sample_id'] = sample_dataframe['Run']
	sample_dataframe['sample_name'] = ['mouse_'+rowData['Cell_type'].replace(' cells', '').replace('cell', 'C')+('_Rep{}'.format(rowData['count']) if rowData['Platform'] == 'ILLUMINA' else '') for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 2. Human PacBio
#############################################

@transform(human_illumina,
		   regex(r'(.*)/human_embryo_ilmn/Groups.txt'),
		   add_inputs(human_pacbio_flnc),
		   r'\1/human-sample_names.csv')

def renameHumanPacbio(infiles, outfile):
	
	# Read data
	illumina_dataframe = pd.read_csv(infiles[0], sep='\t').sort_values('batch').rename(columns={'sample_name': 'sample_id'})
	illumina_dataframe['Platform'] = 'ILLUMINA'

	# Illumina
	illumina_dataframe['count'] = illumina_dataframe.groupby(['batch', 'Group']).cumcount()+1
	illumina_dataframe = illumina_dataframe.sort_values(['batch', 'Group', 'count'])
	illumina_dataframe['sample_name'] = ['human_{Group}_B{batch}_{count}'.format(**rowData).replace('cell', 'C') for index, rowData in illumina_dataframe.iterrows()]
	illumina_dataframe.head()

	# PacBio
	pacbio_dataframe = pd.DataFrame([
		{'sample_id': 'TD01042_1PN', 'Group': '2PN', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_2PN_1'},
		{'sample_id': 'TD01042_4PN', 'Group': '2PN', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_2PN_4'},
		{'sample_id': 'TD01042_EMB1', 'Group': '2PN', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_2PN_EMB1'},
		{'sample_id': 'TD01042_EMB2', 'Group': '2C', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_2C_EMB2'},
		{'sample_id': 'TD01042_EMB3', 'Group': '4C', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_4C_EMB3'},
		{'sample_id': 'TD01042_EMB4', 'Group': '4C', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_4C_EMB4'},
		{'sample_id': 'TD01042_EMB5', 'Group': '8C', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_8C_EMB5'},
		{'sample_id': 'TD01042_EMB6', 'Group': 'blastocyst', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_blastocyst_EMB6'},
		{'sample_id': 'TD01042_MOR', 'Group': 'morula', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_morula_MOR'},
		{'sample_id': 'TD01673_all_embryos_merged', 'Group': 'merged', 'Platform': 'PACBIO_SMRT', 'sample_name': 'human_merged_TD01673'}
	])

	# Merged
	concatenated_dataframe = pd.concat([illumina_dataframe, pacbio_dataframe], axis=0)
	
	# Write
	concatenated_dataframe.to_csv(outfile, index=False)


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