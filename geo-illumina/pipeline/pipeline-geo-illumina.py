#################################################################
#################################################################
############### GEO Illumina RNA-Seq ################
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
r_source = 'pipeline/scripts/geo-illumina.R'
py_source = 'pipeline/scripts/GeoIllumina.py'
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
#import GeoIllumina as P

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
linked_fastq = glob.glob('arion/geo_illumina/s01-datasets.dir/*/fastq/*/*.fastq.gz')
human_filtered_gtf = 'arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf'
human_star_index = 'arion/illumina/s04-alignment.dir/human/all/STAR/index'
human_rsem_index = 'arion/illumina/s04-alignment.dir/human/all/RSEM/index/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.idx.fa'
human_junction_file = 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_junctions.tsv'

#######################################################
#######################################################
########## S1. Fix sample names
#######################################################
#######################################################

#############################################
########## 1. Yan
#############################################

@transform('arion/datasets/yan/yan-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixYanSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'GEO_Accession (exp)': 'geo_sample_id', 'Cell_type': 'geo_cell_type', 'Development_stage': 'geo_development_stage', 'Library Name': 'geo_library_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_development_stage == "early blastomere"')

	# Add cell type
	sample_dataframe['cell_type'] = [x.split(':')[-1].replace('-cell embryo', 'C').replace('Morulae', 'morula').replace('Late blastocyst', 'blastocyst').replace('Zygote', '1C').replace('-cell', 'C').replace('Oocyte', 'oocyte').replace(' ', '') for x in sample_dataframe['geo_library_name']]

	# Sort and add sample number
	sample_dataframe = sample_dataframe.sort_values('cell_type')
	sample_dataframe['sample_number'] = sample_dataframe.groupby('cell_type').cumcount()+1

	# Add sample name
	sample_dataframe['sample_name'] = ['human_{cell_type}_Rep{sample_number}'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 2. Yan
#############################################

@transform('arion/datasets/xue/xue-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixXueSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
		
	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'GEO_Accession (exp)': 'geo_sample_id', 'Cell_type': 'geo_cell_type', 'Organism': 'geo_organism', 'source_name': 'geo_source_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_organism == "Homo sapiens" and geo_cell_type != "primary blood"')

	# Add cell type
	sample_dataframe['cell_type'] = [x.split(':')[-1].replace('-cell blastomere', 'C').replace('zygote', '1C') for x in sample_dataframe['geo_cell_type']]
	sample_dataframe


	# Sort and add sample number
	sample_dataframe = sample_dataframe.sort_values('cell_type')
	sample_dataframe['sample_number'] = sample_dataframe.groupby('cell_type').cumcount()+1

	# Add sample name
	sample_dataframe['sample_name'] = ['human_{cell_type}_Rep{sample_number}'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 3. Liu
#############################################

@transform('arion/datasets/liu/liu-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixLiuSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'Sample Name': 'geo_sample_name', 'Organism': 'geo_organism', 'Library Name': 'geo_library_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_organism == "Homo sapiens" and assay_type == "RNA-Seq"')

	# Add sample name
	sample_dataframe['sample_name'] = ['human_'+x.replace('R_', '').replace('-cell', 'C').replace('Morula', 'morula').replace('Zygote', '1C') for x in sample_dataframe['geo_sample_name']]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 4. Link
#############################################

@subdivide((fixYanSamples, fixXueSamples, fixLiuSamples),
		   regex(r'(.*)/.*-sample_names.csv'),
		   r'\1/fastq/*/*.fastq.gz',
		   r'\1/fastq/{sample_name}/{fastq_basename}')

def linkFASTQ(infile, outfiles, outfileRoot):

	# Read
	sample_dataframe = pd.read_csv(infile)

	# Get dataset name
	dataset = os.path.basename(infile).split('-')[0]

	# Loop through dataframe
	for index, rowData in sample_dataframe.iterrows():

		# Get FASTQ files
		fastq_files = glob.glob('/hpc/users/torred23/pipelines/projects/early-embryo/geo-illumina/arion/datasets/{dataset}/rawdata/{Run}*.fastq.gz'.format(**locals(), **rowData))

		# Loop
		for fastq_file in fastq_files:
			pass

			# # Get basename
			# fastq_basename = os.path.basename(fastq_file)

			# # Get outfile
			# outfile = outfileRoot.format(**locals(), **rowData)
		
			# # Outdir
			# outdir = os.path.dirname(outfile)
			# if not os.path.exists(outdir):
			# 	os.makedirs(outdir)

			# # Link
			# os.system('ln -s {fastq_file} {outfile}'.format(**locals()))

#######################################################
#######################################################
########## S2. FASTQC
#######################################################
#######################################################

#############################################
########## 1. Run
#############################################

# @follows(linkFASTQ)

@transform(linked_fastq,
		   regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/(.*).fastq.gz'),
		   r'\1/s02-fastqc.dir/\2/\4_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

#############################################
########## 2. MultiQC
#############################################

# @follows(runFastQC)

@transform('arion/geo_illumina/s02-fastqc.dir/*',
		   regex(r'(.*)/s02-fastqc.dir/(.*)'),
		   r'\1/multiqc/\2/multiqc_report.html')

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="01:00", GB=10, n=1, print_outfile=False, run_locally=False, stdout=outfile.replace('.html', '.log'))

#######################################################
#######################################################
########## S3. Align
#######################################################
#######################################################

#############################################
########## 1. STAR junctions
#############################################

# @follows(linkFASTQ)

@collate(linked_fastq,
		 regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/.*.fastq.gz'),
		 add_inputs(human_star_index),
		 r'\1/s03-alignment.dir/\2/STAR/pass1/\3/\3-SJ.out.tab')

def getStarJunctions(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]

	# FASTQ string
	fastq_str = ' '.join(fastq_files)

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", print_outfile=True, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass1/*/*_job.log | jsc

#############################################
########## 2. STAR BAM
#############################################

@follows(getStarJunctions)

@collate(linked_fastq,
		 regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/.*.fastq.gz'),
		 add_inputs(human_star_index, r'\1/s03-alignment.dir/\2/STAR/pass1/*/*-SJ.out.tab'),
		 r'\1/s03-alignment.dir/\2/STAR/pass2/\3/\3-Aligned.sortedByCoord.out.bam')

def runStar(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]
	sj_files = infiles[0][2:]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	fastq_str = ' '.join(fastq_files)
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=15, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_outfile=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*.log | jsc
# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*Log.out | lr
# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*.err | xargs wc -l

#############################################
########## 4. Junction counts
#############################################

@follows(runStar)

@collate('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-SJ.out.tab',
		 regex(r'(.*)/(.*)/(.*)/pass2/.*/.*.tab'),
		 add_inputs(human_junction_file),
		 r'\1/\2/\3/\2-junction_counts.tsv')

def getJunctionCounts(infiles, outfile):

	# Infiles
	sj_files = [x[0] for x in infiles]
	junction_file = infiles[0][1]

	# Run
	run_r_job('get_junction_counts', sj_files, outfile, additional_params=junction_file, print_outfile=False, conda_env='env', W='02:00', GB=50, n=3, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# find arion/geo_illumina/s03-alignment.dir -name "*junction_counts*" | xargs rm

#############################################
########## 5. RSEM expression
#############################################

@follows(runStar)

# # @transform('arion/illumina/s04-alignment.dir/human/2C_vs_4C/STAR/pass2/human_4C_B3_9/human_4C_B3_9-Aligned.toTranscriptome.out.bam',
# @transform('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/human_morula_Rep15-Aligned.toTranscriptome.out.bam',
@transform('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
# @transform('arion/geo_illumina/s03-alignment.dir/xue/STAR/pass2/human_1C_Rep1/*-Aligned.toTranscriptome.out.bam',
		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
		   add_inputs(human_rsem_index),
		   r'\1/RSEM/\2/\2.isoforms.results')

def runRsem(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Paired end
	dataset = outfile.split('/')[-4]
	paired_str = '--paired-end' if dataset == 'xue' else ''
	strandedness_str = 'forward' if dataset == 'xue' else 'none'

	# Command
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness {strandedness_str} \
		{paired_str} \
		--estimate-rspd \
		--num-threads 200 \
		--no-bam-output \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log && \
		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())
		# --calc-ci \

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=2, n=25, modules=['rsem/1.3.3'], print_outfile=True, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

# find arion/geo_illumina/s03-alignment.dir/*/RSEM -name "*.log" | grep -v 'rsem' | jsc

#############################################
########## 4. Create BigWig
#############################################

@transform(runStar,
		   suffix('-Aligned.sortedByCoord.out.bam'),
		   '.bw')

def createBigWig(infile, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

#######################################################
#######################################################
########## S5. Expression
#######################################################
#######################################################

#############################################
########## 1. Prepare metadata
#############################################

# @collate('arion/geo_illumina/s03-alignment.dir/*/RSEM/*/*.isoforms.results',
@collate(runRsem,
		 regex(r'(.*)/s03-alignment.dir/(.*)/RSEM/.*.isoforms.results'),
		 r'\1/s04-expression.dir/\2/\2-sample_metadata.txt')

def prepareSampleMetadata(infiles, outfile):

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Prepare dict
	sample_dict = [{'files': x, 'names': x.split('/')[-2], 'cell_type': x.split('/')[-2].split('_')[1]} for x in infiles]

	# Convert to dataframe
	sample_dataframe = pd.DataFrame(sample_dict)

	# Write
	sample_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Aggregate
#############################################

@transform(prepareSampleMetadata,
		   suffix('-sample_metadata.txt'),
		   add_inputs(human_filtered_gtf),
		   '-counts.rda')

def aggregateCounts(infiles, outfile):

	# Run
	run_r_job('aggregate_counts', infiles, outfile, conda_env='env', run_locally=False, W='00:30', GB=30, n=1, wait=False, ow=True, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#######################################################
#######################################################
########## S5. Sashimi
#######################################################
#######################################################

#############################################
########## 1. BAM tables
#############################################

# @collate(runStar,
# @collate('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
# @collate('arion/illumina/s04-alignment.dir/human/all/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
@collate(('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam', 'arion/illumina/s04-alignment.dir/human/all/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam'),
		 regex(r'.*.dir/(.*)/STAR/pass2/.*/.*.bam'),
		 r'arion/geo_illumina/summary.dir/sashimi/settings/\1-bams.txt')

def makeBamTables(infiles, outfile):
		
	# Create dataframe
	bam_dataframe = pd.DataFrame([{'sample_name': x.split('/')[-2], 'bam': os.path.join(os.getcwd(), x), 'cell_type': x.split('/')[-2].split('_')[1]} for x in infiles])

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	bam_dataframe.to_csv(outfile, sep='\t', index=False, header=False)

#############################################
########## 2. Sashimi plots
#############################################

# def sashimiJobs():

# 	# Loop through comparisons
# 	for comparison, comparison_info in comparison_dict.items():
		
# 		# Get BAM paths
# 		bams = {}
# 		for group in comparison_info['contrast']:
# 			with open('arion/illumina/s07-rmats.dir/isoseq/{comparison}/bams/{group}-bams.txt'.format(**locals())) as openfile:
# 				bams[group] = [os.path.join(os.getcwd(), x) for x in openfile.read().split(',')]
				
# 		# Read splicing results
# 		splicing_dataframe = pd.read_table('arion/illumina/summary.dir/{comparison}/{comparison}-isoseq-splicing_summary.tsv'.format(**locals()))
# 		splicing_dataframe['i'] = splicing_dataframe['pval_rmats'].sort_values().index

# 		# Filter
# 		splicing_dataframe = splicing_dataframe.query('significant_rmats==True and significant_suppa==True and event_type == "SE"').sort_values('pval_rmats')
		
# 		# Loop through splicing results
# 		for index, rowData in splicing_dataframe.iterrows():
# 			min_coverage = 5
# 			plot_settings = {
# 				'groups': [{'name': group, 'bams': bams[group], 'psi': round(rowData['PSI'+str(i+1)], ndigits=2), 'color': '#377eb8' if i == 0 else '#e41a1c'} for i, group in enumerate(comparison_info['contrast'])],
# 				'event_id': rowData['Event_id'],
# 				'gene_name': rowData['gene_name'],
# 				'rmats_statistics': 'dPSI={dPSI_rmats:.2f}, FDR={pval_rmats:.1e}'.format(**rowData),
# 				'gtf': reference_dict['isoseq']['gtf_cds'],
# 				'min_coverage': min_coverage
# 			}

# 			# Files
# 			infile = 'arion/illumina/summary.dir/{comparison}/plots/sashimi/{gene_name}-{event_type}-{i}-{min_coverage}.json'.format(**locals(), **rowData)
# 			outfile = infile.replace('.json', '.pdf')
			
# 			# Outdir
# 			outdir = os.path.dirname(outfile)
# 			if not os.path.exists(outdir):
# 				os.makedirs(outdir)

# 			# Write infile
# 			with open(infile, 'w') as openfile:
# 				json.dump(plot_settings, openfile)

# 			yield [infile, outfile]


# # @follows(functionToFollow)

# @files(sashimiJobs)

# def sashimiPlot(infile, outfile):

# 	# Run
# 	run_py_job('ggsashimi', infile, outfile, W="00:10", GB=30, n=1, run_locally=False, conda_env='env')

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