#################################################################
#################################################################
############### Embryo Illumina ################
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
import re
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/embryo-illumina.R'
py_source = 'pipeline/scripts/EmbryoIllumina.py'
P = 'acc_apollo'
q = 'sla'
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
sys.path.append('pipeline/scripts')
import EmbryoIllumina as S

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
# with open('/sc/arion/projects/GuccioneLab/genome-indices/genome-indices.json') as openfile:
# 	genome_indices = json.load(openfile)

##### 3. Variables #####
# Metadata
organism_metadata = {
	'mouse': 'arion/datasets/qiao/qiao-sample_names.csv',
	'human': 'arion/datasets/human/human-sample_names.csv'
}

# Qiao
# qiao_metadata = 'arion/datasets/qiao/qiao-sample_names.csv'
qiao_illumina = glob.glob('arion/datasets/qiao/rawdata/illumina/*.fastq.gz')

# Human
# human_metadata = 'arion/datasets/human/human-sample_names.csv'
human_illumina = glob.glob('arion/datasets/human/human_embryo_ilmn/*/*/*.fastq.gz')

# All
all_illumina_fastq = glob.glob('arion/illumina/s01-fastq.dir/*/*/*/*.f*q.gz')
trimmed_illumina_fastq = glob.glob('arion/illumina/s01-fastq.dir/*/trimmed/*/*.fq.gz')

# Outlier samples
outlier_sample_file = 'pipeline/outlier_samples.json'
with open(outlier_sample_file) as openfile:
	outlier_samples = json.load(openfile)

# References
reference_dict = {
	'human': {
		# 'ensembl': {
		# 	'genome_fasta': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
		# 	'gtf': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.102.gtf',
		# 	'transcript_fasta': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.cdna.all.fa'
		# },
		'isoseq': {
			'genome_fasta': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
			'gtf': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon.gtf',
			'gtf_filtered': 'arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf',
			'talon_abundance': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_abundance_filtered.tsv',
			'talon_junctions': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_junctions.tsv',
			'transcript_fasta': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon.fasta',
			'summary_file': 'arion/isoseq/summary.dir/human-isoseq_summary.tsv',
			'cpat_predictions': 'arion/isoseq/s06-cpat.dir/human/Homo_sapiens.GRCh38.102_talon-cpat.ORF_prob.best_results.tsv',
			'pfam_predictions': 'arion/isoseq/s07-pfam.dir/human/human-translated_pfam.tsv',
			'gtf_cds': 'arion/isoseq/s06-cpat.dir/human/gtf/Homo_sapiens.GRCh38.102_talon.cds.gtf'
		}
	},
	'mouse': {
		# 'ensembl': {
		# 	'genome_fasta': 'arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa',
		# 	'gtf': 'arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.102.gtf',
		# 	'transcript_fasta': 'arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.cdna.all.fa'
		# },
		'isoseq': {
			'genome_fasta': 'arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa',
			'gtf': 'arion/isoseq/s05-talon.dir/mouse/Mus_musculus.GRCm38.102_talon.gtf',
			'gtf_filtered': 'arion/illumina/s04-alignment.dir/mouse/all/gtf/Mus_musculus.GRCm38.102_talon-all-SJ_filtered.gtf',
			'talon_abundance': 'arion/isoseq/s05-talon.dir/mouse/Mus_musculus.GRCm38.102_talon_abundance_filtered.tsv',
			'talon_junctions': 'arion/isoseq/s05-talon.dir/mouse/Mus_musculus.GRCm38.102_talon_junctions.tsv',
			'transcript_fasta': 'arion/isoseq/s05-talon.dir/mouse/Mus_musculus.GRCm38.102_talon.fasta',
			'summary_file': 'arion/isoseq/summary.dir/mouse-isoseq_summary.tsv',
			'cpat_predictions': 'arion/isoseq/s06-cpat.dir/mouse/Mus_musculus.GRCm38.102_talon-cpat.ORF_prob.best_results.tsv',
			'pfam_predictions': 'arion/isoseq/s07-pfam.dir/mouse/mouse-translated_pfam.tsv',
			'gtf_cds': 'arion/isoseq/s06-cpat.dir/mouse/gtf/Mus_musculus.GRCm38.102_talon.cds.gtf'
		}
	}
}

# Comparisons
comparison_file = 'pipeline/comparisons.json'
with open(comparison_file) as openfile:
	comparison_dict = json.load(openfile)

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():

	# Read metadata
	mouse_dataframe = pd.read_csv(organism_metadata['mouse']).query('Platform == "ILLUMINA"')
	human_dataframe = pd.read_csv(organism_metadata['human']).query('Platform == "ILLUMINA"')
	
	# Samples
	sample_dict = {
		'mouse': [{'sample_name': rowData['sample_name'], 'fastq': [x for x in qiao_illumina if rowData['Run'] in x]} for  index, rowData in mouse_dataframe.iterrows()],
		'human': [{'sample_name': rowData['sample_name'], 'fastq': [x for x in human_illumina if x.split('/')[-3] == rowData['sample_id']]} for index, rowData in human_dataframe.iterrows()]
	}

	# Loop
	for organism, samples in sample_dict.items():
		for sample in samples:
			for fastq_file in sample['fastq']:

				# Get FASTQ name
				fastq_basename = os.path.basename(fastq_file)

				# Add infile path
				infile = os.path.join(os.getcwd(), fastq_file)

				# Get outfile
				outfile = 'arion/illumina/s01-fastq.dir/{organism}/raw/{sample_name}/{fastq_basename}'.format(**locals(), **sample)

				# Yield
				yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	# os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/illumina/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

@follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	cmd_str = '''trim_galore --nextera --paired --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='06:00', GB=6, n=6, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

# find arion/atacseq/s01-fastq.dir/*/trimmed -name "job.log" | jsc

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

@follows(linkFASTQ, trimIlluminaAdapters)

@transform(all_illumina_fastq,
		   regex(r'(.*)/s01-fastq.dir/(.*).f.*q.gz'),
		   r'\1/s02-fastqc.dir/\2_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

#############################################
########## 2. MultiQC
#############################################

# @follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/illumina/s02-fastqc.dir/human/raw', 'arion/illumina/multiqc/human_fastqc/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/mouse/raw', 'arion/illumina/multiqc/mouse_fastqc/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/human/trimmed', 'arion/illumina/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/mouse/trimmed', 'arion/illumina/multiqc/mouse_fastqc_trimmed/multiqc_report.html'],
		['arion/illumina/s04-alignment.dir/human', 'arion/illumina/multiqc/human_alignment_trimmed/multiqc_report.html'],
		['arion/illumina/s04-alignment.dir/mouse', 'arion/illumina/multiqc/mouse_alignment_trimmed/multiqc_report.html'],
		['arion/illumina/s04-alignment.dir-old/human/isoseq/all', 'arion/illumina/multiqc/human_test_qc/multiqc_report.html'],
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
########## S3. Splice junctions
#######################################################
#######################################################

#############################################
########## 1. STAR index
#############################################

def starIndexJobs():
	for organism, organism_references in reference_dict.items():
		for source, reference_files in organism_references.items():
			outfile = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/index'.format(**locals())
			yield [reference_files, outfile]

@follows(runFastQC)

@files(starIndexJobs)

def buildStarIndex(infiles, outfile):

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {genome_fasta} --sjdbGTFfile {gtf} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals(), **infiles)

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=False, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

# find arion/illumina/s03-junctions.dir -name "job.log" | js

#############################################
########## 2. STAR junctions
#############################################

def starJunctionJobs():
	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()
	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()
		for sample_name, fastq_files in fastq_dict.items():
			for source in ['isoseq']:
				star_index = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/index'.format(**locals())
				outfile = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/pass1/{sample_name}/{sample_name}-SJ.out.tab'.format(**locals())
				yield [(fastq_files, star_index), outfile]

@follows(linkFASTQ, buildStarIndex)

@files(starJunctionJobs)

def getStarJunctions(infiles, outfile):

	# Split
	fastq_files, star_index = infiles

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# ls arion/illumina/s03-junctions.dir/*/*/STAR/pass1/*/*_job.log | jsc

#############################################
########## 3. STAR BAM
#############################################

def starJobs():
	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()
	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()
		for sample_name, fastq_files in fastq_dict.items():
			for source in ['isoseq']:
				star_index = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/index'.format(**locals())
				sj_files = glob.glob('arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/pass1/*/*-SJ.out.tab'.format(**locals()))
				outfile = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/pass2/{sample_name}/{sample_name}-Aligned.sortedByCoord.out.bam'.format(**locals())
				yield [(fastq_files, star_index, sj_files), outfile]

@follows(getStarJunctions)

@files(starJobs)

def runStar(infiles, outfile):

	# Split
	fastq_files, star_index, sj_files = infiles

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# ls arion/illumina/s03-junctions.dir/*/*/STAR/pass2/*/*Aligned.sortedByCoord.out.log | jsc

#############################################
########## 4. Junction counts
#############################################

def sjCountJobs():
	for organism, organism_references in reference_dict.items():
		for source, reference_files in organism_references.items():
			infiles = [reference_files['talon_junctions']] + glob.glob('arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/pass2/*/*-SJ.out.tab'.format(**locals()))
			outfile = 'arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/{organism}-{source}-junction_counts.tsv'.format(**locals())
			yield [infiles, outfile]

@follows(runStar)

@files(sjCountJobs)

def getJunctionCounts(infiles, outfile):

	# Run
	run_r_job('get_junction_counts', infiles, outfile, run_locally=False, conda_env='env', W='02:00', GB=50, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S4. Alignment
#######################################################
#######################################################

#############################################
########## 1. Filter GTF
#############################################

def filterJobs():
	for organism, comparisons in comparison_dict.items():
		source = 'isoseq'
		reference_files = reference_dict[organism][source]
		gtf = reference_files['gtf']
		talon_abundance_file = reference_files['talon_abundance']
		sj_counts_file = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/{organism}-{source}-junction_counts.tsv'.format(**locals())
		infiles = [gtf, talon_abundance_file, sj_counts_file, outlier_sample_file]
		for comparison in comparisons+['all']:
			if comparison == 'all':
				gtf_prefix = os.path.basename(gtf)[:-len('.gtf')]
				comparison_string = '_vs_'.join(comparison) if comparison != 'all' else comparison
				outfile = 'arion/illumina/s04-alignment.dir/{organism}/{comparison_string}/gtf/{gtf_prefix}-{comparison_string}-SJ_filtered.gtf'.format(**locals())
				yield [infiles, outfile, comparison]

# @follows(getJunctionCounts)

@files(filterJobs)

def filterGTF(infiles, outfile, comparison):

	# Run
	run_r_job('filter_gtf', infiles, outfile, additional_params=comparison, run_locally=False, conda_env='env', W='00:10', GB=10, n=1, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

#############################################
########## 2. STAR index
#############################################

@transform(filterGTF,
		   regex(r'(.*.dir)/(.*?)/(.*)/gtf/.*.gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa'),
		   r'\1/\2/\3/STAR/index')

def buildStarIndexFiltered(infiles, outfile):

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=3, n=15, ow=True, print_cmd=False, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

# find arion/illumina/s04-alignment.dir/*/isoseq/*/STAR/index -name "job.log" | jsc

#############################################
########## 3. RSEM index
#############################################

@transform(filterGTF,
		   regex(r'(arion)/(.*.dir)/(.*?)/(.*)/gtf/(.*).gtf'),
		   add_inputs(r'\1/datasets/reference_genomes/\3/*.dna_sm.primary_assembly.fa'),
		   r'\1/\2/\3/\4/RSEM/index/\5.idx.fa')

def createRsemIndex(infiles, outfile):

	# Command
	basename = outfile[:-len('.idx.fa')]
	cmd_str = ''' rsem-prepare-reference --gtf {infiles[0]} --num-threads 10 {infiles[1]} {basename} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=False, stdout=basename+'.log', stderr=basename+'.err')

# find arion/illumina/s04-alignment.dir/*/*/*/RSEM/index -name "*.log" | jsc

#############################################
########## 4. STAR
#############################################

def starFilteredJobs():

	# Read samples
	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()

	# Loop through organisms
	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()

		# Loop through comparisons
		for comparison in comparison_dict[organism]+['all']:
			comparison_string = '_vs_'.join(comparison) if comparison != 'all' else comparison

			# Loop throush samples
			for sample_name, fastq_files in fastq_dict.items():
				cell_type = sample_name.split('_')[1].replace('2PN', '1C')

				# Select samples for each comparison
				# if cell_type==comparison[0] or cell_type==comparison[1] or comparison == 'all':
				if comparison == 'all' and sample_name not in outlier_samples[organism]:
					star_index = 'arion/illumina/s04-alignment.dir/{organism}/{comparison_string}/STAR/index'.format(**locals())
					sj_files = glob.glob('arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/pass1/*/*-SJ.out.tab'.format(**locals()))
					outfile = 'arion/illumina/s04-alignment.dir/{organism}/{comparison_string}/STAR/pass2/{sample_name}/{sample_name}-Aligned.sortedByCoord.out.bam'.format(**locals())
					yield [(fastq_files, star_index, sj_files), outfile]

@follows(getStarJunctions, buildStarIndexFiltered)

@files(starFilteredJobs)

def runStarFiltered(infiles, outfile):

	# Split
	fastq_files, star_index, sj_files = infiles

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# find arion/illumina/s04-alignment.dir/human/*/STAR/pass2 -name "*.log" | jsc

#############################################
########## 5. RSEM expression
#############################################

@follows(runStarFiltered, createRsemIndex)

# @transform('arion/illumina/s04-alignment.dir/human/2C_vs_4C/STAR/pass2/human_4C_B3_9/human_4C_B3_9-Aligned.toTranscriptome.out.bam',
@transform('arion/illumina/s04-alignment.dir/*/*/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
		   add_inputs(r'\1/RSEM/index/*.idx.fa'),
		   r'\1/RSEM/results/\2/\2.isoforms.results')

def runRsem(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Command
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness none \
		--paired-end \
		--estimate-rspd \
		--num-threads 200 \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log && \
		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())
		# --calc-ci \

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=2, n=25, modules=['rsem/1.3.3'], print_cmd=True, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

#############################################
########## 6. Create BigWig
#############################################

@transform(runStarFiltered,
		   suffix('-Aligned.sortedByCoord.out.bam'),
		   '.bw')

def createBigWig(infile, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

#############################################
########## 7. Filtered FASTA
#############################################

# @transform('arion/illumina/s04-alignment.dir/*/all/gtf/*SJ_filtered.gtf',
# 		   regex(r'(.*.dir)/(.*?)/(.*)/(.*).gtf'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*dna_sm.primary_assembly.fa'),
# 		   r'\1/\2/\3/SQANTI3/\4.fasta')

# def getFilteredFasta(infiles, outfile):

# 	# Command
# 	cmd_str = ''' gffread {infiles[0]} -g {infiles[1]} -w {outfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['gff/2021-02'], W='00:30', GB=10, n=1, stdout=outfile.replace('.fasta', '_fasta.log'), stderr=outfile.replace('.fasta', '_fasta.err'))

#############################################
########## 8. Run SQANTI
#############################################

@transform('arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf',
		   regex(r'(.*.dir)/(.*?)/(.*)/(.*).gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.gtf', r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa', '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt', [r'arion/illumina/s04-alignment.dir/\2/all/STAR/pass2/*/*-SJ.out.tab']),
		   r'\1/\2/\3/SQANTI3_test5/\4-SQANTI3_report.pdf')

def runSqanti(infiles, outfile):

	# Prepare
	sj_files = ','.join(infiles[4])
	dirname = os.path.dirname(outfile)
	basename = os.path.basename(outfile)[:-len('_report.pdf')]

	# Command
			# --polyA_motif_list {infiles[3]} \
			# --cpus 10 \
	cmd_str = ''' export PYTHONPATH=$PYTHONPATH:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/cDNA_Cupcake/ && export PYTHONPATH=$PYTHONPATH:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/cDNA_Cupcake/sequence/ && \
		python /sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/sqanti3_qc.py {infiles[0]} {infiles[1]} {infiles[2]} \
			--skipORF \
			-c {sj_files} \
			--report both \
			-d {dirname} \
			-o {basename}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='SQANTI3_v4.1', W='02:00', n=3, GB=10, print_cmd=True)

# find arion/illumina/s04-alignment.dir/*/all/gtf/SQANTI3 -name "SQANTI3" | xargs rm -r

#######################################################
#######################################################
########## S5. Expression
#######################################################
#######################################################

#############################################
########## 1. Prepare metadata
#############################################

# @collate('arion/illumina/s04-alignment.dir/*/*/RSEM/results/*/*.isoforms.results',
@collate(runRsem,
		 regex(r'(.*)/s04-alignment.dir/(.*)/(.*)/RSEM/.*/.*.isoforms.results'),
		 r'\1/s05-expression.dir/\2/\3/\2_\3-sample_metadata.txt')

def prepareSampleMetadata(infiles, outfile):
	
	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Prepare dict
	sample_dict = [{'files': x, 'names': x.split('/')[-2]} for x in infiles]

	# Convert to dataframe
	sample_dataframe = pd.DataFrame(sample_dict)

	# Add more columns
	organism = outfile.split('/')[-3]
	sample_dataframe['organism'] = [organism for x in sample_dataframe['names']]
	sample_dataframe['cell_type'] = [x.split('_')[1].replace('2PN', '1C') for x in sample_dataframe['names']]

	# Read metadata
	metadata_dataframe = pd.read_csv(organism_metadata[organism]).query('Platform == "ILLUMINA"').rename(columns={'sample_name': 'names'}).drop(['Platform', 'count'], axis=1)

	# Fix organism specific columns
	if organism == 'human':
		metadata_dataframe = metadata_dataframe.drop(['Group', 'batch'], axis=1)
		sample_dataframe['batch'] = [x.split('_')[2] for x in sample_dataframe['names']]
	elif organism == 'mouse':
		metadata_dataframe['strain'] = [x.replace(' ', '').replace('/', '_').replace('x', '_') for x in metadata_dataframe['Strain']]
		metadata_dataframe = metadata_dataframe.drop(['Run', 'Cell_type', 'sample_id', 'Strain'], axis=1)
		sample_dataframe['replicate'] = [x.split('_')[-1] for x in sample_dataframe['names']]

	# Merge
	merged_dataframe = sample_dataframe.merge(metadata_dataframe, on='names')

	# Write
	merged_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Aggregate
#############################################

@follows(filterGTF)

# @transform('arion/illumina/s05-expression.dir/*/*/*-sample_metadata.txt',
@transform(prepareSampleMetadata,
		   regex(r'(.*)/(s05-expression.dir)/(.*)/(.*)/.*-sample_metadata.txt'),
		   add_inputs(r'\1/s04-alignment.dir/\3/\4/gtf/*.gtf'),
		   r'\1/\2/\3/\4/\3_\4-counts.rda')

def aggregateCounts(infiles, outfile):

	# Run
	run_r_job('aggregate_counts', infiles, outfile, conda_env='env', run_locally=False, W='00:30', GB=30, n=1, wait=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# find arion/illumina/s05-expression.dir/*/all -name "*counts*" | xargs rm

#############################################
########## 3. Transcript TPM
#############################################

@transform(aggregateCounts,
		   suffix('counts.rda'),
		   'transcript_tpm.txt')

def getTranscriptTPM(infile, outfile):

	# Run
	run_r_job('get_transcript_tpm', infile, outfile, conda_env='env', run_locally=False, ow=False)

#############################################
########## 4. Gene counts
#############################################

# @transform('arion/illumina/s05-expression.dir/*/all/*-counts.rda',
@transform(aggregateCounts,
		   suffix('counts.rda'),
		   'gene_normalized_counts.tsv')

def getGeneExpression(infile, outfile):

	# Run
	run_r_job('get_gene_expression', infile, outfile, conda_env='env', run_locally=False, ow=False)

#############################################
########## 5. Get size factors
#############################################

@transform('arion/illumina/s05-expression.dir/*/all/*-counts.rda',
		   suffix('-counts.rda'),
		   '-size_factors.tsv')

def getSizeFactors(infile, outfile):

	# Run
	run_r_job('get_size_factors', infile, outfile, conda_env='env', modules=[], W='00:15', GB=20, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 6. Create scaled BigWig
#############################################

@transform('arion/illumina/s04-alignment.dir/human/all/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
		   regex(r'(.*.dir)/(.*?)/(.*)-Aligned.sortedByCoord.out.bam'),
		   add_inputs(r'arion/illumina/s05-expression.dir/\2/all/\2_all-size_factors.tsv'),
		   r'\1/\2/\3-scaled.bw')

def createScaledBigWig(infiles, outfile):

	# Read size factor
	normalization_dict = pd.read_table(infiles[1], index_col='sample_name')['size_factor_reciprocal'].to_dict()
	size_factor = normalization_dict[outfile.split('/')[-2]]

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=1 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_outfile=True)

#############################################
########## 7. Merge scaled BigWig
#############################################

@collate(createScaledBigWig,
		 regex(r'(.*)/s04-alignment.dir/(.*)/all/STAR/pass2/.*?_(.*?)_.*/.*-scaled.bw'),
		 add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chrom.sizes'),
		 r'\1/s05-expression.dir/\2/all/scaled_bw_test6/\2-\3.bw')

def mergeScaledBigWig(infiles, outfile):
	
	# Files
	wig_file = outfile.replace('.bw', '.wig')
	bedgraph_file = outfile.replace('.bw', '.bedgraph')
	infiles_str = ' '.join([x[0] for x in infiles])

	# Command
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} """.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_cmd=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

# wiggletools write_bg - arion/illumina/s05-expression.dir/human/all/scaled_bw/human-morula.wig | bedSort | head

# find arion/illumina/s04-alignment.dir/human/all/STAR/pass2 -name "*scaled.bw"

#######################################################
#######################################################
########## S5. Differential expression
#######################################################
#######################################################

#############################################
########## 1. DESeq2
#############################################

# @follows(aggregateCounts)

@subdivide('arion/illumina/s05-expression.dir/*/*/*-counts.rda',
		   regex(r'(.*)/s05-expression.dir/(.*)/(.*)/.*-counts.rda'),
		   add_inputs(comparison_file),
		   r'\1/s06-differential_expression.dir/\2/\3/\2_\3-*-deseq.tsv',
		   r'\1/s06-differential_expression.dir/\2/\3/\2_\3-{comparison[1]}_vs_{comparison[2]}-{feature}-deseq.tsv')

def runDESeq2(infiles, outfiles, outfileRoot):

	# Loop
	for feature in ['gene', 'transcript']:

		# Get info
		logname = os.path.join(os.path.dirname(outfileRoot), 'deseq-'+feature+'.log')
		jobname = '_'.join(logname[:-len('.log')].split('/')[-3:])
		
		# Run
		run_r_job('run_deseq2', infiles, outfileRoot, conda_env='env', W='24:00', GB=25, n=1, additional_params=feature, ow=False, stdout=logname, jobname=jobname)

# #############################################
# ########## 2. Excel
# #############################################

# # @follows(runDESeq2)

# @collate('arion/illumina/s06-differential_expression.dir/human/isoseq/human_isoseq-*-gene-deseq.tsv',
# 		 regex(r'(.*isoseq)-.*.tsv'),
# 		 r'\1-merged_v7.xlsx')

# def createExcel(infiles, outfile):

# 	# Read differential expression
# 	dataframes = {x.split('-')[-3]: pd.read_table(x).drop(['gene_source', 'lfcSE', 'stat'], axis=1) for x in infiles}

# 	# Read Dalit's genelist
# 	dalit_genelist = pd.read_csv('arion/datasets/dalit/Developmental_genelist_by_Dalit.csv').set_index('Gene_Symbol')['Classification'].to_dict()

# 	# Initialize ExcelWriter
# 	with pd.ExcelWriter(outfile) as writer:

# 		# Loop
# 		for comparison in ['1C_vs_2C', '2C_vs_4C', '4C_vs_8C', '8C_vs_morula', 'morula_vs_blastocyst']:

# 			# Get data
# 			dataframe = dataframes[comparison]

# 			# Split
# 			comparison_groups = comparison.split('_vs_')

# 			# Significance
# 			dataframe['gene_biotype'] = [rowData['gene_biotype'].replace('_', ' ') if rowData['gene_biotype'] != 'not_available' else rowData['gene_category'].replace('_', ' ') for index, rowData in dataframe.iterrows()]

# 			# Significance
# 			dataframe['significant'] = [rowData['padj'] < 0.05 and abs(rowData['log2FoldChange']) > 1 for index, rowData in dataframe.iterrows()]
# 			dataframe['differential_expression'] = ['Not significant' if not rowData['significant'] else 'Up in '+comparison_groups[1] if rowData['log2FoldChange'] > 0 else 'Down in '+comparison_groups[1] for index, rowData in dataframe.iterrows()]

# 			# Round
# 			dataframe['baseMean'] = dataframe['baseMean'].round(1)
# 			dataframe['log2FoldChange'] = dataframe['log2FoldChange'].round(1)
# 			dataframe['pvalue'] = ['{:.2e}'.format(x) for x in dataframe['pvalue']]
# 			dataframe['padj'] = ['{:.2e}'.format(x) for x in dataframe['padj']]

# 			# Dalit geneset
# 			dataframe['Developmental category'] = [dalit_genelist.get(x, '').replace('_', ' ').title() for x in dataframe['gene_name']]
			
# 			# Rename
# 			dataframe = dataframe.rename(columns={'gene_id': 'Gene ID', 'gene_name': 'Gene name', 'gene_biotype': 'Gene biotype', 'baseMean': 'Average expression', 'pvalue': 'P-value', 'padj': 'Adjusted P-value', 'differential_expression': 'Differentially expressed'}).drop(['significant', 'gene_category'], axis=1)

# 			# Write
# 			dataframe.to_excel(writer, sheet_name=comparison, index=False)
			
# 			# Modify
# 			workbook = writer.book
# 			worksheet = writer.sheets[comparison]
# 			format1 = workbook.add_format({'align': 'center'})
# 			worksheet.set_column('A:J', 20, format1)
# 			max_logfc = dataframe['log2FoldChange'].abs().max()
# 			worksheet.conditional_format('E1:E'+str(len(dataframe.index)+1), {'type': '3_color_scale', 'min_value': -max_logfc, 'min_type': 'num', 'mid_value': 0, 'mid_type': 'num', 'max_value': max_logfc, 'max_type': 'num', 'min_color': '#67a9cf', 'mid_color': '#ffffff', 'max_color': '#ef8a62'})
# 			worksheet.autofilter(0, 0, len(dataframe.index), len(dataframe.columns) - 1)

#######################################################
#######################################################
########## S7. Enrichment
#######################################################
#######################################################

#############################################
########## 1. GO
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s06-differential_expression.dir/*/all/*-gene-deseq.tsv',
		   regex(r'(.*)/s06-differential_expression.dir/(.*)/all/(.*)-gene-deseq.tsv'),
		   r'\1/s07-enrichment.dir/\2/go/\3-go_enrichment.tsv')

def runGoEnrichment(infile, outfile):

	# Run
	run_r_job('run_go_enrichment', infile, outfile, conda_env='env', W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'))

#############################################
########## 2. Domain
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s06-differential_expression.dir/*/all/*-transcript-deseq.tsv',
		   regex(r'(.*)/s06-differential_expression.dir/(.*)/all/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s07-pfam.dir/\2/\2-translated_pfam.tsv'),
		   r'\1/s07-enrichment.dir/\2/domain/\3-domain_enrichment.tsv')

def runDomainEnrichment(infiles, outfile):

	# Run
	run_r_job('run_domain_enrichment', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'))

# #############################################
# ########## 3. Repeats
# #############################################

# @follows(runDESeq2)

@transform('arion/illumina/s06-differential_expression.dir/*/all/*-transcript-deseq.tsv',
		   regex(r'(.*)/s06-differential_expression.dir/(.*)/all/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s08-repeatmasker.dir/\2/*_repeatmasker.tsv'),
		   r'\1/s07-enrichment.dir/\2/repeat/\3-repeat_enrichment.tsv')

def runRepeatEnrichment(infiles, outfile):

	# Run
	run_r_job('run_repeat_enrichment', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'))

#######################################################
#######################################################
########## S8. WGCNA
#######################################################
#######################################################

#############################################
########## 1. Pick soft thresholds
#############################################

# @follows(getGeneExpression)

# @transform('arion/illumina/s05-expression.dir/human/all/human_all-gene_normalized_counts.tsv',
@transform('arion/illumina/s05-expression.dir/human/all/human_all-counts.rda',
		   regex(r'(.*)/s05-expression.dir/(.*)/all/.*.rda'),
		   r'\1/s08-wgcna.dir/\2/network/\2-soft_thresholds_signed.rda')

def pickSoftThresholds(infile, outfile):

	# Run
	run_r_job('pick_soft_thresholds', infile, outfile, modules=['R/4.0.3'], W='00:45', GB=25, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 2. Cluster genes
#############################################

# @follows(getGeneExpression)

@transform(pickSoftThresholds,
		   regex(r'(.*)-soft_thresholds_(.*).rda'),
		   r'\1-gene_network_\2.rda')

def clusterGenes(infile, outfile):

	# Run
	run_r_job('cluster_genes', infile, outfile, modules=['R/4.0.3'], W='00:45', GB=50, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Get modules
#############################################

# @follows(getGeneExpression)

@transform(clusterGenes,
		   suffix('.rda'),
		   add_inputs(pickSoftThresholds),
		   '_modules.rda')

def getGeneModules(infiles, outfile):

	# Run
	run_r_job('get_gene_modules', infiles, outfile, modules=['R/4.0.3'], W='00:45', GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 4. Get enrichment
#############################################

# @follows(getGeneExpression)

@transform(getGeneModules,
		   suffix('s.rda'),
		   '_enrichment.rda')

def runModuleEnrichment(infiles, outfile):

	# Run
	run_r_job('run_module_enrichment', infiles, outfile, run_locally=True)# conda env, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 5. Get module preservation
#############################################

# @follows(getGeneExpression)

@transform(('arion/geo_illumina/s04-expression.dir/*/*-counts.rda', 'arion/geo_illumina/s05-primates.dir/*/RSEM/counts/*-counts.rda'),
		   regex(r'.*/(.*)-counts.rda'),
		   add_inputs(pickSoftThresholds, getGeneModules),
		   r'arion/illumina/s08-wgcna.dir/human/module_preservation/\1-module_preservation.rda')

def getModulePreservation(infiles, outfile):

	# Run
	run_r_job('get_module_preservation', infiles, outfile, W='02:00', GB=10, n=5, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 6. Get gene connectivity
#############################################

# @follows(getGeneExpression)

@transform(getGeneModules,
		   suffix('s.rda'),
		   add_inputs(pickSoftThresholds, 'arion/illumina/s08-wgcna.dir/human/network/human-gene_network_signed_adjacency.rda'),
		   '_membership.rda')

def getModuleMembership(infiles, outfile):

	# Run
	run_r_job('get_module_membership', infiles, outfile, W='00:10', GB=30, n=1, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 7. Get gene networks
#############################################

# @follows(getGeneExpression)

@transform(getModuleMembership,
		   regex(r'(.*)/.*/.*.rda'),
		   add_inputs(pickSoftThresholds, getGeneModules, reference_dict['human']['isoseq']['gtf_filtered'], 'arion/datasets/dalit/Developmental_genelist_by_Dalit_NJFedit.csv'),
		   r'\1/gene_networks_top10/network-nodes.tsv')

def getGeneNetworks(infiles, outfile):

	# Run
	run_r_job('get_gene_networks', infiles, outfile, run_locally=True, W='00:10', GB=30, n=1, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# @transform(getGeneModules,
# 		   suffix('s.rda'),
# 		   add_inputs('arion/illumina/s05-expression.dir/human/all/human_all-gene_normalized_counts.tsv'),
# 		   '_correlation.tsv')

# def getModuleCorrelations(infiles, outfile):

# 	# Run
# 	run_r_job('get_module_correlations', infiles, outfile, modules=['R/4.0.3'], stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S9. SUPPA
#######################################################
#######################################################

#############################################
########## 1. Index
#############################################

def suppaIndexJobs():
	for organism, organism_references in reference_dict.items():
		for file_format in ['ioi', 'ioe']:
			infile = organism_references['isoseq']['gtf_filtered']
			outdir = 'arion/illumina/s09-suppa.dir/{organism}/01-events/{file_format}/'.format(**locals())
			yield [infile, outdir, file_format]

# @follows(functionToFollow)

@files(suppaIndexJobs)

def buildSuppaIndex(infile, outdir, file_format):

	# Basename
	basename = outdir.split('/')[-4]

	# Command
	if file_format == 'ioe':
		cmd_str = '''python $SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outdir}{basename} -f {file_format} -e SE SS MX RI FL'''.format(**locals())
	elif file_format == 'ioi':
		cmd_str = '''python $SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outdir}{basename} -f {file_format}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, W='00:15', modules=['suppa/2.3'], GB=10, n=1, stdout=os.path.join(outdir, 'job.log'), jobname='_'.join(outdir.split('/')[-4:]).replace('_01-events', ''), ow=True)

#############################################
########## 2. PSI
#############################################

def psiJobs():
	for organism in ['human', 'mouse']:
		tpm_file = 'arion/illumina/s05-expression.dir/{organism}/all/{organism}_all-transcript_tpm.txt'.format(**locals())
		indices = glob.glob('arion/illumina/s09-suppa.dir/{organism}/01-events/io*/*.io?'.format(**locals()))
		for suppa_index in indices:
			file_format = suppa_index.split('.')[-1]
			if file_format == 'ioe':
				infiles = [tpm_file, suppa_index]
				event_type = suppa_index.split('_')[-2]
			elif file_format == 'ioi':
				infiles = [tpm_file, reference_dict[organism]['isoseq']['gtf_filtered']]
				event_type = 'isoform'
			outfile = 'arion/illumina/s09-suppa.dir/{organism}/02-psi/{organism}_{event_type}.psi'.format(**locals())
			yield [infiles, outfile, event_type] 

@follows(buildSuppaIndex)

@files(psiJobs)

def getSuppaPSI(infiles, outfile, event_type):
	
	# Isoform
	if event_type == 'isoform':
		outname = outfile[:-len('_isoform.psi')]
		cmd_str = 'python $SUPPA_HOME/suppa.py psiPerIsoform -g {infiles[1]} -e {infiles[0]} -o {outname}'.format(**locals())
	# Event
	else:
		outname = outfile[:-len('.psi')]
		cmd_str = 'python $SUPPA_HOME/suppa.py psiPerEvent --ioe-file {infiles[1]} --expression-file {infiles[0]} -o {outname}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="00:30", GB=10, n=1, modules=['suppa/2.3'], ow=False, stderr=outfile.replace('.psi', '.err'))

#############################################
########## 3. Split
#############################################

def splitJobs():
	infiles = glob.glob('arion/illumina/s05-expression.dir/*/*/*-transcript_tpm.txt')+glob.glob('arion/illumina/s09-suppa.dir/*/02-psi/*.psi')
	for infile in infiles:
		infile_split = infile.split('/')
		organism = infile_split[3]
		basename = infile_split[-1].replace('-', '_').replace('.psi', '').replace('.txt', '').replace('all_', '')
		outdir = 'arion/illumina/s09-suppa.dir/{organism}/03-split/{basename}'.format(**locals())
		yield [infile, outdir, basename]

# @follows(getTranscriptTPM, getSuppaPSI)

@files(splitJobs)

def splitData(infile, outdir, basename):

	# Get metadata
	dataframe = pd.read_table(infile)

	# Get groups
	group_dataframe = pd.DataFrame(dataframe.columns).rename(columns={0: 'sample'})
	group_dataframe['cell_type'] = [x.split('_')[1].replace('2PN', '1C') for x in group_dataframe['sample']]
	group_dict = group_dataframe.groupby('cell_type')['sample'].apply(lambda x: list(x))

	# Loop
	for cell_type, samples in group_dict.items():

		# Get outfile
		file_format = 'tsv' if 'tpm' in infile else 'psi'
		outfile = '{outdir}/{basename}_{cell_type}.{file_format}'.format(**locals())

		# Outdir
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		
		# Subset and write
		dataframe[samples].to_csv(outfile, sep='\t', index_label=False, na_rep='nan')

#############################################
########## 4. Differential splicing
#############################################

def suppaRunJobs():
	for organism, comparisons in comparison_dict.items():
		for event_type in ['isoform', 'A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']:
			if event_type == 'isoform':
				suppa_index = 'arion/illumina/s09-suppa.dir/{organism}/01-events/ioi/{organism}.ioi'.format(**locals())
			else:
				suppa_index = 'arion/illumina/s09-suppa.dir/{organism}/01-events/ioe/{organism}_{event_type}_strict.ioe'.format(**locals())
			for comparison in comparisons:
				infiles = []
				for cell_type in comparison:
					infiles.append('arion/illumina/s09-suppa.dir/{organism}/03-split/{organism}_{event_type}/{organism}_{event_type}_{cell_type}.psi'.format(**locals()))
					infiles.append('arion/illumina/s09-suppa.dir/{organism}/03-split/{organism}_transcript_tpm/{organism}_transcript_tpm_{cell_type}.tsv'.format(**locals()))
				infiles.append(suppa_index)
				outfile = 'arion/illumina/s09-suppa.dir/{organism}/04-dpsi/{organism}-{comparison[0]}_vs_{comparison[1]}/{organism}-{comparison[0]}_vs_{comparison[1]}-{event_type}/{organism}-{comparison[0]}_vs_{comparison[1]}-{event_type}.psivec'.format(**locals())
				yield [infiles, outfile]

# @follows(splitData)

@files(suppaRunJobs)

def getDifferentialPSI(infiles, outfile):

	# Full paths
	infiles = [os.path.join(os.getcwd(), x) for x in infiles]

	# Command
	outname = os.path.join(os.getcwd(), outfile.rsplit('.', 1)[0])
	
	# Command
	cmd_str = 'python $SUPPA_HOME/suppa.py diffSplice -m empirical --save_tpm_events --input {infiles[4]} --psi {infiles[0]} {infiles[2]} --tpm {infiles[1]} {infiles[3]} -gc -o {outname}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="03:00", GB=50, n=1, modules=['suppa/2.3', 'python/3.8.2'], stdout=outfile.replace('.psivec', '.log'), stderr=outfile.replace('.psivec', '.err'))

# find arion/illumina/s09-suppa.dir/*/04-dpsi -name "*.log" | jsc
# find arion/illumina/s09-suppa.dir/*/04-dpsi -name "*.log" | xargs wc -l
# find arion/illumina/s09-suppa.dir/*/04-dpsi -name "*.err" | xargs wc -l

# # ls /hpc/users/torred23/pipelines/projects/early-embryo/arion/illumina/s09-suppa.dir/*/*/04-dpsi/*/*/*.log | js | grep -v completed
# find arion/illumina/s09-suppa.dir/mouse/ensembl/04-dpsi -name "*.psivec" | wc -l
# # ls /hpc/users/torred23/pipelines/projects/early-embryo/arion/illumina/s09-suppa.dir/*/*/04-dpsi/*/*/*.err | lr | grep Error > error_samples.txt

#############################################
########## 5. Summary
#############################################

# getDifferentialPSI
@collate('arion/illumina/s09-suppa.dir/*/04-dpsi/*/*/*.psivec',
		 regex(r'(.*)/04-dpsi/(.*)/.*/.*.psivec'),
		 r'\1/05-summaries/\2-suppa_summary.tsv')

def createSuppaSummary(infiles, outfile):

	# Make summary
	summary_dataframe = pd.concat([S.createSuppaSummary(x[:-len('.psivec')]) for x in infiles if 'isoform' not in x])

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	summary_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 6. Summary
#############################################

# getDifferentialPSI
@transform('arion/illumina/s09-suppa.dir/*/04-dpsi/*/*/*isoform.psivec',
		    regex(r'(.*)/04-dpsi/(.*)/.*/.*.psivec'),
		    r'\1/05-summaries/\2-suppa_summary_isoform.tsv')

def createSuppaIsoformSummary(infile, outfile):

	# Make summary
	summary_dataframe = S.createSuppaSummary(infile[:-len('.psivec')], event_type='isoform')

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	summary_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 6. Cluster PSI
#############################################

@collate(createSuppaIsoformSummary,
		 regex(r'(.*)/05-summaries/(human)-.*.tsv'),
		 add_inputs(r'\1/02-psi/\2_isoform.psi'),
		 r'\1/06-psi_clusters/\2_isoform-timepoints.rda')

def clusterPSI(infiles, outfile):

	# Split
	summary_files = [x[0] for x in infiles]
	psi_file = infiles[0][1]
	infiles = [psi_file]+summary_files

	# Run
	run_r_job('cluster_psi', infiles, outfile, modules=['R/4.0.3'], W="02:00", GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 7. Get clusters
#############################################

@transform(clusterPSI,
		   suffix('timepoints.rda'),
		   'clusters.rda')

def getPSIClusters(infile, outfile):

	# Run
	run_r_job('get_psi_clusters', infile, outfile, modules=['R/4.0.3'], W="01:00", GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#######################################################
#######################################################
########## S10. Isoform switching
#######################################################
#######################################################

#############################################
########## 1. Filter
#############################################

def isoformFilterJobs():
	for organism in ['human']: #, 'mouse'
		filtered_gtf = glob.glob('arion/illumina/s04-alignment.dir/{organism}/all/gtf/*102_talon-all-SJ_filtered.gtf'.format(**locals()))[0]
		for file_type in ['gtf_cds', 'cpat_predictions', 'pfam_predictions', 'transcript_fasta']:
			infile = reference_dict[organism]['isoseq'][file_type]
			infile_basename, infile_extension = os.path.basename(infile).rsplit('.', 1)
			infiles = [infile, filtered_gtf]
			outfile = 'arion/illumina/s10-isoform_switching.dir/{organism}/filtered_data/{infile_basename}_filtered.{infile_extension}'.format(**locals())
			yield [infiles, outfile, file_type]

@files(isoformFilterJobs)

def filterIsoformData(infiles, outfile, file_type):

	# Run
	run_r_job('filter_isoform_data', infiles, outfile, additional_params=file_type, conda_env='env')

#############################################
########## 2. Load data
#############################################

@collate('arion/illumina/s10-isoform_switching.dir/human/filtered_data/*',
		 regex(r'(.*)/(.*)/filtered_data/.*'),
		 add_inputs(r'arion/illumina/s05-expression.dir/\2/all/\2_all-sample_metadata.txt', comparison_file),
		 r'\1/\2/\2-isoforms.rda')

def loadIsoformData(infiles, outfile):

	# Split
	infiles = [x[0] for x in infiles]+[infiles[0][1], infiles[0][2]]

	# Run
	run_r_job('load_isoform_data', infiles, outfile, conda_env='env', W="06:00", GB=30, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Run
#############################################

@transform(loadIsoformData,
		  suffix('.rda'),
		  '_results.rda')

def getIsoformSwitching(infile, outfile):

	# Run
	run_r_job('get_isoform_switching', infile, outfile, conda_env='env', W="06:00", GB=30, n=1, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#######################################################
#######################################################
########## S6. rMATS
#######################################################
#######################################################

#############################################
########## 1. Group bams
#############################################

# #runStar
# @collate('arion/illumina/s04-alignment.dir/*/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
# 		 regex(r'(.*)/s04-alignment.dir/(.*)/.*/STAR/pass2/.*/.*?_(.*?)_.*.bam'),
# 		 r'\1/s09-rmats.dir/\2/bams/\2-\3-bams.txt')

# def groupBams(infiles, outfile):

# 	# Create directory
# 	outdir = os.path.dirname(outfile)
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)

# 	# # Write
# 	# with open(outfile, 'w') as openfile: #.replace('2PN', '1C')
# 	# 	openfile.write(','.join(infiles))

# #############################################
# ########## 2. Run prep
# #############################################

# # @follows(groupBams)

# @transform('arion/illumina/s09-rmats.dir/*/bams/*.txt',
# 		  regex(r'(.*)/(s09-rmats.dir)/(.*)/bams/(.*).txt'),
# 		  add_inputs(r'\1/s04-alignment.dir/\3/all/gtf/*.gtf'),
# 		  r'\1/\2/\3/prep/\4')

# def runRmatsPrep(infiles, outfile):

# 	# Read length
# 	read_length = 125 if 'human' in outfile else 150

# 	# Command
# 	outdir = outfile
# 	cmd_str = '''python $RMATS_HOME/rmats.py --b1 {infiles[0]} --gtf {infiles[1]} --readLength {read_length} --variable-read-length --nthread 24 --od {outdir} --tmp {outdir}/tmp --task prep'''.format(**locals()) #--statoff 

# 	# Run
# 	# print(cmd_str)
# 	W = "01:00" if 'mouse' in outfile else "03:00"
# 	run_job(cmd_str, outfile, W=W, GB=5, n=6, modules=['rmats/4.1.0'], stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'), jobname=os.path.basename(outfile)+'_rmats_prep', ow=False, print_outfile=False)

# #############################################
# ########## 3. Run post
# #############################################

# def rmatsJobs():
# 	for organism, comparisons in comparison_dict.items():
# 		if organism == 'human':
# 			for source, reference_files in reference_dict[organism].items():
# 				for comparison in comparisons:
# 					infiles = []
# 					for cell_type in comparison:
# 						infiles.append('arion/illumina/s09-rmats.dir/{organism}/bams/{organism}-{cell_type}-bams.txt'.format(**locals()))
# 						infiles += glob.glob('arion/illumina/s09-rmats.dir/{organism}/prep/{organism}-{cell_type}-bams/tmp/*.rmats'.format(**locals()))
# 					infiles += glob.glob('arion/illumina/s04-alignment.dir/{organism}/all/gtf/*-all-SJ_filtered.gtf'.format(**locals()))
# 					outfile = 'arion/illumina/s09-rmats.dir/{organism}/post/{organism}-{comparison[0]}_vs_{comparison[1]}'.format(**locals())
# 					yield [infiles, outfile]

# # @follows(runStar, linkRmatsPrep)

# @files(rmatsJobs)

# def runRmatsPost(infiles, outfile):

# 	# Split
# 	bam_files = [infiles[0], infiles[2]]
# 	prep_files = [os.path.join(os.getcwd(), x) for x in [infiles[1], infiles[3]]]
# 	gtf_file = infiles[4]

# 	# Get temp dir
# 	tempdir = os.path.join(outfile, 'tmp')

# 	# Read length
# 	read_length = 125 if 'human' in outfile else 150

# 	# Command
# 	cmd_str = '''mkdir -p {tempdir} && ln -s {prep_files[0]} {tempdir} && ln -s {prep_files[1]} {tempdir} &&  \
# 		python $RMATS_HOME/rmats.py \
# 			--b1 {bam_files[0]} \
# 			--b2 {bam_files[1]} \
# 			--gtf {gtf_file} \
# 			--readLength {read_length} \
# 			--variable-read-length \
# 			--nthread 24 \
# 			--od {outfile} \
# 			--tmp {tempdir} \
# 			--task post '''.format(**locals()) #--statoff 

# 	# Run
# 	W = "01:00" if 'mouse' in outfile else "03:00"
# 	if 'morula' in outfile and 'blastocyst' in outfile:
# 		run_job(cmd_str, outfile, W=W, GB=5, n=6, modules=['rmats/4.1.0'], stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'), jobname=os.path.basename(outfile)+'_rmats_post', ow=False, print_cmd=False)

# #############################################
# ########## 2. Run
# #############################################

# def rmatsJobs():
# 	for organism, comparisons in comparison_dict.items():
# 		for source, reference_files in reference_dict[organism].items():
# 			for comparison in comparisons:
# 				infiles = []
# 				for cell_type in comparison:
# 					infiles.append('arion/illumina/s09-rmats.dir/{organism}/bams/{organism}-{cell_type}-bams.txt'.format(**locals()))
# 				infiles.append(reference_files['gtf'])
# 				outfile = 'arion/illumina/s09-rmats.dir/{organism}/results/{organism}-{comparison[0]}_vs_{comparison[1]}_stat'.format(**locals())
# 				yield [infiles, outfile]

# # @follows(runStar)

# @files(rmatsJobs)

# def runRmats(infiles, outfile):

# 	# Read length
# 	read_length = 125 if 'human' in outfile else 150

# 	# Command
# 	cmd_str = '''python $RMATS_HOME/rmats.py --b1 {infiles[0]} --b2 {infiles[1]} --gtf {infiles[2]} --readLength {read_length} --variable-read-length --nthread 24 --od {outfile} --tmp /tmp '''.format(**locals()) #--statoff 

# 	# Run
# 	# print(cmd_str)
# 	W = "01:00" if 'mouse' in outfile else "03:00"
# 	if 'morula' in outfile and 'blastocyst' in outfile:
# 		run_job(cmd_str, outfile, W=W, GB=5, n=6, modules=['rmats/4.1.0'], stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'), jobname=os.path.basename(outfile)+'_rmats', ow=False, print_outfile=False, wait=False)

# python $RMATS_HOME/rMATS_P/prepare_stat_inputs.py --new-output-dir arion/illumina/s09-rmats.dir/human/results/human-morula_vs_blastocyst/stats --old-output-dir arion/illumina/s09-rmats.dir/human/results/human-morula_vs_blastocyst --group-1-indices 0,1,2 --group-2-indices 3,4,5


# find arion/illumina/s09-rmats.dir -name "job.log" | js

#############################################
########## 3. Summary
#############################################

# @follows(runRmats)

# @collate('arion/illumina/s09-rmats.dir/*/*/results/*/*.MATS.JC.txt',
# 		 regex(r'(.*)/results/(.*)/.*.txt'),
# 		 r'\1/summaries/\2-rmats_summary.tsv')

# def createRmatsSummary(infiles, outfile):

# 	# Make summary
# 	summary_dataframe = pd.concat([S.createRmatsSummary(x) for x in infiles])

# 	# Outdir
# 	outdir = os.path.dirname(outfile)
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)

# 	# Write
# 	summary_dataframe.to_csv(outfile, sep='\t', index=False)

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