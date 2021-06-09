#################################################################
#################################################################
############### Embryo IsoSeq ################
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
import sqlite3
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/embryo-isoseq.R'
py_source = 'pipeline/scripts/EmbryoIsoseq.py'
P = 'acc_apollo'
# q = 'express'
# q = 'premium'
q = 'sla'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import EmbryoIsoseq as P

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
with open('/hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/reference_genomes/reference-genomes.json') as openfile:
	genome_indices = json.load(openfile)


##### 3. Variables #####
# Primers
clontech_primer_file = 'arion/isoseq/s01-isoseq.dir/primer.fasta'

# Qiao
qiao_metadata = 'arion/datasets/qiao/qiao-sample_names.csv'
qiao_subreads = 'arion/datasets/qiao/rawdata/pacbio/*/*.subreads.bam'
qiao_flnc = 'arion/isoseq/s01-isoseq.dir/mouse/flnc/*.flnc.bam'

# Human
human_metadata = 'arion/datasets/human/human-sample_names.csv'
human_flnc = 'arion/datasets/human/human_embryo_pacb/*.flnc.bam'

# Illumina
illumina_fastq = 'arion/illumina/s01-fastq.dir/*/trimmed/*/*.fq.gz'

# Pfam
# illumina_dict = {
# 	'human': glob.glob(human_illumina),
# 	'mouse': glob.glob(qiao_illumina)
# }

#######################################################
#######################################################
########## S1. Process data
#######################################################
#######################################################

#############################################
########## 1. Index
#############################################

@transform(qiao_subreads,
		   suffix('.bam'),
		   '.bam.pbi')

def indexSubreads(infile, outfile):

	# Command
	cmd_str = ''' pbindex {infile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 2. CCS
#############################################

@follows(indexSubreads)

@subdivide(qiao_subreads,
		   regex(r'(.*)/datasets/(.*)/rawdata/pacbio/(.*)/.*.subreads.bam'),
		   r'\1/isoseq/s01-isoseq.dir/mouse/ccs/\3/chunk_*/\3_chunk_*.ccs.bam',
		   r'\1/isoseq/s01-isoseq.dir/mouse/ccs/\3/chunk_{i}/\3_chunk_{i}.ccs.bam')

def runCCS(infile, outfiles, outfileRoot):

	# Loop
	n = 50
	for i in range(1, n+1):

		# Get outfile
		outfile = outfileRoot.format(**locals())

		# Command
		cmd_str = 'ccs {infile} {outfile} --min-rq 0.9 --chunk {i}/{n}'.format(**locals())

		# Run
		run_job(cmd_str, outfile, conda_env='isoseq3', W='01:00', GB=4, n=6, print_outfile=False)

#############################################
########## 3. lima
#############################################
# primer.fasta file from https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md#step-2---primer-removal-and-demultiplexing on 2020/11/05

@transform(runCCS,
		   suffix('ccs.bam'),
		   add_inputs(clontech_primer_file),
		   'demux.Clontech_5p--NEB_Clontech_3p.bam')

def runLima(infiles, outfile):

	# Get basename
	basename = outfile.replace('.Clontech_5p--NEB_Clontech_3p', '')

	# Command
	cmd_str = 'lima --isoseq --different --min-passes 1 --split-named --dump-clips --dump-removed --peek-guess {infiles[0]} {infiles[1]} {basename}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_cmd=False)

#############################################
########## 4. Refine
#############################################

@transform(runLima,
		   suffix('.bam'),
		   add_inputs(clontech_primer_file),
		   '.flnc.bam')

def refineReads(infiles, outfile):

	# Command
	cmd_str = 'isoseq3 refine --require-polya {infiles[0]} {infiles[1]} {outfile}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_cmd=False)

#############################################
########## 5. Merge reads
#############################################

@collate(refineReads,
		 regex(r'(.*)/ccs/(.*)/chunk_.*/.*.flnc.bam'),
		 r'\1/flnc/\2.flnc.bam')

def mergeChunks(infiles, outfile):

	# Command
	infiles_str = ' '.join(infiles)
	cmd_str = 'pbmerge -o {outfile} {infiles_str}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:30', GB=15, n=1, print_outfile=False)

#############################################
########## 6. Merge reports
#############################################

# @follows(refineReads)

# @collate('arion/isoseq/s01-isoseq.dir/*/ccs/*/*/*.demux.Clontech_5p--NEB_Clontech_3p.flnc.report.csv',
# 		 regex(r'(.*)/ccs/(.*)/.*/.*.report.csv'),
# 		 r'\1/flnc/\2.flnc.report.csv')

# def mergeReports(infiles, outfile):

# 	# Read
# 	print('Doing {}...'.format(outfile))
# 	dataframe = pd.concat([pd.read_csv(x) for x in infiles])

# 	# Write
# 	dataframe.to_csv(outfile, index=False)

#############################################
########## 7. Convert to FASTA
#############################################

# @transform('arion/isoseq/s01-isoseq.dir/mouse/flnc/SRR10267008.flnc.bam',
@transform((mergeChunks, human_flnc),
		   suffix('.bam'),
		   '.fastq.gz')

def flncToFastq(infile, outfile):

	# Command
	cmd_str = ''' samtools view {infile} | awk '{{printf("@%s\\n%s\\n+\\n%s\\n", $1, $10, $11)}}' | gzip > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.11'], W='01:00', GB=10, n=1, print_outfile=False, print_cmd=False)

#######################################################
#######################################################
########## S2. Align
#######################################################
#######################################################

#############################################
########## 1. Copy FASTQ
#############################################

def fastqJobs():

	# Paths
	paths = {
		'mouse': 'arion/isoseq/s01-isoseq.dir/mouse/flnc/{sample_id}.flnc.fastq.gz',
		'human': 'arion/datasets/human/human_embryo_pacb/{sample_id}.flnc.fastq.gz'
	}

	# Get sample names
	sample_dict = {x.split('/')[-2].replace('qiao', 'mouse'): pd.read_csv(x).query('Platform == "PACBIO_SMRT"').set_index('sample_id')['sample_name'].to_dict() for x in [qiao_metadata, human_metadata]}

	# Loop
	for organism, samples in sample_dict.items():
		for sample_id, sample_name in samples.items():
			infile = os.path.join(os.getcwd(), paths[organism].format(**locals()))
			outfile = 'arion/isoseq/s02-alignment.dir/fastq/{organism}/{sample_name}.flnc.fastq'.format(**locals())
			yield [infile, outfile]

@follows(flncToFastq)

@files(fastqJobs)

def copyFASTQ(infile, outfile):

	# Run
	run_r_job('copy_fastq', infile, outfile, run_locally=True)#conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. minimap2
#############################################

@transform(copyFASTQ,
		   regex(r'(.*)/fastq/(.*)/(.*).flnc.fastq'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.primary_assembly.fa'),
		   r'\1/minimap2/\2/\3/\3-minimap2.sam')

def runMinimap2(infiles, outfile):

	# Command
	outname = outfile.rsplit('.', 1)[0]
	cmd_str = ''' minimap2 -ax splice -uf --secondary=no -C5 -t 30 --MD {infiles[1]} {infiles[0]} | samtools sort -O sam -o {outfile} && \
		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f1,5 > {outname}.mapq '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.9', 'minimap2/2.17'], W='03:00', GB=6, n=6, print_outfile=False, stderr=outfile.replace('.sam', '.log'))

#############################################
########## 3. Filter
#############################################

@transform(runMinimap2,
		   regex(r'(.*).sam'),
		   r'\1.filtered.sam')

def filterSam(infile, outfile):

	# Files
	outname = outfile.replace('.sam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --sam-input --format sam --filter "not unmapped and mapping_quality >= 50" {infile} | samtools sort -O sam -o {outfile} && \
		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f5 -d "	" > {outname}.mapq'''.format(**locals())
	# and not duplicate

	# Run
	run_job(cmd_str, outfile, W='00:15', n=5, GB=1, modules=['samtools/1.9', 'sambamba/0.5.6'], print_cmd=False)

#############################################
########## 4. MultiQC
#############################################

@transform('arion/isoseq/s02-alignment.dir',
		   regex(r'(.*)/s..-(.*).dir'),
		   r'\1/multiqc/\2/multiqc_report.html')

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="00:20", GB=5, n=1, print_outfile=False, run_locally=False)

#######################################################
#######################################################
########## S3. Illumina alignment
#######################################################
#######################################################

#############################################
########## 1. STAR index
#############################################

# @follows(filterSam)

@transform('arion/datasets/reference_genomes/*/*.102.gtf',
		   regex(r'(.*)/(.*)/(.*).102.gtf'),
		   add_inputs(r'\1/\2/\3.dna_sm.primary_assembly.fa'),
		   r'arion/isoseq/s03-illumina_alignment.dir/\2/STAR/index/')

def buildStarIndex(infiles, outfile):

	# Split
	reference_gtf, reference_fasta = infiles

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {reference_fasta} --sjdbGTFfile {reference_gtf} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='01:00', GB=5, n=15, ow=True, print_cmd=False, jobname='_'.join(outfile.split('/')[-4:-1]))

#############################################
########## 2. STAR pass 1
#############################################

@follows(buildStarIndex)

@collate(illumina_fastq,
		regex(r'.*/s01-fastq.dir/(.*)/trimmed/(.*)/.*.fq.gz'),
		add_inputs(r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/index/'),
		r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass1/\2/\2-SJ.out.tab')

def getSJsPass1(infiles, outfile):

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[0][1]} \
		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 50 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=5, n=10, modules=['star/2.7.5b'], q='express')

# ls arion/isoseq/s03-illumina_alignment.dir/*/STAR/pass1/*/*-Log.progress.out | lr

#############################################
########## 3. STAR pass 2
#############################################

@follows(getSJsPass1)

@collate(illumina_fastq,
		regex(r'.*/s01-fastq.dir/(.*)/trimmed/(.*)/.*.fq.gz'),
		add_inputs(r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/index/', r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass1/*/*-SJ.out.tab'),
		r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass2/\2/\2-SJ.out.tab')

def getSJsPass2(infiles, outfile):

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]
	sj_files_str = ' '.join(infiles[0][2:])

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[0][1]} \
		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 50 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="10:00", GB=10, n=10, modules=['star/2.7.5b'], print_cmd=False, stdout=outfile.replace('-SJ.out.tab', '-job.log'), q='premium', wait=False)

# find arion/isoseq/s03-illumina_alignment.dir/*/STAR/pass2 -name "*job.log" | jsc

#############################################
########## 4. Merge SJ files
#############################################

@collate((getSJsPass1, getSJsPass2),
		 regex(r'(.*.dir)/(.*)/STAR/.*.tab'),
		 r'\1/\2/STAR/\2-merged.SJ.out.tab')

def mergeSJs(infiles, outfile):

	# Run
	run_r_job('merge_sjs', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1)

#######################################################
#######################################################
########## S4. TranscriptClean
#######################################################
#######################################################

#############################################
########## 1. Trim genome FASTA
#############################################

@transform('arion/datasets/reference_genomes/*/*.dna_sm.primary_assembly.fa',
		   suffix('.fa'),
		   '_renamed.fa')

def renameFastaChromosomes(infile, outfile):

	# Command
	cmd_str = ''' cut -d ' ' -f1 {infile} > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, W="00:30", GB=15, n=1)

#############################################
########## 2. Clean
#############################################

@follows(mergeSJs, renameFastaChromosomes)

@transform(filterSam,
		   regex(r'(.*)/s02-alignment.dir/minimap2/(.*)/(.*)/.*.sam'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly_renamed.fa', r'arion/isoseq/s03-illumina_alignment.dir/\2/STAR/\2-merged.SJ.out.tab'),
		   r'\1/s04-cleaned_transcripts.dir/\2/\3/\3_clean.sam')

def runTranscriptClean(infiles, outfile):

	# Prefix
	tempdir = os.path.join(os.path.dirname(outfile), 'tmp')
	prefix = outfile[:-len('_clean.sam')]

	# Command
	cmd_str = ''' python /sc/arion/work/torred23/libraries/TranscriptClean/TranscriptClean-2.0.2/TranscriptClean.py \
		--sam {infiles[0]} \
		--genome {infiles[1]} \
		--spliceJns {infiles[2]} \
		--tmpDir {tempdir} \
		--threads 10 \
		--deleteTmp \
		--canonOnly \
		--primaryOnly \
		--outprefix {prefix}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, W="03:00", GB=5, n=10, conda_env='talon', q='sla', stdout=outfile.replace('_clean.sam', '_job.log'), stderr=outfile.replace('_clean.sam', '_job.err'))

# find arion/isoseq/s04-cleaned_transcripts.dir -name "*_job.log" | jsc
# find arion/isoseq/s04-cleaned_transcripts.dir -name "*_job.err" | xargs wc -l

#############################################
########## 3. Flag reads
#############################################

@transform(runTranscriptClean,
		   regex(r'(.*.dir)/(.*)/(.*)/(.*).sam'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly_renamed.fa'),
		   r'\1/\2/\3/flagged/\4_labeled.sam')

def flagReads(infiles, outfile):

	# Prefix
	tempdir = os.path.join(os.path.dirname(outfile), 'tmp')
	prefix = outfile[:-len('_labeled.sam')]

	# Command
	cmd_str = ''' talon_label_reads \
		--f {infiles[0]} \
		--g {infiles[1]} \
		--t 1 \
		--ar 20 \
		--tmpDir {tempdir} \
		--deleteTmp \
		--o {prefix}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, W="02:00", GB=15, n=1, conda_env='talon', stdout=outfile.replace('_labeled.sam', '_job.log'), stderr=outfile.replace('_labeled.sam', '_job.err'))

# find arion/isoseq/s04-cleaned_transcripts.dir/*/* -name "flagged" | xargs rm -r
# find arion/isoseq/s04-cleaned_transcripts.dir/*/*/flagged -name "*job.log" | js
# find arion/isoseq/s04-cleaned_transcripts.dir/*/*/flagged -name "*job.err" | xargs wc -l

#######################################################
#######################################################
########## S5. TALON
#######################################################
#######################################################

#############################################
########## 1. Initialize
#############################################

@transform('arion/datasets/reference_genomes/*/*.102.gtf',
		   regex(r'(.*)/(.*)/(.*)(.102).gtf'),
		   add_inputs(r'\1/\2/\3.dna_sm.primary_assembly.fa'),
		   r'arion/isoseq/s05-talon.dir/\2/\3\4_talon.db')

def initializeTalonDatabase(infiles, outfile):

	# Get parameters
	annotation_name = os.path.basename(outfile)[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	outname = outfile[:-len('.db')]

	# Command
	cmd_str = ''' talon_initialize_database \
		--f {infiles[0]} \
		--a {annotation_name} \
		--g {genome_build} \
		--l 200 \
		--idprefix TALON \
		--5p 1000 \
		--3p 1000 \
		--o {outname} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='01:00', n=1, GB=15, conda_env='talon', print_cmd=False, stdout=outfile.replace('.db', '.log'), stderr=outfile.replace('.db', '.err'), wait=True)

#############################################
########## 2. Config
#############################################

@collate(flagReads,
		 regex(r'(.*)/s04-cleaned_transcripts.dir/(.*?)/.*.sam'),
		 r'\1/s05-talon.dir/\2/\2-talon_config.csv')

def createTalonConfig(infiles, outfile):

	# Create dataframe
	config_dataframe = pd.DataFrame([{
		'dataset_name': x.split('/')[-3],
		'sample_description': x.split('/')[-3],
		'platform': 'PacBio-Sequel2' if 'human' in outfile else 'PacBio-Sequel',
		'sam_file': os.path.join(os.getcwd(), x),
	} for x in infiles])

	# outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	config_dataframe.to_csv(outfile, header=False, index=False)

#############################################
########## 3. Run
#############################################

@follows(createTalonConfig)

@transform(initializeTalonDatabase,
		   regex(r'(.*)/(.*).db'),
		   add_inputs(r'\1/*-talon_config.csv'),
		   r'\1/\2_annotated_QC.log')

def runTalon(infiles, outfile):

	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')

	# Fix paths
	outdir = os.path.dirname(outfile)+'/'
	infiles = [os.path.basename(x) for x in infiles]
	outname = os.path.basename(outfile[:-len('_QC.log')])

	# Command
	cmd_str = ''' cd {outdir} && talon \
		--db {infiles[0]} \
		--f {infiles[1]} \
		--build {genome_build} \
		--threads 10 \
		--cov 0.99 \
		--identity 0.95 \
		--o {outname} '''.format(**locals())

	# Run
	# run_job(cmd_str, outfile, W='06:00', n=10, GB=10, conda_env='talon', print_cmd=False, stdout=outfile.replace('_QC.log', '.log'), stderr=outfile.replace('_QC.log', '.err'))

# wc -l arion/isoseq/s05-talon.dir/*/*_QC.log

# human 9,444,700 reads
# mouse 747,163 reads

# find arion/isoseq/s04-cleaned_transcripts.dir/human -name "*clean_labeled.sam" | xargs wc -l
# find arion/isoseq/s04-cleaned_transcripts.dir/mouse -name "*clean_labeled.sam" | xargs wc -l

#############################################
########## 4. Filter transcripts
#############################################

@follows(runTalon)

@transform(initializeTalonDatabase,
		   regex(r'(.*).db'),
		   r'\1_filtered_transcripts.tsv')

def filterTranscripts(infile, outfile):

	# Open database
	with sqlite3.connect(infile) as conn:

		# Query
		query = """ SELECT DISTINCT gene_ID, transcript_ID 
						FROM transcripts
						WHERE transcript_ID NOT IN (
							SELECT DISTINCT ID
								FROM transcript_annotations
								WHERE (attribute = "ISM_transcript" AND value = "TRUE")
								OR (attribute = "genomic_transcript" AND value = "TRUE")
						) AND transcript_ID NOT IN (
							SELECT DISTINCT transcript_ID
								FROM observed
								GROUP BY transcript_ID
								HAVING MIN(fraction_As) > 0.6
						) AND transcript_ID NOT IN (
							SELECT DISTINCT tx.transcript_ID
								FROM transcripts tx
								LEFT JOIN transcript_annotations ta
									ON tx.transcript_ID=ta.ID
								WHERE n_exons = 1 AND source = 'TALON'
						)
		""".format(**locals())

		# Get dataframe
		id_dataframe = pd.read_sql_query(query, conn)

	# Write
	id_dataframe.to_csv(outfile, header=False, index=False)

# find arion/isoseq/s05-talon.dir -name "*filtered_transcripts.tsv" | xargs wc -l
# find arion/isoseq/s05-talon.dir -name "*filtered_transcripts.tsv" | xargs rm
# find arion/isoseq/s05-talon.dir -name "*reference*" | xargs rm
# find arion/isoseq/s05-talon.dir -name "*.gtf" | xargs rm

#############################################
########## 5. Get GTF
#############################################

@follows(filterTranscripts)

@transform(initializeTalonDatabase,
		   regex(r'(.*).db'),
		   add_inputs(r'\1_filtered_transcripts.tsv'),
		   r'\1.gtf')

def getTalonGTF(infiles, outfile):
	
	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	prefix = outfile[:-len('_talon.gtf')]

	# Command
	cmd_str = ''' talon_create_GTF \
		--db {infiles[0]} \
		--annot {annotation_name} \
		--build {genome_build} \
		--whitelist {infiles[1]} \
		--o {prefix} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='00:10', n=1, GB=10, conda_env='talon', print_cmd=False, stdout=outfile.replace('.gtf', '_reference.log'), stderr=outfile.replace('.gtf', '_reference.err'))

#############################################
########## 6. Get FASTA
#############################################

@transform(getTalonGTF,
		   regex(r'(.*)/(.*)/(.*).gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa'),
		   r'\1/\2/\3.fasta')

def getFASTA(infiles, outfile):

	# Command
	cmd_str = ''' gffread {infiles[0]} -g {infiles[1]} -w {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['gff/2021-02'], W='00:30', GB=10, n=1, stdout=outfile.replace('.fasta', '_fasta.log'), stderr=outfile.replace('.fasta', '_fasta.err'))

#############################################
########## 7. Transcript summary
#############################################

@follows(runTalon)

@transform(filterTranscripts,
		   regex(r'(.*)(_filtered_transcript)s.tsv'),
		   add_inputs(r'\1.db'),
		   r'\1\2_summary.tsv')

def getTalonSummary(infiles, outfile):
	
	# Read IDs
	id_dataframe = pd.read_csv(infiles[0], header=None, names=['gene_ID', 'transcript_ID'])

	# Initialize connection
	conn = sqlite3.connect(infiles[1])

	### Genes
	# Gene annotation
	gene_dataframe = pd.read_sql_query(""" SELECT DISTINCT ID as gene_ID, attribute, value FROM gene_annotations """, conn).pivot(index='gene_ID', columns='attribute', values='value')
	gene_dataframe['gene_category'] = ['known' if rowData['gene_status'] == 'KNOWN' else ','.join([x for x in ['intergenic_novel', 'antisense_gene'] if rowData[x] == 'TRUE']) for index, rowData in gene_dataframe.iterrows()]
	gene_dataframe['gene_source'] = gene_dataframe['gene_source'].fillna('talon')
	gene_dataframe['gene_biotype'] = gene_dataframe['gene_biotype'].fillna('not_available')
	gene_dataframe = gene_dataframe[['gene_id', 'gene_name', 'gene_biotype', 'gene_source', 'gene_category']].reset_index()

	### Transcripts
	# Novel transcript annotation
	query = """ SELECT DISTINCT ID AS transcript_ID, attribute, value FROM transcript_annotations WHERE source == 'TALON' AND attribute LIKE '%transcript%' AND attribute NOT LIKE 'ISM-%' """
	novel_transcript_dataframe = pd.read_sql_query(query, conn).pivot(index='transcript_ID', columns='attribute', values='value')
	novel_transcript_dataframe['category'] = [','.join([key.replace('_transcript', '') for key, value in rowData.items() if value == 'TRUE']) for index, rowData in novel_transcript_dataframe.iterrows()]
	novel_transcript_dataframe = novel_transcript_dataframe[['transcript_id', 'transcript_name', 'category']]
	novel_transcript_dataframe['transcript_source'] = 'talon'
	novel_transcript_dataframe['transcript_biotype'] = 'not_available'

	# Known transcript annotation
	known_transcript_dataframe = pd.read_sql_query(""" SELECT DISTINCT ID AS transcript_ID, attribute, value FROM transcript_annotations WHERE source != 'TALON' """, conn).pivot(index='transcript_ID', columns='attribute', values='value').drop(['ccds_id', 'tag', 'transcript_status', 'transcript_support_level', 'transcript_version', 'source'], axis=1)
	known_transcript_dataframe['category'] = 'FSM'

	# Concatenate
	transcript_dataframe = pd.concat([x[['transcript_id', 'transcript_name', 'category', 'transcript_source', 'transcript_biotype']] for x in (known_transcript_dataframe, novel_transcript_dataframe)], axis=0).rename(columns={'category': 'transcript_category'}).reset_index()

	# Add observed column
	# observed_transcript_ids = pd.read_sql_query(""" SELECT DISTINCT transcript_ID FROM observed """, conn)['transcript_ID']
	# transcript_dataframe['observed_transcript'] = [x in observed_transcript_ids for x in transcript_dataframe['transcript_ID']]

	### Merge
	result_dataframe = id_dataframe.merge(gene_dataframe, on='gene_ID', how='left').merge(transcript_dataframe, on='transcript_ID', how='left').drop(['gene_ID', 'transcript_ID'], axis=1)

	# Close
	conn.close()

	# Write
	result_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 8. Get SJs
#############################################

# @transform(getTalonGTF,
@transform('arion/isoseq/s05-talon.dir/*/*.gtf',
		   regex(r'(.*)/(.*)/(.*).gtf'),
		   r'\1/\2/\3_junctions.tsv')

def getJunctions(infile, outfile):

	# Run
	run_r_job('get_junctions', infile, outfile, conda_env='env', W='02:00', GB=20, n=5, run_locally=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 9. Get abundance
#############################################

@follows(filterTranscripts)

@transform('arion/isoseq/s05-talon.dir/*/*.102_talon.db',
		   regex(r'(.*)/(.*).db'),
		   add_inputs(r'\1/\2_filtered_transcripts.tsv'),
		   r'\1/\2_abundance_filtered.tsv')

def getTalonAbundance(infiles, outfile):
	
	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	prefix = outfile[:-len('_talon_abundance_filtered.tsv')]

	# Command
	cmd_str = ''' talon_abundance \
		--db {infiles[0]} \
		--annot {annotation_name} \
		--build {genome_build} \
		--whitelist {infiles[1]} \
		--o {prefix} '''.format(**locals())
		# --observed \

	# Run
	run_job(cmd_str, outfile, W='01:00', n=1, GB=25, conda_env='talon', print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))


# find arion/isoseq/s05-talon.dir -name "*junctions.*" | xargs rm

#######################################################
#######################################################
########## S6. CPAT
#######################################################
#######################################################

#############################################
########## 1. Run
#############################################

@follows(getFASTA)

@transform('arion/isoseq/s05-talon.dir/*/*_talon.fasta',
		   regex(r'(.*)/s05-talon.dir/(.)(.*)/(.*).fasta'),
		   add_inputs(r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_Hexamer.tsv', r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_logitModel.RData'),
		   r'\1/s06-cpat.dir/\2\3/\4-cpat.ORF_prob.best.tsv')

def runCPAT(infiles, outfile):
	
	# Split
	transcript_fasta, cpat_hexamer, cpat_model = infiles
	basename = outfile.replace('.ORF_prob.best.tsv', '')

	# Command
	cmd_str = ''' ~/.conda/envs/env/bin/cpat.py -g {transcript_fasta} \
		-x {cpat_hexamer} \
		-d {cpat_model} \
		--top-orf=5 \
		--log-file={basename}.log \
		-o {basename}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', GB=25, n=1, print_outfile=False, wait=False)

#############################################
########## 2. Get results
#############################################

@transform(runCPAT,
		   suffix('.tsv'),
		   '_results.tsv')

def formatCPAT(infile, outfile):

	# Rename
	rename_dict = {'seq_ID': 'Sequence Name', 'mRNA': 'RNA size', 'ORF': 'ORF size', 'Fickett': 'Ficket Score', 'Hexamer': 'Hexamer Score', 'Coding_prob': 'Coding Probability'}

	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
	coding_cutoff = 0.364 if 'human' in outfile else 0.44

	# Read
	cpat_dataframe = pd.read_table(infile).rename(columns=rename_dict)[rename_dict.values()]
	cpat_dataframe['Coding Label'] = ['yes' if x >= coding_cutoff else 'no' for x in cpat_dataframe['Coding Probability']]
	cpat_dataframe.index.name = 'Data ID'

	# Write
	cpat_dataframe.to_csv(outfile, sep='\t', index=True)

#############################################
########## 3. Split GTF
#############################################

@subdivide(getTalonGTF,
		   regex(r'(.*)/s05-talon.dir/(.*)/(.*).gtf'),
		   r'\1/s06-cpat.dir/\2/gtf/split/\3_??.gtf',
		   r'\1/s06-cpat.dir/\2/gtf/split/\3_{chunk_nr}.gtf')

def splitGTF(infile, outfiles, outfileRoot):

	# Run
	if not len(outfiles):
		run_r_job('split_gtf', infile, outfileRoot, conda_env='env', W='01:00', GB=15, n=1, run_locally=False, wait=False)

# find arion/isoseq/s06-cpat.dir -name "*talon_*.gtf"| xargs wc -l

#############################################
########## 4. Add CDS
#############################################

@follows(formatCPAT)

@transform(splitGTF,
		   regex(r'(.*)/(gtf/split/.*).gtf'),
		   add_inputs(r'\1/*-cpat.ORF_prob.best.tsv'),
		   r'\1/\2.cds.gtf')

def addCDS(infiles, outfile):

	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
	coding_cutoff = 0.364 if 'human' in outfile else 0.44

	# Run
	run_r_job('add_cds', infiles, outfile, additional_params=coding_cutoff, conda_env='env', W='02:00', GB=25, n=1, run_locally=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

# find arion/isoseq/s06-cpat.dir/ -name "*.cds.log" | jsc
# find arion/isoseq/s06-cpat.dir/ -name "*.cds.err" | lr
# find arion/isoseq/s06-cpat.dir/ -name "*.cds.err" | xargs wc -l

#############################################
########## 5. Merge
#############################################

@collate(addCDS,
		 regex(r'(.*)/split/(.*)_.*.cds.gtf'),
		 r'\1/\2.cds.gtf')

def mergeGTF(infiles, outfile):

	# Run
	run_r_job('merge_gtf', infiles, outfile, conda_env='env', W='00:45', GB=15, n=1, run_locally=False)

#######################################################
#######################################################
########## S10. Pfam
#######################################################
#######################################################

#############################################
########## 1. Translate
#############################################

@transform(runCPAT,
		   regex(r'(.*)/(s06-cpat.dir)/(.*)/(.*).ORF_prob.best.tsv'),
		   add_inputs(r'\1/\2/\3/\4.ORF_seqs.fa'),
		   r'\1/s07-pfam.dir/\3/fasta/\3-translated.fasta')

def translateORFs(infiles, outfile):

	# Run
	run_r_job('translate_orfs', infiles, outfile, conda_env='env', W="01:00", GB=15, n=1, stdout=outfile.replace('.fasta', '.log'), wait=True)

#############################################
########## 2. Split
#############################################

@subdivide(translateORFs,
		   regex(r'(.*).fasta'),
		   r'\1.fasta.*',
		   r'\1.fasta.100')

def splitORFs(infile, outfiles, outfileRoot):

	# Get number
	N = outfileRoot.split('.')[-1]

	# Command
	cmd_str = ''' gt splitfasta -numfiles {N} {infile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitORF_'+outfileRoot.split('/')[-3], wait=False)

#############################################
########## 3. Run
#############################################

@transform(splitORFs,
		   regex(r'(.*)/fasta/(.*).fasta\.(.*)'),
		   add_inputs('/sc/arion/work/torred23/libraries/PfamScan'),
		   r'\1/split/\2_\3_pfam.txt')

def runPfamScan(infiles, outfile):

	# Data directory
	input_fasta, pfam_dir = infiles

	# Command
	cmd_str = ''' {pfam_dir}/pfam_scan.pl \
		-dir {pfam_dir}/data \
		-cpu 50 \
		-fasta {input_fasta} > {outfile}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, conda_env='env', modules=['hmmer/3.3'], W='06:00', GB=2, n=10, print_cmd=False, run_locally=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

# find arion/isoseq/s07-pfam.dir -name "*pfam.log" | jsc

#############################################
########## 3. Format
#############################################

@collate(runPfamScan,
		 regex(r'(.*)/split/(.*)_.*_pfam.txt'),
		 r'\1/\2_pfam.tsv')

def mergePfamResults(infiles, outfile):

	# Initialize
	results = []

	# Loop
	for infile in infiles:

		# Read
		pfam_dataframe = pd.read_csv(infile, comment='#', delim_whitespace=True, header=None)

		# Get column names
		colnames = [x.replace(' ', '_').replace('-', '_') for x in pd.read_csv(infile).iloc[26][0][:-1].replace('# ', '').replace('<', '').split('> ')]

		# Add column names
		pfam_dataframe.columns = colnames

		# Fix sequence ID
		pfam_dataframe['seq_id'] = [x.split('_ORF')[0] for x in pfam_dataframe['seq_id']]

		# Append
		results.append(pfam_dataframe)

	# Concatenate
	result_dataframe = pd.concat(results).query('E_value < 0.1').sort_values('seq_id')

	# Write
	result_dataframe.to_csv(outfile, index=False, sep='\t')
	
#######################################################
#######################################################
########## S8. RepeatMasker
#######################################################
#######################################################

#############################################
########## 1. Split
#############################################

@follows(getFASTA)

@subdivide('arion/isoseq/s05-talon.dir/*/*_talon.fasta',
		   regex(r'(.*)/s05-talon.dir/(.*)/(.*).fasta'),
		   r'\1/s08-repeatmasker.dir/\2/fasta/*.fasta.*',
		   r'\1/s08-repeatmasker.dir/\2/fasta/\3.fasta.100')

def splitFASTA(infile, outfiles, outfileRoot):

	# Get number
	N = outfileRoot.split('.')[-1]

	# Get temp filename
	outdir = os.path.dirname(outfileRoot)
	tempfile = os.path.join(outdir, os.path.basename(infile))

	# Command
	cmd_str = ''' cp {infile} {outdir} && gt splitfasta -numfiles {N} {tempfile} && rm {tempfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitFASTA_'+outfileRoot.split('/')[-3], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'))

#############################################
########## 2. Run
#############################################

@transform(splitFASTA,
		   regex(r'(.*)/fasta/(.*)'),
		   r'\1/split/\2.out')

def runRepeatMasker(infile, outfile):

	# Paths
	infile = os.path.join(os.getcwd(), infile)
	outdir = os.path.dirname(outfile)

	# Species
	species = 'Homo sapiens' if 'human' in outfile else 'Mus musculus'

	# Command
	cmd_str = ''' cd {outdir} && RepeatMasker -species "{species}" -pa 64 -dir . {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['repeatmasker/4.1.1'], W='03:00', GB=3, n=10, print_outfile=False, stdout=outfile.replace('.out', '.log'), stderr=outfile.replace('.out', '.err'))

# find arion/isoseq/s08-repeatmasker.dir/*/split -name "*.log" | jsc

#############################################
########## 3. Filter
#############################################

@collate(runRepeatMasker,
		 regex(r'(.*)/split/(.*talon).*.out'),
		 r'\1/\2_repeatmasker.tsv')

def mergeRepeatMasker(infiles, outfile):

	# Run
	run_r_job('merge_repeatmasker', infiles, outfile, conda_env='env', run_locally=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S9. PhyloP
#######################################################
#######################################################

#############################################
########## 1. Convert GTF
#############################################

# @follows(splitGTF)

# @transform('arion/isoseq/s06-cpat.dir/human/gtf/split/*talon_01.gtf',
@transform('arion/isoseq/s06-cpat.dir/*/gtf/split/*talon_??.gtf',
		   regex(r'(.*)/s06-cpat.dir/(.*)/gtf/split/(.*).gtf'),
		   r'\1/s09-phylop.dir/\2/bed/\3.bed')

def gtfToBed(infile, outfile):

	# Run
	run_r_job('gtf_to_bed', infile, outfile, conda_env='env', W='00:05', GB=10, n=1, stdout=outfile.replace('.bed', '.log'), stderr=outfile.replace('.bed', '.err'))

# find arion/isoseq/s09-phylop.dir/*/bed -name "*.log" | jsc

#############################################
########## 2. Convert GTF
#############################################

# @follows(getTalonGTF)

@transform(gtfToBed,
		   regex(r'(.*)/(.*)/(.*)/bed/(.*).bed'),
		   add_inputs(r'arion/datasets/phylop/\3/*.bw'),
		   r'\1/\2/\3/split/\4_phyloP.tsv')

def getPhyloScores(infiles, outfile):

	# Command
	cmd_str = ''' bigWigAverageOverBed {infiles[1]} {infiles[0]} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:15', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# find arion/isoseq/s09-phylop.dir/*/split -name "*.log" | jsc
# find arion/isoseq/s09-phylop.dir/*/split -name "*.err" | lr

#############################################
########## 3. Merge
#############################################

@collate(getPhyloScores,
# @collate('arion/isoseq/s09-phylop.dir/*/split/*_phyloP.tsv',
		 regex(r'(.*)/(.*)/split/.*.tsv'),
		 r'\1/\2/\2-phyloP.tsv')

def mergePhyloScores(infiles, outfile):

	# Run
	run_r_job('merge_phylo_scores', infiles, outfile, conda_env='env', W='00:05', GB=10, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# bigWigAverageOverBed arion/datasets/phylop/human/hg38.phyloP100way.bw test.bed test.tsv

#######################################################
#######################################################
########## S10. liftOver
#######################################################
#######################################################

#############################################
########## 1. Convert
#############################################

# @follows(addCDS)

# @transform('arion/isoseq/s06-cpat.dir/human/gtf/split/Homo_sapiens.GRCh38.102_talon_01.cds.gtf',
@transform('arion/isoseq/s06-cpat.dir/*/gtf/split/*.cds.gtf',
		   regex(r'(.*)/s06-cpat.dir/(.*)/gtf/split/(.*).gtf'),
		   r'\1/s10-liftover.dir/\2/gp/\3.gp')

def convertToGenePred(infile, outfile):

	# Command
	tempfile = outfile.replace('.gp', '.tmp.gp')
	cmd_str = ''' ldHgGene -gtf -nobin -out={tempfile} ignored "" {infile} && awk 'OFS="\t" {{$2="chr"$2; print}}' {tempfile} > {outfile} && rm {tempfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:05', GB=1, n=5, print_cmd=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

# find arion/isoseq/s10-liftover.dir/*/gp/*.log | jsc

#############################################
########## 2. Lift
#############################################

def liftoverJobs():
	for organism in ['human', 'mouse']:
		gp_files = glob.glob('arion/isoseq/s10-liftover.dir/{organism}/gp/*.gp'.format(**locals()))
		genome_str = organism.replace('human', 'hg38').replace('mouse', 'mm10')
		liftover_files = glob.glob('arion/datasets/liftover/{organism}/{genome_str}*.over.chain.gz'.format(**locals()))
		# liftover_file = 'arion/datasets/liftover/human/hg38ToMm10.over.chain.gz' if organism == 'human' else 'arion/datasets/liftover/mouse/mm10ToHg38.over.chain.gz'
		for liftover_file in liftover_files:
			lift_string = os.path.basename(liftover_file).split('.')[0]
			for gp_file in gp_files:
				basename = os.path.basename(gp_file)[:-len('.gp')]
				infiles = [gp_file, liftover_file]
				outfile = 'arion/isoseq/s10-liftover.dir/{organism}/lift/{lift_string}/{basename}-{lift_string}.gp'.format(**locals())
				yield [infiles, outfile]

@follows(convertToGenePred)

@files(liftoverJobs)

def runLiftOver(infiles, outfile):

	# Unmapped
	unmapped_file = outfile.replace('.gp', '_unmapped.gp')

	# Command
	cmd_str = ''' liftOver -genePred {infiles[0]} {infiles[1]} {outfile} {unmapped_file} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['liftover/09-Jul-2019'], W='00:15', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

# find arion/isoseq/s10-liftover.dir/*/lift/*/*.log | jsc
# find arion/isoseq/s10-liftover.dir/*/lift/*/*.err | lr

#############################################
########## 3. Merge
#############################################

@follows(runLiftOver)

@collate('arion/isoseq/s10-liftover.dir/*/lift/*/*.gp',
		 regex(r'(.*)/lift/(.*)/(.*)_..(.cds)-(.*.gp)'),
		 r'\1/merged/\2/\3\4-\5')

def mergeLiftOver(infiles, outfile):

	# String
	infiles_str = ' '.join(infiles)

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Merge
	os.system('cat {infiles_str} > {outfile}'.format(**locals()))

#############################################
########## 4. Filter
#############################################

# @follows(mergeLiftOver)

@transform('arion/isoseq/s10-liftover.dir/*/merged/*/*.gp',
		   regex(r'(.*)/(.*)/(merged)/(.*)(?!d).{1}.gp'),
		   add_inputs(r'arion/illumina/s04-alignment.dir/\2/all/gtf/*-all-SJ_filtered.gtf'),
		    r'\1/\2/\3/\4_filtered.gp')

def filterGenePred(infiles, outfile):

	# Run
	run_r_job('filter_genepred', infiles, outfile, conda_env='env', W='00:05', GB=10, n=1, run_locally=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

#############################################
########## 5. Convert
#############################################

@transform(filterGenePred,
		   suffix('.gp'),
		   '.gtf')

def convertLiftOver(infile, outfile):

	# Command
	cmd_str = ''' genePredToGtf file {infile} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:15', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

# du -hs /hpc/users/torred23/pipelines/projects/early-embryo/arion/isoseq/s10-liftover.dir/*/merged/*.gtf

#######################################################
#######################################################
########## Summary
#######################################################
#######################################################

#############################################
########## 1. Create summary
#############################################

@collate((getTalonSummary, formatCPAT, mergePfamResults, mergeRepeatMasker),
		 regex(r'(.*)/.*.dir/(.*?)/.*'),
		 r'\1/summary.dir/\2-isoseq_summary.tsv')

def getTranscriptSummary(infiles, outfile):

	# Run
	run_r_job('get_transcript_summary', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1)

# #######################################################
# #######################################################
# ########## S9. CPAT
# #######################################################
# #######################################################

# #############################################
# ########## 1. Run
# #############################################

# @follows(runSqantiFinalPass)

# @transform('arion/isoseq/s08-sqanti_pass2.dir/*/filtered/*_classification.filtered_lite_no_ism_corrected.fasta',
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/.(.*)_classification.filtered_lite_no_ism_corrected.fasta'),
# 		   add_inputs(r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_Hexamer.tsv', r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_logitModel.RData'),
# 		   r'\1/s09-cpat.dir/\2/\2-cpat.ORF_prob.best.tsv')

# def runCPAT(infiles, outfile):
	
# 	# Split
# 	transcript_fasta, cpat_hexamer, cpat_model = infiles
# 	basename = outfile.replace('.ORF_prob.best.tsv', '')

# 	# Command
# 	cmd_str = ''' ~/.conda/envs/env/bin/cpat.py -g {transcript_fasta} \
# 		-x {cpat_hexamer} \
# 		-d {cpat_model} \
# 		--top-orf=5 \
# 		--log-file={basename}.log \
# 		-o {basename}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='06:00', GB=25, n=1, run_locally=False, print_outfile=False, wait=True)

# #############################################
# ########## 2. Get results
# #############################################

# @transform(runCPAT,
# 		   suffix('.tsv'),
# 		   '_results.tsv')

# def formatCPAT(infile, outfile):

# 	# Rename
# 	rename_dict = {'seq_ID': 'Sequence Name', 'mRNA': 'RNA size', 'ORF': 'ORF size', 'Fickett': 'Ficket Score', 'Hexamer': 'Hexamer Score', 'Coding_prob': 'Coding Probability'}

# 	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
# 	coding_cutoff = 0.364 if 'human' in outfile else 0.44

# 	# Read
# 	cpat_dataframe = pd.read_table(infile).rename(columns=rename_dict)[rename_dict.values()]
# 	cpat_dataframe['Coding Label'] = ['yes' if x >= coding_cutoff else 'no' for x in cpat_dataframe['Coding Probability']]
# 	cpat_dataframe.index.name = 'Data ID'

# 	# Write
# 	cpat_dataframe.to_csv(outfile, sep='\t', index=True)

# #############################################
# ########## 3. Split GTF
# #############################################

# @subdivide(runSqantiFinalPass,
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/(.*?)_.*.gtf'),
# 		   r'\1/s09-cpat.dir/\2/gtf/split/\3_??.gtf',
# 		   r'\1/s09-cpat.dir/\2/gtf/split/\3_{chunk_nr}.gtf')

# def splitGTF(infile, outfiles, outfileRoot):

# 	# Run
# 	run_r_job('split_gtf', infile, outfileRoot, conda_env='env', W='01:00', GB=15, n=1, run_locally=False, wait=True)

# #############################################
# ########## 4. Add CDS
# #############################################

# @follows(formatCPAT)

# @transform(splitGTF,
# 		   regex(r'(.*)/(gtf/split/.*).gtf'),
# 		   add_inputs(r'\1/*-cpat.ORF_prob.best.tsv'),
# 		   r'\1/\2.cds.gtf')

# def addCDS(infiles, outfile):

# 	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
# 	coding_cutoff = 0.364 if 'human' in outfile else 0.44

# 	# Run
# 	run_r_job('add_cds', infiles, outfile, additional_params=coding_cutoff, conda_env='env', W='03:00', GB=15, n=1, run_locally=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds
# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds.log | jsc
# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds.err | lr

# #############################################
# ########## 5. Merge
# #############################################

# @collate(addCDS,
# 		 regex(r'(.*)/split/(.*)_.*.cds.gtf'),
# 		 r'\1/\2.cds.gtf')

# def mergeGTF(infiles, outfile):

# 	# Run
# 	run_r_job('merge_gtf', infiles, outfile, conda_env='env', W='00:45', GB=15, n=1, run_locally=False)

# #######################################################
# #######################################################
# ########## S10. Pfam
# #######################################################
# #######################################################

# #############################################
# ########## 1. Translate
# #############################################

# @transform(runCPAT,
# 		   regex(r'(.*)/(s09-cpat.dir)/(.*)/(.*).ORF_prob.best.tsv'),
# 		   add_inputs(r'\1/\2/\3/\4.ORF_seqs.fa'),
# 		   r'\1/s10-pfam.dir/\3/fasta/\3-translated.fasta')

# def translateORFs(infiles, outfile):

# 	# Run
# 	run_r_job('translate_orfs', infiles, outfile, conda_env='env', W="01:00", GB=15, n=1, stdout=outfile.replace('.fasta', '.log'), wait=True)

# #############################################
# ########## 2. Split
# #############################################

# @subdivide(translateORFs,
# 		   regex(r'(.*).fasta'),
# 		   r'\1.fasta.*',
# 		   r'\1.fasta.100')

# def splitORFs(infile, outfiles, outfileRoot):

# 	# Get number
# 	N = outfileRoot.split('.')[-1]

# 	# Command
# 	cmd_str = ''' gt splitfasta -numfiles {N} {infile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitORF_'+outfileRoot.split('/')[-3], wait=True)

# #############################################
# ########## 3. Run
# #############################################

# @transform(splitORFs,
# 		   regex(r'(.*)/fasta/(.*).fasta\.(.*)'),
# 		   add_inputs('/sc/arion/work/torred23/libraries/PfamScan'),
# 		   r'\1/split/\2_\3_pfam.txt')

# def runPfamScan(infiles, outfile):

# 	# Data directory
# 	input_fasta, pfam_dir = infiles

# 	# Command
# 	cmd_str = ''' {pfam_dir}/pfam_scan.pl \
# 		-dir {pfam_dir}/data \
# 		-cpu 50 \
# 		-fasta {input_fasta} > {outfile}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', modules=['hmmer/3.3'], W='06:00', GB=2, n=10, print_cmd=False, run_locally=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

# # find /hpc/users/torred23/pipelines/projects/early-embryo/arion/isoseq/s10-pfam.dir -name "*pfam.log" | jsc

# #############################################
# ########## 3. Format
# #############################################

# @collate(runPfamScan,
# 		 regex(r'(.*)/split/(.*)_.*_pfam.txt'),
# 		 r'\1/\2_pfam.tsv')

# def mergePfamResults(infiles, outfile):

# 	# Initialize
# 	results = []

# 	# Loop
# 	for infile in infiles:

# 		# Read
# 		pfam_dataframe = pd.read_csv(infile, comment='#', delim_whitespace=True, header=None)

# 		# Get column names
# 		colnames = [x.replace(' ', '_').replace('-', '_') for x in pd.read_csv(infile).iloc[26][0][:-1].replace('# ', '').replace('<', '').split('> ')]

# 		# Add column names
# 		pfam_dataframe.columns = colnames

# 		# Fix sequence ID
# 		pfam_dataframe['seq_id'] = [x.split('_ORF')[0] for x in pfam_dataframe['seq_id']]

# 		# Append
# 		results.append(pfam_dataframe)

# 	# Concatenate
# 	result_dataframe = pd.concat(results).query('E_value < 0.1').sort_values('seq_id')

# 	# Write
# 	result_dataframe.to_csv(outfile, index=False, sep='\t')
	
# #######################################################
# #######################################################
# ########## S11. RepeatMasker
# #######################################################
# #######################################################

# #############################################
# ########## 1. Split
# #############################################

# # @follows(runSqantiFinalPass)

# @subdivide('arion/isoseq/s08-sqanti_pass2.dir/*/filtered/*_classification.filtered_lite_no_ism_corrected.fasta',
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/(.*).fasta'),
# 		   r'\1/s11-repeatmasker.dir/\2/fasta/*.fasta.*',
# 		   r'\1/s11-repeatmasker.dir/\2/fasta/\3.fasta.100')

# def splitFASTA(infile, outfiles, outfileRoot):

# 	# Get number
# 	N = outfileRoot.split('.')[-1]

# 	# Get temp filename
# 	outdir = os.path.dirname(outfileRoot)
# 	tempfile = os.path.join(outdir, os.path.basename(infile))

# 	# Command
# 	cmd_str = ''' cp {infile} {outdir} && gt splitfasta -numfiles {N} {tempfile} && rm {tempfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitFASTA_'+outfileRoot.split('/')[-3], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'))

# #############################################
# ########## 2. Run
# #############################################

# @transform(splitFASTA,
# 		   regex(r'(.*)/fasta/(.*)'),
# 		   r'\1/split/\2.out')

# def runRepeatMasker(infile, outfile):

# 	# Paths
# 	infile = os.path.join(os.getcwd(), infile)
# 	outdir = os.path.dirname(outfile)

# 	# Species
# 	species = 'Homo sapiens' if 'human' in outfile else 'Mus musculus'

# 	# Command
# 	cmd_str = ''' cd {outdir} && RepeatMasker -species "{species}" -pa 64 -dir . {infile}'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, modules=['repeatmasker/4.1.1'], W='03:00', GB=3, n=10, print_outfile=False, stdout=outfile.replace('.tbl', '.log'), stderr=outfile.replace('.tbl', '.err'))

# #############################################
# ########## 2. Filter
# #############################################

# @collate(runRepeatMasker,
# 		 regex(r'(.*)/split/(.*?)_.*.out'),
# 		 r'\1/\2_repeatmasker.tsv')

# def mergeRepeatMasker(infiles, outfile):

# 	# Run
# 	run_r_job('merge_repeatmasker', infiles, outfile, conda_env='env', run_locally=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# find arion/isoseq/s05-talon.dir -name "*talon.log" | js

# @transform('arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.102.gtf',
# 		   regex(r'(.*)/(.*)(.102).gtf'),
# 		   add_inputs(r'\1/\2.dna_sm.primary_assembly.fa'),
# 		   r'arion/isoseq/talon/SJs/\2\3_SJs.tsv')
# @files(('arion/isoseq/talon/SJs/gencode.v29.annotation.gtf', 'arion/isoseq/talon/SJs/GRCh38.primary_assembly.genome.fa'),
# 	   'arion/isoseq/talon/SJs/test_SJs.tsv')

# def getReferenceSJs(infiles, outfile):

# 	# Command
# 	cmd_str = ''' python /sc/arion/work/torred23/libraries/TranscriptClean/TranscriptClean-2.0.2/accessory_scripts/get_SJs_from_gtf.py \
# 		--f {infiles[0]} \
# 		--g {infiles[1]} \
# 		--o {outfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W='01:00', n=1, GB=15, conda_env='talon', print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# get SJs and run TranscriptClean


# def collapseJobs():
	
# 	# Get SAM
# 	# sam_files = glob.glob('arion/isoseq/s02-alignment.dir/minimap2/mouse/*1C/*-minimap2.filtered.sam')
# 	sam_files = glob.glob('arion/isoseq/s02-alignment.dir/minimap2/*/*/*-minimap2.filtered.sam')

# 	# Parameters
# 	parameter_list = []
# 	# for max_5_diff in [0, 10, 100, 500, 750, 1000, 1500, 2000, 3000, 5000]:# [1000]:
# 	for max_5_diff in [1000]:
# 		for max_fuzzy_junction in [5]:
# 			# for max_3_diff in [0, 10, 100, 500, 750, 1000, 1500, 2000, 3000, 5000]:# [1000]:
# 			for max_3_diff in [1000]:
# 				parameter_list.append({'max_5_diff': max_5_diff, 'max_fuzzy_junction': max_fuzzy_junction, 'max_3_diff': max_3_diff})
# 	del max_5_diff
# 	del max_fuzzy_junction
# 	del max_3_diff

# 	# Loop
# 	for infile in sam_files:
# 		sample_name = infile.split('/')[-2]
# 		organism = infile.split('/')[-3]
# 		for parameters in parameter_list:
# 			parameters = parameters.copy()
# 			parameters['basename'] = 'arion/isoseq/s03-collapsed.dir/{organism}/{sample_name}/{sample_name}-5p{max_5_diff}-J{max_fuzzy_junction}-3p{max_3_diff}/{sample_name}-5p{max_5_diff}-J{max_fuzzy_junction}-3p{max_3_diff}'.format(**locals(), **parameters)
# 			parameters['fastq'] = 'arion/isoseq/s02-alignment.dir/fastq/{organism}/{sample_name}.flnc.fastq'.format(**locals())
# 			parameters['genome_fasta'] = glob.glob('arion/datasets/reference_genomes/{organism}/*.primary_assembly.fa'.format(**locals()))[0]
# 			outfile = '{basename}.collapsed.filtered.gff'.format(**parameters)
# 			yield [infile, outfile, parameters]

# @follows(filterSam)

# @files(collapseJobs)

# def collapseIsoforms(infile, outfile, parameters):

# 	# SAM file
# 	sam = infile

# 	# Collapse - add --gen_mol_count perhaps
# 	cmd_str = ''' collapse_isoforms_by_sam.py \
# 		--input {fastq} --fq \
# 		--sam {sam} \
# 		--prefix {basename} \
# 		--max_5_diff {max_5_diff} \
# 		--max_fuzzy_junction {max_fuzzy_junction} \
# 		--max_3_diff {max_3_diff} \
# 		--dun-merge-5-shorter \
# 		--gen_mol_count \
# 		--min-coverage 0.99 \
# 		--min-identity 0.95 && \
# 		filter_away_subset.py {basename}.collapsed && \
# 		get_seq_stats.py {basename}.collapsed.filtered.rep.fq
# 	'''.format(**locals(), **parameters) # add deletion of rep.fq files as they are unused

# 	# Run
# 	W = '02:00' if 'mouse' in outfile else '12:00'
# 	GB = 50 if 'human' in outfile else 10
# 	run_job(cmd_str, outfile, modules=['cdna_cupcake/21.0.0'], W=W, GB=GB, n=1, print_outfile=False, stdout=outfile.rsplit('.', 1)[0]+'.log', stderr=outfile.rsplit('.', 1)[0]+'.err')

# #############################################
# ########## 2. SQANTI
# #############################################

# @follows(collapseIsoforms)

# @transform('arion/isoseq/s03-collapsed.dir/*/*/*/*.collapsed.filtered.gff',
# 		   regex(r'(.*)/(.*)/s03-collapsed.dir/(.*)/(.*)/(.*)/.*.gff'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa'),
# 		   r'\1/\2/s03-collapsed.dir/\3/\4/\5/sqanti/\5_sqanti_report.pdf')

# def runSqanti(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if isinstance(x, str) and x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta = infiles

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_sqanti_report.pdf')]

# 	# PYTHONPATH
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PYTHONPATH={PYTHONPATH} && cd {outdir} && sqanti3_qc.py -v && sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 		--gtf \
# 		--dir . \
# 		--output {output}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v1.6.1', W='03:00', GB=10, n=1, print_cmd=True, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'))

# #############################################
# ########## 3. Statistics
# #############################################

# @follows(runSqanti)

# @transform('arion/isoseq/s03-collapsed.dir/*/*/*/sqanti/*sqanti_report.pdf',
# # @transform(runSqanti,
# 		   regex(r'(.*)_sqanti_report.pdf'),
# 		   add_inputs(r'\1_corrected.gtf', r'\1_classification.txt*'),
# 		   r'\1.transcript_stats.tsv')

# def getCollapseStatistics(infiles, outfile):

# 	# Remove PDF
# 	infiles = [x for x in infiles if 'pdf' not in x]
	
# 	# Run
# 	run_r_job('get_collapse_statistics', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1, run_locally=False)

# #############################################
# ########## 4. Make report
# #############################################

# @follows(getCollapseStatistics)

# @collate('arion/isoseq/s03-collapsed.dir/*/*/*/sqanti/*stats.tsv',
# 		 regex(r'(.*)/.*/(.*)/.*/sqanti/.*.tsv'),
# 		 r'\1/reports/\2-collapse_report.pdf')

# def makeCollapseReport(infiles, outfile):

# 	# Run
# 	run_r_job('make_collapse_report', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1, run_locally=False)

# #######################################################
# #######################################################
# ########## S4. Merge
# #######################################################
# #######################################################

# #############################################
# ########## 1. Create links
# #############################################

# def linkJobs():

# 	# Loop
# 	for species in ['mouse', 'human']:

# 		# Get infiles
# 		collapse_settings = '5p1000-J5-3p1000'
# 		infiles = glob.glob('arion/isoseq/s03-collapsed.dir/{species}/*/*-{collapse_settings}/*'.format(**locals()))

# 		# Loop
# 		for infile in infiles:

# 			# Filter
# 			if infile.endswith('.collapsed.filtered.gff') or infile.endswith('.group.txt') or infile.endswith('.filtered.abundance.txt'):# or infile.endswith('.filtered.rep.fq'):

# 				# Get outname
# 				sample_name = os.path.basename(infile).split('-')[0]
# 				file_suffix = os.path.basename(infile).split('.', 1)[-1]

# 				# Get outfile
# 				infile = os.path.join(os.getcwd(), infile)
# 				outfile = os.path.join(os.getcwd(), 'arion/isoseq/s04-merged.dir/{species}/links/{sample_name}/{file_suffix}'.format(**locals()))

# 				# Yield
# 				yield [infile, outfile]

# @follows(collapseIsoforms)

# @files(linkJobs)

# def createChainLinks(infile, outfile):

# 	# Create directory
# 	outdir = os.path.dirname(outfile)
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)

# 	# Create link
# 	os.system('ln -s {infile} {outfile}'.format(**locals()))

# #############################################
# ########## 2. Create config
# #############################################

# @collate(createChainLinks,
# 		 regex(r'(.*)/links/.*'),
# 		 r'\1/results/sample.config')

# def createChainConfig(infiles, outfile):

# 	# Get samples string
# 	sample_dict = {x.split('/')[-2]: os.path.join(os.getcwd(), os.path.dirname(x)) for x in infiles}
# 	sample_str = '\n'.join(['SAMPLE={key};{value}'.format(**locals()) for key, value in sample_dict.items()])

# 	# Get string
# 	config_str = '''{sample_str}
# GROUP_FILENAME=collapsed.group.txt
# GFF_FILENAME=collapsed.filtered.gff
# COUNT_FILENAME=collapsed.filtered.abundance.txt'''.format(**locals()) #FASTQ_FILENAME=collapsed.filtered.rep.fq

# 	# Create directory
# 	outdir = os.path.dirname(outfile)
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)

# 	# Write
# 	with open(outfile, 'w') as openfile:
# 		openfile.write(config_str)
		
# #############################################
# ########## 3. Chain samples
# #############################################

# @transform(createChainConfig,
# 		   suffix('sample.config'),
# 		   'all_samples.chained.gff')

# def chainSamples(infile, outfile):

# 	# Command
# 	cmd_str = ''' cd $(dirname {infile}) && chain_samples.py $(basename {infile}) count_fl && rm tmp_*'''.format(**locals())

# 	# Run
# 	# GB = 50 if 'human' in outfile else 5
# 	# run_job(cmd_str, outfile, modules=['cdna_cupcake/21.0.0'], W='03:00', GB=None, R='himem', q='premium', stdout=outfile.replace('.gff', '.log'), stderr=outfile.replace('.gff', '.err'))
# 	run_job(cmd_str, outfile, modules=['cdna_cupcake/21.0.0'], W='03:00', GB=50, n=5, q = 'premium', stdout=outfile.replace('.gff', '.log'), stderr=outfile.replace('.gff', '.err'))

# #############################################
# ########## 4. Fix gene IDs
# #############################################

# @transform(chainSamples,
# 		   suffix('.gff'),
# 		   '.fixed.gtf')

# def fixChainedGFF(infile, outfile):

# 	# Run
# 	run_r_job('fix_chained_gff', infile, outfile, conda_env='env', W='02:00', GB=15, n=1)

# #######################################################
# #######################################################
# ########## S5. SQANTI First Pass
# #######################################################
# #######################################################

# #############################################
# ########## 1. SQANTI 1
# #############################################

# @transform(fixChainedGFF,
# 		   regex(r'(.*)/(.*)/s04-merged.dir/(.*)/results/.*.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa'),
# 		   r'\1/\2/s05-sqanti_pass1.dir/\3/\3_corrected.gtf')

# def runSqantiPass1(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if isinstance(x, str) and x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta = infiles

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_corrected.gtf')]

# 	# PATHS
# 	PATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PATH={PATH}:$PATH && export PYTHONPATH={PYTHONPATH} && cd {outdir} && \
# 		sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 			--gtf \
# 			--dir . \
# 			--skipORF \
# 			--output {output}
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v1.6.1', W='06:00', GB=50, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'), wait=True)
	
# #######################################################
# #######################################################
# ########## IS5. Alignment indices
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR
# #############################################

# @transform(runSqantiPass1,
# 		   regex(r'(.*)/(.*)/s05-sqanti_pass1.dir/(.*)/.*_corrected.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa'),
# 		   r'\1/\2/s06-indices.dir/\3/STAR/')

# def buildStarIndex(infiles, outfile):

# 	# Split
# 	reference_gtf, reference_fasta = infiles

# 	# Command
# 	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {reference_fasta} --sjdbGTFfile {reference_gtf} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='03:00', GB=5, n=15, ow=True, print_cmd=False, jobname='STAR_index', wait=True)

# #############################################
# ########## 2. Salmon
# #############################################

# @follows(runSqantiPass1)

# @transform('arion/isoseq/s05-sqanti_pass1.dir/*/*_corrected.fasta',
# 		   regex(r'(.*)/s05-sqanti_pass1.dir/(.*)/.*_corrected.fasta'),
# 		   r'\1/s06-indices.dir/\2/salmon')

# def buildSalmonIndex(infile, outfile):

# 	# Command
# 	cmd_str = '''salmon index -t {infile} -i {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['salmon/1.4.0'], W='03:00', GB=8, n=5, print_cmd=False, jobname='salmon_index', wait=True)

# #######################################################
# #######################################################
# ########## S6. Illumina alignment
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR pass 1
# #############################################

# @follows(buildStarIndex)

# @collate(illumina_fastq,
# 		regex(r'.*/s01-fastq.dir/(.*)/(.*)/.*.fastq.gz'),
# 		add_inputs(r'arion/isoseq/s06-indices.dir/\1/STAR/'),
# 		r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass1/\2/\2-SJ.out.tab')

# def getSJsPass1(infiles, outfile):

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[0][1]} \
# 		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 50 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", GB=5, n=10, modules=['star/2.7.5b'], print_cmd=False)

# # find arion/isoseq/s07-illumina_alignment.dir/*/STAR -name "*Log.progress.out" | du -hs
# # ls arion/isoseq/s07-illumina_alignment.dir/*/STAR/pass1/*/*-Log.progress.out | lr
# # ls -lhrt arion/isoseq/s07-illumina_alignment.dir/*/STAR/pass1/*/*-SJ.out.tab

# #############################################
# ########## 2. STAR pass 2
# #############################################

# @follows(getSJsPass1)

# @collate(illumina_fastq,
# 		regex(r'.*/s01-fastq.dir/(.*)/(.*)/.*.fastq.gz'),
# 		add_inputs(r'arion/isoseq/s06-indices.dir/\1/STAR/', r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass1/*/*-SJ.out.tab'),
# 		r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass2/\2/\2-SJ.out.tab')

# def getSJsPass2(infiles, outfile):

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]
# 	sj_files_str = ' '.join(infiles[0][2:])

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[0][1]} \
# 		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 50 \
# 		--sjdbFileChrStartEnd {sj_files_str} \
# 		--limitSjdbInsertNsj 5000000 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="10:00", GB=10, n=10, modules=['star/2.7.5b'], print_cmd=False, stdout=outfile.replace('-SJ.out.tab', '-job.log'), q='premium', wait=False)

# # find arion/isoseq/s07-illumina_alignment.dir/*/STAR/pass2/* -name "*job.log" | jsc
# # ls arion/isoseq/s07-illumina_alignment.dir/*/STAR/pass2/*/*-Log.progress.out | lr

# #############################################
# ########## 3. Salmon
# #############################################

# @follows(buildSalmonIndex)

# @collate(illumina_fastq,
# 		regex(r'.*/s01-fastq.dir/(.*)/(.*)/.*.fastq.gz'),
# 		add_inputs(r'arion/isoseq/s06-indices.dir/\1/salmon'),
# 		r'arion/isoseq/s07-illumina_alignment.dir/\1/salmon/samples/\2/quant.sf')

# def runSalmon(infiles, outfile):

# 	# Split
# 	fastq_files, salmon_index = infiles

# 	# Prefix
# 	outdir = os.path.dirname(outfile)

# 	# Command
# 	cmd_str = '''salmon quant -i {infiles[0][1]} -l A -1 {infiles[0][0]} -2 {infiles[1][0]} --validateMappings -p 200 -o {outdir}'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, modules=['salmon/1.4.0'], W='02:00', GB=2, n=25, print_outfile=False, jobname=outfile.split('/')[-2]+'_salmon', stdout=outfile.replace('quant.sf', 'job.log'))

# # find arion/isoseq/s07-illumina_alignment.dir/*/salmon/samples/ -name "job.log" | jsc

# #############################################
# ########## 4. Aggregate counts
# #############################################

# @collate(runSalmon,
# 		 regex(r'(.*)/(.*)/salmon/.*.quant.sf'),
# 		 r'\1/\2/salmon/\2-transcript_tpm.tsv')

# def aggregateSalmonCounts(infiles, outfile):

# 	# Run
# 	run_r_job('aggregate_salmon_counts', infiles, outfile, conda_env='env', W='00:30', n=1, GB=10, run_locally=True, print_outfile=False)

# #######################################################
# #######################################################
# ########## IS7. SQANTI Second Pass
# #######################################################
# #######################################################

# #############################################
# ########## 1. QC
# #############################################

# @follows(chainSamples, aggregateSalmonCounts, getSJsPass2)

# # @transform('arion/isoseq/s05-sqanti_pass1.dir/*/*_corrected.gtf',
# @transform(runSqantiPass1,
# 		   regex(r'(.*)/(.*)/s05-sqanti_pass1.dir/(.*)/.*.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa', r'arion/isoseq/s07-illumina_alignment.dir/\3/salmon/\3-transcript_tpm.tsv', r'arion/isoseq/s04-merged.dir/\3/results/all_samples.chained_count.txt', r'arion/isoseq/s07-illumina_alignment.dir/\3/STAR/pass2/*/*-SJ.out.tab'),
# 		   r'\1/\2/s08-sqanti_pass2.dir/\3/\3_corrected.gtf')

# def runSqantiPass2(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta, salmon_tpm, fl_counts = infiles[:5]
# 	sj_files = ','.join(infiles[5:])

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_corrected.gtf')]

# 	# PATHS
# 	PATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PATH={PATH}:$PATH && export PYTHONPATH={PYTHONPATH} && cd {outdir} && \
# 		sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 			--gtf \
# 			--skipORF \
# 			--expression {salmon_tpm} \
# 			--fl_count {fl_counts} \
# 			--coverage "{sj_files}" \
# 			--dir . \
# 			--output {output}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v1.6.1', W='06:00', GB=50, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'), wait=True)
	
# #############################################
# ########## 2. Filter
# #############################################

# @transform(runSqantiPass2,
# 		   regex(r'(.*)/(.*)_corrected.gtf'),
# 		   add_inputs(r'\1/\2_corrected.fasta', r'\1/\2_classification.txt'),
# 		   r'\1/\2_classification.filtered_lite.gtf')

# def filterSqanti(infiles, outfile):
	
# 	# Get infiles
# 	infiles = [os.path.join(os.getcwd(), x) if isinstance(x, str) and x.startswith('arion') else x for x in infiles]
# 	gtf_file, fasta_file, classification_file = infiles

# 	# Get output
# 	outdir = os.path.dirname(outfile)

# 	# PATHS
# 	PATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PATH={PATH}:$PATH && export PYTHONPATH={PYTHONPATH} && cd {outdir} && \
# 		sqanti3_RulesFilter.py {classification_file} {fasta_file} {gtf_file} \
# 			--intrapriming 0.6 \
# 			--runAlength 6 \
# 			--max_dist_to_known_end 50 \
# 			--filter_mono_exonic \
# 			--min_cov 3
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v1.6.1', W='03:00', GB=50, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'), wait=True)
	
# # usage: sqanti3_RulesFilter.py [-h] [--sam SAM] [--faa FAA] [-a INTRAPRIMING]
# #                               [-r RUNALENGTH] [-m MAX_DIST_TO_KNOWN_END]
# #                               [-c MIN_COV] [--filter_mono_exonic] [--skipGTF]
# #                               [--skipFaFq] [--skipJunction] [-v]
# #                               sqanti_class isoforms gtf_file

# # Filtering of Isoforms based on SQANTI3 attributes

# # positional arguments:
# #   sqanti_class          SQANTI classification output file.
# #   isoforms              fasta/fastq isoform file to be filtered by SQANTI3
# #   gtf_file              GTF of the input fasta/fastq

# # optional arguments:
# #   -h, --help            show this help message and exit
# #   --sam SAM             (Optional) SAM alignment of the input fasta/fastq
# #   --faa FAA             (Optional) ORF prediction faa file to be filtered by
# #                         SQANTI3
# #   -a INTRAPRIMING, --intrapriming INTRAPRIMING
# #                         Adenine percentage at genomic 3' end to flag an
# #                         isoform as intra-priming (default: 0.6)
# #   -r RUNALENGTH, --runAlength RUNALENGTH
# #                         Continuous run-A length at genomic 3' end to flag an
# #                         isoform as intra-priming (default: 6)
# #   -m MAX_DIST_TO_KNOWN_END, --max_dist_to_known_end MAX_DIST_TO_KNOWN_END
# #                         Maximum distance to an annotated 3' end to preserve as
# #                         a valid 3' end and not filter out (default: 50bp)
# #   -c MIN_COV, --min_cov MIN_COV
# #                         Minimum junction coverage for each isoform (only used
# #                         if min_cov field is not 'NA'), default: 3
# #   --filter_mono_exonic  Filter out all mono-exonic transcripts (default: OFF)
# #   --skipGTF             Skip output of GTF
# #   --skipFaFq            Skip output of isoform fasta/fastq
# #   --skipJunction        Skip output of junctions file
# #   -v, --version         Display program version number.

# #############################################
# ########## 3. ISM filter
# #############################################

# @follows(chainSamples)

# @transform(filterSqanti,
# 		   regex(r'(.*)/(.*).gtf'),
# 		   add_inputs(r'\1/\2_classification.txt'),
# 		   r'\1/filtered/\2_no_ism.gtf')

# def filterISMs(infiles, outfile):

# 	# Run
# 	run_r_job('filter_isms', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1, wait=True)

# #############################################
# ########## 4. SQANTI
# #############################################

# # @transform('arion/isoseq/s08-sqanti_pass2.dir/*/filtered/*_classification.filtered_lite_no_ism.gtf',
# @transform(filterISMs,
# 		   regex(r'(.*)/(.*)/s08-sqanti_pass2.dir/(.*?)/(.*).gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa', r'arion/isoseq/s07-illumina_alignment.dir/\3/salmon/\3-transcript_tpm.tsv', r'arion/isoseq/s04-merged.dir/\3/results/all_samples.chained_count.txt', r'arion/isoseq/s07-illumina_alignment.dir/\3/STAR/pass2/*/*-SJ.out.tab'),
# 		   r'\1/\2/s08-sqanti_pass2.dir/\3/\4_corrected.gtf')

# def runSqantiFinalPass(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta, salmon_tpm, fl_counts = infiles[:5]
# 	sj_files = ','.join(infiles[5:])

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_corrected.gtf')]

# 	# PATHS
# 	PATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3_v1.6.1/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PATH={PATH}:$PATH && export PYTHONPATH={PYTHONPATH} && cd {outdir} && \
# 		sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 			--gtf \
# 			--skipORF \
# 			--expression {salmon_tpm} \
# 			--fl_count {fl_counts} \
# 			--coverage "{sj_files}" \
# 			--dir . \
# 			--output {output}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v1.6.1', W='06:00', GB=50, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'), wait=True)

# #######################################################
# #######################################################
# ########## S9. CPAT
# #######################################################
# #######################################################

# #############################################
# ########## 1. Run
# #############################################

# @follows(runSqantiFinalPass)

# @transform('arion/isoseq/s08-sqanti_pass2.dir/*/filtered/*_classification.filtered_lite_no_ism_corrected.fasta',
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/.(.*)_classification.filtered_lite_no_ism_corrected.fasta'),
# 		   add_inputs(r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_Hexamer.tsv', r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_logitModel.RData'),
# 		   r'\1/s09-cpat.dir/\2/\2-cpat.ORF_prob.best.tsv')

# def runCPAT(infiles, outfile):
	
# 	# Split
# 	transcript_fasta, cpat_hexamer, cpat_model = infiles
# 	basename = outfile.replace('.ORF_prob.best.tsv', '')

# 	# Command
# 	cmd_str = ''' ~/.conda/envs/env/bin/cpat.py -g {transcript_fasta} \
# 		-x {cpat_hexamer} \
# 		-d {cpat_model} \
# 		--top-orf=5 \
# 		--log-file={basename}.log \
# 		-o {basename}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='06:00', GB=25, n=1, run_locally=False, print_outfile=False, wait=True)

# #############################################
# ########## 2. Get results
# #############################################

# @transform(runCPAT,
# 		   suffix('.tsv'),
# 		   '_results.tsv')

# def formatCPAT(infile, outfile):

# 	# Rename
# 	rename_dict = {'seq_ID': 'Sequence Name', 'mRNA': 'RNA size', 'ORF': 'ORF size', 'Fickett': 'Ficket Score', 'Hexamer': 'Hexamer Score', 'Coding_prob': 'Coding Probability'}

# 	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
# 	coding_cutoff = 0.364 if 'human' in outfile else 0.44

# 	# Read
# 	cpat_dataframe = pd.read_table(infile).rename(columns=rename_dict)[rename_dict.values()]
# 	cpat_dataframe['Coding Label'] = ['yes' if x >= coding_cutoff else 'no' for x in cpat_dataframe['Coding Probability']]
# 	cpat_dataframe.index.name = 'Data ID'

# 	# Write
# 	cpat_dataframe.to_csv(outfile, sep='\t', index=True)

# #############################################
# ########## 3. Split GTF
# #############################################

# @subdivide(runSqantiFinalPass,
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/(.*?)_.*.gtf'),
# 		   r'\1/s09-cpat.dir/\2/gtf/split/\3_??.gtf',
# 		   r'\1/s09-cpat.dir/\2/gtf/split/\3_{chunk_nr}.gtf')

# def splitGTF(infile, outfiles, outfileRoot):

# 	# Run
# 	run_r_job('split_gtf', infile, outfileRoot, conda_env='env', W='01:00', GB=15, n=1, run_locally=False, wait=True)

# #############################################
# ########## 4. Add CDS
# #############################################

# @follows(formatCPAT)

# @transform(splitGTF,
# 		   regex(r'(.*)/(gtf/split/.*).gtf'),
# 		   add_inputs(r'\1/*-cpat.ORF_prob.best.tsv'),
# 		   r'\1/\2.cds.gtf')

# def addCDS(infiles, outfile):

# 	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
# 	coding_cutoff = 0.364 if 'human' in outfile else 0.44

# 	# Run
# 	run_r_job('add_cds', infiles, outfile, additional_params=coding_cutoff, conda_env='env', W='03:00', GB=15, n=1, run_locally=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds
# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds.log | jsc
# # ls arion/isoseq/s09-cpat.dir/*/gtf/split/*.cds.err | lr

# #############################################
# ########## 5. Merge
# #############################################

# @collate(addCDS,
# 		 regex(r'(.*)/split/(.*)_.*.cds.gtf'),
# 		 r'\1/\2.cds.gtf')

# def mergeGTF(infiles, outfile):

# 	# Run
# 	run_r_job('merge_gtf', infiles, outfile, conda_env='env', W='00:45', GB=15, n=1, run_locally=False)

# #######################################################
# #######################################################
# ########## S10. Pfam
# #######################################################
# #######################################################

# #############################################
# ########## 1. Translate
# #############################################

# @transform(runCPAT,
# 		   regex(r'(.*)/(s09-cpat.dir)/(.*)/(.*).ORF_prob.best.tsv'),
# 		   add_inputs(r'\1/\2/\3/\4.ORF_seqs.fa'),
# 		   r'\1/s10-pfam.dir/\3/fasta/\3-translated.fasta')

# def translateORFs(infiles, outfile):

# 	# Run
# 	run_r_job('translate_orfs', infiles, outfile, conda_env='env', W="01:00", GB=15, n=1, stdout=outfile.replace('.fasta', '.log'), wait=True)

# #############################################
# ########## 2. Split
# #############################################

# @subdivide(translateORFs,
# 		   regex(r'(.*).fasta'),
# 		   r'\1.fasta.*',
# 		   r'\1.fasta.100')

# def splitORFs(infile, outfiles, outfileRoot):

# 	# Get number
# 	N = outfileRoot.split('.')[-1]

# 	# Command
# 	cmd_str = ''' gt splitfasta -numfiles {N} {infile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitORF_'+outfileRoot.split('/')[-3], wait=True)

# #############################################
# ########## 3. Run
# #############################################

# @transform(splitORFs,
# 		   regex(r'(.*)/fasta/(.*).fasta\.(.*)'),
# 		   add_inputs('/sc/arion/work/torred23/libraries/PfamScan'),
# 		   r'\1/split/\2_\3_pfam.txt')

# def runPfamScan(infiles, outfile):

# 	# Data directory
# 	input_fasta, pfam_dir = infiles

# 	# Command
# 	cmd_str = ''' {pfam_dir}/pfam_scan.pl \
# 		-dir {pfam_dir}/data \
# 		-cpu 50 \
# 		-fasta {input_fasta} > {outfile}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', modules=['hmmer/3.3'], W='06:00', GB=2, n=10, print_cmd=False, run_locally=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

# # find /hpc/users/torred23/pipelines/projects/early-embryo/arion/isoseq/s10-pfam.dir -name "*pfam.log" | jsc

# #############################################
# ########## 3. Format
# #############################################

# @collate(runPfamScan,
# 		 regex(r'(.*)/split/(.*)_.*_pfam.txt'),
# 		 r'\1/\2_pfam.tsv')

# def mergePfamResults(infiles, outfile):

# 	# Initialize
# 	results = []

# 	# Loop
# 	for infile in infiles:

# 		# Read
# 		pfam_dataframe = pd.read_csv(infile, comment='#', delim_whitespace=True, header=None)

# 		# Get column names
# 		colnames = [x.replace(' ', '_').replace('-', '_') for x in pd.read_csv(infile).iloc[26][0][:-1].replace('# ', '').replace('<', '').split('> ')]

# 		# Add column names
# 		pfam_dataframe.columns = colnames

# 		# Fix sequence ID
# 		pfam_dataframe['seq_id'] = [x.split('_ORF')[0] for x in pfam_dataframe['seq_id']]

# 		# Append
# 		results.append(pfam_dataframe)

# 	# Concatenate
# 	result_dataframe = pd.concat(results).query('E_value < 0.1').sort_values('seq_id')

# 	# Write
# 	result_dataframe.to_csv(outfile, index=False, sep='\t')
	
# #######################################################
# #######################################################
# ########## S11. RepeatMasker
# #######################################################
# #######################################################

# #############################################
# ########## 1. Split
# #############################################

# # @follows(runSqantiFinalPass)

# @subdivide('arion/isoseq/s08-sqanti_pass2.dir/*/filtered/*_classification.filtered_lite_no_ism_corrected.fasta',
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)/filtered/(.*).fasta'),
# 		   r'\1/s11-repeatmasker.dir/\2/fasta/*.fasta.*',
# 		   r'\1/s11-repeatmasker.dir/\2/fasta/\3.fasta.100')

# def splitFASTA(infile, outfiles, outfileRoot):

# 	# Get number
# 	N = outfileRoot.split('.')[-1]

# 	# Get temp filename
# 	outdir = os.path.dirname(outfileRoot)
# 	tempfile = os.path.join(outdir, os.path.basename(infile))

# 	# Command
# 	cmd_str = ''' cp {infile} {outdir} && gt splitfasta -numfiles {N} {tempfile} && rm {tempfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitFASTA_'+outfileRoot.split('/')[-3], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'))

# #############################################
# ########## 2. Run
# #############################################

# @transform(splitFASTA,
# 		   regex(r'(.*)/fasta/(.*)'),
# 		   r'\1/split/\2.out')

# def runRepeatMasker(infile, outfile):

# 	# Paths
# 	infile = os.path.join(os.getcwd(), infile)
# 	outdir = os.path.dirname(outfile)

# 	# Species
# 	species = 'Homo sapiens' if 'human' in outfile else 'Mus musculus'

# 	# Command
# 	cmd_str = ''' cd {outdir} && RepeatMasker -species "{species}" -pa 64 -dir . {infile}'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, modules=['repeatmasker/4.1.1'], W='03:00', GB=3, n=10, print_outfile=False, stdout=outfile.replace('.tbl', '.log'), stderr=outfile.replace('.tbl', '.err'))

# #############################################
# ########## 2. Filter
# #############################################

# @collate(runRepeatMasker,
# 		 regex(r'(.*)/split/(.*?)_.*.out'),
# 		 r'\1/\2_repeatmasker.tsv')

# def mergeRepeatMasker(infiles, outfile):

# 	# Run
# 	run_r_job('merge_repeatmasker', infiles, outfile, conda_env='env', run_locally=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# # find arion/isoseq/s11-repeatmasker.dir -name "*.fasta.*.log" | jsc

# #######################################################
# #######################################################
# ########## Summary
# #######################################################
# #######################################################

# #############################################
# ########## 1. Create
# #############################################

# def summaryJobs():
# 	for organism in ['human', 'mouse']:
# 		sqanti_file = 'arion/isoseq/s08-sqanti_pass2.dir/{organism}/filtered/{organism}_classification.filtered_lite_no_ism_classification.txt'.format(**locals())
# 		cpat_file = 'arion/isoseq/s09-cpat.dir/{organism}/{organism}-cpat.ORF_prob.best_results.tsv'.format(**locals())
# 		pfam_file = 'arion/isoseq/s10-pfam.dir/{organism}/{organism}-translated_pfam.tsv'.format(**locals())
# 		repeatmasker_file = 'arion/isoseq/s11-repeatmasker.dir/{organism}/{organism}_repeatmasker.tsv'.format(**locals())
# 		infiles = [sqanti_file, cpat_file, pfam_file, repeatmasker_file]
# 		outfile = 'arion/isoseq/summary.dir/{organism}-isoseq_summary.tsv'.format(**locals())
# 		yield [infiles, outfile]

# @files(summaryJobs)

# def createIsoseqSummary(infiles, outfile):

# 	# Run
# 	run_r_job('create_isoseq_summary', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1)

# #######################################################
# #######################################################
# ########## TALON
# #######################################################
# #######################################################

# #############################################
# ########## 1. minimap2
# #############################################

# @transform('arion/isoseq/s02-alignment.dir/fastq/mouse/*.fastq',
# 		   regex(r'(.*)/fastq/(.*)/(mouse.*).flnc.fastq'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*.primary_assembly.fa'),
# 		   r'arion/isoseq/talon/s1-alignment.dir/\3/\3-minimap2.sam')

# def runMinimap2_v2(infiles, outfile):

# 	# Command
# 	outname = outfile.rsplit('.', 1)[0]
# 	cmd_str = ''' minimap2 -ax splice -uf --secondary=no -C5 -t 30 --MD {infiles[1]} {infiles[0]} | samtools sort -O sam -o {outfile} && \
# 		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f1,5 > {outname}.mapq '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['samtools/1.9', 'minimap2/2.17'], W='03:00', GB=6, n=6, print_outfile=False, stderr=outfile.replace('.sam', '.log'))

# #############################################
# ########## 2. Filter
# #############################################

# @transform(runMinimap2_v2,
# 		   regex(r'(.*).sam'),
# 		   r'\1.filtered.sam')

# def filterSam_v2(infile, outfile):

# 	# Files
# 	outname = outfile.replace('.sam', '')

# 	# Command
# 	cmd_str = ''' sambamba view --with-header --nthreads 30 --sam-input --format sam --filter "not unmapped and mapping_quality >= 50" {infile} | samtools sort -O sam -o {outfile} && \
# 		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f5 -d "	" > {outname}.mapq'''.format(**locals())
# 	# and not duplicate

# 	# Run
# 	run_job(cmd_str, outfile, W='00:15', n=5, GB=1, modules=['samtools/1.9', 'sambamba/0.5.6'], print_cmd=False)

# #############################################
# ########## 2. Reference SJs
# #############################################

# @transform('arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.102.gtf',
# 		   regex(r'(.*)/(.*)(.102).gtf'),
# 		   add_inputs(r'\1/\2.dna_sm.primary_assembly.fa'),
# 		   r'arion/isoseq/talon/s2-transcript_clean.dir/\2\3_SJs.tsv')

# def getReferenceSJs(infile, outfile):

# 	# Command
# 	cmd_str = ''' python ${TranscriptCleanPath}/accessory_scripts/get_SJs_from_gtf.py \
# 		--f {infiles[0]} \
# 		--g {infiles[1]} \
# 		--o {outfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W='01:00', n=1, GB=15, conda_env='talon', print_cmd=False)

# # get SJs and run TranscriptClean

# ############################################
# ######### 3. Initialize database
# ############################################

# @transform('arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.102.gtf',
# 		   regex(r'(.*)/(.*)(.102).gtf'),
# 		   add_inputs(r'\1/\2.dna_sm.primary_assembly.fa'),
# 		   r'arion/isoseq/talon/s2-talon.dir/\2\3_talon.db')

# def initializeTalonDatabase(infiles, outfile):

# 	# Get parameters
# 	annotation_name = '_'.join(infiles[0].split('.')[-3:-1])
# 	genome_name = 'mm10' if 'mouse' in infiles[0] else 'hg38'
# 	outname = outfile[:-len('.db')]

# 	# Command
# 	cmd_str = ''' talon_initialize_database \
# 		--f {infiles[0]} \
# 		--a {annotation_name} \
# 		--g {genome_name} \
# 		--l 0 \
# 		--idprefix TALON \
# 		--5p 500 \
# 		--3p 300 \
# 		--o {outname} '''.format(**locals())

# # 1000bp each
# # -l 200
# # keep monoexons that are in annotation

	# # Run
	# run_job(cmd_str, outfile, W='01:00', n=1, GB=15, conda_env='talon', print_cmd=False, stdout=outfile.replace('.db', '.log'), stderr=outfile.replace('.db', '.err'))

# ############################################
# ######### 4. Run TALON
# ############################################

# @transform(initializeTalonDatabase,
# 		   regex(r'(.*)/(.*).db'),
# 		   add_inputs('arion/isoseq/talon/s2-talon.dir/mouse-samples.csv'),
# 		   r'\1/annotated/\2_annotated_QC.log')

# def annotateTalonDatabase(infiles, outfile):
	
# 	# Get parameters
# 	build = 'mm10' if 'mouse' in infiles[1] else 'hg38'
# 	outname = outfile[:-len('_QC.log')]

# 	# Command
# 	cmd_str = ''' talon \
# 		--db {infiles[0]} \
# 		--f {infiles[1]} \
# 		--build {build} \
# 		--threads 50 \
# 		--cov 0.9 \
# 		--identity 0.9 \
# 		--o {outname} '''.format(**locals())
# # 99 coverage
# # 95 identity
# 	# Run
# 	# run_job(cmd_str, outfile, W='06:00', n=6, GB=10, conda_env='talon', print_cmd=False, stdout=outfile.replace('.db', '.log'), stderr=outfile.replace('.db', '.err'))

# ############################################
# ######### 5. Create GTF
# ############################################

# @follows(annotateTalonDatabase)

# @transform(initializeTalonDatabase,
# 		   suffix('.db'),
# 		   '_annotated_v2.gtf')

# def getTalonGTF(infile, outfile):
	
# 	# Get parameters
# 	build = 'mm10' if 'mouse' in infile else 'mm10'
# 	annotation_name = 'GRCm38_102' if 'mouse' in infile else 'GRCm38_102'

# 	# Command
# 	cmd_str = ''' talon_create_GTF \
# 		--db {infile} \
# 		--build {build} \
# 		--observed \
# 		--annot {annotation_name} \
# 		--o {outfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W='06:00', n=6, GB=10, conda_env='talon', print_cmd=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

#############################################
########## 2. Get abundance
#############################################

# mouse_isoforms = glob.glob('')
# human_isoforms = glob.glob('')


	# elif parameters['method'] == 'tama':

	# 	# Command
	# 	outname = outfile[:-len('.bed')]
	# 	cmd_str = ''' python $TAMA_HOME/tama_collapse.py \
	# 		-s {sam} \
	# 		-f {genome_fasta} \
	# 		-x no_cap \
	# 		-a {max_5_diff} \
	# 		-m {max_fuzzy_junction} \
	# 		-z {max_3_diff} \
	# 		-c 99 \
	# 		-i 95 \
	# 		-rm low_mem \
	# 		-p {basename} && \
	# 		python $TAMA_HOME/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {basename}.bed {basename}.gtf
	# 	'''.format(**locals(), **parameters)

	# 	# Run
	# 	run_job(cmd_str, outfile, modules=['python/2.7.17'], W='01:00', GB=30, n=1, print_outfile=False, stdout=outfile.rsplit('.', 1)[0]+'.log', stderr=outfile.rsplit('.', 1)[0]+'.err', q='premium')

# #############################################
# ########## 1. Split SAM
# #############################################

# @subdivide(filterSam,
# 		   regex(r'(.*)/s02-alignment.dir/minimap2/(.*)/(.*)/.*.sam'),
# 		   r'\1/s03-collapsed.dir/\2/split/\3/\3.split.*.sam',
# 		   r'\1/s03-collapsed.dir/\2/split/\3/\3.split.00.sam')

# def splitSam(infile, outfiles, outfileRoot):

# 	# Parameters
# 	header_file = os.path.join(os.path.dirname(outfileRoot), os.path.basename(infile).replace('.sam', '_header.txt'))
# 	split_basename = outfileRoot[:-len('00.sam')]

# 	# Command
# 	cmd_str = ''' samtools view -H {infile} > {header_file} && samtools view {infile} | split - {split_basename} --numeric-suffixes -l 50000 --filter='cat {header_file} - > $FILE.sam' && rm {header_file} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['samtools/1.9'], W='00:30', GB=25, n=1, print_outfile=False, print_cmd=False)

# #############################################
# ########## 2. Collapse isoforms
# #############################################

# def collapseJobs():
	
# 	# Get SAM
# 	# sam_files = glob.glob('arion/isoseq/s02-alignment.dir/minimap2/mouse/*1C/*-minimap2.filtered.sam')
# 	# sam_files = glob.glob('arion/isoseq/s02-alignment.dir/minimap2/*/*/*-minimap2.filtered.sam')
# 	# sam_files = glob.glob('arion/isoseq/s03-collapsed.dir/*/split/*/*.split.*.sam')
# 	sam_files = glob.glob('arion/isoseq/s03-collapsed.dir/*/split/human_2C_EMB2/*.split.*.sam')

# 	# Parameters
# 	parameter_dict = {
# 		'A': {'max_fuzzy_junction': 0, 'max_5_diff': 10, 'max_3_diff': 10},
# 		'B': {'max_fuzzy_junction': 5, 'max_5_diff': 10, 'max_3_diff': 10},
# 		'C': {'max_fuzzy_junction': 10, 'max_5_diff': 10, 'max_3_diff': 10},
# 		'D': {'max_fuzzy_junction': 5, 'max_5_diff': 100, 'max_3_diff': 10},
# 		'E': {'max_fuzzy_junction': 5, 'max_5_diff': 10, 'max_3_diff': 1000},
# 		'F': {'max_fuzzy_junction': 5, 'max_5_diff': 100, 'max_3_diff': 1000},
# 	}

# 	# Write
# 	with open('arion/isoseq/s03-collapsed.dir/collapse_parameters.json', 'w') as openfile:
# 		json.dump(parameter_dict, openfile)

# 	# Loop
# 	for infile in sam_files:
# 		split_name = os.path.basename(infile)[:-len('.sam')]
# 		sample_name = split_name.split('.split')[0]
# 		organism = infile.split('/')[3]
# 		for letter, parameters in parameter_dict.items():
# 			for method in ['cdna_cupcake', 'tama']:
# 				parameters = parameters.copy()
# 				parameters['method'] = method
# 				parameters['basename'] = 'arion/isoseq/s03-collapsed.dir/{organism}/collapsed/{sample_name}/{split_name}-{method}-{letter}/{split_name}-{method}-{letter}'.format(**locals())
# 				parameters['fastq'] = 'arion/isoseq/s02-alignment.dir/fastq/{organism}/{sample_name}.flnc.fastq'.format(**locals())
# 				parameters['genome_fasta'] = glob.glob('arion/datasets/reference_genomes/{organism}/*.primary_assembly.fa'.format(**locals()))[0]
# 				outfile = parameters['basename']+('.gtf' if method == 'tama' else '.collapsed.gff')
# 				yield [infile, outfile, parameters]

# # @follows(filterSam)

# @files(collapseJobs)

# def collapseIsoforms(infile, outfile, parameters):

# 	# SAM file
# 	sam = infile

# 	# cDNA Cupcake
# 	if parameters['method'] == 'cdna_cupcake':

# 		# Collapse
# 		cmd_str = ''' collapse_isoforms_by_sam.py \
# 			--input {fastq} --fq \
# 			--sam {sam} \
# 			--prefix {basename} \
# 			--max_5_diff {max_5_diff} \
# 			--max_fuzzy_junction {max_fuzzy_junction} \
# 			--max_3_diff {max_3_diff} \
# 			--dun-merge-5-shorter \
# 			--min-coverage 0.99 \
# 			--min-identity 0.95 \
# 		'''.format(**locals(), **parameters)

# 		# Run
# 		# run_job(cmd_str, outfile, conda_env='cdna_cupcake', W='00:15', GB=15, n=1, print_outfile=False, stdout=outfile.rsplit('.', 1)[0]+'.log', stderr=outfile.rsplit('.', 1)[0]+'.err')

# 	elif parameters['method'] == 'tama':

# 		# Command
# 		outname = outfile[:-len('.bed')]
# 		cmd_str = ''' python $TAMA_HOME/tama_collapse.py \
# 			-s {sam} \
# 			-f {genome_fasta} \
# 			-x no_cap \
# 			-a {max_5_diff} \
# 			-m {max_fuzzy_junction} \
# 			-z {max_3_diff} \
# 			-c 99 \
# 			-i 95 \
# 			-rm low_mem \
# 			-p {basename} && \
# 			python $TAMA_HOME/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {basename}.bed {basename}.gtf
# 		'''.format(**locals(), **parameters)

# 		# Run
# 		run_job(cmd_str, outfile, modules=['python/2.7.17'], W='01:00', GB=30, n=1, print_outfile=False, stdout=outfile.rsplit('.', 1)[0]+'.log', stderr=outfile.rsplit('.', 1)[0]+'.err', q='premium')


#######################################################
#######################################################
########## S3. TAMA
#######################################################
#######################################################

# #############################################
# ########## 1. Split SAM
# #############################################

# @subdivide(filterSam,
# 		   regex(r'(.*)/s02-alignment.dir/minimap2/(.*)/(.*)/.*.sam'),
# 		   r'\1/s03-tama.dir/\2/split/\3/\3.split.*.sam',
# 		   r'\1/s03-tama.dir/\2/split/\3/\3.split.00.sam')

# def splitSam(infile, outfiles, outfileRoot):

# 	# Parameters
# 	header_file = os.path.join(os.path.dirname(outfileRoot), os.path.basename(infile).replace('.sam', '_header.txt'))
# 	split_basename = outfileRoot[:-len('00.sam')]

# 	# Command
# 	cmd_str = ''' samtools view -H {infile} > {header_file} && samtools view {infile} | split - {split_basename} --numeric-suffixes -l 50000 --filter='cat {header_file} - > $FILE.sam' && rm {header_file} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['samtools/1.9'], W='00:30', GB=25, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 2. Collapse
#############################################

# @transform('*.sam',
# # @transform(splitSam,
# 		   regex(r'(.*.dir)/(.*)/split/(.*)/(.*).sam'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*.primary_assembly.fa'),
# 		   r'\1/\2/collapse/\3/\4/\4.tama.bed')

# def tamaCollapse(infiles, outfile):

# 	# Command
# 	outname = outfile[:-len('.bed')]
# 	cmd_str = ''' python $TAMA_HOME/tama_collapse.py \
# 		-s {infiles[0]} \
# 		-f {infiles[1]} \
# 		-x no_cap -i 95 -z 1000 -rm low_mem \
# 		-p {outname} && bedtools getfasta -name -fi {infiles[1]} -bed {outfile} -fo {outname}.fasta
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['python/2.7.17', 'BEDTools/2.29.0'], W='10:00', GB=30, n=1, print_outfile=False, stdout=outfile.replace('.bed', '.log'))

# jj | grep tama | grep RUN | wc -l; jj | grep tama | grep PEND | wc -l
# ls arion/isoseq/s03-tama.dir/*/collapse/*/*/*.tama.fasta | wc -l
# ls arion/isoseq/s03-tama.dir/*/collapse/*/*/*.tama.bed | wc -l
# ls arion/isoseq/s03-tama.dir/*/split/*/*.sam | wc -l
# for i in $(ls arion/isoseq/s03-tama.dir/human/collapse/*/*/*.tama.log); do seconds=$(cat $i | grep 'Run time'  | cut -f2 -d : | sed 's/ //g' | sed 's/sec.//g'); status=$(cat $i | sed -n '35,35p'); echo $(basename $i | sed 's/.tama.log//g') $seconds $status; done > arion/isoseq/s03-tama.dir/job_status.txt; cat arion/isoseq/s03-tama.dir/job_status.txt | sort -k2,2n
# cat arion/isoseq/s03-tama.dir/job_status.txt | grep -v completed | cut -f1 -d ' ' > arion/isoseq/s03-tama.dir/failed_jobs.txt

# This script collapses mapped transcript models

# optional arguments:
#   -h, --help  show this help message and exit
#   -s S        Sorted sam file (required)
#   -f F        Genome fasta file (required)
#   -p P        Output prefix (required)
#   -x X        Capped flag: capped or no_cap
#   -e E        Collapse exon ends flag: common_ends or longest_ends (default
#               common_ends)
#   -c C        Coverage (default 99)
#   -i I        Identity (default 85)
#   -icm ICM    Identity calculation method (default ident_cov for including
#               coverage) (alternate is ident_map for excluding hard and soft
#               clipping)
#   -a A        5 prime threshold (default 10)
#   -m M        Exon/Splice junction threshold (default 10)
#   -z Z        3 prime threshold (default 10)
#   -d D        Flag for merging duplicate transcript groups (default is
#               merge_dup will merge duplicates ,no_merge quits when duplicates
#               are found)
#   -sj SJ      Use error threshold to prioritize the use of splice junction
#               information from collapsing transcripts(default no_priority,
#               activate with sj_priority)
#   -sjt SJT    Threshold for detecting errors near splice junctions (default is
#               10bp)
#   -lde LDE    Threshold for amount of local density error near splice
#               junctions that is allowed (default is 1000 errors which
#               practically means no threshold is applied)
#   -ses SES    Simple error symbol. Use this to pick the symbol used to
#               represent matches in the simple error string for LDE output.
#   -b B        Use BAM instead of SAM
#   -log LOG    Turns on/off output of collapsing process. (default on, use
#               log_off to turn off)
#   -v V        Prints out version date and exits.
#   -rm RM      Run mode allows you to use original or low_mem mode, default is
#               original
#   -vc VC      Variation covwerage threshold: Default 5 reads

#############################################
########## 3. Config
#############################################

# # @collate('arion/isoseq/s03-tama.dir/mouse/collapse/*/*.tama.bed',
# @collate(tamaCollapse,
# 		 regex(r'(.*)/(.*)/collapse/.*.bed'),
# 		 r'\1/\2/merge/\2-filelist.txt')

# def tamaConfig(infiles, outfile):

# 	# Get dict
# 	file_dict = [{'file_name': x, 'cap_flag': 'no_cap', 'merge_priority': '1,1,1', 'source_name': x.split('/')[-2]} for x in infiles]

# 	# Convert to dataframe
# 	file_dataframe = pd.DataFrame(file_dict)

# 	# Create directory
# 	outdir = os.path.dirname(outfile)
# 	if not os.path.exists(outdir):
# 		os.makedirs(outdir)

# 	# Write
# 	file_dataframe.to_csv(outfile, sep='\t', header=None, index=False)

# #############################################
# ########## 4. Merge
# #############################################

# @transform(tamaConfig,
# 		   suffix('filelist.txt'),
# 		   'tama_merged.gtf')

# def tamaMerge(infile, outfile):

# 	# Command
# 	outname = outfile[:-len('.gtf')]
# 	cmd_str = ''' python $TAMA_HOME/tama_merge.py -f {infile} -d merge_dup -p {outname} && \
# 		python $TAMA_HOME/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {outname}.bed {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['python/2.7.17'], W='01:00', GB=10, n=1, print_cmd=False, stdout=outfile.replace('.gtf', '.log'))

# # usage: tama_merge.py [-h] [-f F] [-p P] [-e E] [-a A] [-m M] [-z Z] [-d D]
# #                  [-s S] [-cds CDS]

# # This script merges transcriptomes.

# # optional arguments:
# #   -h, --help  show this help message and exit
# #   -f F        File list
# #   -p P        Output prefix
# #   -e E        Collapse exon ends flag: common_ends or longest_ends (Default is
# #               common_ends)
# #   -a A        5 prime threshold (Default is 10)
# #   -m M        Exon ends threshold/ splice junction threshold (Default is 10)
# #   -z Z        3 prime threshold (Default is 10)
# #   -d D        Flag for merging duplicate transcript groups (default no_merge
# #               quits when duplicates are found, merge_dup will merge
# #               duplicates)
# #   -s S        Use gene and transcript ID from a merge source. Specify source
# #               name from filelist file here.
# #   -cds CDS    Use CDS from a merge source. Specify source name from filelist
# #               file here.

# #############################################
# ########## 5. Get FL counts
# #############################################

# @follows(tamaCollapse)

# @transform('arion/isoseq/s03-tama.dir/*/collapse/*/*/*tama_trans_report.txt',
# 		   suffix('_trans_report.txt'),
# 		   '.abundance.txt')

# def getFlCounts(infile, outfile):
	
# 	# Read
# 	report_dataframe = pd.read_table(infile)

# 	# Get counts
# 	report_dataframe['norm_fl'] = ['{:.4e}'.format(x) for x in report_dataframe['num_clusters']/report_dataframe['num_clusters'].sum()]

# 	# Subset
# 	report_dataframe = report_dataframe.rename(columns={'transcript_id': 'pbid', 'num_clusters': 'count_fl'})[['pbid', 'count_fl', 'norm_fl']]
	
# 	# Get values
# 	report_dataframe_string = report_dataframe.to_csv(sep='\t', index=False)
# 	fl_total = report_dataframe['count_fl'].sum()

# 	# String
# 	output_str = '''#
# # -----------------
# # Field explanation
# # -----------------
# # count_fl: Number of associated FL reads
# # norm_fl: count_fl / total number of FL reads, mapped or unmapped
# # Total Number of FL reads:  {fl_total} 
# #
# {report_dataframe_string}'''.format(**locals())

# 	# Write
# 	with open(outfile, 'w') as openfile:
# 		openfile.write(output_str)

#######################################################
#######################################################
########## S3. cDNA Cupcake
#######################################################
#######################################################

#############################################
########## 1. Regular (cupcake)
#############################################

# @transform(mergeIsoseqAlignments,
# 		   regex(r'(.*)/(s3-merged.dir)/(.*)/(.*).sorted.uniq'),
# 		   add_inputs(r'\1/\2/\3/*.fastq', r'\1/\2/\3/*.cluster.csv'),
# 		   r'\1/s4-transcripts.dir/\3/regular/\4.collapsed.min_fl_2.filtered.gff')

# def cupcakeCollapse(infiles, outfile):

# 	# Split infiles
# 	sam_file, fastq_file, cluster_file = infiles

# 	# Variables
# 	min_fl = outfile[:-len('.filtered.gff')].split('min_fl_')[-1]
# 	basename = outfile[:-len('.collapsed.min_fl_2.filtered.gff')]

# 	# Command (collapse and get abundance)
# 	cmd_str = '''
# 		collapse_isoforms_by_sam.py --input {fastq_file} --fq -s {sam_file} --dun-merge-5-shorter -o {basename} && \
# 		get_abundance_post_collapse.py {basename}.collapsed {cluster_file} && \
# 		filter_by_count.py --min_count {min_fl} --dun_use_group_count {basename}.collapsed && \
# 		filter_away_subset.py {basename}.collapsed.min_fl_{min_fl} && \
# 		get_seq_stats.py {basename}.collapsed.min_fl_{min_fl}.filtered.rep.fq
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3.env', W='00:30', GB=10, n=1, print_outfile=False, stdout=outfile.replace('.gff', '.log'))

# usage:	 collapse_isoforms_by_sam.py [-h] [--input INPUT] [--fq] -s SAM -o
#                                    PREFIX [-c MIN_ALN_COVERAGE]
#                                    [-i MIN_ALN_IDENTITY]
#                                    [--max_fuzzy_junction MAX_FUZZY_JUNCTION]
#                                    [--max_5_diff MAX_5_DIFF]
#                                    [--max_3_diff MAX_3_DIFF]
#                                    [--flnc_coverage FLNC_COVERAGE]
#                                    [--dun-merge-5-shorter]

# optional arguments:
#   -h, --help            show this help message and exit
#   --input INPUT         Input FA/FQ filename
#   --fq                  Input is a fastq file (default is fasta)
#   -s SAM, --sam SAM     Sorted GMAP SAM filename
#   -o PREFIX, --prefix PREFIX
#                         Output filename prefix
#   -c MIN_ALN_COVERAGE, --min-coverage MIN_ALN_COVERAGE
#                         Minimum alignment coverage (default: 0.99)
#   -i MIN_ALN_IDENTITY, --min-identity MIN_ALN_IDENTITY
#                         Minimum alignment identity (default: 0.95)
#   --max_fuzzy_junction MAX_FUZZY_JUNCTION
#                         Max fuzzy junction dist (default: 5 bp)
#   --max_5_diff MAX_5_DIFF
#                         Maximum allowed 5' difference if on same exon
#                         (default: 1000 bp)
#   --max_3_diff MAX_3_DIFF
#                         Maximum allowed 3' difference if on same exon
#                         (default: 100 bp)
#   --flnc_coverage FLNC_COVERAGE
#                         Minimum # of FLNC reads, only use this for aligned
#                         FLNC reads, otherwise results undefined!
#   --dun-merge-5-shorter
#                         Don't collapse shorter 5' transcripts (default: turned
#                         off)
#######################################################
#######################################################
########## IS4. SQANTI First Pass
#######################################################
#######################################################

#############################################
########## 1. SQANTI 1
#############################################

# @transform(tamaMerge,
# 		   regex(r'(.*)/(.*)/s03-tama.dir/.*/(.*)-tama_merged.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa'),
# 		   r'\1/\2/s04-sqanti_pass1.dir/\3/\3_corrected.gtf')

# def runSqantiPass1(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if isinstance(x, str) and x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta = infiles

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_corrected.gtf')]

# 	# PYTHONPATH
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/cDNA_Cupcake:/sc/arion/work/torred23/libraries/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PYTHONPATH={PYTHONPATH} && cd {outdir} && sqanti3_qc.py -v && sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 		--gtf \
# 		--force_id_ignore \
# 		--dir . \
# 		--output {output}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3.env', W='02:00', GB=5, n=3, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.log'))
	
# #######################################################
# #######################################################
# ########## IS5. Alignment indices
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR
# #############################################

# @transform(runSqantiPass1,
# 		   regex(r'(.*)/(isoseq)/.*/(.*)/.*.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa'),
# 		   r'\1/\2/s05-indices.dir/\3/STAR/')

# def buildStarIndex(infiles, outfile):

# 	# Split
# 	reference_gtf, reference_fasta = infiles

# 	# Command
# 	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {reference_fasta} --sjdbGTFfile {reference_gtf} --runThreadN 32 --outFileNamePrefix {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='06:00', GB=10, n=4, ow=True, print_cmd=False, jobname='-'.join(outfile.split('/')[-3:-1]), stdout=os.path.join(outfile, 'job.log'))

# #############################################
# ########## 2. Salmon
# #############################################

# @follows(runSqantiPass1)

# @transform('arion/isoseq/s04-sqanti_pass1.dir/*/*_corrected.fasta',
# 		   regex(r'(.*)/s04-sqanti_pass1.dir/(.*)/.*.fasta'),
# 		   r'\1/s05-indices.dir/\2/salmon/')

# def buildSalmonIndex(infile, outfile):

# 	# Command
# 	cmd_str = '''salmon index -t {infile} -i {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['salmon/1.2.1'], W='06:00', GB=8, n=4, print_cmd=False, ow=True, jobname='-'.join(outfile.split('/')[-3:-1]), stdout=os.path.join(outfile, 'job.log'))

# #######################################################
# #######################################################
# ########## S6. Illumina alignment
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR pass 1
# #############################################

# @follows(buildStarIndex)

# @collate(illumina_fastq,
# 		 regex(r'.*/(.*)/(.*)/.*.fastq.gz'),
# 		 add_inputs(r'arion/isoseq/s05-indices.dir/\1/STAR'),
# 		 r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass1/\2/\2-SJ.out.tab')

# def getSJsPass1(infiles, outfile):

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[0][1]} \
# 		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 48 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", GB=5, n=8, modules=['star/2.7.5b'], print_outfile=False, stdout=outfile.replace('-SJ.out.tab', '.log'))
	
# #############################################
# ########## 2. STAR pass 2
# #############################################

# @follows(getSJsPass1)

# @collate(illumina_fastq,
# 		 regex(r'.*/(.*)/(.*)/.*.fastq.gz'),
# 		 add_inputs(r'arion/isoseq/s05-indices.dir/\1/STAR', r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass1/*/*-SJ.out.tab'),
# 		 r'arion/isoseq/s07-illumina_alignment.dir/\1/STAR/pass2/\2/\2-SJ.out.tab')

# def getSJsPass2(infiles, outfile):

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]
# 	sj_files_str = ' '.join(infiles[0][2:])

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[0][1]} \
# 		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 48 \
# 		--sjdbFileChrStartEnd {sj_files_str} \
# 		--limitSjdbInsertNsj 5000000 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	# run_job(cmd_str, outfile, W="09:00", GB=10, n=8, modules=['star/2.7.5b'], print_outfile=True, q='express')

# #############################################
# ########## 3. Salmon
# #############################################

# @follows(buildSalmonIndex)

# @collate(illumina_fastq,
# 		 regex(r'.*/(.*)/(.*)/.*.fastq.gz'),
# 		 add_inputs(r'arion/isoseq/s05-indices.dir/\1/salmon'),
# 		 r'arion/isoseq/s07-illumina_alignment.dir/\1/salmon/samples/\2/quant.sf')

# def runSalmon(infiles, outfile):

# 	# Split
# 	fastq_files, salmon_index = infiles

# 	# Prefix
# 	outdir = os.path.dirname(outfile)

# 	# Command
# 	cmd_str = '''salmon quant -i {infiles[0][1]} -l A -1 {infiles[0][0]} -2 {infiles[1][0]} --validateMappings -p 100 -o {outdir}'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, modules=['salmon/1.2.1'], W='10:00', GB=3, n=10, print_outfile=False, jobname=outfile.split('dir/')[-1].replace('/quant.sf', '').replace('/', '_'), stdout=outfile.replace('quant.sf', 'job.log'), q='express')

# # ls arion/isoseq/s07-illumina_alignment.dir/*/salmon/samples/*/job.log | wc -l
# # ls arion/isoseq/s07-illumina_alignment.dir/*/salmon/samples/*/job.log | js

# #############################################
# ########## 4. Aggregate counts
# #############################################

# # @collate(runSalmon,
# @collate('arion/isoseq/s07-illumina_alignment.dir/*/salmon/samples/*/quant.sf',
# 		 regex(r'(.*.dir)/(.*)/salmon/samples/.*.quant.sf'),
# 		 r'\1/\2/salmon/\2-transcript_tpm.tsv')

# def aggregateSalmonCounts(infiles, outfile):

# 	# Run
# 	run_r_job('aggregate_salmon_counts', infiles, outfile, conda_env='env', W='00:30', n=1, GB=10, run_locally=True, print_outfile=False)

# #######################################################
# #######################################################
# ########## IS7. SQANTI Second Pass
# #######################################################
# #######################################################

# #############################################
# ########## 1. QC
# #############################################

# # @follows(aggregateSalmonCounts, getSJsPass2)

# # @transform(runSqantiPass1,
# @transform('arion/isoseq/s04-sqanti_pass1.dir/*/*_corrected.gtf',
# 		   regex(r'(.*)/(.*)/s04-sqanti_pass1.dir/(.*)/.*.gtf'),
# 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.gtf', r'\1/datasets/reference_genomes/\3/*.primary_assembly.fa', r'arion/isoseq/s07-illumina_alignment.dir/\3/salmon/\3-transcript_tpm.tsv', r'arion/isoseq/s07-illumina_alignment.dir/\3/STAR/pass2/*/*-SJ.out.tab'),
# 		   r'\1/\2/s08-sqanti_pass2.dir/\3/\3_corrected.gtf')

# def runSqantiPass2(infiles, outfile):
	
# 	# Add paths
# 	infiles = [os.path.join(os.getcwd(), x) if x.startswith('arion') else x for x in infiles]
# 	transcript_gtf, ensembl_gtf, ensembl_fasta, salmon_tpm = infiles[:4]
# 	sj_files = ','.join(infiles[4:]) # same result, even though in params.txt file only last file is shown

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_corrected.gtf')]

# 	# PYTHONPATH
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/cDNA_Cupcake:/sc/arion/work/torred23/libraries/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PYTHONPATH={PYTHONPATH} && cd {outdir} && sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
# 		--gtf \
# 		--force_id_ignore \
# 		--expression {salmon_tpm} \
# 		--coverage "{sj_files}" \
# 		--dir . \
# 		--output {output}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3.env', W='06:00', GB=15, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.std'))
	
# # usage: sqanti3_qc.py [-h] [--min_ref_len MIN_REF_LEN] [--force_id_ignore]
# #                      [--aligner_choice {minimap2,deSALT,gmap}]
# #                      [--cage_peak CAGE_PEAK]
# #                      [--polyA_motif_list POLYA_MOTIF_LIST]
# #                      [--polyA_peak POLYA_PEAK] [--phyloP_bed PHYLOP_BED]
# #                      [--skipORF] [--is_fusion] [-g] [-e EXPRESSION]
# #                      [-x GMAP_INDEX] [-t CPUS] [-n CHUNKS] [-o OUTPUT]
# #                      [-d DIR] [-c COVERAGE] [-s SITES] [-w WINDOW]
# #                      [--genename] [-fl FL_COUNT] [-v] [--isoAnnotLite]
# #                      [--gff3 GFF3]
# #                      isoforms annotation genome

# # Structural and Quality Annotation of Novel Transcript Isoforms

# # positional arguments:
# #   isoforms              Isoforms (FASTA/FASTQ or gtf format; By default
# #                         "FASTA/FASTQ". For GTF, use --gtf
# #   annotation            Reference annotation file (GTF format)
# #   genome                Reference genome (Fasta format)

# # optional arguments:
# #   -h, --help            show this help message and exit
# #   --min_ref_len MIN_REF_LEN
# #                         Minimum reference transcript length (default: 200 bp)
# #   --force_id_ignore     Allow the usage of transcript IDs non related with
# #                         PacBio's nomenclature (PB.X.Y)
# #   --aligner_choice {minimap2,deSALT,gmap}
# #   --cage_peak CAGE_PEAK
# #                         FANTOM5 Cage Peak (BED format, optional)
# #   --polyA_motif_list POLYA_MOTIF_LIST
# #                         Ranked list of polyA motifs (text, optional)
# #   --polyA_peak POLYA_PEAK
# #                         PolyA Peak (BED format, optional)
# #   --phyloP_bed PHYLOP_BED
# #                         PhyloP BED for conservation score (BED, optional)
# #   --skipORF             Skip ORF prediction (to save time)
# #   --is_fusion           Input are fusion isoforms, must supply GTF as input
# #                         using --gtf
# #   -g, --gtf             Use when running SQANTI by using as input a gtf of
# #                         isoforms
# #   -e EXPRESSION, --expression EXPRESSION
# #                         Expression matrix (supported: Kallisto tsv)
# #   -x GMAP_INDEX, --gmap_index GMAP_INDEX
# #                         Path and prefix of the reference index created by
# #                         gmap_build. Mandatory if using GMAP unless -g option
# #                         is specified.
# #   -t CPUS, --cpus CPUS  Number of threads used during alignment by aligners.
# #                         (default: 10)
# #   -n CHUNKS, --chunks CHUNKS
# #                         Number of chunks to split SQANTI3 analysis in for
# #                         speed up (default: 1).
# #   -o OUTPUT, --output OUTPUT
# #                         Prefix for output files.
# #   -d DIR, --dir DIR     Directory for output files. Default: Directory where
# #                         the script was run.
# #   -c COVERAGE, --coverage COVERAGE
# #                         Junction coverage files (provide a single file or a
# #                         file pattern, ex: "mydir/*.junctions").
# #   -s SITES, --sites SITES
# #                         Set of splice sites to be considered as canonical
# #                         (comma-separated list of splice sites). Default:
# #                         GTAG,GCAG,ATAC.
# #   -w WINDOW, --window WINDOW
# #                         Size of the window in the genomic DNA screened for
# #                         Adenine content downstream of TTS
# #   --genename            Use gene_name tag from GTF to define genes. Default:
# #                         gene_id used to define genes
# #   -fl FL_COUNT, --fl_count FL_COUNT
# #                         Full-length PacBio abundance file
# #   -v, --version         Display program version number.
# #   --isoAnnotLite        Run isoAnnot Lite to output a tappAS-compatible gff3
# #                         file
# #   --gff3 GFF3           Precomputed tappAS species specific GFF3 file. It will
# #                         serve as reference to transfer functional attributes

# #############################################
# ########## 2. Filter
# #############################################

# @transform(runSqantiPass2,
# 		   regex(r'(.*)/(.*)_corrected.gtf'),
# 		   add_inputs(r'\1/\2_corrected.fasta', r'\1/\2_classification.txt'),
# 		   r'\1/\2_classification.filtered_lite.fasta')

# def filterSqanti(infiles, outfile):
	
# 	# Get infiles
# 	infiles = [os.path.join(os.getcwd(), x) if isinstance(x, str) and x.startswith('arion') else x for x in infiles]
# 	gtf_file, fasta_file, classification_file = infiles

# 	# Get output
# 	outdir = os.path.dirname(outfile)
# 	output = os.path.basename(outfile)[:-len('_classification.filtered_lite.gtf')]

# 	# PYTHONPATH
# 	PYTHONPATH = '/sc/arion/work/torred23/libraries/cDNA_Cupcake:/sc/arion/work/torred23/libraries/cDNA_Cupcake/sequence'

# 	# Command
# 	cmd_str = ''' export PYTHONPATH={PYTHONPATH} && cd {outdir} && sqanti3_RulesFilter.py {classification_file} {fasta_file} {gtf_file} \
# 		--intrapriming 0.6 \
# 		--runAlength 6 \
# 		--max_dist_to_known_end 50 \
# 		--min_cov 3
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3.env', W='00:30', GB=10, n=1, print_cmd=False, stderr=outfile.replace('.gtf', '.err'), stdout=outfile.replace('.gtf', '.std'))
	
# # usage: sqanti3_RulesFilter.py [-h] [--sam SAM] [--faa FAA] [-a INTRAPRIMING]
# #                               [-r RUNALENGTH] [-m MAX_DIST_TO_KNOWN_END]
# #                               [-c MIN_COV] [--filter_mono_exonic] [--skipGTF]
# #                               [--skipFaFq] [--skipJunction] [-v]
# #                               sqanti_class isoforms gtf_file

# # Filtering of Isoforms based on SQANTI3 attributes

# # positional arguments:
# #   sqanti_class          SQANTI classification output file.
# #   isoforms              fasta/fastq isoform file to be filtered by SQANTI3
# #   gtf_file              GTF of the input fasta/fastq

# # optional arguments:
# #   -h, --help            show this help message and exit
# #   --sam SAM             (Optional) SAM alignment of the input fasta/fastq
# #   --faa FAA             (Optional) ORF prediction faa file to be filtered by
# #                         SQANTI3
# #   -a INTRAPRIMING, --intrapriming INTRAPRIMING
# #                         Adenine percentage at genomic 3' end to flag an
# #                         isoform as intra-priming (default: 0.6)
# #   -r RUNALENGTH, --runAlength RUNALENGTH
# #                         Continuous run-A length at genomic 3' end to flag an
# #                         isoform as intra-priming (default: 6)
# #   -m MAX_DIST_TO_KNOWN_END, --max_dist_to_known_end MAX_DIST_TO_KNOWN_END
# #                         Maximum distance to an annotated 3' end to preserve as
# #                         a valid 3' end and not filter out (default: 50bp)
# #   -c MIN_COV, --min_cov MIN_COV
# #                         Minimum junction coverage for each isoform (only used
# #                         if min_cov field is not 'NA'), default: 3
# #   --filter_mono_exonic  Filter out all mono-exonic transcripts (default: OFF)
# #   --skipGTF             Skip output of GTF
# #   --skipFaFq            Skip output of isoform fasta/fastq
# #   --skipJunction        Skip output of junctions file
# #   -v, --version         Display program version number.

# #######################################################
# #######################################################
# ########## S8. CPAT
# #######################################################
# #######################################################

# # split fasta

# #############################################
# ########## 1. Run
# #############################################

# @transform(filterSqanti,
# 		   regex(r'(.*)/s08-sqanti_pass2.dir/(.*)_classification.filtered_lite.fasta'),
# 		   add_inputs('/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/Human_Hexamer.tsv', '/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/Human_logitModel.RData'),
# 		   r'\1/s08-cpat.dir/\2-cpat.ORF_prob.best.tsv')

# def runCPAT(infiles, outfile):
	
# 	# Split
# 	transcript_fasta, cpat_hexamer, cpat_model = infiles
# 	basename = outfile.replace('.ORF_prob.best.tsv', '')

# 	# Command
# 	cmd_str = ''' ~/.conda/envs/env/bin/cpat.py -g {transcript_fasta} \
# 		-x {cpat_hexamer} \
# 		-d {cpat_model} \
# 		--top-orf=5 \
# 		-o {basename}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='00:30', GB=10, n=1, run_locally=False, print_outfile=False)
	

# # merge

# #############################################
# ########## 2. Get results
# #############################################

# @transform(runCPAT,
# 		   suffix('.tsv'),
# 		   '_results.tsv')

# def formatCPAT(infile, outfile):

# 	# Rename
# 	rename_dict = {'seq_ID': 'Sequence Name', 'mRNA': 'RNA size', 'ORF': 'ORF size', 'Fickett': 'Ficket Score', 'Hexamer': 'Hexamer Score', 'Coding_prob': 'Coding Probability'}

# 	# Read
# 	cpat_dataframe = pd.read_table(infile).rename(columns=rename_dict)[rename_dict.values()]
# 	cpat_dataframe['Coding Label'] = ['yes' if x >= 0.364 else 'no' for x in cpat_dataframe['Coding Probability']]
# 	cpat_dataframe.index.name = 'Data ID'

# 	# Write
# 	cpat_dataframe.to_csv(outfile, sep='\t', index=True)

# #######################################################
# #######################################################
# ########## S9. Pfam
# #######################################################
# #######################################################

# #############################################
# ########## 1. Translate
# #############################################

# @transform(runCPAT,
# 		   regex(r'(.*)/(s08-cpat.dir)/(.*).ORF_prob.best.tsv'),
# 		   add_inputs(r'\1/\2/\3.ORF_seqs.fa'),
# 		   r'\1/s09-pfam.dir/fasta/Exosc3-translated.fasta')

# def translateORFs(infiles, outfile):

# 	# Run
# 	run_r_job('translate_orfs', infiles, outfile, conda_env='env', W="01:00", GB=10, n=1, run_locally=False, print_outfile=False)

# #############################################
# ########## 2. Split
# #############################################

# @subdivide(translateORFs,
# 		   regex(r'(.*).fasta'),
# 		   r'\1.fasta.*',
# 		   r'\1.fasta.50')

# def splitORFs(infile, outfiles, outfileRoot):

# 	# Get number
# 	N = outfileRoot.split('.')[-1]

# 	# Command
# 	cmd_str = ''' gt splitfasta -numfiles {N} {infile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitORF')

# #############################################
# ########## 3. Run
# #############################################

# @transform(splitORFs,
# 		   regex(r'(.*)/fasta/(.*).fasta\.(.*)'),
# 		   add_inputs(pfam_dir),
# 		   r'\1/split/\2_\3_pfam.txt')

# def runPfamScan(infiles, outfile):

# 	# Data directory
# 	input_fasta, pfam_dir = infiles

# 	# Command
# 	cmd_str = ''' {pfam_dir}/pfam_scan.pl \
# 		-dir {pfam_dir}/data \
# 		-cpu 50 \
# 		-fasta {input_fasta} > {outfile}
# 	'''.format(**locals())
	
# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', modules=['hmmer/3.3'], W='02:00', GB=2, n=10, print_cmd=False, run_locally=False)

# #############################################
# ########## 3. Format
# #############################################

# @transform(runPfamScan,
# 		   suffix('.txt'),
# 		   '_results.tsv')

# def formatPfamResults(infile, outfile):

# 	# Read
# 	pfam_dataframe = pd.read_csv(infile, comment='#', delim_whitespace=True, header=None)

# 	# Get column names
# 	colnames = [x.replace(' ', '_').replace('-', '_') for x in pd.read_csv(infile).iloc[26][0][:-1].replace('# ', '').replace('<', '').split('> ')]

# 	# Add column names
# 	pfam_dataframe.columns = colnames

# 	# Fix sequence ID
# 	pfam_dataframe['seq_id'] = [x.split('_ORF')[0] for x in pfam_dataframe['seq_id']]

# 	# Filter
# 	pfam_dataframe = pfam_dataframe[pfam_dataframe['E_value'] < 0.1]

# 	# Write
# 	pfam_dataframe.to_csv(outfile, index=False, sep='\t')
	
#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

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