#################################################################
#################################################################
############### Embryo ATAC-Seq ################
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
r_source = 'pipeline/scripts/embryo-atacseq.R'
py_source = 'pipeline/scripts/EmbryoAtacseq.py'
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
#import EmbryoAtacseq as P

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
### Samples
# Read samples
liu_dataframe = pd.read_csv('/hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/liu/liu-samples.csv', comment='#')
liu_dataframe = liu_dataframe.rename(columns={x: x.replace(' ', '_') for x in liu_dataframe.columns})
liu_dataframe['fastq'] = ['/hpc/users/torred23/pipelines/projects/early-embryo/embryo-atacseq/arion/datasets/liu/rawdata/{x}.fastq.gz'.format(**locals()) for x in liu_dataframe['Run']]
liu_dataframe['sample_name'] = ['human_'+x.replace('_A', '').replace('Zygote', '1C').replace('Morula', 'morula').replace('-cell', 'C') for x in liu_dataframe['Sample_Name']]

# Samples
sample_dict = {
    'human': liu_dataframe.query('Assay_Type == "ATAC-seq" and BioSampleModel == "Human"')[['Run', 'sample_name', 'fastq']].rename(columns={'Run': 'run'}).to_dict(orient='records')
}

# All FASTQ
all_illumina_fastq = glob.glob('arion/atacseq/s01-fastq.dir/*/*/*/*.f*q.gz')

# Reference
reference_dict = {
	'human': {
		'filtered_gtf': 'arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf',
		'talon_abundance': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_abundance_filtered.tsv'		
	}
}

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def qcJobs():
	for organism, samples in sample_dict.items():
		for sample_info in samples:
			infile = sample_info['fastq']
			outfile = 'arion/atacseq/s01-fastq.dir/{organism}/raw/{sample_name}/{sample_name}.fastq.gz'.format(**locals(), **sample_info)
			yield [infile, outfile]

@files(qcJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/atacseq/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

@follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	if len(infiles) == 1:
		cmd_str = '''trim_galore --nextera --cores 6 --output_dir {outdir} {infiles[0]}'''.format(**locals())
	elif len(infiles) == 2:
		cmd_str = '''trim_galore --nextera --paired --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='02:00', GB=6, n=6, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

# find arion/illumina/s01-fastq.dir/*/trimmed -name "job.log" | jsc

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

@follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/atacseq/s02-fastqc.dir/human/raw', 'arion/atacseq/multiqc/human_fastqc/multiqc_report.html'],
		['arion/atacseq/s02-fastqc.dir/human/trimmed', 'arion/atacseq/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/atacseq/s03-alignment.dir/human/bowtie2', 'arion/atacseq/multiqc/human_alignment/multiqc_report.html'],
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
########## 1. Bowtie index
#############################################

@transform('arion/datasets/reference_genomes/*/*.primary_assembly.fa',
		   regex(r'.*/(.*)/(.*).fa'),
		   r'arion/atacseq/s03-alignment.dir/\1/bowtie2/index/\2.1.bt2')

def buildBowtieIndex(infile, outfile):

	# Command
	outname = outfile.replace('.1.bt2', '')
	cmd_str = '''bowtie2-build --threads 30 {infile} {outname}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bowtie2/2.4.1'], W='02:00', GB=6, n=6, stdout=outfile.replace('.1.bt2', '.log'), stderr=outfile.replace('.1.bt2', '.err'))

#############################################
########## 2. Bowtie
#############################################

@follows(trimIlluminaAdapters, buildBowtieIndex)

@collate('arion/atacseq/s01-fastq.dir/human/trimmed/*/*_trimmed.fq.gz',
 		 regex(r'(.*)/s01-fastq.dir/(.*)/trimmed/(.*)/.*.fq.gz'),
		 add_inputs(r'\1/s03-alignment.dir/\2/bowtie2/index/*primary_assembly.1.bt2'),
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

# find arion/atacseq/s03-alignment.dir -name "*.log" | js

#############################################
########## 3. Filter
#############################################
# Removes duplicates, and any reads not mapping to chromosomes 1-19, X, Y.
# "[XS] == null" - multimappers
# "not unmapped"
# "not duplicate"
# "ref_id <= 21 and ref_id != 19" - chromosomes 1-19,X,Y (no chrM)
# "mapping_quality >= 30"
# fragment lengths is empty for single end data

@transform(runBowtie,
		   suffix('.bam'),
		   '_filtered.bam')

def filterBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "ref_id <= 21 and ref_id != 19 and not unmapped and not duplicate and mapping_quality >= 30" {infile} > {outfile} && \
		samtools index -b {outfile} && samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | awk '$9>0' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > {outname}.fragment_lengths.txt'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'], print_outfile=False)

#############################################
########## 4. Create BigWig
#############################################

@transform((runBowtie, filterBam),
		   suffix('.bam'),
		   add_inputs('/sc/arion/projects/GuccioneLab/genome-indices/hg38/blacklists/hg38-blacklist.v2.bed'),
		   '.bw')

def createBigWig(infiles, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=10 --skipNonCoveredRegions --numberOfProcessors=48 --normalizeUsing RPKM --blackListFileName {infiles[1]} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=8, GB=4, print_outfile=False)

#############################################
########## 5. Sort BAM by QNAME (for Genrich)
#############################################

@transform(filterBam,
		   suffix('.bam'),
		   '_nsorted.bam')

def sortBAM(infile, outfile):

	# Command
	cmd_str = ''' samtools sort -n --threads 30 -o {outfile} {infile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.11'], W='02:00', n=6, GB=6, print_outfile=False)

#######################################################
#######################################################
########## S4. Peaks
#######################################################
#######################################################

#############################################
########## 1. MACS2
#############################################

@collate(filterBam,
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/(.*)_Rep.*/.*_filtered.bam'),
		 r'\1/s04-peaks.dir/\2/macs2/\3/\3_peaks.xls')

def runMacs2(infiles, outfile):

	# Infiles
	infiles_str = ' '.join(infiles)
	basename = outfile.replace('_peaks.xls', '')
	organism = outfile.split('/')[-4].replace('human', 'hs').replace('mouse', 'mm')

	# Command
	cmd_str = ''' macs2 callpeak \
		-t {infiles_str}\
		-f BAM \
		--nomodel \
		--nolambda \
		--extsize 250 \
		-q 0.01 \
		-n {basename} \
		-g {organism}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['macs/2.1.0'], W='06:00', n=1, GB=30, print_cmd=False, stdout=outfile.replace('.xls', '.log'), stderr=outfile.replace('.xls', '.err'))

#############################################
########## 2. gEnrich
#############################################

@collate(sortBAM,
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/(.*)_Rep.*/.*_filtered_nsorted.bam'),
		 r'\1/s04-peaks.dir/\2/genrich/\3/\3-genrich.narrowPeak')

def runGenrich(infiles, outfile):

	# Infiles
	bam_str = ','.join(infiles)
	basename = outfile.replace('.narrowPeak', '')

	# Command
	cmd_str = ''' Genrich -t {bam_str} \
		-o {outfile} \
		-q 0.01 \
		-j -y -v '''.format(**locals())
		# -b {basename}.bed \
		# -R {basename}.pcr_duplicates.txt \
		# -r

	# Run
	run_job(cmd_str, outfile, modules=['genrich/0.6'], W='03:00', n=1, GB=30, print_cmd=False, stdout=outfile.replace('.narrowPeak', '.log'), stderr=outfile.replace('.narrowPeak', '.err'))

#######################################################
#######################################################
########## S5. TSS coverage
#######################################################
#######################################################

#############################################
########## 1. Get TSS BED
#############################################

def transcriptJobs():
	for organism, reference_info in reference_dict.items():
		infile = reference_info['filtered_gtf']
		outfile = 'arion/atacseq/s05-tss_coverage.dir/{organism}/{organism}-tss.bed'.format(**locals())
		yield [infile, outfile]

@files(transcriptJobs)

def getTssRange(infile, outfile):

	# Run
	run_r_job('get_tss', infile, outfile, run_locally=True)#, conda_env='env')#, modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. Get counts
#############################################

# @follows(filterBam, getTssRange)

@collate('arion/atacseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bam',
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/.*/.*.bam'),
		 add_inputs(r'\1/s05-tss_coverage.dir/\2/\2-tss.bed'),
		 r'\1/s05-tss_coverage.dir/\2/\2-tss_coverage.tsv')

def getTssCoverage(infiles, outfile):

	# Split files
	bam_files = [x[0] for x in infiles]
	bed_file = infiles[0][1]

	# Prepare
	bam_str = ' '.join(bam_files)
	bam_str_tsv = '\t'.join([x.split('/')[-2] for x in bam_files])

	# Command
	cmd_str = ''' echo "chr\tstart\tend\ttranscript_id\tscore\tstrand\t{bam_str_tsv}" > {outfile} && bedtools multicov -bams {bam_str} -bed {bed_file} >> {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='01:00', n=1, GB=25)#, print_cmd=False, stdout=outfile.replace('.narrowPeak', '.log'), stderr=outfile.replace('.narrowPeak', '.err'))

	# # Run
	# run_r_job('get_tss_coverage', bam_files, outfile, additional_params=bed_file, run_locally=True)#conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 3. Get peaks
#############################################

# @follows(runGenrich, getTssRange)

@transform('arion/atacseq/s04-peaks.dir/human/genrich/*/*-genrich.narrowPeak',
		   regex(r'(.*)/s04-peaks.dir/(.*)/genrich/(.*)/.*.narrowPeak'),
		   add_inputs(r'\1/s05-tss_coverage.dir/\2/\2-tss.bed'),
		   r'\1/s05-tss_coverage.dir/\2/peaks/\2-\3-tss_peaks_intersect.bed')

def getTssPeaks(infiles, outfile):

	# Command
	cmd_str = ''' bedtools intersect -wa -a {infiles[1]} -b {infiles[0]} > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='00:30', n=1, GB=25, print_cmd=False)#, stdout=outfile.replace('.narrowPeak', '.log'), stderr=outfile.replace('.narrowPeak', '.err'))

#######################################################
#######################################################
########## Plots
#######################################################
#######################################################

#############################################
#############################################
########## 1. Isoform heatmap
#############################################
#############################################

#############################################
########## 1.1 Split GTF by isoform class
#############################################

def splitJobs():
	for organism, reference_info in reference_dict.items():
		infiles = list(reference_info.values())
		outfile = 'arion/atacseq/summary_plots.dir/isoform_heatmaps_5000/{organism}/gtf/Known.gtf'.format(**locals())
		yield [infiles, outfile]

@files(splitJobs)

def splitGTF(infiles, outfile):

	# Run
	run_r_job('split_gtf', infiles, outfile, run_locally=True)#conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 1.2 Matrix
#############################################

# @follows(createBigWig, splitGTF)

@collate('arion/atacseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bw',
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/.*/.*.bw'),
		 add_inputs(r'\1/summary_plots.dir/isoform_heatmaps_5000/\2/gtf/*.gtf'),
		 r'\1/summary_plots.dir/isoform_heatmaps_5000/\2/\2-matrix.gz')

def computeIsoformMatrix(infiles, outfile):

	# Get order
	order_dict = {
		'bigwig': {'1C': 1, '2C': 2, '4C': 3, '8C': 4, 'morula': 5, 'ICM': 6, 'TE': 7},
		'gtf': {'Known': 1, 'NIC': 2, 'NNC': 3, 'Antisense': 4, 'Intergenic': 5}
	}

	# Split
	bigwigs = [x[0] for x in infiles]
	gtfs = infiles[0][1:]

	# Get bigwig order
	bigwig_str = ' '.join(pd.DataFrame({
		'bigwig': bigwigs,
		'order': [order_dict['bigwig'][os.path.basename(x).split('_')[1]] for x in bigwigs]
	}).sort_values('order')['bigwig'])

	# Get GTF order
	gtf_str = ' '.join(pd.DataFrame({
		'gtf': gtfs,
		'order': [order_dict['gtf'][os.path.basename(x).split('.')[0]] for x in gtfs]
	}).sort_values('order')['gtf'])

	# Command
	cmd_str = ''' computeMatrix scale-regions -S {bigwig_str} \
					-R {gtf_str} \
					--beforeRegionStartLength 3000 \
					--regionBodyLength 5000 \
					--afterRegionStartLength 3000 \
					--numberOfProcessors 48 \
					--skipZeros -o {outfile}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=6, GB=4, print_cmd=False, ow=False, wait=True)

	# Write samples
	bigwig_names = [x.split('/')[-2].replace('_', ' ') for x in bigwig_str.split(' ')]
	jsonfile = outfile.replace('.gz', '.json')
	if not os.path.exists(jsonfile):
		with open(jsonfile, 'w') as openfile:
			openfile.write(json.dumps(bigwig_names))

#############################################
########## 1.3 Plot
#############################################

@transform(computeIsoformMatrix,
		   regex(r'(.*).gz'),
		   add_inputs(r'\1.json'),
		   r'\1_heatmap.png')

def plotIsoformHeatmap(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		samples_label = '" "'.join(json.load(openfile))

	# Command
	cmd_str = ''' plotHeatmap -m {infiles[0]} \
					--samplesLabel "{samples_label}" \
					--heatmapWidth 7 \
					-out {outfile}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=True)

#############################################
########## 1.4 Plot
#############################################

# @transform('arion/atacseq/summary_plots.dir/isoform_heatmaps_500/human/human-matrix.gz',
@transform(computeIsoformMatrix,
		   regex(r'(.*).gz'),
		   add_inputs(r'\1.json'),
		   r'\1_profile.png')

def plotIsoformProfile(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		samples_label = '" "'.join(json.load(openfile))

	# Command
	cmd_str = ''' plotProfile -m {infiles[0]} \
		--numPlotsPerRow 4 \
		--plotWidth 10 \
		--colors "#6BAED6" "#78C679" "#EE6A50" "#66C2A4" "#E9967A"\
		--samplesLabel "{samples_label}" \
		-out {outfile}
	'''.format(**locals())
		# --perGroup \

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=True)

#############################################
########## 1.4 Plot
#############################################

@transform(computeIsoformMatrix,
		   regex(r'(.*).gz'),
		   add_inputs(r'\1.json'),
		   r'\1_profile_v2.png')

def plotIsoformProfile2(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		# samples_label = '" "'.join(json.load(openfile))
		samples_label = [' ' for x in json.load(openfile)]

	# Command
	cmd_str = ''' plotProfile -m {infiles[0]} \
		--numPlotsPerRow 3 \
		--perGroup \
		--plotWidth 10 \
		--samplesLabel {samples_label} \
		-out {outfile}
	'''.format(**locals())
		# --colors "#6BAED6" "#78C679" "#EE6A50" "#66C2A4" "#E9967A"\

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=True)

#############################################
#############################################
########## 2. TSS heatmap
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

def tssJobs():
	for organism, reference_info in reference_dict.items():
		infiles = list(reference_info.values())
		outfile = 'arion/atacseq/summary_plots.dir/tss_heatmaps/{organism}/bed_5k/Known_TSS.bed'.format(**locals())
		yield [infiles, outfile]

@files(tssJobs)

def splitTSS(infiles, outfile):

	# Run
	run_r_job('split_tss', infiles, outfile, run_locally=True)#conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 1.2 Matrix
#############################################

# @follows(createBigWig, splitGTF)

@collate('arion/atacseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bw',
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/.*/.*.bw'),
		 add_inputs(r'\1/summary_plots.dir/tss_heatmaps/\2/bed_5k/*.bed'),
		 r'\1/summary_plots.dir/tss_heatmaps/\2/\2-matrix_5k.gz')

def computeTssMatrix(infiles, outfile):

	# Get order
	order_dict = {
		'bigwig': {'1C': 1, '2C': 2, '4C': 3, '8C': 4, 'morula': 5, 'ICM': 6, 'TE': 7},
		'bed': {'Known_TSS': 1, 'Novel_TSS': 2, 'Antisense_TSS': 3, 'Intergenic_TSS': 4}
	}

	# Split
	bigwigs = [x[0] for x in infiles]
	bigwigs = [x for x in bigwigs if 'Rep1' in x]
	beds = infiles[0][1:]

	# Get bigwig order
	bigwig_str = ' '.join(pd.DataFrame({
		'bigwig': bigwigs,
		'order': [order_dict['bigwig'][os.path.basename(x).split('_')[1]] for x in bigwigs]
	}).sort_values('order')['bigwig'])

	# Get GTF order
	bed_str = ' '.join(pd.DataFrame({
		'bed': beds,
		'order': [order_dict['bed'][os.path.basename(x).split('.')[0]] for x in beds]
	}).sort_values('order')['bed'])

	# Command
	cmd_str = ''' computeMatrix reference-point -S {bigwig_str} \
					-R {bed_str} \
					--referencePoint center \
					--beforeRegionStartLength 500 \
					--afterRegionStartLength 500 \
					--numberOfProcessors 48 \
					--skipZeros -o {outfile}
	'''.format(**locals())

	# Run
	# run_job(cmd_str, outfile, conda_env='env', W='10:00', n=6, GB=4, print_cmd=False, ow=False, wait=False)

	# Write samples
	bigwig_names = [x.split('/')[-2].replace('_', ' ') for x in bigwig_str.split(' ')]
	jsonfile = outfile.replace('.gz', '.json')
	if not os.path.exists(jsonfile):
		with open(jsonfile, 'w') as openfile:
			openfile.write(json.dumps(bigwig_names))

#############################################
########## 1.3 Plot
#############################################

@transform(computeTssMatrix,
		   regex(r'(.*)_(.*).gz'),
		   add_inputs(r'\1_\2.json'),
		   r'\1_heatmap_\2.png')

def plotTssHeatmap(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		samples_label = '" "'.join(json.load(openfile))

	# Command
	cmd_str = ''' plotHeatmap -m {infiles[0]} \
					--samplesLabel "{samples_label}" \
					--heatmapWidth 7 \
					--refPointLabel TSS \
					-out {outfile}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=True)

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