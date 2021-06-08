#################################################################
#################################################################
############### Embryo ATAC-Seq - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tibble))
suppressPackageStartupMessages(require(glue))
suppressPackageStartupMessages(require(tidyr))

##### 2. Other libraries #####

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

split_gtf <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read abundance
    abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(transcript_novelty_v2=ifelse(gene_novelty=='Known', transcript_novelty, gene_novelty)) %>% select(transcript_id, transcript_novelty_v2)

    # Merge
    merged_gtf <- gtf %>% left_join(abundance_dataframe, by='transcript_id') %>% replace_na(list(transcript_novelty_v2='Known'))

    # Split
    gtf_split <- split(merged_gtf, merged_gtf$transcript_novelty_v2)

    # Get directory
    outdir <- dirname(outfile)

    # Loop
    for (transcript_class in names(gtf_split)) {
        
        # Get outfile
        gtf_outfile <- glue('{outdir}/{transcript_class}.gtf')
        
        # Export
        rtracklayer::export(gtf_split[[transcript_class]], gtf_outfile, format='gtf') # %>% head(5000)
        
    }

}

#######################################################
#######################################################
########## S5. TSS coverage
#######################################################
#######################################################

#############################################
########## 1. Get TSS BED
#############################################

get_tss <- function(infile, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infile)

    # Extract TSSs
    tss_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end))) %>% distinct %>% mutate(start=tss, end=tss+1, score=0) %>% select(chr, start, end, transcript_id, score, strand)

    # Write
    fwrite(tss_dataframe, file=outfile, sep='\t', col.names=FALSE)

}

#############################################
########## 2. Get counts
#############################################

get_tss_coverage <- function(infiles, outfile, bed_file) {

}

#############################################
#############################################
########## 2. TSS heatmap
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

split_tss <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read abundance
    abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(transcript_novelty_v2=ifelse(gene_novelty=='Known', transcript_novelty, gene_novelty)) %>% select(transcript_id, gene_novelty, transcript_novelty_v2)

    # Extract TSSs
    transcript_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end)), tss_coordinates=paste0('chr', chr, ':', tss)) %>% distinct %>% select(transcript_id, tss_coordinates, strand)

    # Merge
    merged_dataframe <- transcript_dataframe %>% left_join(abundance_dataframe, by='transcript_id') %>% replace_na(list(gene_novelty='Known', transcript_novelty='Known')) %>% replace_na(list(transcript_novelty_v2='Known'))
    head(merged_dataframe)

    # Pivot
    tss_dataframe <- merged_dataframe %>% group_by(tss_coordinates, strand) %>% summarize(transcript_types=paste(unique(transcript_novelty_v2), collapse=',')) %>% 
        mutate(tss_category=ifelse(grepl('Known', transcript_types), 
                                'Known_TSS',
                                ifelse(grepl('Intergenic', transcript_types), 
                                        'Intergenic_TSS', 
                                        ifelse(grepl('Antisense', transcript_types), 
                                                'Antisense_TSS',
                                                'Novel_TSS')))) %>% group_by(tss_category)

    # Split
    dataframes <- setNames(tss_dataframe %>% group_split, tss_dataframe %>% group_keys %>% pull(tss_category))

    # Get directory
    outdir <- dirname(outfile)

    # Loop
    for (tss_class in names(dataframes)) {
        
        # Get bed
        bed_dataframe <- dataframes[[tss_class]] %>% mutate(
            chr=gsub('chr(.*):.*', '\\1', tss_coordinates),
            start=as.numeric(gsub('.*:(.*)', '\\1', tss_coordinates))-1, #-1 is actually unnecessary, because BED is 0-based
            end=start+2, score=0) %>% select(chr, start, end, score, transcript_types, strand)

        # Subset
        # nr_rows <- min(nrow(bed_dataframe), 500)
        # row_idx <- sample(1:nrow(bed_dataframe), nr_rows)
        # bed_dataframe <- as.data.frame(bed_dataframe)[row_idx,]
        
        # Get outfile
        bed_outfile <- glue('{outdir}/{tss_class}.bed')

        # Export
        fwrite(bed_dataframe, bed_outfile, sep='\t', col.names = FALSE)

    }

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################