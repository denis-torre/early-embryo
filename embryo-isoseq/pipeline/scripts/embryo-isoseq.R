#################################################################
#################################################################
############### Embryo IsoSeq - R Support #################
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
########## S2. Align
#######################################################
#######################################################

#############################################
########## 1. Copy FASTQ
#############################################

copy_fastq <- function(infile, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(ShortRead))
    
    # Read
    message(glue('Reading {infile}...'))
    fastq <- readFastq(infile)

    # Find duplicate names
    duplicate_bool <- duplicated(as.character(id(fastq)))

    # Check duplicate names
    if (any(duplicate_bool)) {
        
        # Get indices
        duplicate_idx <- which(duplicate_bool)
        message('Following duplicate sequence names found:')
        message(paste(unique((as.character(id(fastq)[duplicate_bool]))), collapse=', '))
        
        # Get name
        name_dataframe <- id(fastq)[duplicate_idx] %>% as.data.frame %>% group_by(x) %>% mutate(count=row_number()) %>% ungroup %>% mutate(index=duplicate_idx) %>% mutate(new_name=paste0(x, '_', count))    
        
        # Loop
        for (i in 1:nrow(name_dataframe)) {
        
            # Get sequence index
            sequence_index <- name_dataframe[i,'index'][[1]]
            
            # Replace sequence with new name
            fastq[sequence_index] <- ShortReadQ(sread=sread(fastq[sequence_index]), quality=quality(fastq[sequence_index]), id=BStringSet(name_dataframe[i,'new_name'][[1]]))
            
        }

    } else {
        message('No duplicate sequence names found.')
    }

    # Write
    writeFastq(fastq, outfile, compress=FALSE)
}

#######################################################
#######################################################
########## S3. Illumina alignment
#######################################################
#######################################################

#############################################
########## 4. Merge SJ files
#############################################

merge_sjs <- function(infiles, outfile) {
    
    # Split infiles
    infiles_split <- sapply(c('pass1', 'pass2'), function(x) infiles[grepl(x, infiles)], simplify=FALSE)

    # Get colnames
    sj_colnames <- c('chr', 'intron_start', 'intron_end', 'strand', 'intron_motif', 'annotated', 'unique_junction_reads', 'multimapping_junction_reads', 'max_spliced_alignment_overhang')

    # Get annotations
    annotation_dataframe <- lapply(infiles_split[['pass1']], function(x) {
        df <- fread(x);
        colnames(df) <- sj_colnames;
        df <- df %>% select(chr, intron_start, intron_end, strand, intron_motif, annotated)
        return(df)
    }) %>% bind_rows %>% distinct %>% rename('annotated_ensembl'='annotated')

    # Read
    sj_dataframe <- lapply(infiles_split[['pass2']], function(x) {
        df <- fread(x);
        colnames(df) <- sj_colnames;
        return(df)
    }) %>% bind_rows %>% filter(strand != 0) %>% group_by(chr, intron_start, intron_end, strand, intron_motif, annotated) %>% summarize(tot_unique_junction_reads=sum(unique_junction_reads), tot_multimapping_junction_reads=sum(multimapping_junction_reads), max_spliced_alignment_overhang=max(max_spliced_alignment_overhang))

    # Merge and filter
    min_reads <- 3
    merged_dataframe <- sj_dataframe %>% left_join(annotation_dataframe, by=c('chr', 'intron_start', 'intron_end', 'strand', 'intron_motif')) %>% filter((tot_unique_junction_reads >= min_reads | tot_multimapping_junction_reads >= min_reads) | annotated_ensembl == 1) %>% select(-annotated_ensembl)

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t', col.names=FALSE)

}

#######################################################
#######################################################
########## S6. Illumina alignment
#######################################################
#######################################################

#############################################
########## 4. Aggregate counts
#############################################

# aggregate_salmon_counts <- function(infiles, outfile) {

#     # Load
#     require(tximport)

#     # Fix names
#     names(infiles) <- gsub('.*/(.*)/quant.sf', '\\1', infiles)

#     # Import
#     txi <- tximport(infiles, type = "salmon", txOut = TRUE, countsFromAbundance='scaledTPM')

#     # Get dataframe
#     tpm_dataframe <- txi$abundance %>% as.data.frame %>% rownames_to_column('ID')

#     # Write
#     fwrite(tpm_dataframe, file=outfile, sep='\t')

# }

#######################################################
#######################################################
########## S3. Collapse
#######################################################
#######################################################

#############################################
########## 3. Statistics
#############################################

# get_collapse_statistics <- function(infiles, outfile) {
    
#     # Get parameters
#     sample <- gsub('.*/(.*)-5p.*_corrected.gtf', '\\1', infiles[1])
#     max_5_diff <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\1', infiles[1])
#     max_fuzzy_junction <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\2', infiles[1])
#     max_3_diff <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\3', infiles[1])

#     # Read GTF
#     gtf <- rtracklayer::readGFF(infiles[1])

#     # Read SQANTI
#     sqanti_dataframe <- fread(infiles[2]) %>% select(isoform, structural_category, associated_gene) %>% rename('transcript_id'='isoform') %>% mutate(gene_id=gsub('(.*)\\..*', '\\1', transcript_id))

#     # Get genes
#     pb_gene_dataframe <- sqanti_dataframe %>% group_by(gene_id) %>% summarize(ensembl_gene_id=paste(unique(associated_gene), collapse=','))

#     ### Statistics
#     # Number of transcripts per gene
#     transcript_dataframe <- gtf %>% filter(type=='transcript') %>% group_by(gene_id) %>% tally(name = 'nr_transcripts') %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff) %>% merge(pb_gene_dataframe, by='gene_id')

#     # Number of exons per transcript
#     exon_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id) %>% tally(name = 'nr_exons') %>% merge(sqanti_dataframe, by='transcript_id') %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff)

#     # Number of junction-matched transcripts
#     matching_dataframe <- gtf %>% filter(type=='exon') %>% mutate(coordinates=paste0(start, '-', end)) %>% group_by(seqid, gene_id, transcript_id) %>% 
#         summarize(nr_exons=length(coordinates), coordinates_merged=paste0(sort(coordinates), collapse='-')) %>% filter(nr_exons > 1) %>%
#         mutate(junctions=paste0(glue('chr{seqid}:'), '5p-', gsub('.*?-(.*)-.*', '\\1', coordinates_merged), '-3p')) %>% group_by(gene_id, junctions) %>% tally(name='nr_matched_transcripts') %>%
#         merge(pb_gene_dataframe, by='gene_id') %>% arrange(-nr_matched_transcripts) %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff)

#     ### Write
#     fwrite(transcript_dataframe, file=outfile, sep='\t', row.names=FALSE)
#     fwrite(exon_dataframe, file=gsub('transcript', 'exon', outfile), sep='\t', row.names=FALSE)
#     fwrite(matching_dataframe, file=gsub('transcript', 'matching_transcript', outfile), sep='\t', row.names=FALSE)

# }

#############################################
########## 4. Make report
#############################################

# make_collapse_report <- function(infiles, outfile) {

#     # Library
#     suppressPackageStartupMessages(require(ggplot2))
#     suppressPackageStartupMessages(require(ggrepel))
    
#     # Transcripts
#     transcript_dataframe <- lapply(infiles[grepl('\\.transcript_stats.tsv', infiles)], fread) %>% bind_rows %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_transcripts=mean(nr_transcripts))

#     # Exons
#     exon_dataframe <- lapply(infiles[grepl('exon_stats.tsv', infiles)], fread) %>% bind_rows %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_exons=mean(nr_exons), nr_monoexons=sum(nr_exons==1))

#     # Matching transcripts
#     matching_dataframe <- lapply(infiles[grepl('matching_transcript_stats.tsv', infiles)], fread) %>% bind_rows %>% filter(nr_matched_transcripts > 1) %>% mutate(nr_matched_transcripts=nr_matched_transcripts-1) %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_matching_transcripts=mean(nr_matched_transcripts), nr_matching_transcripts=sum(nr_matched_transcripts))

#     # Merge
#     merged_dataframe <- transcript_dataframe %>% merge(exon_dataframe, by=c('max_5_diff', 'max_3_diff')) %>% merge(matching_dataframe, by=c('max_5_diff', 'max_3_diff'))

#     # Plot
#     gp <- ggplot(merged_dataframe, aes(x=nr_monoexons, y=nr_matching_transcripts, color=as.factor(max_5_diff), group=max_3_diff)) +
#         geom_point() +
#         geom_line() +
#         geom_label_repel(data=merged_dataframe %>% filter(max_5_diff==max(merged_dataframe$max_5_diff)) %>% mutate(label=glue('{formatC(max_3_diff, big.mark=",")}bp 3\'')), aes(label=label), nudge_y=10000, size=3, color='black') +
#         scale_color_brewer(type='div', palette=7, direction=-1) +
#         scale_x_continuous(labels = scales::comma) +
#         scale_y_continuous(labels = scales::comma) +
#         labs(x='Mono-exonic transcripts', y='Transcripts with alternate start/end', color='5\' collapse\nthreshold (bp)', title=glue('Transcript collapsing parameters for {strsplit(basename(outfile), "-")[[1]][1]}')) +
#         theme_classic() + theme(plot.title = element_text(hjust = 0.5))

#     # Save
#     ggsave(outfile, plot=gp, height=5, width=9)

# }

#######################################################
#######################################################
########## S4. Merged
#######################################################
#######################################################

#############################################
########## 4. Fix gene IDs
#############################################

# fix_chained_gff <- function(infile, outfile) {

#     # Library
#     require(rtracklayer)

#     # Read
#     gtf <- import(infile)

#     # Get transcript coordinates
#     coordinate_dataframe <- gtf %>% as.data.frame %>% filter(type=='exon') %>% mutate(boundaries=paste0(start, '-', end)) %>% group_by(gene_id, transcript_id) %>% summarize(coordinates=paste0(boundaries, collapse='-'))

#     # Duplicated transcripts
#     duplicate_transcripts <- coordinate_dataframe$transcript_id[which(duplicated(coordinate_dataframe$coordinates))]

#     # Filter
#     gtf_filtered <- gtf[!gtf$transcript_id %in% duplicate_transcripts,]

#     # Write
#     export(gtf_filtered, outfile, format='gtf')
# }

#######################################################
#######################################################
########## S8. SQANTI Second Pass
#######################################################
#######################################################

#############################################
########## 3. ISM filter
#############################################

# filter_isms <- function(infiles, outfile) {

#     # Library
#     require(rtracklayer)

#     # Read GTF
#     gtf <- readGFF(infiles[1])

#     # Get ISMs
#     ism_ids <- fread(infiles[2]) %>% filter(structural_category=='incomplete-splice_match') %>% pull(isoform)

#     # Remove
#     gtf_filtered <- gtf %>% filter(!transcript_id %in% ism_ids)

#     # Write
#     export(gtf_filtered, outfile)
# }

#######################################################
#######################################################
########## S5. TALON
#######################################################
#######################################################

#############################################
########## 8. Get SJs
#############################################

get_junctions <- function(infile, outfile) {
    
    # Library
    require(parallel)

    # Read GTF
    gtf <- rtracklayer::readGFF(infile) %>% select(-source) %>% filter(type=='exon')# %>% head(500)# %>% select(gene_id, transcript_id, seqid, start, end)

    # Get multiexonic transcripts
    exon_counts <- table(gtf$transcript_id)
    multiexon_transcripts <- names(exon_counts)[exon_counts > 1]

    # Filter
    gtf_filtered <- gtf %>% filter(transcript_id %in% multiexon_transcripts)

    # Split
    gtf_split <- split(gtf_filtered, gtf_filtered$transcript_id)

    # Get cores
    cores <- 30

    # Make cluster
    cluster <- makeCluster(cores)

    # Load libraries
    clusterEvalQ(cluster, library("dplyr"));

    # Loop
    junction_dataframe <- parSapply(cluster, gtf_split, function(x) {
        
        # Get GTF
        transcript_gtf <- x %>% arrange(start) %>% select(gene_id, transcript_id, exon_id, seqid, start, end, strand)
        
        # Get strand
        strand <- unique(transcript_gtf$strand)
        
        # Get junctions
        junctions <- transcript_gtf %>% group_by(gene_id, transcript_id, seqid, strand) %>% mutate(junction_start=end, junction_end=lead(start), exon_ids=paste0(exon_id, '-', lead(exon_id))) %>% ungroup %>% head(-1) %>% select(gene_id, transcript_id, exon_ids, seqid, junction_start, junction_end, strand)
        
    }, simplify=FALSE) %>% bind_rows

    # Write
    fwrite(junction_dataframe, outfile, sep='\t')#, compress='gzip')

}

#######################################################
#######################################################
########## S6. CPAT
#######################################################
#######################################################

#############################################
########## 3. Split GTF
#############################################

split_gtf <- function(infile, outfileRoot) {
    
    # Library
    suppressPackageStartupMessages(require(GenomicFeatures))
    suppressPackageStartupMessages(require(rtracklayer))

    # Read GTF
    gtf <- import(infile)

    # Get gene chunks
    gene_ids <- unique(gtf$gene_id)
    gene_chunks <- split(gene_ids, cut(seq_along(1:length(gene_ids)), 50, labels=FALSE))

    # Loop
    for (i in names(gene_chunks)) {

        # Get subset
        message(glue('Doing chunk {i}...'))
        gtf_subset <- gtf[gtf$gene_id %in% gene_chunks[[i]]]
        
        # Get outfile
        chunk_nr <- stringr::str_pad(i, width=2, pad='0')
        outfile <- glue(outfileRoot)
        
        # Write
        export(gtf_subset, outfile, format='gtf')
        
    }
}

#############################################
########## 4. Add CDS
#############################################

add_cds <- function(infiles, outfile, coding_cutoff) {
    
    # Library
    suppressPackageStartupMessages(require(ensembldb))
    suppressPackageStartupMessages(require(rtracklayer))

    ### 1. Create EnsemblDB
    # Read GTF
    gtf <- import(infiles[1])

    # Remove gene boundaries (some novel transcripts are outside of Ensembl-defined ones)
    gtf <- gtf[!gtf$type == 'gene',]

    # Add fake CDS from first exon of first gene
    exon <- gtf[gtf$type=='exon'][1,]
    gtf_cds <- c(gtf, GRanges(seqnames(exon), ranges(exon), strand=strand(exon), type='CDS', source='PacBio', gene_id=exon$gene_id, transcript_id=exon$transcript_id))

    # Get genome info
    organism <- ifelse(grepl('mouse', outfile), 'Mus_musculus', 'Homo_sapiens')
    genomeVersion <- ifelse(grepl('mouse', outfile), 'GRCm38', 'GRCh38')

    # Create database
    db_outfile <- gsub('.gtf', '.sqlite', outfile, fixed=TRUE)
    db_path <- ensDbFromGRanges(gtf_cds, organism = organism, genomeVersion = genomeVersion, version = 102, outfile=db_outfile)
    ensdb <- EnsDb(db_path)

    ### 2. Read CPAT results
    # Read CPAT
    cpat_dataframe <- fread(infiles[2]) %>% dplyr::filter(seq_ID %in% unique(gtf$transcript_id) & Coding_prob >= as.numeric(coding_cutoff))

    # ORF ranges
    orf_ranges <- IRanges(start=cpat_dataframe$ORF_start, width=cpat_dataframe$ORF, names=cpat_dataframe$seq_ID)

    ### 3. Add CDS coordinates
    # Get CDS genomic coordinates
    message('Mapping ORFS to genomic coordinates...')
    cds_ranges <- transcriptToGenome(orf_ranges, ensdb)
    cds_ranges

    # Add phase
    message('Adding phase info...')
    cds_ranges_phase <- lapply(cds_ranges, function(x) {

        # Get CDS mod
        cds_phase <- x %>% as.data.frame %>% arrange(ifelse(strand=='+', start, -start)) %>% mutate(
            cds_mod = (end-start+1)%%3,
            phase = 0,
        )

        # Add phase
        if (nrow(cds_phase) > 1) {
            for (i in 2:nrow(cds_phase)) {
                cds_phase[i,'phase'] <- (3-cds_phase[i-1, 'cds_mod']+cds_phase[i-1, 'phase']) %% 3
            }
        }
        
        # Convert
        result <- makeGRangesFromDataFrame(cds_phase, keep.extra.columns=TRUE);

        # Return
        return(result)

    })

    # Split GTF by transcript to add CDS
    gtf_split <- split(gtf, gtf$transcript_id)

    # Add CDS to GTF if available, merging CDS coordinates by exon ID to preserve exon metadata
    message('Merging...')
    result_gtf <- do.call('c', lapply(names(gtf_split), function(x) {
        result <- gtf_split[[x]];
        if (x %in% names(cds_ranges_phase)) {
            cds_ranges_annotated <- as.data.frame(cds_ranges_phase[[x]]) %>% dplyr::select(start, end, width, phase, exon_id) %>% left_join(as.data.frame(result) %>% dplyr::select(-start, -end, -width, -phase), by='exon_id') %>% mutate(type='CDS') %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE);
            result <- c(result, cds_ranges_annotated)
        }
        return(result)
    }))

    # Subset columns
    result_gtf_subset <- result_gtf %>% as.data.frame %>% dplyr::select(seqnames, start, end, width, strand, source, type, score, phase, gene_id, transcript_id, exon_id, gene_name, transcript_name) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    # Export
    export(result_gtf_subset, outfile, format='gtf')
    
}

#############################################
########## 5. Merge
#############################################

merge_gtf <- function(infiles, outfile) {
    
    # Load
    suppressPackageStartupMessages(require(rtracklayer))

    # Read, concatenate and sort
    merged_gtf <- do.call('c', lapply(infiles, import)) %>% sortSeqlevels %>% sort

    # Write
    export(merged_gtf, outfile, format='gtf')

}

#######################################################
#######################################################
########## S7. Pfam
#######################################################
#######################################################

#############################################
########## 1. Translate
#############################################

translate_orfs <- function(infiles, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(Biostrings))

    # Read ORFs
    orf_dataframe <- fread(infiles[1])

    # Read FASTA
    nucleotide_fasta <- readDNAStringSet(infiles[2], format="fasta")

    # Fix names
    names(nucleotide_fasta) <- gsub('(.*ORF_.).*', '\\1', names(nucleotide_fasta))

    # Subset
    nucleotide_fasta <- nucleotide_fasta[orf_dataframe$ID]

    # Translate
    aa_fasta <- translate(nucleotide_fasta, if.fuzzy.codon='solve')

    # Write
    writeXStringSet(aa_fasta, file=outfile)

}

#######################################################
#######################################################
########## S8. RepeatMasker
#######################################################
#######################################################

#############################################
########## 3. Merge
#############################################

merge_repeatmasker <- function(infiles, outfile) {

    # Read
    repeatmasker_dataframe <- lapply(infiles, function(x) {

        # Read
        df_string <- readChar(x, file.info(x)$size)

        # Replace asterisk
        df_string <- gsub(' *', '*', df_string, fixed=TRUE)

        # Read
        con <- textConnection(df_string)
        df <- read.table(con, fill=TRUE, skip=2)
        close(con)

        # Add column names
        colnames(df) <- c('sw_score', 'pct_div', 'pct_del', 'pct_ins', 'query_sequence', 'query_begin', 'query_end', 'query_left', 'match', 'matching_repeat', 'repeat_class', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID')

        # Return
        return(df)

    }) %>% bind_rows

    # Write
    fwrite(repeatmasker_dataframe, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)

}

#######################################################
#######################################################
########## Summary
#######################################################
#######################################################

#############################################
########## 1. Create
#############################################

get_transcript_summary <- function(infiles, outfile) {

    # Read TALON summary
    talon_dataframe <- fread(infiles[1])

    # Read CPAT
    cpat_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='Sequence Name', 'ORF_length'='ORF size', 'coding_probability'='Coding Probability', 'coding'='Coding Label') %>% select(transcript_id, ORF_length, coding_probability, coding)

    # Read Pfam
    pfam_dataframe <- fread(infiles[3]) %>% rename('transcript_id'='seq_id') %>% group_by(transcript_id) %>% summarize(nr_domains=length(hmm_name), domains=paste0(unique(hmm_name), collapse=','))

    # Read RepeatMasker
    repeat_dataframe <- fread(infiles[4]) %>% rename('transcript_id'='query_sequence') %>% group_by(transcript_id) %>% summarize(nr_repeats=length(matching_repeat), repeats=paste0(unique(matching_repeat), collapse=','))

    # Merge
    merged_dataframe <- list(talon_dataframe, cpat_dataframe, pfam_dataframe, repeat_dataframe) %>% purrr::reduce(left_join, by = "transcript_id") %>% replace_na(list(nr_repeats=0, nr_domains=0, coding='no'))

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t')
    
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################