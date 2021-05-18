#################################################################
#################################################################
############### Embryo Illumina - R Support #################
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
########## S3. Splice junctions
#######################################################
#######################################################

#############################################
########## 3. Junction counts
#############################################

get_junction_counts <- function(infiles, outfile) {

    # Read junctions
    junction_dataframe <- fread(infiles[1]) %>% mutate(junction_start=junction_start+1, junction_end=junction_end-1) %>% rename('chrom'='seqid') %>% select(-exon_ids)

    # Read SJ counts
    count_dataframe <- lapply(infiles[2:length(infiles)], function(x) {
        
        # Read STAR junctions
        sj_dataframe <- fread(x, col.names = c('chrom', 'junction_start', 'junction_end', 'strand', 'intron_motif', 'annotated', 'uniquely_mapped_reads', 'multimapping_reads', 'max_spliced_alignment_overhang')) %>% 
            filter(strand !=0 ) %>% mutate(strand=recode(strand, '1'='+', '2'='-')) %>% select(-intron_motif, -annotated, -max_spliced_alignment_overhang)
        
        # Merge to introns
        merged_dataframe <- junction_dataframe %>% left_join(sj_dataframe, by=c('chrom', 'junction_start', 'junction_end', 'strand')) %>% replace_na(list(uniquely_mapped_reads=0, multimapping_reads=0)) %>% mutate(sample=gsub('.*/(.*)/.*', '\\1', x))

    }) %>% bind_rows

    # Count
    summary_dataframe <- count_dataframe %>% group_by(gene_id, transcript_id, sample) %>% summarize(
        min_unique_sj=min(uniquely_mapped_reads),
        mean_unique_sj=mean(uniquely_mapped_reads),
        max_unique_sj=max(uniquely_mapped_reads),
        min_multi_sj=min(multimapping_reads),
        mean_multi_sj=mean(multimapping_reads),
        max_multi_sj=max(multimapping_reads)
    )# %>% drop_na

    # Write
    fwrite(summary_dataframe, file=outfile, sep='\t')

}

#######################################################
#######################################################
########## S4. Alignment
#######################################################
#######################################################

#############################################
########## 1. Filter GTF
#############################################

filter_gtf <- function(infiles, outfile, comparison) {
    
    # Read GTF
    gtf <- rtracklayer::import(infiles[1])

    # Get outlier samples
    organism <- gsub('.*.dir/(.*?)/.*', '\\1', outfile)
    outlier_samples <- rjson::fromJSON(file=infiles[4])[[organism]]

    # Read and filter abundance
    abundance_dataframe <- fread(infiles[2]) %>% as.data.frame
    abundance_dataframe <- abundance_dataframe[,!colnames(abundance_dataframe) %in% outlier_samples]
    abundance_dataframe$fl_counts <- apply(abundance_dataframe[,12:ncol(abundance_dataframe)], 1, sum)

    # Read and filter junctions
    jc_dataframe <- fread(infiles[3]) %>% mutate(cell_type=gsub('.*?_(.*?)_.*', '\\1', sample))
    if (length(outlier_samples) > 0) {
        jc_dataframe <- jc_dataframe %>% filter(!sample %in% outlier_samples)
    }
    print('Samples used:')
    print(unique(jc_dataframe$sample))
    print('Outlier samples:')
    print(outlier_samples)

    # Filter samples within each comparison
    if (comparison != 'all') {
        jc_dataframe <- jc_dataframe %>% filter(cell_type %in% comparison)
    }

    # Collapse SJ counts
    count_dataframe <- jc_dataframe %>% group_by(transcript_id) %>% summarize(max_min_unique_sj=max(min_unique_sj), max_min_multi_sj=max(min_multi_sj))

    # Get SJ transcripts
    junction_transcripts <- count_dataframe %>% filter(max_min_unique_sj >= 3) %>% pull(transcript_id)

    # Get low FL count isoforms (min 3 for human, min 2 for mouse)
    min_fl_threshold <- ifelse(grepl('human', outfile), 3, 2)
    transcripts_to_remove <- abundance_dataframe %>% filter(fl_counts < min_fl_threshold & transcript_novelty != 'Known') %>% pull(annot_transcript_id)

    # Get SJ-supported, FL-filtered transcripts
    filtered_transcripts <- setdiff(junction_transcripts, transcripts_to_remove)

    # Filter GTF
    gtf_filtered <- gtf[gtf$transcript_id %in% filtered_transcripts,]

    # Export
    rtracklayer::export(gtf_filtered, outfile, format='gtf')
    
}

#######################################################
#######################################################
########## S5. Expression
#######################################################
#######################################################

#############################################
########## 2. Aggregate
#############################################

aggregate_counts <- function(infiles, outfile) {
    
    # Load
    suppressPackageStartupMessages(require(tximeta)) #tximeta_1.8.2
    suppressPackageStartupMessages(require(SummarizedExperiment))
    suppressPackageStartupMessages(require(DESeq2))

    # Get organism
    organism <- gsub('.*/(.*?)_.*.rda', '\\1', outfile)

    # Get sample dataframe
    sample_dataframe <- fread(infiles[1])

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[2]) %>% select(-source) %>% filter(type=='transcript')

    # Get transcript information
    transcript_dataframe <- gtf %>% select(gene_id, gene_name, gene_status, transcript_id, transcript_name, transcript_status, transcript_biotype, NNC_transcript, NIC_transcript, intergenic_transcript, antisense_transcript) %>% distinct
    rownames(transcript_dataframe) <- transcript_dataframe$transcript_id

    # Get gene information
    gene_dataframe <- transcript_dataframe %>% group_by(gene_id, gene_name, gene_status) %>% summarize(
        nr_transcripts=length(unique(transcript_id)),
        novel_transcripts=sum(transcript_status=='NOVEL'),
        NNC_transcripts=sum(NNC_transcript=='TRUE', na.rm=TRUE),
        NIC_transcripts=sum(NIC_transcript=='TRUE', na.rm=TRUE),
        intergenic_transcripts=sum(intergenic_transcript=='TRUE', na.rm=TRUE),
        antisense_transcripts=sum(antisense_transcript=='TRUE', na.rm=TRUE)
    ) %>% as.data.frame
    rownames(gene_dataframe) <- gene_dataframe$gene_id

    # Read isoform counts
    se <- tximeta(sample_dataframe, type='rsem', txIn=TRUE, txOut=TRUE)
    rowData(se) <- transcript_dataframe[rownames(rowData(se)),]

    # Read gene counts
    gse <- tximeta(sample_dataframe, type='rsem', txIn=TRUE, txOut=FALSE, tx2gene=transcript_dataframe %>% select(transcript_id, gene_id))
    rowData(gse) <- gene_dataframe[rownames(rowData(gse)),]

    # Create DESeq dataset
    if (organism == 'human') {
        dds_list <- list(
            'transcript' = DESeqDataSet(se, design=~cell_type+batch), #+quality
            'gene' = DESeqDataSet(gse, design=~cell_type+batch) #+quality
        )
    } else if (organism == 'mouse') {
        dds_list <- list(
            'transcript' =  DESeqDataSet(se, design=~cell_type),
            'gene' = DESeqDataSet(gse, design=~cell_type)
        )
    }

    # Save
    save(se, gse, dds_list, file=outfile)

}

#############################################
########## 3. TPM
#############################################

get_transcript_tpm <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(library(DESeq2))

    # Load
    load(infile)
    
    # Get tpm
    tpm_dataframe <- assay(se, 'abundance') %>% as.data.frame

    # Write
    write.table(tpm_dataframe, file=outfile, quote=FALSE, sep='\t')

}

#######################################################
#######################################################
########## S5. Differential expression
#######################################################
#######################################################

#############################################
########## 1. DESeq2
#############################################

run_deseq2 <- function(infiles, outfileRoot, feature) {
    
    # Library
    suppressPackageStartupMessages(require(rjson))
    suppressPackageStartupMessages(require(DESeq2))

    # Load
    load(infiles[1])

    # Read comparisons
    organism <- gsub('.*.dir/(.*?)/(.*?)/.*', '\\1', outfileRoot)
    comparisons <- fromJSON(file = infiles[2])[[organism]]

    # Get DDS
    dds <- dds_list[[feature]]

    # Filter lowly expressed transcripts or genes
    # keep_sample <- colData(dds)$cell_type %in% comparison
    # keep_gene <- rowSums(counts(dds[,keep_sample])) >= 5
    # dds <- dds[keep_gene,]

    # Run DESeq
    dds <- DESeq(dds)

    # Loop through comparisons
    for (comparison in comparisons) {

        # Get results
        deseq_dataframe <- results(dds, contrast = c('cell_type', comparison[2], comparison[1]), alpha=0.05) %>% as.data.frame %>% rownames_to_column(ifelse(feature == 'transcript', 'transcript_id', 'gene_id')) %>% filter(baseMean > 0)
        
        # Add annotation
        if (feature == 'gene') {
            gene_dataframe <- as.data.frame(rowData(dds))[,1:5]
            result_dataframe <- gene_dataframe %>% right_join(deseq_dataframe, by='gene_id')
        } else if (feature == 'transcript') {
            transcript_dataframe <- as.data.frame(rowData(dds)) %>% select(transcript_id, transcript_name, transcript_status, transcript_biotype, gene_id, gene_name, gene_status, NNC_transcript, NIC_transcript, intergenic_transcript, antisense_transcript)
            result_dataframe <- deseq_dataframe %>% left_join(transcript_dataframe, by='transcript_id')
        }

        # Sort
        result_dataframe <- result_dataframe %>% arrange(pvalue)

        # Get outfile
        outfile <- glue(outfileRoot)

        # Write
        fwrite(result_dataframe, file=outfile, sep='\t')

    }

        # Check if SQANTI
        # if ('sqanti_dataframe' %in% ls()) {

        #     # Select columns
        #     sqanti_dataframe <- sqanti_dataframe %>% mutate(gene=gsub('(.*)\\..*', '\\1', isoform)) %>% select(isoform, gene, structural_category, associated_gene, associated_transcript, subcategory)
            
        #     # Feature level
        #     if (feature == 'transcript') {

        #         # Merge
        #         deseq_dataframe <- deseq_dataframe %>% merge(sqanti_dataframe, by='isoform')

        #     } else if (feature == 'gene') {

        #         # Merge
        #         deseq_dataframe <- deseq_dataframe %>% merge(sqanti_dataframe %>% select(gene, associated_gene) %>% distinct, by='gene')

        #     }

        #     # Read genes
        #     gtf <- rtracklayer::import(infiles[2])

        #     # Make dataframe
        #     gene_dataframe <- gtf %>% as.data.frame %>% select(gene_id, gene_name) %>% distinct %>% rename('gene_id'='associated_gene')

        #     # Merge
        #     deseq_dataframe <- gene_dataframe %>% merge(deseq_dataframe, by='associated_gene') %>% arrange(padj)

        # }

        # # Write
        # fwrite(deseq_dataframe, file=outfile, sep='\t')

}

#######################################################
#######################################################
########## S9. Enrichment
#######################################################
#######################################################

#############################################
########## 1. GO
#############################################

run_go_enrichment <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(msigdbr))
    suppressPackageStartupMessages(require(fgsea))

    # Read signature
    deseq_dataframe <- fread(infile)
    gene_signature <- deseq_dataframe %>% filter(!grepl('TALON', gene_name)) %>% group_by(gene_name) %>% slice_min(order_by=padj, n=1) %>% ungroup %>% pull(stat, gene_name)

    
    # Get geneset
    organism <- ifelse(grepl('human', infile), 'Homo sapiens', 'Mus musculus')
    ms <- msigdbr(species = organism, category='C5', subcategory='BP')
    genesets <- ms %>% group_by(gs_name) %>% summarize(genes=list(unique(gene_symbol))) %>% pull(genes, gs_name)

    # Run
    gsea_dataframe <- fgsea(pathways = genesets, stats = gene_signature, minSize  = 15, maxSize  = 5000, eps=0) %>% select(pathway, pval, padj, NES, size) %>% arrange(padj)

    # Write
    fwrite(gsea_dataframe, file=outfile, sep='\t')

}


#############################################
########## 2. Domain
#############################################

run_domain_enrichment <- function(infiles, outfile) {

    # Library
    suppressPackageStartupMessages(require(msigdbr))
    suppressPackageStartupMessages(require(fgsea))

    # Read signature
    deseq_dataframe <- fread(infiles[1])
    gene_signature <- deseq_dataframe %>% pull(stat, transcript_id)
    
    # Get geneset
    domains <- fread(infiles[2]) %>% group_by(hmm_name) %>% summarize(genes=list(unique(seq_id))) %>% pull(genes, hmm_name)

    # Run
    gsea_dataframe <- fgsea(pathways = domains, stats = gene_signature, minSize  = 5, maxSize  = 5000, eps=0) %>% select(pathway, pval, padj, NES, size) %>% arrange(padj)

    # Write
    fwrite(gsea_dataframe, file=outfile, sep='\t')
}

#############################################
########## 3. Repeats
#############################################

run_repeat_enrichment <- function(infiles, outfile) {

    # Library
    suppressPackageStartupMessages(require(msigdbr))
    suppressPackageStartupMessages(require(fgsea))

    # Read signature
    deseq_dataframe <- fread(infiles[1])
    gene_signature <- deseq_dataframe %>% pull(stat, transcript_id)
    
    # Get geneset
    repeats <- fread(infiles[2]) %>% group_by(matching_repeat) %>% summarize(genes=list(unique(query_sequence))) %>% pull(genes, matching_repeat)

    # Run - increase maxSize
    gsea_dataframe <- fgsea(pathways = repeats, stats = gene_signature, minSize  = 5, maxSize  = 10000, eps=0) %>% select(pathway, pval, padj, NES, size) %>% arrange(padj)

    # Write
    fwrite(gsea_dataframe, file=outfile, sep='\t')
}

#######################################################
#######################################################
########## S10. Isoform switching
#######################################################
#######################################################

#############################################
########## 1. Load data
#############################################

load_isoform_data <- function(infiles, outfile) {

    # Library
    suppressPackageStartupMessages(require(IsoformSwitchAnalyzeR))
    suppressPackageStartupMessages(require(rjson))

    # Read metadata
    metadata_dataframe <- fread(infiles[1])

    # Read expression
    salmon_infiles <- metadata_dataframe %>% pull(files, names)
    salmonQuant <- importIsoformExpression(sampleVector = salmon_infiles)

    # Get organism
    organism <- gsub('.*/(.*?)_.*', '\\1', infiles[1])

    # Get metadata
    if (organism == 'mouse') {
        design_dataframe <- metadata_dataframe %>% dplyr::rename('sampleID'='names', 'condition'='cell_type') %>% dplyr::select(sampleID, condition)
    } else if (organism == 'human') {
        design_dataframe <- metadata_dataframe %>% dplyr::rename('sampleID'='names', 'condition'='cell_type') %>% dplyr::select(sampleID, condition, batch)
    }
    design_dataframe <- design_dataframe %>% mutate(condition=paste0('embryo_', condition))

    # Read comparisons
    comparison_dataframe <- fromJSON(file = infiles[6])[[organism]]%>% as.data.frame %>% apply(1, function(x) paste0('embryo_', x)) %>% as.data.frame %>% dplyr::rename('condition_1'='V1', 'condition_2'='V2')
    rownames(comparison_dataframe) <- NULL

    # Create switchAnalyzeRlist
    aSwitchList <- importRdata(
        isoformCountMatrix   = salmonQuant$counts,
        isoformRepExpression = salmonQuant$abundance,
        designMatrix         = design_dataframe,
        comparisonsToMake    = comparison_dataframe,
        isoformExonAnnoation = infiles[2],
        isoformNtFasta       = infiles[3],
        ignoreAfterPeriod    = FALSE
    )

    # Add CPAT
    aSwitchList <- analyzeCPAT(
        switchAnalyzeRlist   = aSwitchList,
        pathToCPATresultFile = infiles[4],
        codingCutoff         = ifelse(organism=='human', 0.364, 0.44),
        removeNoncodinORFs   = FALSE   # because ORF was added from CPAT
    )

    # Add Pfam
    aSwitchList <- analyzePFAM(
        switchAnalyzeRlist   = aSwitchList,
        pathToPFAMresultFile = infiles[5],
        showProgress=FALSE
    )

    # Save
    save(aSwitchList, file=outfile)
}

#############################################
########## 2. Run
#############################################

get_isoform_switching <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(IsoformSwitchAnalyzeR))

    # Load
    load(infile)

    # Filter
    switchListFiltered <- preFilter(
        aSwitchList,
        geneExpressionCutoff = 10, # default 1
        isoformExpressionCutoff = 3, # default 0
        IFcutoff=0.01,
        acceptedGeneBiotype = NULL,
        acceptedIsoformClassCode = NULL,
        removeSingleIsoformGenes = TRUE,
        reduceToSwitchingGenes=FALSE,
        reduceFurtherToGenesWithConsequencePotential = FALSE,
        onlySigIsoforms = FALSE,
        keepIsoformInAllConditions=FALSE,
        alpha=0.05,
        dIFcutoff = 0.1,
        quiet=FALSE
    )

    # Run DEXSeq
    switchListAnalyzed <- isoformSwitchTestDEXSeq(
        switchAnalyzeRlist = switchListFiltered,
        reduceToSwitchingGenes=TRUE
    )

    # Alternative splicing
    switchListAnalyzed <- analyzeAlternativeSplicing(
        switchListAnalyzed,
        onlySwitchingGenes=TRUE,
        alpha=0.05,
        dIFcutoff = 0.1
    )

    # Consequences
    consequencesToAnalyze <- c('intron_retention', 'coding_potential', 'NMD_status', 'domains_identified')
    
    # Switch consequences
    switchListAnalyzed_SwitchConseq <- analyzeSwitchConsequences(
        switchListAnalyzed,
        consequencesToAnalyze=consequencesToAnalyze,
        alpha=0.05,
        dIFcutoff=0.1,
        onlySigIsoforms=FALSE,
        ntCutoff=50,
        ntFracCutoff=NULL,
        ntJCsimCutoff=0.8,
        AaCutoff=10,
        AaFracCutoff=0.5,
        AaJCsimCutoff=0.9,
        removeNonConseqSwitches=TRUE
    )

    # Save
    save(switchListAnalyzed, switchListAnalyzed_SwitchConseq, file=outfile)

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################