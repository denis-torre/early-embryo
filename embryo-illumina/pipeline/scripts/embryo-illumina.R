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
    filtered_multiexonic_transcripts <- setdiff(junction_transcripts, transcripts_to_remove)

    # Get FL-filtered, monoexonic transcripts (won't have SJ support)
    filtered_monoexonic_transcripts <- abundance_dataframe %>% filter(n_exons == 1 & transcript_novelty == 'Known' & fl_counts >= min_fl_threshold) %>% pull(annot_transcript_id)

    # Concatenate
    filtered_transcripts <- c(filtered_multiexonic_transcripts, filtered_monoexonic_transcripts)

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

#############################################
########## 3. Gene counts
#############################################

get_gene_expression <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(library(DESeq2))

    # Load
    load(infile)
    
    # Get expression
    # expression_type <- gsub('.*gene_(.*).tsv', '\\1', outfile)
    # abundance_dataframe <- assay(gse, expression_type) %>% as.data.frame %>% rownames_to_column('gene_id')
    dds <- estimateSizeFactors(dds_list[['gene']])
    expression_dataframe <- counts(dds, normalized=TRUE) %>% as.data.frame %>% rownames_to_column('gene_id')

    # Write
    write.table(expression_dataframe, file=outfile, quote=FALSE, sep='\t', row.names=FALSE)

}

#############################################
########## 5. Get size factors
#############################################

get_size_factors <- function(infile, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(DESeq2))

    # Load
    load(infile)

    # Get counts
    count_matrix <- counts(dds_list[['gene']])

    # Get factors
    normalization_factors <- edgeR::calcNormFactors(object = count_matrix, method = "TMM")

    # Get library sizes
    library_sizes <- colSums(count_matrix)

    # Get dataframe
    normalization_dataframe <- data.frame(normalization_factor=normalization_factors, library_size=library_sizes) %>% rownames_to_column('sample_name') %>% mutate(size_factor=normalization_factor*library_size/1000000, size_factor_reciprocal=1/size_factor)

    # Write
    fwrite(normalization_dataframe, file=outfile, sep='\t')

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
########## S7. Enrichment
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
########## S8. WGCNA
#######################################################
#######################################################

#############################################
########## 1. Pick soft thresholds
#############################################

pick_soft_thresholds <- function(infile, outfile) {

    # Library
    require(DESeq2)
    require(WGCNA)

    # Load
    load(infile)
    
    # Get dds
    dds <- estimateSizeFactors(dds_list[['gene']])

    # Get counts
    count_dataframe <- counts(dds, normalized=FALSE)

    # Get gene counts
    gene_counts <- apply(count_dataframe, 1, sum)
    
    # Filter genes
    good_genes <- names(gene_counts)[gene_counts > 10]

    # Subset
    dds_subset <- dds[good_genes,]
    
    # Get expression
    expression_dataframe <- counts(dds_subset, normalized=TRUE)

    # Filter expression
    log1p_dataframe <- log10(expression_dataframe[good_genes,]+1) %>% t

    # Pick soft threshold
    powers <- 1:20
    sft <- pickSoftThreshold(log1p_dataframe, powerVector=powers, networkType='signed')

    # Save
    save(sft, powers, log1p_dataframe, file=outfile)

}

#############################################
########## 2. Cluster genes
#############################################

cluster_genes <- function(infile, outfile) {

    # Library
    require(WGCNA)

    # Load
    load(infile)

    # Get adjacency
    adjacency <- adjacency(log1p_dataframe, power = 13, type='signed')

    # Turn adjacency into topological overlap
    TOM <- TOMsimilarity(adjacency, TOMType = "signed");
    dissTOM <- 1-TOM # ??? correct, between 1 and 0 even for signed
    rownames(dissTOM) <- colnames(log1p_dataframe)
    colnames(dissTOM) <- colnames(log1p_dataframe)

    # Cluster genes
    geneTree <- hclust(as.dist(dissTOM), method = "average")

    # Save
    save(geneTree, dissTOM, file=outfile)
    save(adjacency, file=gsub('.rda', '_adjacency.rda', outfile))

}

#############################################
########## 3. Get modules
#############################################

get_gene_modules <- function(infiles, outfile) {

    # Library
    require(WGCNA)

    # Load
    load(infiles[1])
    load(infiles[2])

    # Find modules
    initial_module_numbers <- cutreeDynamic(dendro = geneTree, distM = dissTOM, minClusterSize = 30, method = 'tree')
    initial_module_colors <- labels2colors(initial_module_numbers)

    # Calculate and cluster eigengenes
    MEList <- moduleEigengenes(log1p_dataframe, colors = initial_module_colors)
    MEDiss <- 1-cor(MEList$eigengenes)
    METree <- hclust(as.dist(MEDiss), method = "average")

    # Merge eigengenes
    MEDissThres <- 0.05
    merge <- mergeCloseModules(log1p_dataframe, initial_module_colors, cutHeight = MEDissThres, verbose = 3)
    colorOrder <- c("grey", standardColors(50));
    merged_module_colors <- merge$colors
    merged_module_numbers <- match(merged_module_colors, colorOrder)-1;
 
    # Modules
    module_dataframe <- data.frame(gene_id=rownames(dissTOM), module_name=paste0('module_', merged_module_numbers), module_color=merged_module_colors)

    # Modules
    color_dataframe <- module_dataframe %>% select(module_name, module_color) %>% mutate(module_number=as.numeric(gsub('.*_(.*)', '\\1', module_name))) %>% distinct %>% arrange(module_number)

    # Eigengenes
    me_dataframe <- merge$newMEs
    new_module_names <- color_dataframe %>% mutate(old_name=paste0('ME', module_color)) %>% pull(module_name, old_name)
    colnames(me_dataframe) <- new_module_names[colnames(me_dataframe)]
    me_dataframe <- me_dataframe[,new_module_names] %>% rownames_to_column('sample_name')

    # Save
    save(module_dataframe, color_dataframe, me_dataframe, merge, initial_module_colors, merged_module_colors, file=outfile)
   
    # Plot
    pdf(gsub('.rda', '.pdf', outfile), onefile=TRUE, height=4, width=13)

    # Initial clusters
    plotDendroAndColors(geneTree, initial_module_colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
    
    # Plot module clustering
    plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
    abline(h=MEDissThres, col = "red")
    
    # Merged clusters
    plotDendroAndColors(geneTree, cbind(initial_module_colors, merged_module_colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
    dev.off()

}

#############################################
########## 4. Get enrichment
#############################################

run_module_enrichment <- function(infile, outfile) {

    # Load
    load(infile)

    # Get modules
    module_list <- split(module_dataframe$gene_id, module_dataframe$module_name)

    # Run enrichment
    gprofiler_results <- gprofiler2::gost(query = module_list, multi_query = FALSE)

    # Save
    save(gprofiler_results, file=outfile)

    # Log
    writeLines(capture.output(sessionInfo()), gsub('.rda', '.log', outfile))

}

#############################################
########## 5. Get module preservation
#############################################

get_module_preservation <- function(infiles, outfile) {

    # Library
    suppressPackageStartupMessages(require(WGCNA))
    suppressPackageStartupMessages(require(DESeq2))
    
    # Load data
    load(infiles[1])
    load(infiles[2])
    load(infiles[3])

    # Get size factors
    dds <- estimateSizeFactors(dds_list[['gene']])

    # Get expression
    test_log1p_dataframe <- log10(counts(dds, normalized=TRUE)+1) %>% t

    # Remove genes with no variance in test dataset
    gene_variance <- apply(test_log1p_dataframe, 2, var)
    variable_genes <- names(gene_variance)[gene_variance > 0]

    # Intersect with network expression data
    filtered_genes <- intersect(colnames(log1p_dataframe), variable_genes)

    # Get matching colors
    modules <- module_dataframe %>% pull(module_name, gene_id) 
    filtered_modules <- modules[filtered_genes]

    # Create expression list
    multiExpr <- list(
        'network_data' = list('data' = log1p_dataframe[,filtered_genes]),
        'test_data' = list('data' = test_log1p_dataframe[,filtered_genes])
    )

    # Enable threads
    enableWGCNAThreads(10)

    # Get preservation
    module_preservation <- modulePreservation(
        multiData = multiExpr,
        multiColor = list('network_data'=filtered_modules),
        referenceNetworks = 1,
        networkType = 'signed',
        maxModuleSize = 1000,
        nPermutations = 1000,
        parallelCalculation = TRUE,
        verbose = 3
    )

    # Save
    save(module_preservation, file=outfile)

}

#############################################
########## 6. Get gene connectivity
#############################################

get_module_membership <- function(infiles, outfile) {
    
    # Library
    require(WGCNA)

    # Load
    load(infiles[1])
    load(infiles[2])
    load(infiles[3])

    # Get module membership
    membership_dataframe <- signedKME(log1p_dataframe,  me_dataframe %>% column_to_rownames('sample_name'))#, corOptions="use = 'p', method = 'spearman'")

    # Get connectivity
    module_assignments <- module_dataframe %>% pull(module_name, gene_id)
    connectivity_dataframe <- intramodularConnectivity(adjacency, module_assignments[colnames(adjacency)])

    # Save
    save(membership_dataframe, connectivity_dataframe, file=outfile)

}

#############################################
########## 7. Get gene networks
#############################################

get_gene_networks <- function(infiles, outfile) {
    
    # Load
    load(infiles[1])
    load(infiles[2])
    load(infiles[3])

    # Get gene symbol
    name_dataframe <- rtracklayer::readGFF(infiles[4]) %>% select(gene_id, gene_name) %>% distinct

    # Read genes
    developmental_dataframe <- fread(infiles[5], header=TRUE)

    # Sort
    colnames(membership_dataframe) <- gsub('kME', 'mo', colnames(membership_dataframe))
    gene_dataframe <- membership_dataframe %>% rownames_to_column('gene_id') %>% pivot_longer(-gene_id, names_to = 'module_name', values_to = 'module_membership') %>% inner_join(module_dataframe, by=c('gene_id', 'module_name'))
    id2symbol <- name_dataframe %>% pull(gene_name, gene_id)

    # # Top N genes by connectivity
    selected_modules <- paste0('module_', c(2,6,3,8,13,32,40))
    # nr_genes <- as.numeric(gsub('.*top(.*?)/.*', '\\1', outfile))
    # top_dataframe <- gene_dataframe %>% group_by(module_name) %>% slice_max(order_by = module_membership, n=nr_genes) %>% filter(module_name %in% selected_modules)
    
    # Select top 3 novel genes per cluster
    novel_dataframe <- gene_dataframe %>% filter(grepl('TALON', gene_id)) %>% group_by(module_name) %>% slice_max(order_by = module_membership, n = 3) %>% filter(module_name %in% selected_modules)

    # Select 10 known genes per cluster, prioritizing developmental ones
    known_dataframe <- gene_dataframe %>% filter(grepl('ENS', gene_id)) %>% mutate(gene_symbol=id2symbol[gene_id], developmental=ifelse(gene_symbol %in% developmental_dataframe$Gene_Symbol, 100, 1)) %>% group_by(module_name) %>% 
        slice_max(order_by = developmental*module_membership, n = 15) %>% filter(module_name %in% selected_modules) %>% arrange(module_name)

    # Merge
    top_dataframe <- rbind(novel_dataframe, known_dataframe) %>% select(-gene_symbol, -developmental)

    # Split
    gene_ids <- setNames(top_dataframe %>% group_split, top_dataframe %>% group_keys %>% pull(module_name)) %>% lapply(function(x) x %>% pull(gene_id))

    # Get correlations
    correlations <- lapply(gene_ids, function(x) {
        correlation_matrix <- cor(log1p_dataframe[,x], method='spearman')
        colnames(correlation_matrix) <- id2symbol[x]
        rownames(correlation_matrix) <- id2symbol[x]
        correlation_matrix[!upper.tri(correlation_matrix)] <- NA
        edge_dataframe <- correlation_matrix %>% as.data.frame %>% rownames_to_column('source_gene_id') %>% pivot_longer(-source_gene_id, names_to = 'target_gene_id', values_to = 'correlation') %>% drop_na
        colnames(edge_dataframe) <- gsub('_gene_id', '', colnames(edge_dataframe))
        edge_dataframe
    })
                                                                                                                    
    # Write edges
    for (module_name in names(correlations)) {
        module_outfile <- gsub('network-nodes', paste0(module_name, '-network'), outfile)
        fwrite(correlations[[module_name]], file=module_outfile, sep='\t')
    }

    # Write nodes
    fwrite(name_dataframe %>% filter(gene_id %in% do.call('c', gene_ids)) %>% rename('shared_name'='gene_id', 'name'='gene_name'), file=outfile, sep='\t')

}

# #############################################
# ########## 5. Get module correlations
# #############################################

# get_module_correlations <- function(infiles, outfile) {

#     # Load
#     load(infiles[1])

#     # Read expression
#     normalized_count_dataframe <- fread(infiles[2]) %>% filter(grepl('TALON', gene_id)) %>% column_to_rownames('gene_id') %>% as.matrix

#     # Normalize expression
#     log1p_matrix <- log10(normalized_count_dataframe+1)

#     # Get eigengenes
#     eigengene_matrix <- me_dataframe %>% column_to_rownames('sample_name') %>% as.matrix %>% t

#     # Get common samples
#     common_samples <- intersect(colnames(log1p_matrix), colnames(eigengene_matrix))
#     log1p_matrix <- log1p_matrix[,common_samples]
#     eigengene_matrix <- eigengene_matrix[,common_samples]

#     # Correlate
#     correlation_results <- list()
#     for (i in 1:nrow(log1p_matrix)) {
#         for (j in 1:nrow(eigengene_matrix)) {
#             spearman_results <- suppressWarnings(cor.test(log1p_matrix[i,], eigengene_matrix[j,], method='spearman'))
#             correlation_results[[length(correlation_results)+1]] <- data.frame(gene_id=rownames(log1p_matrix)[i], module_name=rownames(eigengene_matrix)[j], rho=spearman_results$estimate[[1]], pvalue=spearman_results$p.value)
#         }
#     }

#     # Merge
#     correlation_dataframe <- correlation_results %>% bind_rows %>% arrange(pvalue)

#     # Write
#     fwrite(correlation_dataframe, file=outfile, sep='\t')

# }

#######################################################
#######################################################
########## S9. SUPPA
#######################################################
#######################################################

#############################################
########## 6. Cluster PSI
#############################################

cluster_psi <- function(infiles, outfile) {
    
    # Library
    require(Mfuzz)
    
    # Read PSI
    psi_dataframe <- fread(infiles[1]) %>% column_to_rownames('V1')
    
    # Read differential results
    suppa_dataframe <- lapply(infiles[2:length(infiles)], function(x) {
        fread(x) %>% mutate(comparison=gsub('.*/human-(.*)-.*.tsv', '\\1', x))
    }) %>% bind_rows
    
    # Get significant events
    significant_events <- suppa_dataframe %>% filter(pval < 0.05 & abs(dPSI) > 0.1) %>% pull(Event_id) %>% unique

    # Filter
    filtered_dataframe <- psi_dataframe[significant_events,]

    # Sample dataframe
    sample_dataframe <- data.frame(sample_name=colnames(filtered_dataframe)) %>% mutate(cell_type=gsub('2PN', '1C', gsub('human_(.*?)_.*', '\\1', sample_name)))

    # Get average
    average_dataframe <- filtered_dataframe %>% rownames_to_column('event_id') %>% pivot_longer(-event_id, names_to = 'sample_name', values_to = 'psi') %>% left_join(sample_dataframe, by='sample_name') %>% 
        group_by(event_id, cell_type) %>% summarize(mean_psi=mean(psi, na.rm=TRUE)) %>% pivot_wider(id_cols = 'event_id', names_from = 'cell_type', values_from = 'mean_psi') %>% column_to_rownames('event_id')

    # Get times
    times <- c('1C'=0, '2C'=1, '4C'=2, '8C'=3, 'morula'=4, 'blastocyst'=5)

    # Timepoint dataframe
    timepoint_dataframe <- data.frame(sample_name=colnames(average_dataframe)) %>% rowwise %>% mutate(time=times[sample_name]) %>% column_to_rownames('sample_name')

    # Convert
    timepoint_annotation <- AnnotatedDataFrame(data=timepoint_dataframe)

    # Variable metadata
    varMetadata <- data.frame(labelDescription='Time')
    rownames(varMetadata) <- 'time'

    # Create expression set
    eset <- ExpressionSet(assayData=as.matrix(average_dataframe), phenoData = timepoint_annotation, varMetadata = varMetadata)

    # Filter
    eset.r <- filter.NA(eset, thres=0.05)

    # Replace NA
    eset.f <- fill.NA(eset.r, mode="mean") #knnw

    # Filter
    eset.f2 <- filter.std(eset.f,min.std=0)

    # Standardize
    eset.s <- standardise(eset.f2)

    # Estimate
    m1 <- mestimate(eset.s)

    # Set cluster range
    cluster_ranges <- seq(5, 30, by=1)

    # Selection
    pdf(gsub('.rda', '.pdf', outfile), height=5, width=9)
    dmin_results <- Dmin(eset.s, m=m1, crange=cluster_ranges, repeats=5, visu=TRUE)
    cselection_results <- cselection(eset.s, m=m1, crange=cluster_ranges, repeats=5, visu=TRUE)
    dev.off()

    # Save
    save(eset.s, m1, cluster_ranges, dmin_results, cselection_results, times, file=outfile)
}

#############################################
########## 7. Get PSI clusters
#############################################

get_psi_clusters <- function(infile, outfile) {
    
    # Library
    require(Mfuzz)

    # Load
    load(infile)

    # Cluster
    cluster_results <- mfuzz(eset.s, c=10, m=m1)
    
    # Selection
    pdf(gsub('.rda', '.pdf', outfile), height=9, width=15)
    mfuzz.plot(eset.s,cl=cluster_results,mfrow=c(4,4), new.window = FALSE, time.labels=names(times))
    dev.off()

    # Save
    save(cluster_results, file=outfile)
}

#######################################################
#######################################################
########## S10. Isoform switching
#######################################################
#######################################################

#############################################
########## 1. Filter
#############################################

filter_isoform_data <- function(infiles, outfile, file_type) {

    # Get transcript ids
    transcript_ids <- rtracklayer::readGFF(infiles[2]) %>% pull(transcript_id) %>% unique

    # Read and filter
    if (file_type == 'gtf_cds') {
        
        # Read GTF
        gtf <- rtracklayer::import(infiles[1])
        
        # Filter
        gtf_filtered <- gtf[gtf$transcript_id %in% transcript_ids]

        # Export
        rtracklayer::export(gtf_filtered, outfile, format='gtf')

    } else if (file_type == 'cpat_predictions') {
        
        # Read CPAT
        cpat_dataframe <- fread(infiles[1]) %>% filter(`Sequence Name` %in% transcript_ids) %>% mutate(`Data ID`=1:n()-1)
        
        # Write
        fwrite(cpat_dataframe, file=outfile, sep='\t')
        
    } else if (file_type == 'pfam_predictions') {
            
        # Read Pfam
        pfam_dataframe <- fread(infiles[1]) %>% filter(seq_id %in% transcript_ids)
        
        # Write
        fwrite(pfam_dataframe, file=outfile, sep='\t')
        
    } else if (file_type == 'transcript_fasta') {

        # Read FASTA
        talon_fasta <- Biostrings::readDNAStringSet('arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon.fasta', format="fasta")

        # Subset
        filtered_fasta <- talon_fasta[transcript_ids]

        # Write
        Biostrings::writeXStringSet(filtered_fasta, file=outfile)

    }
}


#############################################
########## 1. Load data
#############################################

load_isoform_data <- function(infiles, outfile) {

    # Library
    suppressPackageStartupMessages(require(IsoformSwitchAnalyzeR))
    suppressPackageStartupMessages(require(rjson))

    # Read metadata
    metadata_dataframe <- fread(infiles[5])

    # Read expression
    salmon_infiles <- metadata_dataframe %>% pull(files, names)
    salmonQuant <- importIsoformExpression(sampleVector = salmon_infiles)

    # Get organism
    organism <- gsub('.*.dir/(.*?)/.*', '\\1', infiles[1])

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
        pathToCPATresultFile = infiles[1],
        codingCutoff         = ifelse(organism=='human', 0.364, 0.44),
        removeNoncodinORFs   = FALSE   # because ORF was added from CPAT
    )

    # Add Pfam
    aSwitchList <- analyzePFAM(
        switchAnalyzeRlist   = aSwitchList,
        pathToPFAMresultFile = infiles[4],
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