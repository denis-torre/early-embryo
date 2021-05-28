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
########## 1. Split GTF by isoform class
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
        rtracklayer::export(gtf_split[[transcript_class]] %>% head(5000), gtf_outfile, format='gtf')
        
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