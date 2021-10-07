#################################################################
#################################################################
############### Embryo methylation - R Support #################
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
########## S3. Bismark
#######################################################
#######################################################

#############################################
########## 1. Average replicates
#############################################

get_average_methylation <- function(infiles, outfile) {

    # Read data
    bedgraph_dataframe <- lapply(infiles, function(x) {
        fread(x, skip = 1, col.names = c('chr', 'start', 'end', 'pct')) %>% mutate(sample=basename(x))
    }) %>% bind_rows %>% pivot_wider(id_cols = c(chr, start, end), names_from = sample, values_from = pct, values_fill = 0)

    # Get average
    bedgraph_dataframe$average <- rowMeans(bedgraph_dataframe[,4:ncol(bedgraph_dataframe)])

    # Subset
    average_dataframe <- bedgraph_dataframe[,c('chr', 'start', 'end', 'average')]

    # Write
    fwrite(average_dataframe, file=outfile, sep='\t', col.names = FALSE)

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################