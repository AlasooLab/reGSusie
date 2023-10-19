library(data.table)
library(optparse)
library(dplyr)
library(qqman)


# Argument parser
option_list <- list( 
    make_option('--phenotype_id', type='character',
    help='ID of the phenotype to be analysed.'),
    make_option('--regenie_file', type='character',
    help='List of the regenie step 2 outputs')
);

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
phenotype_id <- options$phenotype_id

full_regenie <- read.table(options$regenie_file, header=T)

full_regenie <- full_regenie[order(full_regenie$CHROM, full_regenie$GENPOS), ]
rownames(full_regenie) <- 1:nrow(full_regenie)  

write.table(full_regenie, gzfile(paste0(phenotype_id, '_full_regenie.tsv.gz')), row.names=F, sep='\t', quote=F) 

set.seed(42)  
threshold <- -log10(5e-08)

filtered_df <- full_regenie[full_regenie$LOG10P >= threshold, ]
num_rows_to_remove <- floor(0.9 * nrow(full_regenie))

poss_remove <- full_regenie[full_regenie$LOG10P < threshold, ]
rows_to_remove <- poss_remove[-sample(1:nrow(poss_remove), num_rows_to_remove), ]

final_df <- rbind(filtered_df, rows_to_remove)

final_df <- final_df[order(final_df$CHROM, final_df$GENPOS), ]

png(paste0(phenotype_id,'_Manhattan_plot.png'), width=1200, height=800, units='px')
manhattan(final_df, chr='CHROM', bp='GENPOS', p='LOG10P', snp='ID', logp=F)
dev.off()