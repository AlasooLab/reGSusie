library(data.table)
library(optparse)
library(dplyr)
library(stringr) 
library(GenomicRanges)

find_lead_SNPs <- function(data, win_size, MHC_start, MHC_end, remove_MHC) {
  final <- data[-c(1:nrow(data)), ]
  
  while (nrow(data) > 0) {  
    lead_snp <- data[abs(data$log10P) == max(data$log10P), ] 
    final <- rbind(final, lead_snp)
    data <- data[!(data$chromosome == lead_snp$chromosome &
                   data$position > lead_snp$position-win_size &
                   data$position < lead_snp$position+win_size), ]
  }
  final$start <- final$position - win_size
  final$end <- final$position + win_size 
  
  if (remove_MHC == T) {
    final <- final[!(final$chromosome == 6 & ((final$start > MHC_start & final$start  < MHC_end) |
                  (final$end > MHC_start & final$end < MHC_end))), ]

  }
   
  ifelse (final$start < 0, 
          final$region <- paste0(final$chromosome, ':', 0, '-', final$end), 
          final$region <- paste0(final$chromosome, ':', final$start, '-', final$end))

  #final$region <- str_replace(final$region, '23:', 'X:')
  final <- final %>% arrange(chromosome, position)

  return(final)
}

# Argument parser
option_list <- list( 
    make_option('--regenie_out', type='character',
    help='Regenie step 2 output: GWAS summary statistics file.'),
    make_option('--phenotype_id', type='character',
    help='ID of the phenotype to be analysed.'),
    make_option('--chr', type='character',
    help='Chromosome to be analysed.'),
    make_option('--p_threshold', type='double', default=5e-08,
    help='GWAS p-value threshold. Default is 5e-8'),
    make_option('--window_size', type='integer',default=1500000, 
    help='Genomic window size for defining loci. Default is +/- 1.5 Mb from lead SNP.'),
    make_option('--MHC_start', type='integer', default=28510120,
    help='MHC region starting point. Defaultly hg38 assembly coordinates are used.'),
    make_option('--MHC_end', type='integer', default=33480577,
    help='MHC region ending point. Defaultly hg38 assembly coordinates are used.'),
    make_option('--remove_MHC', type='logical', default=T,
    help='Option to remove MHC region. Default is TRUE.')
);

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
p_thresh <- options$p_threshold
win_size <- options$window_size
MHC_start <- options$MHC_start
MHC_end <- options$MHC_end
remove_MHC <- options$remove_MHC

regenie_output <- read.table(gzfile(options$regenie_out), header=T)
phenotype_id <- options$phenotype_id
chr <- options$chr

formatted_file <- regenie_output[, c('ID', 'CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N', 'BETA', 'SE', 'LOG10P')]
colnames(formatted_file) <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'aaf', 'info', 'n', 'beta', 'se', 'log10P')
formatted_file <- formatted_file %>% mutate(maf=ifelse(aaf > 0.5, 1-aaf, aaf))
formatted_file <- formatted_file %>% mutate(an=2*n)
formatted_file <- formatted_file %>% mutate(ac=2*n*aaf)

#fwrite(formatted_file, file=paste0('standard_', phenotype_id, '_chr', chr), sep='\t', quote=F, row.names=F)

filter <- formatted_file %>% filter(log10P > -log10(p_thresh))
message(paste0('Number of statistically significant variants: ', nrow(filter)))

regions <- find_lead_SNPs(filter, win_size, MHC_start, MHC_end, remove_MHC)
phenotype <- ifelse(nrow(filter) == 0, '', phenotype_id)

gr_object <- GRanges(seqnames=regions$chromosome, ranges=IRanges(start=as.numeric(regions$start), end=as.numeric(regions$end)))
reduced_gr <- reduce(gr_object)
printable_regions <- data.frame(matrix(ncol=4, nrow=0))
colnames(printable_regions) <- c('chromosome', 'start', 'end', 'region')

if (length(reduced_gr) > 0) {
   for (i in 1:length(reduced_gr)) {
      printable_regions[nrow(printable_regions) + 1,] <- c(chr, start(ranges(reduced_gr))[i], end(ranges(reduced_gr))[i], 
                                                           paste0(chr, ':', start(ranges(reduced_gr))[i], '-', end(ranges(reduced_gr))[i]))
   }
}

#fwrite(data.table(phenotype=phenotype, region=printable_regions$region), 
#       file=paste0('chr', chr, '_', phenotype_id, '_regions.txt'), sep='\t', quote=F, row.names=F)


if (nrow(printable_regions) == 0) {
  fwrite(data.table(), paste0(phenotype_id, '_NULL_region.txt'), sep='\t', quote=F, row.names=F)
}

if (nrow(printable_regions) > 0) {
  for (i in 1:nrow(printable_regions)) {
    sub_reg <- formatted_file %>% filter(chromosome == as.numeric(printable_regions$chromosome[i]) & 
                                         position > as.numeric(printable_regions$start[i]) & 
                                         position < as.numeric(printable_regions$end[i]))
    
    if (nrow(sub_reg) <= 50) {
      fwrite(data.table(), paste0(phenotype_id, '_NULL_region.txt'), sep='\t', quote=F, row.names=F)
    }

    if (nrow(sub_reg) > 50) {     
      sub_reg <- sub_reg[!duplicated(sub_reg[, 1]), ]
      sub_reg <- sub_reg %>% mutate(chromosome = ifelse(chromosome < 10, paste0('0', chromosome), chromosome))
      fwrite(sub_reg, paste0(phenotype_id, '_', printable_regions$region[i], '_region.txt'), sep=' ', quote=F, row.names=F)

    }
  }
}



