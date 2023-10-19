library(data.table)
library(optparse)
library(dplyr)
library(stringr) 
library(GenomicRanges)

find_regions <- function(data, win_size, MHC_start, MHC_end, remove_MHC, win_shrink_ratio, max_region_width) {
  final <- data[-c(1:nrow(data)), ]
  data_copy <- data
  
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
  
  final$region <-  ifelse (final$start < 0, paste0(final$chromosome, ':', 0, '-', final$end),
                           paste0(final$chromosome, ':', final$start, '-', final$end))

  final <- final %>% arrange(chromosome, position)

  gr_object <- GRanges(seqnames=final$chromosome, ranges=IRanges(start=final$start, end=final$end))
  reduced_gr <- reduce(gr_object)

  regions_to_keep <- reduced_gr[width(reduced_gr) <= max_region_width]
  regions_to_shrink <- reduced_gr[width(reduced_gr) > max_region_width]

  data_gr <- with(data_copy, GRanges(seqnames=chromosome, ranges=IRanges(start=position, end=position)))
  overlap_indexes <- findOverlaps(data_gr, regions_to_shrink)
  subset_data <- data_copy[queryHits(overlap_indexes), ]

  if (nrow(subset_data) > 0) {
    regions_to_keep_new <- find_regions(subset_data, win_size=win_size * win_shrink_ratio,
                                          MHC_start=MHC_start, MHC_end=MHC_end, remove_MHC=remove_MHC,
                                          win_shrink_ratio=win_shrink_ratio, max_region_width=max_region_width)
    final_regions <- c(regions_to_keep, regions_to_keep_new)
  } else {
    final_regions <- regions_to_keep
  }

  return(sort(final_regions))
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
    make_option('--window_shrink_ratio', type='double',default=0.95, 
    help='Ratio by which the window size is reduced in each iteration. Default is 0.95.'),
    make_option('--max_region_width', type='integer',default=4500000, 
    help='Maximum allowed width for a genomic region. Regions wider than this value will be further reduced. Default is 4.5 Mb.'),
    make_option('--MHC_start', type='integer', default=28510120,
    help='MHC region starting point. Defaultly hg38 assembly coordinates are used.'),
    make_option('--MHC_end', type='integer', default=33480577,
    help='MHC region ending point. Defaultly hg38 assembly coordinates are used.'),
    make_option('--remove_MHC', type='logical', default=T,
    help='Option to remove MHC region. Default is TRUE.'),
    make_option('--bgen_chr_has_zero', type='logical', default=F,
    help='Specify if chromosome numbers (1...9) in BGEN file have leading zeros (i.e. 01...09). Default is FALSE.')
);

option_parser <- OptionParser(option_list=option_list)
options <- parse_args(option_parser)
p_thresh <- options$p_threshold
win_size <- options$window_size
win_shrink_ratio <- options$window_shrink_ratio
max_region_width <- options$max_region_width
MHC_start <- options$MHC_start
MHC_end <- options$MHC_end
remove_MHC <- options$remove_MHC
bgen_chr_has_zero <- options$bgen_chr_has_zero

regenie_output <- read.table(gzfile(options$regenie_out), header=T)
phenotype_id <- options$phenotype_id
chr <- options$chr

formatted_file <- regenie_output[, c('ID', 'CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N', 'BETA', 'SE', 'LOG10P')]
colnames(formatted_file) <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'aaf', 'info', 'n', 'beta', 'se', 'log10P')
formatted_file <- formatted_file %>% mutate(maf=ifelse(aaf > 0.5, 1-aaf, aaf))
formatted_file <- formatted_file %>% mutate(an=2*n)
formatted_file <- formatted_file %>% mutate(ac=2*n*aaf)

filter <- formatted_file %>% filter(log10P > -log10(p_thresh))
message(paste0('Number of statistically significant variants: ', nrow(filter)))

regions <- find_regions(filter, win_size, MHC_start, MHC_end, remove_MHC, win_shrink_ratio, max_region_width)
phenotype <- ifelse(nrow(filter) == 0, '', phenotype_id)

printable_regions <- data.frame(matrix(ncol=4, nrow=0))
colnames(printable_regions) <- c('chromosome', 'start', 'end', 'region')

if (length(regions) > 0) {
   for (i in 1:length(regions)) {
      printable_regions[nrow(printable_regions) + 1,] <- c(chr, ifelse(start(ranges(regions))[i] < 0, 0, start(ranges(regions))[i]),
                                                           end(ranges(regions))[i], 
                                                           paste0(chr, ':', 
                                                                  ifelse(start(ranges(regions))[i] < 0, 0, start(ranges(regions))[i]),
                                                                  '-', end(ranges(regions))[i]))
   }
}

if (nrow(printable_regions) == 0) {
  fwrite(data.table(), paste0(phenotype_id, '_NULL_region.txt'), sep='\t', quote=F, row.names=F)
} else {
  for (i in 1:nrow(printable_regions)) {
    if(printable_regions$chromosome[i] == 'X') {
      sub_reg <- formatted_file %>% filter(chromosome == 23 & 
                                         position > as.numeric(printable_regions$start[i]) & 
                                         position < as.numeric(printable_regions$end[i]))
      sub_reg$chromosome <- 'X'
      
    } else {
      sub_reg <- formatted_file %>% filter(chromosome == as.numeric(printable_regions$chromosome[i]) & 
                                         position > as.numeric(printable_regions$start[i]) & 
                                         position < as.numeric(printable_regions$end[i]))
    }
    if (nrow(sub_reg) <= 50) {
      fwrite(data.table(), paste0(phenotype_id, '_NULL_region.txt'), sep='\t', quote=F, row.names=F)
    } else {     
      sub_reg <- sub_reg[!duplicated(sub_reg[, 1]), ]
      if (bgen_chr_has_zero) {
        sub_reg <- sub_reg %>% mutate(chromosome = ifelse(chromosome != 'X' & chromosome < 10, paste0('0', chromosome), chromosome))
      }
      fwrite(sub_reg, paste0(phenotype_id, '_', printable_regions$region[i], '_region.txt'), sep=' ', quote=F, row.names=F) 
    }
  }
}

