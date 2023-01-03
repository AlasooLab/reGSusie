#!/usr/bin/env Rscript

# Installing the needed packages
packages <- c('optparse', 'dplyr', 'ggplot2', 'ggrepel')

for (package in packages) {
  if(!require(package, character.only = T)) {
    install.packages(package, dependencies = T, repos = c(CRAN = "https://cran.r-project.org"))
  }
}

# Defining the plotting style
theme_signature <- function() {
  theme_minimal(base_size=14, base_family="Rubik") +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line=element_line(size=.13, color='#343434'),
      axis.text=element_text(color='#343434'),
      axis.ticks.x=element_line(size=0.5, color='#343434'),
      axis.ticks.y=element_line(size=0.5, color='#343434'),
      axis.title=element_text(color="#343434"),
      axis.title.y=element_text(hjust=1, margin=margin(0, 6, 0, 15, "pt")),
      axis.title.x=element_text(hjust=0, margin=margin(6, 0, 15, 0, "pt")),
      plot.subtitle=element_text(color="#343434", size=11),
      plot.title=element_text(color="#343434", size=20),
      plot.title.position="plot",
      plot.caption=element_text(hjust=0, color="#343434"),
      plot.caption.position="plot",
      plot.margin=margin(.5, .5, .5, .5, "cm"),
      plot.background=element_rect(fill="#f7f7f7"),
      legend.text=element_text(color="#343434"),
      legend.title=element_text(color="#343434"),
      strip.text=element_text(color="black")) 
}

# Reading the data and modifying it for ploting
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--phenotype_id"), type="character", default=NULL,
              help="Map file for plot title in Manhattan plot [default= %default]", metavar="character")
);

option_parser = OptionParser(option_list=option_list);
option = parse_args(option_parser, positional_arguments=0);

file <- option$options$file
print(paste("reading file:", file))
data <- read.table(gzfile(file), header=T)

phenotype_id <- option$options$phenotype_id

output_prefix = file
if(!is.null(option$options$out)) {
  output_prefix = option$options$out
}


plot_data <- data %>% 
	     group_by(CHROM) %>% 
  	     summarise(chr_len=max(GENPOS)) %>% 
   	     mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
             select(-chr_len) %>%
	     left_join(data, ., by=c("CHROM"="CHROM")) %>%
  	     arrange(CHROM, GENPOS) %>%
             mutate(BPcum=GENPOS+tot) %>%
	     mutate(annot=ifelse(LOG10P>-log10(5e-08), 'T', 'F'))

axisdf = plot_data %>% 
	 group_by(CHROM) %>%
	 summarize(center=(max(BPcum)+min(BPcum))/2)


options(bitmapType='cairo')
png(paste(output_prefix, '_manhattan.png', sep=''), width=2000, height=1000)
ggplot(plot_data, aes(x=BPcum, y=LOG10P)) +
    geom_point(aes(color=as.factor(CHROM)), alpha=0.7, size=1.5) +
    scale_color_manual(values = rep(c('#26abff', '#929497'), 22 )) +
    scale_x_continuous(label = axisdf$CHROM, breaks=axisdf$center) +
    xlab('Chromosome') +
    scale_y_continuous(expand = c(0, 0), limits=c(0, max(plot_data$LOG10P))) +
    ylab(expression(-log[10](p))) +
    geom_text_repel(data=subset(plot_data, annot=='T'), 
                    aes(label=ID, color=as.factor(CHROM)), size=2) +
    geom_hline(yintercept=-log10(5e-08), color="#dd1c77") +
    geom_hline(yintercept=-log10(1e-05), color="#dd1c77", linetype='dashed') +
    theme_signature() +
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

dev.off()
