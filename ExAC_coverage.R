## This scripts aims to assess ExAC coverage from public sources.
## It uses only canonical transcript UTRs and also subtracts overlapping coding regions.

#load packages
library(tidyverse)

# Import data
# From Ensembl API using script canonical.pl

canonical_transcripts <- read_tsv("canonical_transcripts/human_canonical_transcripts_v93.fa", 
                                  col_names = c("transcript"), 
                                  col_types = "c")


#import ExAC coverage
#files are downloaded from ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/coverage

chrs <-c(as.character(1:22), "X", "Y")
ExAC_list = list()

for (i in 1:length(chrs)){
  file=gzfile(sprintf("ExAC_coverage/Panel.chr%s.coverage.txt.gz", chrs[i]))
  ExAC_list[[i]] <- read_tsv(file, skip = 1, col_names = c("chr", "pos", "mean", "median"), 
                             col_types = "cinn_________")
}

ExACdf <- do.call("bind_rows", ExAC_list)

rm(ExAC_list)


## Import position dataframe from Position_dfs.R script

posdf <- read_tsv('position_dataframe_UTRs_CDS.tsv', col_types = cols(chr = 'c'))


## Annotate positions with Exac coverage

pos_cov <- posdf %>% left_join(ExACdf)

## generate position dataframe for all CDS regions

gtf <- as_tibble(rtracklayer::import('refs/gencode.v19.annotation.gtf'))

allCDS <- gtf %>% filter(type == 'CDS') %>% distinct(seqnames, start, end, .keep_all = T)

seqfun <- function (a, b, c){
  pos <- b:c
  chr <- rep(as.character(a), length(pos))
  data.frame(chr, pos, stringsAsFactors = F)
}

allCDSposdf <- bind_rows(mapply(seqfun, 
                                allCDS$seqnames, 
                                allCDS$start, 
                                allCDS$end, 
                                SIMPLIFY = F)) 

allCDSposdf <- allCDSposdf %>% distinct(chr, pos)

allCDSposdf$chr <- sub("chr", "", allCDSposdf$chr)

#extraxt all CDS positions from

pos_cov_woCDS <- pos_cov %>% anti_join(allCDSposdf)

#plots
##Plots
# ggplot(filter(pos_cov_woCDS, !is.na(mean)), aes(x=pos_relative)) +
#   geom_bar() +
#   coord_cartesian(xlim = c(-200,200)) +
#   ggtitle("ExAC coverage for UTR regions by relative positions") +
#   ylab("Count of UTRs for each relative position") +
#   xlab("Relative position in UTR from the first/last exon") +
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5))


possum <- summarize(group_by(filter(pos_cov_woCDS, !is.na(mean)), pos_relative), mean_coverage = mean(mean))

ExAC_covplot <- ggplot(data=possum, aes(x=pos_relative, y=mean_coverage)) +
  geom_bar(stat="identity") +
  ylab("Mean coverage") +
  xlab("Distance from start (UTR5) or stop (UTR3) codon") +
  coord_cartesian(xlim = c(-200,200), ylim = c(0,50)) +
  theme_bw() +
  theme(axis.text=element_text(size=10), 
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        legend.key.size = unit(2, 'line' ),
        #legend.position = c(0.75,0.9),
        #legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line')) +  annotate('text', x = -150, y = 45, label = 'UTR5', size = 8) +
  annotate('text', x = 150, y = 45, label = 'UTR3', size = 8)

ggsave("ExAC_coverage.pdf", plot = ExAC_covplot, width = unit(6, 'in'), height  = unit(4, 'in'))
