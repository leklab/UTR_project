library(tidyverse)

#import position dataframe

posdf <- read_tsv('position_dataframe_UTRs_CDS.tsv', col_types = cols(chr = 'c'))

#Extract Kozak as defined by first 4 nucleotides from CDS and last 6 from UTR5

KOZposdf <- posdf %>% filter(Region %in% c('UTR5', 'CDS') & pos_relative %in% -6:4 & chr != 'M')

rm(posdf)

#import fasta

library(seqinr)

fasta <- read.fasta('/Users/sander/data_references/Homo_sapiens_assembly19.fasta', 'DNA')

pos_nuc <- function(chr, pos){
  fasta[[chr]][pos]
}

KOZposdf$nucleotide <- toupper(mapply(pos_nuc, KOZposdf$chr, KOZposdf$pos))

KOZposdf <- KOZposdf %>% mutate(mRNA = case_when(
  strand == '+' ~ nucleotide,
  nucleotide == 'A' ~ 'T',
  nucleotide == 'T' ~ 'A',
  nucleotide == 'C' ~ 'G',
  nucleotide == 'G' ~ 'C',
  TRUE ~ 'NA'
))

library(seqLogo)

nuc_tbl <- table(KOZposdf$mRNA, KOZposdf$pos_relative)

nuc_mat <- matrix(prop.table(nuc_tbl, 2), nrow = 4)

p <- makePWM(nuc_mat)

seqLogo(p, ic.scale=T,  xaxis = F, yaxis = F)



