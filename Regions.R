#Get all protein coding and canonical transcripts lines from gtf and separate into UTR3s, UTR5s and CDS regions

## Load libraries

#install.packages("rtracklayer")
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")

library(tidyverse)

#Import gtf from gencode.

gtf <- as_tibble(rtracklayer::import('refs/gencode.v19.annotation.gtf'))

#Import canonical transcripts

# Output from canonical.pl script.

canonical_transcripts <- read_tsv("canonical_transcripts/human_canonical_transcripts_v93.fa", 
                                  col_names = c("transcript"), 
                                  col_types = "c")

# select protein coding canonical trancsripyts from gtf
# seprate transcript version number for easier matching

gtf_v2 <- gtf %>% separate(transcript_id, c("transcript", "version"), sep = "\\.")
canonical_transcripts_v2 <- canonical_transcripts %>% separate(transcript, c("transcript", "version"), sep = "\\.")

#write_tsv(canonical_transcripts_v2$transcript, canonical_tr_wo_versions.txt)

prot_coding_canonical_v2 <- gtf_v2 %>% 
  filter(transcript_type == "protein_coding") %>% 
  filter(transcript %in% canonical_transcripts_v2$transcript)

#Select unique transcripts

transcripts <- unlist(distinct(prot_coding_canonical_v2, transcript))

# get a data frame indicating first and last CDS position

cdslist <- list()

for (i in 1:length(transcripts)) {
  tr <- transcripts[[i]]
  cds <- prot_coding_canonical_v2 %>% filter(type == 'CDS' & transcript == tr) %>% arrange(start)
  first_cds <- if (cds$strand[1] == '+') first(cds$start) else last(cds$end)
  last_cds <- if (cds$strand[1] == '+') last(cds$end) else first(cds$start)
  cdslist[[i]] <- c(transcript = tr, first_cds = first_cds, last_cds = last_cds)
}

cds_df <- do.call("bind_rows", cdslist)
cds_df$first_cds <- as.integer(cds_df$first_cds)
cds_df$last_cds <- as.integer(cds_df$last_cds)

# joni first and last CDS positions to gtf dataframe and annotate UTR3s and UTR5s
prot_coding_canonical_v3 <- prot_coding_canonical_v2 %>% left_join(cds_df) %>% mutate(UTRtype = case_when(
  type == 'UTR' & strand == '+' & end < first_cds ~ 'UTR5',
  type == 'UTR' & strand == '+' & start > last_cds ~ 'UTR3',
  type == 'UTR' & strand == '-' & end < last_cds ~ 'UTR3',
  type == 'UTR' & strand == '-' & start > first_cds ~ 'UTR5',
  TRUE ~ as.character(type)))


#### Generate GATK style intervals for every region of interest.

types <- rep(c("CDS", "start_codon", "stop_codon", "UTR3", "UTR5"), each = 3)
strands <- rep(c('+', '-', 'all'), length(unique(types)))

create_interval_files <- function(type, strand, df) {
  #discard chrM and filter for type, and select first cols
  type_df <- df %>% 
    filter(seqnames != 'chrM')  %>% 
    filter(UTRtype == !!type) 
  
  # remove 'chr' for b37 ref compatability
  type_df$seqnames <- sub("chr", "", type_df$seqnames)
  
  #filter for pluss and minus strands
  
  if (strand == '+') {
    type_df <- type_df %>% filter(strand == '+') %>% select(seqnames, start, end)
    prefix=paste(type, '_plus.intervals', sep = "")
  } else if (strand == '-') {
    type_df <- type_df %>% filter(strand == '-') %>% select(seqnames, start, end)
    prefix=paste(type, '_minus.intervals', sep = "")
  } else {
    type_df <- type_df %>% select(seqnames, start, end)
    prefix=paste(type, '.intervals', sep = "")
  }

  # format as intervals
  
  type.intervals <- type_df %>% 
    unite(chrstart, seqnames, start, sep = ":") %>%
    unite(chrstartend, chrstart, end, sep = "-")
  
  #write intervals file
  write_tsv(type.intervals, prefix, col_names = F)

}

for (i in 1:length(types)) {
  create_interval_files(type = types[i], strand = strands[i], df = prot_coding_canonical_v3)
}

## Now you have 15 files in your directory which serve as intervals for GATK select variants tool.

# save gtf dataframe to be used as an input for other scripts.

# write_tsv(prot_coding_canonical_v3, 'gtf_df_wUTRs_separated.tsv')


