# Load packages

library(tidyverse)

# Import modified gtf with UTR5s and UTR3s separated in UTRtype column. It was generated in Regions.R script

gtf_utrs <- read_tsv('gtf_df_wUTRs_separated.tsv')

# define function for generating relative positions inside the region (e.g. last UTR5 position is -1, and first CDS position 1)

pos_df_fun <- function(region, df, isUTR5 = F) {
  
  region_df <- df %>% filter(UTRtype == !!region)

  uniq_tr <- list() 
  transcripts <- unlist(distinct(region_df, transcript))
  
  for (i in 1:length(transcripts)) {
    t <- region_df %>% filter(transcript == transcripts[i]) %>% arrange(start)
    pos <- unlist(mapply(seq, t$start, t$end))
    
    if (isUTR5) {
      if (t$strand[1] == "+") {
          pos_relative <- -(length(pos)):-1
        } else {
          pos_relative <- -1:-(length(pos))
        }
    } else {
        if (t$strand[1] == "+") {
          pos_relative <- 1:length(pos)
        } else {
          pos_relative <- length(pos):1
        }
    } 
     
      uniq_tr[[i]] <- tibble(chr = unlist(t[1,1]), 
                                   start = unlist(t[1,2]), 
                                   end = tail(unlist(t[,3]),1), 
                                   strand = unlist(t[1,5]), 
                                   transcript = unlist(t[1,11]),
                                   version = unlist(t[1,12]),
                                   exons = dim(t)[1],
                                   genomic_posvec = list(pos),
                                   relative_posvec = list(pos_relative))
  }                
  
  uniq_tr <- do.call("bind_rows", uniq_tr)
  
  #create position data frames (has all positions in the dataset )
  
  poslist <- list()
  
  for (i in 1:dim(uniq_tr)[1]) {
    poslist[[i]] <- tibble(chr = uniq_tr$chr[i],
                               pos = unlist(uniq_tr$genomic_posvec[i]),
                               pos_relative = unlist(uniq_tr$relative_posvec[i]),
                               transcript = uniq_tr$transcript[i],
                               transcript_ver = uniq_tr$version[i],
                               strand = uniq_tr$strand[i],
                               Region = region)
  }
  
  do.call("bind_rows", poslist)
  
}

UTR5posdf <- pos_df_fun('UTR5', gtf_utrs, isUTR5 = T)
CDSposdf <- pos_df_fun('CDS', gtf_utrs, isUTR5 = F)
UTR3posdf <- pos_df_fun('UTR3', gtf_utrs, isUTR5 = F)

posdf <- bind_rows(UTR5posdf, CDSposdf, UTR3posdf) 

posdf$chr <- sub("chr", "", posdf$chr)

write_tsv(posdf, 'position_dataframe_UTRs_CDS.tsv')
