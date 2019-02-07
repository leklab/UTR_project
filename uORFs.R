library(tidyverse)

# Input was sup table 4 from Gerstein lab paper https://doi.org/10.1093/nar/gky188 

# The input to this lab was produced using following grep command.
# cat complete_active_uORF_predictions_hg19.tsv | grep 'uORF_ATG' | grep -Ff ../canonical_tr_wo_versions.txt > complete_active_uORF_predictions_hg19_canonicalATGonly.tsv
# canonical_tr_wo_versions.txt was producsed in Regions.R script, contains just 1 column with all canonical transcripts without version numbers.

#in this analysis only active uORFs (active defined in original paper) on canonical transcript and with ATG start codons are used.

uorfs <- read_tsv('Gerstein_uORFs/complete_active_uORF_predictions_hg19_canonicalATGonly.tsv', col_names = c('uORF_ID', 'uORF_score', 'start_codon', 'chromosome', 'strand', 'start_coordinate', 'end_coordinate', 'gene', 'additional_transcripts', 'type', 'peptide_score'))

uorfs$chromosome <- sub("chr", "", uorfs$chromosome)

uorfs_plus <- uorfs %>% filter(strand == '+')
uorfs_minus <- uorfs %>% filter(strand == '-')


seqfun <- function (a, b, c, d, e, strand){
  pos <- b:c
  chr <- rep(as.character(a), length(pos))
  codon <- rep(as.character(d), length(pos))
  id <- rep(as.character(e), length(pos))
  
  data.frame(chr, pos, codon, id, stringsAsFactors = F)
}

plus_uorfs_pos_df <- bind_rows(mapply(seqfun, 
                                      uorfs_plus$chromosome, 
                                      uorfs_plus$start_coordinate, 
                                      uorfs_plus$end_coordinate, 
                                      uorfs_plus$start_codon, 
                                      uorfs_plus$uORF_ID, 
                                      strand = "+", SIMPLIFY = F)) %>% 
  mutate(strand = '+')



minus_uorfs_pos_df <- bind_rows(mapply(seqfun, 
                                       uorfs_minus$chromosome, 
                                       uorfs_minus$end_coordinate, 
                                       uorfs_minus$start_coordinate, 
                                       uorfs_minus$start_codon, 
                                       uorfs_minus$uORF_ID, strand = "-", SIMPLIFY = F)) %>% 
  mutate(strand = '-') 

uorf_pos_df <- bind_rows(plus_uorfs_pos_df, minus_uorfs_pos_df) %>% arrange(chr, pos)


posdf <- read_tsv('position_dataframe_UTRs_CDS.tsv', col_types = cols(chr = 'c'))


# as uorf dataframe includes only intronic regions inside coordinates, need to extract those and assign relative positions in respect of mRNA sequence.
uorf_pos_df_exonic <- uorf_pos_df %>% semi_join(filter(posdf, Region %in% c('UTR5', 'CDS')))

#assign relative position
unique_uorf <- unlist(uorf_pos_df_exonic %>% distinct(id))

uorf_poslist <- list()

for (i in 1:length(unique_uorf)) {
  u <- uorf_pos_df_exonic %>% filter(id == unique_uorf[i])
  if (u$strand[1] == '+') {
    u$uorf_pos <- 1:dim(u)[1]
  } else {
    u$uorf_pos <- dim(u)[1]:1
  }
  uorf_poslist[[i]] <- u
}

uorf_pos_df_exonic_wpos <- do.call("bind_rows", uorf_poslist)

write_tsv(uorf_pos_df_exonic_wpos, 'uorf_canonical_ATGonly_pos_df.tsv')

# uorf_pos_df_exonic_wpos <- read_tsv('uorf_canonical_ATGonly_pos_df.tsv', col_types = cols(chr = 'c'))

# make selection of UTR5 uorf locations

UTR5_uorf_pos <- filter(posdf, Region == 'UTR5') %>% left_join(uorf_pos_df_exonic_wpos, by = c('chr', 'pos')) %>% 
  mutate(uORF = ifelse(is.na(uorf_pos), FALSE, TRUE))

UTR5_uorf_pos_count <- UTR5_uorf_pos[UTR5_uorf_pos$uORF,] %>% count(pos_relative, uORF)

# ggplot(UTR5_uorf_pos_count, aes(pos_relative, n)) + geom_bar(stat = 'identity') +
#   coord_cartesian(xlim = c(-1000, 0)) +
#   theme_bw() + ggtitle('Count of uORFs per each UTR position relative to main ATG')
# 
# ggsave('uORF_positions_in_UTRs.png')

# Plot for figure 4

cumulative <- ggplot(UTR5_uorf_pos[UTR5_uorf_pos$uORF,], aes(x=abs(pos_relative))) + stat_ecdf() +
  coord_cartesian(xlim = c(0, 1500)) +
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
        legend.spacing.x = unit(0.2, 'line')) +
  ylab('Fraction') + xlab('Distance from start codon') +
  scale_x_continuous(breaks = seq(0, 1500, by = 100))

ggsave('uORF_positions_in_UTR_cumulative.pdf', plot = cumulative, width = unit(6, 'in'), height  = unit(4, 'in'))




# UTR5cds_uorf_pos <- filter(posdf, Region %in% c('UTR5', 'CDS')) %>% left_join(uorf_pos_df_exonic_wpos, by = c('chr', 'pos')) %>%
#   mutate(uORF = ifelse(is.na(uorf_pos), FALSE, TRUE))
# 
# UTR5cds_uorf_pos_count <- UTR5cds_uorf_pos[UTR5cds_uorf_pos$uORF,] %>% count(pos_relative, uORF)
# 
# uORF_positions_in_UTR_CDS <- ggplot(UTR5cds_uorf_pos_count, aes(pos_relative, n)) + geom_bar(stat = 'identity') +
#   coord_cartesian(xlim = c(-1000, 500)) +
#   theme_bw() + ggtitle('Count of uORFs per each UTR position relative to main ATG')
# 
# 
# ggsave('uORF_positions_in_UTR_CDS.pdf', plot = uORF_positions_in_UTR_CDS, width = unit(6, 'in'), height  = unit(4, 'in'))
# 



