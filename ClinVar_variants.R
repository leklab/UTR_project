# A script to produce panels B and C for Figure 1.

# Load libraries
library(tidyverse)

######################################################################################
### Clinvar Data

# Download clinvar vcf file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
# In this script the relese clinvar_20180930.vcf.gz is used

# First I prepared the vcf file with GATK:

# java -jar /Users/sander/data_tools/GenomeAnalysisTK-3.6-0/GenomeAnalysisTK.jar \
#  -R /Users/sander/data_references/Homo_sapiens_assembly19.fasta \
#  -T VariantsToTable \
#  -V clinvar_20180930.vcf.gz \
#  --allowMissingData \
#  -F CHROM -F POS -F REF -F ALT -F ID -F QUAL -F TYPE -F CLNSIG -F CLNSIGINCL -F CLNVI -F CLNVC -F MC -F GENEINFO \
#  -o clinvar_20180930.table


# Read in the clinvar

clinvar <- read_tsv('UTR_gnomad/clinvar_20180930.table', col_types = cols('CHROM'  = 'c', 'CLNSIGINCL' = 'c', 'CLNVI' = 'c')) 

######################################################################################
### gnomAD data

# Intervals for extracting gnomAD variants were made by script Reions.R , 
# the commands for extracting variants are listed in region.select.sh

# Read in gnomAD

# Files from UTR and cds
files <- c('gnomad.genomes.UTR5.table', 'gnomad.genomes.UTR3.table', 'gnomad.genomes.CDS.table')
for (file in files) {
  prefix <- str_remove(file, 'gnomad.genomes.')
  prefix <- str_remove(prefix, '.table')
  
  df <- read_tsv(file, col_types = 'cicccdiid') %>%
    mutate(Region = prefix)
  
  assign(prefix, df)
}

# Combine utr and cds datasets, and add frequency categories

utr_cds <- bind_rows(UTR5,CDS,UTR3) %>% 
  mutate(TYPE = ifelse((nchar(REF) == 1 & nchar(ALT) == 1), 'snps', 'indels'),
  FREQ = case_when(
  AC == 1 ~ 'AC = 1',
  AC %in% 2:10 ~ 'AC = 2..10',
  AC > 10 & AF < 0.001 ~ 'AC > 10, AF < 0.1%',
  AF < 0.01 ~ 'AF = 0.1%..1%',
  AF < 0.05 ~ 'AF = 1%..5%',
  TRUE ~ 'AF > 5%'
)) 

utr_cds$UTR <- factor(utr_cds$Region, levels = c('UTR5', 'CDS', 'UTR3'))
utr_cds$FREQ <- factor(utr_cds$FREQ, levels = c('AF > 5%', 'AF = 1%..5%', 'AF = 0.1%..1%', 'AC > 10, AF < 0.1%', 'AC = 2..10', 'AC = 1'))

#clean up

rm(UTR5,CDS,UTR3,files,df)

#######################################################################################################################

### Make the plots for Figure 1

# Join gnomad and clinvar variants

utr_cds_inclinvar <- utr_cds %>% inner_join(clinvar, by = c('CHROM', 'POS', 'REF', 'ALT'))

gnomad_pct_in_clinvar <- tibble(Region = factor(c('UTR5', 'CDS', 'UTR3'), levels =  c('UTR5', 'CDS', 'UTR3')), 
                                GnomAD = as.numeric(table(utr_cds$UTR)), 
                                ClinVar = as.numeric(table(utr_cds_inclinvar$UTR))) %>% 
                        mutate(Pct_in_ClinVar = ClinVar / GnomAD)

# Plot for fig 1A

gnomad_clinvar <- ggplot(gnomad_pct_in_clinvar, aes(x=Region, y=Pct_in_ClinVar)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text=element_text(size=15), 
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        legend.key.size = unit(2, 'line' ),
        legend.position = c(0.75,0.9),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line' )) +
  ylab('Fraction')

ggsave('gnomad_pct_in_clinvar.pdf', gnomad_clinvar)

########
# Classify each ClinVar variant according to the region affected with hierarchy in case of multiple transcripts, exonic > 5'UTR > 3'UTR > intronic

exonic <- c('splice_acceptor_variant', 'splice_donor_variant', 'missense_variant', 'nonsense', 'frameshift_variant', 'synonymous_variant') 

clinvar2 <- clinvar %>% mutate(Region = case_when(
  str_detect(MC, paste(exonic, collapse="|")) ~ 'Coding',
  str_detect(MC, '5_prime_UTR_variant') ~ 'UTR5',
  str_detect(MC, '3_prime_UTR_variant') ~ 'UTR3',
  str_detect(MC, 'intron_variant') ~ 'Intron',
  TRUE ~ 'Upstream_downstream'
), Class = case_when(
  CLNSIG %in% c('Likely_pathogenic', 'Pathogenic', 'Pathogenic/Likely_pathogenic') ~ '(Likely) pathogenic',
  CLNSIG %in% c('Likely_benign', 'Benign', 'Benign/Likely_benign') ~ '(Likely) benign',
  TRUE ~ 'VUS'
)) %>% filter(Region != 'Upstream_downstream')


# Create a summary counts table for annotation
clinvar_counts <- clinvar2 %>% count(Region) %>% mutate(label = paste('n=', n, sep = ''))

# Plot for Fig 1B
sumplot <- ggplot(clinvar2, aes(x ="", fill = Region)) +
  geom_bar(position = "fill") + 
  coord_polar("y") +
  theme_void() +
  theme(axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid  = element_blank(), 
        axis.title = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        legend.key.size = unit(2, 'line' ),
        legend.position = c(0.5,0.05),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line' ))

ggsave('sumplot.pdf', sumplot)


# Plot for Fig 1C
p <- ggplot(clinvar2, aes(x=1, fill = Class)) +
  geom_bar(position = "fill") + 
  facet_wrap(Region ~ .) + 
  coord_polar("y", direction = -1) +
  scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF")) +
  geom_text(data=clinvar_counts, mapping = aes(x = 1.6, y = 1, label = label), size = 5, inherit.aes = F) +
  theme_minimal() + 
  theme(axis.line = element_blank(),
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid  = element_blank(), 
        axis.title = element_blank(),
        strip.text = element_text(size=15),
        strip.placement = 'inside',
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        legend.key.size = unit(2, 'line' ),
        legend.position = c(0.5,0.01),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line' ))

ggsave('Fig1B.pdf', plot = p)

