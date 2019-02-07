#Script for producing Figure 2 in grant application

library(tidyverse)

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

utr_cds$Region <- factor(utr_cds$Region, levels = c('UTR5', 'CDS', 'UTR3'))
utr_cds$FREQ <- factor(utr_cds$FREQ, levels = c('AF > 5%', 'AF = 1%..5%', 'AF = 0.1%..1%', 'AC > 10, AF < 0.1%', 'AC = 2..10', 'AC = 1'))

#clean up

rm(UTR5,CDS,UTR3,files,df)


# Generate allele frequency plot

AFplot <- ggplot(utr_cds, aes(Region, fill = FREQ)) +
  geom_bar(position = 'fill') +
  theme_classic() +
  theme(axis.text=element_text(size=15), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.key.size = unit(2, 'line' ),
        #legend.position = c(0.75,0.9),
        #legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line')) +
  scale_fill_brewer(aesthetics = c('color', 'fill'),direction = -1) +
  ylab('Fraction')

ggsave('Fig2A.pdf', plot = AFplot)

# Generate region variant counts plot

#read_in position dataframe
posdf <- read_tsv('position_dataframe_UTRs_CDS.tsv', col_types = cols(chr = 'c'))
uorf_pos_df<- read_tsv('uorf_canonical_ATGonly_pos_df.tsv', col_types = cols(chr = 'c'))

#eliminate duplicate variants
utr_cds_unique <- utr_cds %>% distinct(CHROM, POS, REF, ALT, .keep_all = T)

#join variants with regions and uorfs
regs <- utr_cds %>% left_join(posdf, by = c( 'CHROM' = 'chr', 'POS' = 'pos', 'Region'))  %>% 
  left_join(uorf_pos_df,by = c( 'CHROM' = 'chr', 'POS' = 'pos'), suffix = c('.var', '.uorf')) 

# add specific regions
regs <- regs %>% mutate(kozak = ifelse(Region %in% c('UTR5', 'CDS') & pos_relative %in% -6:4, 'in_kozak', 'not_in_kozak'),
                        start_codon = ifelse(Region == 'CDS' & pos_relative %in% 1:3, 'start_codon', 'not_start_codon'),
                        kozak_but_not_start = ifelse(Region %in% c('UTR5', 'CDS') & pos_relative %in% c(-6:-1,4), 'in_kozak_but_not_start', 'not_in_kozak_or_start'),
                        uORF = ifelse(!is.na(uorf_pos), 'in_uORF', 'not_in_uORF'))

# eliminate duplicates

regs_unique <- regs %>% distinct(CHROM, POS, REF, ALT, .keep_all = T)

#write_tsv(regs_unique, 'regs.tsv')

categories <- regs_unique  %>% filter(Region %in% c('UTR5', 'CDS')) %>% filter(TYPE != 'indels') %>%
  select(CHROM, POS, REF, ALT, Region, pos_relative, uorf_pos, FREQ, uORF, kozak,  start_codon, kozak_but_not_start)

categories_long <- categories %>% gather(key, value, -CHROM, -POS, -REF, -ALT, -Region, -pos_relative, -uorf_pos, -FREQ)

cat_long_in <- categories_long %>% filter(value %in% c('in_kozak', 'in_kozak_but_not_start', 'start_codon', 'in_uORF'))

# make dataframe showing total width of each region and number of variants in each region to calculate variant density.
region <- c('uORF',  'mATG', 'Kozak', 'UTR5', 'CDS', 'UTR3')

width <- c(dim(semi_join(uorf_pos_df, filter(posdf, Region %in% c('UTR5', 'CDS'))))[1], 
           dim(filter(posdf, Region == 'CDS' & pos_relative %in% 1:3))[1], 
           dim(filter(posdf, Region %in% c('UTR5', 'CDS') & pos_relative %in% c(-6:-1,4)))[1],
           dim(filter(posdf, Region == 'UTR5'))[1],
           dim(filter(posdf, Region == 'CDS'))[1],
           dim(filter(posdf, Region == 'UTR3'))[1])

count <- c(table(cat_long_in$value)['in_uORF'],
           table(cat_long_in$value)['start_codon'],
           table(cat_long_in$value)['in_kozak_but_not_start'],
           table(utr_cds_unique$Region)['UTR5'],
           table(utr_cds_unique$Region)['CDS'],
           table(utr_cds_unique$Region)['UTR3'])

var_per_kb <- count / (width / 1000)

cat_df <- tibble(region, width, count, var_per_kb)

cat_df$region <- factor(cat_df$region, levels = c('UTR5', 'uORF', 'Kozak', 'mATG', 'CDS', 'UTR3'))

# Plot figure 2B
fig2b <- ggplot(cat_df, aes(x=region, y=var_per_kb)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = count, y = var_per_kb + 5 ), size = 6) +
  theme_classic() +
  theme(axis.text=element_text(size=15), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        legend.key.size = unit(2, 'line' ),
        legend.position = c(0.75,0.9),
        legend.direction = 'horizontal',
        legend.spacing.x = unit(0.2, 'line' )) +
  ylab('Variants per kilobase')

ggsave('variant_densities.pdf', plot = fig2b)

