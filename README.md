# UTR project
Scripts used for Lek Labs UTR project.

## General workflow

_This analysis relies on canonical intervals, as defined in Ensembl ensembl.org ; gnomAD genomes dataset https://gnomad.broadinstitute.org/; gencode annotation v19 https://www.gencodegenes.org/, uORFs are extracted from from Gerstein lab paper https://doi.org/10.1093/nar/gky188. See comments in scripts for details_

1) Import gtf file, and separate UTR3s and UTR5s, then make GATK-compatible intervals files for each region of interest (UTR5, CDS, UTR3) `Regions.R`
2) Download gnomad vcf file for genomes dataset, and select variants based on the regions defined in #1. `region_select.sh`
3) Make position dataframe: each position is marked with region and relative position in that region (UTR5 - position from ATG, negative integers; CDS - position from ATG, positive integers; UTR3 position from stop codon, positive integers). `Position_dfs.R`
4) For ClinVar variants analysis run (Figure 1) `ClinVar_variants.R`, dependent on #2
5) For uORF analysis run `uORFs.R`, dependent on #3.
6) For Kozak sequence logio plot run `Kozak.R`, dependent on #3
7) For ExAC coverage run `ExAC_coverage.R`.
8) For analysis on allele frequencies and counts of variants in different regions run `Figure2.R`, dependent on #2 and #3.
