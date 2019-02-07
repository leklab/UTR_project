# downloaded the gnomad genomes vcf file
# run the following command for each region from Region.R script. No need to separate for different strands at this point.
# region=UTR5, region=UTR3, region=CDS
# GATK version used The Genome Analysis Toolkit (GATK) v3.7-0-gcfedb67, Compiled 2016/12/12 11:21:18

java -Xmx34G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
        -T SelectVariants  \
        -R ~/project/resources/Homo_sapiens_assembly19.fasta \
        -L ${region}.intervals \
        -V gnomad.genomes.r2.1.sites.vcf.bgz \
        -o gnomad.genomes.${region}.vcf.gz
        
java -Xmx20G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
     -R ~/project/resources/Homo_sapiens_assembly19.fasta \
     -T VariantsToTable \
     -V gnomad.genomes.${region}.vcf.gz \
     -F CHROM -F POS -F REF -F ALT -F ID -F QUAL -F AC -F AN -F AF \
     -o gnomad.genomes.${region}.table
     