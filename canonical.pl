use strict;
use warnings;
 
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";
 
$reg->load_registry_from_db(
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
   -dbname  => 'homo_sapiens_core_93_37'
);
 
open(OUTFILE, ">human_canonical_transcripts_v93.fa");
 
my $gene_adaptor = $reg->get_adaptor('human', 'core', 'gene');
my @genes = @{$gene_adaptor->fetch_all};
 
while(my $gene = shift @genes) {
  print
    (OUTFILE
        $gene->canonical_transcript->stable_id, ".", $gene->canonical_transcript->version, "\n");
}
