#!/usr/bin/env perl 
use strict;
use warnings;
use Data::Dumper;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use FindBin qw($RealBin);
use lib "$RealBin/lib/vcfhacks/lib/";
use VcfReader;

if (@ARGV != 1){
    die <<EOT

Converts ClinVar TSV file (as available from https://github.com/macarthur-lab/clinvar) into a VCF file with the ClinVar information placed into the INFO field. 

Usage: $0 clinvar.tsv > clinvar.vcf

EOT
;
}
 
my $in = shift;
#OPEN AND CHECK ClinVar INPUT
my $FH;
if ($in =~ /\.gz$/){
    $FH = new IO::Uncompress::Gunzip $in or die "IO::Uncompress::Gunzip failed while opening $in for reading:\n$GunzipError";
}else{
    open ($FH, $in ) or die "Can't open $in for reading: $!\n";
}
my %clinvar_columns = ();
my $header = <$FH>;
chomp (my @head = split("\t", $header)); 
my $n = 0;
%clinvar_columns = map { lc($_) => $n++ } map { (my $trim = $_) =~ s/^#+//; $trim } @head;
my @clinvar_fields = qw(
    mut
    measureset_id
    symbol
    clinical_significance
    review_status
    hgvs_c
    hgvs_p
    all_submitters
    all_traits
    all_pmids
    pathogenic
    conflicted
);
foreach my $req (@clinvar_fields){
    if (not exists $clinvar_columns{$req}){
        die "Required column ($req) not found in HGMD file $in\nFound the following columns:\n" 
            .join("\n", @head) . "\n";
    }
}

writeHeader();

while (my $line = <$FH>){
    chomp $line;
    my @fields = split("\t", $line); 
    my $vcf_record = convertToVcf(\@fields);
    print "$vcf_record\n";
}
close $FH;


###########################################################
sub writeHeader{
    print "##fileformat=VCFv4.2\n";
    foreach my $f (@clinvar_fields){
        (my $tag = $f) =~ s/\s/_/g; 
        print "##INFO=<ID=$tag,Number=.,Type=String,Description=\"$f field from MaCarthur Lab's ClinVar table (https://github.com/macarthur-lab/clinvar), source: \"$in\">\n"
    }
    print '#' . join("\t",  ( qw / CHROM POS ID REF ALT QUAL FILTER INFO / ) ) . "\n";
}

###########################################################
sub convertToVcf{
    my $clinvar = shift;
    my @out = ();
    my @details = ();
    my $chrom = $clinvar->[ $clinvar_columns{"chrom"} ];
    my $pos = $clinvar->[ $clinvar_columns{"pos"} ];
    my $ref = $clinvar->[ $clinvar_columns{"ref"} ];
    my $alt = $clinvar->[ $clinvar_columns{"alt"} ];
    push @out, 
        $chrom,
        $pos,
        $clinvar->[ $clinvar_columns{"measureset_id"} ],
        $ref,
        $alt,
        ".",
        ".", #blank QUAL and FILTER fields
    ;
    foreach my $d (@clinvar_fields){ 
        (my $tag = $d) =~ s/\s/_/g;
        (my $value = $clinvar->[ $clinvar_columns{$d} ]) =~ s/\s/_/g;
        push @details, "$tag=$value";
    }
    
    push @out, join(";", @details);
    return join("\t", @out) ;
}
