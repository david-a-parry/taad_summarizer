#!/usr/bin/env perl 
#David A. Parry, October 2015

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Term::ProgressBar;
use POSIX qw/strftime/;
use FindBin;
use lib "$FindBin::Bin/lib";
use IdParser;
use EnsemblRestQuery;

my @gene_ids = ();
my %opts = ( i => \@gene_ids);
GetOptions(
    \%opts,
    'l|list=s',
    'i|id=s{,}',#  => \@gene_ids,
    's|species=s',
    'o|output=s',
    'q|quiet',
    'h|?|help',
) or usage("Syntax error");

usage() if $opts{h};
usage("Error: a gene list or gene ID must be provided") if not $opts{i} and not $opts{l};
$opts{s} = "human" if not $opts{s};
my $parser = new IdParser();
my $restQuery = new EnsemblRestQuery();
my %id_mapping = (); #keep track of which genes came from which ID input
my %genes = (); #key is gene ID, value is hashses of gene info retrieved from Ensembl
my %transcript_ranks = (); #key is gene, value is hash of transcript IDs to scores (higher being more preferential)
my %transcript_xref = (); #key is transcript ID, values are hashes of databases and respective IDs
if ($opts{l}){
    push @gene_ids, readList();
}

die "No gene IDs provided - nothing to do!\n" if not @gene_ids;

getGenesFromIds();
if (not %genes){
    die "Nothing to do!\n";
}
rankTranscriptsAndCrossRef();


print Dumper %transcript_ranks;

print Dumper %transcript_xref;

#########################################################
sub rankTranscriptsAndCrossRef{
     foreach my $ensg (keys %genes){
        if (not $genes{$ensg}->{Transcript}){
            print STDERR "WARNING: No transcripts identified for \"$ensg}}->{id}\"/\"" 
                . $genes{$ensg}->{display_name} . "\"\n";
            next;
        }
        my @longest_trans = ();#
        my $longest = 0; 
        foreach my $t (@{$genes{$ensg}->{Transcript}}){
            $transcript_ranks{$ensg}->{$t->{id}} = 0;
            $transcript_ranks{$ensg}->{$t->{id}}++ if $t->{is_canonical};
            if ($t->{Translation}){
                $transcript_ranks{$ensg}->{$t->{id}}++ ;
                if ($t->{Translation}->{length} > $longest){
                    @longest_trans = ($t->{id}); 
                    $longest = $t->{Translation}->{length};
                }elsif($t->{Translation}->{length} == $longest){
                    push @longest_trans, $t->{id}; 
                }
                
            }
            addXrefs($t->{id}); 

       foreach my $l (@longest_trans){
            $transcript_ranks{$ensg}->{$l}++;
        }
    }
}

#########################################################
sub addXrefs{
    my $id = shift;
    my $external = $restQuery->getXrefs(id => $id, all_levels => 1);
    if (ref $external eq 'ARRAY'){
        foreach my $ext (@$external){
            if (grep { $ext->{dbname} eq $_ }  
                qw ( 
                    RefSeq_mRNA 
                    CCDS 
                    Uniprot/SWISSPROT
                ) 
            ){
                #NOTE that Uniprot entries may include partial/isoform matches
                push @{$transcript_xref{$id}->{$ext->{dbname}}} , $ext->{primary_id};
            }
        }
    }
}

#########################################################
sub getGenesFromIds{
    foreach my $g (@gene_ids){
        $parser->parseId($g); 
        if (not $opts{q}){
            print STDERR "Interpretting ID \"$g\" as of type \"" . $parser->get_identifierType() . "\"...\n";
        }
        my $gene_hash; 
        if ($parser->get_isEnsemblId()){
            if ( $parser->get_isTranscript() ){
                $gene_hash = geneFromEnst($g);
            }elsif( $parser->get_isProtein() ) {
                $gene_hash = geneFromEnsp($g);
            }else{
                $gene_hash = $restQuery->lookUpEnsId($g, 1);
            }
        }elsif($parser->get_isTranscript()  or $parser->get_isProtein() ) {
            my $transcript = $restQuery->getTranscriptViaXreg($g, $opts{s});
            if ($transcript){
                if (exists $transcript->{id}){
                    $gene_hash = geneFromEnst($transcript->{id});
                }
            }
        }else{
            my $gene = $restQuery->getGeneViaXreg($g, $opts{s});
            if (ref $gene eq 'ARRAY'){
                if ($gene->[0]->{id}){
                    $gene_hash = $restQuery->lookUpEnsId($gene->[0]->{id}, 1);
                }
            }
        }
        print Dumper $gene_hash;#DEBUG
        if (not $gene_hash){
            print STDERR "WARNING: Could not identify gene for ID \"$g\"\n";
        }else{
            if (exists $id_mapping{$gene_hash->{id}}){
                print STDERR "WARNING: Duplicated gene identified: $gene_hash->{id} identified from $g and from ".
                    join(", ", @{$id_mapping{$gene_hash->{id}}}) . "\n";
            }
            push @{$id_mapping{$gene_hash->{id}}}, $g;
            $genes{$gene_hash->{id}} = $gene_hash;
        }
    }
}


#########################################################
sub geneFromEnst{
    my $id = shift;
    return $restQuery->getParent($id, 1);
}
#########################################################
sub geneFromEnsp{
    my $id = shift;
    my $par = $restQuery->getParent($id);
    if ($par){
        if (exists $par->{id}){
            return geneFromEnst($par->{id});
        }
    }
}
#########################################################
sub readList{
    open (my $LIST, $opts{l}); 
    my @ids = (); 
    while (my $line = <$LIST>){
        chomp $line; 
        my @s = split(/\s+/, $line); 
        push @ids, $s[0];
    }
    if (not @ids){
        print STDERR "WARNING: No ids found in file $opts{l}.\n";
    }
    return @ids;
}


#########################################################
sub usage{
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print <<EOT
    
    usage: $0 [-i gene_list.txt | -s gene_symbol | -g gene_id | -t transcript_id ] [options]

    Options:

    -l, --list FILE
        Input file of gene/transcript/protein IDs one per line. Any whitespace following the first word of each line will be ignored.
    -i, --ids
        One or more gene/transcript/protein identifiers to look up. May be used instead or as well as --list file.
    -s, --species
        Species to search (only applies to non-Ensembl identifiers). Default = human.
    -o, --ouptut FILE
        Output file. Optional. Default is STDOUT.
    -?, -h, --help
        Show this help message.


    Information:
    
    Note that even when using transcript or protein IDs, the parent gene will be identified and all transcripts processed.

EOT
;
    exit 1 if $msg;
    exit;
}
