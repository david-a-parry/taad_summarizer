#!/usr/bin/env perl 
#David A. Parry, October 2015

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Term::ProgressBar;
use POSIX qw/strftime/;
use HTTP::Tiny;
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
    't|transcripts=s',
    'u|uniprot_info=s',
    'q|quiet',
    'h|?|help',
) or usage("Syntax error");

usage() if $opts{h};
usage("Error: a gene list or gene ID must be provided") if not $opts{i} and not $opts{l};
#usage("Error: both --transcripts and --uniprot_info output filenames must be provided") if not $opts{t} or not $opts{u};
$opts{s} = "human" if not $opts{s};
my $parser = new IdParser();
my $restQuery = new EnsemblRestQuery();
my $http = HTTP::Tiny -> new(); 

my %id_mapping = (); #keep track of which genes came from which ID input
my %genes = (); #key is gene ID, value is hashses of gene info retrieved from Ensembl
my %transcript_ranks = (); #key is gene, value is hash of transcript IDs to scores (higher being more preferential)
my %transcript_xref = (); #key is transcript ID, values are hashes of databases and respective IDs
my %uniprot_info = (); 
my %enst_to_uniprot = ();

if ($opts{l}){
    push @gene_ids, readList();
}

die "No gene IDs provided - nothing to do!\n" if not @gene_ids;
=cut
open (my $TRANSC, ">$opts{t}") or die "Could not open $opts{t} for writing: $!\n";
open (my $UNIPRO, ">$opts{u}") or die "Could not open $opts{u} for writing: $!\n";
=cut
getGenesFromIds();
if (not %genes){
    die "Nothing to do!\n";
}
rankTranscriptsAndCrossRef();

=cut
print Dumper %transcript_ranks;
print Dumper %transcript_xref;
print Dumper %uniprot_info;
=cut 

foreach my $ensg (keys %genes){
    my $symbol = $genes{$ensg}->{display_name};
    foreach my $t (sort keys 
            {
                $transcript_ranks{$ensg}->{$a} <=> 
                $transcript_ranks{$ensg}->{$b} 
            } 
    %{$transcript_ranks{$ensg}}){
        
        
        
    }
}


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
        }
        foreach my $t (@longest_trans){
            $transcript_ranks{$ensg}->{$t}++;
        }
        foreach my $t (@{$genes{$ensg}->{Transcript}}){
            $transcript_ranks{$ensg}->{$t}++ if exists $enst_to_uniprot{$t};
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
                if ($ext->{dbname} eq 'Uniprot/SWISSPROT'){
                    if (not exists $uniprot_info{$ext->{primary_id}}){
                        getUniprotData($ext->{primary_id});
                    }
                    if (exists $enst_to_uniprot{$id}){
                        push @{$transcript_xref{$id}->{$ext->{dbname}}} , $ext->{primary_id};
                    }
                }else{
                    push @{$transcript_xref{$id}->{$ext->{dbname}}} , $ext->{primary_id};
                }
            }
        }
    }
}

#########################################################
sub getUniprotData{
    my $id = shift;
    my $site = "http://www.uniprot.org/uniprot";
    my $url = "$site/$id.txt";
    my $txt = getHttpData($url);
    die "Failed to retrieve Uniprot info for $id.\n"
        ."Tried URL: $url\nExiting\n" unless $txt;
    $url = "$site/$id.gff";
    my $gff = getHttpData($url);
    die "Failed to retrieve Uniprot info for $id.\n"
        ."Tried URL: $url\nExiting\n" unless $gff;
    # below collects only the transcript that code for
    # the canonical uniprot isoform.
    my @t = parseUniprotFlat($txt);
    if (not @t){
        die "ERROR: Could not identify Ensembl transcript for canonical isoform of $id!\n";
    }
    foreach my $t (@t){
        $enst_to_uniprot{$t} = $id;
    }
    $uniprot_info{$id} = parseUniprotGff($gff);
}

#########################################################
sub parseUniprotFlat{
    my $txt = shift;
    my @lines = split("\n", $txt);
    my $canonical; 
    my @transcripts = ();
    foreach my $l (@lines){
        if ($l =~ /^CC\s+IsoId=(\S+);\s+Sequence=Displayed;/){
            $canonical = $1;
            next;
        }
        if ($canonical){
            if($l =~ /DR\s+Ensembl;\s+(ENST\d+);\s+ENSP\d+;\s+ENSG\d+\. \[$canonical\]/){
                push @transcripts, $1;
            }
        }
    }
    return @transcripts;
}
#########################################################
sub parseUniprotGff{
    my $gff = shift;
    my @lines = split("\n", $gff);
    my %features = ();
    foreach my $l (@lines){
        next if $l =~ /^#/;
        chomp $l;
        my @g = split("\t", $l);
        my ($f, $start, $end, $details) = @g[2..4,8];
        next if $f eq 'Chain';
        if ($details =~ /Note=(.*)[\n\r;]/){
            $f .= "|$1";
        }
        push @{$features{"$start-$end"}}, $f;
    }
    return \%features;
}   


#########################################################
sub getHttpData{
    my $url = shift;
    my $response = $http->get($url);
    return if not $response->{success};
    return $response->{content};
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
    
    usage: $0 -l gene_list.txt -t transcript_output.txt -u uniprot_output.txt 
           $0 -i gene_symbol/gene_id/transcript_id -t transcript_output.txt -u uniprot_output.txt

    Options:

    -l, --list FILE
        Input file of gene/transcript/protein IDs one per line. Any whitespace following the first word of each line will be ignored.
    -i, --ids
        One or more gene/transcript/protein identifiers to look up. May be used instead or as well as --list file.
    -s, --species
        Species to search (only applies to non-Ensembl identifiers). Default = human.
    -t, --transcripts FILE
        Output file for transcripts, their ranks and cross refs. Required.
    -u, --uniprot-info FILE
        Output file for uniprot information for corresponding proteins. Required.
    -?, -h, --help
        Show this help message.


    Information:
    
    Note that even when using transcript or protein IDs, the parent gene will be identified and all transcripts processed.

EOT
;
    exit 1 if $msg;
    exit;
}
