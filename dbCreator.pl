#!/usr/bin/env perl 
#David A. Parry, October 2015

use strict;
use warnings;
use POSIX qw/strftime/;
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
    'u|uniprot-info=s',
    'c|cdd-features=s',
    'q|quiet',
    'h|?|help',
) or usage("Syntax error");

usage() if $opts{h};
usage("Error: a gene list or gene ID must be provided") if not $opts{i} and not $opts{l};
$opts{s} = "human" if not $opts{s};
my $parser = new IdParser();
my $restQuery = new EnsemblRestQuery();
my $http = HTTP::Tiny -> new(); 

my %id_mapping = (); #keep track of which genes came from which ID input
my %genes = (); #key is gene ID, value is hashses of gene info retrieved from Ensembl
my %transcript_ranks = (); #key is gene, value is hash of transcript IDs to scores (higher being more preferential)
my %transcript_xref = (); #key is transcript ID, values are hashes of databases and respective IDs
my %uniprot_info = (); 
my %uniprot_to_genename = ();
my %enst_to_uniprot = ();

if ($opts{l}){
    push @gene_ids, readList();
}

die "No gene IDs provided - nothing to do!\n" if not @gene_ids;
my $TRANSC = \*STDOUT; 
if ($opts{t}){
    open ($TRANSC, ">$opts{t}") or die "Could not open $opts{t} for writing: $!\n";
}
open (my $UNIPRO, ">$opts{u}") or die "Could not open $opts{u} for writing: $!\n" if $opts{u};
open (my $CDD,    ">$opts{c}") or die "Could not open $opts{c} for writing: $!\n" if $opts{c};

getGenesFromIds();
if (not %genes){
    die "Nothing to do!\n";
}
rankTranscriptsAndCrossRef();

printTranscriptInfo();
close $TRANSC if $opts{t} and $TRANSC;

printUniprotInfo() if $opts{u};
close $UNIPRO if $UNIPRO;

retrieveAndPrintCddInfo() if $opts{c};
close $CDD if $CDD;

#########################################################
sub parseCddResult{
    my $data = shift;
    my @lines = split("\n", $data);
    foreach my $l (@lines){
        next if $l =~ /^#/;
        chomp $l;
        next if not $l;
        next if $l =~ /^Query/;
        my @s = split("\t", $l); 
        if ($s[0] =~ /Q\#\d+ - (\S+)/){
            my $u = $1;
            my $name = $uniprot_to_genename{$u};
            print $CDD join("\t", 
              (
                $name, 
                $u,
                $s[2],
                $s[3],
              ) 
            ) . "\n";
        }else{
            die "Could not parse CDD result:\n$l\n";
        }
    }
}

#########################################################
sub retrieveAndPrintCddInfo{
    my $result = retrieveCddFeatures([keys %uniprot_info]);
    #print $result;#DEBUG
    my $cdd_data = parseCddResult($result);
    print $CDD $cdd_data;
}


#########################################################
sub retrieveCddFeatures{
    my $queries = shift; 
    my $url = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";
    my %opts = (
      evalue   => 0.001,
      tdata    => 'feats',
      queries  => $queries
                    
    );
    my $rid;
    my $response = $http->post_form(
        $url,
        \%opts,
    );
    die "Error: ", $response->{status}
    unless $response->{success};

    if($response->{content} =~ /^#cdsid\s+([a-zA-Z0-9-]+)/m) {
        $rid =$1;
        informUser( "Search with Request-ID $rid started.\n");
    }else{
        die "Submitting the search failed,\n can't make sense of response: $response->content\n";
    }
 
    # checking for completion, wait 2 seconds between checks

    my $done = 0;
    my $status = -1;
    while ($done == 0) {
        sleep(2);
        $response = $http->post_form(
            $url,
            [
              'tdata' => "feats",
              'cdsid' => $rid
            ],
        );
        die "Error: ", $response->{status}
          unless $response->{success};

        if ($response->{content} =~ /^#status\s+([\d])/m) {
            $status = $1;
            if($status == 0) {
                $done = 1;
                informUser("CDD search has completed, retrieving results ..\n");
            }elsif($status == 3) {
                my $time = strftime( "%H:%M:%S", localtime );
                informUser("[$time] CDD search still running...\n");
            }elsif($status == 1) {
                die "Invalid request ID\n";
            }elsif($status == 2) {
                die "Invalid input - missing query information or search ID\n";
            }elsif($status == 4) {
                die "Queue Manager Service error\n";
            }elsif($status == 5) {
                die "Data corrupted or no longer available\n";
            }
        }else{
            die "Checking search status failed,\ncan't make sense of response: $response->{content}\n";
        }
    }
    # retrieve and display results
    $response = $http->post_form(
        $url,
        [
            'cdsid'  => $rid
        ],
    );
    die "Error: ", $response->{status}
      unless $response->{success};

    return $response->{content};
}


#########################################################
sub printUniprotInfo{
    informUser("Writing Uniprot info to $opts{u}...\n");
    print $UNIPRO '#' . join("\t", 
                qw / 
                    GeneName
                    UniprotId
                    Start
                    End
                    Feature
                    Note
                    /) . "\n";
    foreach my $id (sort keys %uniprot_info){
        my $name = $uniprot_to_genename{$id};
        foreach my $k ( sort by_unipro_coords keys %{$uniprot_info{$id}}){
            my ($start, $end) = split('-', $k);
            foreach my $f (@{$uniprot_info{$id}->{$k}}){
                my ($feature, $note) = split(/\|/, $f );
                $note ||= '.';
                print $UNIPRO join("\t", 
                    (
                        $name,
                        $id,
                        $start,
                        $end,
                        $feature,
                        $note,
                    ) 
                ) . "\n";
            }
        }
    }
}

#########################################################
sub by_unipro_coords{
    my @asplit = split('-', $a);
    my @bsplit = split('-', $b);
    if (my $diff = $asplit[0] - $bsplit[0]){
        return $diff;
    }
    return $asplit[1] - $bsplit[1];
}

#########################################################
sub printTranscriptInfo{
    informUser("Writing transcript info to $opts{t}...\n");
    print $TRANSC '#' . join("\t", 
                qw /
                    Symbol
                    EnsemblGeneID
                    EnsemblTranscriptID
                    EnsemblProteinID
                    RefSeq_mRNA
                    CCDS
                    Uniprot
                    TranscriptScore
                /) . "\n";

    foreach my $ensg (keys %genes){
        foreach my $t (sort 
        {
            $transcript_ranks{$ensg}->{$b} <=> 
            $transcript_ranks{$ensg}->{$a} ||
            exists $enst_to_uniprot{$b} <=> 
            exists $enst_to_uniprot{$a}  
        } 
        keys %{$transcript_ranks{$ensg}} ){
            my $symbol = $genes{$ensg}->{display_name};
            my $score = $transcript_ranks{$ensg}->{$t};
            my $ensp = ".";
            my $refseq = ".";
            my $ccds = ".";
            my $uniprot = ".";
            if (exists $transcript_xref{$t}->{Translation}){
                $ensp = $transcript_xref{$t}->{Translation};
            }
            if (exists $transcript_xref{$t}->{RefSeq_mRNA}){
                $refseq = join("/", @{ $transcript_xref{$t}->{RefSeq_mRNA} } );
            } 
            if (exists $transcript_xref{$t}->{CCDS}){
                $ccds = join("/", @{ $transcript_xref{$t}->{CCDS} } );
            } 
            if (exists $enst_to_uniprot{$t}){
                $uniprot =  $enst_to_uniprot{$t};
            }
            
            print $TRANSC join("\t", 
                (
                 $symbol, 
                 $ensg,
                 $t, 
                 $ensp,
                 $refseq,
                 $ccds,
                 $uniprot,
                 $score,
                ) 
            ) . "\n";
                 
        }
    }
}
#########################################################
sub rankTranscriptsAndCrossRef{
     foreach my $ensg (keys %genes){
        if (not $genes{$ensg}->{Transcript}){
            print STDERR "WARNING: No transcripts identified for \"$ensg\"/\"" 
                . $genes{$ensg}->{display_name} . "\"\n";
            next;
        }
        my @longest_trans = ();#
        my $longest = 0; 
        foreach my $t (@{$genes{$ensg}->{Transcript}}){
            $transcript_ranks{$ensg}->{$t->{id}} = 0;
            $transcript_ranks{$ensg}->{$t->{id}}++ if $t->{is_canonical};
            if ($t->{Translation}){
                $transcript_xref{$t->{id}}->{Translation} = $t->{Translation}->{id};
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
        foreach my $id (@longest_trans){
            $transcript_ranks{$ensg}->{$id}++;
        }
        foreach my $t (@{$genes{$ensg}->{Transcript}}){
            $transcript_ranks{$ensg}->{$t->{id}} += 2 if exists $enst_to_uniprot{$t->{id}};
            $transcript_ranks{$ensg}->{$t->{id}}++ if exists $transcript_xref{$t->{id}}->{RefSeq_mRNA};
            $transcript_ranks{$ensg}->{$t->{id}}++ if exists $transcript_xref{$t->{id}}->{CCDS};
        }
    }
}

#########################################################
sub addXrefs{
    my $id = shift;
    informUser( "Retrieving external IDs for transcript $id...\n");
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
    informUser( "Retrieving Uniprot text data for $id...\n");
    my $site = "http://www.uniprot.org/uniprot";
    my $url = "$site/$id.txt";
    my $txt = getHttpData($url);
    die "Failed to retrieve Uniprot info for $id.\n"
        ."Tried URL: $url\nExiting\n" unless $txt;
    $url = "$site/$id.gff";
    informUser( "Retrieving Uniprot GFF data for $id...\n");
    my $gff = getHttpData($url);
    die "Failed to retrieve Uniprot info for $id.\n"
        ."Tried URL: $url\nExiting\n" unless $gff;
    # below collects only the transcript that code for
    # the canonical uniprot isoform.
    my ($name, @t)  = parseUniprotFlat($txt);
    if (not @t){
        die "ERROR: Could not identify Ensembl transcript for canonical isoform of $id!\n";
    }
    foreach my $t (@t){
        $enst_to_uniprot{$t} = $id;
    }
    $uniprot_info{$id} = parseUniprotGff($gff);
    $uniprot_to_genename{$id} = $name;
}

#########################################################
sub parseUniprotFlat{
    my $txt = shift;
    my @lines = split("\n", $txt);
    my $canonical; 
    my $name = '';
    my @transcripts = ();
    foreach my $l (@lines){
        if ($l =~ /^GN\s+Name=(\S+);/){
            $name = $1;
            next;
        }
        if ($l =~ /^CC\s+IsoId=(\S+);\s+Sequence=Displayed;/){
            $canonical = $1;
            next;
        }
        if ($canonical){
            if($l =~ /^DR\s+Ensembl;\s+(ENST\d+);\s+ENSP\d+;\s+ENSG\d+\. \[$canonical\]/){
                push @transcripts, $1;
            }
        }elsif($l =~ /^DR\s/){#if only one isoform there won't be IsoId's
            if($l =~ /^DR\s+Ensembl;\s+(ENST\d+);\s+ENSP\d+;\s+ENSG\d+/){
                push @transcripts, $1;
            }
        }
    }
    return ($name, @transcripts);
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
        next if $f eq 'Alternative sequence';
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
            informUser( "Interpretting ID \"$g\" as of type \"" . $parser->get_identifierType() . "\"...\n");
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
            informUser("Identifying Ensembl gene via transcript cross-reference...\n");
            my $transcript = $restQuery->getTranscriptViaXreg($g, $opts{s});
            if ($transcript){
                if (exists $transcript->{id}){
                    $gene_hash = geneFromEnst($transcript->{id});
                }
            }
        }else{
            informUser("Identifying Ensembl gene via gene cross-reference...\n");
            my $gene = $restQuery->getGeneViaXreg($g, $opts{s});
            if (ref $gene eq 'ARRAY'){
                if ($gene->[0]->{id}){
                    $gene_hash = $restQuery->lookUpEnsId($gene->[0]->{id}, 1);
                }
            }
        }
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
    informUser("Identifying parent gene from Ensembl transcript...\n");
    my $id = shift;
    return $restQuery->getParent($id, 1);
}
#########################################################
sub geneFromEnsp{
    informUser("Identifying parent gene from Ensembl protein...\n");
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
sub informUser{
    return if $opts{q};
    my $msg = shift;
    print STDERR $msg;
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
        Output file for transcripts, their ranks and cross refs. Will print to STDOUT if not provided.
    -u, --uniprot-info FILE
        Output file for uniprot information for corresponding proteins. Optional.
    -c, --cdd-features FILE
        Output file for Conserved Domain Database feature residue information for corresponding proteins. Optional.
    -?, -h, --help
        Show this help message.


    Information:
    
    Note that even when using transcript or protein IDs, the parent gene will be identified and all transcripts processed.

EOT
;
    exit 1 if $msg;
    exit;
}
