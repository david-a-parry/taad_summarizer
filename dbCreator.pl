#!/usr/bin/env perl 
#David A. Parry, October 2015

use strict;
use warnings;
use DBI;
use POSIX qw/strftime/;
use Getopt::Long;
use Data::Dumper;
use Term::ProgressBar;
use POSIX qw/strftime/;
use HTTP::Tiny;
use URI::Encode qw (uri_decode);
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
    'd|db=s', #sqlite database output
    's|species=s',
    'e|et_folder=s',#optional evolutionary trace folder containing protein predictions
    'q|quiet',
    'h|?|help',
) or usage("Syntax error");

usage() if $opts{h};
usage("Error: a gene list or gene ID must be provided") if not $opts{i} and not $opts{l};
usage("Error: please specify a filename for your database using the --db option!\n") 
    if not $opts{d};
$opts{s} = "human" if not $opts{s};


#set up database
if (-e $opts{d}){
    informUser("INFO: --db file '$opts{d}' already exists - additional data will be appended to this database file\n");
    #should really amend this so that appending is possible, not a priority yet though
}
my $driver   = "SQLite";
my $max_commit = 5000;
my $dbh = DBI->connect("DBI:$driver:$opts{d}", {RaiseError => 1}) 
    or die "Could not create sqlite database '$opts{d}': " . DBI->errstr . "\n";
#$dbh->do("create database $opts{d}") or die "Could not create database '$opts{d}': " . $dbh->errstr . "\n";
#$dbh->do("use $opts{d}") or die "Could not execute 'use $opts{d}' in sqlite: " . $dbh->errstr . "\n"; 

#set up our query objects
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

getGenesFromIds();
if (not %genes){
    die "Nothing to do!\n";
}
rankTranscriptsAndCrossRef();

outputTranscriptInfo();

outputEvolutionaryTraceInfo();

outputUniprotInfo() ; 

retrieveAndOutputCddInfo() ; 

$dbh->disconnect(); 


#########################################################
sub createTable{
    my ($table, $fields, $field_to_prop) = @_;
    my $create_string = join(", ", map { "$_ $field_to_prop->{$_}" } @$fields);
    my $stmt = "CREATE TABLE IF NOT EXISTS $table ($create_string)";
    informUser("Creating/appending to '$table' table in $opts{d}.\n");
    $dbh->do($stmt) or die "Error creating '$table' table: " . $dbh->errstr . "\n";
}

#########################################################
sub getInsertionQuery{
    my $table = shift;
    my $fields = shift;
    my $field_list = join(", ", @$fields);
    my $placeholders = join ", ", map {'?'} @$fields;
    return "INSERT OR REPLACE INTO $table ($field_list) VALUES ($placeholders)";
}


#########################################################
sub getSelectQuery{
    my $table = shift;
    my $fields = shift;
    my $field_list = join(", ", @$fields);
    my $queries = join " and ", map {"$_ == ?"} @$fields;
    #my $queries = join " and ", map {"($_ == ? or ($_ is null and ? == 1))"} @$fields;
    return "SELECT $field_list from $table where $queries";
}

#########################################################
sub addRow{
    my ($insth, $vals, $selth, ) = @_;
    my @multi_vals = map {($_, defined($_) ? 0 : 1) } @$vals;
    if ($selth){
        $selth->execute(@$vals) or die "Could not execute duplicate check query: " 
       # $selth->execute(@multi_vals) or die "Could not execute duplicate check query: " 
          . $selth->errstr;
        if ($selth->fetchrow_arrayref){
            informUser("Ignoring duplicate entry " . join("|", @$vals) . "\n");
            return 0;
        }
    }
    my $n = $insth->execute(@$vals) or die "Could not insert values: " . $insth->errstr;
    return $n;
}

#########################################################
sub parseCddFeats{
    my $data = shift;
    my $f = shift; 
    my @lines = split("\n", $data);
    my $insert_query = getInsertionQuery("cdd", $f); 
    my $select_query = getSelectQuery("cdd", $f); 
    $dbh->do('begin');
    my $insth= $dbh->prepare( $insert_query );
    my $selth= $dbh->prepare( $select_query );
    my $n = 0;
    foreach my $l (@lines){
        next if $l =~ /^#/;
        chomp $l;
        next if not $l;
        next if $l =~ /^Query/;
        my @s = split("\t", $l); 
        if ($s[0] =~ /Q\#\d+ - (\S+)/){
            my $u = $1;
            my $name = $uniprot_to_genename{$u};
            my @coords = sort { $a <=> $b } 
                         map { s/^[A-Z]+//; $_ } 
                         split(",", $s[3]);
            my @values = (
                $u,
                $name, 
                "Feature",
                $s[2],
                $s[3],
                $coords[0],
                $coords[-1],
            ) ;
            $n += addRow($insth, \@values, $selth);
            if ($n == $max_commit ){
                $dbh->do('commit');
                $dbh->do('begin');
                $n = 0;
            }
        }else{
            informUser("WARNING: Could not parse CDD result:\n$l\n");
        }
    }
    $dbh->do('commit');
}

#########################################################
sub parseCddHits{
    my $data = shift;
    my $f = shift; 
    my @lines = split("\n", $data);
    my $insert_query = getInsertionQuery("cdd", $f); 
    my $select_query = getSelectQuery("cdd", $f); 
    $dbh->do('begin');
    my $insth= $dbh->prepare( $insert_query );
    my $selth= $dbh->prepare( $select_query );
    my $n = 0;
    foreach my $l (@lines){
        next if $l =~ /^#/;
        chomp $l;
        next if not $l;
        next if $l =~ /^Query/;
        my @s = split("\t", $l); 
        if ($s[0] =~ /Q\#\d+ - (\S+)/){
            my $u = $1;
            my $name = $uniprot_to_genename{$u};

            my @values = (
                $u,
                $name, 
                "Hit",
                $s[8],
                undef,
                $s[3],
                $s[4],
            ) ;
            $n += addRow($insth, \@values, $selth);
            if ($n == $max_commit ){
                $dbh->do('commit');
                $dbh->do('begin');
                $n = 0;
            }
        }else{
            informUser("WARNING: Could not parse CDD result:\n$l\n");
        }
    }
    $dbh->do('commit');
}

#########################################################
sub retrieveAndOutputCddInfo{
    my $feats = retrieveCddFeatures([keys %uniprot_info], 'feats');
    my $hits  = retrieveCddFeatures([keys %uniprot_info], 'hits');
    my @fields = qw / 
                 UniprotId   
                 symbol 
                 ResultType 
                 Feature 
                 Residues 
                 Start 
                 End
                 /; 
    my %f_to_prop = ( 
                 UniprotId   => "TEXT not null",  
                 symbol      => "TEXT",
                 ResultType  => "TEXT",
                 Feature     => "TEXT",
                 Residues    => "TEXT",
                 Start       => "INT not null",
                 End       => "INT not null",
                 );
    createTable('cdd', \@fields, \%f_to_prop);
    informUser("Adding data retrieved from CDD to 'cdd' table of $opts{d}.\n");
    parseCddFeats($feats, \@fields);
    parseCddHits($hits, \@fields);
}


#########################################################
sub retrieveCddFeatures{
    my $queries = shift; 
    my $tdata = shift;
    $tdata ||= 'feats';
    my $url = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";
    informUser( "Setting up CDD '$tdata' search...\n");
    my %opts = (
      evalue   => 0.001,
      tdata    => $tdata,
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
        informUser( "CDD  '$tdata' search with Request-ID $rid started.\n");
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
                informUser("CDD search still running...\n");
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
sub outputUniprotInfo{
    my @fields = qw / 
                GeneName
                UniprotId
                Start
                End
                Feature
                Note
                /; 
    my %f_to_prop = ( 
                GeneName    => "TEXT",
                UniprotId   => "TEXT",
                Start		=> "INT not null",
                End		    => "INT not null",
                Feature		=> "TEXT not null",
                Note		=> "TEXT",
                 );
    createTable('uniprot', \@fields, \%f_to_prop);
    my @lol = ();
    my $insert_query = getInsertionQuery("uniprot", \@fields); 
    my $select_query = getSelectQuery("uniprot", \@fields); 
    $dbh->do('begin');
    my $insth= $dbh->prepare( $insert_query );
    my $selth= $dbh->prepare( $select_query );
    informUser("Adding data retrieved from Uniprot to 'uniprot' table of $opts{d}.\n");
    my $n = 0; 
    foreach my $id (sort keys %uniprot_info){
        my $name = $uniprot_to_genename{$id};
        foreach my $k ( sort by_unipro_coords keys %{$uniprot_info{$id}}){
            my ($start, $end) = split('-', $k);
            foreach my $f (@{$uniprot_info{$id}->{$k}}){
                my ($feature, $note) = split(/\|/, $f );
                my @values = (
                        $name,
                        $id,
                        $start,
                        $end,
                        $feature,
                        $note,
                    );
                $n += addRow($insth, \@values, $selth);
                if ($n == $max_commit ){
                    $dbh->do('commit');
                    $dbh->do('begin');
                    $n = 0;
                }
            }
        }
    }
    $dbh->do('commit');
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
sub outputEvolutionaryTraceInfo{
    return if not $opts{e};
    my @fields = qw /
                    EnsemblTranscriptID
                    RefSeq_peptide
                    Position
                    WildTypeResidue     
                    MutantResidue   
                    Score    
                /;
    my %f_to_prop = ( 
                EnsemblTranscriptID => "TEXT not null",
                RefSeq_peptide      => "TEXT",
                Position            => "INT not null",
                WildTypeResidue     => "TEXT",
                MutantResidue       => "TEXT",
                Score               => "INT not null",
                );
    createTable('EvolutionaryTrace', \@fields, \%f_to_prop);
    informUser("Adding local Evolutionary Trace data to 'evolutionarytrace' table of $opts{d}.\n");
    my $insert_query = getInsertionQuery("evolutionarytrace", \@fields); 
    my $select_query = getSelectQuery("evolutionarytrace", \@fields); 
    my $insth= $dbh->prepare( $insert_query );
    my $selth= $dbh->prepare( $select_query );
    $dbh->do('begin');
    my $n = 0;
    foreach my $t (keys %transcript_xref){
        if (exists $transcript_xref{$t}->{RefSeq_peptide}){
            foreach my $pep (@{$transcript_xref{$t}->{RefSeq_peptide}}){ 
                my $et_file = "$opts{e}/$pep.pred";
                if (not -e $et_file){
                    informUser("WARNING: Could not find Evolutionary Trace " . 
                      "prediction file for $pep ".
                      "($et_file). No ET values will be calculated for this protein.\n");
                    next;
                }
                open (my $ET, $et_file) or die "Could not open $et_file: $!\n";
                while (my $line = <$ET>){
                    chomp $line;
                    next if not $line;
                    my ($res, $score) = split(/\s+/, $line); 
                    if ($res =~ /([A-Z])(\d+)([A-Z])/){
                        my ($wt, $pos, $mut)  = ($1, $2, $3);
                        my @values = (
                             $t, 
                             $pep,
                             $pos,
                             $wt,
                             $mut,
                             $score,
                        ); 
                        $n += addRow($insth, \@values);
                        if ($n == $max_commit ){
                            $dbh->do('commit');
                            $dbh->do('begin');
                            $n = 0;
                        }
                    }else{
                        informUser("ERROR: Could not parse residue '$res' for $pep! Format not understood.\n");
                    }
                }
            }
        }
    }
    $dbh->do('commit');
}


#########################################################
sub outputTranscriptInfo{
    my @fields = qw /
                    Symbol
                    EnsemblGeneID
                    EnsemblTranscriptID
                    EnsemblProteinID
                    RefSeq_mRNA
                    RefSeq_peptide
                    CCDS
                    Uniprot
                    TranscriptScore
                    TranscriptRank	
                /;
    my %f_to_prop = ( 
                Symbol			    => "TEXT",
                EnsemblGeneID	    => "TEXT not null",
                EnsemblTranscriptID => "TEXT primary key not null",
                EnsemblProteinID    => "TEXT",
                RefSeq_mRNA			=> "TEXT",
                RefSeq_peptide      => "TEXT",
                CCDS			    => "TEXT",
                Uniprot			    => "TEXT",
                TranscriptScore		=> "INT not null",
                TranscriptRank		=> "INT not null",
                );
    createTable('transcripts', \@fields, \%f_to_prop);
    informUser("Adding transcript data retrieved from Ensembl to 'transcripts' table of $opts{d}.\n");
    my $insert_query = getInsertionQuery("transcripts", \@fields); 
    my $select_query = getSelectQuery("transcripts", \@fields); 
    my $insth= $dbh->prepare( $insert_query );
    my $selth= $dbh->prepare( $select_query );
    $dbh->do('begin');
    my $n = 0;
    foreach my $ensg (keys %genes){
        my $rank = 0;
        foreach my $t (sort 
        {
            $transcript_ranks{$ensg}->{$b} <=> 
            $transcript_ranks{$ensg}->{$a} ||
            exists $enst_to_uniprot{$b} <=> 
            exists $enst_to_uniprot{$a}  
        } 
        keys %{$transcript_ranks{$ensg}} ){
            $rank++;
            my $symbol = $genes{$ensg}->{display_name};
            my $score = $transcript_ranks{$ensg}->{$t};
            my $ensp = undef;
            my $ref_mrna = undef;
            my $ref_pep = undef;
            my $ccds = undef;
            my $uniprot = undef;
            if (exists $transcript_xref{$t}->{Translation}){
                $ensp = $transcript_xref{$t}->{Translation};
            }
            if (exists $transcript_xref{$t}->{RefSeq_mRNA}){
                $ref_mrna = join("/", @{ $transcript_xref{$t}->{RefSeq_mRNA} } );
            } 
            if (exists $transcript_xref{$t}->{RefSeq_peptide}){
                $ref_pep = join("/", @{ $transcript_xref{$t}->{RefSeq_peptide} } );
            } 
            if (exists $transcript_xref{$t}->{CCDS}){
                $ccds = join("/", @{ $transcript_xref{$t}->{CCDS} } );
            } 
            if (exists $enst_to_uniprot{$t}){
                $uniprot =  $enst_to_uniprot{$t};
            }
            
            my @values = (
                 $symbol, 
                 $ensg,
                 $t, 
                 $ensp,
                 $ref_mrna,
                 $ref_pep,
                 $ccds,
                 $uniprot,
                 $score,
                 $rank,
            ); 
            $n += addRow($insth, \@values, $selth);
            if ($n == $max_commit ){
                $dbh->do('commit');
                $dbh->do('begin');
                $n = 0;
            }
        }
    }
    $dbh->do('commit');
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
                    RefSeq_peptide
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
        next if $f eq 'Sequence conflict';
        my @dets = split(";", $details);
        foreach my $d (@dets){
            if ($d =~ /Note=(.+)/){
                $f .= "|$1";
            }
        }
        $f = uri_decode($f); 
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
    my $time = strftime( "%H:%M:%S", localtime );
    my $msg = shift;
    print STDERR "[$time] $msg";
}
        

#########################################################
sub usage{
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print <<EOT
    
    usage: $0 -l gene_list.txt -d output.db
           $0 -i gene_symbol/gene_id/transcript_id --d output.db

    Options:

    -l, --list FILE
        Input file of gene/transcript/protein IDs one per line. Any whitespace following the first word of each line will be ignored.
    -i, --ids
        One or more gene/transcript/protein identifiers to look up. May be used instead or as well as --list file.
    -d, --db FILE
        Output file for database. Required.
    -s, --species
        Species to search (only applies to non-Ensembl identifiers). Default = human.
    -?, -h, --help
        Show this help message.


    Information:
    
    Note that even when using transcript or protein IDs, the parent gene will be identified and all transcripts processed.

EOT
;
    exit 1 if $msg;
    exit;
}
