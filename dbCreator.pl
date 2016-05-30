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
use File::Tee qw(tee);
use FindBin;
use lib "$FindBin::Bin/lib/vcfhacks/lib";
use IdParser;
use EnsemblRestQuery;
use VcfReader;

my @gene_ids = ();
my %opts = ( i => \@gene_ids);
GetOptions(
    \%opts,
    'l|list=s',
    'i|id=s{,}',#  => \@gene_ids,
    'd|db=s', #sqlite database output
    's|species=s',
    'e|et_folder=s',#optional evolutionary trace folder containing protein predictions
    'm|hgmd=s', #optional HGMD VEP annotated VCF for creation of HGMD variant table
    'c|clinvar=s', #optional ClinVar VEP annotated VCF for creation of ClinVar variant table
    'a|assembly=s', #optional - specify assembly for HGMD consequences - default = GRCh37
    'x|error_log=s',#optional log file for printing messages
    'q|quiet',
    'h|?|help',
) or usage("Syntax error");

usage() if $opts{h};
usage("Error: a gene list or gene ID must be provided") if not $opts{i} and not $opts{l};
usage("Error: please specify a filename for your database using the --db option!\n") 
    if not $opts{d};
$opts{s} = "human" if not $opts{s};
$opts{a} = "GRCh37" if not $opts{a};

if ($opts{x}){#write to log if specified
    tee(STDERR, '>>', $opts{x});
}
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
my %uniprot_to_ensp = ();
my %enst_to_uniprot = ();
my %mappings = (); #protein pos mappings to prevent redundant lookups
                #key = protein id -> assembly -> start -> end = coordinate 

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

outputClinvarInfo();

outputHgmdInfo();

outputEvolutionaryTraceInfo();

outputUniprotInfo() ; 

retrieveAndOutputCddInfo() ; 

$dbh->disconnect(); 

#########################################################
sub outputHgmdInfo{
    my @hgmd_fields =  qw /
            hgmd_id
            disease
            variant_class
            gene_symbol
            hgvs
    /;
   
    my @get_vep =  qw /
            symbol
            gene
            consequence
            allele
            feature
            hgvsc
            hgvsp
            exon
            intron
            cdna_position
            cds_position
            amino_acids
            protein_position
            lof
            lof_filter
            lof_info
            lof_flags
      /;
     
    my @clin_fields = qw /
                hgmd_id
                disease
                variant_class
                gene_symbol
                hgvs
                assembly           
                chrom              
                pos
                ref                
                alt                
                id             
    /;

    my @csq_fields = qw /
                hgmd_id
                feature
                consequence 
                cdna_position
                cds_position
                protein_position
                amino_acids
                hgvsc              
                hgvsp              
                lof
                lof_filter
                lof_info
                lof_flags
    /;
    my %f_to_prop = ( 
                feature             => "TEXT not null",
                hgmd_id             => "TEXT",
                disease             => "TEXT",
                variant_class       => "TEXT",
                gene_symbol         => "TEXT",
                hgvs                => "TEXT",
                assembly            => "TEXT not null",
                chrom               => "TEXT not null",
                pos                 => "INT not null",
                ref                 => "TEXT not null",
                alt                 => "TEXT not null",
                id                  => "TEXT",
                consequence         => "TEXT not null",
                cdna_position       => "int",
                cds_position        => "int",
                protein_position    => "int",
                amino_acids         => "TEXT",
                hgvsc               => "TEXT",
                hgvsp               => "TEXT",
                lof                 => "TEXT",
                lof_filter          => "TEXT",
                lof_info            => "TEXT",
                lof_flags           => "TEXT",
    );
    informUser("Adding local HGMD data to 'HGMD' and 'HGMD_VEP' table of $opts{d}.\n");
    createTable('HGMD', \@clin_fields, \%f_to_prop);
    createTable('HGMD_VEP', \@csq_fields, \%f_to_prop);
    return if not $opts{m};#if not specified create the empty table and return
    my $FH = VcfReader::openVcf($opts{m}); 
    my @vhead = VcfReader::getHeader($opts{m});
    die "VCF header not OK for $opts{m}\n" if not VcfReader::checkHeader(header => \@vhead);
    my %vep_fields = VcfReader::readVepHeader(header => \@vhead);
    my $clin_insert_query = getInsertionQuery("HGMD", \@clin_fields); 
    my $clin_select_query = getSelectQuery("HGMD", \@clin_fields); 
    my $clin_insth= $dbh->prepare( $clin_insert_query );
    my $clin_selth= $dbh->prepare( $clin_select_query );
    my $csq_insert_query = getInsertionQuery("HGMD_VEP", \@csq_fields); 
    my $csq_select_query = getSelectQuery("HGMD_VEP", \@csq_fields); 
    my $csq_insth= $dbh->prepare( $csq_insert_query );
    my $csq_selth= $dbh->prepare( $csq_select_query );
    $dbh->do('begin');
    my $n = 0;
    while (my $l = <$FH>){
        next if $l =~ /^#/;
        my %var_fields = (assembly => $opts{a}); 
        my @split = split("\t", $l);
        my @vep_csq = VcfReader::getVepFields
        (
            line       => \@split,
            vep_header => \%vep_fields,
            field      => \@get_vep,
        );
        #we shouldn't have multiallelic sites in our HGMD VCF (right?)
        # - so no need to check alleles in VEP CSQ
        my $transcript_found = 0;
        foreach my $csq (@vep_csq){ 
            next if not $csq->{gene};
            next if not (exists $transcript_ranks{$csq->{gene}});
            next if not (exists $transcript_ranks{$csq->{gene}}->{$csq->{feature}});
            if (not $transcript_found){
                #GET HGMD ANNOTATIONS FROM INFO FIELD
                foreach my $f (@hgmd_fields){
                    $var_fields{$f} = VcfReader::getVariantInfoField(\@split, $f);
                }
                $transcript_found = 1; 
            }
            my %fields_for_csq = (hgmd_id => $var_fields{hgmd_id}, assembly => $opts{a});
            foreach my $f (@get_vep){
                $fields_for_csq{$f} = $csq->{$f};
            }
            my @csq_values = map { $fields_for_csq{$_} } @csq_fields;
            $n += addRow($csq_insth, \@csq_values);
            if ($n == $max_commit ){
                $dbh->do('commit');
                $dbh->do('begin');
                $n = 0;
            }
        }
        if ($transcript_found){
            #COLLECT VCF FIELDS (CHROM POS etc)
            foreach my $f (qw /chrom pos ref alt id / ){ 
                $var_fields{$f} =  VcfReader::getVariantField(\@split, uc($f));
            }
            my @clin_values = map { $var_fields{$_} } @clin_fields;
            $n += addRow($clin_insth, \@clin_values);
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
sub outputClinvarInfo{
    my @clinvar_fields =  qw /
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
    /;
   
    my @get_vep =  qw /
            symbol
            gene
            consequence
            allele
            feature
            hgvsc
            hgvsp
            exon
            intron
            cdna_position
            cds_position
            amino_acids
            protein_position
            lof
            lof_filter
            lof_info
            lof_flags
      /;
    my @clin_fields = qw / 
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
                assembly
                chrom              
                pos
                ref                
                alt 
    /;
    my @csq_fields = qw /
                feature
                measureset_id
                consequence 
                cdna_position
                cds_position
                protein_position
                amino_acids
                hgvsc              
                hgvsp              
                lof
                lof_filter
                lof_info
                lof_flags
    /;
    my %f_to_prop = ( 
                feature               => "TEXT not null",
                measureset_id         => "TEXT not null",
                symbol                => "TEXT",
                clinical_significance => "TEXT", 
                review_status         => "TEXT",
                hgvs_c                => "TEXT",
                hgvs_p                => "TEXT",
                all_submitters        => "TEXT",
                all_traits            => "TEXT",
                all_pmids             => "TEXT",
                pathogenic            => "INT not null",
                conflicted            => "INT not null", 
                assembly              => "TEXT not null",
                chrom                 => "TEXT not null",
                pos                   => "INT not null",
                ref                   => "TEXT not null",
                alt                   => "TEXT not null",
                consequence           => "TEXT not null",
                cdna_position         => "int",
                cds_position          => "int",
                protein_position      => "int",
                amino_acids           => "TEXT",
                hgvsc                 => "TEXT",
                hgvsp                 => "TEXT",
                lof                   => "TEXT",
                lof_filter            => "TEXT",
                lof_info              => "TEXT",
                lof_flags             => "TEXT",
    );
    informUser("Adding local ClinVar data to 'ClinVar' and 'ClinVar_VEP' tables of $opts{d}.\n");
    createTable('ClinVar', \@clin_fields, \%f_to_prop);
    return if not $opts{c};
    my $FH = VcfReader::openVcf($opts{c}); 
    my @vhead = VcfReader::getHeader($opts{c});
    die "VCF header not OK for $opts{m}\n" if not VcfReader::checkHeader(header => \@vhead);
    my %vep_fields = VcfReader::readVepHeader(header => \@vhead);
    my $clin_insert_query = getInsertionQuery("ClinVar", \@clin_fields); 
    my $clin_select_query = getSelectQuery("ClinVar", \@clin_fields); 
    my $clin_insth= $dbh->prepare( $clin_insert_query );
    my $clin_selth= $dbh->prepare( $clin_select_query );

    createTable('ClinVar_VEP', \@csq_fields, \%f_to_prop);
    my $csq_insert_query = getInsertionQuery("ClinVar_VEP", \@csq_fields); 
    my $csq_select_query = getSelectQuery("ClinVar_VEP", \@csq_fields); 
    my $csq_insth= $dbh->prepare( $csq_insert_query );
    my $csq_selth= $dbh->prepare( $csq_select_query );
    $dbh->do('begin');
    my $n = 0;
    while (my $l = <$FH>){
        next if $l =~ /^#/;
        my %var_fields = (assembly => $opts{a}); 
        my @split = split("\t", $l);
        my @vep_csq = VcfReader::getVepFields
        (
            line       => \@split,
            vep_header => \%vep_fields,
            field      => \@get_vep,
        );
        #we shouldn't have multiallelic sites in our ClinVar VCF (right?)
        # - so no need to check alleles in VEP CSQ
        my $transcript_found = 0;
        foreach my $csq (@vep_csq){ 
            next if not $csq->{gene};
            next if not (exists $transcript_ranks{$csq->{gene}});
            next if not (exists $transcript_ranks{$csq->{gene}}->{$csq->{feature}});
            #only collect INFO fields if we haven't already
            if (not $transcript_found){
                #GET ClinVar ANNOTATIONS FROM INFO FIELD
                foreach my $f (@clinvar_fields){
                    $var_fields{$f} = VcfReader::getVariantInfoField(\@split, $f);
                } 
                $transcript_found = 1;
            }
            my %fields_for_csq = (measureset_id => $var_fields{measureset_id});
            foreach my $f (@get_vep){
                $fields_for_csq{$f} = $csq->{$f};
            }
            my @csq_values = map { $fields_for_csq{$_} } @csq_fields;
            $n += addRow($csq_insth, \@csq_values);
            if ($n == $max_commit ){
                $dbh->do('commit');
                $dbh->do('begin');
                $n = 0;
            }
        }
        if ($transcript_found){
            #COLLECT VCF FIELDS (CHROM POS etc)
            foreach my $f (qw /chrom pos ref alt / ){ 
                $var_fields{$f} =  VcfReader::getVariantField(\@split, uc($f));
            }
            my @clin_values = map { $var_fields{$_} } @clin_fields;
            $n += addRow($clin_insth, \@clin_values);
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
                         map { s/^[A-Z]+//g; $_ } 
                         split(/[\,\-]/, $s[3]);
            my ($grch37_pos, $grch38_pos ) = ("", "");
            if ( $uniprot_to_ensp{$u} and ref $uniprot_to_ensp{$u} eq 'ARRAY'){
                ($grch37_pos, $grch38_pos ) = genomicPosFromEnsp
                (
                    ids   => $uniprot_to_ensp{$u},
                    start => $coords[0],
                    end   => $coords[-1],
                );
            } 
            my @values = (
                $u,
                $name, 
                "Feature",
                $s[2],
                $s[3],
                $coords[0],
                $coords[-1],
                $grch37_pos, 
                $grch38_pos, 
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
            my ($grch37_pos, $grch38_pos ) = ("", "");
            if ( $uniprot_to_ensp{$u} and ref $uniprot_to_ensp{$u} eq 'ARRAY'){
                ($grch37_pos, $grch38_pos ) = genomicPosFromEnsp
                (
                    ids   => $uniprot_to_ensp{$u},
                    start => $s[3],
                    end   => $s[4],
                ); 
            }
            my @values = (
                $u,
                $name, 
                "Hit",
                $s[8],
                undef,
                $s[3],
                $s[4],
                $grch37_pos, 
                $grch38_pos, 
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
                 GRCh37Pos   
                 GRCh38Pos  
                 /; 
    my %f_to_prop = ( 
                 UniprotId   => "TEXT not null",  
                 symbol      => "TEXT",
                 ResultType  => "TEXT",
                 Feature     => "TEXT",
                 Residues    => "TEXT",
                 Start       => "INT not null",
                 End       => "INT not null",
                 GRCh37Pos   => "TEXT",
                 GRCh38Pos   => "TEXT",
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
                GRCh37Pos
                GRCh38Pos
                /; 
    my %f_to_prop = 
    ( 
        GeneName    => "TEXT",
        UniprotId   => "TEXT",
        Start		=> "INT not null",
        End		    => "INT not null",
        Feature		=> "TEXT not null",
        Note		=> "TEXT",
        GRCh37Pos   => "TEXT",
        GRCh38Pos   => "TEXT",
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
                my ($grch37_pos, $grch38_pos ) = ("", "");
                if ($uniprot_to_ensp{$id} and ref $uniprot_to_ensp{$id} eq 'ARRAY'){
                    ($grch37_pos, $grch38_pos ) = genomicPosFromEnsp
                    (
                        ids   => $uniprot_to_ensp{$id},
                        start => $start,
                        end   => $end,
                    ); 
                }
                my @values = 
                (
                    $name,
                    $id,
                    $start,
                    $end,
                    $feature,
                    $note,
                    $grch37_pos,
                    $grch38_pos,
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
sub genomicPosFromEnsp{
    my %args = @_;
    my (%grch37_pos, %grch38_pos );
    informUser
    (
        "Mapping GRCh37 and GRCh38 coordinates from protein positons ".
        "$args{start}-$args{end} for ". join(",", @{$args{ids}}) . "...\n"
    );
    $restQuery->useGRCh37Server(); 
    foreach my $id (@{$args{ids}}){
        if (exists $mappings{$id}->{GRCh37}->{$args{start}}->{$args{end}}){
            my $regions = $mappings{$id}->{GRCh37}->{$args{start}}->{$args{end}};
            $grch37_pos{$regions} = undef;
        }elsif (my $regions = getGenomicRegionsFromEnsp
            (
                id       => $id,
                start    => $args{start},
                end      => $args{end},
                assembly => 'GRCh37',
            )
        ){
            $grch37_pos{$regions} = undef;
            $mappings{$id}->{GRCh37}->{$args{start}}->{$args{end}} = $regions;
        }
    }
    $restQuery->useDefaultServer(); 
    foreach my $id (@{$args{ids}}){
        if (exists $mappings{$id}->{GRCh38}->{$args{start}}->{$args{end}}){
            my $regions = $mappings{$id}->{GRCh38}->{$args{start}}->{$args{end}};
            $grch38_pos{$regions} = undef;
        }elsif (my $regions = getGenomicRegionsFromEnsp
            (
                id       => $id,
                start    => $args{start},
                end      => $args{end},
                assembly => 'GRCh38',
            )
        ){
            $grch38_pos{$regions} = undef;
            $mappings{$id}->{GRCh38}->{$args{start}}->{$args{end}} = $regions;
        }
    }
    if (keys %grch37_pos > 1){
        informUser(
            "WARNING: More than one set of genomic positions found for ".
            join ("/", @{$args{ids}} ) . " $args{start}-$args{end} for GRCh37. Found :\n".
            join("\n", keys %grch37_pos) . "\n"
        );
    }
    if (keys %grch38_pos > 1){
        informUser(
            "WARNING: More than one set of genomic positions found for ".
            join ("/", @{$args{ids}} ) . " $args{start}-$args{end} for GRCh38. Found :\n".
            join("\n", keys %grch38_pos) . "\n"
        );
    }
    if (not keys(%grch37_pos)){
        informUser(
            "WARNING: No genomic positions found for ".
            join ("/", @{$args{ids}} ) . " $args{start}-$args{end} for GRCh37.\n"
        );
        if (keys %grch38_pos){
            informUser(
                "WARNING: Attempting to map found GRCh38 coordinates to GRCh37\n"
            );
            if (my $regions = convertAssemblyCoordinates
                (
                    "GRCh38",
                    "GRCh37",
                    [keys %grch38_pos],
                )
            ){
                $grch37_pos{$regions} = undef;
            }
        }
    }
    if (not keys(%grch38_pos)){
        informUser(
            "WARNING: No genomic positions found for ".
            join ("/", @{$args{ids}} ) . " $args{start}-$args{end} for GRCh38.\n"
        );
        if (keys %grch37_pos){
            informUser(
                "WARNING: Attempting to map found GRCh37 coordinates to GRCh38\n"
            );
            if (my $regions =convertAssemblyCoordinates 
                (
                    "GRCh37",
                    "GRCh38",
                    [keys %grch38_pos],
                )
            ){
                $grch38_pos{$regions} = undef;
            }
        }
    }
    return 
    (
        join(";", keys %grch37_pos), 
        join(";", keys %grch38_pos), 
    );
}

#########################################################
sub convertAssemblyCoordinates{
    my ($original, $new, $regions) = @_;
    my @pos = (); 
    foreach my $r (@$regions){
        my $endpoint = "map/human/$original/$r:1/$new";
        my $coords = $restQuery->queryEndpoint($endpoint);
        if (ref $coords eq 'HASH'){
            foreach my $map (@{$coords->{mappings}}){
                push @pos, sprintf
                (   "%s:%s-%s",
                    $map->{mapped}->{seq_region_name},
                    $map->{mapped}->{start},
                    $map->{mapped}->{end},
                ); 
            }
        }
    }
    return join(";", @pos);
}
#########################################################
sub getGenomicRegionsFromEnsp{
    my %args = @_;
    my @pos = (); 
    my $coords = $restQuery->proteinPosToGenomicPos(%args);
    if (not exists $coords->{mappings} or 
        ref $coords->{mappings} ne 'ARRAY'){
        informUser("WARNING: Could not map genomic coordinates ".
                   "for $args{id} $args{start}-$args{end} ($args{assembly}) - server ".
                   "did not return an ARRAY reference.\n"
        );
    }else{
        foreach my $mapping (@{$coords->{mappings}}){
            if ($mapping->{assembly_name} ne $args{assembly}){
                informUser("WARNING: Unexpected assembly name (".
                           $mapping->{assembly_name} . ") for coordinates ".
                           "retrieved for $args{id} $args{start}-$args{end}. ".
                           "Expected assembly '$args{assembly}'.\n"
                );
            }else{
                my $reg = sprintf("%s:%s-%s",
                                  $mapping->{seq_region_name},
                                  $mapping->{start},
                                  $mapping->{end},
                ); 
                push @pos, $reg;
            }
        }
    }
    return join(";", @pos);
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
            addXrefs($t, $genes{$ensg}); 
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
    my $enst = shift; 
    my $ensg_ts = shift;
    my $id = $enst->{id};
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
                        getUniprotData($ext->{primary_id}, $ensg_ts);
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
    my $ensg_ts = shift;
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
    my ($name, $length, @t)  = parseUniprotFlat($txt, $id);
    if (not @t){
        my $matched_length = 0;
        if ($length){
            foreach my $enst (@{$ensg_ts->{Transcript}}){
                if ($enst->{Translation}){
                    if ($enst->{Translation}->{length} == $length){
                        $matched_length++;
                        $enst_to_uniprot{$enst->{id}} = $id;
                        push @{$uniprot_to_ensp{$id}},  $enst->{Translation}->{id};
                        informUser(
"WARNING: Could not identify Ensembl transcript for canonical isoform of $id " .
"- basing match of $enst->{id} on basis of shared translation length.\n");
                    }
                }
            }
        }
        if (not $matched_length){
            informUser ("WARNING: Could not identify Ensembl transcript for canonical isoform of $id!\n");
        }
    }
    foreach my $t (@t){
        $enst_to_uniprot{$t} = $id;
        my @ensts = grep { $_->{id} eq $t } @{$ensg_ts->{Transcript}}; 
        if (@ensts){
            push @{$uniprot_to_ensp{$id}},  map { $_->{Translation}->{id} }  @ensts;
        }
    }
    $uniprot_info{$id} = parseUniprotGff($gff);
    $uniprot_to_genename{$id} = $name;
}

#########################################################
sub parseUniprotFlat{
    my $txt = shift;
    my $u_id = shift;
    my @lines = split("\n", $txt);
    my $canonical; 
    my $name = '';
    my $u_length = 0;
    my @transcripts = ();
    for (my $i = 0; $i < @lines; $i++){
        if ($lines[$i] =~ /^ID\s+\S+\s+.*\s(\d+)\s+AA\./){
            $u_length = $1;
            next;
        }
        if ($lines[$i] =~ /^GN\s+Name=(\S+);/){
            $name = $1;
            next;
        }
        if ($lines[$i] =~ /^CC\s+IsoId=((\S+)(,\s+\S+)*);/ ){
            my @iso_ids = split(/\,\s+/, $1);
            if ($lines[$i] =~ /^CC\s+IsoId=(\S+)(,\s+\S+)*;\s+Sequence=Displayed;/){
                ($canonical) = grep { $_ =~ /^$u_id-\d+$/ } @iso_ids;
                next;
            }elsif($lines[$i + 1] =~ /Sequence=Displayed/){ 
                ($canonical) = grep { $_ =~ /^$u_id-\d+$/ } @iso_ids;
                next;
            }
        }
        if ($canonical){
            if($lines[$i] =~ /^DR\s+Ensembl;\s+(ENST\d+);\s+ENSP\d+;\s+ENSG\d+\. \[$canonical\]/){
                push @transcripts, $1;
            }
        }elsif($lines[$i] =~ /^DR\s/){#if only one isoform there won't be IsoId's
            if($lines[$i] =~ /^DR\s+Ensembl;\s+(ENST\d+);\s+ENSP\d+;\s+ENSG\d+/){
                push @transcripts, $1;
                
            }
        }
    }
    return ($name, $u_length, @transcripts);
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
        my @lookups = ();
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
            if ($transcript and ref $transcript eq 'HASH'){
                if (exists $transcript->{id}){
                    $gene_hash = geneFromEnst($transcript->{id});
                }
            }else{
                informUser( "WARNING: No transcript identified for ID \"$g\"\n");
            }
        }else{
            informUser("Identifying Ensembl gene via gene cross-reference...\n");
            my $gene = $restQuery->getGeneViaXreg($g, $opts{s});
            if (ref $gene eq 'ARRAY'){
                #go through each item to see if it matches the display-name
                foreach my $ge (@$gene){
                    if ($ge->{id}){
                        my $ge_hash = $restQuery->lookUpEnsId($ge->{id}, 1);
                        if (uc($ge_hash->{display_name}) eq uc($g)){
                        #if gene symbol matches then we use this entry
                            $gene_hash = $ge_hash;
                            last;
                        }else{
                            push @lookups, $ge_hash;
                        }
                    }
                }
                if (not $gene_hash){
                    if (@lookups == 1){
                        $gene_hash = $lookups[0];
                    }
                }
            }
        }
        if (not $gene_hash){
            informUser("WARNING: Could not identify gene for ID \"$g\"\n");
            if (@lookups){
                my $idstring = join("\n", map { $_->{display_name} } @lookups );
                informUser
                (
                    "Identified the following non-matching display names:\n".
                    "$idstring\n"
                );
            }
        }else{
            informUser
            (
                "Found gene with display name " . 
                $gene_hash->{display_name} . " and Ensmbl gene ID ".
                $gene_hash->{id} . " for input '$g'.\n"
            );
            if (exists $id_mapping{$gene_hash->{id}}){
                informUser
                (
                    "WARNING: Duplicated gene identified: $gene_hash->{id} ".
                    "identified from $g and from ".
                    join(", ", @{$id_mapping{$gene_hash->{id}}}) . "\n"
                );
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
