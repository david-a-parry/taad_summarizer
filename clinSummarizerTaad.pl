#!/usr/bin/env perl 
#David A. Parry, June 2015

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Term::ProgressBar;
use POSIX qw/strftime/;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use List::Util qw(first sum);
use Pod::Usage;
use File::Basename;
use FindBin;
use Tabix; 
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use HTTP::Tiny;
use JSON;
use lib "$FindBin::Bin/lib";
use VcfReader;
use TextToExcel;

my $progressbar;
my @allele_balance = ();
my %seq_cache = ();
my %opts = (b => \@allele_balance);
GetOptions(
    \%opts,
    'i|input=s',        #vcf input
    'o|output=s',       #xlsx output
    'n|no_blanks',      #no blank samples
    'm|hgmd=s',         #vcf of HGMD variations converted with hgmdMartToVcf.pl
    't|transcripts=s',  #optional tsv list of transcript IDs in order of priority
    'c|clinvar=s',      #optional ClinVar VCF to add ClinVar CLINSIG annotations
    'b|allele_balance=f{,}', #min and optional max alt allele ratio per sample call
    'd|depth=i',         #optional min depth for sample call
    'a|allele_cutoff=f', #remove if allele is present in this proportion or more calls
    'y|frequency=f',     #filter on this allele frequency if dbSNP/EVS/ExAC annot. found
    'f|filter_output=s', #optional output file for calls filtered on allele_cutoff
    'r|rest_server=s',   #URL of REST server to use if not the default (http://grch37.rest.ensembl.org)
    'g|gq=f',            #min GQ quality for calls
    's|scan_gxy',       #scan for GXY sequences throughout sequence rather than relying on VEP annotation
    'primer_file=s',      #read a tab delimited file of primers and coordinates to assign primer IDs to variants
    't|transcripts=s',  #optional tsv list of transcript IDs plus links to uniprot IDs in order of priority
    'rules=s',        #optional tsv file of mutation rules
    'u|uniprot_data=s', #optional uniprot data file
    'c|cdd_data=s',     #optional CDD data file
    'phenotype_file=s',   #read a csv file of patient IDs and phenotype information
    'p|progress',       #show a progress bar?
    'x|do_not_merge',    #do not merge cells
    '1|single_sheet',
    'v|verbose',        #print extra progress information
    'h|?|help',
    'manual',
) or pod2usage(-exitval => 2, -message => "Syntax error.\n"); 

pod2usage( -verbose => 1 ) if $opts{h};
pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -exitval => 2, -message => "-i/--input is required" ) if (not $opts{i});
pod2usage( -exitval => 2, -message => "-m/--hgmd is required" ) if (not $opts{m});

#create our http client for ensembl REST queries 
my $http = HTTP::Tiny->new();
my $server = $opts{rest_server} ? $opts{rest_server} : 'http://grch37.rest.ensembl.org';
my $min_gq = defined $opts{g} ? 0 : $opts{g}; #default min GQ of 0
my $af = defined $opts{y} ? 0.001 : $opts{y}; #default of 0.1 % cutoff in external databases

#open VCF, get samples and get VEP annotations
informUser("Checking input VCF\n");
my @vhead              = VcfReader::getHeader($opts{i});
die "VCF header not OK for $opts{i}\n" if not VcfReader::checkHeader(header => \@vhead);
my %samples_to_columns = VcfReader::getSamples(header => \@vhead, get_columns => 1);
my %vep_fields         = VcfReader::readVepHeader(header => \@vhead);
my %info_fields = VcfReader::getInfoFields(header => \@vhead);
my %af_info_fields     = getAfAnnotations(\@vhead);
my $total_var;

#check primer_file if supplied
my @primers = (); #create an array of hash refs with id, chr, start, end, f and r keys
if ($opts{primer_file}){
    informUser("Checking primer file\n");
    @primers = readPrimerFile();
}

#check validations_file if supplied
my %validated = (); #key is variant unique ID, value is 1 for validated, 2 for did not validate, 0 for not done
if ($opts{validations_file}){
    informUser("Checking validations file\n");
    %validated = readValidationsFile();
}
#check phenotype_file if supplied
my $phenotypes ; # hash ref key is sample ID, value is anon hash of phenotype fields to values
my @pheno_fields = (); #list of phenotype fields found
if ($opts{phenotype_file}){
    informUser("Checking phenotype file\n");
    ($phenotypes, @pheno_fields) = readPhenotypeFile();
    @pheno_fields = sort @pheno_fields;
}
#check HGMD VCF and retrieve search arguments 

informUser("Checking HGMD vcf\n");
my %search_args = getSearchArgs();

# check ClinVar TSV file (https://github.com/macarthur-lab/clinvar) 
# if supplied and retrieve search arguments 

informUser("Checking ClinVar\n");
my %clinvar_sargs = getClinVarSearchArgs();
my %clnsig_codes = getClnSigCodes();

#set consequence ranks;
my %so_ranks = ();
setConsequenceRanks();

#open and check transcripts list

my %transcript_ranks = ();
my %enst_to_uniprot = ();
setTranscriptsRanks();

#read uniprot data
my %uniprot_features = ();
readUniprotFile();

#read CDD data
my %cdd_features = ();
readCddFile();

#check rules

my %rules = ();
setMutationRules();


#setup our hash of headers

my %headers = getHeaders();

#setup output
#Excel related global variables
my %sheets = ();
my $header_formatting;    
my $std_formatting;     
my $url_format;
my $xl_obj;
#filehandle for --allele_cutoff filtered variants
informUser("Preparing output file\n");
my $FILTER_OUT;
setupOutput();

#set up counts for variant indexes of each sheet
my %variant_counts = map { $_ => 0 } qw / HGMD LOF DamagingMissense BenignMissense Other / ;

#store variants per sample in this hash for writing sample sheet
my %sample_vars = (); 

#set up progress bar if user wants one
my $next_update = 0;
my $n = 0;
if ($opts{p}){
    informUser("Calculating file length for progress bar...\n");
    $total_var = VcfReader::countVariants( $opts{i} );
    informUser("$opts{i} has $total_var variants.\n");
    $progressbar = Term::ProgressBar->new(
        { name => "Analyzing", 
          count => ($total_var), 
          ETA => "linear" 
        } 
    );
    $progressbar->rbrack('>');
}

#get filehandle and start reading variants
my $FH = VcfReader::_openFileHandle($opts{i}); 
informUser("Commencing variant analysis\n");
while (my $l = <$FH>){
    next if $l =~ /^#/;
    assessAndWriteVariant($l);
    $n++;
    if ($progressbar){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
}
close $FH;
if ($progressbar){
    $progressbar->update($total_var) if $total_var >= $next_update;
}
if ($FILTER_OUT){
    close $FILTER_OUT;
}
$progressbar = undef;
informUser("Done writing variant sheets - writing sample summary...\n"); 
writeSampleSummary();
informUser("Done.\n");
if ($xl_obj){
    $xl_obj->DESTROY();
}

###########################################################
sub getAfAnnotations{
    my $h = shift;
    return if not $af;
    my %af_found = ();
    my @af_fields =  qw ( 
        AS_CAF
        AS_G5A
        AS_G5
        AS_COMMON
        EVS_EA_AF
        EVS_AA_AF
        EVS_ALL_AF
    );
    foreach my $key (keys %info_fields){
        my $warning = <<EOT
[WARNING] Found expected frequency annotation ($key) in INFO fields, but 'Number' field is $info_fields{$key}->{Number}, expected 'A'. Ignoring this field.
EOT
;
        my $info = <<EOT
[INFO] Found allele frequency annotation: $key. This will be used for filtering on allele frequency.
EOT
;
        if (grep { $key eq $_ } @af_fields){
            if ($info_fields{$key}->{Number} ne 'A'){
                print STDERR $warning;
            }else{
                print STDERR $info;
                $af_found{$key} = $info_fields{$key};
            }
        }else{
            if ($key =~ /^FVOV_AF_\S+$/){
                if ($info_fields{$key}->{Number} ne 'A'){
                    print STDERR $warning;
                }else{
                    print STDERR $info;
                    $af_found{$key} = $info_fields{$key};
                }
            }
        }
    }
    return %af_found;
}
 
###########################################################
sub writeSampleSummary{
    my $sample_sheet = $xl_obj->get_worksheets()->[$sheets{SampleSummary}];
    my $row = 1;
    my %most_damaging = ();#keys are categories, values are array of sample variant array refs
    foreach my $s (
            sort {$samples_to_columns{$a} <=> $samples_to_columns{$b}} 
            keys %samples_to_columns
    ){
        #write out numbers of each variant class per sample
        my $col = 0;
        foreach my $field ( @{$headers{sample}} ){
            if ($field eq 'Sample'){
                $sample_sheet->write($row, $col, $s)
            }elsif (exists $sample_vars{$s}->{$field} ){
                $sample_sheet->write($row, $col, scalar @{$sample_vars{$s}->{$field}});
                #$sample_sheet->write_comment($row, $col,  join("\n", @{$sample_vars{$s}->{$field}}), visible => 0);
            }else{
                $sample_sheet->write($row, $col, 0);
            }
            $col++;
        }
        $row++;

        #determine most damaging consequence category for this sample
        # and store in most damaging to allow output in order of most 
        # likely pathogenic 
        my @csq = sort byCsqClass keys %{$sample_vars{$s}}; 
        if (not @csq){ #no variants for this sample
            push @{ $most_damaging{NONE} }, [0, $s]; 
        }else{
            if (@{$sample_vars{$s}->{$csq[0]}} > 1){ 
            #if multiple variants with most damaging category chose the highest CADD score
            #CADD score will be first entry in array ref
            
                @{$sample_vars{$s}->{$csq[0]}} = sort {$b->[0] <=> $a->[0]} @{$sample_vars{$s}->{$csq[0]}};
            }
            push @{ $most_damaging{$csq[0]} }, $sample_vars{$s}->{$csq[0]}->[0]; 

        }
    }
    foreach my $k (sort byCsqClass keys %most_damaging){
        #order each category by CADD score
        #CADD score will be first entry in array ref
        foreach my $var ( sort {$b->[0] <=> $a->[0]} @{ $most_damaging{$k} } ){
            #remove cadd score from beginning of line
            shift @$var; 
            $xl_obj->writeLine
            (
                line => $var,
                worksheet  => $sheets{MostDamaging},
            );
        }        
    }
}

###########################################################
sub byCsqClass{
    if ($a eq $b){
        return 0;
    }
    foreach my $class (qw / 
            HGMD_DM 
            ClinVarPathogenic 
            LOF 
            CollagenGlyXY 
            HGMD_other 
            DamagingMissense 
            BenignMissense
            Other
            /){

        if ($a eq "$class"){
            return -1;
        }elsif ($b eq "$class"){
            return 1;
        }
    }
    return 0;
}
        



###########################################################
sub getClinVarSearchArgs{
    return if not $opts{c};
    #use bgzip compressed version of clinvar file from https://github.com/macarthur-lab/clinvar
    my ($bgz, $clinVarCols) = checkClinvarFile($opts{c});
    my $index = "$bgz.tbi";
    my $iterator = Tabix->new(-data =>  $bgz, -index => $index) ;
    my %sargs = ( tabix_iterator => $iterator, columns => $clinVarCols ); 
    return %sargs;
}

###########################################################
sub readPhenotypeFile{
    my %pheno = (); #key is variant unique ID, value is anon hash of phenotype info
    open (my $PHENO, $opts{phenotype_file}) 
      or die "Can't open phenotype file $opts{phenotype_file}: $!\n";
    chomp (my $header= <$PHENO>);
    $header =~ s/^#+//;
    my %cols = getColumns($header, 1, ',');
    if (not exists $cols{sample_id}){
        die "Could not find required column 'sample_id' in phenotype file!\n";
    }
    while ( my $line = <$PHENO> ){
        chomp $line;
        my @split = split(",", $line); 
        foreach my $f (sort keys %cols){
            next if $f eq 'sample_id';
            $pheno{$split[ $cols{sample_id} ] } ->{$f} = $split[ $cols{$f} ] ;
        }
    }
    close $PHENO;
    delete $cols{sample_id}; 
    return \%pheno, keys %cols;
}

###########################################################
sub readValidationsFile{
    my %val = (); #key is variant unique ID, value is 1 for validated, 2 for did not validate, 0 for not done
    open (my $VAL, $opts{validations_file}) 
      or die "Can't open validations file $opts{validations_file}: $!\n";
    my $header= <$VAL>;
    my %cols = getColumns($header);
    my @fields =  qw / UID VALIDATED / ;
    foreach my $f (@fields){
        if (not exists $cols{$f}){
            die "Could not find required column '$f' in validation file!\n";
        }
    }
    while ( my $line = <$VAL> ){
        chomp $line;
        my @split = split("\t", $line); 
        $val{$split[ $cols{UID} ] } = $split[ $cols{VALIDATED} ] ;
    }
    close $VAL;
    return %val;
}
###########################################################
sub readPrimerFile{
    open (my $PRIMERS, $opts{primer_file}) 
      or die "Can't open primer file $opts{primer_file}: $!\n";
    my @prime = (); #create an array of hash refs with id, chr, start, end, f and r keys
    my $header = <$PRIMERS>;
    my %cols = getColumns($header); 
    my @fields =  qw / ID CHROM START END F R / ;
    foreach my $f (@fields){
        if (not exists $cols{$f}){
            die "Could not find required column '$f' in primers file!\n";
        }
    }
    while ( my $line = <$PRIMERS> ){
        chomp $line;
        my @split = split("\t", $line); 
        my %pr = map { $_ => $split[ $cols{$_} ] } @fields;
        push @prime, \%pr; 
    }
    close $PRIMERS;
    @prime = sort {   
            $a->{CHROM} cmp $b->{CHROM} ||
            $a->{START} <=> $b->{START} || 
            $a->{END} <=> $b->{END}     || 
            $a->{ID} cmp $b->{ID}
    } @prime;
    return @prime;
}
###########################################################
sub checkClinvarFile{
#returns columns hash of col name to col number
    my $cv = shift;
    if ($cv !~ /\.(b)*gz$/){
        print STDERR "ClinVar file ($cv) is missing a .gz extension - " .
            "attempting to compress with bgzip.\n";
        return compressClinVarTsv($cv); 
    }elsif( not -e "$cv.tbi"){
        print STDERR "No index found for ClinVar file ($cv). Attempting to index with tabix...\n";
        return indexClinVar($cv);
    }elsif( -M "$cv.tbi" > -M $cv){
        print STDERR "Tabix index is older than ClinVar file ($cv). Attempting to re-index with tabix...\n";
        return indexClinVar($cv);
    }else{
        my %col = getClinVarColumns($cv);
        return $cv, \%col;
    }
}
###########################################################
sub compressClinVarTsv{
    my $cv = shift;
    if (not `which bgzip`){
        die "Could not find bgzip executable - please install bgzip and ensure ".
            "it is in your PATH or compress and index $cv manually.\n";
    }
    my %columns = getClinVarColumns($cv);
    print STDERR "Compressing $cv with bgzip...\n";
    system("bgzip -c $cv > $cv.gz"); 
    checkExit($?);
    $cv = "$cv.gz";
    indexClinVar($cv, \%columns); 
    return ($cv, \%columns);
}
###########################################################
sub indexClinVar{
    my $cv = shift;
    my %columns = ();
    if (my $col = shift){
        %columns = %$col;
    }else{
        %columns = getClinVarColumns($cv);
    }
    my $seqcol = $columns{chrom} + 1;
    my $poscol = $columns{pos}   + 1;
    system("tabix -S 1 -s $seqcol -b $poscol -e $poscol $cv"); 
    checkExit($?);
    return %columns;
}
###########################################################
sub checkExit{
    my $exit = shift;
    return if not $exit;
    if ($exit == -1) {
        print "failed to execute: $!\n";
    }elsif ($exit & 127) {
        printf "command died with signal %d, %s coredump\n",
        ($exit & 127),  ($exit & 128) ? 'with' : 'without';
    }elsif ($exit) {
        printf "command exited with value %d\n", $exit >> 8;
    }
    die "Error executing command. Exiting.\n";
}
    
###########################################################
sub getClinVarColumns{
    my $cv = shift;
    my $CVAR; 
    if ($cv =~ /\.(b)*gz$/){
        $CVAR = new IO::Uncompress::Gunzip $cv, MultiStream => 1 
          or die "IO::Uncompress::Gunzip failed while opening $cv for reading: \n$GunzipError";
    }else{
        open ($CVAR, $cv) or die "Can't open $cv for reading header: $!\n";
    }
    my $header = <$CVAR>;
    my %columns = getColumns($header, 1);
    for my $c ( qw / 
                    chrom
                    pos
                    ref
                    alt
                    mut
                    clinical_significance
                    pathogenic
                    all_traits
                    all_pmids
                / ) { 
        if (not exists $columns{$c}){
            die "Required column ($c) not found in ClinVar file ($cv) header.\n";
        }
    }
    close $CVAR;
    return %columns;
}
###########################################################
sub getSearchArgs{
    die "Header not ok for input ($opts{m}) "
      if not VcfReader::checkHeader( vcf => $opts{m} );
    return VcfReader::getSearchArguments($opts{m});
}
    
###########################################################
sub assessAndWriteVariant{
#identify transcript to report (make sure variant is functional and choose highest ranked transcript)
#write variant details to HGMD Excel sheet
#add var info to sample summary sheet and link to variant in main worksheet;
    my $line = shift;
    chomp ($line); 
    my @split = split("\t", $line); 
    #see if variant overlaps any variants in HGMD file
    my %min= VcfReader::minimizeAlleles(\@split);
    foreach my $al (sort {$a<=>$b} keys %min){
        if (%af_info_fields){ #check for annotateSnps.pl etc. frequencies
           my %af_info_values = getAfInfoValues(\@split);
           next if ( alleleAboveMaf($al - 1, \%af_info_values) );
        }
        my @hgmd_matches = searchHgmd($min{$al});
        my @clinvar_matches = searchClinVar($min{$al});
        writeToSheet(\@split, $min{$al}, \@hgmd_matches, \@clinvar_matches);
    }
}

###########################################################
sub getAfInfoValues{
    my $l = shift;
    my %values = ();
    foreach my $k (keys %af_info_fields){
        $values{$k} = VcfReader::getVariantInfoField($l, $k);
    }
    return %values;
}

###########################################################
sub alleleAboveMaf{
    my $i = shift; #1-based index of alt allele to assess
    my $af_values = shift; #hash ref of INFO fields to their values
    foreach my $k (keys %{$af_values}){
        next if not $af_values->{$k};
        next if $af_values->{$k} eq '.';
        my @split = split(",", $af_values->{$k}); 
        next if $split[$i] eq '.';
        if ($k eq "AS_G5" or $k eq "AS_G5A"){
            if ($af <= 0.05 and $split[$i] > 0){
                return 1;
            }
        }elsif ($k eq "AS_COMMON"){
            if ($af <= 0.01 and $split[$i] > 0){
                return 1;
            }
        }else{#should be a standard allele freq now
            if ($info_fields{$k}->{Type} eq 'Float' or $k eq 'AS_CAF'){
                return 1 if $split[$i] >= $af;
            }else{
                print STDERR "WARNING: Don't know how to parse INFO field: $k.\n";
            }
        }
    }
    return 0;
}

###########################################################
sub writeToSheet{
    my ($l, $var, $matches, $clinvar_matches) = @_;
    #split line, minimized variant hash, HGMD matching lines, ClinVar matching lines
    
    #filter line if allele is too common
    if ($opts{a}){
        my %allele_counts = VcfReader::countAlleles
            (
                  line => $l,
                  minGQ => $min_gq,
            );
        my $total_allele_count = 0;
        foreach my $k (keys %allele_counts){
            $total_allele_count += $allele_counts{$k}; 
        }
        if ($total_allele_count > 0){
            if ($allele_counts{$var->{ALT_INDEX}}/$total_allele_count >= $opts{a}){
                print $FILTER_OUT join("\t", @$l) . "\n" if $opts{f};
                return;
            }
        }
    }
    my @row = (); #values to write out to excel sheet
    my @split_cells = (); #values to write in split cells spanning row
    my @primer_hits = (); #primers that could PCR variant if $opts{primer_file}
    my %hgmd = ();#keys are HGMD fields, values are array refs of values
    #collect HGMD database annotations if found
    my @hgmd_fields = 
    ( qw /
            hgmd_id
            disease
            variant_class
            gene_symbol
            hgvs
        /
    );
    if (@$matches){
        foreach my $h (@$matches){
            my @match = split("\t", $h);
            foreach my $f (@hgmd_fields){
                push @{$hgmd{$f}}, 
                        VcfReader::getVariantInfoField(\@match, $f);
            }
        }
        foreach my $f (@hgmd_fields){
            push @row, join(",", @{$hgmd{$f}});
        }
    }else{
        foreach my $f (@hgmd_fields){
            push @row, join("-");
        }
    }
        
    #collect ClinVar database annotations if found
    my ($clnSig, $clnTraits, $cvarPathognic, $cvarConflicted) 
                         = getClinSig($clinvar_matches, $var);
    push @row, $clnSig, $clnTraits; 
    push @row, $cvarConflicted;
    
    #add standard VCF fields from this allele only
    foreach my $f ( qw / CHROM ORIGINAL_POS ORIGINAL_REF ORIGINAL_ALT / ){
        push @row, $var->{$f};
    }
    push @row, VcfReader::getVariantField($l, 'ID');
    push @row, VcfReader::getVariantField($l, 'QUAL');
    push @row, VcfReader::getVariantField($l, 'FILTER');
    #deal with VEP annotations
    my @get_vep =  
      qw /
            Symbol
            Consequence
            Allele
            feature
            canonical
            hgvsc
            hgvsp
            polyphen
            sift
            DOMAINS
            Amino_acids
            protein_position
            ensp
            LoF
            LoF_Filter
            LoF_info
            LoF_flags

      /;

    my @vep_csq = VcfReader::getVepFields
    (
        line       => $l,
        vep_header => \%vep_fields,
        field      => \@get_vep,
    );
    my $ref      = VcfReader::getVariantField($l, 'REF');
    my @alts     = split(",", VcfReader::getVariantField($l, 'ALT'));
    my @vep_alts = VcfReader::altsToVepAllele
    (
        ref  => $ref,
        alts => \@alts,
    );
    my %alt_to_vep = ();
    @alt_to_vep{@alts} = @vep_alts;
    my @csq_to_rank = ();
    foreach my $csq (@vep_csq){ 
        #collect info for each transcript before selecting which to use in report
        #skip consequences for different alleles if variant site is multiallelic
        next if $csq->{allele} ne $alt_to_vep{ $var->{ORIGINAL_ALT} };
        #if we've provided a list of RefSeq ranks skip if this gene symbol is not listed
        if (keys %transcript_ranks){
            my $sym = uc($csq->{symbol});
            next if not exists $transcript_ranks{$sym};
        }
        push @csq_to_rank , $csq;
    }
    my $csq_to_report = rankTranscriptsAndConsequences(\@csq_to_rank); 
    
    my @vep_fields = (qw / symbol feature consequence hgvsc hgvsp / );
    my $most_damaging_csq = getMostDamagingConsequence($csq_to_report);
    foreach my $f (@vep_fields){
        push @row, $csq_to_report->{$f};
    }
    #choose sheet to write to
    my $s_name;
    my %lofs = map {$_ => undef} qw /
        transcript_ablation  
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        frameshift_variant
        stop_lost
        start_lost
        transcript_amplification
    /;
    my @cadd = split(",", VcfReader::getVariantInfoField($l, "CaddPhredScore"));
    my $allele_cadd = $cadd[ ($var->{ALT_INDEX} -1) ];
    if (@$matches){
        $s_name = "HGMD";
    }elsif($cvarPathognic){
        $s_name = "ClinVarPathogenic";
    }elsif (exists $lofs{$most_damaging_csq}){
        $s_name = "LOF";
    }elsif ($most_damaging_csq eq 'missense_variant'){
        #check if in collagen domain
        if (exists $csq_to_report->{glyxy} ){
            $s_name = "CollagenGlyXY";
        #check if damaging or benign...
        }elsif ($allele_cadd  >= 10 and 
            $csq_to_report->{polyphen} =~ /damaging/i and 
            $csq_to_report->{sift} =~ /deleterious/i
        ){
            $s_name = "DamagingMissense";
        }else{
            $s_name = "BenignMissense";
        }
    }elsif( $most_damaging_csq =~  /^inframe_(inser|dele)tion$/ or
            $most_damaging_csq eq 'protein_altering_variant'){
        #check if in collagen domain
        if (exists $csq_to_report->{glyxy}){
            $s_name = "CollagenGlyXY";
        }elsif ($allele_cadd >= 10){
            $s_name = "DamagingMissense";
        }else{
            $s_name = "BenignMissense";
        }
    }else{
        $s_name = "Other";
    }

    push @row, $csq_to_report->{domains};
    push @row, $csq_to_report->{glyxy_pos};
    push @row, $csq_to_report->{domain_coords};

    foreach my $f (qw / lof lof_filter lof_info lof_flags polyphen sift/){
        push @row, $csq_to_report->{$f};
    }

    push @row, $allele_cadd;

    my $uid_base = sprintf
        (
            "%s:%s-%s/%s-", 
            $var->{CHROM},  
            $var->{ORIGINAL_POS},
            $var->{ORIGINAL_REF},
            $var->{ORIGINAL_ALT},
        ); #we add sample IDS to this string to create UIDs for sample variants

    #TO DO
    #search primers
    if (@primers){
        @primer_hits = searchPrimers($var);
    }
    #add sample info to new array of array refs to be written 
    # in split cells alongside global variant fields
    my %samp_gts = VcfReader::getSampleActualGenotypes
        (
              line => $l, 
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    my %samp_gqs = VcfReader::getSampleGenotypeField
        (
              line => $l, 
              field => "GQ",
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    my %samp_ads = VcfReader::getSampleGenotypeField
         (
              line => $l,
              field => "AD",
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    
    my $variant_has_valid_sample = 0;
    foreach my $s (
            sort {$samples_to_columns{$a} <=> $samples_to_columns{$b}} 
            keys %samples_to_columns
    ){
        my @alts = split(/[\/\|]/, $samp_gts{$s});
        if (grep { $_ eq $var->{ORIGINAL_ALT} } @alts ){ 
            my @ads = split(",", $samp_ads{$s}); 
            my $depth = sum(@ads);
            if ($opts{d}){
                next if $opts{d} > $depth;
            }
            if ($opts{g}){
                eval "$opts{g} <= $samp_gqs{$s}" or next;
            }
            my $ab = 0;
            if ( $depth > 0){
                $ab = $ads[$var->{ALT_INDEX}]/$depth;
            }
            if (@{$opts{b}}){
                next if $ab < $opts{b}->[0];
                if (scalar @{$opts{b}} > 1){
                    next if $ab > $opts{b}->[1];
                }
            }
            $variant_has_valid_sample++;
            my $uid = "$uid_base$s";
            my @sample_cells = ( $s, $samp_gts{$s}, $samp_ads{$s}, $ab, $samp_gqs{$s}, $uid );
            if (@primers){
                if (@primer_hits){
                    push @sample_cells, join("/", @primer_hits); 
                }else{
                    push @sample_cells, '-';
                }
            }
            if (%validated){
                #search validations
                if (exists $validated{$uid}){
                    push @sample_cells, $validated{$uid};
                }else{
                    push @sample_cells, 0;
                }
            }
            if ($opts{phenotype_file}){
                #search phenotypes
                foreach my $ph (@pheno_fields){
                    if (exists $phenotypes->{$s}){
                        push @sample_cells, $phenotypes->{$s}->{$ph};
                    }else{
                        push @sample_cells, '-';
                    }
                }
            }
            push @split_cells, \@sample_cells;
            my $var_class = $s_name; 
            if ($s_name eq 'HGMD'){
                if (grep {$_ eq 'DM'} @{$hgmd{variant_class}}){
                    $var_class = "HGMD_DM";
                }else{
                    $var_class = "HGMD_other";
                }
            }
            {
                no warnings 'uninitialized';
                push @{$sample_vars{$s}->{$var_class}}, [$allele_cadd, @sample_cells, $var_class, @row, ];
            }
        }
    }
    if (not $opts{n} or $variant_has_valid_sample){
        my $category = $s_name;
        $variant_counts{$category}++;
        if ($opts{1}){
            $s_name = 'Variants';
        }
        unshift @row, $category;
        unshift @row, "$category.$variant_counts{$category}";
        if ($opts{x}){ #do not merge
            foreach my $s (@split_cells){
                my @full_line = (@row, @$s); 
                $xl_obj->writeLine
                (
                    line => \@full_line,
                    worksheet  => $sheets{$s_name},
                );
            }
        }else{
            $xl_obj->writeLine
            (
                line       => \@row, 
                worksheet  => $sheets{$s_name},
                succeeding => \@split_cells,
            );
        }
    }
}
###########################################################
sub searchPrimers{
#returns array of matching primer IDs
    my $var = shift;
    my @hits = (); 
    my $u = $#primers;
    my $l = 0;
    my $var_end = $var->{POS} + length($var->{REF}); 
    while ( $l <= $u ){
        my $i = int ( ( $u + $l ) / 2 ) ;
        if ($var->{CHROM} lt $primers[$i]->{CHROM}){
            $u = $i - 1;
        }elsif ($var->{CHROM} gt $primers[$i]->{CHROM}){
            $l = $i + 1;
        }else{#same chrom
            my $target_start = $primers[$i]->{START} + length( $primers[$i]->{F} ) + 50;
            my $target_end = $primers[$i]->{END} - length( $primers[$i]->{R} ) - 50;
                #require primers to be at least 50 bp away from target nucleotide
            if ($var->{POS} < $target_start){
                $u = $i -1;
            }elsif($var_end > $target_end){
                $l = $i + 1;
            }else{
                push @hits, $primers[$i]->{ID};
                #search for hits either site;
                for (my $j = $i - 1; $j >= 0; $j--){
                    my $t_start = $primers[$j]->{START} + length( $primers[$j]->{F} ) + 50;
                    my $t_end = $primers[$j]->{END} - length( $primers[$j]->{R} ) - 50;
                    if ($var->{POS} >= $t_start && $var_end <= $t_end){
                        push @hits, $primers[$j]->{ID};
                    }else{
                        last if $var_end < $t_start;
                    }
                } 
                for (my $j = $i + 1; $j < @primers; $j++){
                    my $t_start = $primers[$j]->{START} + length( $primers[$j]->{F} ) + 50;
                    my $t_end = $primers[$j]->{END} - length( $primers[$j]->{R} ) - 50;
                    if ($var->{POS} >= $t_start && $var_end <= $t_end){
                        push @hits, $primers[$j]->{ID};
                    }else{
                        last if $var->{POS} > $t_end;
                    }
                }
                return @hits;
            }
        }
    }
    return @hits;
}

###########################################################
#kept for legacy in case need to switch back to ClinVar VCF
sub getClinVarCodeVcf{
    my ($clinvars, $var) = @_;
    #array ref of clinvar VCF lines and a minimized variant hash from input
    my @annots = ();
    foreach my $cl (@$clinvars){
        my @match = split("\t", $cl);
        my @alts = split(",", VcfReader::getVariantField(\@match, "ALT"));
        my @clnsigs = split(",", VcfReader::getVariantInfoField(\@match, "CLNSIG"));
        my @codes = ();
        if (@alts == @clnsigs){#CLNSIG codes are a bit sketchy, sometimes one per allele sometimes not
            if (my $i = alleleMatches($var, $cl)){
                push @codes, split(/\|/, $clnsigs[$i-1]);
            }
        }else{
            foreach my $c (@clnsigs){
                push @codes, split(/\|/, $c);
            }
        }
        foreach my $code (@codes){
            push @annots, "$code ($clnsig_codes{$code})";
        }
    }
    return join(", ", @annots); 
}

###########################################################
sub getClinSig{
    my ($clinvars, $var) = @_;
    #array ref of clinvar lines and a minimized variant hash from input
    my @clnsig = ();
    my @trait  = ();
    my $isPathogenic = 0;
    my $conflicted   = 0;
    foreach my $cl (@$clinvars){
        my @match = split("\t", $cl);
        push @clnsig, $match[$clinvar_sargs{columns}->{clinical_significance}] ;
        push @trait, $match[$clinvar_sargs{columns}->{all_traits}] ;
        $isPathogenic +=  $match[$clinvar_sargs{columns}->{pathogenic}] ;
        $conflicted   +=  $match[$clinvar_sargs{columns}->{conflicted}] ;
    }
    return join(", ", @clnsig), join(", ", @trait), $isPathogenic, $conflicted; 
}
    
###########################################################
#kept for legacy in case need to switch back to ClinVar VCF
sub searchClinVarVcf{
    return if not $opts{c};
    my $var = shift;#$var is a ref to a single entrey from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    #below is a hack because ClinVar file should be tsv.gz not VCF
    my @hits = VcfReader::searchByRegion(
        %clinvar_sargs,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatchesClinVar($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}

###########################################################

sub searchClinVar{
    return if not $opts{c};
    my $var = shift;#$var is a ref to a single entrey from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    #below is a hack because ClinVar file should be tsv.gz not VCF
    my @hits = VcfReader::searchByRegion(
        %clinvar_sargs,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatchesClinVar($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}

###########################################################
sub searchHgmd{
    my $var = shift;#$var is a ref to a single entry from minimized alleles hash
    #simplify alleles and check if there's a match in HGMD file
    my @matches = ();
    my @hits = VcfReader::searchByRegion(
        %search_args,
        chrom => $var->{CHROM},
        start => $var->{POS},
        end   => $var->{POS} + length($var->{REF}) - 1,
    );
    foreach my $h (@hits){
        if (alleleMatches($var, $h)){
            push @matches, $h;
        }
    }
    return @matches;
}


##########################################################
sub alleleMatches{
    my ($var, $line) = @_;
    #$var is a ref to a single entrey from minimized alleles hash
    #$line is a VCF entry
    #returns allele code for matching alt
    my @split = split ("\t", $line);
    my $chrom = VcfReader::getVariantField(\@split, 'CHROM');
    next if $chrom ne $var->{CHROM}; 
    my %l_min= VcfReader::minimizeAlleles(\@split);
    foreach my $k (keys %l_min){
        next if $l_min{$k}->{POS} != $var->{POS};
        next if $l_min{$k}->{REF} ne $var->{REF};
        next if $l_min{$k}->{ALT} ne $var->{ALT};
        return $k;
    }
    return 0;
}
##########################################################
sub alleleMatchesClinVar{
    my ($var, $line) = @_;
    #$var is a ref to a single entrey from minimized alleles hash
    #$line is a clinvar entry
    #returns 1 if it matches
    my @split = split ("\t", $line);
    my $chrom = $split[$clinvar_sargs{columns}->{chrom}] ;
    next if $chrom ne $var->{CHROM}; 
    my $pos = $split[$clinvar_sargs{columns}->{pos}] ;
    #alleles should already be minimized
    my $ref = $split[$clinvar_sargs{columns}->{ref}] ;
    my $alt = $split[$clinvar_sargs{columns}->{alt}] ;
    return 0 if $pos != $var->{POS};
    return 0 if $ref ne $var->{REF};
    return 0 if $alt ne $var->{ALT};
    return 1;
}
###########################################################
sub setupOutput{
    my @suffixes = (".vcf", ".txt");
    if (not $opts{o}){
        my ($out, $dir, $extension) = fileparse($opts{i}, @suffixes);
        $opts{o} = "$dir/$out.xlsx";
    }else{
        $opts{o} .= ".xlsx" if $opts{o} !~ /\.xlsx$/;
    }
#    $url_format = $xl_obj->createFormat(
#        color     => 'blue',
#        underline => 1,
#    );  
    if ($opts{1}){
        $xl_obj = TextToExcel->new( file=> $opts{o}, name => "Variants");
        $sheets{Variants} = 0;
        $header_formatting = $xl_obj->createFormat(bold => 1);
        writeHeader($sheets{Variants}, $headers{Variants});
    }else{
        $xl_obj = TextToExcel->new( file=> $opts{o}, name => "HGMD");
        $header_formatting = $xl_obj->createFormat(bold => 1);
        writeHeader($sheets{HGMD}, $headers{Variants});
        $sheets{HGMD} = 0;
        $sheets{ClinVarPathogenic} = addSheet("ClinVarPathogenic", $headers{Variants});
        $sheets{LOF} = addSheet("LOF", $headers{Variants});
        $sheets{CollagenGlyXY} = addSheet("CollagenGlyXY", $headers{Variants});
        $sheets{DamagingMissense} = addSheet("DamagingMissense", $headers{Variants});
        $sheets{BenignMissense} = addSheet("BenignMissense", $headers{Variants});
        $sheets{Other} = addSheet("Other", $headers{Variants});
    }

    $sheets{SampleSummary} = addSheet("Sample Summary", $headers{sample});
    $sheets{MostDamaging} = addSheet("MostDamaging", $headers{mostdamaging});
        
    if ($opts{f}){
        open ($FILTER_OUT, ">$opts{f}") or die "Could not create filter output file \"$opts{f}\": $!\n";
        print $FILTER_OUT join("\n", @vhead) . "\n";
    }
}

###########################################################
sub getHeaders{
    my %h = ();

    @{$h{Variants}} =  ( #header for all variant sheets
        qw/
            index
            Category
            Hgmd_ID
            Disease
            variant_class
            HGMD_Symbol
            HGVS
            ClinVarSig
            ClinVarTrait
            ClinVarConflicted
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            DomainPosition
            DomainCoordinates
            LoF
            LoF_Filter
            LoF_info
            LoF_flags
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
            UID
        /
    );
=cut
    @{$h{HGMD}} =  ( 
        qw / 
            index
            Hgmd_ID
            Disease
            variant_class
            HGMD_Symbol
            HGVS
            ClinVarSig
            ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
            UID
     /);
    @{$h{LOF}} =  ( 
        qw / 
            index
            ClinVarSig
            ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            LoF
            LoF_Filter
            LoF_info
            LoF_flags
            CADD
            Sample
            GT
            AD
            AB
            GQ
            UID
     /);
    
    @{$h{clinvar}} =  ( 
        qw / 
            index
            ClinVarSig
            ClinVarTrait
            Conflicted
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
            UID
     /);

    @{$h{glyxy}} =  ( 
        qw / 
            index
            ClinVarSig
            ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            DomainPosition
            DomainCoordinates
            Polyphen
            SIFT
            CADD
            Sample
            GT
            AD
            AB
            GQ
            UID
     /);




    @{$h{missense}} =  ( 
        qw / 
            index
            ClinVarSig
            ClinVarTrait
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature 
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            Polyphen
            SIFT
            CADD    
            Sample
            GT
            AD
            AB
            GQ
            UID
     /);
=cut

    @{$h{sample}} =  (  
        qw / 
            Sample
            HGMD_DM
            ClinVarPathogenic
            CollagenGlyXY
            HGMD_other
            LOF
            DamagingMissense
            BenignMissense
            Other
            /
    );
    @{$h{mostdamaging}} =  ( #header for all variant sheets
        qw/
            Sample
            GT
            AD
            AB
            GQ
            UID
            Category
            Hgmd_ID
            Disease
            variant_class
            HGMD_Symbol
            HGVS
            ClinVarSig
            ClinVarTrait
            ClinVarConflicted
            Chrom
            Pos
            Ref
            Alt
            ID
            Qual
            Filter
            Symbol
            Feature
            Consequence 
            HGVSc 
            HGVSp 
            Domains
            DomainPosition
            DomainCoordinates
            LoF
            LoF_Filter
            LoF_info
            LoF_flags
            Polyphen
            SIFT
            CADD
        /
    );

    my $md_spl = 6;#colu no. to add primers or validated columns to mostdamaging sheet
    my $v_spl = 9; #column no. to add uniprot/cdd domain info to variant sheets
    if ($opts{primer_file}){
        push @{$h{Variants}}, "primers";
        splice (@{$h{mostdamaging}}, $md_spl++, 0, "primers");
    }
    if ($opts{validations_file}){
        push @{$h{Variants}}, "validated?(1=yes,2=no,0=not_done)";
        splice (@{$h{mostdamaging}}, $md_spl++, 0, "validated?(1=yes,2=no,0=not_done)");
    }
    if ($opts{phenotype_file}){
        foreach my $ph (@pheno_fields){
            push @{$h{Variants}}, $ph;
            splice (@{$h{mostdamaging}}, $md_spl++, 0, $ph); 
        }
    }
    return %h;
}
###########################################################
sub writeHeader{
    my $sheet = shift; 
    my $header = shift;
    $xl_obj->writeLine(line => $header, format => $header_formatting, worksheet => $sheet);
}

###########################################################
sub addSheet{
    my ($name, $header) = @_; 
    my $sheet = $xl_obj->addWorksheet($name);
    writeHeader($sheet, $header);
    return $sheet;
}

###########################################################
sub rankTranscriptsAndConsequences{
    my $csq_array = shift;#ref to array of VEP consequences 
    @$csq_array = rankConsequences(@$csq_array); 
    my $most_damaging = $csq_array->[0]->{consequence} ;
    @$csq_array = rankTranscripts(@$csq_array); 
    if ($most_damaging eq 'missense_variant' or 
        $most_damaging eq 'protein_altering_variant' or  
        $most_damaging =~ /^inframe_(inser|dele)tion$/
    ){
        #check if is mutation of glycine in collagen triple helix and return
        if (my $glyxy_csq = checkCollagenDomain($csq_array)){
            return $glyxy_csq;
        }
    }
    return first { $_->{consequence} eq $most_damaging } @$csq_array;
}
###########################################################
sub checkCollagenDomain{
# only applicable to missense and inframe indels
# returns first VEP consequence with a missense altering a Gly 
# in a collagen triple helix or an inframe indel that is alters
# a number of amino acids not divisible by 3
    my $csq_array = shift;#ref to array of VEP consequences 
    foreach my $c (@$csq_array){
        my $dom_hash;
        my $gxy_frame = 0;
        my $gxy_length;
        if (not $opts{s}){#if not using the --scan_gxy option
            next if not $c->{domains};
            my @domains = split ("&", $c->{domains}); #domains variant overlaps
            if (grep { /Pfam_domain:PF01391/i } @domains){
                if ($opts{v}){
                    informUser("Found variant in collagen triple helix\n");
                    informUser("Checking domain sequence using Ensembl REST...\n");
                }
        #get coordinates of domains to find Gly-X-Y pattern using ensembl REST API
                #my $seq_hash  = ensRestQuery("$server/sequence/id/$c->{ensp}?"); 
                my $rest_url = "$server/overlap/translation/$c->{ensp}?";
                my $pfam_ar = ensRestQuery($rest_url);
                if (ref $pfam_ar ne 'ARRAY'){
                    die "Required array reference from Ensembl REST query: $rest_url\n" 
                }
                $dom_hash = findOverlappingDomain($pfam_ar, "PF01391", $c);
                
                # CHECK SEQ TO MAKE FIND 'FRAME' of GLY-X-Y repeat 
                my $seq_hash = getProteinSequence($c->{ensp}); 
                $dom_hash->{seq} = substr(
                                    $seq_hash->{seq}, 
                                    $dom_hash->{start} -1,  
                                    $dom_hash->{end} - $dom_hash->{start} + 1, 
                ); 
                #find the first repeat of 5 or more Gly-X-Ys 
                # -- TO DO - more sophisticated HMM perhaps?
                if ($opts{v}){
                    informUser("domain $dom_hash->{id}: $dom_hash->{start}-$dom_hash->{end} - $dom_hash->{seq}\n");
                    informUser("frame: $gxy_frame\n");
                }
            }
        }else{
            #scan sequence ourselves for GXY repeats
            next if not $c->{ensp};#skip non-coding
            next if $c->{symbol} !~ /^COL/;#skip non-collagens
            next if not $c->{amino_acids};#skip any non-coding variant
            my $seq_hash = getProteinSequence($c->{ensp}); 
            $dom_hash = findGlyXY($seq_hash->{seq}, $c); 
            if ($dom_hash and $opts{v}){
                informUser("Found putative Gly-X-Y domain in $c->{symbol} ($c->{ensp}) $dom_hash->{start}-$dom_hash->{end}:\n$dom_hash->{seq}\n");
            }
        }#limitation at the moment - does not allow for deletion spanning more than one GlyXY domains
        if ($dom_hash){
            $gxy_length =  1 + $dom_hash->{end} - $dom_hash->{start}; 
            if ($dom_hash->{seq} =~ /(G[A-Z]{2}){5,}/){
                $gxy_frame = $-[0] % 3;
                $gxy_length -= $gxy_frame;
            }else{
                die "Could not find Gly-X-Y repeat in PF01391 domain seq:\n$dom_hash->{seq}\n\n";
            }
            my @aa = split ("/", $c->{amino_acids} ) ;
            s/-// for @aa; #turn deleted AAs to empty strings
            if (length($aa[0]) == length($aa[1])){
                #check whether an essential glycine has been altered
                if (essentialGlyAltered(\@aa, $c, $dom_hash->{start} + $gxy_frame, $dom_hash->{end})){
                    my $gxy_pos = 1 + $c->{protein_position} - $dom_hash->{start}; 
                    $c->{glyxy_pos} = "$gxy_pos/$gxy_length";
                    $c->{glyxy} = 1;
                    $c->{domain_coords} = "$dom_hash->{id}: $dom_hash->{start}-$dom_hash->{end}"; 
                    return $c;
                }
            }else{
                my ($p_start, $p_end) = split("-", $c->{protein_position}); 
                my ($ref, $alt);
                ($ref, $alt, $p_start) = simplifyIndel($aa[0], $aa[1], $p_start);
                $p_end = $p_start + length($ref) - 1;
                $p_end = $p_end >= $p_start ? $p_end : $p_start;#fix in case length of $ref is 0
                my $diff = length($aa[0]) - length($aa[1]); 
                if ($diff > 0){#deletion
                    #check which portion of an indel lies within domain in case only partial overlap
                    my $dp_start = $p_start >= $dom_hash->{start} ? $p_start : $dom_hash->{start};
                    my $dp_end   = $p_end   <= $dom_hash->{end}   ? $p_end   : $dom_hash->{end};
                    if ($dp_start != $p_start){
                        my $trim = $dp_start - $p_start;
                        substr($aa[0], 0, $trim, ""); #remove beginning portion that does not lie in domain
                        substr($aa[1], 0, $trim, ""); #remove beginning portion that does not lie in domain
                    }
                    if ($dp_end != $p_end){
                        my $trim = $dp_end - $p_end;
                        $aa[0] = substr($aa[0], 0, $trim ); #remove end portion that does not lie in domain
                        $aa[1] = substr($aa[1], 0, $trim ); #remove end portion that does not lie in domain
                    }
                    $diff = length($aa[0]) - length($aa[1]); 
                }
                #check if insertion/deletion is divisible by 3 
                #if not divisible by 3 then triple helix broken
                if ($diff % 3){
                    $c->{glyxy} = 1;
                    $c->{domain_coords} = "$dom_hash->{id}: $dom_hash->{start}-$dom_hash->{end}"; 
                    return $c;
                }
                # if divisible by 3 check if a pos that should be a 
                # glycine has been altered
                if (essentialGlyAltered(\@aa, $c, $dom_hash->{start} + $gxy_frame, $dom_hash->{end})){
                    $c->{glyxy} = 1;
                    my $start_distance = $p_start > $dom_hash->{start} ? 1 + $p_start - $dom_hash->{start} : 1;
                    my $end_distance   = $p_end   < $dom_hash->{end}   ? 1 +  $p_end  - $dom_hash->{start} : "";
                    $c->{glyxy_pos} = "$start_distance-$end_distance/$gxy_length";
                    $c->{domain_coords} = "$dom_hash->{id}: $dom_hash->{start}-$dom_hash->{end}"; 
                    return $c;
                }
            }
        }
    }
}

###########################################################
sub findGlyXY{
    my ($seq, $csq) = @_; 
    my @gxys = scanForGlyXY($seq);
    my ($p_start, $p_end) = split("-", $csq->{protein_position}); 
    $p_end ||= $p_start; 
    foreach my $hash (@gxys){
        if ($hash->{start} <= $p_start and $hash->{end} >= $p_start){
            #variant lies within this domain
            return $hash;
        }elsif ($hash->{start} <= $p_end and $hash->{end} >= $p_end){
            return $hash;
        }
    }
    return undef;
}

###########################################################
sub scanForGlyXY{
    my $seq = shift;
    my @glyxys = ();
    my $n = 0;
    while ($seq =~ /((G[A-Z]{2}){5,})/g){
        my %dom = (); 
        $n++;
        $dom{start} = $-[0] + 1;
        $dom{end} = $-[0] + length($1);
        $dom{id} = "glyxy-$n";
        push @glyxys, \%dom;
    }
    #merge hashes within 6 aa of each other to allow for
    #occasional non-gly-x-y triplets
    my $prev_hash = $glyxys[0];
    my @merged = ();
    for (my $i = 1; $i < @glyxys; $i++){
        my $gap = $glyxys[$i]->{start} - $prev_hash->{end};
        if ($gap == 7 or $gap == 4){
            $prev_hash->{end} = $glyxys[$i]->{end};
        }else{
            $prev_hash->{seq} = substr(
                              $seq,
                              $prev_hash->{start} -1,  
                              $prev_hash->{end} - $prev_hash->{start} + 1, 
            );
            push @merged, $prev_hash;
            $prev_hash = $glyxys[$i];
        }
    }
    if ($prev_hash){
        $prev_hash->{seq} = substr(
                          $seq,
                          $prev_hash->{start} -1,  
                          $prev_hash->{end} - $prev_hash->{start} + 1, 
        );
        push @merged, $prev_hash ;
    }
    return @merged;
}
        


###########################################################
sub essentialGlyAltered{
    #return 1 if Gly at pos 0, 3, 6, etc. of domain altered
    my ($aa, $c, $dom_start, $dom_end) = @_;
    # ref to amino acid array (ref and alt), VEP csq hash ref 
    # and domain details hash ref
    my ($p_start, $p_end) = split("-", $c->{protein_position}); 
    my @gly_idxs = grep { substr($aa->[0], $_, 1) eq "G" } 0..(length($aa->[0]));
    foreach my $i (@gly_idxs){
        my $pos = $i + $p_start ;#actual protein pos of this Gly
        next if $pos < $dom_start;
        last if $pos > $dom_end;
        my $dist = $pos - $dom_start;#distance of Gly from domain start
        next if $dist % 3; #not an essential Gly if not at pos 0 of 3 aa repeat
        if ($i > length($aa->[1])){#Essential Gly deleted(?)
            return 1;
        }
        if (substr($aa->[1], $i, 1) ne 'G'){ # Essential Gly altered
            return 1;
        }
    }
    #if insertion, check we have Glys in the right positions
    if (length($aa->[0]) < length($aa->[1])){
        my $frame = ($p_start - $dom_start) % 3 ;
        for (my $i = 0; $i < length($aa->[1]); $i += 3){
            next if $i < length($aa->[0]); #only check the remaining portion of seq not in the Ref AAs
            my $res = substr($aa->[1], $i + $frame, 1);
            if ($res ne 'G'){
                return 1;
            }
        }
    }
    return 0;
}

###########################################################
sub simplifyIndel{
    my ($ref, $alt, $p_start) = @_;
    while (length($ref) > 1 and length($alt) > 1){
        while ( (substr($ref, -1, 1)) eq (substr($alt, -1, 1)) 
                and (length($ref) > 0) and (length($alt) > 0)
        ){
            substr($ref, -1, 1, "");
            substr($alt, -1, 1, "");
        }
        while ( (substr($ref, 0, 1)) eq (substr($alt, 0, 1)) 
                and (length($ref) > 0) and (length($alt) > 0)
        ){
            substr($ref, 0, 1, "");
            substr($alt, 0, 1, "");
            $p_start++;
        }
    }
    return ($ref, $alt, $p_start);
}
###########################################################
sub findOverlappingDomain{
    my ($pfam_ar, $id, $csq) = @_;
    my ($p_start, $p_end) = split("-", $csq->{protein_position}); 
    $p_end ||= $p_start; 
    foreach my $hash (@$pfam_ar){
        next if $hash->{id} ne $id;
        next if $hash->{seq_region_name} ne $csq->{ensp};
        if ($hash->{start} <= $p_start and $hash->{end} >= $p_start){
            #variant lies within this domain
            return $hash;
        }elsif ($hash->{start} <= $p_end and $hash->{end} >= $p_end){
            #variant lies within this domain
            return $hash;
        }
    }
    die "No overlap found with $id domain for $csq->{hgvsc}/$csq->{hgvsp} in REST database\n";
}
###########################################################
sub getProteinSequence{
    my $ensp = shift;
    if (not exists $seq_cache{$ensp}){
        my $rest_url = "$server/sequence/id/$ensp";
        $seq_cache{$ensp} = ensRestQuery($rest_url);
        die "No sequence returned from REST query: $rest_url\n" if not $seq_cache{$ensp}->{seq};
    }
    return $seq_cache{$ensp};

}
###########################################################
sub ensRestQuery{
    my $url = shift;
    my $response = $http->get($url, {
          headers => { 'Content-type' => 'application/json' }
    });
    die "Ensembl REST query ('$url') failed!\n" unless $response->{success};
    if(length $response->{content}) {
      return decode_json($response->{content});
    }
    die "No content for Ensembl REST query ('$url')!\n";
}
        
    
###########################################################
sub rankTranscripts{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { getTranscriptsRanks( uc($a->{feature}) ) <=> getTranscriptsRanks( uc($b->{feature}) ) } @vars;
}
###########################################################
sub getTranscriptsRanks{
    my $symbol = shift;
    my $transcript = shift; 
    return -1 if not exists $transcript_ranks{$symbol};
    my $idx = first { $transcript_ranks{$symbol}->[$_] eq $transcript } 0 ..  $#{$transcript_ranks{$symbol}}; #also returns undef if transcript not found
    return -1 if not defined $idx;
    return $idx;
}
###########################################################
sub setTranscriptsRanks{
    return if ( not $opts{t} );
    open (my $TR, $opts{t}) or die "Can't open --transcripts file ($opts{t}) for reading: $!\n";
    my $header = <$TR>;
    $header =~ s/^#+//;
    my %tr_columns = getColumns($header);
    foreach my $req ( qw / symbol transcript uniprot / ){
        if (not exists  $tr_columns{$req}){
            die "Could not find required column '$req' ".
                "in header of --transcripts file $opts{t}. ".
                "Found the following columns:\n" . 
                join("\n", sort { $tr_columns{$a} <=> $tr_columns{$b} } keys %tr_columns) . 
                "\n";
        }
    }

    while (my $line = <$TR>){
        chomp $line;
        my @split = split("\t", $line); 
        my $symbol = uc($split[ $tr_columns{symbol} ]);
        my $transcript = uc ($split[ $tr_columns{transcript} ]); 
        push @{$transcript_ranks{$symbol} }, $transcript; 
        my $uniprot = uc($split[ $tr_columns{uniprot} ]);
        if ($uniprot ne '.'){
             $enst_to_uniprot{$transcript} = $uniprot;
        }
    }
    close $TR;
}

###########################################################
sub readUniprotFile{
    return if ( not $opts{u} );
    open (my $UNI, $opts{u}) or die "Can't open --uniprot file ($opts{u}) for reading: $!\n";
    my $header = <$UNI>;
    $header =~ s/^#+//;
    my %tr_columns = getColumns($header);
    foreach my $req ( qw / UniprotId Start End Feature Note / ){
        if (not exists  $tr_columns{$req}){
            ##TODO
        }
    }
    ##TODO
    close $UNI;
}

###########################################################
sub readCddFile{
    ...
}

###########################################################

sub setMutationRules{
    ...
}

###########################################################
sub getColumns{
    #returns a hash of column names to 0-based column number
    my ($header, $to_lower, $delimiter) = @_;
    $delimiter ||= "\t"; 
    chomp (my @head = split($delimiter, $header)); 
    my $n = 0;
    my %columns = (); 
    if ($to_lower){
        %columns = map { lc($_) => $n++ } @head;
    }else{
        %columns = map { $_ => $n++ } @head;
    }
    return %columns;
}


###########################################################
sub rankConsequences{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { 
        $so_ranks{getMostDamagingConsequence($a)} <=> 
        $so_ranks{getMostDamagingConsequence($b)} 
    } @vars;
    
}

###########################################################
sub getMostDamagingConsequence{
    my $csq = shift;#hash ref to VEP consequences for single transcript/feature
    my @s_csq = split("&", $csq->{consequence} );
    @s_csq = sort { $so_ranks{$a} <=> $so_ranks{$b} } @s_csq;
    return $s_csq[0];
}
    
###########################################################
sub setConsequenceRanks{
    my @so_terms = qw /
        transcript_ablation  
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        frameshift_variant
        stop_lost
        start_lost
        transcript_amplification
        inframe_insertion
        inframe_deletion
        missense_variant
        protein_altering_variant
        splice_region_variant
        incomplete_terminal_codon_variant
        stop_retained_variant
        synonymous_variant
        coding_sequence_variant
        mature_miRNA_variant
        5_prime_UTR_variant
        3_prime_UTR_variant
        non_coding_transcript_exon_variant
        intron_variant
        NMD_transcript_variant
        non_coding_transcript_variant
        upstream_gene_variant
        downstream_gene_variant
        TFBS_ablation
        TFBS_amplification
        TF_binding_site_variant
        regulatory_region_ablation
        regulatory_region_amplification
        feature_elongation
        regulatory_region_variant
        feature_truncation
        intergenic_variant 
    /;
    my $n = 0;
    %so_ranks = map { $_ => $n++ } @so_terms; 
}

###########################################################
sub getClnSigCodes{
    return (
        0 => "Uncertain significance",
        1 => "not provided",
        2 => "Benign",
        3 => "Likely benign",
        4 => "Likely pathogenic",
        5 => "Pathogenic",
        6 => "drug response",
        7 => "histocompatibility",
        255 => "other"
    );
}

#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    if ($progressbar){
        $progressbar->message( "[INFO - $time] $msg" );
    }else{
        print STDERR "[INFO - $time] $msg";
    }
}
#################################################

=head1 NAME

clinSummarizerTaad.pl - assess known/potentially pathogenic mutations for TAAD

=head1 SYNOPSIS

        clinSummarizerTaad.pl -i [vcf file] -m [HGMD mart converted vcf file] [options]
        clinSummarizerTaad.pl -h (display help message)
        clinSummarizerTaad.pl --manual (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file (prefiltered on public databases for allele frequency).

=item B<-o    --output>

Output xlsx file. Defaults to input file with .xlsx extension added.

=item B<-n    --no_blanks>

Do not output variants which do not have at least one sample passing the variant criteria.

=item B<-m    --hgmd>

VCF of HGMD variations converted with hgmdMartToVcf.pl.

=item B<-c    --clinvar>

Optional ClinVar VCF to add ClinVar CLINSIG annotations.

=item B<-t    --transcripts>

Optional tsv list of transcript IDs in order of priority.

=item B<-b    --allele_balance>

Min and optional max alt allele ratio per sample call.

=item B<-d    --depth>

Optional min depth for sample call.

=item B<-a    --allele_cutoff>

Remove if allele is present in this proportion or more calls.

=item B<-f    --filter_output>

Optional output file for calls filtered on allele_cutoff.

=item B<-r    --rest_server>

URL of REST server to use if not the default (http://grch37.rest.ensembl.org).

=item B<-g    --gq>

Min GQ quality for calls.

=item B<-p    --progress>

Display a progress bar.

=item B<-v    --verbose>

Show verbose progress information.

=item B<-h    --help>

Display help message.

=item B<--manual>

Show manual page.

=back 

=cut

=head1 DESCRIPTION

Creates an xlsx spreadsheet with worksheets for TAAD mutations based on categorisations such as known HGMD mutations, ClinVar pathogenic/likely pathogenic mutations, loss of function mutations, mutations of Glycines in Collagen triple helices etc.

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
