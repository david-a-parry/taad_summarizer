#!/usr/bin/perl
#David Parry August 2011

#TO DO - implement an allele frequency filter which uses information from ped files
# to only count per family

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Term::ProgressBar;
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin";
use ParseVCF;
my $vcf;
my $out;
my @samples;
my $check_presence_only;
my $ignore_non_existing;#don't exit if a sample is not found
my @reject = ();#reject if allele is present in these samples
my @reject_except = (); #reject all except these samples
my $af = 0; #filter if allele frequency equal to or greater than this
my $threshold;
my $quality = 20;
my $allele_depth_cutoff = 0;#even if genotype isn't called use this value to filter on reported allele depth
my $allele_ratio_cutoff = 0;#even if genotype isn't called use this value to filter on relative reported allele depth
my $aff_genotype_quality ;#will convert to $genotype_quality value if not specified
my $unaff_genotype_quality ;#will convert to $genotype_quality value if not specified
my $confirm_missing; #only consider variants where there is sufficient genotype information from all --reject samples to exclude
my $num_matching;
my $help;
my $manual;
my $progress;
my %opts = ('existing' => \$ignore_non_existing,
        'input' =>\$vcf,
        "output" => \$out,
        'samples' => \@samples,
        'reject' => \@reject,
        'reject_all_except' => \@reject_except,
        'frequency' => \$af,
        'threshold' => \$threshold,
        'presence' => \$check_presence_only,
        'confirm_missing' => \$confirm_missing,
        'quality' => \$quality,
        'aff_quality' => \$aff_genotype_quality,
        'un_quality' => \$unaff_genotype_quality,
        'depth_allele_cutoff' => \$allele_depth_cutoff,
        'allele_ratio_cutoff' => \$allele_ratio_cutoff,
        'num_matching' => \$num_matching,
        'help' => \$help,
        "manual" => \$manual,
        'progress' => \$progress);

GetOptions(\%opts,
        'existing' => \$ignore_non_existing,
        'input=s' =>\$vcf,
        "output=s" => \$out,
        'samples=s{,}' => \@samples,
        'r|reject=s{,}' => \@reject,
        'x|reject_all_except:s{,}' => \@reject_except,
        'frequency=f' => \$af,
        'threshold=i' => \$threshold,
        'p|presence' => \$check_presence_only,
        'confirm_missing' => \$confirm_missing,
        'quality=i' => \$quality,
        'a|aff_quality=i',
        'un_quality=i',
        'depth_allele_cutoff=f' => \$allele_depth_cutoff,
        'z|allele_ratio_cutoff=f' => \$allele_ratio_cutoff,
        'num_matching=i' => \$num_matching,
        'help' => \$help,
        "manual" => \$manual,
        'b|progress' => \$progress) or pod2usage(-message=> "syntax error.\n");
pod2usage(-verbose => 2) if $manual;
pod2usage(-verbose => 1) if $help;
pod2usage(-message=> "syntax error: --input (-i) argument is required.\n") if not $vcf;
pod2usage(-message=> "syntax error: you must specify samples using at least one of the arguments --samples (-s), --reject (-r) or --reject_all_except (-x).\n") if not @samples and not @reject and not @reject_except;
pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) if ($quality < 0 );
pod2usage(-message => "--depth_allele_cutoff must be a value between 0.00 and 1.00.\n", -exitval => 2) if ($allele_depth_cutoff < 0 or $allele_depth_cutoff > 1.0 );
pod2usage(-message => "--allele_ratio_cutoff must be greater than 0.00.\n", -exitval => 2) if ($allele_ratio_cutoff < 0 );
pod2usage(-message => "--frequency must be a value between 0.00 and 1.00.\n", -exitval => 2) if ($af < 0 or $af > 1.0 );
if (defined $aff_genotype_quality){
    pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) 
        if $aff_genotype_quality < 0;
}else{
    $aff_genotype_quality = $quality;
}
if (defined $unaff_genotype_quality){
    pod2usage(-message => "Genotype quality scores must be 0 or greater.\n", -exitval => 2) 
        if $unaff_genotype_quality < 0;
}else{
    $unaff_genotype_quality = $quality;
}
print STDERR "Warning - --num_matching has no effect when --presence flag is set.\n" if $check_presence_only and $num_matching;

print STDERR "Initializing VCF input ($vcf)\n";
my $vcf_obj = ParseVCF->new(file=> $vcf);
my @not_found = ();
my @samples_checked = ();
my @reject_checked = ();
   
my ($samples_found, $samples_not_found) = check_samples(\@samples);
my ($reject_found, $reject_not_found) = check_samples(\@reject);
if (@$samples_not_found or @$reject_not_found){
    my @not_found = (@$samples_not_found, @$reject_not_found);
    if (not $ignore_non_existing){
        print STDERR "Warning - could not find the following samples in VCF:\n".join("\n", @not_found)."\n";
        if (@samples and not @$samples_found){
            print STDERR "Warning - no samples specified by --samples identified in VCF.\n";
        }
        if (@reject and not @$reject_found){
            print STDERR "Warning - no samples specified by --reject identified in VCF.\n";
        }
        @samples = @$samples_found;
        @reject = @$reject_found;
    }else{
        die "Could not find the following samples in VCF:\n".join("\n", @not_found)."\n";
    }
}
if (@reject_except){
    my @all = $vcf_obj->getSampleNames();
    push @reject_except, @samples; 
    my %subtract = map {$_ => undef} @reject_except;
    @all = grep {!exists $subtract{$_} } @all;
    push @reject, @all;
    my %seen = ();
    @reject = grep { ! $seen{$_}++} @reject;
}

if (not @reject and not @samples){
    print STDERR "Warning - no samples from --samples (-s), --reject (-r) or --reject_all_except (-x) argument found to filter. Your output will remain unchanged.\n";
}

my $OUT;
if ($out){
    open ($OUT, ">$out") || die "Can't open $out for writing: $!\n";
}else{
    $OUT = \*STDOUT;
}
my $progressbar;
if ($progress){
    if ($vcf eq "-"){
        print STDERR "Can't use --progress option when input is from STDIN\n";
        $progress = 0;
    }else{
        $progressbar = Term::ProgressBar->new({name => "Filtering", count => $vcf_obj->countLines("variants"), ETA => "linear", });
    }
}
my $next_update = 0;
print $OUT  $vcf_obj->getHeader(0) ."##filterOnSample.pl\"";
my @opt_string = ();
foreach my $k (sort keys %opts){
    if (not ref $opts{$k}){
        push @opt_string, "$k=$opts{$k}";
    }elsif (ref $opts{$k} eq 'SCALAR'){
        if (defined ${$opts{$k}}){
            push @opt_string, "$k=${$opts{$k}}";
        }else{
            push @opt_string, "$k=undef";
        }
    }elsif (ref $opts{$k} eq 'ARRAY'){
        if (@{$opts{$k}}){
            push @opt_string, "$k=" .join(",", @{$opts{$k}});
        }else{
            push @opt_string, "$k=undef";
        }
    }
}
print $OUT join(" ", @opt_string) . "\"\n" .  $vcf_obj->getHeader(1);
my $line_count = 0;
LINE: while (my $line = $vcf_obj->readLine){
    $line_count++;
    if ($progress){
        $next_update = $progressbar->update($line_count) if $line_count >= $next_update;
    }
    my %alleles = ();
    my %min_allele_ratios = ();#collect the minimum allele depth ratio in called samples for comparison with --reject samples
    #do samples first for efficiency (if they don't have a variant allele)
    if (@samples){
SAMPLE: foreach my $sample (@samples){
            my $call = $vcf_obj->getSampleCall(sample => $sample, minGQ => $aff_genotype_quality);
            if ($call =~ /(\d+)[\/\|](\d+)/){
                if ($1 == 0 and $2 == 0){
                    if ($check_presence_only or $num_matching){
                        next SAMPLE;
                    }else{
                        next LINE;#by default only keep variants present in all samples
                    }
                }else{
                    #samples only get added to %alleles hash if they are called and not reference
                    #because we compare the samples in %alleles hash only this allows variants 
                    #with no call to go through
                    push (@{$alleles{$sample}}, $1, $2);
                    
                    if ($allele_ratio_cutoff){#find the min ad ratios for --samples
                        my $ad = $vcf_obj->getSampleGenotypeField(sample => $sample, field => 'AD');
                        if ($ad){ 
                            my @ads = split(",", $ad);
                            @ads = grep {! /\./ } @ads; #sometimes with no reads an AD of '.' is given
                            my $dp = sum(@ads);
                            if ($dp){
                                for (my $i = 0; $i < @ads; $i++){
                                    my $ratio = $ads[$i]/$dp;
                                    if (exists $min_allele_ratios{$i}){
                                        $min_allele_ratios{$i} = 
                                            ($ratio >= $min_allele_ratios{$i} ? $ratio : $min_allele_ratios{$i});
                                    }else{
                                        $min_allele_ratios{$i} = $ratio;
                                    }
                                }
                            }
                        }
                    }
                }
            }else{
                next LINE unless $check_presence_only or $num_matching;
            }
        }
        next LINE if (keys%alleles < 1);#i.e. if only reference (0) or no calls were found
    }#otherwise we'll collect --reject alleles and see if there are any alts that aren't in %reject_alleles

    my %reject_alleles;
    my $total_reject = 0;
    if (@reject){
        foreach my $reject (@reject){
            my $call = $vcf_obj->getSampleCall(sample => $reject, minGQ => $unaff_genotype_quality);
            if ($call =~ /(\d+)[\/\|](\d+)/){
                $reject_alleles{$1}++; #store alleles from rejection samples as keys of %reject_alleles
                $reject_alleles{$2}++;
                $total_reject += 2;
            }elsif($confirm_missing){
                next LINE;
            }
            if($allele_depth_cutoff){
                #reject even if uncalled if proportion of alt allele >= $allele_depth_cutoff
                my $ad = $vcf_obj->getSampleGenotypeField(sample => $reject, field => 'AD');
                if ($ad){ 
                    my @ads = split(",", $ad);
                    @ads = grep {! /\./ } @ads; #sometimes with no reads an AD of '.' is given
                    my $dp = sum(@ads);
                    if ($dp){
                        for (my $i = 0; $i < @ads; $i++){
                            $reject_alleles{$i}++ if $ads[$i]/$dp >= $allele_depth_cutoff;
                        }
                    }
                }
            } 
            if ($allele_ratio_cutoff){#find the max ad ratios for --reject
                my $ad = $vcf_obj->getSampleGenotypeField(sample => $reject, field => 'AD');
                if ($ad){ 
                    my @ads = split(",", $ad);
                    @ads = grep {! /\./ } @ads; #sometimes with no reads an AD of '.' is given
                    my $dp = sum(@ads);
                    if ($dp){
                        for (my $i = 0; $i < @ads; $i++){
                            my $ratio = $ads[$i]/$dp;
                            if (exists $min_allele_ratios{$i}){
                                #compare the ratio of $rejects ad/dp ratio to the minimum ad/dp 
                                #for this allele in @samples
                                if ( $min_allele_ratios{$i} == 0){
                                    #don't div by 0, but if $ratio is greater than 0 filter
                                    $reject_alleles{$i}++ if $ratio;
                                }elsif ( $ratio/$min_allele_ratios{$i} >= $allele_ratio_cutoff){
                                    $reject_alleles{$i}++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    my %r_allele_counts = ();
    if ($af and @reject){
        #if using allele frequency filtering delete any allele from %reject_alleles that has an
        #allele frequency less than $af
        if ($total_reject > 0){
            foreach my $k (keys %reject_alleles){
                if ($reject_alleles{$k}/$total_reject < $af){
                    delete $reject_alleles{$k};
                }
            }
        }
    }elsif($af){
        %r_allele_counts = $vcf_obj->countAlleles(minGQ => $unaff_genotype_quality);
    }
    if (not @samples){
        my @ref_alt = $vcf_obj->readAlleles();
        for (my $i = 1; $i < @ref_alt; $i++){
            push @{$alleles{alt}}, $i if not exists $reject_alleles{$i};
        }
        next LINE if (keys%alleles < 1);#i.e. if only reference (0) or no calls were found
    }
    my %var_call;
    my %breached;
    my %count;
    my %genotypes = $vcf_obj->countGenotypes();
    foreach my $samp (keys%alleles){
    #if we're looking for alleles that match in ALL samples than we only need to check a single hash entry
ALLELE: foreach my $allele (@{$alleles{$samp}}){
            next ALLELE if ($allele == 0);
            next ALLELE if (exists $reject_alleles{$allele});
            $count{$allele}++; #count number of unrejected alleles for comparison with threshold by storing as key of %count hash
            if ($threshold){
                my $t = 0;
                foreach my $k (keys %genotypes){
                    $t += $genotypes{$k} if ($k =~ /(^$allele[\/\|][\.\d+]|[\.\d+][\/\|]$allele)$/);
                }
                if ($t > $threshold){
                    $breached{$allele}++;
                    next ALLELE;
                }
            }
            if ($af){
                if (exists $r_allele_counts{$allele} and $total_reject > 0){
                    next ALLELE if $r_allele_counts{$allele}/$total_reject >= $af;
                }
            }
            my $allele_matches = 0;
            foreach my $sample (keys %alleles){
                $allele_matches++ if (grep (/^$allele$/, @{$alleles{$sample}})); #we're counting the no. of samples with matching allele, not no. of occcurences (zygosity) of allele
            }
            if ($check_presence_only){#don't worry if all samples have variant or not if checking presence only
                $var_call{$allele}++ ;
            }elsif($num_matching){
                $var_call{$allele}++ if $allele_matches >= $num_matching;
            }else{
                $var_call{$allele}++ if ($allele_matches == keys%alleles); #i.e. if all of our called '--keep' sample genotypes match this allele in either het or homo state
            }
        }
    }
    next LINE if (keys %var_call < 1 );#if we don't have at least one valid variant allele
    if (keys %count and keys %breached){
        next LINE if (keys %count == keys %breached);
    }
    print $OUT "$line\n";
}

if ($progressbar){
        $progressbar->update($vcf_obj->countLines("variants")) if $vcf_obj->countLines("variants") >= $next_update;
}
$vcf_obj->DESTROY();
close $OUT; 

#################################################
sub check_samples{
    my ($sample_ref) = @_;
    my @found;
    my @not_found;
    foreach my $s (@$sample_ref){
        if (not $vcf_obj->checkSampleInVcf($s)){
            push @not_found, $s;
        }else{
            push @found, $s;
        }
    }
    return (\@found, \@not_found);
}
 
#################################################

=head1 NAME

filterOnSample.pl - filter variants in vcf that belong to specific samples.

=cut

=head1 SYNOPSIS

    filterOnSample.pl --input [var.vcf] --samples [samples to keep variants if present in all] --reject [samples to reject variants from if present in any] 

=head1 ARGUMENTS

=over 8 

=item B<-i    --input>

vcf file input.

=item B<-o    --output>

output filename.

=item B<-s    --samples>

IDs of samples to keep variants from. Variants will be kept only if present in ALL of these samples in either heterozygous or homozygous states unless --presence or --num_matching flags are set.  Samples must be entered as contained in the vcf header.

=item B<-p    --presence>

Use this flag to print variants present in any sample specified by the --samples option rather than variants present in all.

=item B<-n    --num_matching>

Use this flag to print variants present in at least this many samples rather than only variants present in all.

=item B<-r    --reject>

IDs of samples to reject variants from. Variants will be rejected if present in ANY of these samples unless --allele_frequency is set.

=item B<-x    --reject_all_except>

Reject variants present in all samples except these. If used without an argument all samples in VCF that are not specified by --samples argument will be used to reject variants. If one or more samples are given as argument to this option then all samples in VCF that are not specified by --samples argument or this argument will be used to reject variants.

=item B<-f    --frequency>

Reject variants if the allele frequency (decimal value between 0.00 and 1.00) in the VCF is equal to or greater than this value. If --reject or --reject_all_except arguments are used only the relevent samples will be counted when calculating allele frequency. Otherwise, the allele frequency will be calculated for all samples with a variant call.

=item B<-t    --threshold>

Reject variants present in more than this number of samples in the VCF. Counts all samples in the VCF irrespective of whether they are specified by the --samples or any other argument.

=item B<-c    --confirm_missing>

Use this flag to look only for variants that are present only in --samples and are confirmed absent from all --reject samples. This means that as well as filtering variants with alleles present in --reject samples, variants that contain no calls (or genotype qualities below the --un_quality threshold) in --reject samples will also be filtered. In this way you may identify likely de novo variants in a sample specified by --samples by specifying parental samples with the --reject option, thus avoiding variants where there is not sufficient information to confirm de novo inheritance.

=item B<-d    --depth_allele_cutoff>

Fraction cut-off for allele depth. When specified, the allele depth will be assessed for all samples specified using the --reject argument and any allele with a proportion of reads greater than or equal to this value will be rejected even if the sample genotypes are not called by the genotyper. For example, a value of 0.1 would reject any allele making up 10 % of reads for a sample specified by --reject even if the sample is called as homozygous reference. Must be a value between 0.0 and 1.0. 

=item B<-z    --allele_ratio_cutoff>

Relative fraction cut-off for allele depth between --samples and --reject.  When specified, for each allele the proportion allele depth will be calculated for all samples specified by --samples argument (if they have a genotype call for the given allele) and the minimum allele proportion identified amongst those samples. Similarly, the fraction allele depth for any of the samples specified using the --reject argument (even where the genotype has not been called) will be calculated. If the ratio of the fraction allele depth for any of the --reject samples compared to the minimum allele proportion in the --samples samples is equal to or higher than the value specified for this argument the allele will be filtered.

For example, a variant allele might have a value for [allele depth]/[depth] in one --samples sample of 0.2 (i.e. the allele makes up 20 % of reads at this site for one sample). A --reject sample might have a a value for [allele depth]/[depth] value of 0.1 (i.e. the allele makes up 10 % of reads at this site for this sample). If the value for --allele_ratio_cutoff is given as 0.5 the allele will be filtered as the proportion of reads for the --reject sample is half that of a --samples sample. If the value for --allele_ratio_cutoff is given as 1.0 alleles will only be filtered if the proportion of reads for a --reject sample is equal to or higher than that of a samples --sample. If the value for --allele_ratio_cutoff is given as 0.1 alleles will be filtered if the proportion of reads for a --reject sample is 1/10th or higher than that of a samples --sample, and so on. In this manner you can filter variants if there is a similar level of evidence for the presence of the same variant in a --reject sample even if the genotype hasn't been called, in order to remove variants which are either likely spurious calls in --samples samples or spurious no-calls in --reject samples.

=item B<-a    --aff_quality>

Minimum genotype qualities to consider for samples specified by --samples argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 20. Overrides any value specified by --quality argument.

=item B<-u    --un_quality>

Minimum genotype qualities to consider for samples specified by --reject argument only. Any sample call with a genotype quality below this threshold will be considered a no call. Default is 20. Overrides any value specified by --quality argument.

=item B<-q    --quality>

Minimum genotype quality for all samples - sets both --aff_quality and --un_quality to this value.

=item B<-e    --existing>

Use this flag to cause the program to ignore non-existing samples rather than exiting.

=item B<-b    --progress>

Show a progress bar.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page

=back
=cut

=head1 EXAMPLES

 filterOnSample.pl -i [var.vcf] -s Sample1
 (only print variants if Sample1's genotype has a variant allele)

 filterOnSample.pl -i [var.vcf] -s Sample1 -r Sample2 Sample3
 (look for variants in Sample1 but reject variants also present in Sample2 or Sample3)
 
 filterOnSample.pl -i [var.vcf] -s Sample1 -r Sample2 Sample3 -a 20 -u 30
 (as above but variants are only kept if Sample1 has a genotype quality score of 20 or higher)
 
 filterOnSample.pl -i [var.vcf] -s Sample1 -r Sample2 Sample3 -a 20 -u 30
 (as above but variants are only rejected when alleles are present in Sample2 or Sample3 with a genotype quality score of 30 or higher)
 
 filterOnSample.pl -i [var.vcf] -s Sample1 -r Sample2 Sample3 -t 4
 (same but also reject variants present in 4 or more samples)
 
 filterOnSample.pl -i [var.vcf] -s Sample1 -x 
 (look for variants in Sample1 and reject variants if present in any other sample in the VCF)
 
 filterOnSample.pl -i [var.vcf] -s Sample1 -x Sample2 
 (look for variants in Sample1 and reject variants if present in any other sample in the VCF except for Sample2)
 
 filterOnSample.pl -i [var.vcf] -s child -r mum dad -c -q 30
 (look for apparent de novo variants in child, ignoring variants for which mum or dad have no calls or genotype qualities below threshold genotype quality of 30)

=head1 DESCRIPTION

This program reads a VCF file and filters variants depending on which samples contain variant. Samples to keep variants from can be specified using --samples (-s) and samples to reject variants from can be specified with --reject (-r). 

=cut

=head1 AUTHOR

David A. Parry

University of Leeds

=head1 COPYRIGHT AND LICENSE

Copyright 2013  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


=cut
