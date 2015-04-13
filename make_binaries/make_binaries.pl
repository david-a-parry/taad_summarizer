#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy::Recursive qw(rcopy);
use FindBin qw($Bin);

my @needs_tabix = qw (
        annotateSnps.pl
        filterOnEvsMaf.pl
        filterOnSample.pl
        filterVcfOnVcf.pl
        getVariantsByLocation.pl
        sampleCallsToInfo.pl
);
my $dir = "$Bin/../";
opendir (my $DIR, $dir) or die "Cannot read current directory $dir: $!\n";
my @files = readdir($DIR); 
close $DIR;
@files = grep {$_ !~ /for_binaries|make_binaries/ } @files;
chdir $dir or die "Cannot move to directory $dir: $!\n";
my $bin_dir = "for_binaries_$^O";
mkdir($bin_dir) or die "$!\n";
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f =~ /\.pl$/){
        print STDERR "Making refactored copy of $f...\n";
        (my $exe = $f) =~ s/\.pl$//;
        my $out = "$bin_dir/$f";
        open (my $IN, $f) or die "Can't read file $f: $!\n";
        open (my $OUT, ">$out") or die "Can't open $out for writing: $!\n";
        while (my $line = <$IN>){
            $line =~ s/$f/$exe/g; 
            print $OUT $line;
        }
        close $IN;
        close $OUT;
    }else{
        print STDERR "Copying $f...\n";
        rcopy($f, "$bin_dir/$f") or die "error copying $f: $!\n"; 
    }
}
#we only make binaries now because we need to be sure that 
#the libs folder has already been copied
foreach my $f (@files){
    next if $f =~ /^\./;
    if ($f =~ /\.pl$/){
        (my $exe = $f) =~ s/\.pl$//;
        chdir($bin_dir);
        my $pp_cmd = "pp -c $f -o $exe";
        if (grep {$_ eq $f} @needs_tabix){
            $pp_cmd .= " -M Tabix";
        }
        print STDERR "Making binary with command: $pp_cmd\n";
        system($pp_cmd); 
        if ($?){
            print STDERR "WARNING - $pp_cmd exited with status $?\n";
        }else{
            print STDERR "Done.\n"; 
        }
        chdir("..");
    }
}
 