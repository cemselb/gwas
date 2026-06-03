#!/usr/bin/perl
use strict;
use warnings;

my $cohort_name   = "cohortName";
my $phenotype     = "phenotype.txt";
my $analysis_name = "case_control";
my $out_dir_base  = "logit";
my $chr_start     = 1;
my $chr_end       = 22;

my $output_script = "${cohort_name}_logit_chrs_cmds.sh";

open my $OUT, '>', $output_script
    or die "Could not open $output_script for writing: $!";

foreach my $chr ($chr_start .. $chr_end) {
    my @cmd = (
        "mlogit",
        "-i ${cohort_name}_${chr}.out.gz",
        "-p $phenotype",
        "-o ${out_dir_base}/${chr}/${cohort_name}_${analysis_name}_${chr}.output.txt"
    );
    
    print $OUT join(" ", @cmd) . "\n";
}

close $OUT;
print "Cmds written to: $output_script\n"; 
