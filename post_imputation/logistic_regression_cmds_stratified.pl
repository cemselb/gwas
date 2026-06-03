#!/usr/bin/perl
use strict;
use warnings;

my $cohort_name   = "cohortName";
my $analysis_name = "case_control";
my $input_dir     = "..";
my $out_dir_base  = "logit";
my $chr_start     = 1;
my $chr_end       = 22;

##### stratification
# for standard analysis, leave empty: my @subgroups = (); for sex-stratified analysis, use:   my @subgroups = ("males", "females");
my @subgroups = ("males", "females");

#base phenotype file only used if @subgroups is empty
my $base_phenotype = "phenotype.txt";

my $suffix = @subgroups ? join("_", @subgroups) : "chrs";
my $output_script = "${cohort_name}_logit_${suffix}_cmds.sh";


open my $OUT, '>', $output_script
    or die "Could not open $output_script for writing: $!";

foreach my $chr ($chr_start .. $chr_end) {
    
    if (@subgroups) {
        foreach my $group (@subgroups) {
            my $pheno = "${cohort_name}_${group}_phenotype.txt";
            
            my @cmd = (
                "mlogit",
                "-i ${input_dir}/${cohort_name}_${chr}.out.gz",
                "-p $pheno",
                "-o ${out_dir_base}/${chr}/${cohort_name}_${analysis_name}_${chr}_${group}.output.txt"
            );
            print $OUT join(" ", @cmd) . "\n";
        }
    } 
    else {
        my @cmd = (
            "mlogit",
            "-i ${input_dir}/${cohort_name}_${chr}.out.gz",
            "-p $base_phenotype",
            "-o ${out_dir_base}/${chr}/${cohort_name}_${analysis_name}_${chr}.output.txt"
        );
        print $OUT join(" ", @cmd) . "\n";
    }
}

close $OUT;
print "Commands written to: $output_script\n";
