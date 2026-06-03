#!/usr/bin/perl
use strict;
use warnings;

my $cohort_name   = "cohortName";
my $analysis_name = "case_control";
my $region_type   = "nonPAR";        # nonPAR, PAR
my $impute_state  = "phased";        # phased, unphased
my $chr           = 23;

my $base_dir = "impute/";

my $output_script = "${cohort_name}_logit_commands_chr${chr}.sh";

my %chr_int = ( 
    chr23 => [
        "2699555 7700000",   "7700001 12700001",  "12700002 17700002",
        "17700003 22700003", "22700004 27700004", "27700005 32700005",
        "32700006 37700006", "37700007 42700007", "42700008 47700008",
        "47700009 52700009", "52700010 57700010", "57700011 62700011",
        "62700012 66021550", "68021551 72700013", "72700014 77700014",
        "77700015 82700015", "82700016 87700016", "87700017 92700017",
        "92700018 97700018", "97700019 102700019", "102700020 107700020",
        "107700021 112700021", "112700022 117700022", "117700023 122700023",
        "122700024 127700024", "127700025 132700025", "132700026 137700026",
        "137700027 142700027", "142700028 147700028", "147700029 152700029",
        "152700030 154930230", "66021551 68021550"
    ]
);

open my $OUT, '>', $output_script
    or die "Could not open $output_script for writing: $!";
  
foreach my $range (@{ $chr_int{"chr${chr}"} }) {
    
    my $range_str = $range;
    $range_str =~ s{ }{_};

    my $input_file = "${base_dir}chr${chr}/${cohort_name}_chr${chr}_${region_type}_${impute_state}_${range_str}.txt.gz";
    my $pheno_file = "${base_dir}chr${chr}/${cohort_name}_phenotype_chr${chr}.txt";
    my $out_file   = "${base_dir}chr${chr}/results/${analysis_name}_chr${chr}_${region_type}.${range_str}.output.txt";

    my @cmd = (
        "mlogit", 
        "-i $input_file", 
        "-p $pheno_file",
        "-o $out_file"
    );

    print $OUT join(" ", @cmd) . "\n";
}

close $OUT;
print "Commands written to: $output_script\n";
