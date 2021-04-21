#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# usage: ./parse_bitscores.pl <orthlist> <path to .search files>

print "Pulling together bitscores for all orthologous genes. This could take some time...\n";

my $orthlist = shift @ARGV;
my %orthologs;
my @strains;

open ORTH, $orthlist;
while (<ORTH>) {
    chomp;
    my @split = split /\",\"/;
    for (@split) {s/\"//g};
    for (@split) {s/___[0-9]*//g};
    for (@split) {s/\r//g};
    # skip the header and any genes only present in one strain
    if ($split[0] eq "Gene") {
        @strains = @split[14..$#split];
    } else {
        next if ($split[3] ==1);
        foreach my $gene (14..$#split) {
            @{$orthologs{$split[$gene]}} = ($split[0], $gene-14);
        }
    }
}
# print Dumper(\%orthologs);
close ORTH;

my $folder = shift @ARGV;
my @files = `ls $folder/*.search`;
my %bitscores;

foreach my $file (@files) {
	chomp $file;
	print $file;
    open IN, $file;
    while (<IN>) {
        chomp;
        next if ($_ =~ /^#/);
        my @split = split /\s+/;
        for (@split) {s/\"//g};
        if (defined($orthologs{$split[0]}) && $split[6]<0.001) {
            $bitscores{$orthologs{$split[0]}[0]}{$split[3]}[$orthologs{$split[0]}[1]] = $split[7];
        }
    }
}
# print Dumper (\%{$bitscores{mreB}});
close IN;

open OUT, "> bitscores.tsv";
open MODELS, "> models_used.tsv";
print OUT "\t", join("\t", @strains), "\n";
foreach my $gene (keys(%bitscores)) {
    my $bestscore = 0;
    my $bestmodel;
    foreach my $model (keys(%{$bitscores{$gene}})) {
        my $sum;
        map { $sum += $_ } grep {defined} @{$bitscores{$gene}{$model}};
       if ($sum > $bestscore) {
            $bestmodel = $model;
           $bestscore = $sum;
        }
    }
    print OUT $gene, "\t", join("\t", map { $_ // '' } @{$bitscores{$gene}{$bestmodel}}), "\n";
    print MODELS $gene, "\t", $bestmodel, "\n";
}
