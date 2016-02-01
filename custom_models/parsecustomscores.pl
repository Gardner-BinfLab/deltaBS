#!/usr/bin/perl

use warnings;
use strict;

# input: ./parsecustomscores.pl [path to lookup table of models] [path to ortholog list] [output file name]

# read in models and their orthologs for scoring then print out a table of scores for each representative

my $lookup = shift @ARGV;
my $orthlist = shift @ARGV;
my $outfile = shift @ARGV;

my %queries;

open IN, $lookup;
while (<IN>) {
	chomp;
	my @split = split /\t/;
	my $model = shift @split;
	if ($model =~ /^(\S+)/) {
		$model = $1;
	}
	@{$queries{$model}} = @split;
	
}

my %lookup;

my @searchfiles = `ls */*.search`;
foreach my $search (@searchfiles) {
	open SEARCH, $search;
	while (<SEARCH>) {
		chomp;
		next if ($_ =~ /#/);
		my @split = split /\s+/;
		$lookup{$split[3]}{$split[0]} = $split[7];
	}
}

open ORTHS, $orthlist;
while (<ORTHS>) {
	chomp;
	my @split = split /\t/;
	foreach my $gene (@split) {
		$orthologs{$split[0]} = $split[1];
	}
}

open OUT, "> $comp/customscores.tsv";
open MISS, "> $comp/no_cust_score.txt";
my $index2;
foreach my $gene (keys(%orthlist)) {
	if (defined($lookup{$gene}{gene}) && defined($lookup{$gene}{$orthologs{$gene}})) {
		my $dbs = $lookup{$gene}{$gene} - $lookup{$gene}{$orthologs{$gene}};
		print OUT "$gene\t$orthologs{$gene}\t$lookup{$gene}{$gene}\t$lookup{$gene}{$orthologs{$gene}}\t$dbs\n";
	}
	else {
		if (defined($orthologs{$gene}) && $orthologs{$gene} ne "NA") {
			print MISS "$1\t$orthologs{$gene}\n";
		}
	}
	
	}
}
close OUT;
close MISS;
