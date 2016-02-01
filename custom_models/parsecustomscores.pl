#!/usr/bin/perl

use warnings;
use strict;

# input: ./parsecustomscores.pl [path to ortholog list] [directory with hmmsearch files in it]
# output: file containing custom scores, and ortholog list of those genes which didn't receive custom scores, for Pfam analysis

# read in models and their orthologs for scoring then print out a table of scores for each representative

my $orthlist = shift @ARGV;
my $dir = shift @ARGV;

my %lookup;
my @searchfiles = `ls $dir*.search`;
foreach my $search (@searchfiles) {
	open SEARCH, $search;
	while (<SEARCH>) {
		chomp;
		next if ($_ =~ /#/);
		my @split = split /\s+/;
		$lookup{$split[3]}{$split[0]} = $split[7];
	}
}

my %orthologs;
open ORTHS, $orthlist or die "Couldn't open orthlist";
while (<ORTHS>) {
	chomp;
	my @split = split /\t/;
	if ($#split >0) {
		$orthologs{$split[0]} = $split[1];
	}
}

open OUT, "> customscores.tsv";
open MISS, "> no_cust_score.txt";
my $index2;
foreach my $gene (keys(%orthologs)) {
	next if ($orthologs{$gene} eq "NA");
	if (defined($lookup{$gene}{$gene}) && defined($lookup{$gene}{$orthologs{$gene}})) {
		my $dbs = $lookup{$gene}{$gene} - $lookup{$gene}{$orthologs{$gene}};
		print OUT "$gene\t$orthologs{$gene}\t$lookup{$gene}{$gene}\t$lookup{$gene}{$orthologs{$gene}}\t$dbs\n";
	}
	else {
		if (defined($orthologs{$gene}) && $orthologs{$gene} ne "NA") {
			print MISS "$gene\t$orthologs{$gene}\n";
		}
	}
}
close OUT;
close MISS;
