#!/usr/bin/perl

# input: ./runphmmer.abacus.pl uniref_directory

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my $uniref_dir = shift @ARGV;

# split Uniref90

open UNIREF, "$uniref_dir/uniref90.fasta" or die "Couldn't open Uniref file for splitting";
system "mkdir $uniref_dir/uniref_pieces";
my $count = 0;
my $block = 0;
# go through file and split the sequence database up into chunks of 100,000 sequences
while (<UNIREF>) {
	if ($_ =~ />/) {
		if ($count == 0) {
			$block++;
			open OUT, "> $uniref_dir/uniref_pieces/uniref$block.fasta";
		}
		if ($count == 1000001) {
			$count = 0;
			close OUT;
			$block++;
			open OUT, "> $uniref_dir/uniref_pieces/uniref$block.fasta";
		}
		$count++;
	}
	print OUT $_;
}

