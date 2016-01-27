#!/usr/bin/perl

# input: ./runphmmer.abacus.pl uniref_directory query_proteome

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my $uniref_dir = shift @ARGV;
my $proteome = shift @ARGV;

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
		elsif ($count == 100001) {
			$count = 0;
			close OUT;
			$block++;
			open OUT, "> $uniref_dir/uniref_pieces/uniref$block.fasta";
		}
		else {
			$count++;
		}
	}
	print OUT $_;
}

# counts how many pieces Uniref has been split into
my $chunks = `ls $uniref_dir/uniref_pieces/* | wc -l`;
chomp $chunks;
if ($chunks =~ /\s(\S+)$/) {
	$chunks=$1;
}

system "mkdir alignments";

# run query proteome against each chunk of Uniref using qsub
system "mkdir -p $HOME/$$";
open OUT, "> runjob.sh";
print OUT "#!/bin/bash
phmmer -A $HOME/alignments/phmmer\$SGE_TASK_ID $proteome$uniref_dir/uniref_pieces/uniref\$SGE_TASK_ID\.fasta
";
close OUT;
system "chmod +x runjob.sh";
system "qsub -pe multi_thread 8 -l h_vmem=2G -N ph$$ -q \42all.q\42 -v PATH -v PERL5LIB -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.sh";
