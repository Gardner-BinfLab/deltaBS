#!/usr/bin/perl

# input: ./runjackhmmer.abacus.pl uniref_directory query_proteome
# phmmer would make more sense since we're only doing one iteration, but phmmer results files don't have the name of the query sequence in them so there's no way to trace back to the original sequence

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my $uniref_dir = shift @ARGV;
my $proteome = shift @ARGV;

my $chunks = `ls $uniref_dir/uniref_pieces/ | wc -l`;
chomp $chunks;
if ($chunks =~ /\s(\S+)$/) {
	$chunks=$1;
}

system "mkdir alignments";

# run query proteome against each chunk of Uniref using qsub
system "mkdir -p $HOME/$$";
open OUT, "> runjob.sh";
print OUT "#!/bin/bash
jackhmmer -N 1 -A $HOME/alignments/jackhmmer\$SGE_TASK_ID $proteome $uniref_dir/uniref_pieces/uniref\$SGE_TASK_ID\.fasta
";
close OUT;
system "chmod +x runjob.sh";
system "qsub -pe multi_thread 8 -l h_vmem=2G -N jh$$ -q \42all.q\42 -v PATH -v PERL5LIB -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.sh";
