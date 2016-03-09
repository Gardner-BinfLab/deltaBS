#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my @files = `ls $HOME/filtered/*`;
open OUT, "> jobFiles.list";
foreach my $file (@files) {
	chomp $file;
	if ($file =~ /$HOME\/filtered\/(\S+).fasta/) {
		print OUT "$1\n";
	}
}
close OUT;
my $chunks = `ls $HOME/filtered/* | wc -l`;
chomp $chunks;
if ($chunks =~ /\s(\S+)$/) {
	$chunks=$1;
}

system "mkdir filteredandaligned";

system "mkdir -p /home/new26/$$";
system "qsub -pe multi_thread 8 -N mafft$$ -S /usr/bin/perl -q \42all.q\42 -V -v PATH -v PERL5LIB -v MAFFT_BINARIES -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.mafft.pl";

