#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my @files = `ls $HOME/filteredandaligned/*`;
open OUT, "> jobFiles.list";
foreach my $file (@files) {
	chomp $file;
	if ($file =~ /$HOME\/filteredandaligned\/(\S+).afa/) {
		print OUT "$1\n";
	}
}
close OUT;
my $chunks = `ls $HOME/filteredandaligned/* | wc -l`;
chomp $chunks;
if ($chunks =~ /\s(\S+)$/) {
	$chunks=$1;
}
system "mkdir models";

system "mkdir -p /home/new26/$$";
system "qsub -pe multi_thread 8 -N hmmbuild$$ -S /usr/bin/perl -q \42all.q\42 -V -v PATH -v PERL5LIB -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.hmmbuild.pl";

