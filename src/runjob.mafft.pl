#!/usr/bin/perl
use warnings;
use strict;
my $HOME = $ENV{"HOME"};

open LIST, "jobFiles.list" or die "Error: $!\n";
my @labels;
while (<LIST>) {
	chomp;
	my $name = $_;
	push(@labels, $name);
}
my $whichFile = $ENV{SGE_TASK_ID} - 1;
my $filename = $labels[$whichFile];
system "$HOME/bin/mafft $HOME/filtered/$filename.fasta > $HOME/filteredandaligned/$filename.afa";


