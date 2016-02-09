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
system "/opt/bio/hmmer/bin/hmmbuild --informat afa $HOME/models/$filename.hmm $HOME/filteredandaligned/$filename.afa";

