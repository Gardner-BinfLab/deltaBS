#!/usr/bin/perl

use warnings;
use strict;

#input: ./customscustomscoring.pl [path to custom model hmm file] [path to ortholog list] [path to lookup table of models]

# read in orthlist, make a list of all custom models in database
# check each line of orthlist to see if there are any custom models for it
# print orthlist for each model, with the model in the first column

my $modelfile = shift @ARGV;
my $orthlist = shift @ARGV;
my $lookup = shift @ARGV;

my @models;
my %lookup;

open MODELS, $modelfile;
while (<MODELS>) {
	if ($_ =~ /NAME\s+(\S+)$/) {
		push @models, $1;
	}
}
close MODELS;
print "$#models\n";

open ORTHLIST, $orthlist;
while (<ORTHLIST>) {
	my @split = split /\s+/;
	foreach my $gene (@split) {
		@{$lookup{$gene}} = @split;
	}
}

open HITS, "> $lookup";

foreach my $model (@models) {
	if (defined($lookup{$model})) {
		print HITS $model;
		foreach my $gene (@{$lookup{$model}}) {
			print HITS "\t", $gene;
		}
		print HITS "\n";
	}
	else {
		print "Model not found in lookup: $model";
	}
}
close HITS;