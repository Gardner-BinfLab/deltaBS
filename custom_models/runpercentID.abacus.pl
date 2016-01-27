#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $HOME = $ENV{"HOME"};

my $chunks = `ls $HOME/alignments/* | wc -l`;
chomp $chunks;
if ($chunks =~ /\s(\S+)$/) {
	$chunks=$1;
}

system "mkdir -p /home/new26/$$";
open OUT, "> runjob.sh";
print OUT "#!/bin/bash
./percentID.pl $HOME/alignments/jackhmmer\$SGE_TASK_ID
";
close OUT;
system "chmod +x runjob.sh";
system "qsub -pe multi_thread 8 -l h_vmem=2G -N pid$$ -q \42all.q\42 -v PATH -v PERL5LIB -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.sh";
