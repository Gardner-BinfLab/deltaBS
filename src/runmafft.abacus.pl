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

#		$$ is process ID		$SGE_TASK_ID is the task ID within the process
#echo \42SGE_TASK_ID:[\$SGE_TASK_ID]\42 >> /home/new26/$$/SGE_TASK_IDs


open OUT, "> runjob.sh";
print OUT "#!/usr/bin/perl
use warnings;
use strict;

open LIST, "jobFiles.list" or die "Error: $!\n";
my @labels;
while (<LIST>) {
	chomp;
	my $name = $_;
	push(@labels, $name);
}
my $whichFile = $ENV{SGE_TASK_ID} - 1;
my $filename = $labels[$whichFile];
system "/home/new26/bin/mafft /home/new26/filtered_noquery/40percent/$filename.fasta > /home/new26/filteredandaligned_noquery/40percent/$filename.afa";

";
close OUT;


system "mkdir -p /home/new26/$$";
system "qsub -pe multi_thread 8 -N mafft$$ -S /usr/bin/perl -q \42all.q\42 -V -v PATH -v PERL5LIB -v MAFFT_BINARIES -e $HOME/$$ -o $HOME/$$/$$\_qsub.out -t 1-$chunks runjob.pl";

#/home/new26/bin/mafft /home/new26/filtered/30percent/$genes[$SGE_TASK_ID-1].fasta > /home/new26/filteredandaligned/30percent/$genes[$SGE_TASK_ID-1].afa
