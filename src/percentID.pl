#!/usr/bin/perl

# takes output from each phmmer search and prints filtered sequences to file

use warnings;
use strict;
use Data::Dumper;

# arguments: filename 

my $filename = shift @ARGV;		# jackhmmer chunk
my $outdir = "filtered";		# filtered_sequences
my $HOME = $ENV{"HOME"};
system "mkdir $HOME/$outdir";
system "mkdir $HOME/troubleshooting";

my %sequences;
my %gaps;
my %percentids;
my $seqname;
my $original = "";
my %overthreshold;
my %matchcounts;
my @genes;

open IN, $filename;
while(<IN>) {
    chomp;
	if ($_ =~ /=GF\ ID\ (\S+)-i1/) {
		$original = $1;
		push @genes, $1;
		print $1, "\n";
	}
	elsif ($_ =~ /^#/) {
		next;
	}
	# goes through each aligned sequence and records whether each position is a gap or not, as well as the sequence
	elsif ($_ =~ /(\S+)\s+(\S+)$/){
		my @line = split("", $2);
		foreach my $res (@line) {
			if (isAA ($res)) {
				push @{$gaps{$original}{$1}}, 0;
			}
			else {
				push @{$gaps{$original}{$1}}, 1;
			}
		}
		push @{$sequences{$original}{$1}}, @line;
	}

}

close IN;

print Dumper (\@genes);

# go through each hit and determine the percentage identity
foreach my $gene (@genes) {
	open PIDS, ">> $HOME/troubleshooting/$gene.pids.txt";
	foreach my $seq (keys(%{$sequences{$gene}})){
		next if ($seq eq "$gene");
		my $alignmentlength = $#{$sequences{$gene}{$seq}} + 1;
		$matchcounts{$seq}=0;
		foreach my $residue (0..$#{$sequences{$gene}{$seq}}) {
			if ($gaps{$gene}{$gene}[$residue] == 1 && $gaps{$gene}{$seq}[$residue] == 1) {
				$alignmentlength--;
			}
			elsif (uc $sequences{$gene}{$gene}[$residue] eq uc $sequences{$gene}{$seq}[$residue]) {
				$matchcounts{$seq}++;
			}
		}
		my $percentid = $matchcounts{$seq}/$alignmentlength;
		$percentids{$gene}{$seq} = $percentid;
		printf PIDS "%0.2f\t$seq\n", ($matchcounts{$seq}/$alignmentlength*100);
	}
#	close PIDS;
	#produce output file once all PIDs have been calculated
	my $outfile = "$HOME/filtered/$gene.fasta";
	print $outfile, "\n";
	# print query sequence when file is created - DON'T!!!
#	if (-e $outfile) {
		open OUT, ">> $outfile";
#	}
#	else {
#		open OUT, "> $outfile";
#		print OUT ">$gene\n";
#		foreach my $res (0..$#{$sequences{$gene}{$gene}}) {
#			if (isAA($sequences{$gene}{$gene}[$res])) {
#				print OUT $sequences{$gene}{$gene}[$res];
#			}
#		}
#		print OUT "\n";
#	}
#	if (-e $outfile2) {
#	}
#	else {
#		open OUT2, "> $outfile2";
#		print OUT2 ">$gene\n";
#		foreach my $res (0..$#{$sequences{$gene}{$gene}}) {
#			if (isAA($sequences{$gene}{$gene}[$res])) {
#				print OUT2 $sequences{$gene}{$gene}[$res];
#			}
#		}
#		print OUT2 "\n";
#	}
	# then print all qualifying sequences
	foreach my $seq (keys(%{$percentids{$gene}})) {
		next if ($seq eq $original);
		if($percentids{$gene}{$seq}>0.1) {				# 40% ID cutoff is hard-coded in here - edit at own risk
			print OUT ">$seq\n";
			foreach my $res (0..$#{$sequences{$gene}{$seq}}) {
				if (isAA($sequences{$gene}{$seq}[$res])) {
					print OUT $sequences{$gene}{$seq}[$res];
				}
			}
			print OUT "\n";
		}
	}
}


############################
#isAA: check if single character belongs to the IUPAC amino-acid code.
sub isAA {
    
    my $aa = shift;
    return 0 if (not defined($aa) or length($aa) != 1);
    my $iupac = 'ABCDEFGHIKLMNPQRSTVWXY';
    
    if ($aa=~/[$iupac]/i){
        return 1;
    }
	else {
		return 0;
	}
}
