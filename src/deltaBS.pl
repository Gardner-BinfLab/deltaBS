#!/usr/bin/env perl 

use warnings;
use strict;

use Bio::Perl;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
use Bio::DB::Fasta;
use Statistics::Distributions;
use Data::Dumper;

use Getopt::Long;

my ($filetype, $file1, $file2, $pfamannot1, $pfamannot2, $phmmerannot1, $phmmerannot2, $outdir, $orthlist, $hmmer_path, $cpus, $hmm_lib_path, $tmp_dir, $verbose, $post, $dirty, $help);
my $identifier = int(rand(10000));

&GetOptions(
	"f|filetype=s"	     => \$filetype,
	"f1|file1=s"	     => \$file1,
	"f2|file2=s"	     => \$file2,
	"o|outdir=s"	     => \$outdir,
	"ol|orthlist=s"	     => \$orthlist,
        "pa1|pfamannot1=s"   => \$pfamannot1,
        "pa2|pfamannot2=s"   => \$pfamannot2,
        "ph1|phmmerannot1=s" => \$phmmerannot1,
        "ph2|phmmerannot2=s" => \$phmmerannot2,
	"hp|hmmerpath=s"     => \$hmmer_path,
	"c|cpus=i"	     => \$cpus,
	"hd|hmmlibdir=s"     => \$hmm_lib_path,
	"t|tempdir=s"	     => \$tmp_dir,
	"v|verbose"	     => \$verbose,
        "d|dirty"            => \$dirty,
	"p|post"	     => \$post,
	"h|help"	     => \$help
	);


if ($help) {
	&help();
	exit(1);
}

if(not defined($filetype) or ($filetype ne "embl" & $filetype ne "genbank" & $filetype ne "fasta")){
	print "ERROR: filetype must be specified as either 'embl', 'genbank', or 'fasta'\n";
	&help();
	exit(1);
} 
elsif ($verbose) {
	print STDERR "filetype set to $filetype\n";
}

if (not defined($file1) or not defined($file2)){
	print "ERROR: need $filetype files to work with!\n";
	&help();
	exit(1);
}

elsif (not -s $file1 or not -s $file2){
	print "ERROR: [$file1] &/or [$file2] is missing or empty!\n";
	&help();
	exit(1);
}

elsif (not defined($outdir)){
	print "ERROR: output directory required!\n";
	&help();
	exit(1);
}

elsif (not defined($hmmer_path) && (not defined $ENV{'HMMER'})){
	print "ERROR: HMMER path not provided and HMMER environment variable not set!\n";
	&help();
	exit(1);
}

#elsif (not defined($easel_path) && (not defined $ENV{'EASEL'})){
#	print "ERROR: EASEL path not provided and EASEL environment variable not set!\n";
#	&help();
#	exit(1);
#}


elsif (not defined($hmm_lib_path) && (not defined $ENV{'DELTABS'})){
	print "ERROR: HMM library path not set and DELTABS environment variable not set!\n";
	&help();
	exit(1);
}

if(defined $ENV{'HMMER'} && not defined($hmmer_path)){
	warn "WARNING: HMMER path not set, using ",$ENV{'HMMER'},"\n";
	$hmmer_path = $ENV{'HMMER'};
}

if(defined $ENV{'DELTABS'} && not defined($hmm_lib_path)){
	warn "WARNING: HMM library path not set, using ",$ENV{'DELTABS'},"/deltaBS.hmmlib\n";
	$hmm_lib_path = $ENV{'DELTABS'};
}

if (not defined($tmp_dir)){
	warn "WARNING: no temporary directory set, will attempt to use /tmp\n";
	$tmp_dir = "/tmp/DBS.$identifier";
	system("mkdir -p $tmp_dir") == 0 or die "Unable to create $tmp_dir: $!";
} else {
	$tmp_dir = "$tmp_dir/DBS.$identifier";
	print STDERR "Using $tmp_dir as temporary directory.\n" if $verbose;
	system("mkdir -p $tmp_dir") == 0 or die "Unable to create $tmp_dir: $!";
}

if(not defined($cpus)){
	warn "WARNING: number CPUs for hmmer not set, using 1 as default\n";
	$cpus = 1;
}

if (not defined($orthlist)){
	warn "WARNING: no ortholog list provided, will attempt automatic orthology determination with phmmer.\n";
}

if (defined($post)){
	print STDERR "Post-processing enabled, will produce Genome Properties and GO annotations\n" if $verbose;
}

system("mkdir -p $outdir") == 0 or die "Couldn't create $outdir: $!";

#generate fasta files, get annotation hash refs
my ($fasta1, $fasta2);
if($filetype ne "fasta"){
	print STDERR "Generating proteome fasta files in $tmp_dir...\n" if $verbose;
	$fasta1 = "$tmp_dir/DELTABS_embl1.fasta";
	$fasta2 = "$tmp_dir/DELTABS_embl2.fasta";
	my $embl1_genes = &write_CDS_fasta($file1,$filetype,$fasta1);
	my $embl2_genes = &write_CDS_fasta($file2,$filetype, $fasta2);
	print STDERR "done.\n" if $verbose;
} else {
	print STDERR "Input in fasta format, no need to generate proteome files.\n" if $verbose;
	$fasta1 = $file1;
	$fasta2 = $file2;
}
######################################################################
#run hmmscan on the proteomes with Pfam HMMs:

if(not defined($pfamannot1)){
    $pfamannot1 = "$file1-pfam_hmmscan1.tbl";
    print STDERR "Running hmmscan on [$file1] sequences with HMM database...\n" if $verbose;
    system("$hmmer_path/hmmscan -o /dev/null --noali --cpu $cpus --domtblout $pfamannot1 --cut_tc $hmm_lib_path/deltaBS.hmmlib $fasta1 1>&2") == 0 or die "hmmscan failed: $!";
}
else{
    print STDERR"skipping hmmscan\'ing [$file1]. Using [$pfamannot1] instead.\n" if $verbose;
}

if(not defined($pfamannot2)){
    $pfamannot2 = "$file2-pfam_hmmscan1.tbl";
    print STDERR "Running hmmscan on [$file2] sequences with HMM database...\n" if $verbose;
    system("$hmmer_path/hmmscan -o /dev/null --noali --cpu $cpus --domtblout $pfamannot2 --cut_tc $hmm_lib_path/deltaBS.hmmlib $fasta2 1>&2") == 0 or die "hmmscan failed: $!";
}
else{
    print STDERR "skipping hmmscan\'ing [$file2]. Using [$pfamannot2] instead.\n" if $verbose;
}

print STDERR "done hmmscan\47ing.\n" if $verbose;

## want to add something in here that looks at whether the position in the hmm is about the same (or at least output these tables)

#my $pfamannot1 = "tmp/DBS.8392/DELTABS_hmmscan1.tbl";
#my $pfamannot2 = "tmp/DBS.8392/DELTABS_hmmscan2.tbl";
#parse hmmscan results
print STDERR "Parsing hmmscan results...\n" if $verbose;
die "hmmscan results $pfamannot1 empty for $fasta1!" if (-z $pfamannot1);
my $scan1 = &parse_hmmscan_tbl($pfamannot1);
die "hmmscan results $pfamannot2 empty for $fasta2!" if (-z $pfamannot2);
my $scan2 = &parse_hmmscan_tbl($pfamannot2);
print STDERR "done.\n" if $verbose;

######################################################################

#get ortholog list
my %orths;
if(not defined($orthlist)){
	print STDERR "Since no ortholog list provided, predicting orthologs with phmmer...\n" if $verbose;
	my $ref;
	if(defined($phmmerannot1) and defined($phmmerannot2)){
	    $ref = &predict_orths_phmmer($fasta1, $fasta2, $phmmerannot1, $phmmerannot2);
	}
	else{
	    $ref = &predict_orths_phmmer($fasta1, $fasta2);
	}
	%orths = %$ref;
	open OUT, ">", "$outdir/orthlist.dbs";
	foreach my $key (keys(%orths)){
		print OUT $key,"\t",$orths{$key},"\n";
	}
	print STDERR "done. Ad hoc ortholog list printed to $outdir/orthlist.dbs\n" if $verbose;
	close OUT;
} else {
	print STDERR "Reading ortholog list: $orthlist...\n" if $verbose;
	open ORTHS, "<$orthlist";
	while(<ORTHS>){
		next if ($_ =~ /^#/);
		die "Ortholog list must take form ID1\\tID2\\n" if ($_ !~ /.+\t.+\n/);
		chomp;
		my @splat = split /\t/;
		$orths{$splat[0]} = $splat[1];
	}
	close ORTHS;
	print STDERR "Done.\n" if $verbose;
}



#filter orthologs on identical architecture
print STDERR "Determining and filtering domain architecture, calculating delta BS...\n" if $verbose;
my %dbs;
my @score_dist;
my $inc = 0;
my $inc_file = "$outdir/inc_archs.dbs";
unlink $inc_file if ( -e $inc_file);
foreach my $key (keys(%orths)){
	next if(!defined($scan1->{$key}) || !defined($scan2->{$orths{$key}}));
	my $arch1 = &get_domain_arch($scan1->{$key});
	my $arch2 = &get_domain_arch($scan2->{$orths{$key}});
	if(&comp_archs($arch1, $arch2)){
		#domain, eval, score, start, end
		foreach my $i (0..scalar($#$arch1)){
		    my $score = $arch1->[$i][2] - $arch2->[$i][2];#####COMPUTE DELTA BITSCORE HERE! 
			push @score_dist, $score;
			#domain, start1, end1, start2, end2, score
			@{$dbs{$key}[$i]} = ($arch1->[$i][0], $arch1->[$i][3], $arch1->[$i][4],$arch1->[$i][2], $arch2->[$i][3], $arch2->[$i][4],$arch2->[$i][2], $score);
		}
	} else {
		$inc = 1;
		open OUT, ">>", $inc_file;
		print OUT $key,";",$orths{$key},";";
		my $first = 1;
		foreach my $dom (@$arch1){
			if($first){
				print OUT $dom->[0];
				$first = 0;
			}else{
				print OUT ",",$dom->[0];
			}
		}
		print OUT ";";
		$first = 1;
		foreach my $dom (@$arch2){
			if($first){
				print OUT $dom->[0];
				$first = 0;
			}else{
				print OUT ",",$dom->[0];
			}
		}
		print OUT "\n";
	}
}
close OUT;
print STDERR "done.\n" if $verbose;
print STDERR "Incompatible architectures detected in orthologous CDSes, printed to $outdir/inc_archs.dbs\n" if $inc;

die "No compatible architectures found; something horribly wrong - check input files\n" if (scalar @score_dist == 0); #really shouldn't happen...unless you give it the wrong Pfam annotation file

print STDERR "Calculating empirical cutoff:\n" if $verbose;
my $emp_cut = &calc_emp_cutoff(@score_dist);
print STDERR $emp_cut, "\n";

#get population mean
my $mean = &calc_mean(@score_dist);
#get population SD
my $sd = &calc_sd(\@score_dist, $mean);
#shrink mean & SD
print STDERR "Shrinking mean and standard deviation:\n" if $verbose;
print STDERR "original mean: $mean; original sd: $sd\n" if $verbose;
my $shrunken = &shrink_dist(\@score_dist, $sd);
$sd = $shrunken->[0];
$mean = $shrunken->[1];
print STDERR "final mean: $mean; final sd: $sd\n";
#calculate Z-scores and stat sig
#note, using mean of 0 for Z-score, calculated SD
print STDERR "Calculating Z-scores, p-values...\n" if $verbose;
my @outtable;
foreach my $key (keys(%dbs)){
	foreach my $dom (@{$dbs{$key}}){
		my $ex;
		if(abs($dom->[7]) > $emp_cut){
			$ex = 1;
		} else {
			$ex = 0;
		}
		my $z = &calc_z($dom->[7], 0, $sd);
		my $pval;
		if ($z > 0){
			$pval = 2*(Statistics::Distributions::uprob($z));
		} else {
			$pval = 2*(Statistics::Distributions::uprob(abs($z)));
		}
		push @outtable, [$key,$orths{$key}, $dom->[0],$dom->[1],$dom->[2],$dom->[3],$dom->[4],$dom->[5],$dom->[6],$dom->[7],$z,$pval, $ex];
	}
}
print STDERR "done.\n" if $verbose;

print STDERR "Sorting, printing results...\n" if $verbose;
#sort out table, print results!
my @sortedtable = sort{$b->[9] <=> $a->[9]} @outtable;
open OUT, ">$outdir/results.dbs";
print OUT "# MEAN: ",$mean," SD: ", $sd," EMP_CUT: ", $emp_cut,"\n";
print OUT "gene_1\tgene_2\tHMM\tstart_gene_1\tend_gene_1\tbitscore_gene_1\tstart_gene_2\tend_gene_2\tbitscore_gene_2\tdelta-bitscore\tz-score\tp-value\tloss_of_function\n";
foreach my $row (@sortedtable){
	my $first = 1;
	foreach my $entry (@$row){
		if($first){
			print OUT $entry;
			$first = 0;
			next;
		}
		print OUT "\t",$entry;
	}
	print OUT "\n";
}
close OUT;
print STDERR "done.\n" if $verbose;

print STDERR "Cleaning up temp files...\n" if $verbose;
#clean up
if(not defined($dirty)){
    system("rm -rf $tmp_dir") == 0 or die "Couldn't rm $tmp_dir: $!";
    print STDERR "done.\n" if $verbose;
}

if(defined($post)){
	print STDERR "Post-processing results...\n" if $verbose;
	#build HMM to step link hash
	my $step = &build_step_ev_hash($hmm_lib_path);
	#build HMM to GO-term link hash
	my $go = &build_go_hash($hmm_lib_path);
	#build step to genome property link hash
	my $s2p = &build_s2gp_hash($hmm_lib_path);
	#build genome property def hash
	my $gp = &build_gp_hash($hmm_lib_path);
	
	#build HMM GO and genome property output
	open GO, ">$outdir/gomapping.dbs";
	foreach my $row (@sortedtable){
		my $dom = $row->[2];
		$dom =~ s/\.\d+//; #strip off Pfam versioning	
		my $goterms = $go->{$dom};
		foreach my $term (@$goterms){
			#print $dom,"\t",$term,"\n";
			print GO $row->[0],"\t",$row->[1],"\t",$row->[2],"\t",$term,"\t",$row->[9],"\t",$row->[10],"\t",$row->[11],"\n";
		}
	}
	close GO;
	open GP, ">$outdir/genprops.dbs";
	foreach my $row (@sortedtable){
		my $dom = $row->[2];
		$dom =~ s/\.\d+//; #strip off Pfam versioning	
		my $steps = $step->{$dom};
		foreach my $step (@$steps){
			#prop, desc, req?
			my $prop = $s2p->{$step};
			#genprop acc, type, desc 
			my $desc = $gp->{$prop->[0]};
			print GP $row->[0],"\t",$row->[1],"\t",$row->[2],"\t",$desc->[0],"\t",$desc->[1],"\t",$prop->[1],"\t",$desc->[2],"\t",$prop->[2],"\t",$row->[9],"\t",$row->[10],"\t",$row->[11],"\n";
		}
	}
	close GP;
}

print STDERR "Done! Results in $outdir\n";

###################################################
#build_s2gp_hash: parse PROP_STEP.TABLE, return
#hash mapping steps to genome properties
###################################################
sub build_s2gp_hash {
	my ($path) = @_;
	open IN, "<$path/PROP_STEP.TABLE" or die "Couldn't open $path/PROP_STEP.TABLE:$!";
	my %s2gphash;
	while(<IN>){
		chomp;
		my @line = split "\t";
		#step = prop, desc, required?
		$s2gphash{$line[0]} = [$line[1], $line[3],$line[4]];
	}
	close IN;

	return \%s2gphash;
}

###################################################
#build_gp_hash: parse PROP_DEF_WITHOUT_DESCRIPTION_FIELD.TABLE,
#return hash containing genome property description
###################################################
sub build_gp_hash {
	my ($path) = @_;
	open IN, "<$path/PROP_DEF_WITHOUT_DESCRIPTION_FIELD.TABLE" or die "Couldn't open $path/PROP_DEF_WITHOUT_DESCRIPTION_FIELD.TABLE:$!";
	my %gphash;
	while(<IN>){
		chomp;
		my @line = split "\t";
		#prop = genprop acc, type, desc
		$gphash{$line[0]} = [$line[3],$line[2],$line[1]];
	}
	close IN;
	return \%gphash;
}
###################################################
#build_go_hash: parse pfam2go and TIGR_GO_LINK,
#return hash containing domain to GO term mapping
###################################################
sub build_go_hash {
	my($path) = @_;
	open IN, "<$path/pfam2go" or die "Couldn't open $path/pfam2go:$!";
	my %gohash;
	while(<IN>){
		next if ($_ =~ /^!/);
		chomp;
		$_ =~ m/(PF\d\d\d\d\d)/;
		my $dom = $1;
		$_ =~ m/(GO:\d\d\d\d\d\d\d)/;
		my $term = $1;
		push @{$gohash{$dom}}, $term;
	}
	close IN;
	open IN, "<$path/TIGR_GO_LINK" or die "Couldn't open $path/TIGR_GO_LINK:$!";
	while(<IN>){
		chomp;
		my @line = split "\t";
		push @{$gohash{$line[0]}},$line[1];
	}
	close IN;

	return \%gohash;
}

###################################################
#build_step_ev_hash: parse TIGRfam/Genome
#Properties STEP_EV_LINK table, return hash mapping
#domains to steps
###################################################
sub build_step_ev_hash {
	my($path) = @_;
	open IN, "<$path/STEP_EV_LINK.TABLE" or die "Couldn't open $path/STEP_EV_LINK.TABLE:$!";
	my %gphash;
	while(<IN>){
		chomp;
		my @line = split "\t";
		next if ($line[3] !~ /HMM/);
		push @{$gphash{$line[2]}}, $line[1];
	}
	close IN;
	return \%gphash;
}

###################################################
#calc_z: Calculate z-score given value, mean, and
#standard deviation
###################################################
sub calc_z {
	my($value, $mean, $sd) = @_;
	return (($value - $mean)/$sd);
}
###################################################
#shrink_dist: Attempt to remove outliers by
#iteratively removing elements > 5sd's from the
#mean; return shrunk sd and mean
###################################################
sub shrink_dist {
	my ($values, $sd, $iter) = @_;
	if(!defined($iter)){
		$iter = 0;
	} else {
		$iter ++;
	}
	my @new_values;
	my $mean = &calc_mean(@$values);
	my $threshold = 5*$sd;
	foreach my $value (@$values){
		push @new_values, $value if (abs($value - $mean) < $threshold);
	}
	$mean = &calc_mean(@new_values);
	my $new_sd = &calc_sd(\@new_values, $mean);
	print STDERR "new mean: $mean; new sd: $new_sd\n" if $verbose;
	my $returned;
	if($new_sd < $sd - 0.0001){
		die "shrink_dist failed to converge in 1000 iterations\n" if $iter > 1000;
	 	$returned = &shrink_dist(\@new_values, $new_sd, $iter);
	} else {
		return([$new_sd, $mean]);
	}
	return($returned);
}


###################################################
#calc_mean: Calculate mean given an array of values
###################################################
sub calc_mean {
	my (@values) = @_;
	my $sum;
	my $count = 0;
	foreach (@values) {
		if ($_ != 0){
			$sum += $_;
		}
		$count++;
	}
        if ($count > 0){
	    my $mean = $sum / $count;
        } else {
            my $mean = 0;
        }
	return $mean;
}
###################################################
#calc_emp_cutoff:  Identify narrower side of dist,
#return 2.5% cutoff
###################################################
sub calc_emp_cutoff {
	my (@values) = @_;
	my @pos;
	my @neg;

	foreach (@values) {
		push(@pos, $_) if($_ >= 0);
		push(@neg, -$_) if($_ <= 0);
	}
	@pos = sort{$a <=> $b} @pos;
	@neg = sort{$a <=> $b} @neg;
	my($p_cut, $n_cut);
	if($#pos > 0){
		$p_cut = $pos[int($#pos*.975)];
	} else {
		$p_cut = 0;
	}
	if($#neg > 0) {
		$n_cut = $neg[int($#neg*.975)];
	} else {
		$n_cut = 0;
	} 
	return($p_cut) if $p_cut < $n_cut;
	return($n_cut);
}
	
###################################################
#calc_sd: Calculate standard deviation given an
#array of values and the population mean
###################################################
sub calc_sd {
	my ($values, $mean) = @_;
	my ($x, $sd);
	#print "WTF!!\n";
	#print $mean, "MEAN\n";
	#foreach my $value (@$values){
		#print $value, "\n";
	#}
	my $count = 0;
	foreach (@$values) {
		if($_ != 0){
			$x += ($_ - $mean) ** 2;
			$count++;
		}
	}
        if ($counts > 0){
	    $sd = sqrt((1/scalar($count)) * $x);
	    $sd = 1 if $sd < 1;
	    return $sd;
        }
        else {
            return 1;
        }
}
###################################################
#comp_archs: given two domain architectures,
#determine if they are equivalent
###################################################
sub comp_archs {
	my($arch1,$arch2) = @_;
	#print "REF: $arch1\n";
	#print Dumper $arch1;
	#print "REF: $arch2\n";
	#print Dumper $arch2;
	#print "ARCH $arch1 $arch2\n";
	if($#$arch1 != $#$arch2){
		return 0;
	}

	foreach my $i (0..$#$arch1){
		#print "ARCH: $i ", $arch1->[$i][0],"\t",$arch2->[$i][0],"\n";
		if($arch1->[$i][0] ne $arch2->[$i][0]){
		#	print "inc!\n";
			#foreach my $j ($i+1..$#$arch1){
			#	print "ARCH: $j ", $arch1->[$j][0],"\t",$arch2->[$j][0],"\n";
			#}
			return 0;
		}
	}

	return 1;
}

###################################################
#get_domain_arch: given domain hits, resolve 
#overlaps, preferentially use TIGRfams over Pfams
###################################################
sub get_domain_arch {
	my ($hit) = @_;
	my @arch;
	my @sorted = sort {$a->[3] <=> $b->[3]} @$hit; #sort on start coords
	my %segments;
	my $seg_count = 0;
	my @cont;
	
	#define segments: range of overlapping domain hits
	push @{$cont[0]}, 0;
	$segments{0} = 0;
	foreach my $i1 (0..($#sorted - 1)){
		foreach my $i2 (($i1 + 1)..$#sorted){
			if($sorted[$i1][4] > $sorted[$i2][3]){
				$segments{$i2} = $seg_count;
				push @{$cont[$seg_count]}, $i2;
				last;
			} else {
				$seg_count++;
				push @{$cont[$seg_count]}, $i2;
				last;
			}
		}
	}
	
	#resolve segments to single domains
	#domain, eval, score, start, end
	foreach my $seg (0..$seg_count) {
		if($#{$cont[$seg]} == 0){ #trivial case, n domains = 1
			$arch[$seg] = [@{$sorted[$cont[$seg][0]]}];
			next;
		}
		
		#old code to prioritize TIGRfams	
		#my @TIGRs;
		#my @evals;
		#foreach my $i (0..$#{$cont[$seg]}){
		#	push(@TIGRs,$cont[$seg][$i]) if ($sorted[$cont[$seg][$i]][0] =~ /^TIGR/); 
		#	push(@evals, $sorted[$cont[$seg][$i]][1]);
		#}
		#if($#TIGRs == 0){ #if only one TIGRfam, take it
		#	$arch[$seg] = [@{$sorted[$TIGRs[0]]}];
		#	next;
		#}
		#elsif($#TIGRs > 0){ #if multiple TIGRfams, compete them
		#	my $min = 1;
		#	my $winner;
		#	foreach my $dom (@TIGRs){
		#		if($sorted[$dom][1] < $min){
		#			$min = $sorted[$dom][1];
		#			$winner = $dom;
		#		}
		#	}
		#	$arch[$seg] = [@{$sorted[$winner]}];
		#}

		my @NEWs;
		foreach my $i (0..$#{$cont[$seg]}){
			push(@NEWs,$cont[$seg][$i]) if ($sorted[$cont[$seg][$i]][0] =~ /^newfam/); 
			#push(@evals, $sorted[$cont[$seg][$i]][1]);
		}
		if($#NEWs == 0){ #if only one newfam, take it
			$arch[$seg] = [@{$sorted[$NEWs[0]]}];
			next;
		}
		elsif($#NEWs > 0){ #if multiple newfams, compete them
			my $min = 1;
			my $winner;
			foreach my $dom (@NEWs){
				if($sorted[$dom][1] < $min){
					$min = $sorted[$dom][1];
					$winner = $dom;
				}
			}
			$arch[$seg] = [@{$sorted[$winner]}];
		}
		else { #all Pfams, compete them
			my $min = 1;
			my $winner;
			foreach my $dom (@{$cont[$seg]}){
			        if($sorted[$dom][1] < $min){
				    $min = $sorted[$dom][1];
				    $winner = $dom;
				}
			}
			$arch[$seg] = [@{$sorted[$winner]}];
		}
	}
	return \@arch;
}
	

###################################################
#predict_orths_phmmer: given two fasta files run 
#reciprocal phmmer searches; take lowest average
#evalue as indicative of orthology.
#ADD SUM OVER MULTIPLE LOCAL HITS AND A COVERAGE THRESHOLD!!!
###################################################
sub predict_orths_phmmer {
	my($fasta1, $fasta2) = ($_[0], $_[1]);
	my ($line, $seq, $name);
	my ($f1, $f2, $l1, $l2);
	my $f1_tmp = "$tmp_dir/DELTABS_f1.fa";
	my $f2_tmp = "$tmp_dir/DELTABS_f2.fa";
	my $f1_tbl = "$tmp_dir/DELTABS_fa1.tbl";
	my $f2_tbl = "$tmp_dir/DELTABS_fa2.tbl";
	my $f1_domtbl = "$tmp_dir/DELTABS_fa1.domtbl";
	my $f2_domtbl = "$tmp_dir/DELTABS_fa2.domtbl";
	my %orths;
	my %pairs;
	my %coverage;
	
	#read in fasta files
	#$f1 = &read_fasta($fasta1);
	#$f2 = &read_fasta($fasta2);	
	$l1 = &find_seq_lengths($fasta1); 
	$l2 = &find_seq_lengths($fasta2); 
	
	my $l1_no = scalar keys(%$l1); 
	
	if (defined($_[2]) && defined($_[3])){
	    ($f1_domtbl, $f2_domtbl)=($_[2], $_[3]); 
	}
	else {
	#loop through f1 searching each seq against fasta2
	#reciprocally search seqs from f2 which are hit at
	#E < 0.1
	    system("$hmmer_path/phmmer -o $f1_tbl.phmmer --noali --cpu $cpus --domtblout $f1_domtbl --tblout $f1_tbl -E 0.000001 $fasta1 $fasta2") == 0 or die "phmmer failed: $!";
	    system("$hmmer_path/phmmer -o $f2_tbl.phmmer --noali --cpu $cpus --domtblout $f2_domtbl --tblout $f2_tbl -E 0.000001 $fasta2 $fasta1") == 0 or die "phmmer failed: $!";
	}
	
#	my $f1_hits = &parse_phmmer_tbl2($f1_tbl);
#	my $f2_hits = &parse_phmmer_tbl2($f2_tbl);
	my ($f1_hits, $f1_lengths) = &parse_phmmer_tbl3($f1_domtbl);
	my ($f2_hits, $f2_lengths) = &parse_phmmer_tbl3($f2_domtbl);
	
	foreach my $id1 (keys(%$l1)){
	    
	    my $max = 1;
	    #print "predict_orths_phmmer:[$id1]\t";
	    #print "\t2:1[$f2_hits->{$id1}]\n";
	    
	    next if(!$f2_hits->{$id1});
	    my $f2_orth;
	    foreach my $id2 (keys(%$l2)){
		next if(!$f1_hits->{$id2});
		next if(!$f1_hits->{$id2}{$id1});
		next if(!$f2_hits->{$id1}{$id2});
		#print "\t[$id2]";
		#print "HERE!\t";
		$pairs{"$id1:$id2"}=$f1_hits->{$id2}{$id1} + $f2_hits->{$id1}{$id2};
		my @coverage = (($f1_lengths->{$id2}{$id1}{$id2}/$l2->{$id2}), ($f1_lengths->{$id2}{$id1}{$id1}/$l1->{$id1}), ($f2_lengths->{$id1}{$id2}{$id2}/$l2->{$id2}), ($f2_lengths->{$id1}{$id2}{$id1}/$l1->{$id1}));
		
		#print "coverage:[@coverage]\t";
		$coverage{"$id1:$id2"}=minA( @coverage  ); 
		#print "\t\t$id1:$id2 -> pairs:[" . $pairs{"$id1:$id2"} . "]\tcoverage:[" . $coverage{"$id1:$id2"} . "]"; 
	    }
	    #print "\n";
	}
	
	
	my %seen;
	#Sort on sum of reciprocal phmmer bitscores:
	foreach my $p ( sort {$pairs{$b} <=> $pairs{$a}} keys %pairs){
	    
	    my ($id1, $id2) = split(/:/, $p); 
	    next if defined $seen{$id1};
	    next if defined $seen{$id2};
	    if ($pairs{$p} > 20 && $coverage{"$id1:$id2"}>0.75){
		$orths{$id1} = $id2;
		($seen{$id2},$seen{$id1})=(1,1);
		#print "$id1\t$id2\t$pairs{$p}\t" . $coverage{"$id1:$id2"} . "\n";
		#print "YES!\n"; 
	    }
	}
	
	return \%orths;
}


###################################################
#parse_phmmer_tbl: parse phmmer table output
###################################################
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
sub parse_phmmer_tbl {
	my($tbl) = @_;
	my %hits;
	open TBL, "<", $tbl or die "Cannot open $tbl: $!";

	while(<TBL>){
		next if($_ =~ /^#/);
		chomp;
		my @splat = split /\s+/;
		$hits{$splat[0]} = $splat[5];
	}
	close TBL;
	return 0 if(! scalar(keys(%hits)));
	return \%hits;
}

###################################################
#parse_phmmer_tbl2: parse phmmer table output
###################################################
#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
sub parse_phmmer_tbl2 {
	my($tbl) = @_;
	my %hits;
	open TBL, "<", $tbl or die "Cannot open $tbl: $!";
	
	while(<TBL>){
		next if($_ =~ /^#/);
		chomp;
		my @splat = split /\s+/;
		$hits{$splat[0]}{$splat[2]} = 0 if (not defined($hits{$splat[0]}{$splat[2]})); 
		$hits{$splat[0]}{$splat[2]} += $splat[5]; #really should check if hits overlap
	}
	close TBL;
	return 0 if(! scalar(keys(%hits)));
	return \%hits;
}

###################################################
#parse_phmmer_tbl3: parse phmmer domain table output
###################################################
##                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
## target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
##------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#t0001                -             21 SL1344_0001          -             28   5.6e-09   31.9   5.3   1   1   1.3e-12   5.7e-09   31.9   5.3     8    28     1    21     1    21 0.98 -
sub parse_phmmer_tbl3 {
	my($tbl) = @_;
	my (%hits, %lengths);
	open TBL, "<", $tbl or die "Cannot open $tbl: $!";
	
	while(<TBL>){
		next if($_ =~ /^#/);
		chomp;
		my @splat = split /\s+/;
		$hits{$splat[0]}{$splat[3]} = 0 if (not defined($hits{$splat[0]}{$splat[3]})); 
		$hits{$splat[0]}{$splat[3]} += $splat[7]; #really should check if hits overlap
		
		$lengths{$splat[0]}{$splat[3]}{$splat[0]} = 0 if (not defined($lengths{$splat[0]}{$splat[3]}{$splat[0]})); 
		$lengths{$splat[0]}{$splat[3]}{$splat[3]} = 0 if (not defined($lengths{$splat[0]}{$splat[3]}{$splat[3]})); 
		$lengths{$splat[0]}{$splat[3]}{$splat[3]} += ($splat[16] - $splat[15] + 1); #hmm coord -- query
		$lengths{$splat[0]}{$splat[3]}{$splat[0]} += ($splat[18] - $splat[17] + 1); #ali coord -- target 
		
		#printf "parse_phmmer_tbl3: [$splat[0]]:[$splat[3]]\tbits[$splat[7]]\tl1:[%d]\tl2:[%d]\n", ($splat[16] - $splat[15] + 1), ($splat[18] - $splat[17] + 1);
		
	}
	close TBL;
	return 0 if(! scalar(keys(%hits)));
	return (\%hits, \%lengths);
}

###################################################
#parse_hmmscan_tbl: parse hmmscan domain table
###################################################
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from to  from    to  acc description of target
##DO SOMETHING CLEVER WITH SPLIT DOMAINS????
sub parse_hmmscan_tbl {
	my ($tbl) = @_;
	my %hits;
	open TBL, "<", $tbl;
	
	while(<TBL>){
		next if($_ =~ /^#/);
		chomp;
		my @splat = split /\s+/;
		#domain, eval, score, start, end
		push @{$hits{$splat[3]}}, [$splat[1], $splat[6], $splat[7], $splat[17], $splat[18]];
	}
	close TBL;
	return \%hits;
}

###################################################
#read_fasta: given a fasta file name, read it in to
#a hash with ID => sequence
###################################################
sub read_fasta {
	my($fasta) = @_;
	my ($line, $name, $seq);
	my %f;

	open FA, "<$fasta";
	$line = readline(FA);
	$line =~ /^>(.+)\n/;
	$name = $1;
	while(<FA>){
		if ($_ =~ /^>(.+)\n/){
			$f{$name} = $seq;
			$name = $1;
			$seq = '';
			next;
		}
		chomp;
		$seq .= $_;
	}
	$f{$name} = $seq;
	close FA;

	return(\%f);
}


###################################################
#find_seq_lengths: given a fasta file name, read it in to
#a hash with ID => seq length
###################################################
sub find_seq_lengths {
    
	my($fasta) = @_;
	my ($name);
	my %l;
	
	my $fh = Bio::DB::Fasta->newFh($fasta);
  	while (my $seq = <$fh>) {
   		$l{$seq->display_id()} = $seq->length;
	#	print "find_seq_lengths:\t[".$seq->display_id()."]\t[".$seq->length."]\n"; 
  	}

	return(\%l);	
}

###################################################
#write_CDS_fasta: given an embl file, dump fasta
#file containing all coding annotations; return
#reference to hash linking identifiers to gene 
#names and products
###################################################
sub write_CDS_fasta {
	my($embl, $format, $outfile) = @_;
	my %gene_info;
	
	my $cds_coordinates  = cds_locations($embl, $format);
    	my $annotation_file =  Bio::SeqIO->new(-file => $embl, -format => $format) or die "Error: Couldnt open $format file: $!\n";
	#my $annotation_file =  Bio::SeqIO->new(-file => $embl, -format => 'GENBANK') or die "Error: Couldnt open GENBANK file: $!\n";

	open OUT, ">$outfile" or die "Couldn't open file for writing: $!\n";

	while (my $sequence_annotation = $annotation_file->next_seq()) 
	{
  		for my $feature ($sequence_annotation->get_SeqFeatures())
  		{
    			next if !($feature->primary_tag eq 'CDS' || $feature->primary_tag eq 'polypeptide');
			next if $feature->has_tag("pseudo");
    			my $feature_id    = get_feature_id($feature);
    			my $gene_name     = get_gene_name($feature);
    			my $product_value = get_product_value($feature);
			my $start = $feature->start;
			my $end = $feature->end;
			my $sequence = '';

			@{$gene_info{$feature_id}} = ($gene_name, $product_value);
			if($feature->has_tag("translation")){
				my @value = $feature->get_tag_values("translation");
				$sequence = $value[0];
			}
			elsif($feature->strand == 1){
				$sequence = translate_as_string($feature->entire_seq->subseq($start, $end), -complete=>1,
											    -codontable_id=>11);	
			} else {
				$sequence = translate_as_string(revcom($feature->entire_seq->subseq($start, $end))->seq(), -complete=>1,
													   -codontable_id=>11);
			}

			print OUT ">$feature_id\n";
			print OUT substr($sequence, 0, 80, '')."\n" while (length($sequence));
		}
	}
	close OUT;
	return(\%gene_info)
}
 

###################################################
#is_gene_within_cds
###################################################
sub is_gene_within_cds
{
   my($cds_coordinates, $gene_feature) = @_;
   for my $current_coords(@{$cds_coordinates})
   {
     next if( $current_coords->[0] > $gene_feature->start);
     next if( $current_coords->[1] < $gene_feature->end);
     return 1;
   }
   
   return 0;
}
###################################################
#get_product_value: Get gene product
###################################################
sub get_product_value
{
  my($feature) = @_;
  my $product = "";
  my @junk;
  if($feature->has_tag('product'))
  {
    ($product, @junk) = $feature->get_tag_values('product');
  }
  
  return $product;
}

###################################################
#get_gene_name: Get gene name of genomic feature
###################################################
sub get_gene_name
{
  my($feature) = @_;
  my $gene_name;
  my @junk;
  if($feature->has_tag('gene'))
  {
    ($gene_name, @junk) = $feature->get_tag_values('gene');
  }
  else
  {
    $gene_name = get_feature_id($feature);
  }
  $gene_name =~ s/\W//g;
  return $gene_name;
}

###################################################
#get_feature_id: Get id of genomic feature
###################################################
sub get_feature_id
{
  my($feature) = @_;
  my $feature_id = int(rand(10000));
  my @junk;
  if($feature->has_tag('locus_tag'))
  {
    ($feature_id, @junk) = $feature->get_tag_values('locus_tag');
  }
  elsif($feature->has_tag('ID'))
  {
    ($feature_id, @junk) = $feature->get_tag_values('ID');
  }
  elsif($feature->has_tag('systematic_id'))
  {
    ($feature_id, @junk) = $feature->get_tag_values('systematic_id');
  }
  else
  {
    $feature_id = join("_",($feature->seq_id(), $feature->strand, $feature->start, $feature->end ));
  }
  $feature_id =~ s/^"|"$//g;
  return $feature_id ;
}

###################################################
#cds_locations: determine locations of "CDS"
#annotations in embl file
###################################################
sub cds_locations
{
  my($embl_file, $format) = @_;
  my @cds_coordinates;
  
  my $annotation_file =  Bio::SeqIO->new(-file => $embl_file, -format => $format) or die "Error: Couldnt open the annotation file\n";
#  my $annotation_file =  Bio::SeqIO->new(-file => $embl_file, -format => 'GENBANK') or die "Error: Couldnt open the annotation file\n";
  while (my $sequence_annotation = $annotation_file->next_seq()) 
  {
    for my $feature ($sequence_annotation->get_SeqFeatures())
    {
      next if !($feature->primary_tag eq 'CDS');
      push(@cds_coordinates, [$feature->start,$feature->end]);
    }
  }
  return \@cds_coordinates;
}

######################################################################
#Max and Min
#max
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

#min
sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}
######################################################################
#Max and Min for arrays:
#max
sub maxA {
    my $max = $_[0];
    foreach my $a (@_){
	$max = max($max, $a) if isNumeric($a);
    }
    return $max;
}

#min
sub minA {
    my $min = $_[0];
    foreach my $a (@_){
	$min = min($min, $a) if isNumeric($a);
    }
    return $min;
}

######################################################################
sub isNumeric {
    my $num = shift;
    if ($num=~/^-?\d+\.?\d*$/) { 
	return 1; 
    }
    else {
	return 0;
    }
}

###################################################
#help
###################################################
sub help {
	print STDERR <<EOF;
####################################################################
DeltaBS: monitoring functional changes in protein domains using HMMs
Version 0.1 5/4/2013; Author: Lars Barquist
####################################################################
Usage: $0 -f <filetype> -f1 <reference file> -f2 <comparator file> -o <output directory> -hp <path to hmmer> -hd <path to hmm libraries> -t <temp directory> [-options]

Options:
	-h  / --help		:	This screen
	-f  / --filetype	:	Specify filetype, must be 'embl', 'genbank', or 'fasta'
	-f1 / --file1		:	Reference genome/proteome in filetype format
	-f2 / --file2		:	Comparator genome/proteome in filetype format
	-pa1/ --pfamannot1      :       Pfam annotations of proteome1
	-pa2/ --pfamannot2      :       Pfam annotations of proteome2
	-ph1/ --phmmerannot1    :       phmmer domtblout of proteome1 vs proteome2
	-ph2/ --phmmerannot2    :       phmmer domtblout of proteome2 vs proteome1
	-o  / --outdir		:	Output directory
	-ol / --orthlist        :       Ortholog list file
	-hp / --hmmerpath	:	Path to hmmer installation
	-hd / --hmmlibdir	:	Path to hmm libraries (and annotation files for post-processing)
	-c  / --cpus            :       Number of CPUs for hmmsearch, phmmer etc to use. 
	-t  / --tempdir		:	Path to temporary directory
	-p  / --post		:	Enable post-processing (pathways, etc. EXPERIMENTAL)
	-d  / --dirty           :       Do not delete /tmp file
	-v  / --verbose		:	Turn on verbose messaging

EOF
}
