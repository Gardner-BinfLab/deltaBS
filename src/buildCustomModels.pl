#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use LWP::Simple;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

use Getopt::Long;

my ($proteome, $database, $testMode, $verbose, $help);
my $datadir = ".";
&GetOptions(
        "t|testmode"         => \$testMode, 
	"p|proteome=s"       => \$proteome,
	"d|datadir=s"        => \$datadir,
        "db|database=s"      => \$database,
	"v|verbose"	     => \$verbose,
	"h|help"	     => \$help
	);

if ($help) {
	&help();
	exit(1);
}

if (not defined $proteome and not defined $testMode) {
    print "PROTEOME FILE NOT DEFINED! Use the [-p] option.\n";
    &help();
    exit(1);
}

#fetch uniref90 if it's not there
if ( not defined($database) and ((not -s "$datadir/uniref90.fasta.gz") or (not -s "$datadir/uniref90.fasta")) ){
    
    print "fetching [ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz], writing to [$datadir/uniref90.fasta.gz]\n" if (defined $verbose);
    
    open OUT, "> $datadir/uniref90.fasta.gz";
    select OUT; $| = 1;
    my $cnt=0;
    while(!is_success(getprint("ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"))){
	print STDERR "\t request failed, retrying\n";
	sleep(10);
	die "ERROR: failed to fetch uniref90.fasta.gz!" if ($cnt > 20); 
	$cnt++;
    }
}
else {
    print "Skipping fetch of uniref90...\n" if (defined $verbose);
}

#unzip uniref90 file:
if( not defined($database) and (-s "$datadir/uniref90.fasta.gz" and not -s "$datadir/uniref90.fasta")){
    print "unzipping [$datadir/uniref90.fasta.gz]\n" if (defined $verbose);
    my $status = gunzip "$datadir/uniref90.fasta.gz" => "$datadir/uniref90.fasta" 
        or die "gunzip failed: $GunzipError\n";
    $database = "$datadir/uniref90.fasta";
}
elsif(not defined($database)){
    $database = "$datadir/uniref90.fasta";
}

#generate a small test DB and sample test proteome 
if(defined $testMode){
    print "Generating test database [$datadir/buildCustomModels-test-proteome]\n" if (defined $verbose);
    $proteome = "$datadir/buildCustomModels-test-proteome" if (not defined($proteome));
    ($database,$proteome) = generateTestDatasets($database,$proteome);
    print "Generated test data [($database,$proteome)]\n\n\n" if (defined $verbose);
}

######################################################################
#run jackhmmer searches:
#split uniref90.fasta and paralellise if you want to run on a cluster:
my $jackOut  = "$datadir/jackhmmer-alignments.stk";
my $jackOut2 = "$datadir/jackhmmer-alignments.jackhmmer";

my $jackExe = "jackhmmer -N 1 -A $jackOut $proteome $database > $jackOut2";
print "Running jackhmmer! [$jackExe]\n" if (defined $verbose);
system("$jackExe") and 
    die "FATAL: failed to execute [$jackExe]\n[$!]";

######################################################################
#Filter sequences that are less than 40% PID similar to reference
print "Running [filterDivSeqs($jackOut,40)]!\n" if (defined $verbose);
my $filteredAlignment = filterDivSeqs($jackOut,40); 

######################################################################
#build new hmms:
my $buildExe = "hmmbuild $filteredAlignment\.hmm $filteredAlignment";
print "Running hmmbuild! [$buildExe]\n" if (defined $verbose);
system("$buildExe") and 
    die "FATAL: failed to execute [$buildExe]\n[$!]";

print "New custom models are available from: [$filteredAlignment\.hmm]!\n" if defined($verbose);

exit(0);
######################################################################



sub generateTestDatasets {
    ($database,$proteome) = @_;
    
    my ($dbSize, $proteomeSize)  = (0,0);
    
    my $printMe = 0;
    open(DBIN, "< $database") or die "FATAL: failed to open [$database] for reading!\n[$!]";
    open(DBUT, "> " . $database . ".testdb") or die "FATAL: failed to open [$database\.testdb] for reading!\n[$!]";
    open(PTUT, "> " . $proteome . ".test"  ) or die "FATAL: failed to open [$proteome\.test] for reading!\n[$!]";
    
    while(my $dbin = <DBIN>){
	if($dbin =~ /^>.*Chorismate synthase/ or $dbin =~ /^>.*Pyruvate carboxylase/){
	    $printMe=1;
	    $dbSize++;
	}
	elsif($dbin =~ /^>/){
	    $printMe=0;
	}
	
	print DBUT $dbin if ($printMe);
	
	if($proteomeSize < 10 && $printMe){
	    print PTUT $dbin;
	    $proteomeSize++;
	}
	
	last if ($dbSize > 1000 && $printMe == 0);
    }
    
    close(PTUT);
    close(DBUT);
    close(DBIN);
    
    $database = $database . ".testdb";
    $proteome = $proteome . ".test";
    
    return ($database,$proteome);
}


######################################################################
#filterDivSeqs: Filter sequences that are less than  $pidThresh% PID similar to reference
sub filterDivSeqs {
    my ($inFile, $pidThresh) = @_;
    
    #MAY NEED TO SPLIT ALIGNMENTS INTO SEPERATE FILES:

    #1. fetch alignment ID's 
    my $alnIDs = `grep ^"#=GF ID" $inFile`; 
#=GF ID UniRef90_Q81FQ3-i1
    my @alnIDsTemp = split(/\n/, $alnIDs); 
    my @alnIDs     = ();
    foreach my $alnID (@alnIDsTemp){
	if($alnID =~ /^#=GF ID\s+(\S+)/){
	    push(@alnIDs, $1);
	}
    }
    
    my $outFile = $inFile . ".filtered." . $pidThresh . ".stk";    
    open(UT, "> $outFile") or die "FATAL: failed to open [$outFile]\n[$!]";
    #2. fetch alignments 
    foreach my $alnID (@alnIDs){
	
	my $fetchExe = "esl-afetch $inFile $alnID | egrep -v \'^#=GS|^#=GC|^#=GR\' > /tmp/$$.tempaln && esl-reformat pfam /tmp/$$.tempaln";
	print "Running: [$fetchExe]\n\n" if (defined($verbose));
	open(AIN, "$fetchExe |" ) or die "FATAL: unable to open pipe for [$fetchExe].\n[$!]";
	
	my $originalQuery = "";
	my %seqsToFilter  = ();
	my %seqsToPadding = ();
	
	print UT "# STOCKHOLM 1.0\n";
	while(my $ain = <AIN>){
	    if ($ain =~ /^#=GF\s+ID\s+(\S+)-i1/){
		$originalQuery = $1;
		print UT "#=GF ID newfams-$originalQuery\-filtered$pidThresh\n\n";
		#print    "Original query:[$originalQuery]\n[#=GF ID newfams-$originalQuery\-filtered$pidThresh]\n\n\n";
	    }
	    elsif($ain =~ /^(\S+)(\s+)(\S+)$/){
		$seqsToPadding{$1}=$2;
		$seqsToFilter{$1} =$3;
	    }
	}
	
	if(-s "/tmp/$$.tempaln" && not defined($verbose)){
	    system("rm /tmp/$$.tempaln");
	}
	
	
	if( not defined $seqsToFilter{$originalQuery} ){
	    print STDERR "WARNING: failed to find originalQuery:[$originalQuery] in alignment [$alnIDs] in file [$inFile]!";
	    print UT "//\n";
	    next; 
	}
	
	#print "org:[$seqsToFilter{$originalQuery}]\n";
	print UT $originalQuery . $seqsToPadding{$originalQuery} . $seqsToFilter{$originalQuery} . "\n";

	foreach my $seqId (keys %seqsToFilter){
	    
	    next if ($seqId eq $originalQuery);
	    my $pid = computePID($seqsToFilter{$originalQuery}, $seqsToFilter{$seqId} );
	    
	    #print "[$pid]\n";
	    #print "[$seqId] new:[$seqsToFilter{$seqId}] [$pid < $pidThresh]\n";
	    if($pid < $pidThresh){
		#print "deleting [$seqId]\n";
		delete $seqsToFilter{$seqId};
	    } else {
		#print "printing [$seqId]\n";
		print UT $seqId . $seqsToPadding{$seqId} . $seqsToFilter{$seqId} . "\n";
	    }
	    
	}
	print UT "//\n";
    }
    close(UT);

    return $outFile;
}

######################################################################
#computePID: compute the identity between two sequences.
sub computePID {
    my $a = shift;
    my $b = shift;
    my @a = split(//, $a);
    my @b = split(//, $b);
    my ($sim, $len) = (0, 0);
    
    if (scalar(@a) != scalar(@b)){
	return 0;
    }
    else {
	
 	for (my $i=0; $i<scalar(@b); $i++){
	    $a[$i] = uc($a[$i]);
	    $b[$i] = uc($b[$i]);

	    if ( (isAA($a[$i]) || isAA($b[$i])) ){
		$len++; #don't count gap-gaps
		if($a[$i] eq $b[$i]){
		    $sim++;
		}
	    }
	}

    }
    
    if ($len>0){
	return 100*$sim/$len;
    }
    else {
	return 0;
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


######################################################################
#help
###################################################
sub help {
	print STDERR <<EOF;
####################################################################
buildCustomModels.pl: build custom profile HMMs for a proteome using UniRef90 as a source DB 
    Version 0.1 2016/04/22; Author: Nicole Wheeler, Lars Barquist, Paul Gardner
####################################################################
Usage: $0 -d <data directory for writing> -p <proteome file in fasta format> [-options]

Options:
	
	-d  / --datadir  <s>    :       Data directory for writing output to
	-db / --database <s>    :       Database file default:[datadir/uniref90.fasta]
	-p  / --proteome <s>	:	A proteome file in fasta format 
	-t  / --testmode	:	Generate small database and proteome files for testing 
	-v  / --verbose		:	Turn on verbose messaging
	-h  / --help		:	This screen

EOF
}




