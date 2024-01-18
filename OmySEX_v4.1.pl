#!/usr/bin/perl

#OmySEX_V4.1.pl
#V4.1 made the cntrl_ratio an argument to the function.
#V4 update added a flag to end of the .genos concatenation and made output an equal length to standard .genos row entries
#V4 update parallelized and takes argument for number of jobs/processors to use. 
# OmySEX_test_v3.pl
# by Nate 
# Sex fish by GTseq
# This version of the script uses information from the version 2 genotyping output as a control for the sex marker.
# It then generates a genotypic sex and appends the sex marker to the end of the .genos file.
# The sex marker is then included in the summary file generated by the GTseq_GenoCompile_v2.pl script.
# update 5/2016: Incorporates fuzzy match to male specific probe allowing for a single base difference in probe.  Makes
# sex marker data compatible with GTseq_ErrorReport.pl script.

use strict; 
use warnings;
use Parallel::ForkManager;
use String::Approx 'amatch';
use PerlIO::gzip; #Kim added gzip reading capability
use Getopt::Long qw(GetOptions);

die("Usage: $0 --jobs=[integer: number of processors] --ctrl_ratio=[real number : ratio X OTreads = X-counts]\nExecute in directory with .genos and .fastq files\n") unless @ARGV == 2;

my $ctrl_ratio;
my $max_processors;

# Get and validate input arguments
GetOptions(
	'jobs=i' => \$max_processors,
	'ctrl_ratio=f'=>\$ctrl_ratio,
) or die("Bad arguments\nUsage: $0 --jobs=[integer: number of processors] --ctrl_ratio=[real number : ratio X OTreads = X-counts]\nExecute in directory with .genos and .fastq files\n");




sub ProcessFastqXY{
	my $samples = $_;
	open (FILE, "<:gzip", "$samples") or die;	
	$samples =~ s/.fastq.gz//;
	my $genos = "$samples.genos";
	my $OT_reads = 0;
	my $primer_counts = 0;
	my $primerOT = 0;
	my $perofallOTreads = 0;
	open (READ, "<$genos") or die;
		while (<READ>) {
			if ($. == 1) {
				my @info = split ",", $_;
				my @info2 = split ":", $info[2];
				$OT_reads = $info2[1];
					}
				}
	close READ;
	open (OUT, '>>', $genos) or die "Error opening $genos\n";
	my $counts = 0;
	my $cntrl_counts = 0;
	my $sex_geno = "00";
	my $geno_class = "NA";
	while (<FILE>) {
		chomp;
		my $info_line = $_;
		my $seq_line = <FILE>;
		my $info_line2 = <FILE>;
		my $qual_line = <FILE>;
		if($seq_line =~ m/^GCGCATTTGTATGGTGAAAA/) {$primer_counts++}
		if(($seq_line =~ m/^GCGCATTTGTATGGTGAAAA/) && (amatch("ATGTGTTCATATGCCAG", ["S1"], $seq_line))) {$counts++}
	}
	if($primer_counts == 0) {$primer_counts = 1}
	$primerOT = $counts/$primer_counts*100;
	$primerOT = sprintf("%.3f", $primerOT);
	#Kim added if/else statement to deal with divide by zero error
	if ($OT_reads == 0) {$perofallOTreads = 0}
	else {$perofallOTreads = $counts/$OT_reads*100}
	$perofallOTreads = sprintf("%.3f", $perofallOTreads);
	$cntrl_counts = $OT_reads * $ctrl_ratio; # ctrl_ratio = 0.001 as of 01/2024
	$cntrl_counts = int ( $cntrl_counts );
	if ($cntrl_counts == 0) {$cntrl_counts = 1}
	my $real_counts = $counts; # added by Zak to correctly count OT reads
	if ($counts == 0) {$counts = 1}
	my $ratio = $cntrl_counts/$counts;
	$ratio = sprintf("%.3f", $ratio);
	if ($cntrl_counts + $counts < 10) {$sex_geno = "00", $geno_class = "NA"}
	elsif ($ratio >= 10) {$sex_geno = "XX", $geno_class = "A1HOM"}
	elsif ($ratio <= 0.1) {$sex_geno = "XY", $geno_class = "A2HOM"}
	elsif ($ratio <= 0.2) {$sex_geno = "00", $geno_class = "NA"}
	elsif ($ratio <= 5) {$sex_geno = "XY", $geno_class = "HET"}
	
	print OUT "OmyY1_2SEXY,X=$cntrl_counts,Y=$counts,$ratio,$sex_geno,$geno_class,0,0,$real_counts,$primerOT,$perofallOTreads,0.00,Pass,0\n";
close FILE;
close OUT;
}



my @Files = `ls *.fastq.gz`;
chomp ( @Files );


my $fork= new Parallel::ForkManager($max_processors);


foreach (@Files) {
$fork->start and next; # do the fork
ProcessFastqXY($_);
$fork->finish; # do the exit in the child process
}
$fork->wait_all_children;


