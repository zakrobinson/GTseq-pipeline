#!/usr/bin/perl
# V5 update: (Z.Robinson; October 2022) Made the regular expression more general to tolerate differences in plate names and accurately generate sample name in output. 
# V5 update: Process additional flag in genos file that specifies marker type and correctly calculate OT percentage. Also adds ability to print species diagnostics. 
# Modification (Kim Vertacnik; May 2022) of CRITFC_GTseq_GenoCompile_v2_IFI_colon.pl (Nate Campbell)
# v4 update: Getopt arguments; combined GTseq_GenoCompile_v2 with CRITFC_GTseq_GenoCompile_v2
# Usage: GTseq_GenoCompile_v4.pl <path to directory with .genos files> --name=['index' or 'sample'] --output=['count' or 'genotype'] --filter=[0-100] > libraryName.csv
# Compiles genotypes for multiple samples into a single .csv file
#  Sample-specific output fields: sampleName,onTargetReadCount,onTargetPercent,genotypingPercentage,IFIscore,markerGenotypes
#   IFI (individual fuzziness index): measure of DNA cross contamination calculated from background signal (read counts from homozygous and no-call loci)
#    Low scores are better than high scores
#    See GTseq_Genotyper_v4
# Input is .genos files from GTseq_Genotyper_v4
#  Files must be in the same directory
# '--name' output formats
#  index: i131_79_TOSS1107_OtsPBT20-PRHFa-1624
#  sample: OtsPBT20-PRHFa-1624
# '--output' formats
#  count: A1reads + A2reads
#  genotype: A1A2 bases
# '--filter' is the genotyping percentage cutoff (threshhold) value
#   Genotyping percentage: for each sample, how many panel markers were successfully genotyped
#   If a sample's genotyping percentage is below the cutoff, all its genotypes are converted to "0:0"

use strict; 
use warnings;
use Getopt::Long qw(GetOptions);

die("Usage: $0 <path to directory with .genos files> --name=['index' or 'sample'] --output=['count' or 'genotype'] --filter=[0-100] > outputFile.csv\n") unless @ARGV == 4;

my $flag;
my $genoThresh = 0;
my $nameFormat;

# Assign current working directory
my $cwdPath = shift @ARGV;

# Get and validate input arguments
GetOptions(
	'filter=i' => \$genoThresh,
	'name=s'   => \$nameFormat,
	'output=s' => \$flag,
) or die("Bad arguments\nUsage: $0 <path to directory with .genos files> --name=['index' or 'sample'] --output=['count' or 'genotype'] --filter=[0-100]\n");

if ($genoThresh < 0 || $genoThresh > 100) {die("Genotype percentage filter out of range: --filter=[0-100]\n")}

$nameFormat eq 'index' || $nameFormat eq 'sample' or die("Sample naming format: --name=['index' or 'sample']\n");

$flag eq 'count' || $flag eq 'genotype' or die("Genotype format: --output=['count' or 'genotype']\n");

# Get all input file names (all the *.genos files in the current working directory)
my @FileList = glob("$cwdPath/*.genos");
chomp ( @FileList );

print "Sample,Raw_Reads,On-Target_Reads,%On-Target,%GT,IFI";    #Beginning of the header line

# Open the first input .genos file and appends the marker names (first item of each row) to the header line
# Assumes that all the .genos files in the current working directory were genotyped with the same species panel
open (FILE1, "<$FileList[0]") or die("Could not open $!");
while (<FILE1>) {
	if ($. > 1){
		my @marker = split(",", $_);
		my $assay1 = $marker[0];
		print ",$assay1";
	}
}
print "\n";
close FILE1;

# Add each .genos file's summary and genotype data
foreach my $samples (@FileList) {
	my $rawRead_count  = 0;
	my $onTarget_count = 0;
	my $GT_pct         = 0;
	my $noGT_count     = 0;	
	my $marker_count   = 0;
	my $IFI            = 0;
	my $sample_name    = $samples;
	
	
	# Derive the sample name from the input filename
	$sample_name =~ s/$cwdPath\///;
	$sample_name =~ s/.genos//;

	if ($nameFormat eq 'index') {
		print "$sample_name,";
	}
	elsif ($nameFormat eq 'sample') {
		$sample_name =~ s/^i[0-9][0-9][0-9]_[0-9][0-9]_([A-Z]|[a-z])*[0-9]*[A-Z]?_//;
		print "$sample_name,";
	}
	
	# For each sample (.genos) file
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		# From the header line get the raw read count and IFI score
		if ($. == 1) {
			chomp;
			my @summary = split(",", $_);    
			$summary[1] =~ s/.+\://;
			print "$summary[1],";
			$rawRead_count = $summary[1];
			$IFI = $summary[4];
			$IFI =~ s/.+\://;
		}
		# From subsequent lines get the A1 and A2 read counts for each marker
		elsif ($. > 1) {
			chomp;
			# Get a1, a2, and total on-target read counts
			my @info  = split(",", $_);    # marker,A1=count,A2=count,ratio,genotype
			my $diagFlag; # diagFlag ensures that secondary SNPs on an amplicon do not inflate OT%
			if(exists $info[13]){$diagFlag = $info[13];} # If marker flag isn't present treat markers as primary
			else				{$diagFlag = 1;}
			
							
			unless (($diagFlag > 0) && ($diagFlag % 2 != 1)){ # Do not count any markers that have the diagFlag == 2 or 4. These "non-primary" variant targets on an amplicon that has a primary target (diagFlag==1[polymorphic] or 3[species]). 
					$marker_count++; # At present, the exclusion of diagFlag==2 or 4 makes this an amplicon count ZLR. 
					my $geno  = $info[4];
					if ($geno =~ m/NA|00/) {$noGT_count++}    # Count of how many markers with no genotype call : Only Considers diagFlag= 0[dominant sex marker],1[primary snp],3[primary species]
            		my $OTcount = $info[8];
					$onTarget_count = $onTarget_count + $OTcount;  # Cumulative for all markers
					
			}
		}
	}
	close FILE;

    # Calculate the sample's summary stats	
	my $OT_pct;
	
	
	if ($rawRead_count == 0) {$OT_pct = 0.00;}
	else                     {$OT_pct = $onTarget_count / $rawRead_count * 100;}

	$OT_pct          = sprintf("%.2f", $OT_pct);
	$GT_pct          = 100 - ( $noGT_count / $marker_count * 100 );
	my $Print_GT_pct = sprintf("%.2f", $GT_pct);
	
	print "$onTarget_count,$OT_pct,$Print_GT_pct,$IFI,";    # Appended to 'print "$sample_name,";'
	
    # Get the genotypes for each marker	
	open (FILE, "<$samples") or die;
	while (<FILE>) {
		if ($. > 1) {
		chomp;
		my @marker = split(",", $_);    # marker,A1=count,A2=count,ratio,genotype
		my @Gbases = split("", $marker[4]);    # split A1A2 genotype into A1 and A2 alleles
		my $geno = "$Gbases[0]:$Gbases[1]";
		my $alleleRead_count = 0;
		my $diagFlag;
		
		if(exists $marker[13]){$diagFlag = $marker[13];}
		else				  {$diagFlag = 1;}
		
		$marker[1] =~ s/.*=//;
		$marker[2] =~ s/.*=//;
		$alleleRead_count = $marker[1] + $marker[2];
		
        # If the sample's %GT is above the filtering threshhold, print genotypes in desired format
        # ZLR added the diagFlag, which will renders species markers (diaFlags= 3[primary] or 4[secondary]) are always printed
        # if no longer desired Species markers can be made subject to GT filter by changing 3 and 4 to 1 and 2. Without effecting accuracy of OT%
		# If below print 0:0 for all markers
		
		if(($flag eq 'genotype') && ($GT_pct >= $genoThresh || $diagFlag > 2 )) {print "$geno,";}
			elsif(($flag eq 'count') && ($GT_pct >= $genoThresh || $diagFlag > 2)) {print "$alleleRead_count,";}
			else {print "0:0,";}
		} 
	}
	print "\n";
	close FILE;
}

exit(0);
