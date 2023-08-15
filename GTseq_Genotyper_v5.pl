#!/usr/bin/perl



# Modification (Zak Robinson; October 2022) of GTseq_Genotyper_v4.pl (Kim Vertacnik; Nate Campbell)
# v5 update: Address total On-target percentages greater than 100%. Process each line of fastq as it is read to avoid double counting reads. Multidimensional hash addresses multiple snps per amplicon. 
# v4 update: new structure for microhaps; removes fuzzy match sequence printing
# Usage: GTseq_Genotyper_v4.pl probeSeq.csv sample.fastq.gz
# Genotypes .fastq file using GTseq probes. Outputs .genos file with genotype calls for each GTseq marker
#  Input undemultiplexed fastq file must be gzip compressed
# probeSeq.csv lists [Locus Name,Allele1,Allele2,ProbeSeq1,ProbeSeq2,FWD_Primer,A1_correction,A2_correction] for all loci
# Output files start with a summary header line followed by locus-specific data 
#  Locus-specific output fields: LocusName, Allele1_counts, Allele2_counts, A1/A2-ratio, Genotype, Genotype_Class, A1_Correction_value, A2_correction_value, On_target (probe match) reads, Locus OT_percentage (probe match of all reads with fwd primer), Locus On-target reads as percentage of total on-target reads,IFI score,null allele pass/fail, ProbeSeq Flag
#  Locus OT_percentage is defined as: #reads beginning with forward primer sequence and containing a probe sequence / #All reads beginning with forward primer sequence * 100
#  Locus On-target reads is defined as: total number of on-target reads for locus / Overall on-target reads in panel * 100
#  IFI (individual fuzziness index) score: measure of DNA cross contamination calculated from background signal (read counts from homozygous and no-call loci)
#   Low scores are better than high scores
#  Pass/fail check for possible null alleles. Flagged null alleles could also be from the co-amplification of a paralogous locus

use strict;
use warnings;
use PerlIO::gzip;
use String::Approx 'amatch';

die "Usage: $0 probeSeq.csv sample.fastq.gz\n" unless @ARGV == 2;

# Primer and probe data preliminaries
my %probedata;
my %refProbeData;
my @info;
# Allele and read counting preliminaries
my $onTarget_total_count = 0;
my $rawRead_count        = 0;
my %onTarget_count;
my %a1_probe_count;
my %a2_probe_count;
my %offTarget_count;
my %fuzzyMatchCount;
my %probeLength;
my %probe_on_read_start;
# Genotyping Preliminaries
my %nullTest;
my %printLine;
my $background_count     = 0;
my $arraySize            = 0;
my $homozygous_count     = 0;
my $ifi                  = 0;


# Process probeSeq file and generate multidimensional hash for read/allele counting.
open(PROBES, "<$ARGV[0]") or die "error reading $ARGV[0]\n";
while (<PROBES>) {
    chomp;
    @info = split(/,/, $_);
    $info[3] =~ s/R|Y|S|W|K|M|B|D|H|V|N/./g; # This is suboptimal replacing IUPAC with generic wildcard
    $info[4] =~ s/R|Y|S|W|K|M|B|D|H|V|N/./g; # This is suboptimal replacing IUPAC with generic wildcard
	$info[5] = substr $info[5], 0, 14;
	# Add Correction factors and Primary target flag if absent from files
	unless (exists $info[6]) {push(@info, 0);}
	unless (exists $info[7]) {push(@info, 0);} #Assumes no corrections if none provided.
	unless (exists $info[8]) {push(@info, 1);} # Assumes an unflagged marker is primary target
	# Add reverse complements of Probseqs
	my $revcomp1;
    my $revcomp2;
    $revcomp1 = reverse $info[3];
    $revcomp1 =~ tr/[ACGT]/[TGCA]/;
    $revcomp1 =~ tr/\]\[/\[\]/;
    $revcomp2 = reverse $info[4];
    $revcomp2 =~ tr/[ACGT]/[TGCA]/;
    $revcomp2 =~ tr/\]\[/\[\]/;
    push(@info,$revcomp1,$revcomp2); # Add's RCprobes to @info array at position 9, 10 
	#populate multidimensional hash
	$probedata{$info[5]}{$info[0]} = [@info]; # Probe data stored in multidimensional hashes in the following format data{FWD_primer}{Locus_Name}=> probeSeq data array
	#Basic hash of loci and details
	$refProbeData{$info[0]}=[@info]; #Basic hash for data{locus} => probeSeq data; Convenient for genotyping step with minor memory cost.
	#Initializing locus specific data hashes
	$a1_probe_count{$info[0]}  = 0;
	$a2_probe_count{$info[0]}  = 0;
	$fuzzyMatchCount{$info[0]} = 0;
	$onTarget_count{$info[0]}  = 0;
	$offTarget_count{$info[0]} = 0;
	$probeLength{$info[0]} = 20;    # I don't know why the defaults are 20 and not 0
	$probe_on_read_start{$info[0]} = 20; 
	$nullTest{$info[0]} = "Pass";	
}
close PROBES;

# Process fastq.gz file and count reads, alleles, fuzzy matches, etc. 

my $fastq_start = -1;    # (-1) Because the sequence line is line 2 of 4 and we want modulus==1
open FASTQ, "<:gzip", $ARGV[1] or die "error reading $ARGV[1]\n";
while (<FASTQ>) {
    chomp;
    if (++$fastq_start % 4 == 1) {
     $rawRead_count++;
     my $R1_seq = $_;
     my $FP_seq = substr $R1_seq, 0, 14;
     my $read_OTcounted = 0;
     if (exists $probedata{$FP_seq}){ # Does the read have FWD primer sequence first 15 BP
     	for my $locus (sort keys %{$probedata{$FP_seq}}){ #
    		if ($R1_seq =~ m/$probedata{$FP_seq}{$locus}[3]|$probedata{$FP_seq}{$locus}[9]/) {
                $probe_on_read_start{$locus} = $-[0];    # The read bp position where the probe starts aligning
              	$a1_probe_count{$locus}++;
               	$onTarget_count{$locus}++;
                $read_OTcounted++;
            }# if 
            elsif ($R1_seq =~ m/$probedata{$FP_seq}{$locus}[4]|$probedata{$FP_seq}{$locus}[10]/) {
             	$probe_on_read_start{$locus} = $-[0];
              	$a2_probe_count{$locus}++;
               	$onTarget_count{$locus}++;
              	$read_OTcounted++;
            }
            else {
                $offTarget_count{$locus}++;
                my $fuzz      = 0;
                my $a1_fuzzy  = $probedata{$FP_seq}{$locus}[3];
                my $subnumber = 0;

                #Set up amatch parameters
                my $fuzzy_set  = "[ \"I1\", \"S1\", \"D1\" ]";
                if ($a1_fuzzy  =~ m/\[|\.|\?|\+/) {
                    $a1_fuzzy  =~ s/\[..\]/N/g;
                    $a1_fuzzy  =~ tr/[.?+]/[NNN]/;
                    $subnumber = $a1_fuzzy =~ tr/N/N/ + 1;
                    $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";
                }
                my $a1_fuzzy_rc  = $probedata{$FP_seq}{$locus}[9];
                if ($a1_fuzzy_rc =~ m/\[|\.|\?|\+/) {
                    $a1_fuzzy_rc =~ s/\[..\]/N/g;
                    $a1_fuzzy_rc =~ tr/[.?+]/[NNN]/;
                    $subnumber   = $a1_fuzzy_rc =~ tr/N/N/ + 1;
                    $fuzzy_set   = "[ \"I1\", \"S$subnumber\", \"D1\" ]";
                }
                my $a2_fuzzy   = $probedata{$FP_seq}{$locus}[4];
                if ($a2_fuzzy  =~ m/\[|\.|\?|\+/) {
                    $a2_fuzzy  =~ s/\[..\]/N/g;
                    $a2_fuzzy  =~ tr/[.?+]/[NNN]/;
                    $subnumber = $a2_fuzzy =~ tr/N/N/ + 1;
                    $fuzzy_set = "[ \"I1\", \"S$subnumber\", \"D1\" ]";
                }
                my $a2_fuzzy_rc  = $probedata{$FP_seq}{$locus}[10];
                if ($a2_fuzzy_rc =~ m/\[|\.|\?|\+/) {
                    $a2_fuzzy_rc =~ s/\[..\]/N/g;
                    $a2_fuzzy_rc =~ tr/[.?+]/[NNN]/;
                    $subnumber   = $a2_fuzzy_rc =~ tr/N/N/ + 1;
                    $fuzzy_set   = "[ \"I1\", \"S$subnumber\", \"D1\" ]";
                }

                if (length $a1_fuzzy >= length $a2_fuzzy) {$probeLength{$locus} = length $a1_fuzzy}
				else                                      {$probeLength{$locus} = length $a2_fuzzy;}

                #Search for fuzzy matches
                if (amatch($a1_fuzzy, $fuzzy_set, $R1_seq)) {
                    my $fuzzy_str = substr $R1_seq, $probe_on_read_start{$locus}, $probeLength{$locus};
                    if (length $fuzzy_str > 10) {$fuzz++;}

                }
                elsif (amatch($a1_fuzzy_rc, $fuzzy_set, $R1_seq)) {
                    my $fuzzy_str = substr $R1_seq, $probe_on_read_start{$locus}, $probeLength{$locus};
					$fuzzy_str = reverse $fuzzy_str; $fuzzy_str =~ tr/ACGT/TGCA/;
					if (length $fuzzy_str > 10) {$fuzz++;}
                }
                elsif (amatch($a2_fuzzy, $fuzzy_set, $R1_seq)) {
                    my $fuzzy_str = substr $R1_seq, $probe_on_read_start{$locus}, $probeLength{$locus};
					if (length $fuzzy_str > 10) {$fuzz++;}
                }
                elsif (amatch($a2_fuzzy_rc, $fuzzy_set, $R1_seq)) {
                    my $fuzzy_str = substr $R1_seq, $probe_on_read_start{$locus}, $probeLength{$locus};
					$fuzzy_str = reverse $fuzzy_str; $fuzzy_str =~ tr/ACGT/TGCA/;
					if (length $fuzzy_str > 10) {$fuzz++;}
                }
                if ($fuzz > 0) {$fuzzyMatchCount{$locus}++}
            }
            
    	
    	}# close for loop
    	if ($read_OTcounted > 0) {$onTarget_total_count++;} # Addresses double counting on reads for total on target; If the read was counted as OT for any probe for the given amplicon (FWD primer) it's OT. 
     }# close If has fwd primer    
}# close if it is the sequence line in fastq
 } # close while loop    
close FASTQ;

# For each marker, process allele counts and call genotypes

foreach my $loci (sort keys %refProbeData){
    my $A1fix  = 0;
    my $A2fix  = 0;
    my $sum_xy = $a1_probe_count{$loci} + $a2_probe_count{$loci};
    # per marker IFI prelims
    my $marker_homo_ct  = 0;
    my $marker_bkgrd_ct = 0;
    my $marker_ifi      = 0;
    #
    $a1_probe_count{$loci} =  $a1_probe_count{$loci} - ($sum_xy / 4 * $refProbeData{$loci}[6]);    # I don't know why 4
    if ($a1_probe_count{$loci} < 0) {$a1_probe_count{$loci} = 0}
    $a1_probe_count{$loci} = int($a1_probe_count{$loci});
    
    $a2_probe_count{$loci} = $a2_probe_count{$loci} - ($sum_xy / 4 * $refProbeData{$loci}[7]);
    if ($a2_probe_count{$loci} < 0) {$a2_probe_count{$loci} = 0}
    $a2_probe_count{$loci} = int($a2_probe_count{$loci});
    
    my $geno      = "00";
    my $genoclass = "NA";

    # Fix allele counts to non-zero number for ratio calculation
    if ($a1_probe_count{$loci} == 0) {$A1fix = 0.1}
    else                             {$A1fix = $a1_probe_count{$loci}}
    if ($a2_probe_count{$loci} == 0) {$A2fix = 0.1}
    else                             {$A2fix = $a2_probe_count{$loci}}
    
    my $ratio = $A1fix / $A2fix;
    $ratio = sprintf("%.3f", $ratio);

    # Genotype call
    if ($a1_probe_count{$loci} + $a2_probe_count{$loci} < 10) {
        $geno      = "00";                                               # Low allele count loci (less than 10 matching reads) get "00" genotype
        $genoclass = "NA";
    }
    elsif ($ratio >= 10) {
        $geno      = "$refProbeData{$loci}[1]$refProbeData{$loci}[1]";
        $genoclass = "A1HOM";                                            # Allele1 Homozygotes
        $homozygous_count = $homozygous_count + $a1_probe_count{$loci};
        $background_count = $background_count + $a2_probe_count{$loci};
        $marker_homo_ct = $marker_homo_ct + $a1_probe_count{$loci};
        $marker_bkgrd_ct = $marker_bkgrd_ct + $a2_probe_count{$loci};
    }
    elsif (($ratio < 10) && ($ratio > 2)) {
        $geno      = "00";
        $genoclass = "NA";                                               # In-betweeners
        $homozygous_count = $homozygous_count + $a1_probe_count{$loci};
        $background_count = $background_count + $a2_probe_count{$loci};
        $marker_homo_ct = $marker_homo_ct + $a1_probe_count{$loci};
        $marker_bkgrd_ct = $marker_bkgrd_ct + $a2_probe_count{$loci};
    }
    elsif ($ratio <= 0.1) {
        $geno      = "$refProbeData{$loci}[2]$refProbeData{$loci}[2]";
        $genoclass = "A2HOM";                                            # Allele2 Homozygotes
        $homozygous_count = $homozygous_count + $a2_probe_count{$loci};
        $background_count = $background_count + $a1_probe_count{$loci};
        $marker_homo_ct = $marker_homo_ct + $a2_probe_count{$loci};
        $marker_bkgrd_ct = $marker_bkgrd_ct + $a1_probe_count{$loci};
    }
    elsif ($ratio < 0.5) {
        $geno      = "00";
        $genoclass = "NA";                                               # In-betweeners
        $homozygous_count = $homozygous_count + $a2_probe_count{$loci};
        $background_count = $background_count + $a1_probe_count{$loci};
        $marker_homo_ct = $marker_homo_ct + $a2_probe_count{$loci};
        $marker_bkgrd_ct = $marker_bkgrd_ct + $a1_probe_count{$loci};
    }
    elsif ($ratio <= 2) {
        $geno      = "$refProbeData{$loci}[1]$refProbeData{$loci}[2]";
        $genoclass = "HET";                                              # Heterozygotes
    }
	

    
    # Null allele check
    if ($sum_xy == 0) {$sum_xy = 0.1};
	my $fuzz_ratio = $fuzzyMatchCount{$loci} / $sum_xy;
	if (($genoclass !~ "HET") && ($fuzz_ratio > 0.8) && ($sum_xy > 10)) {$nullTest{$loci} = "Fail"}
    
    # Calculate summary stats
    if ($marker_homo_ct == 0) {$marker_homo_ct = 1}
    $marker_ifi = $marker_bkgrd_ct / $marker_homo_ct * 100;
    $marker_ifi = sprintf("%.2f", $marker_ifi);

    if ($offTarget_count{$loci} == 0) {$offTarget_count{$loci} = 0.01}
    my $onTarget_of_markerFwdPrimer_per = ($onTarget_count{$loci} / ($offTarget_count{$loci} + $onTarget_count{$loci})) * 100;
    $onTarget_of_markerFwdPrimer_per = sprintf("%.1f", $onTarget_of_markerFwdPrimer_per);
    
    my $onTarget_of_allFwdPrimer_per = 0;
    if ($onTarget_total_count > 0) {$onTarget_of_allFwdPrimer_per= $onTarget_count{$loci} / $onTarget_total_count * 100;}
    $onTarget_of_allFwdPrimer_per= sprintf("%.3f", $onTarget_of_allFwdPrimer_per);

    $printLine{$loci} = "$loci,$refProbeData{$loci}[1]=$a1_probe_count{$loci},$refProbeData{$loci}[2]=$a2_probe_count{$loci},$ratio,$geno,$genoclass,$refProbeData{$loci}[6],$refProbeData{$loci}[7],$onTarget_count{$loci},$onTarget_of_markerFwdPrimer_per,$onTarget_of_allFwdPrimer_per,$marker_ifi,$nullTest{$loci},$refProbeData{$loci}[8]";
}

if ($onTarget_total_count == 0) {$onTarget_total_count = 1}
my $onTarget_percent = $onTarget_total_count / $rawRead_count * 100 unless $rawRead_count == 0;
if ($rawRead_count == 0) {$onTarget_percent = 0.0}
$onTarget_percent = sprintf("%.1f", $onTarget_percent);

if ($homozygous_count == 0) {$homozygous_count = 1}
$ifi = $background_count / $homozygous_count * 100;
$ifi = sprintf("%.2f", $ifi);

# Print summary line (line 1 of the .genos file)
print "$ARGV[1],Raw_reads:$rawRead_count,On-Target_reads:$onTarget_total_count,%On-Target:$onTarget_percent,IFI_score:$ifi\n";

# Print marker lines
foreach my $loci2 (sort keys %refProbeData) {print "$printLine{$loci2}\n";}

exit(0);

