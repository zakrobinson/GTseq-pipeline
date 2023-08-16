# GTseq-pipeline
A slightly modified version of Perl-based GTseq genotyping pipeline presented in Campbell et al. (2015)\
\
Originally authored by Nate Campbell and later modified by Kim Vertacnik and Zak Robinson. 

### Bug Fixes and Modifications
1) <i>GTseq_Genotyper_v5.pl</i> now allows for multiple snps for amplicons. At present, each SNP is called independently and phased haplotypes are not provided.
2) A typo resulting in asymetrical heterozygote calling based on the ratio of allelic counts has been corrected.
3) An additional column, termed a "probeSeq flag", has been added to probeSeq files. These flags allow for accurate counting of on-target reads in <i>GTseq_GenoCompile_v5.1.pl</i> and to omit certain loci from genotyping success \(%GT\) filters. If the probeSeq flag is omitted from file, all loci assume a value of 1.

|ProbeSeq Flag|Explanation|
|-----------|-----------|
|0|Sex marker with stand-alone script|
|1|Flagship SNP of amplicon|
|2|Secondary SNP|
|3|Flagship SNP of amplicon. Not subject to %GT filter|
|4|Secondary SNP. Not subject to %GT filter|

## Example Usage On A Linux System
### Notes on generating a raw FASTQ file:
<i>GTseq_Demultiplex.py</i> operates on a non-demultiplexed FASTQ file from an Illumina sequencing run.\
The provided, example file <i>Undetermined_S0_R1_001.fastq.gz</i> was generated using the following command using a generic sample sheet:

bcl-convert --output-directory \<PATH\> --force --bcl-input-directory \<PATH-TO-RUN\> --sample-sheet SampleSheetV2_GTseq.csv --strict-mode true --no-lane-splitting true

### Using <i>GTseq_Demultiplex.py</i> 

Using the provided example:
```
cd ExampleData
../GTseq_Demultiplex.py --barcodeFile NS2K_85-bc.csv --fastq Undetermined_S0_R1_001.fastq.gz --i5rc Y --seqType GTseq --gzipOutput N
```
This should result in three fastq files i202_03_TOSS2145_Ots_Harvest1.fastq, i202_09_TOSS2145_Ots_Harvest2.fastq, and i211_21_TOSS2154_Ots_Harvest3.fastq

Now GZIP and move them to a new directory:
```
gzip i*
mkdir READS
mv i* ./READS
cd READS
```
### Genotype with <i>GTseq_Genotyper_v5.pl</i> 
With one FASTQ:
```
../../GTseq_Genotyper_v5.pl ../OtsGTseqV9.2_363-probeSeqs.csv i202_03_TOSS2145_Ots_Harvest1.fastq.gz > i202_03_TOSS2145_Ots_Harvest1.genos
```
Parallelized example:

```
ls | grep fastq | while read -r LINE; do
genos_out=$(echo $LINE | sed -e s/fastq.gz/genos/g);
echo "../../GTseq_Genotyper_v5.pl ../OtsGTseqV9.2_363-probeSeqs.csv $LINE > $genos_out" >> GenotypeCommands.cmds;
done  
parallel -j 3 < GenotypeCommands.cmds > Genotyper.log 2>&1 

```

Move '.genos' files to thier own directory

```
mkdir ../GENOS
mv *genos ../GENOS
cd ../GENOS
```

### Compile Genotypes into wide-format CSV file using <i>GTseq_GenoCompile_v5.1.pl</i>
```
../../GTseq_GenoCompile_v5.1.pl . --name=index --output=genotype --filter=0 > NS2K_85-GenosNF.csv
../../GTseq_GenoCompile_v5.1.pl . --name=sample --output=genotype --filter=90 > NS2K_85-ProgGenos90.csv

```
### Optional Analyses
Look for off-target species among these putative chinook:\
CallSpecies.py is available from https://github.com/zakrobinson/CallSpecies
```
./CallSpecies.py --SpeciesSeq SpeciesSeq.csv --inGENO NS2K_85-ProgGenos90.csv --outGENO NS2K_85-ProgGenos90_CALLED.csv
```





