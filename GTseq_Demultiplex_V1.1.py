#!/usr/bin/python3

# Modification (Kim Vertacnik; May 2022) of GTseq_BarcodeSplit_MP.py (Nate Campbell)
# v1 update: use of defaultdict; argparse arguments; combined GTseq, poolSeq, i5 reverse complement (previously separate scripts)
#v1.1 added sequencerID to argparse arguments.
# Usage: GTseq_Demultiplex.py --barcodeFile= --fastq= --i5rc=Y/N --seqType=GTseq/poolSeq --gzipOutput=Y/(N)
# barcodeFile is a .csv file with each sample's i5 and i7 barcode sequences
#  SampleID,PlateID,i7name,i7sequence,i5name,i5sequence
#   A header line is required
#   The input i5 sequences do not need to be in reverse complement
#    --i5rc for reverse complementing the i5 barcodes (depends on Illumina instrument)
# Input undemultiplexed fastq file must be gzip compressed
# Output fastq files are uncompressed by default
#  --gzipOutput argument is optional (default is N)
#   Gzip in python is really slow: It's faster to output .fastq and use linux gzip to compress
#  --seqType determines the output file name format
#   GTseq: i7name_i5name_PlateID_SampleID.fastq
#   poolSeq: sampleID_Library.fastq
# Update sequencerID (ie, if not the CRITFC NextSeq 2000)
# maxSample and batchSize determine how many processors may be used
# Important! Reads are appended to the output files
#  If this script is run more than once, the original output files must be deleted to avoid duplicate reads

import argparse
import gzip
import os
import shutil
import sys
from multiprocessing import Process
from time import sleep
from datetime import datetime
from collections import defaultdict

maxSample = 4000    # Currently capped at 40 processing cores
batchSize = 100

parser = argparse.ArgumentParser()
parser.add_argument('--barcodeFile',
                    type=str,
                    required=True,
                    metavar='<sample i5 and i7 barcode sequences>'
                    )
parser.add_argument('--fastq',
                    type=str,
                    required=True,
                    metavar='<undemultiplexed fastq.gz file>'
                    )
parser.add_argument('--i5rc', 
                    type=str, 
                    required=True,
                    choices=['Y', 'N'],
                    help='Reverse complement i5 barcode sequences (Y/N)'
                    )
parser.add_argument('--seqType', 
                    type=str, 
                    required=True,
                    choices=['GTseq', 'poolSeq'],
                    )
parser.add_argument('--sequencerID', 
                    type=str, 
                    required=False,
                    default='@VH00647',
                    metavar='',
                    help="sequencerID in header lines default: %(default)s"
                    )
parser.add_argument('--gzipOutput',
                    nargs='?',
                    type=str,
                    choices=['Y', 'N'],
                    const=1,
                    default='N',
                    help='Compress output fastq files (optional, default is N)'
                    )
args = parser.parse_args()

sequencerID = args.sequencerID

list_dict = defaultdict(list)

individuals = 0
start = 1
end = 0
csv_lineNo = 0

f = open(args.barcodeFile, 'r')
for line in f:
    csv_lineNo = csv_lineNo + 1
file_end = csv_lineNo
individuals = csv_lineNo - 1
if individuals > batchSize + start:
    end = batchSize + start
else:
     end = file_end
f.close()

sets = 0
Samples = True
starting_ind = 1

while Samples == True:    # While loop is broken with "end == file_end' or "start > file_end" conditions
    csv_lineNo2 = 0
    
    f = open(args.barcodeFile, 'r')
    list0 = []
    for line in f:
        csv_lineNo2 = csv_lineNo2 + 1
        if csv_lineNo2 > start and csv_lineNo2 <= end:
            list0.append(line)
    sets = sets + 1
    size = len(list0)
    list_dict[sets] = list0
    
    ending_ind = end - 1

    if end == file_end:
        Samples = False
        break
    elif start > file_end:
        Samples = False
        break

    start = end
    end = start + batchSize
    starting_ind = starting_ind + batchSize

    if end > file_end:
        end = file_end
        ending_ind = file_end - 1

    f.close()

def Main():
    startTime = datetime.now()
    print("Start: ", datetime.today().strftime('%Y-%m-%d %H:%M:%S'))

    bcLinesN = len(open(args.barcodeFile).readlines())    # Exit if barcodeFile exceedes maxSample value
    if bcLinesN > maxSample:
            sys.exit("Sample sheet exceeds limit of " + str(maxSample) + " samples")

    if args.i5rc == 'Y':
        print("Sample sheet i5 barcodes will be reverse complemented")
    else:
        print("Sample sheet i5 barcodes will NOT be reverse complemented")
    
    if args.gzipOutput == 'Y':
        try:
            os.mkdir('output')
        except:
            pass

    process_dict = {}

    for batch in list_dict:
        process_dict[batch] = Process(target=split_file, args=(list_dict[batch],))
        #print(f"Starting batch {batch}")    # redrum python3.6.5 doesn't like f-strings?
        print("Starting batch {}".format(batch))
        process_dict[batch].start()
    
    print("Keepin it 100! Your files will be ready shortly...")
    while process_dict[batch].is_alive():
        sleep(5)

    if args.gzipOutput == 'Y':
        compress_output()

    print("End: ", datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
    endTime = datetime.now()
    elapsedTime = endTime-startTime
    print("Run time: %02d:%02d:%02d:%02d" % (elapsedTime.days, elapsedTime.seconds // 3600, elapsedTime.seconds // 60 % 60, elapsedTime.seconds % 60))
    return

def rev_comp(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return "".join(complement[base] for base in reversed(seq))

def split_file(individual_list):
    fqInGz = gzip.open(args.fastq, mode="rt")
    BC_Dict = {}
    File_Dict = {}
    handle_dict = {}

    for lines in individual_list:
        stuff = lines.rsplit("\n")[0].split(',')    # sampleName,PlateID,i7name,i7sequence,i5name,i5sequence
        if args.i5rc == 'Y':
            stuff[5] = rev_comp(stuff[5])
        if args.seqType == 'GTseq':
            name = stuff[2] + '_' + stuff[4] + '_' + stuff[1] + '_' + stuff[0] + '.fastq'    # i7name_i5name_PlateID_sampleName.fastq
        elif args.seqType == 'poolSeq':
            name = stuff[0] + '_' + stuff[1] + '.fastq'    # sampleName_LibraryID.fastq
        BC_Dict[name] = stuff[3] + '+' + stuff[5]    # i7sequence+i5sequence
        File_Dict[stuff[3] + '+' + stuff[5]] = name
        if args.gzipOutput == 'N':
            handle_dict[stuff[3] + '+' + stuff[5]] = open(File_Dict[stuff[3] + '+' + stuff[5]], 'a')    # reads are appended
        elif args.gzipOutput == 'Y':
            handle_dict[stuff[3] + '+' + stuff[5]] = open(os.path.join('output',File_Dict[stuff[3] + '+' + stuff[5]]), 'a')    # reads are appended

    writelines = 0
    for line in fqInGz:
        if writelines < 5 and writelines > 0:
            f_out.write(line)
            writelines = writelines + 1
            if writelines == 4:
                writelines = 0
        if sequencerID in line:
            info = line.rsplit("\n")[0].split(':')    # @VH00647:7:AAAN5NVM5:1:1101:64813:1114 1:N:0:GTGGTT+GTGCGC
            BC = info[9]    # read i7sequence+i5sequence
            if BC in File_Dict:
                f_out = handle_dict[BC]
                f_out.write(line)
                writelines = 1
    fqInGz.close()
    return

def compress_output():
    for filename in os.listdir('output'):
        if filename.endswith('.fastq'):
            infile = open(os.path.join('output',filename),'rb')
            data = bytearray(infile.read())
            infile.close()
            outfile = gzip.open(filename + '.gz','wb')
            outfile.write(data)
            outfile.close()
    try:
        shutil.rmtree('output')
    except:
        pass
    return

if __name__ == '__main__':
    Main()
