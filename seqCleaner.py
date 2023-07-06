#!/usr/bin/env python3

import sys
import os
import pgzip
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('-fmt', '--format', help='File format', choices=["fastq", "fq", "fasta","fa"], default="fasta")
parser.add_argument('-in', '--input_file', help='Input file', required='True')
parser.add_argument('-out', '--output_file', help='Output file')
parser.add_argument('-min', '--min_seqlen', help='Desired minmum sequence length', type=int, default=200)
parser.add_argument('-f', '--force', help='Force overwrite target file if it exists', action='store_true')
parser.add_argument('-t', '--threads', help='Number of threads. Useful if processing compressed files', type=int)


args = parser.parse_args()

infile = args.input_file
outfile = args.output_file
minlen = args.min_seqlen
forceOverwrite = args.force
ft = args.format

if args.threads:
    proc = args.threads
else:
    proc = 1

if ft in ['fastq', 'fq']:
    fmt = 'fastq'
else:
    fmt = 'fasta'

temp_out = f'temporary.{fmt}'


# Testing

def getREcNumLens(file, fmt, proc): # Get the number of records, lenths of the smallest and largest records.
    L = []
    n = 0

    if file.endswith('gz'):
        with pgzip.open(infile, 'rt', thread=proc, blocksize=4*10**8) as handle:
            for rec in SeqIO.parse(handle, fmt):
                n += 1
                l = len(rec.seq)
                L.append(l)
                L.sort()
                m = L[0]
                p = L[-1]
    else:

        for rec in SeqIO.parse(file, fmt):
            n += 1
            l = len(rec.seq)
            L.append(l)
            L.sort()
            m = L[0]
            p = L[-1]

    return([f'{n:,}', f'{p:,}', f'{m:,}'])


if not outfile:
    outfile = os.path.join(os.path.dirname(infile), str(minlen) + '_' + os.path.basename(infile))

if os.path.exists(outfile):
    if not forceOverwrite:
        userIn = ''
        print('\nThe target file exists. Do you want to overwrite it?\n --> y or yes to confirm; n or no to abort.\n')
        while userIn != 'okay':
            userConf = input()
            if userConf.strip().lower() in 'no':
                userIn = 'okay'
                sys.exit()
            elif userConf.strip().lower() in 'yes':
                userIn = 'okay'
            else:
                print('\nResponse not recognized.\n --> y or yes to confirm; n or no to abort.\n')

records = []


if infile.endswith('gz'):
    with pgzip.open(infile, 'rt') as handle:
        for record in SeqIO.parse(handle, fmt):
            if len(record.seq) >= minlen:
                records.append(record)
else:

    for record in SeqIO.parse(infile, fmt):
        if len(record.seq) >= minlen:
            records.append(record)


SeqIO.write(records, temp_out, fmt)

dataIn = getREcNumLens(infile, fmt, proc)
dataOut = getREcNumLens(temp_out, fmt, proc)


if outfile.endswith('gz'):
    with open(temp_out, 'rt') as fp:
        data = fp.read()
        with pgzip.open(outfile, 'wt', thread=proc, blocksize=4*10**8) as f:
            f.write(data)
    os.remove(temp_out)

else:
    os.rename(temp_out, outfile)


print('')
print(f'Original file ({infile}): {dataIn[0]} records --Largest: {dataIn[1]} --Smallest: {dataIn[2]}')
print(f'Output file ({outfile}): {dataOut[0]} records --Largest: {dataOut[1]} --Smallest: {dataOut[2]}')
print('')
