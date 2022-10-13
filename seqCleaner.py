#!/usr/bin/env python3

import sys, os
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('-in', '--input_fasta', help='Input fasta', required='True')
parser.add_argument('-out', '--output_fasta', help='Output file')
parser.add_argument('-min', '--min_seqlen', help='Desired minmum sequence length', type=int, default=200)
parser.add_argument('-f', '--force', help='Force overwrite target file if it exists', action='store_true')


args = parser.parse_args()

infile = args.input_fasta
outfile = args.output_fasta
minlen = args.min_seqlen
forceOverwrite = args.force


def getREcNumLens(file): # Get the number of records, lenths of the smallest and largest records.
    first_record = next(SeqIO.parse(file, 'fasta'))
    m = p = len(first_record.seq)
   # p = len(first_record.seq)
    n = 0
    recs = []
    for rec in SeqIO.parse(file, 'fasta'):
        n += 1
        if len(rec.seq) < m:
            m = len(rec.seq)
        else:
            p = len(rec.seq)
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
for record in SeqIO.parse(infile, 'fasta'):
    if len(record.seq) >= minlen:
        records.append(record)

SeqIO.write(records, outfile, 'fasta')

dataIn = getREcNumLens(infile)
dataOut = getREcNumLens(outfile)

print('')
print(f'Original file ({infile}): {dataIn[0]} records --Largest: {dataIn[1]} --Smallest: {dataIn[2]}')
print(f'Output file ({outfile}): {dataOut[0]} records --Largest: {dataOut[1]} --Smallest: {dataOut[2]}')
print('')

