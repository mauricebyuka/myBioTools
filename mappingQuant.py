#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 22:16:27 2022

@author: mpb5554
"""

import os
import subprocess
import sys
import csv
import argparse


parser = argparse.ArgumentParser()

required = parser.add_argument_group('Required arguments')

required.add_argument("-f", "--bamDir",
                    help="Folder containing bam files.", required=True)
required.add_argument("-o", "--outFile",
                    help=f"Output file", required=True)
args = parser.parse_args()


def quantifyMap(bam):
     
    cmd1 = f'samtools view -c {bam}'
    execute1 = subprocess.run(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if execute1.returncode == 0:
        all_reads = int(execute1.stdout.decode('utf-8').strip())
    else:
        print(f'Something is wrong with the bam file: {bam}')
        sys.exit()

    cmd2 = f'samtools view -c -F 260 {bam}'
    execute2 = subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if execute2.returncode == 0:
        map_reads = int(execute2.stdout.decode('utf-8').strip())#.split('\n')
    else:
        print(f'Something is wrong with the bam file: {bam}')
        sys.exit()    
    
    percentMap = round((map_reads/all_reads)*100, 2)


    return [all_reads, map_reads, percentMap]

dataDir = args.bamDir
mappReport = args.outFile

samples = []
for file in os.listdir(dataDir):
    if file.endswith('_raw.bam'):
        sample = file.strip('_raw.bam')
        samples.append(sample)

with open(mappReport, 'w') as report:
    print('\n Generating mapping summary report\n')
    writer = csv.writer(report, dialect = csv.excel_tab)
    writer.writerow(['sample', 'rtot', 'rmapped', 'perc_rmapped', 'dtot', 'dmapped', 'perc_dmapped', 'perc_dtot', 'perc_dmapped'])
    for sample in sorted(samples):
        data = [sample]
        raw = os.path.join(dataDir, f'{sample}_raw.bam')
        dedup = os.path.join(dataDir, f'{sample}_dedupped.bam')
        L = quantifyMap(raw)
        M = quantifyMap(dedup)
        N = [round(M[0]/L[0]*100, 2), round(M[1]/L[1]*100, 2)]
        P = L + M + N
        data.extend(P)
        writer.writerow(data)
        print(f'--- Added data for sample: {sample} ---')

print('\nDone!\n')

