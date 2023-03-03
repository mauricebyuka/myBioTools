#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import sys
import argparse
from Bio import SearchIO

parser = argparse.ArgumentParser(add_help=False)

required = parser.add_argument_group('Required arguments')
required.add_argument('-q', '--query', help='Input fasta', required=True)
required.add_argument('-ht', '--hit_out', help='Output file for hits', required=True)
required.add_argument('-db', '--blastdb', help='Path to the blast database', required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
optional.add_argument('-nh', '--no_hit', help='Output file for no hits', default='no_hit.tsv')
optional.add_argument('-p', '--blastProgram', help='Blast program to be used', choices=['blastn', 'blastp'], default='blastn')
optional.add_argument('-s', '--short', help='Short query sequences', action='store_true')
optional.add_argument('-c', '--clean', help='Delete intermediate xml files', action='store_true')
optional.add_argument('-t', '--threads', help='Threads to be used by blastn [Default: 3/4 of available cpus]')

filtering = parser.add_argument_group('Arguments for filtering')

filtering.add_argument('-g', '--gi', help='restrict search to specified domain gilist ', choices=['archae', 'bacteria', 'eukaryota', 'virus'], default='all')
filtering.add_argument('-ng', '--neg_gi', help='exclude specified domain gilist from search', choices=['archae', 'bacteria', 'eukaryota', 'virus'], default='none')

args = parser.parse_args()

infile = args.query
blast_results = args.hit_out
base_dir = os.path.dirname(os.path.abspath(blast_results))
base_name = os.path.basename(blast_results).split('.')[0]
if not os.path.exists(base_dir):
    os.makedirs(base_dir)

nh_path = args.no_hit

if not nh_path:
    no_hits = os.path.join(base_dir, 'no_hit.tsv')
else:
    no_hits = nh_path

blast_out = os.path.join(base_dir, f'{base_name}_blast_results.xml')

Threads = args.threads
if Threads:
    cpus = Threads
else:
    cpus = int(int(os.cpu_count()) * 3 / 4)

######
gil = args.gi
ngil = args.neg_gi
prog = args.blastProgram


def cmdgen(program, gil,ngil):
    if gil == 'all' and ngil == 'none':
        cmd = f'{program} -num_threads {cpus} -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
        cmd_s = f'{program} -num_threads {cpus} -task blastn-short -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
        print('\nNo gi list provided: The entire database will be searched.\n')
    elif gil != 'all' and ngil == 'none':
        gilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{gil}')
        cmd = f'{program} -num_threads {cpus} -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
        cmd_s = f'{program} -num_threads {cpus} -task blastn-short -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
        print(f'\nA gi list was provided. The search will be limited to the {gil} entries.\n')

    elif ngil != 'none' and gil == 'all':
        ngilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{ngil}')
        cmd = f'{program} -num_threads {cpus} -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
        cmd_s = f'{program} -num_threads {cpus} -task blastn-short -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
        print(f'A negative gi list was provided. {ngil} will be excluded from search.')
    else:
        print('\nOnly one of -g (--gi) or -ng (--neg_gi) can be used. Not both at the same time.\n')
        sys.exit()

    return [cmd, cmd_s]

commands = cmdgen(prog, gil,ngil)
cmd = commands[0]
cmd_s = commands[1]

# if gil == 'all' and ngil == 'none':
#     cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
#     cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
#     print('\nNo gi list provided: The entire database will be searched.\n')
# elif gil != 'all' and ngil == 'none':
#     gilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{gil}')
#     cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
#     cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
#     print(f'\nA gi list was provided. The search will be limited to the {gil} entries.\n')

# elif ngil != 'none' and gil == 'all':
#     ngilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{ngil}')
#     cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
#     cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
#     print(f'A negative gi list was provided. {ngil} will be excluded from search.')
# else:
#     print('\nOnly one of -g (--gi) or -ng (--neg_gi) can be used. Not both at the same time.\n')
#     sys.exit()
######

if not os.path.exists(blast_out):
    print("\n Runing BLAST. May take some time ......\n")
    if args.short:
        print(f'\n\tRunning--->: {cmd_s}\n')
        os.system(cmd_s)
    else:
        print(f'\n\tRunning--->: {cmd}\n')
        os.system(cmd)
else:
    print("\n Blast results file exists. Extracting data from existing file.\n")

if os.path.exists(blast_out):
    with open(blast_results, "w") as results, open(no_hits, "w") as alt:
        res_writer = csv.writer(results, dialect=csv.excel_tab)
        alt_writer = csv.writer(alt, dialect=csv.excel_tab)
        res_writer.writerow(['Query_id', 'Best hit', 'Accession #', 'Percent ID', 'Hit length', 'Hit start', 'Hit end', 'Query length', 'Query span', 'Query cov', 'Query start', 'Query end', 'Query frame', 'Hit frame'])
        alt_writer.writerow(['Query_id', 'Best hit', 'Accession #', 'Percent ID', 'Hit length', 'Query length', 'HSP length', 'Hit start', 'Hit end', 'Query cov', 'q_start', 'q_end'])
        for record in SearchIO.parse(blast_out, "blast-xml"):

            if not len(record.hits) == 0:
                best_hit = record[0]
                best_hsp = best_hit[0]
                query_id = best_hsp.query_id
                accession = best_hit.accession
                os.system(f'blastdbcmd -db {args.blastdb} -entry {accession} -outfmt "%t" > log.txt')
                with open('log.txt') as log:
                    name = log.readline().strip()
                hsp_len = best_hsp.aln_span
                qstart = best_hsp.query_start + 1
                qend = best_hsp.query_end
                percent_id = str(round(((best_hsp.ident_num / best_hsp.aln_span) * 100), 2)) + "%"
                HSPperQuery = str(best_hsp.aln_span) + "/" + str(record.seq_len)
                hstart = best_hsp.hit_start + 1
                hend = best_hsp.hit_end
                qlen = record.seq_len
                qcov = str(round(((best_hsp.query_span / qlen) * 100), 2)) + "%"
                hframe = best_hsp.hit_frame
                qframe = best_hsp.query_frame
                print(f' {query_id}::{accession}::{name}::{percent_id}::{qcov}')
                print('\tHit start: ', best_hsp.hit_start_all)
                print('\tHit end: ', best_hsp.hit_end_all)
                print('\t--------------------------------')
                print('\tQuery length: ', qlen)
                print('\tQuery start: ', best_hsp.query_start_all)
                print('\tQuery end: ', best_hsp.query_end_all)
                print('\n')

                res_writer.writerow([query_id, name, accession, percent_id, best_hit.seq_len, hstart, hend, qlen, best_hsp.query_span, qcov, qstart, qend, qframe, hframe])
            else:
                alt_writer.writerow([record.id, "None found", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"])

    cleaning = args.clean
    if cleaning:
        os.system(f'rm {blast_out}')
    print(" Done!!\n")
else:
    print(" Blast results file was not found!\n")
