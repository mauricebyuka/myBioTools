#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
import sys
import argparse
from Bio import SearchIO

parser = argparse.ArgumentParser()

parser.add_argument('-q', '--query', help='Input fasta', required='True')
parser.add_argument('-ht', '--hit_out', help='Output file for hits', required=True)
parser.add_argument('-nh', '--no_hit', help='Output file for no hits', default='no_hit.tsv')
parser.add_argument('-db', '--blastdb', help='Path to the blast database', default='/media/mpb5554/data2/blastdb/nt')
parser.add_argument('-s', '--short', help='short query sequences', action='store_true')
parser.add_argument('-g', '--gi', help='restrict search to gilist ', choices=['archae', 'bacteria', 'eukaryota', 'virus', 'all'], default='all')
parser.add_argument('-ng', '--neg_gi', help='exclude gilist from search', choices=['archae', 'bacteria', 'eukaryota', 'virus', 'none'], default='none')

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
cpus = int(int(os.cpu_count()) * 3 / 4)

######
gil = args.gi

if gil == 'all':
    cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
    cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
else:
    gilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{gil}')
    cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
    cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -gilist {gilis} -query {infile} -out {blast_out} -outfmt 5'
######

######
ngil = args.neg_gi

if ngil == 'none':
    cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
    cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -query {infile} -out {blast_out} -outfmt 5'
else:
    ngilis = os.path.join(os.path.dirname(args.blastdb), f'gi_lists/{ngil}')
    cmd = f'blastn -num_threads {cpus} -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
    cmd_s = f'blastn -num_threads {cpus} -task blastn-short -db {args.blastdb} -negative_gilist {ngilis} -query {infile} -out {blast_out} -outfmt 5'
######

if not os.path.exists(blast_out):
    print("\n Runing BLAST. May take some time ......\n")
    if args.short:
        os.system(cmd_s)
    else:
        os.system(cmd)   
else:
    print("\n Blast results file exits. Extracting data from existing file.\n")

if os.path.exists(blast_out):
    with open(blast_results, "w") as results, open(no_hits, "w") as alt:
        res_writer = csv.writer(results, dialect=csv.excel_tab)
        alt_writer = csv.writer(alt, dialect=csv.excel_tab)
        res_writer.writerow(['Query_id', 'Best hit', 'Accession #', 'Percent ID', 'Hit length', 'Hit start', 'Hit end', 'Query length', 'Query span',  'Query cov', 'Query start', 'Query end'])
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
                print(f' {query_id}::{accession}::{name}::{percent_id}::{qcov}')
                print('\tHit start: ', best_hsp.hit_start_all)
                print('\tHit end: ', best_hsp.hit_end_all)
                print('\t--------------------------------')
                print('\tQuery length: ', qlen)
                print('\tQuery start: ', best_hsp.query_start_all)
                print('\tQuery end: ', best_hsp.query_end_all)
                print('\n')

                res_writer.writerow([query_id, name, accession, percent_id, best_hit.seq_len, hstart, hend, qlen, best_hsp.query_span, qcov, qstart, qend])
            else:
                alt_writer.writerow([record.id, "None found", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"])
    print(" Done!!\n")
else:
    print(" Blast results file was not found!\n")

# print(f'{acession}:: {name}:: {percent_id}')
# print('\tHit start: ', best_hsp.hit_start)
# print('\tHit end: ', best_hsp.hit_end)
# print('\t--------------------------------')
# print('\tQuery start: ', best_hsp.hit_start)
# print('\tQuery end: ', best_hsp.hit_all)
# print('')
