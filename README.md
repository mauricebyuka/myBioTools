# myBioTools
This repository contains the scripts that automate some of the tasks that I use on a regular basis. They were not designed with a general case usage in mind but rather for perfoming specific tasks, although some of them may allow some level of flexibility. Below is a list of the scripts and their usage. 

### quickBlaster

This is used for quick blast screening of a set of sequences in a fast file. It is a simple blastn wrapper that outputs the best alignment for each sequence in the query file. This may not be useful if one wants to look at multiple hits and compare them.

**Dependancies:** [Command line Blast](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [Biopython](https://biopython.org/wiki/Download).

```
usage: quickBlaster.py -q QUERY -ht HIT_OUT [-db BLASTDB] [-h] [-nh NO_HIT] [-s] [-g {archae,bacteria,eukaryota,virus,all}]
                       [-ng {archae,bacteria,eukaryota,virus,none}]

Required arguments:
  -q QUERY, --query QUERY
                        Input fasta
  -ht HIT_OUT, --hit_out HIT_OUT
                        Output file for hits
  -db BLASTDB, --blastdb BLASTDB
                        Path to the blast database

Optional arguments:
  -h, --help            show this help message and exit
  -nh NO_HIT, --no_hit NO_HIT
                        Output file for no hits
  -s, --short           short query sequences

Arguments for filtering:
  -g {archae,bacteria,eukaryota,virus,all}, --gi {archae,bacteria,eukaryota,virus,all}
                        restrict search to gilist
  -ng {archae,bacteria,eukaryota,virus,none}, --neg_gi {archae,bacteria,eukaryota,virus,none}
                        exclude gilist from searc
```

### seqCleaner

The script will take a fasta file and remove sequences that are shorter than a specified minimum length. It will report the number of records, the largest and shorted records for the input and output files. There are multiple tools to do the same task.

**Dependancies:**  [Biopython](https://biopython.org/wiki/Download).

```
usage: seqCleaner.py [-h] -in INPUT_FASTA [-out OUTPUT_FASTA] [-min MIN_SEQLEN] [-f]

optional arguments:
  -h, --help            show this help message and exit
  -in INPUT_FASTA, --input_fasta INPUT_FASTA
                        Input fasta
  -out OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        Output file
  -min MIN_SEQLEN, --min_seqlen MIN_SEQLEN
                        Desired minmum sequence length
  -f, --force           Force overwrite target file if it exists

```

### mappingQuant

Teh script works on a collection of pairs of alignment (bam) files : raw files and dedupped files (after removal of aptical and PCR duplicates). The script outputs a table with mapping statistics. At this point, I use the script to works with the output files from the [GenMapViz](https://github.com/mauricebyuka/GenMapViz) program. Raw bam files have to end in _raw.bam and dedupped files in _dedupped.bam.

Working on making it more flexible.

**Dependancies:**  [SAMtools](http://www.htslib.org/).

```
optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -f BAMDIR, --bamDir BAMDIR
                        Folder containing bam files.
  -o OUTFILE, --outFile OUTFILE
                        Output file

```


