# myBioTools
Tools I use for tasks I do repetively.


### QuickBlaster

This is used for quick blast screening of a set of sequences in a fast file. It is a simple blastn wrapper that outputs the best alignment for each sequence in the query file.
This may not be useful if one wants to look at multiple hits and compare them.

```
usage: quickBlaster.py [-h] -q QUERY -ht HIT_OUT [-nh NO_HIT] [-db BLASTDB] [-s] [-g {archae,bacteria,eukaryota,virus,all}]
                       [-ng {archae,bacteria,eukaryota,virus,none}]

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Input fasta
  -ht HIT_OUT, --hit_out HIT_OUT
                        Output file for hits
  -nh NO_HIT, --no_hit NO_HIT
                        Output file for no hits
  -db BLASTDB, --blastdb BLASTDB
                        Path to the blast database
  -s, --short           short query sequences
  -g {archae,bacteria,eukaryota,virus,all}, --gi {archae,bacteria,eukaryota,virus,all}
                        restrict search to gilist
  -ng {archae,bacteria,eukaryota,virus,none}, --neg_gi {archae,bacteria,eukaryota,virus,none}
                        exclude gilist from search
```

### QquickBlaster





