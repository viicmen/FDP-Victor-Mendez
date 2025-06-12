#!bin/bash/

module load HMMER

'''
Pfam-A.hmm is the whole Pfam hmm database free to download in the InterPro website.
tf_dbd_list.txt is a list of different TF DNA binding domains.
'''

hmmfetch -o tfdbd.hmm -f Pfam-A.hmm tf_dbd_list.txt
hmmpress tfdbd.hmm
hmmscan --domtblout output -o dummy -E 0.001 tfdbd.hmm sequences.fasta
