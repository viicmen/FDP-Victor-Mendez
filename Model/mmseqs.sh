#!/bin/bash

module load MMseqs2

'''
modcre_pdbs.fasta is a multifasta with all the sequences that ModCRE utilise for modelling structures
filtered_sequences.fasta is the multifasta (for each database) with the sequences of the motifs in the non-redundant dataset.
all_sequences.fasta is similar to "filtered_sequences.fasta" but for the redundant dataset.
'''

mmseqs createdb modcre_pdbs.fasta pdb_fasta_db
mmseqs createindex pdb_fasta_db tmp

for db in 'JASPAR' 'cis-bp' 'hocomoco';
do

mmseqs createdb $db/filtered_sequences.fasta $db/filtered_fasta_db
mmseqs createindex $db/filtered_fasta_db tmp

mmseqs search -s 7.5 --max-seq-id 1.0 --num-iterations 4 --max-seqs 10000 $db/filtered_fasta_db pdb_fasta_db $db/search_modcre tmp
mmseqs convertalis $db/filtered_fasta_db pdb_fasta_db $db/search_modcre $db/search_modcre.m8

mmseqs createdb $db/all_sequences.fasta $db/all_fasta_db
mmseqs createindex $db/all_fasta_db tmp

mmseqs search -s 7.5 --max-seq-id 1.0 --num-iterations 4 --max-seqs 10000 $db/filtered_fasta_db $db/all_fasta_db $db/search_nn tmp
mmseqs convertalis $db/filtered_fasta_db $db/all_fasta_db $db/search_nn $db/search_nn.m8
done
