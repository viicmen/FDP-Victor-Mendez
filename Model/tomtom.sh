#!/bin/bash

module load MEME

'''
all.memes is a file with the PWMs of the motifs in the non-redundant dataset (for each database) in MEME format.
'''

for db in 'Cisbp' 'Jaspar' 'Hocomoco';
do
cat $db/pwms/* > $db/all_memes.meme
tomtom -thresh 1 -oc $db/tomtom $db/all_memes.meme $db/all_memes.meme
done
