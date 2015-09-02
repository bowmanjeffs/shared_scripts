#!/bin/bash
## run as run_pplacer_single.sh [query].fasta [reference package]
## if reference alignment is large this script could take a long time!

## get the basename from query input name
b=$(basename "$1" .fasta) && 
echo $b &&
echo $2 &&

## remove old files, if they exist
rm $b.clean.joined.aligned.jplace
rm $b.clean.joined.aligned.csv
rm $b.clean.joined.aligned.tre

## get rid of all the punctuation that breaks pplacer
## there could be more bad punctuation not listed here!
tr " " "_" < $b.fasta | tr -d "%" | tr -d "\'" | tr -d "," | tr -d ";" | tr -d "(" | tr -d ")" | tr ":" "_" | tr "=" "_" | tr "." "_" | tr -d "\"" | tr -d "\*" | tr -d "[" | tr -d "]" | tr "-" "_" > $b.clean.fasta &&

## join the reference alignment (found in the reference package) with the cleaned up query
cat $b.clean.fasta $2/*.fasta > $b.clean.joined.fasta &&

## use clustalo to align
clustalo_1.2 --dealign --force --full --iterations=3 -i $b.clean.joined.fasta -o $b.clean.joined.aligned.fasta &&

## run pplacer
pplacer -p --keep-at-most 10 -c $2 $b.clean.joined.aligned.fasta &&

## run guppy to create a fat tree and a csv file of placements
guppy fat --point-mass --pp --edpl 0.03 $b.clean.joined.aligned.jplace -o $b.clean.joined.aligned.phyloxml &&
guppy to_csv --point-mass --pp $b.clean.joined.aligned.jplace > $b.clean.joined.aligned.csv
