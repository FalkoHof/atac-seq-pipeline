#!/bin/bash
############
# Merge bedfiles with bedtools merge and check for invalid regions (0 bp)
# Parameters:
# 1 input globbing pattern
# 2 output file
############
#command line arguments
f_in=$1
f_out=$2
TMPDIR=$3

#zero index based
start_col=1
stop_col=2

#get random value for unique tmp_id
tmp_id=$RANDOM

#merge all peaks
cat $f_in | sort -k1,1V -k2,2n -T $TMPDIR | bedtools merge > \
  $TMPDIR/$tmp_id.bed

#run a small python script to kick out peaks with the size of 0 if present
python check_bed.py $TMPDIR/$tmp_id.bed 1 2 \
  2> $f_out.log > $f_out

#remove the tmp file
rm $TMPDIR/$tmp_id.bed
