#!/bin/bash
input_dir=$1 #uses folder supplied by command line args
x=1
find $input_dir -type f -name "*.bam" | while read line
do
  b=$(basename $line)
  printf "%d %s\n" $x $b
  ((x++))
done
