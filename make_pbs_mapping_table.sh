#!/bin/bash
input_dir=$1 #uses folder supplied by command line args
x=1
find $input_dir -type f -name "*.bam" | while read line
do
  b=$(basename $line)
  y=${b%.bam}
  printf "%d %s\n" $x $y
  ((x++))
done
