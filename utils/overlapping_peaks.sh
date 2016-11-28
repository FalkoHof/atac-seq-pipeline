#!/bin/bash
############
# Perform bedtools substract and intersect on a set of bed files.
# Parameters:
# 1 globbing pattern for the path conatining the bed files
# 2 output dir for bedtools
############
# make sure beedtools is in your path
# ml BEDTools/2.26.0-foss-2016a

bed_files=$1
bedtools_out=$2

mkdir -p $bedtools_out
mkdir $deeptools_out

#1. substract and intersect bed files
for f1 in $bed_files
do
  for f2 in $bed_files
  do
    if [ "$f1" != "$f2" ]; then
      #TODO:implement sample name getter
      name1=${f1%%.*}
      name2=${f2%%.*}
      bedtools substract -a $f1 - b $f2 > $bedtools_out/$name1$name2_substract.bed
      bedtools intersect -a $f1 - b $f2 -wa > $bedtools_out/$name1$name2_intersect.bed
    fi
  done
done
