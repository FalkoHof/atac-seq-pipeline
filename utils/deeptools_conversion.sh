#!/bin/bash
############
# Convert bam to bw and run deeptools multi bam summary
# Parameters:
# 1 input dir containing the bam files
# 2 bed file
# 3 output dir for big wig files
# 4 output dir for deeptools results
# 5 temp_dir
############
# make sure the following modules in your path
# ml deepTools/2.2.4-foss-2015a-Python-2.7.9
# ml OpenSSL/1.0.1p-foss-2015a

bam_dir=$1
bed_file=$2
wig_dir=$3
dt_out=$4
TMPDIR=$5
#arabidopsis nuclear genome size
tair10_size=119146348

if [ ! -d "$bam_dir" ]; then
  echo "Input directory $bam_dir does not exist!" 1>&2
  exit 1
fi

dir_arr=($wig_dir $dt_out $TMPDIR)

for dir in ${dir_arr}; do
  if [ ! -d "$dir" ]; then
    mkdir -pv $dir
  fi
done

#the the filenames
bam_files=$(ls $bam_dir | grep .bam)

#do the conversion
for f in ${bam_files}; do
  bamCoverage \
    -b $bam_dir/$f  \
    -o $wig_files/${f%.*}.bw \
    --binSize=1 \
    --normalizeTo1x $tair10_size \
    --ignoreDuplicates \
    --maxFragmentLength=100
done

#function to run for making the different bam summeries
dt_multi_bam_summmary ()
{
  multiBamSummary \
    $1
    -b $2 \
    -out $3
}

#parameters on what to run
param_1=()
param_1+=("BED-file --BED $bed_file")
param_1+=("bins -bs 100")
param_1+=("bins -bs 1000")

param_2=$bam_files

param_3=()
param_3+=("$dt_out/bw_compare_${bed_file%%.*}.npz")
param_3+=("$dt_out/bw_compare_bins_100bp.npz")
param_3+=("$dt_out/bw_compare_bins_1000bp.npz")

#some nested loops that execture the respective commands for the comparison
for p1 in $param_1
do
  for p3 in param_3
  do
    dt_multi_bam_summmary $p1 $p2 $p3
  done
done
