#!/bin/bash
base_dir=/Volumes/lab/members/falko/atac-seq/data/atac-seq
output_dir=/Volumes/lab/members/falko/atac-seq/tdf
genome_path=/Users/falko.hofmann/igv/genomes/tair_Nodine.genome

cd $base_dir
file_lis_bt=bam_files_paths_bt.txt
file_lis_bt_2=bam_files_paths_bt2.txt

echo "Running bowtie files..."
while read f; do
  echo "Converting: "$f
  igvtools count -z 10 -w 10 "$base_dir/$f" "$output_dir/$(basename $f).tdf" "$genome_path"
done < $base_dir/$file_lis_bt
echo "Running bowtie files... - Done"

echo "Running bowtie2 files..."
while read f; do
  echo "Converting: "$f
  igvtools count -z 10 -w 10 "$base_dir/$f" "$output_dir/$(basename $f).tdf" "$genome_path"
done < $base_dir/$file_lis_bt_2
echo "Running bowtie2 files... - Done"
