#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_peak-calling
#PBS -J 1-17
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/161124_atac-seq/peak-calling/161124_atac-seq_^array_index^_peak_calling.log
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=48gb

# === begin ENVIRONMENT SETUP ===
#some parameters on what to execute.
#TODO:
split_files=1
create_bed=1
run_fseq=1
run_macs2=1
create_wig=1
clean=1
#set to the number of available cores
threads=8
#path from which the script is exectuted
pipe_dir=/lustre/scratch/users/$USER/pipelines/atac-seq-pipeline
#set ouput base dir
base_dir=/lustre/scratch/users/$USER/atac-seq
#location of the mapping file for the array job
pbs_mapping_file=$pipe_dir/pbs_mapping_file.txt
#super folder of the temp dir, script will create subfolders with $sample_name
temp_dir_base=$base_dir/temp
#json filters for bamtools
subnucl_filter=$pipe_dir/bamtools_filter/bamtools_subnucl.json
nucl_filter=$pipe_dir/bamtools_filter/bamtools_polynucl.json
#download and compile fseq from git & set max heap space in fseq exectuable
fseq=/lustre/scratch/users/$USER/software/F-seq/dist~/fseq/bin/fseq
#igvtools=/lustre/scratch/users/$USER/software/IGVTools/igvtools
tair10_size=119146348

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as line number
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $pbs_mapping_file` #read mapping file
names_mapped=($input_mapper)
sample_dir=${names_mapped[1]} # get the sample dir
sample_name=`basename $sample_dir` #get the base name of the dir as sample name

temp_dir=$temp_dir_base/$sample_name

mkdir -p $temp_dir

#print some output for logging
echo '#########################################################################'
echo 'Starting ATAC-seq peak-calling pipeline for: ' $sample_name
echo 'Sample directory: ' $sample_dir
echo 'Mapping file: ' $pbs_mapping_file
#IDEA: maybe add some diagnostic output the chosen modes.
echo '#########################################################################'

# Load the required modules each loop due to deeptools loading differnt libs
ml BamTools/2.4.0-foss-2016a
ml BEDTools/2.26.0-foss-2016a
ml Java/1.8.0_66


#some error handling function
function error_exit
{
  echo "$1" 1>&2
  exit 1
}

#TODO: add conditional initilization..
#set other temp dir location
TMPDIR=$temp_dir

bam_files=$sample_dir/alignments
split_bam=$bam_files/split_bam
bed_files=$sample_dir/bed_files
fseq_files=$sample_dir/fseq_peaks
macs2_files=$sample_dir/macs2_peaks
wig_files=$sample_dir/wig_files

f=($(ls $bam_files | grep -e "unique\.bam$"))

if [[ "${#f[@]}" -ne "1" ]]; then
  error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
fi

if [ $split_files -eq 1 ]; then
  echo "Splitting bam files..."
  mkdir -p $split_bam
  #get the subnucleosomal reads and sort them
  bamtools filter -in $bam_files/$f -script $subnucl_filter | \
    samtools sort -m 3G -@ $threads - -o $split_bam/${f%.*}.subnucl.bam
  samtools index $split_bam/${f%.*}.subnucl.bam
  #get the nucleosomal reads and sort them
  bamtools filter -in $bam_files/$f -script $nucl_filter | \
    samtools sort -m 3G -@ $threads - -o $split_bam/${f%.*}.nucl.bam
  samtools index $split_bam/${f%.*}.nucl.bam
  echo "Splitting bam files... - Done"
fi

if [ $create_bed -eq 1 ]; then
  echo "Converting bam files to bed..."
  mkdir -p $bed_files
  #convert to bed, keep reads from nuclear chromosomes, then sort and store them
  bedtools bamtobed -i $bam_files/$f | grep "^Ath_chr[1-5]" |\
    sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.bed
  bedtools bamtobed -i $split_bam/${f%.*}.subnucl.bam | \
    grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.subnucl.bed
  bedtools bamtobed -i $split_bam/${f%.*}.nucl.bam | \
    grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.nucl.bed
  echo "Converting bam files to bed... - Done"
fi

if [ $run_fseq -eq 1 ]; then
  echo "Peak-calling with f-seq..."
  mkdir -p $fseq_files/combined
  mkdir -p $fseq_files/subnucl
  mkdir -p $fseq_files/nucl
  #run fseq on the different files...
  sh $fseq -v -f 0 -of npf -t 2.0 -o $fseq_files/combined \
    $bed_files/${f%.*}.bed
  sh $fseq -v -f 0 -of npf -t 2.0 -o $fseq_files/subnucl \
    $bed_files/${f%.*}.subnucl.bed
  sh $fseq -v -f 0 -of npf -t 2.0 -o $fseq_files/nucl \
    $bed_files/${f%.*}.nucl.bed
  echo "Peak-calling with f-seq... - Done"

  echo "Merging f-seq files..."
  #put the fseq output file names in some arrays
  fseq_combined=($(ls $fseq_files/combined | grep -e "Ath_chr[1-5].npf"))
  fseq_subnucl=($(ls $fseq_files/subnucl | grep -e "Ath_chr[1-5].npf"))
  fseq_nucl=($(ls $fseq_files/nucl | grep -e "Ath_chr[1-5].npf"))
  #concatenate them to one file and sort
  cd $fseq_files/combined
  cat ${fseq_combined[@]} | sort -k 1,1 -k2,2n > \
    $fseq_files/combined/${f%.*}_combined_fseq.npf
  rm -v ${fseq_combined[@]}

  cd $fseq_files/subnucl
  cat ${fseq_subnucl[@]} | sort -k 1,1 -k2,2n > \
    $fseq_files/subnucl/${f%.*}_subnucl_fseq.npf
  rm -v ${fseq_subnucl[@]}

  cd $fseq_files/nucl
  cat ${fseq_nucl[@]} | sort -k 1,1 -k2,2n > \
    $fseq_files/nucl/${f%.*}_nucl_fseq.npf
  rm -v ${fseq_nucl[@]}
  echo "Merging f-seq files... - Done"
fi

if [ $create_wig -eq 1 ]; then
  echo "Creating normalized bigwig files..."
  mkdir -p $wig_files
  #load module
  ml deepTools/2.2.4-foss-2015a-Python-2.7.9
  ml OpenSSL/1.0.1p-foss-2015a

  bamCoverage \
    -b $bam_files/$f \
    -o $wig_files/${f%.*}.bw \
    --binSize=1 \
    --normalizeTo1x $tair10_size \
    --ignoreDuplicates \
    --numberOfProcessors=$threads \
    --ignoreForNormalization Ath_chrm Ath_chrc

  bamCoverage \
    -b $split_bam/${f%.*}.subnucl.bam  \
    -o $wig_files/${f%.*}.subnucl.bw \
    --binSize=1 \
    --normalizeTo1x $tair10_size \
    --ignoreDuplicates \
    --numberOfProcessors=$threads \
    --ignoreForNormalization Ath_chrm Ath_chrc

  bamCoverage \
    -b $split_bam/${f%.*}.nucl.bam  \
    -o $wig_files/${f%.*}.nucl.bw \
    --binSize=1 \
    --normalizeTo1x $tair10_size \
    --ignoreDuplicates \
    --numberOfProcessors=$threads \
    --ignoreForNormalization Ath_chrm Ath_chrc
  echo "Creating normalized bigwig files... - Done"
fi

if [ $run_macs2 -eq 1 ]; then
  echo "Peak-calling with MACS2"
  mkdir -p $macs2_files
  ml reset
  ml MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5

  mkdir -p $macs2_files/combined
  mkdir -p $macs2_files/subnucl
  mkdir -p $macs2_files/nucl

  macs2 callpeak \
    -t $bam_files/$f \
    -f BAMPE \
    -g $tair10_size \
    -n ${f%.*} \
    -B \
    -m 5 50 \
    --fix-bimodal \
    -q 0.05 \
    --call-summits \
    --outdir $macs2_files/combined

  macs2 callpeak \
    -t $split_bam/${f%.*}.subnucl.bam \
    -f BAMPE \
    -g $tair10_size \
    -n ${f%.*}_subnucl \
    -B \
    -m 5 50 \
    --fix-bimodal \
    -q 0.05 \
    --call-summits \
    --outdir $macs2_files/subnucl

  #call peaks for the nucleosomal fraction
  macs2 callpeak \
    --broad \
    -t $split_bam/${f%.*}.nucl.bam \
    -f BAMPE \
    -g $tair10_size \
    -n ${f%.*}_nucl \
    --outdir $macs2_files/nucl \
    --nomodel \
    --shift -37 \
    --extsize 73 \
    -B \
    -q 0.05
  echo "Peak-calling with MACS2 - Done"
fi

if [ $clean -eq 1 ]; then
  echo "Cleaning up..."
  rm -rv $bed_files
  rm -rv $temp_dir
  echo "Cleaning up... - Done"
fi

echo 'Peak-calling pipeline complete!'
