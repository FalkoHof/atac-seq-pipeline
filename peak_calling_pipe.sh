#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_peak-calling
#PBS -J 1-6
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160526_atac-seq/160526_atac-seq_^array_index^_peak_calling.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=34gb

# === begin ENVIRONMENT SETUP ===
pipe_dir=/lustre/scratch/users/$USER/pipelines/atac-seq-pipeline
#set ouput base dir
base_dir=/lustre/scratch/users/$USER/atac-seq
#TODO: specify paths here...
threads=8
#download and compile fseq from git & set max heap space in fseq exectuable
fseq=/lustre/scratch/users/$USER/software/F-seq/dist~/fseq/bin/fseq
#location of the mapping file for the array job
pbs_mapping_file=$pipe_dir/pbs_mapping_file.txt
#super folder of the temp dir, script will create subfolders with $sample_name
temp_dir_base=$base_dir/temp
#json filters for bamtools
subnucl_filter=$pipe_dir/bamtools_filter/bamtools_subnucl.json
nucl_filter=$pipe_dir/bamtools_filter/bamtools_polynucl.json

#TODO: clean up modules
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load BamTools/2.4.0-goolf-1.4.10
module load BEDTools/v2.17.0-goolf-1.4.10
module load Java/1.8.0_66

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as line number
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $pbs_mapping_file` #read mapping file
names_mapped=($input_mapper)
sample_dir=${names_mapped[1]} # get the sample dir
sample_name=`basename $sample_dir` #get the base name of the dir as sample name

temp_dir=$temp_dir_base/sample_name

mkdir -p $temp_dir

#print some output for logging
echo '#########################################################################'
echo 'Starting ATAC-seq peak-calling pipeline for: ' $sample_name
echo 'Sample directory: ' $sample_dir
echo 'Mapping file: ' $pbs_mapping_file
#IDEA: maybe add some diagnostic output the chosen modes.
echo '#########################################################################'

#some error handling function
function error_exit
{
  echo "$1" 1>&2
  exit 1
}

#TODO and conditional initilization..
aligner_dirs=()
split_dirs=()
bed_dirs=()
peaks_dirs=()

if [ $bowtie_1 -eq 1 ]; then
  bt_1_files=$sample_dir/bowtie
  bt_1_split=$bt_1_files/split_bam
  bt_1_bed=$bt_1_files/bed_files
  bt_1_peaks=$bt_1_files/peak_calling

  mkdir -p bt_1_split
  mkdir -p bt_1_peaks
  mkdir -p bt_1_bed

  aligner_dirs+=(bt_1_files)
  split_dirs+=(bt_1_split)
  bed_dirs+=(bt_1_bed)
  peaks_dirs+=(bt_1_peaks)
fi

if [ $bowtie_2 -eq 1 ]; then
  bt_2_files=$sample_dir/bowtie2
  bt_2_split=$bt_2_files_split/split_bam
  bt_2_bed=$bt_1_files/bed_files
  bt_2_peaks=$bt_2_files/peak_calling

  mkdir -p bt_2_split
  mkdir -p bt_2_peaks
  mkdir -p bt_2_bed

  aligner_dirs+=(bt_2_files)
  split_dirs+=(bt_2_split)
  bed_dirs+=(bt_2_bed)
  peaks_dirs+=(bt_2_peaks)
fi
# === end ENVIRONMENT SETUP ===
##### Starting the ATAC-seq peak calling pipeline #####
for (( i = 0 ; i < ${#aligner_dirs[@]} ; i++ )); do

  f=($(ls ${aligner_dirs[$i]} | grep -e ".sorted.bam"))
  
  if [[ "${#f[@]}" -ne "1" ]]; then
    error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
  fi

  #FIXME: switch to a loop based system..
  echo "Splitting bam files..."
  #get the subnucleosomal reads and sort them
  bamtools filter -in ${aligner_dirs[$i]}/$f -script $subnucl_filter | \
    samtools sort -m 3G -@ $threads - -o ${split_dirs[$i]}/${f%.*}.subnucl.bam
  samtools index ${split_dirs[$i]}/${f%.*}.subnucl.bam
  #get the nucleosomal reads and sort them
  bamtools filter -in ${aligner_dirs[$i]}/$f -script $nucl_filter | \
    samtools sort -m 3G -@ $threads - -o ${split_dirs[$i]}/${f%.*}.nucl.bam
  samtools index ${split_dirs[$i]}/${f%.*}.nucl.bam
  #calculate some stats...
  reads_subnucl=$(samtools view -c -f 1 ${split_dirs[$i]}/${f%.*}.subnucl.bam)
  reads_nucl=$(samtools view -c -f 1 ${split_dirs[$i]}/${f%.*}.nucl.bam)
  enrichment_subncl = $((reads_subnucl / reads_nucl))
  #some string assignments
  read_stats="Subnucleosomal reads: $reads_subnucl\n"
  read_stats+="Nucleosomal reads: $reads_nucl\n"
  read_stats+="Ratio: $enrichment_subncl"
  #print and save them
  printf $read_stats | tee ${split_dirs[$i]}/ratios.txt
  echo "Splitting bam files... - Done"

  echo "Converting bam files to bed..."
  #convert to bed, keep reads from nuclear chromosomes, then sort and store them
  bedtools bamtobed -i ${aligner_dirs[$i]}/$f | grep "^Ath_chr[1-5]" |\
    sort -k1,1V -k2,2n > ${bed_dirs[$i]}/${f%.*}.bed
  bedtools bamtobed -i ${split_dirs[$i]}/${f%.*}.subnucl.bam | \
    grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n > ${bed_dirs[$i]}/${f%.*}.subnucl.bed
  bedtools bamtobed -i ${split_dirs[$i]}/${f%.*}.nucl.bam | \
    grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n > ${bed_dirs[$i]}/${f%.*}.nucl.bed
  echo "Converting bam files to bed... - Done"
  #make some output folders
  mkdir -p ${peaks_dirs[$i]}/combined
  mkdir -p ${peaks_dirs[$i]}/subnucl
  mkdir -p ${peaks_dirs[$i]}/nucl

  echo "Peak-calling with f-seq..."
  #run fseq on the different files...
  sh $fseq -v -f 0 -of npf -t 1.0 -o ${peaks_dirs[$i]}/combined \
    ${bed_dirs[$i]}/${f%.*}.bed
  sh $fseq -v -f 0 -of npf -t 1.0 -o ${peaks_dirs[$i]}/subnucl \
    ${bed_dirs[$i]}/${f%.*}.subnucl.bed
  sh $fseq -v -f 0 -of npf -t 1.0 -o ${peaks_dirs[$i]}/nucl \
    ${bed_dirs[$i]}/${f%.*}.nucl.bed
  echo "Peak-calling with f-seq... - Done"

  echo "Merging f-seq files..."
  #put the fseq output file names in some arrays
  fseq_combined=($(ls ${peaks_dirs[$i]}/combined | grep -e "Ath_chr[1-5].npf"))
  fseq_subnucl=($(ls ${peaks_dirs[$i]}/subnucl | grep -e "Ath_chr[1-5].npf"))
  fseq_nucl=($(ls ${peaks_dirs[$i]}/nucl | grep -e "Ath_chr[1-5].npf"))
  #concatenate them to one file and sort
  cat ${fseq_combined[@]} | sort -k 1,1 -k2,2n > \
    ${peaks_dirs[$i]}/combined/${f%.*}_fseq.np
  cat ${fseq_subnucl[@]} | sort -k 1,1 -k2,2n > \
    ${peaks_dirs[$i]}/subnucl/${f%.*}_fseq.np
  cat ${fseq_nucl[@]} | sort -k 1,1 -k2,2n > \
    ${peaks_dirs[$i]}/nucl/${f%.*}_fseq.np
  echo "Merging f-seq files... - Done"

  #TODO: implement removal of the old per chromosome files.
  #rm -vr ${fseq_combined[@]}
  #rm -vr ${fseq_subnucl[@]}
  #rm -vr ${fseq_nucl[@]}
done
#TODO: implement clean loop that removes or zips the bed files.
