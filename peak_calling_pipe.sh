#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_peak-calling
#PBS -J 1-6
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160526_atac-seq/160526_atac-seq_^array_index^_peak_calling.log
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=8:mem=34gb

# === begin ENVIRONMENT SETUP ===
####set to 0 (false) or 1 (true) to let the repsective code block run
#1.split and filter bam files
filter_bam=0
#2. peak calling via fseq
run_fseq=1
#3. peak calling via macs2
run_macs2=0

#pbs array for different samples
case ${PBS_ARRAY_INDEX} in
        1) NAME='bc_50_proto_160416';;
        2) NAME='bc_50_proto_160516';;
        3) NAME='mg_50_proto_160516';;
        4) NAME='tp_100_proto_160516';;
        5) NAME='test_1';;
        6) NAME='test_2';;
esac

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as the line number #####
#mapping_file=$base_dir/pbs_mapping_file.txt
#input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $mapping_file`
#names_mapped=($input_mapper)
#NAME=${names_mapped[1]}

#set base dirs
base_dir=/lustre/scratch/users/$USER/atac_seq
pipe_dir=/lustre/scratch/users/$USER/pipes/atac-seq-pipeline

bam_files=$base_dir/bam_files/aligned/name_sorted
bam_files_coordinate_sorted=$base_dir/bam_files/aligned/coordinate_sorted

peak_dir=$base_dir/peak_calling
output_macs2=$peak_dir/macs2
output_fseq=$peak_dir/fseq

split_bam_files=$bam_files/split_bam
split_bam_coordiante=$bam_files_coordinate_sorted/split_bam

split_bed_files=$base_dir/bed_files/split_bed

#set background control path and effective genome length for macs2
effective_genome_size=1.2e8
#control=$bam_files/50k_bc_background_160326.bam

# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load BamTools/2.4.0-goolf-1.4.10
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
module load BEDTools/v2.17.0-goolf-1.4.10
module load Java/1.8.0_66
#download and compile fseq from git & set max heap space in fseq binary
fseq=/lustre/scratch/users/$USER/software/F-seq/dist~/fseq/bin
# === end ENVIRONMENT SETUP ===

##### Startint the ATAC-seq peak calling pipeline #####
echo 'Starting peak calling ATAC-seq pipeline for '${NAME}
#1.split and filter bam files
if [ $filter_bam -eq 1 ]; then
  echo "1 - Preparing bam files..."
  #make output folders
  mkdir -p $split_bam_files
  mkdir -p $split_bam_coordiante
  mkdir -p $bam_files_coordinate_sorted/split_bam
  #1.1 split bam files into by size ( <100, and nucleosome sizes)
  echo "1.1 - Splitting bam files to nucleosomal & nucleosomal..."
  bamtools filter -in $bam_files/${NAME}.bam \
    -out $split_bam_files/${NAME}.subnucl.bam \
    -script $pipe_dir/bamtools_filter/bamtools_subnucl.json

  bamtools filter -in $bam_files/${NAME}.bam \
    -out $split_bam_files/${NAME}.polynucl.bam \
    -script $pipe_dir/bamtools_filter/bamtools_polynucl.json
  echo "1.1 - Splitting bam files to nucleosomal & nucleosomal... - Done"
  #1.2 sort the bam files
  echo "1.2 - Sorting bam files..."
  samtools sort -m 4G -@ 8 -o $split_bam_coordiante/${NAME}.subnucl.bam \
    $split_bam_files/${NAME}.subnucl.bam

  samtools sort -m 4G -@ 8 -o $split_bam_coordiante/${NAME}.polynucl.bam \
    $split_bam_files/${NAME}.polynucl.bam
  echo "1.2 - Sorting bam files... - Done"
fi

#2. peak calling via fseq
##  use fseq to call peaks for the subnucleosomal reads, as it is sensitiver
##  than macs2
if [ $run_fseq -eq 1 ]; then
  echo "2 - Starting fseq peak calling..."
  #make ouput folders
  mkdir -p $split_bed_files
  mkdir -p $output_fseq
  #2.1 convert subnucleosomal bam to bed file for fseq
  echo "2.1 - Converting bam to bed..."
  bedtools bamtobed -i $split_bam_coordiante/${NAME}.subnucl.bam \
    > $split_bed_files/${NAME}.subnucl.bed
  echo "2.1 - Converting bam to bed... - Done"
  #2.2 call peaks with fseq
  echo "2.2 - Calling peaks with fseq"
  sh fseq -v -f 0 -of npf -t 2.0 -o $output_fseq \
    $split_bed_files/${NAME}.subnucl.bed
  echo "2.2 - Calling peaks with fseq... - Done"
fi

#3. peak calling via macs2
##  call narrow and broad peaks with macs2. The strength of macs2 is however
##  to call the broad peak regions of the polynucleosomal tracks
if [ $run_macs2 -eq 1 ]; then
  echo "3 - Starting MACS2 peak calling..."
  # make folders and set temp dir to prevent memory issues
  export TMPDIR=$WORKDIR/macs_tmp
  mkdir -p $TMPDIR
  mkdir -p $output_macs2/narrow
  mkdir -p $output_macs2/broad
  #2.1 call narrow peaks with macs2 in the subnucleosomal fraction
  echo "2.1 - Calling narrow peaks..."
  macs2 callpeak -t $split_bam_files/${NAME}.subnucl.bam -f BAMPE \
    -g $effective_genome_size -n ${NAME}_subnucl --outdir $output_macs2/narrow \
    --nomodel --shift -50 --extsize 100 -B -q 0.01
  echo "2.1 - Calling narrow peaks... - Done"
  #2.2 call broad peaks with macs2 in the polynucleosomal fraction
  echo "2.2 - Calling broad peaks..."
  macs2 callpeak --broad -t $split_bam_files/${NAME}.polynucl.bam -c $control -f BAMPE \
    -g $effective_genome_size -n ${NAME}_background --outdir $output_macs2/broad \
    --nomodel --shift -37 --extsize 73 -B -q 0.01
  echo "2.2 - Calling broad peaks... - Done"
fi
echo "Peak calling pipeline for ${NAME} complete."
