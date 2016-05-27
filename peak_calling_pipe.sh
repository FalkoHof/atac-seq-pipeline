#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_peak-calling
#PBS -J 1-6
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160526_atac-seq/160526_atac-seq_^array_index^_peak_calling.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb

# === begin ENVIRONMENT SETUP ===

#pbs array for different samples
case ${PBS_ARRAY_INDEX} in
        1) NAME='bc_50_proto_160416';;
        2) NAME='bc_50_proto_160516';;
        3) NAME='mg_50_proto_160516';;
        4) NAME='tp_100_proto_160516';;
        5) NAME='test_1';;
        6) NAME='test_2';;
#esac

#set base dirs
base_dir=/lustre/scratch/users/$USER/atac_seq
pipe_dir=/lustre/scratch/users/$USER/pipes/atac-seq-pipeline

bam_files=$base_dir/bam_files/aligned/name_sorted
bam_files_coordinate_sorted=$base_dir/bam_files/aligned/coordinate_sorted

output_macs2=$base_dir/peak_calling
split_bam_files=$bam_files/split_bam

#mapping_file=$base_dir/pbs_mapping_file.txt

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as the line number #####
#input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $mapping_file`
#names_mapped=($input_mapper)
#NAME=${names_mapped[1]}

#set background control path and effective genome length for macs2
#control=$bam_files/50k_bc_background_160326.bam
effective_genome_size=1.2e8

#make folders
mkdir -p $output_macs2
mkdir -p $split_bam_files

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load BamTools/2.4.0-goolf-1.4.10
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
# === end ENVIRONMENT SETUP ===

#1.filter bam
bamtools filter -in $bam_files/${NAME}.bam -out $split_bam_files/${NAME}.subnucl.bam \
  -script $pipe_dir/bamtools_filter/filter_subnucleo.json
bamtools filter -in $bam_files/${NAME}.bam -out $split_bam_files/${NAME}.polynucl.bam \
  -script $pipe_dir/bamtools_filter/filter_nucleo.json

#samtools sort -n -m 4G -@ 8 -o $split_bam_files/${NAME}.subnucl.bam \
#  $split_bam_files/${NAME}.subnucl.bam

#samtools sort -n -m 4G -@ 8 -o \
#  $split_bam_files/${NAME}.polynucl.bam


#2. peak_call via macs2
mkdir -p $output_macs2/narrow/

macs2 callpeak -t $split_bam_files/${NAME}.subnucl.bam -f BAMPE \
  -g $effective_genome_size -n ${NAME}_subnucl --outdir $output_macs2/no_control \
  --nomodel --shift -100 --extsize 200 -B -q 0.01

#mkdir -p $output_macs2/control/broad
#macs2 callpeak --broad -t $bam_files/${NAME}.bam -c $control -f BAMPE \
#  -g $effective_genome_size -n ${NAME}_background --outdir $output_macs2/control/broad \
#  --nomodel --shift -100 --extsize 200 -B -q 0.01
