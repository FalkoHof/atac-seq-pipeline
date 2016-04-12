#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_peak-calling
#PBS -J 1-6
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160410_atac-seq/160410_atac-seq_^array_index^_peak_calling.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb

effective_genome_size=116e6
control=50k_bc_background_160326
output_dir=

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load BamTools/2.4.0-goolf-1.4.10
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
# === end ENVIRONMENT SETUP ===

case ${PBS_ARRAY_INDEX} in
        1) NAME='100k_bc_160326';;
        2) NAME='75k_bc_160326';;
        3) NAME='50k_bc_160326';;
        4) NAME='50k_bc_160328';;
        5) NAME='test_1';;
        6) NAME='test_2';;
esac


#2. peak_call via macs2
mkdir -p $output_dir/no_control

macs2 callpeak -t ${NAME}.subnucl.bam -f BAMPE -g $effective_genome_size\
  -n ${NAME}_offset --outdir $output_dir/no_control \
  --nomodel --shift -100 --extsize 200 -B -q 0.01


mkdir -p $output_dir/control/narrow
macs2 callpeak -t ${NAME}.bam -c $control -f BAMPE -g $effective_genome_size\
    -n ${NAME}_offset --outdir $output_dir/no_model \
    --nomodel --shift -100 --extsize 200 -B -q 0.01


mkdir -p $output_dir/broad
macs2 callpeak --broad -t ${NAME}.bam -c $control -f BAM -g $effective_genome_size \
  -n ${NAME} --broad-cutoff 0.1 –outdir $output_dir/with_model

macs2 callpeak --broad -t ${NAME}.bam -c $control -f BAMPE -g $effective_genome_size\
      -n ${NAME}_offset --outdir $output_dir/no_model \
      --nomodel --shift -100 --extsize 200 -B -q 0.01

#2. peak_call via fseq
fseq=/lustre/scratch/users/$USER/software/F-seq/dist~/fseq/bin/fseq

#2.1 split file into multiple bed files

#2.2 build wig files



fseq
