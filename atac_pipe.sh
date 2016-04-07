#!/bin/bash
#PBS -P atac-seq
#PBS -N atac-seq-pipe
#PBS -J 1-69
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/$USER/log/160308_htseq-count/160308_htseq-count_^array_index^.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=4:mem=18gb


##### specify folders and variables #####
#TODO: complete
#set base dir
base_dir=/lustre/scratch/users/$USER/atac_seq
#folders for bam files
bam_files=$base_dir/bam_files
bam_files_unmapped=$bam_files/unmapped
bam_files_aligned=$bam_files/aligned
#folders for temp files
temp_dir=$base_dir/temp
#folder for fastq files
fastq_files=$base_dir/fastq
#folder for bowtie2 aligment logs
log_files=$base_dir/logs
#folder for fastqc files
fastqc_output=$base_dir/fastqc
#set bowtie2 index location
bt2_index=
#file that maps input file base names to pbs array number
mapping_file=/lustre/scratch/users/falko.hofmann/rna_seq/pbs_mapping_file.txt

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as the line number #####
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $mapping_file`
names_mapped=($input_mapper)
NAME=${names_mapped[1]}

##### load modules and assign local repos#####
#TODO: complete
module load SAMtools/1.3-goolf-1.4.10
module load BEDTools/v2.17.0-goolf-1.4.10
module load Bowtie2/2.1.0-goolf-1.4.10
module load Java/1.8.0_77
module load Python/2.7.9-foss-2015a
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5

picard=/lustre/scratch/users/$USER/software/picard/2.2.1

##### Make folders before starting the pipeline#####
#TODO: complete
mkdir -p $bam_files_unmapped
mkdir -p $bam_files_aligned
mkdir -p $fastq_files
mkdir -p $fastqc_output/${NAME}

mkdir -p $temp_dir
mkdir -p $log_files

##### Startint the ATAC-seq pipeline #####
#1.mapping
echo "1 - Start mapping part..."
#1.1 sort bam file
echo "1.1 - Name sorting bam file..."
samtools sort -n -m 4G -@ 12 -o $bam_files_unmapped/${NAME}.sorted.bam \
  $bam_files_unmapped/${NAME}.bam
echo "1.1 - Name sorting bam file... - Done"

#1.2 run fastqc
echo "1.2 - Running fastqc..."
fastqc -o $fastqc_output/${NAME} $bam_files_unmapped/${NAME}.sorted.bam
echo "1.2 - Running fastqc... - Done"

#1.3 convert bam to fq
echo "1.3 - Converting bam to fastq..."
bedtools bamtofastq -i $bam_files_unmapped/${NAME}.sorted.bam  \
  -fq $fastq_files/${NAME}.end1.fq  \
  -fq2 $fastq_files/${NAME}.end2.fq
echo "1.3 - Converting bam to fastq... - Done"

#1.4 align to genome
echo "1.4 - Starting alignment with bowtie2..."
bowtie2 --threads 12 \
  --very-sensitive \
  --maxins 2000 \
  --no-discordant \
  --no-mixed \
  -x $bt2_index \
  -1 $fastq_files/${NAME}.end1.fq \
  -2 $fastq_files/${NAME}.end2.fq \
  -S $bam_files_aligned/${NAME}.sam \
  2> $log_files/${NAME}_bt2_summary.txt
echo "1.4 - Starting alignment with bowtie2... - Done"

#1.5 sort mapped reads and convert to bam
samtools view -bS $bam_files_aligned/${NAME}.sam \
  | samtools sort -n -m 4G -@ 12 -o $bam_files_aligned/${NAME}.bam -
echo "1 - Finished mapping part..."

#2.file conversions
#TODO: implement
echo "2 - Starting post processing..."

#2.1 remove duplicates
#java -Djava.io.tmpdir=$temp_dir -jar $picard/picard.jar \
java -jar $picard/picard.jar MarkDuplicates I=$bam_files_aligned/${NAME}.bam \
  O=$bam_files_aligned/${NAME}_unique.bam M=$log_files/${NAME}_dup_metrics.txt \
  AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR

#2.2 convert to bed file
#2.3 sort bed file by mate id
#2.4 offset data
python add_offset_for_fp.py $fin $fout
#2.5 convert back to bam file


#3.analysis
#TODO: implement
#3.1 extract read length
#3.2 peakcalling
#3.2.1 macs2
#3.2.2 fseq
#3.3 estimate library complexity via preseq
#3.4 do footprinting
