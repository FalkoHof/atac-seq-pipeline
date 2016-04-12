#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_mapping
#PBS -J 1-7
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160411_atac-seq/160411_atac-seq_^array_index^_mapping.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb


# === begin ENVIRONMENT SETUP ===
##### specify folders and variables #####
#set base dir
pipe_dir=/lustre/scratch/users/$USER/pipes/atac-seq-pipeline
base_dir=/lustre/scratch/users/$USER/atac_seq
#folders for bam files
bam_files=$base_dir/bam_files
bam_files_unmapped=$bam_files/unmapped
bam_files_aligned=$bam_files/aligned
bam_files_name_sorted=$bam_files_aligned/name_sorted
bam_files_coordinate_sorted=$bam_files_aligned/coordinate_sorted
bam_files_offsetted=$bam_files_aligned/offsetted
#folders for bed files
bed_files=$base_dir/bed_files
#folders for temp files
temp_dir=$base_dir/temp
#folder for fastq files
fastq_files=$base_dir/fastq
fastq_pre=$fastq_files/pre_alignment
fastq_post=$fastq_files/post_alignment
#folder for bowtie2 aligment logs
log_files=$base_dir/logs
#folder for fastqc files
fastqc_output=$base_dir/fastqc
#set bowtie2 index location
bt2_index=/lustre/scratch/users/$USER/indices/bowtie2/Col/Col
#genome sizes for converting back to bam
tair10_genome_size=$pipe_dir/chromLength.txt
#file that maps input file base names to pbs array number
mapping_file=$base_dir/pbs_mapping_file.txt
read_length_dir=$base_dir/read_length

####set to 0 (false) or 1 (true) to let the repsective code block run
#1.2 align to genome
mapping=0
#2. convert sam to bam, keep only concordant mapped reads
converting=1
#3. do some post processing
remove_duplicates=1
post_processing=1
#4. delete unecessary files from temp_dir
clean=0
##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as the line number #####
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $mapping_file`
names_mapped=($input_mapper)
NAME=${names_mapped[1]}

##### load modules and assign local repos#####
module load SAMtools/1.3-goolf-1.4.10
module load BEDTools/v2.17.0-goolf-1.4.10
module load Bowtie2/2.1.0-goolf-1.4.10
module load Java/1.8.0_77
module load Python/2.7.9-foss-2015a
module load FastQC/0.11.5-foss-2015a

picard=/lustre/scratch/users/$USER/software/picard/picard-tools-2.2.1

##### Make folders before starting the pipeline#####
mkdir -p $bam_files_unmapped
mkdir -p $bam_files_aligned
mkdir -p $bam_files_name_sorted
mkdir -p $bam_files_coordinate_sorted
mkdir -p $bam_files_offsetted

mkdir -p $fastq_files
#mkdir -p $fastqc_output/${NAME}
mkdir -p $fastq_pre/${NAME}
mkdir -p $fastq_post/${NAME}

mkdir -p $temp_dir
mkdir -p $log_files
mkdir -p $read_length_dir

mkdir -p $bed_files
# === end ENVIRONMENT SETUP ===

##### Startint the ATAC-seq pipeline #####
echo 'Starting ATAC-seq pipeline for '${NAME}
##1.mapping
if [ $mapping -eq 1 ]; then
  echo "1 - Start mapping part..."
  #1.1 sort bam file
  echo "1.1 - Name sorting bam file..."
  samtools sort -n -m 4G -@ 8 -o $bam_files_unmapped/${NAME}.sorted.bam \
  $bam_files_unmapped/${NAME}.bam
  echo "1.1 - Name sorting bam file... - Done"
  #1.2 run fastqc
  echo "1.2 - Running fastqc..."
  fastqc -o $fastqc_pre/${NAME} $bam_files_unmapped/${NAME}.sorted.bam
  echo "1.2 - Running fastqc... - Done"
  #1.3 convert bam to fq
  echo "1.3 - Converting bam to fastq..."
  bedtools bamtofastq -i $bam_files_unmapped/${NAME}.sorted.bam  \
    -fq $fastq_files/${NAME}.end1.fq  \
    -fq2 $fastq_files/${NAME}.end2.fq
  echo "1.3 - Converting bam to fastq... - Done"
  echo "1.4 - Starting alignment with bowtie2..."
  bowtie2 --threads 8 \
    --very-sensitive \
    --maxins 2000 \
    --no-discordant \
    --no-mixed \
    -x $bt2_index \
    -1 $fastq_files/${NAME}.end1.fq \
    -2 $fastq_files/${NAME}.end2.fq \
    -S $temp_dir/${NAME}.sam \
    2> $log_files/${NAME}_bt2_summary.txt
    echo "1.4 - Starting alignment with bowtie2... - Done"
fi
if [ $converting -eq 1 ]; then
    #1.5 convert to bam, get mapped concordant mapped reads and and convert to bam
    samtools view -bS $temp_dir/${NAME}.sam > $temp_dir/${NAME}.bam
    samtools view -bhf 0x2 $temp_dir/${NAME}.bam > $temp_dir/${NAME}_concordant_only.bam
    samtools sort -n -m 4G -@ 8 -o $bam_files_name_sorted/${NAME}.bam \
      $temp_dir/${NAME}_concordant_only.bam
    samtools sort -m 4G -@ 8 -o $bam_files_coordinate_sorted/${NAME}.bam \
      $temp_dir/${NAME}_concordant_only.bam
    #rm $temp_dir/${NAME}.sam
    #fastqc -o $fastqc_post/${NAME} $bam_files_coordinate_sorted/${NAME}.bam
    echo "1 - Finished mapping part."
fi

##2.file conversions, duplicates removal and offsetting
if [ $remove_duplicates  -eq 1 ]; then
  echo "2 - Starting post processing..."
  #2.1 remove duplicates
  echo "2.1 - Removing duplicates..."
  java -jar -Xmx60g $picard/picard.jar MarkDuplicates \
    I=$bam_files_coordinate_sorted/${NAME}.bam \
    O=$bam_files_coordinate_sorted/${NAME}_unique.bam M=$log_files/${NAME}_dup_metrics.txt \
    AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
  echo "2.1 - Removing duplicates... - Done"
  #2.2 sort bam file
  echo "2.2 - Sorting unique reads bam..."
  samtools sort -n -m 4G -@ 8 -o $bam_files_name_sorted/${NAME}_unique.bam \
    $bam_files_coordinate_sorted/${NAME}_unique.bam
  echo "2.2 - Sorting unique reads bam - Done..."
fi

if [ $post_processing  -eq 1 ]; then
  echo "2.3 - Converting bam to bed..."
  bedtools bamtobed -i $bam_files_name_sorted/${NAME}_unique.bam > $bed_files/${NAME}_unique.bed
  bedtools bamtobed -i $bam_files_name_sorted/${NAME}.bam > $bed_files/${NAME}.bed
  echo "2.3 - Converting bam to bed... - Done"

  #2.4 offset data
  echo "2.4 - Offsetting data..."
  python $pipe_dir/add_offset_for_fp.py $bed_files/${NAME}_unique.bed \
    $bed_files/${NAME}_unique_offset.bed
  python $pipe_dir/add_offset_for_fp.py $bed_files/${NAME}.bed \
    $bed_files/${NAME}_offset.bed

  #2.5 convert back to bam file
  bedToBam -i $bed_files/${NAME}_unique_offset.bed -g $tair10_genome_size \
    | samtools sort -m 4G -@ 8 -o $bam_files_offsetted/${NAME}_unique_offset.bam
  bedToBam -i $bed_files/${NAME}_offset.bed -g $tair10_genome_size \
    | samtools sort -m 4G -@ 8 -o $bam_files_offsetted/${NAME}_offset.bam
  echo "2.4 - Offsetting data... - Done"

  echo "2.5 - Extracting read lenght..."
  python $pipe_dir/extract_read_length.py -g -v -o $read_length_dir \
    $bed_files/${NAME}.bed
  echo "2.5 - Extracting read length... - Done"
  echo "2 - Finished post processing."
fi
#3. clean up and delete unecessary files
if [ $clean  -eq 1 ]; then
  echo "Cleaning up..."
  rm $temp_dir/${NAME}.sam
  rm $temp_dir/*.bam
  rm $temp_dir/*.bed
  rm $fastq_files/${NAME}.end1.fq
  rm $fastq_files/${NAME}.end2.fq
  echo "Cleaning up... - Done"
fi
echo 'ATAC-seq pipeline complete.'
