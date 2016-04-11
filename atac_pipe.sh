#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq-pipe
#PBS -J 1-7
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/160410_atac-seq/160410_atac-seq_^array_index^.log
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=64gb

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
bam_files_uniqe=$bam_files_aligned/unique
bam_files_offsetted=$bam_files_aligned/offsetted

#folders for bed files
bed_files=$base_dir/bed_files
#folders for temp files
temp_dir=$base_dir/temp
#folder for fastq files
fastq_files=$base_dir/fastq
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

####set to 0 (false) or 1 (true) to let the repsective code block run
mapping=0
post_processing=1

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
module load FastQC/0.11.5-foss-2015a
#module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5

picard=/lustre/scratch/users/$USER/software/picard/picard-tools-2.2.1

##### Make folders before starting the pipeline#####
#TODO: complete
mkdir -p $bam_files_unmapped
mkdir -p $bam_files_aligned
mkdir -p $bam_files_name_sorted
mkdir -p $bam_files_coordinate_sorted
mkdir -p $bam_files_uniqe
mkdir -p $bam_files_offsetted

mkdir -p $fastq_files
mkdir -p $fastqc_output/${NAME}

mkdir -p $temp_dir
mkdir -p $log_files

mkdir -p $bed_files

##### Startint the ATAC-seq pipeline #####
if [ $mapping -eq 1 ]; then
  echo 'Starting ATAC-seq pipeline for '${NAME}
  ##1.mapping
  echo "1 - Start mapping part..."
  #1.1 sort bam file
  echo "1.1 - Name sorting bam file..."
  samtools sort -n -m 4G -@ 8 -o $bam_files_unmapped/${NAME}.sorted.bam \
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
    #1.5 convert to bam, get mapped concordant mapped reads and and convert to bam
    samtools view -bS $temp_dir/${NAME}.sam > $temp_dir/${NAME}.bam
    samtools view -hf 0x2 $temp_dir/${NAME}.bam > $temp_dir/${NAME}_concordant_only.bam
    samtools sort -n -m 4G -@ 8 -o $bam_files_name_sorted/${NAME}.bam \
      $temp_dir/${NAME}_concordant_only.bam
    samtools sort -m 4G -@ 8 -o $bam_files_coordinate_sorted/${NAME}.bam \
      $temp_dir/${NAME}_concordant_only.bam
    #rm $temp_dir/${NAME}.sam
    echo "1 - Finished mapping part."
fi

if [ $post_processing  -eq 1 ]; then
  ##2.file conversions
  echo "2 - Starting post processing..."
  #2.1 remove duplicates
  echo "2.1 - Removing duplicates..."
  java -jar -Xmx60g $picard/picard.jar MarkDuplicates \
    I=$bam_files_coordinate_sorted/${NAME}.bam \
    O=$bam_files_uniqe/${NAME}_unique.bam M=$log_files/${NAME}_dup_metrics.txt \
    AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
  echo "2.1 - Removing duplicates... - Done"
 #2.2 convert to bed file
  echo "2.2 - Converting bam to bed..."
  bedtools bamtobed -i $bam_files_uniqe/${NAME}_unique.bam > $bed_files/${NAME}_unique.bed
  bedtools bamtobed -i $bam_files_name_sorted/${NAME}.bam > $bed_files/${NAME}.bed
  echo "2.2 - Converting bam to bed... - Done"
  #2.3 sort bed file by mate id
  sort -k4,4 -t $'\t' $bed_files/${NAME}_unique.bed > $temp_dir/${NAME}_unique_mate_sorted.bed
  sort -k4,4 -t $'\t' $bed_files/${NAME}.bed > $temp_dir/${NAME}_mate_sorted.bed

  #2.4 offset data
  echo "2.4 - Offsetting data..."
  python $pipe_dir/add_offset_for_fp.py $temp_dir/${NAME}_unique_mate_sorted.bed \
    $bed_files/${NAME}_unique_offset.bed
  python $pipe_dir/add_offset_for_fp.py $temp_dir/${NAME}_mate_sorted.bed \
    $bed_files/${NAME}_offset.bed

  rm $temp_dir/${NAME}_mate_sorted.bed
  rm $temp_dir/${NAME}_unique_mate_sorted.bed

  #2.5 convert back to bam file
  bedToBam -i $bed_files/${NAME}_unique_offset.bed -g $tair10_genome_size \
    | samtools sort -m 4G -@ 8 -o $bam_files_offsetted/${NAME}_unique_offset.bam
  bedToBam -i $bed_files/${NAME}_offset.bed -g $tair10_genome_size \
    | samtools sort -m 4G -@ 8 -o $bam_files_offsetted/${NAME}_offset.bam
  echo "2.4 - Offsetting data... - Done"

  #echo "2.5 - Extracting read length..."
  #TODO:implement for new script version
  #run script to extract the read length per mate pair
  #python extract_read_length.py.py $bed_files/${NAME}_unique.bed \
  # $log_files/${NAME}readLength_bt2.tab
  #sort -n -k2,2 -t $'\t' $WSEQ/${NAME}/bowtie/readLength_${NAME}_bt.tab > $WSEQ/${NAME}/bowtie/readLength_${NAME}_bt_sorted.tab
  #echo "2.5 - Extracting read length... - Done"
  echo "2 - Finished post processing."
fi
#TODO: implement
#3.analysis
#3.1 extract read length
#3.2 peakcalling
#3.2.1 macs2

#3.2.2 fseq
#3.3 estimate library complexity via preseq
#3.4 do footprinting
echo 'ATAC-seq pipeline complete.'
