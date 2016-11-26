#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N atac-seq_mapping
#PBS -J 1-17
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/falko.hofmann/log/161124_atac-seq/161124_atac-seq_^array_index^_mapping.log
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=16:mem=64gb

pipe_dir=/lustre/scratch/users/$USER/pipelines/atac-seq-pipeline
#set ouput base dir
base_dir=/lustre/scratch/users/$USER/atac-seq
#folder for bowtie indices
bt_2_index=/lustre/scratch/users/$USER/indices/bowtie2/Col-0
#location of the mapping file for the array job
pbs_mapping_file=$pipe_dir/pbs_mapping_file.txt
#super folder of the temp dir, script will create subfolders with $sample_name
temp_dir_base=$base_dir/temp

#convert bam to fastq
convert_bam=1
#delete unecessary files
clean=1
#specfify if alignment should be run
align=1
#specify number of threads
threads=16 #set this to the number of available cores

##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as line number
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $pbs_mapping_file` #read mapping file
names_mapped=($input_mapper)
sample_dir=${names_mapped[1]} # get the sample dir
sample_name=`basename $sample_dir` #get the base name of the dir as sample name

temp_dir=$temp_dir_base/$sample_name

mkdir -p $temp_dir

#print some output for logging
echo '#########################################################################'
echo 'Starting ATAC-seq pipeline for: ' $sample_name
echo 'Sample directory: ' $sample_dir
echo 'Mapping file: ' $pbs_mapping_file
#IDEA: maybe add some diagnostic output about bt1 and bt2
echo '#########################################################################'

#some error handling function
function error_exit
{
  echo "$1" 1>&2
  exit 1
}

#load some modules...
ml SAMtools/1.3.1-foss-2015b
ml BEDTools/2.26.0-foss-2015b
ml picard/2.3.0
ml Bowtie2/2.2.7-foss-2015b

#specify some variables for adaptor trimming...
skewer=/lustre/scratch/users/falko.hofmann/software/skewer/skewer
adaptor_1=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
adaptor_2=CTGTCTCTTATACACATCTGACGCTGCCGACGA

# get all bam files in folder
f=($(ls $sample_dir | grep -e ".bam"))
#throw error if more or less than 1 file is present
if [[ "${#f[@]}" -ne "1" ]]; then
  error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
fi

fastq_dir=$sample_dir/fastq
mkdir -p $fastq_dir
mkdir -p $sample_dir/logs

if [ $convert_bam -eq 1 ]; then
  echo "Converting bam to fastq..."

  mkdir -p $sample_dir/fastq
  #do the sorting....
  samtools sort -n -m 3G -@ $threads -o $fastq_dir/${f%.*}.sorted.bam \
    $sample_dir/$f

  bedtools bamtofastq -i $fastq_dir/${f%.*}.sorted.bam  \
    -fq $fastq_dir/${f%.*}.1.fq  \
    -fq2 $fastq_dir/${f%.*}.2.fq
  echo "Converting bam to fastq... - Done"

  echo "Trimming adaptors with skewer..."
  #run skewer
  $skewer \
    -x $adaptor_1 \
    -y $adaptor_2 \
    $fastq_dir/${f%.*}.1.fq $fastq_dir/${f%.*}.2.fq \
    -o $fastq_dir/${f%.*} \
    --quiet

  mv -v $fastq_dir/*.log $sample_dir/logs
  echo "Trimming adaptors with skewer... - Done"
fi

if [ $align -eq 1 ]; then

    echo "Aligning with bowtie2..."
    mkdir -p $sample_dir/alignments

    bowtie2 --threads $threads \
      --very-sensitive \
      --maxins 2000 \
      -x $bt_2_index \
      -1 $fastq_dir/${f%.*}-trimmed-pair1.fastq \
      -2 $fastq_dir/${f%.*}-trimmed-pair2.fastq \
      -S $sample_dir/alignments/${f%.*}.sam \
      2> $sample_dir/logs/${f%.*}_summary.txt
    echo "Aligning with bowtie2... - Done"

    echo "Converting to bam..."
    samtools sort -m 3G -@ $threads $sample_dir/alignments/${f%.*}.sam \
      -o $sample_dir/alignments/${f%.*}.sorted.bam
    echo "Converting to bam... - Done"

    echo "Removing duplicates..."
    java -jar -Xmx60g ${EBROOTPICARD}/picard.jar MarkDuplicates \
      I=$sample_dir/alignments/${f%.*}.sorted.bam \
      O=$sample_dir/alignments/${f%.*}.no_dups.bam \
      M=$sample_dir/logs/${f%.*}dup_metrics.txt \
      AS=true REMOVE_DUPLICATES=true TMP_DIR=$temp_dir
    echo "Removing duplicates... - Done"

    echo "Quality filtering..."
    samtools view -bhf 0x2 -q 30 $sample_dir/alignments/${f%.*}.no_dups.bam | \
      samtools sort -m 3G -@ $threads - -o $sample_dir/alignments/${f%.*}.bam

    samtools view -bhf 0x2 -q 40 $sample_dir/alignments/${f%.*}.no_dups.bam | \
      samtools sort -m 3G -@ $threads - -o $sample_dir/alignments/${f%.*}.unique.bam
    echo "Quality filtering... - Done"

    echo "Indexing bam files..."
    #indes the bam file for downstream applications..
    samtools index $sample_dir/alignments/${f%.*}.bam
    samtools index $sample_dir/alignments/${f%.*}.unique.bam
    echo "Indexing bam files... - Done"
fi

if [ $clean -eq 1 ]; then
  echo "Cleaning up..."
  rm -v $sample_dir/alignments/${f%.*}.sam
  rm -v $sample_dir/alignments/${f%.*}.sorted.bam
  rm -v $sample_dir/alignments/${f%.*}.no_dups.bam
  rm -rv $sample_dir/fastq
  rm -rv $temp_dir
  echo "Cleaning up... - Done"
fi
