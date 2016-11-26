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

# #TODO: set to the appropriate file
# macs2_control=
# effective_genome_size=1.2e8
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

temp_dir=$temp_dir_base/sample_name

mkdir -p $temp_dir

#print some output for logging
echo '#########################################################################'
echo 'Starting ATAC-seq peak-calling pipeline for: ' $sample_name
echo 'Sample directory: ' $sample_dir
echo 'Mapping file: ' $pbs_mapping_file
#IDEA: maybe add some diagnostic output the chosen modes.
echo '#########################################################################'

# Load the required modules each loop due to deeptools loading differnt libs
ml SAMtools/1.3.1-foss-2016a
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

mkdir -p $split_bam
mkdir -p $bed_files
mkdir -p $fseq_files
mkdir -p $wig_files

f=($(ls $bam_files | grep -e "unique\.bam$"))

if [[ "${#f[@]}" -ne "1" ]]; then
  error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
fi

echo "Splitting bam files..."
#get the subnucleosomal reads and sort them
bamtools filter -in $bam_files/$f -script $subnucl_filter | \
  samtools sort -m 3G -@ $threads - -o $split_bam/${f%.*}.subnucl.bam
samtools index $split_bam/${f%.*}.subnucl.bam
#get the nucleosomal reads and sort them
bamtools filter -in $bam_files/$f -script $nucl_filter | \
  samtools sort -m 3G -@ $threads - -o $split_bam/${f%.*}.nucl.bam
samtools index $split_bam/${f%.*}.nucl.bam
echo "Splitting bam files... - Done"

echo "Converting bam files to bed..."
#convert to bed, keep reads from nuclear chromosomes, then sort and store them
bedtools bamtobed -i $bam_files/$f | grep "^Ath_chr[1-5]" |\
  sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.bed
bedtools bamtobed -i $bam_files/${f%.*}.subnucl.bam | \
  grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.subnucl.bed
bedtools bamtobed -i $bam_files/${f%.*}.nucl.bam | \
  grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > $bed_files/${f%.*}.nucl.bed
echo "Converting bam files to bed... - Done"
#make some output folders
mkdir -p $fseq_files/combined
mkdir -p $fseq_files/subnucl
mkdir -p $fseq_files/nucl

echo "Peak-calling with f-seq..."
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
fseq_combined=($(ls ${peaks_dirs[$i]}/combined | grep -e "Ath_chr[1-5].npf"))
fseq_subnucl=($(ls ${peaks_dirs[$i]}/subnucl | grep -e "Ath_chr[1-5].npf"))
fseq_nucl=($(ls ${peaks_dirs[$i]}/nucl | grep -e "Ath_chr[1-5].npf"))
#concatenate them to one file and sort
cd ${peaks_dirs[$i]}/combined
cat ${fseq_combined[@]} | sort -k 1,1 -k2,2n > \
  ${peaks_dirs[$i]}/combined/${f%.*}_combined_fseq.npf

rm -v ${fseq_combined[@]}

cd ${peaks_dirs[$i]}/subnucl
cat ${fseq_subnucl[@]} | sort -k 1,1 -k2,2n > \
  ${peaks_dirs[$i]}/subnucl/${f%.*}_subnucl_fseq.npf
rm -v ${fseq_subnucl[@]}

cd ${peaks_dirs[$i]}/nucl
cat ${fseq_nucl[@]} | sort -k 1,1 -k2,2n > \
  ${peaks_dirs[$i]}/nucl/${f%.*}_nucl_fseq.npf
rm -v ${fseq_nucl[@]}
echo "Merging f-seq files... - Done"

echo "Creating normalized bigwig files..."
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
  -b $bam_files/${f%.*}.subnucl.bam  \
  -o $wig_files/${f%.*}.subnucl.bw \
  --binSize=1 \
  --normalizeTo1x $tair10_size \
  --ignoreDuplicates \
  --numberOfProcessors=$threads \
  --ignoreForNormalization Ath_chrm Ath_chrc

bamCoverage \
  -b $bam_files/${f%.*}.nucl.bam  \
  -o $wig_files/${f%.*}.nucl.bw \
  --binSize=1 \
  --normalizeTo1x $tair10_size \
  --ignoreDuplicates \
  --numberOfProcessors=$threads \
  --ignoreForNormalization Ath_chrm Ath_chrc
echo "Creating normalized bigwig files... - Done"


echo "Peak-calling with MACS2"
ml MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5

macs2 callpeak \
  -t $bam_files/$f \
  -f BAMPE \
  -g $effective_genome_size \
  -n ${f%.*} \
  -B \
  -m 5 50 \
  --fix-bimodal \
  -q 0.05 \
  --call-summits \
  --outdir $macs2_files/combined


macs2 callpeak \
  -t $bam_files/${f%.*}.subnucl.bam \
  -f BAMPE \
  -g $effective_genome_size \
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
  -t $bam_files/${f%.*}.nucl.bam \
  -f BAMPE \
  -g $effective_genome_size \
  -n ${f%.*}_nucl \
  --outdir $macs2_files/nucl \
  --nomodel \
  --shift -37 \
  --extsize 73 \
  -B \
  -q 0.05
echo "Peak-calling with MACS2 - Done"


if [ $clean -eq 1 ]; then
  echo "Cleaning up..."
  rm -rv $bed_files
  rm -rv $temp_dir
  echo "Cleaning up... - Done"
fi



# aligner_dirs=()
# split_dirs=()
# bed_dirs=()
# peaks_dirs=()
# bw_dirs=()
#some if statements that allow control over what is run and make sure the
#right stuff is added to the appropriate arrays
# if [ $bowtie_1 -eq 1 ]; then
#   bt_1_files=$sample_dir/bowtie
#   bt_1_split=$bt_1_files/split_bam
#   bt_1_bed=$bt_1_files/bed_files
#   bt_1_peaks=$bt_1_files/peak_calling
#   bt_1_wg=$bt_1_files/wig_files
#
#   mkdir -p $bt_1_split
#   mkdir -p $bt_1_peaks
#   mkdir -p $bt_1_bed
#   mkdir -p $bt_1_wg
#
#   aligner_dirs+=($bt_1_files)
#   split_dirs+=($bt_1_split)
#   bed_dirs+=($bt_1_bed)
#   peaks_dirs+=($bt_1_peaks)
#   bw_dirs+=($bt_1_wg)
# fi
#
# if [ $bowtie_2 -eq 1 ]; then
#   bt_2_files=$sample_dir/bowtie2
#   bt_2_split=$bt_2_files/split_bam
#   bt_2_bed=$bt_2_files/bed_files
#   bt_2_peaks=$bt_2_files/peak_calling
#   bt_2_wg=$bt_2_files/wig_files
#
#   mkdir -p $bt_2_split
#   mkdir -p $bt_2_peaks
#   mkdir -p $bt_2_bed
#   mkdir -p $bt_2_wg
#
#   aligner_dirs+=($bt_2_files)
#   split_dirs+=($bt_2_split)
#   bed_dirs+=($bt_2_bed)
#   peaks_dirs+=($bt_2_peaks)
#   bw_dirs+=($bt_2_wg)
# fi
# # === end ENVIRONMENT SETUP ===
# ##### Starting the ATAC-seq peak calling pipeline #####
# for (( i = 0 ; i < ${#aligner_dirs[@]} ; i++ )); do
#
#   # Load the required modules each loop due to deeptools loading differnt libs
#   ml SAMtools/1.3.1-foss-2016a
#   ml BamTools/2.4.0-foss-2016a
#   ml BEDTools/2.26.0-foss-2016a
#   ml Java/1.8.0_66
#
#   echo "Running pipeline for aligner: ${aligner_dirs[$i]##/*/}"
#
#   f=($(ls ${aligner_dirs[$i]} | grep -e "\.bam$"))
#
#   if [[ "${#f[@]}" -ne "1" ]]; then
#     error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
#   fi
#
#   #FIXME: switch to a loop based system..
#   echo "Splitting bam files..."
#   #get the subnucleosomal reads and sort them
#   bamtools filter -in ${aligner_dirs[$i]}/$f -script $subnucl_filter | \
#     samtools sort -m 3G -@ $threads - -o ${split_dirs[$i]}/${f%.*}.subnucl.bam
#   samtools index ${split_dirs[$i]}/${f%.*}.subnucl.bam
#   #get the nucleosomal reads and sort them
#   bamtools filter -in ${aligner_dirs[$i]}/$f -script $nucl_filter | \
#     samtools sort -m 3G -@ $threads - -o ${split_dirs[$i]}/${f%.*}.nucl.bam
#   samtools index ${split_dirs[$i]}/${f%.*}.nucl.bam
#   #calculate some stats...
#   #reads_subnucl=$(samtools view -c -f 1 ${split_dirs[$i]}/${f%.*}.subnucl.bam)
#   #reads_nucl=$(samtools view -c -f 1 ${split_dirs[$i]}/${f%.*}.nucl.bam)
#   #enrichment_subncl=$(($reads_subnucl / $reads_nucl))
#   #some string assignments
#   #read_stats="Subnucleosomal reads: $reads_subnucl\n"
#   #read_stats+="Nucleosomal reads: $reads_nucl\n"
#   #read_stats+="Ratio: $enrichment_subncl"
#   #print and save them
#   #printf $read_stats | tee ${split_dirs[$i]}/ratios.txt
#   echo "Splitting bam files... - Done"
#
#   echo "Converting bam files to bed..."
#   #convert to bed, keep reads from nuclear chromosomes, then sort and store them
#   bedtools bamtobed -i ${aligner_dirs[$i]}/$f | grep "^Ath_chr[1-5]" |\
#     sort -k1,1V -k2,2n -T $temp_dir > ${bed_dirs[$i]}/${f%.*}.bed
#   bedtools bamtobed -i ${split_dirs[$i]}/${f%.*}.subnucl.bam | \
#     grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > ${bed_dirs[$i]}/${f%.*}.subnucl.bed
#   bedtools bamtobed -i ${split_dirs[$i]}/${f%.*}.nucl.bam | \
#     grep "^Ath_chr[1-5]" | sort -k1,1V -k2,2n -T $temp_dir > ${bed_dirs[$i]}/${f%.*}.nucl.bed
#   echo "Converting bam files to bed... - Done"
#   #make some output folders
#   mkdir -p ${peaks_dirs[$i]}/combined
#   mkdir -p ${peaks_dirs[$i]}/subnucl
#   mkdir -p ${peaks_dirs[$i]}/nucl
#
#   echo "Peak-calling with f-seq..."
#   #run fseq on the different files...
#   sh $fseq -v -f 0 -of npf -t 2.0 -o ${peaks_dirs[$i]}/combined \
#     ${bed_dirs[$i]}/${f%.*}.bed
#   sh $fseq -v -f 0 -of npf -t 2.0 -o ${peaks_dirs[$i]}/subnucl \
#     ${bed_dirs[$i]}/${f%.*}.subnucl.bed
#   sh $fseq -v -f 0 -of npf -t 2.0 -o ${peaks_dirs[$i]}/nucl \
#     ${bed_dirs[$i]}/${f%.*}.nucl.bed
#   echo "Peak-calling with f-seq... - Done"
#
#   echo "Merging f-seq files..."
#   #put the fseq output file names in some arrays
#   fseq_combined=($(ls ${peaks_dirs[$i]}/combined | grep -e "Ath_chr[1-5].npf"))
#   fseq_subnucl=($(ls ${peaks_dirs[$i]}/subnucl | grep -e "Ath_chr[1-5].npf"))
#   fseq_nucl=($(ls ${peaks_dirs[$i]}/nucl | grep -e "Ath_chr[1-5].npf"))
#   #concatenate them to one file and sort
#   cd ${peaks_dirs[$i]}/combined
#   cat ${fseq_combined[@]} | sort -k 1,1 -k2,2n > \
#     ${peaks_dirs[$i]}/combined/${f%.*}_combined_fseq.npf
#     #${peaks_dirs[$i]}/combined/${f%.*}_fseq_combined.npf
#
#   rm -v ${fseq_combined[@]}
#
#   cd ${peaks_dirs[$i]}/subnucl
#   cat ${fseq_subnucl[@]} | sort -k 1,1 -k2,2n > \
#     ${peaks_dirs[$i]}/subnucl/${f%.*}_subnucl_fseq.npf
#   rm -v ${fseq_subnucl[@]}
#
#   cd ${peaks_dirs[$i]}/nucl
#   cat ${fseq_nucl[@]} | sort -k 1,1 -k2,2n > \
#     ${peaks_dirs[$i]}/nucl/${f%.*}_nucl_fseq.npf
#   rm -v ${fseq_nucl[@]}
#   echo "Merging f-seq files... - Done"
#
#   echo "Creating normalized bigwig files..."
#   #load module
#   ml deepTools/2.2.4-foss-2015a-Python-2.7.9
#   ml OpenSSL/1.0.1p-foss-2015a
#
#   bamCoverage \
#     -b ${aligner_dirs[$i]}/$f \
#     -o ${bw_dirs[$i]}/${f%.*}.bw \
#     --binSize=1 \
#     --normalizeTo1x $tair10_size \
#     --ignoreDuplicates \
#     --numberOfProcessors=$threads \
#     --ignoreForNormalization Ath_chrm Ath_chrc
#
#   bamCoverage \
#     -b ${split_dirs[$i]}/${f%.*}.subnucl.bam  \
#     -o ${bw_dirs[$i]}/${f%.*}.subnucl.bw \
#     --binSize=1 \
#     --normalizeTo1x $tair10_size \
#     --ignoreDuplicates \
#     --numberOfProcessors=$threads \
#     --ignoreForNormalization Ath_chrm Ath_chrc
#
#   bamCoverage \
#     -b ${split_dirs[$i]}/${f%.*}.nucl.bam  \
#     -o ${bw_dirs[$i]}/${f%.*}.nucl.bw \
#     --binSize=1 \
#     --normalizeTo1x $tair10_size \
#     --ignoreDuplicates \
#     --numberOfProcessors=$threads \
#     --ignoreForNormalization Ath_chrm Ath_chrc
#   echo "Creating normalized bigwig files... - Done"
#
#   #TODO: implement removal of the old per chromosome files.
#   #rm -vr ${fseq_combined[@]}
#   #rm -vr ${fseq_subnucl[@]}
#   #rm -vr ${fseq_nucl[@]}
# done
#
# if [ $clean -eq 1 ]; then
#   echo "Cleaning up..."
#   for dir in "${bed_dirs[@]}"
#   do
#     rm -rv $dir
#   done
#   echo "Cleaning up... - Done"
# fi
#
#

#
#
# if [ $macs2 -eq 1 ]; then
#   bt_2_files=$sample_dir/bowtie2
#   bt_2_split=$bt_2_files_split/split_bam
#   bt_2_bed=$bt_1_files/bed_files
#   bt_2_peaks=$bt_2_files/peak_calling
#
#   macs2 callpeak \
#       -t $treatment \
#       -c $control \
#       -f BAMPE \
#       -g $effective_genome_size \
#       -n $name \
#       -B \
#       -m 5 50 \
#       --fix-bimodal \
#       -q 0.05 \
#       --call-summits \
#
#       effective_genome_size=1.2e8
#
#
#   #TODO: load moduleß
#
#   mkdir -p $bt_2_split
#   mkdir -p $bt_2_peaks
#   mkdir -p $bt_2_bed
#
#   aligner_dirs+=($bt_2_files)
#   split_dirs+=($bt_2_split)
#   bed_dirs+=($bt_2_bed)
#   peaks_dirs+=($bt_2_peaks)
# fi
#
#
#
#
#
# #TODO: implement clean loop that removes or zips the bed files.
