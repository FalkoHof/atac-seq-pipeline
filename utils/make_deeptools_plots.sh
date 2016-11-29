#!/bin/bash
############
# Script that takes an gz numpy array from deeptools and invokes different
# plotting functions
# Parameters:
# 1 npz file
############
# make sure the following modules in your path
# ml deepTools/2.2.4-foss-2015a-Python-2.7.9
# ml OpenSSL/1.0.1p-foss-2015a

dt_corr_f=$1
out_dir=$2



#1. make all the correlation plots
params_corMethod=("spearman" "pearson")
params_whatToPlot=("heatmap" "scatterplot")
params_correction=("--skipZeros" "--removeOutliers")

params_plot=()
file_names=()

#nested loops to generate different combinations of parameters and file_names
params=""
name=""
for p1 in ${params_corMethod[@]}; do
  name=""
  params=""
  name=$name"$p1"
  params="$params--corMethod $p1"
  #some temp variabels that allow resetting to the previous state
  temp_name=$name
  temp_param=$params
  for p2 in ${params_whatToPlot[@]}; do
     params="$params --whatToPlot $p2"
     name=$name"_$p2"
     #add to array containing all the parameters
     params_plot+=("$params")
     file_names+=("$name")
     #echo "$params"
     #echo "$name"
     for p3 in ${params_correction[@]}; do
       params="$params $p3"
       name=$name"_${p3#--}"
       params_plot+=("$params")
       file_names+=("$name")
      #  echo "$params"
      #  echo "$name"
     done
     #reset to pre-loop conditions
     params=$temp_param
     name=$temp_name
  done
done

for i in "${!params_plot[@]}"; do
  plotCorrelation \
    ${params_plot[i]}\
    --corData $dt_corr_f \
    --plotFile $out_dir/${dt_corr_f%.*}_${file_names[i]}.pdf \
    --plotNumbers
done

#2. plot pca

plotPCA -in $dt_corr_f \
  -o ${dt_corr_f%.*}_pca.pdf
