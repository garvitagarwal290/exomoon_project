#!/bin/bash
###### Job name ###### 
#PBS -N 918_SNR_planet_only_pandora_modeling
###### Output files ######
#PBS -o /tiara/home/gagarwal/modeling_results/output_918_SNR_july13_nores_planet_only
#PBS -e /tiara/home/gagarwal/modeling_results/error_918_SNR_july13_nores_planet_only
###### Number of nodes and cores ######
#PBS -l walltime=96:00:00
#PBS -l select=1:ncpus=1
###### Queue name ######
#PBS -q serial
###### Sends mail to yourself when the job begins and ends ######
#PBS -M garvit.agarwal@students.iiserpune.ac.in
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd /tiara/home/gagarwal/modeling_results/
source /tiara/home/gagarwal/miniconda3/etc/profile.d/conda.sh
conda activate exomoon_project
python /tiara/home/gagarwal/code/modeling_SNR.py 56 planet_only 1.0 resume
