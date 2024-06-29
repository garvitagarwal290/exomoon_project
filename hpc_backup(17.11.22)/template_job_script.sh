#!/bin/bash
###### Job name ###### 
#PBS -N pandora_modeling
###### Output files ######
#PBS -o /tiara/home/gagarwal/modeling_results/out_file
#PBS -e /tiara/home/gagarwal/modeling_results/err_file
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
python /tiara/home/gagarwal/code/modeling.py 74 planet_only
