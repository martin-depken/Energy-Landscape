#!/bin/bash
#PBS -l nodes=1:ppn=1
 
parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/sfdejong/11_03_2019/jobs.txt)
parameterArray=($parameters)
 
Filenames=${parameterArray[0]}
python3  genomicpipeline.py  $Filenames