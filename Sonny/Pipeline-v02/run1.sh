#!/bin/bash
#PBS -l nodes=1:ppn=1
 
parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/sfdejong/04_01_2019/jobs1.txt)
parameterArray=($parameters)
 
Filenames=${parameterArray[0]}
Cas=${parameterArray[1]}
PAM=${parameterArray[2]}
Rev=${parameterArray[3]}
python2  /home/sfdejong/04_01_2019/genomicpipeline.py  $Filenames $Cas $PAM $Rev