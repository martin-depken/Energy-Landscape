#!/bin/bash
#PBS -l nodes=1:ppn=20
parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/mklein1/Diewertje/1_3_2019/jobs_1_3_2019_1_50.txt)
parameterArray=($parameters)

ModelID=${parameterArray[0]}
OutputMonitor=${parameterArray[1]}
OutputResults=${parameterArray[2]}
OutputInitMonitor=${parameterArray[3]}
UseHPC05=${parameterArray[4]}
python2  /home/mklein1/Diewertje/Pipeline_fit_Finkelstein.py $ModelID $OutputMonitor $OutputResults $OutputInitMonitor $UseHPC05

