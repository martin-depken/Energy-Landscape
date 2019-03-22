#!/bin/bash
#PBS -l nodes=1:ppn=20

parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/mklein1/Energy_Landscape_dCas9/22_3_2019/jobs_22_3_2019_1_50.txt)
parameterArray=($parameters)

ReplicaID=${parameterArray[0]}
ModelID=${parameterArray[1]}
OutputMonitor=${parameterArray[2]}
OutputResults=${parameterArray[3]}
OutputInitMonitor=${parameterArray[4]}
UseHPC05=${parameterArray[5]}
python2  /home/mklein1/Energy_Landscape_dCas9/fit_landscape_Boyle.py $ReplicaID $ModelID $OutputMonitor $OutputResults $OutputInitMonitor $UseHPC05


