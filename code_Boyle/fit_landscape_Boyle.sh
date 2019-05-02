#!/bin/bash
#PBS -l nodes=1:ppn=20

parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/mklein1/Energy_Landscape_dCas9/30_4_2019/jobs_30_4_2019_26_38.txt)
parameterArray=($parameters)

ReplicaID=${parameterArray[0]}
ModelID=${parameterArray[1]}
OutputResults=${parameterArray[2]}
OutputMonitor=${parameterArray[3]}
OutputInitMonitor=${parameterArray[4]}
UseHPC05=${parameterArray[5]}
python2  /home/mklein1/Energy_Landscape_dCas9/fit_landscape_Boyle.py $ReplicaID $ModelID $OutputResults $OutputMonitor $OutputInitMonitor $UseHPC05


