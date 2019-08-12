#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -m e
#PBS -M diewertje@kpnplanet.nl
parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/dddekker/BEP/12_6_2019_cas12_conc8/jobs_12_6_2019_cas12_conc8_1_25.txt)
parameterArray=($parameters) 

ModelID=${parameterArray[0]}
OutputMonitor=${parameterArray[1]}
OutputResults=${parameterArray[2]}
OutputInitMonitor=${parameterArray[3]}
UseHPC05=${parameterArray[4]}

python2  /home/dddekker/BEP/Pipeline_fit_Finkelstein_cas12.py $ModelID $OutputMonitor $OutputResults $OutputInitMonitor $UseHPC05


