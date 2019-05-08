#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -m e
#PBS -M diewertje@kpnplanet.nl
parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/dddekker/BEP/2_5_2019_var_rates_conc3/jobs_2_5_2019_var_rates_conc3_1_25.txt)
parameterArray=($parameters) 

ModelID=${parameterArray[0]}
OutputMonitor=${parameterArray[1]}
OutputResults=${parameterArray[2]}
OutputInitMonitor=${parameterArray[3]}
UseHPC05=${parameterArray[4]}

python2  /home/dddekker/BEP/Pipeline_fit_Finkelstein_2.py $ModelID $OutputMonitor $OutputResults $OutputInitMonitor $UseHPC05


