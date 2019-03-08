#!/bin/bash
#PBS -l nodes=1:ppn=20
/home/dddekker/miniconda3/bin/python3.7  /home/dddekker/BEP/Pipeline_fit_Finkelstein.py 'init_limit_general_energies_v2' 'monitor.txt' 'fit_results.txt' 'init_monitor.txt' 1

