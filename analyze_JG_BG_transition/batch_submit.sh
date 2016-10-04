#!/bin/bash

#BSUB -P P93300601
#BSUB -W 12:00
#BSUB -n 1
#BSUB -q geyser

module load python
module load all-python-libs

python analyze_JG_BG_transition.py -j JG_iteration_5 
