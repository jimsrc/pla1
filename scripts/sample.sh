#!/bin/bash 

export PLA1=$HOME/pla1          # repo-dir on lamp.iafe.uba.ar
THIS_DIR=$PLA1/scripts
EXEC=${THIS_DIR}/rr_hydra.py
NPROCS="-np 24"

# runlog
LOGFILE=${THIS_DIR}/run_r0.9_NmS128_Nm2d128_eps4.64e-6.log

CONDAENV=work2      # conda environment where we work
echo -e "\n [+] Switching to $CONDAENV environment in Anaconda framework."
# load the mpirun we need:
source activate $CONDAENV

echo -e "\n [+] Starting run @ $(date)"
echo -e " [+] runlog: $LOGFILE\n"
# NOTE: the $(which mpirun) belongs to Anaconda!
mpirun $NPROCS $EXEC > $LOGFILE 2>&1

echo -e "\n [+] FINISHED WITH STATUS: $? @ $(date)\n"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# NOTE:
# To kill the processes, you may use this. 
# WARNING: Careful because you might have other processes with the 
# same <pattern> (containing the run script name) that DON'T BELONG 
# TO THIS RUN.
# Normally that's not a problem if you are using the script <pattern>
# for this run only.
# kill $(ps -ef | grep -i <pattern> | grep -v "grep\|vim " | awk '{print $2}')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#EOF
