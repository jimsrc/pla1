#!/bin/bash 

export PLA1=$HOME/pla1
THIS_DIR=$PLA1/scripts
EXEC=${THIS_DIR}/rr_hydra.py
cd $THIS_DIR
PY_SCRIPT=rr_hydra.py
NPROCS="-np 32"

# corrida
LOGFILE=${THIS_DIR}/run_r0.9_NmS256_Nm2d256.log

# load the mpirun we need:
source activate work2

# note that the $(which mpirun) belongs to $ANACONDA
mpirun $NPROCS $EXEC > $LOGFILE 2>&1

# NOTE:
# to kill the processes, you may use this (carefull 
# with kill other processes NO BELONGIN TO THIS RUN):
# kill $(ps -ef | grep -i <pattern> | grep -v "grep\|vim " | awk '{print $2}')


#EOF
