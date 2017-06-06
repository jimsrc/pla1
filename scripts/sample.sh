#!/bin/bash 

#THIS_DIR=/home/jimmy.meza/pla1/scripts
export PLA1=/home/masiasmj/pla1
THIS_DIR=$PLA1/scripts
EXEC=${THIS_DIR}/rr_hydra.py
cd $THIS_DIR
PY_SCRIPT=rr_hydra.py
INPUTS="$PY_SCRIPT dummy_str"
NPROCS="-np 4"

# corrida
LOGFILE=${THIS_DIR}/test.log

# load the mpirun we need:
source activate work2

#mpirun $NPROCS $EXEC $INPUTS > $LOGFILE 2>&1
mpirun $NPROCS $EXEC > $LOGFILE 2>&1

# NOTE:
# to kill the processes, you may use this (carefull 
# with kill other processes NO BELONGIN TO THIS RUN):
# kill $(ps -ef | grep -i <pattern> | grep -v "grep\|vim " | awk '{print $2}')


#EOF
