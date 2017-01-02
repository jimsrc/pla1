#!/bin/bash 
#PBS -l nodes=1:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N pla1_r=2AU
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/pla1/scripts
EXEC=/home/jimmy.meza/pla1/scripts/ifc.x
cd $THIS_DIR
PY_SCRIPT=rr_hydra.py
INPUTS="$PY_SCRIPT dummy_str"
NPROCS="-np 24"

# corrida
mpirun $NPROCS $EXEC $INPUTS 1> ${THIS_DIR}/mon1.log 2> ${THIS_DIR}/mon2.log
#EOF
