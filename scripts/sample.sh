#!/bin/bash 
#PBS -l nodes=1:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N tau_Ek=1e7
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/pla1/scripts
EXEC=/home/jimmy.meza/pla1/scripts/ifc.x
cd $THIS_DIR
PY_SCRIPT=rr_hydra.py
INPUTS="$PY_SCRIPT dummy_str"
NPROCS="-np 24"

# corrida
mpirun $NPROCS $EXEC $INPUTS
#EOF
