#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N run.sh
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/pla1/scripts
EXEC=/home/jimmy.meza/pla1/scripts/ifc.x
cd $THIS_DIR
PY_SCRIPT=rr_hydra.py
INPUTS="$PY_SCRIPT dummy_str"
NPROCS="-np 48"

# corrida
mpirun $NPROCS $EXEC $INPUTS
#EOF
