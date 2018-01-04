#!/bin/bash 

export REPO=$HOME/pla1__Lc          # Git work-tree on lamp.iafe.uba.ar
EXEC=$REPO/scripts/rr_hydra.py
NPROCS="-np 24"                 # number of processors


CONDAENV=work2      # conda environment where we work
echo -e "\n [+] Switching to $CONDAENV environment in Anaconda framework."
# load the mpirun we need:
source activate $CONDAENV || exit 1

# check we are in a committed state
if [[ $(git status | grep "nothing to commit") =~ "nothing to commit" ]]; then
    echo " [+] clean Git repo."
else
    echo -e "\n [-] dirty Git repo. Commit input changes && commit!\n" && exit 1
fi

ro="0.90"  #$(printf "%.2f" 0.5)                            # [AU] heliodistance
ofname="$REPO/out/newLc/r_${ro}__$(git rev-parse HEAD).h5"  # output HDF5

# check if there's already a .h5 with this name!
[[ -f $ofname ]] && {
    echo -e " [-] HDF5 file already exists! $(ls -lh $ofname)\n CHECK before running!!\n";
    exit 1;
}
LOGFILE=$REPO/out/newLc/$(basename $ofname).log             # runlog

#-- report on screen
echo ""
echo -e " [+] Starting run @ $(date)"
echo -e " [+] runlog: $LOGFILE"
echo -e " [*] output: $ofname"
echo ""

# NOTE: the $(which mpirun) belongs to Anaconda!
mpirun $NPROCS $EXEC -- \
    --fname_out $ofname \
    --ro $ro \
    --Nm_slab 128 \
    --Nm_2d 128 \
    --tmax 4e4 \
    --eps 4.64e-5 \
    --sigma 0.3 \
    --Nph 8 \
    --Nth 16 \
    --trails \
    --taubd 0.8 1.1 1.4 1.6 \
    > $LOGFILE 2>&1

echo -e "\n [*] FINISHED WITH STATUS: $? @ $(date)\n"

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
