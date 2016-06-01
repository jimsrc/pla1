#!/bin/bash

function get_sources(){
    local l=(
    control.h
    #main_mpi.cc
    general.cc
    general.h
    defs_turb.cc
    defs_turb.h
    funcs.cc
    funcs.h
    nr3.h
    stepperbs.cc
    stepperbs.h
    odeintt.h
    odeintt.cc
    )
    echo "${l[@]}"
}


function check_diffs(){
    echo -e "\n ---> comparing with:\n " ${dir_comp}"\n"
    for fname in ${fname_list[@]}; do
        fname2=${dir_comp}/$fname
        diff_result=$(diff $fname $fname2) # compare text
        #if cmp -s $fname $fname2; then
        if [ "$diff_result" == "" ]; then
            cc="\033[32m" # synced
            ok=1
        else
            cc="\033[31m" # not synced
            ok=0
        fi
        echo -e $cc "--> [$ok] $fname"
    done
}

#get_sources
fname_list=$(get_sources)
echo $fname_list
dir_comp=$PLAS/src # directorio con el q vamos a comparar

#--------------------- execute
if [[ -n $1 ]]; then
    if [[ "$1"=="sync" ]]; then
        # copy file from ${dir_comp} to $PWD
        /bin/cp -p ${dir_comp}/$2 $PWD/.
        echo " ----> SYNCED WITH: " ${dir_comp}/$2
    fi
fi

check_diffs
#-----------------------------

#EOF
