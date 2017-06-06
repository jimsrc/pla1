#!/bin/bash

CONTROLFILE=./control.h
CYTHONCFILE=./macros.pyx
# build
# NOTE: this deletes those lines that:
# - start with "//"
# - has an empty line
# - define a function macro such as "#define NORM(x) ..."
# - define the "PRIVATE_OR_PUBLIC" 'cause it's meaningless to us
# - define the "CONTROL" macro 'cause that's just a guard,
# and set to "0" those definitions that are commented (i.e. start with "//#define").
# So this will translate the $CONTROLFILE definitions into lines like 
#  DEF SOMETHNG = "algo"
#  DEF OTHER = "0"
# That's all.

#--- header
echo "\
\"\"\"
 preprocessor for Cython, taken out from 
 the $CONTROLFILE definitions. This file was
 generated automatically by ./build_cython_macros.sh
 DATE: $(date)
\"\"\"\
" > $CYTHONCFILE

#--- contents
cat $CONTROLFILE | awk '{if($1=="#define" && $1!="") {print $2, $3} else if($1=="//#define") {print $2, 0} }' | awk '{if($1 !~ /\(/ && $1!="PRIVATE_OR_PUBLIC" && $1!="CONTROL") print "DEF " $1 " = \"" $2"\""}' >> $CYTHONCFILE

#--- EOF
echo -e "\n#EOF" >> $CYTHONCFILE


#EOF
