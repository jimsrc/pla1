SHELL=/bin/bash
.PHONY: default help object executable all clean

#SOURCE_C = $(wildcard *.cc)
SOURCE_C = general.cc defs_turb.cc funcs.cc odeintt.cc stepperbs.cc
OBJECTS_C = $(patsubst %.cc, %.o, $(SOURCE_C))
OBJECTS_C_GPROF = $(patsubst %.cc, %_gp.o, $(SOURCE_C))
CFLAGS = -O3
CXX = /usr/local/bin/mpic++
#EXE = ../CRs.diff.x
EXE = ./main.x
# gprof stuff
EXE_WPROF=./main_wgprof.x 	# profiling / debugging
OPT_GP=-O0 #-O3 #-O0 	# for use w/ gprof

cython:
	python setup.py build_ext --inplace

all: $(OBJECTS_C)
	${CXX} -o $(EXE)  $^  ${CFLAGS}

objects: $(OBJECTS_C)

%.o: %.cc
	${CXX} ${CFLAGS} -c $^ 

%_gp.o: %.cc
	$(CXX) $(OPT_GP) -c $^ -o $@

gprof_exe: $(OBJECTS_C_GPROF)
	${CXX} -pg ${OPT_GP} $^ -o $(EXE_WPROF)

clean:
	rm *.[xo] *.so *.so
