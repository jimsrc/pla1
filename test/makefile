# we should test:
# * compiling w/ different macros in control.h
# * run w/ the different macros (just to check to it runs OK)
# * validate w/ known quantivative results (e.g. Giacalone and Shalchi)
# * generate figs to check that we obtain "reasonable" stuff
SHELL=/bin/bash
fControl=../src/control.h
MACRO__trail="WATCH_TRAIL"

default:
	make help

test_Bmodel:
	ipython test_Bmodel.py

test_compile__watch_trails:
	# check if WATCH_TRAIL macro is activated
	#NotActivated := $(shell grep "//#define WATCH_TRAIL" ../src/control_copy.h > /dev/null; echo $$?)
	@grep "//#define ${MACRO__trail}" ${fControl} > /dev/null; \
		if [ $$? -eq 0 ]; then \
			echo " [*] activating macro..."; \
			sed -Ei "s/\/\/#define ${MACRO__trail}/#define ${MACRO__trail}/g" ${fControl}; \
		else \
			echo " macro already activated."; \
		fi
	# now compile
	@cd ../src && make clean && make \
		&& echo -e "\n [+] compiled ok.\n" \
		|| { \
			echo -e "\n [-] compilation error.\n"; \
			exit 1; \
		}
	# run test
	#

#EOF
