#include ../../../include.mk
binPath = ../bin
#libPath = ../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard runVelvet.py mhc_runVelvet.py))

#modules = pSLs
modules = 
.PHONY: %.all

targets = ${progs:%=${binPath}/%} ${modules:%=%.all}

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@


%.all:
	cd $* && make all

clean :
	rm -f ${progs:%=${binPath}/%}
	cd ${modules:%=%} && make clean
