#include ../../../include.mk
binPath = ../bin
#libPath = ../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard cactus_writeRunCactusMakefile.py cactus_mafToReferenceSeq.py cactus_writeExperimentXml.py eval_mafMapCoor.py eval_bedMapCoor.py eval_MAFToBed.py eval_MAFFilter.py mhc_mafExtract.py mhc_mafToGenomeSeqs.py mhc_writeRepeatMaskingMakefiles.py))

modules = pSLs
.PHONY: %.all

targets = ${progs:%=${binPath}/%} ${modules:%=%.all}

all: ${targets} ${binPath}/cactus_print

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

#${binPath}/cactus_print : *.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
#	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_print cactus_print.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}


%.all:
	cd $* && make all

clean :
	rm -f ${progs:%=${binPath}/%}
	#rm -f ${binPath}/cactus_print
	cd ${modules:%=%} && make clean
