include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard eval_pslFilter.py eval_pslSwitchQueryTarget.py))

targets = ${progs:%=${binPath}/%}
ifneq ($(wildcard ${kentLibWeb}),)
  targets += ${binPath}/eval_pslMapCoor ${binPath}/eval_pslConvertCoor
  # FIXME: tmp hack, both sonlib commonC and kent libraries have logInfo and logDebug functions
  # kent ones are not actually used, but get dragged in.  As a tmp work around, we force
  # the sonlib ones to be brought in first by the linker.
  jkLibHack = -Wl,-ulogDebug ${libPath}/sonLib.a 
endif

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/eval_pslMapCoor : eval_pslMapCoor.c
	${cxx} ${cflags} eval_pslMapCoor.c -lpthread -lz -I ${kentInc} -o ${binPath}/eval_pslMapCoor ${jkLibHack} ${kentLibWeb} 

${binPath}/eval_pslConvertCoor : eval_pslConvertCoor.c
	${cxx} ${cflags} eval_pslConvertCoor.c -lpthread -lz -I ${kentInc} -o ${binPath}/eval_pslConvertCoor ${jkLibHack} ${kentLibWeb} 

clean :
	rm -f ${binPath}/eval_pslMapCoor
	rm -f ${progs:%=${binPath}/%}
