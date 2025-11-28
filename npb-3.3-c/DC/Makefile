SHELL=/bin/sh
BENCHMARK=dc
BENCHMARKU=DC

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h

# List of files which can be merged into a single AST.
MERGE_SRC = adc.c dc.c extbuild.c rbt.c jobcntl.c \
            ${COMMON}/c_print_results.c ${COMMON}/wtime.c \

# List of files which can not be imported by TSAR,
# for example due to symbol name redifinition.
SINGLE_SRC = ${COMMON}/c_timers.c

#------------------------------------------------------------------------------#

include ../sys/make.vars

# npbparams.h is provided for backward compatibility with NPB compilation
# header.h: npbparams.h

${PROGRAM}: config
	${MAKE} $(ACTION);

include ../sys/make.tsar

adc.ast adc.o:           adc.c npbparams.h
dc.ast dc.o:             dc.c adcc.h adc.h macrodef.h npbparams.h
extbuild.ast extbuild.o: extbuild.c adcc.h adc.h macrodef.h npbparams.h
rbt.ast rbt.o:           rbt.c adcc.h adc.h rbt.h macrodef.h npbparams.h
jobcntl.ast jobcntl.o:   jobcntl.c adcc.h adc.h macrodef.h npbparams.h

clean:
	- rm -f *.o *.ast *.ll
	- rm -f npbparams.h core
	- rm -f {../,}ADC.{logf,view,dat,viewsz,groupby,chunks}.* 

