SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = globals.h npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h ${COMMON}/randdp.h

# List of files which can be merged into a single AST.
MERGE_SRC = mg.o \
            ${COMMON}/print_results.c ${COMMON}/wtime.c \
			${COMMON}/${RAND}.c

# List of files which can not be imported by TSAR,
# for example due to symbol name redifinition.
SINGLE_SRC = ${COMMON}/c_timers.c

#------------------------------------------------------------------------------#

include ../sys/make.vars

${PROGRAM}: config
	${MAKE} $(ACTION);

include ../sys/make.tsar

mg.ast mg.o:		mg.c globals.h npbparams.h

clean:
	- rm -f *.o *~ *.ast *.ll
	- rm -f npbparams.h core

