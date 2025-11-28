SHELL=/bin/sh
BENCHMARK=ep
BENCHMARKU=EP

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h ${COMMON}/randdp.h

# List of files which can be merged into a single AST.
MERGE_SRC = ep.c \
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

ep.ast ep.o: ep.c npbparams.h

clean:
	- rm -f *.o *~ *.ast *.ll
	- rm -f npbparams.h core

