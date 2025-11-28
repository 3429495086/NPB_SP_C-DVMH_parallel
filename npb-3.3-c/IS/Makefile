SHELL=/bin/sh
BENCHMARK=is
BENCHMARKU=IS

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h

# List of files which can be merged into a single AST.
MERGE_SRC = is.c \
            ${COMMON}/c_print_results.c ${COMMON}/wtime.c \

# List of files which can not be imported by TSAR,
# for example due to symbol name redifinition.
SINGLE_SRC = ${COMMON}/c_timers.c

#------------------------------------------------------------------------------#

include ../sys/make.vars

${PROGRAM}: config
	${MAKE} $(ACTION);

include ../sys/make.tsar

is.ast is.o:             is.c  npbparams.h


clean:
	- rm -f *.o *~ mputil* *.ast *.ll
	- rm -f npbparams.h core
