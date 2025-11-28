SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = globals.h npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h ${COMMON}/randdp.h

# List of files which can be merged into a single AST.
MERGE_SRC = appft.c auxfnct.c fft3d.c mainft.c verify.c \
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

appft.ast appft.o:		appft.c  global.h npbparams.h
auxfnct.ast auxfnct.o:	auxfnct.c  global.h npbparams.h
fft3d.ast fft3d.o:		fft3d.c  global.h npbparams.h
mainft.ast mainft.o:	mainft.c  global.h npbparams.h
verify.ast verify.o:	verify.c  global.h npbparams.h

clean:
	- rm -f *.o *~ mputil* *.ast *.ll
	- rm -f ft npbparams.h core
