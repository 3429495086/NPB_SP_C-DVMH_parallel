SHELL=/bin/sh
BENCHMARK=ua
BENCHMARKU=UA

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = header.h npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h

# List of files which can be merged into a single AST.
MERGE_SRC = ua.c convect.c diffuse.c adapt.c move.c mason.c \
            precond.c utils.c transfer.c verify.c setup.c\
            ${COMMON}/print_results.c ${COMMON}/wtime.c \

# List of files which can not be imported by TSAR,
# for example due to symbol name redifinition.
SINGLE_SRC = ${COMMON}/c_timers.c

#------------------------------------------------------------------------------#

include ../sys/make.vars

# npbparams.h is included by header.h
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to 
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect). 

# header.h: npbparams.h

${PROGRAM}: config
	${MAKE} $(ACTION);

include ../sys/make.tsar

ua.ast ua.o:             ua.c       header.h npbparams.h
setup.ast setup.o:       setup.c    header.h npbparams.h
convect.ast convect.o:   convect.c  header.h npbparams.h
adapt.ast adapt.o:       adapt.c    header.h npbparams.h
move.ast move.o:         move.c     header.h npbparams.h
diffuse.ast diffuse.o:   diffuse.c  header.h npbparams.h
mason.ast mason.o:       mason.c    header.h npbparams.h
precond.ast precond.o:   precond.c  header.h npbparams.h
transfer.ast transfer.o: transfer.c header.h npbparams.h
utils.ast utils.o:       utils.c    header.h npbparams.h
verify.ast verify.o:     verify.c   header.h npbparams.h

clean:
	- rm -f *.o *~ mputil* *.ast *.ll
	- rm -f npbparams.h core
