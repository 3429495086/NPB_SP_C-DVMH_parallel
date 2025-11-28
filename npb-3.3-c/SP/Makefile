SHELL=/bin/sh
BENCHMARK=sp
BENCHMARKU=SP

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = header.h npbparams.h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h

# List of files which can be merged into a single AST.
MERGE_SRC = sp.c initialize.c exact_solution.c exact_rhs.c \
            set_constants.c adi.c rhs.c x_solve.c ninvr.c \
            y_solve.c pinvr.c z_solve.c tzetar.c add.c txinvr.c \
            error.c verify.c \
            ${COMMON}/print_results.c ${COMMON}/wtime.c \
            #lhsx.c lhsy.c lhsz.c


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

sp.ast sp.o:                         sp.c  header.h npbparams.h
initialize.ast initialize.o:         initialize.c  header.h npbparams.h
exact_solution.ast exact_solution.o: exact_solution.c  header.h npbparams.h
exact_rhs.ast exact_rhs.o:           exact_rhs.c  header.h npbparams.h
set_constants.ast set_constants.o:   set_constants.c  header.h npbparams.h
adi.ast adi.o:                       adi.c  header.h npbparams.h
rhs.ast rhs.o:                       rhs.c  header.h npbparams.h
#lhsx.ast lhsx.o:                    lhsx.c  header.h npbparams.h
#lhsy.ast lhsy.o:                    lhsy.c  header.h npbparams.h
#lhsz.ast lhsz.o:                    lhsz.c  header.h npbparams.h
x_solve.ast x_solve.o:               x_solve.c  header.h npbparams.h
ninvr.ast ninvr.o:                   ninvr.c  header.h npbparams.h
y_solve.ast y_solve.o:               y_solve.c  header.h npbparams.h
pinvr.ast pinvr.o:                   pinvr.c  header.h npbparams.h
z_solve.ast z_solve.o:               z_solve.c  header.h npbparams.h
tzetar.ast tzetar.o:                 tzetar.c  header.h npbparams.h
add.ast add.o:                       add.c  header.h npbparams.h
txinvr.ast txinvr.o:                 txinvr.c  header.h npbparams.h
error.ast error.o:                   error.c  header.h npbparams.h
verify.ast verify.o:                 verify.c  header.h npbparams.h

clean:
	- rm -f *.o *~ mputil* *.ast *.ll
	- rm -f npbparams.h core
