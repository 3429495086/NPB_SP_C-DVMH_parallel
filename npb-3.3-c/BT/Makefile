SHELL=/bin/sh
BENCHMARK=bt
BENCHMARKU=BT
VEC=
ACTION=exec

include ../config/make.def

include ../sys/make.common

#------------------- Lists of sources and headers -----------------------------#

HEAD = header.h npbparams.h work_lhs$(VEC).h \
       ${COMMON}/print_results.h ${COMMON}/timers.h \
       ${COMMON}/wtime.h ${COMMON}/type.h

# List of files which can be merged into a single AST.
MERGE_SRC = exact_solution.c exact_rhs.c \
          set_constants.c adi.c  rhs.c \
          x_solve$(VEC).c y_solve$(VEC).c solve_subs.c \
          z_solve$(VEC).c add.c error.c initialize.c \
          bt.c verify.c \
          ${COMMON}/print_results.c ${COMMON}/wtime.c

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
	@if [ x$(VERSION) = xvec ] ; then	\
		${MAKE} VEC=_vec $(ACTION);		\
	elif [ x$(VERSION) = xVEC ] ; then	\
		${MAKE} VEC=_vec $(ACTION);		\
	else					\
		${MAKE} $(ACTION);			\
	fi

include ../sys/make.tsar

bt.ast bt.o:                         bt.c  header.h npbparams.h
initialize.ast initialize.o:         initialize.c  header.h npbparams.h
exact_solution.ast exact_solution.o: exact_solution.c  header.h npbparams.h
exact_rhs.ast exact_rhs.o:           exact_rhs.c  header.h npbparams.h
set_constants.ast set_constants.o:   set_constants.c  header.h npbparams.h
ast.ast adi.o:                       adi.c  header.h npbparams.h
rhs.ast rhs.o:                       rhs.c  header.h npbparams.h
x_solve$(VEC).ast x_solve$(VEC).o:   x_solve$(VEC).c  header.h work_lhs$(VEC).h npbparams.h
y_solve$(VEC).ast y_solve$(VEC).o:   y_solve$(VEC).c  header.h work_lhs$(VEC).h npbparams.h
z_solve$(VEC).ast z_solve$(VEC).o:   z_solve$(VEC).c  header.h work_lhs$(VEC).h npbparams.h
solve_subs.ast solve_subs.o:         solve_subs.c  npbparams.h
add.ast add.o:                       add.c  header.h npbparams.h
error.ast error.o:                   error.c  header.h npbparams.h
verify.ast verify.o:                 verify.c  header.h npbparams.h

clean:
	- rm -f *.o *~ mputil* *.ast *.ll
	- rm -f  npbparams.h core
